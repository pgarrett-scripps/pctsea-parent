#!/usr/bin/env python3
"""Prepare protein spectral-count files as canonical pCTSEA gene queries.

This is a validation data tool, not a runtime dependency of pCTSEA. It joins
UniProt accessions to gene names and protein lengths from a pinned FASTA file.
It writes both raw gene-level spectrum counts and NSAF-style values.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


GENE_PATTERN = re.compile(r"(?:^|\s)GN=([^\s]+)")


@dataclass(frozen=True)
class Protein:
    accession: str
    gene: str
    length: int


@dataclass(frozen=True)
class Observation:
    locus: str
    spectrum_count: float
    accession: str
    gene: str
    length: int


@dataclass
class FileSummary:
    file: str
    input_rows: int = 0
    target_species_rows: int = 0
    foreign_species_rows: int = 0
    contaminant_rows: int = 0
    decoy_rows: int = 0
    unsupported_rows: int = 0
    missing_accession_rows: int = 0
    mapped_rows: int = 0
    output_genes: int = 0
    collapsed_rows: int = 0
    total_spectrum_count: float = 0.0
    target_spectrum_count: float = 0.0
    mapped_spectrum_count: float = 0.0

    def as_row(self) -> list[str]:
        row_rate = self.mapped_rows / self.target_species_rows if self.target_species_rows else 0.0
        count_rate = (
            self.mapped_spectrum_count / self.target_spectrum_count
            if self.target_spectrum_count
            else 0.0
        )
        return [
            self.file,
            str(self.input_rows),
            str(self.target_species_rows),
            str(self.foreign_species_rows),
            str(self.contaminant_rows),
            str(self.decoy_rows),
            str(self.unsupported_rows),
            str(self.missing_accession_rows),
            str(self.mapped_rows),
            str(self.output_genes),
            str(self.collapsed_rows),
            format_number(self.total_spectrum_count),
            format_number(self.target_spectrum_count),
            format_number(self.mapped_spectrum_count),
            f"{row_rate:.6f}",
            f"{count_rate:.6f}",
        ]


SUMMARY_HEADER = [
    "file",
    "input_rows",
    "target_species_rows",
    "foreign_species_rows",
    "contaminant_rows",
    "decoy_rows",
    "unsupported_rows",
    "missing_accession_rows",
    "mapped_rows",
    "output_genes",
    "collapsed_rows",
    "total_spectrum_count",
    "target_spectrum_count",
    "mapped_spectrum_count",
    "mapped_target_row_fraction",
    "mapped_spectrum_count_fraction",
]


def fasta_records(path: Path) -> Iterable[tuple[str, str]]:
    header: str | None = None
    sequence: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence)
                header = line[1:]
                sequence = []
            elif header is None:
                raise ValueError(f"{path}: sequence occurs before the first FASTA header")
            else:
                sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def load_uniprot_fasta(path: Path) -> dict[str, Protein]:
    proteins: dict[str, Protein] = {}
    for header, sequence in fasta_records(path):
        identifier = header.split(maxsplit=1)[0]
        parts = identifier.split("|")
        if len(parts) < 3 or parts[0] not in {"sp", "tr"}:
            continue
        gene_match = GENE_PATTERN.search(header)
        if gene_match is None or not sequence:
            continue
        accession = parts[1].upper()
        proteins[accession] = Protein(
            accession=accession,
            gene=gene_match.group(1).upper(),
            length=len(sequence),
        )
    if not proteins:
        raise ValueError(f"{path}: no UniProt records with gene names were found")
    return proteins


def parse_locus(locus: str) -> tuple[str, str] | None:
    parts = locus.strip().split("|")
    if len(parts) < 3 or parts[0] not in {"sp", "tr"}:
        return None
    accession = parts[1].strip().upper()
    entry_name = parts[2].strip().upper()
    if not accession or not entry_name:
        return None
    return accession, entry_name


def read_observations(
    path: Path,
    proteins: dict[str, Protein],
    target_suffix: str,
) -> tuple[list[Observation], FileSummary, list[list[str]]]:
    observations: list[Observation] = []
    unresolved: list[list[str]] = []
    summary = FileSummary(file=path.name)
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 2:
            raise ValueError(f"{path}: expected a two-column tab-separated file")
        if header[0].strip().lower() != "locus" or header[1].strip().lower() != "spectrum count":
            raise ValueError(f"{path}: expected 'Locus<TAB>Spectrum Count'")
        for line_number, fields in enumerate(reader, start=2):
            if not fields or not any(field.strip() for field in fields):
                continue
            if len(fields) < 2:
                raise ValueError(f"{path}:{line_number}: expected two tab-separated columns")
            locus = fields[0].strip()
            try:
                spectrum_count = float(fields[1].strip())
            except ValueError as error:
                raise ValueError(f"{path}:{line_number}: invalid spectrum count") from error
            if spectrum_count < 0:
                raise ValueError(f"{path}:{line_number}: spectrum count cannot be negative")
            summary.input_rows += 1
            summary.total_spectrum_count += spectrum_count
            if locus.lower().startswith("reverse_"):
                summary.decoy_rows += 1
                unresolved.append([path.name, str(line_number), locus, "decoy"])
                continue
            if "contaminant" in locus.lower():
                summary.contaminant_rows += 1
                unresolved.append([path.name, str(line_number), locus, "contaminant"])
                continue
            parsed = parse_locus(locus)
            if parsed is None:
                summary.unsupported_rows += 1
                unresolved.append([path.name, str(line_number), locus, "unsupported_locus"])
                continue
            accession, entry_name = parsed
            if not entry_name.endswith(f"_{target_suffix}"):
                summary.foreign_species_rows += 1
                unresolved.append([path.name, str(line_number), locus, "foreign_species"])
                continue
            summary.target_species_rows += 1
            summary.target_spectrum_count += spectrum_count
            protein = proteins.get(accession)
            if protein is None:
                summary.missing_accession_rows += 1
                unresolved.append([path.name, str(line_number), locus, "accession_not_in_fasta"])
                continue
            observations.append(
                Observation(
                    locus=locus,
                    spectrum_count=spectrum_count,
                    accession=accession,
                    gene=protein.gene,
                    length=protein.length,
                )
            )
            summary.mapped_rows += 1
            summary.mapped_spectrum_count += spectrum_count
    return observations, summary, unresolved


def aggregate_observations(
    observations: Iterable[Observation],
    aggregation: str,
    nsaf_scale: float,
) -> tuple[dict[str, float], dict[str, float], int]:
    raw_by_gene: dict[str, list[float]] = defaultdict(list)
    saf_by_gene: dict[str, list[float]] = defaultdict(list)
    row_count = 0
    for observation in observations:
        row_count += 1
        raw_by_gene[observation.gene].append(observation.spectrum_count)
        saf_by_gene[observation.gene].append(observation.spectrum_count / observation.length)
    reducer = max if aggregation == "max" else sum
    raw = {gene: reducer(values) for gene, values in raw_by_gene.items()}
    saf = {gene: reducer(values) for gene, values in saf_by_gene.items()}
    saf_total = sum(saf.values())
    nsaf = {
        gene: value / saf_total * nsaf_scale
        for gene, value in saf.items()
        if saf_total > 0
    }
    return raw, nsaf, row_count - len(raw)


def format_number(value: float) -> str:
    return f"{value:.12g}"


def write_query(path: Path, values: dict[str, float], top_genes: int | None = None) -> None:
    selected = values
    if top_genes is not None:
        selected = dict(
            sorted(values.items(), key=lambda item: (-item[1], item[0]))[:top_genes]
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene", "value"])
        for gene in sorted(selected):
            writer.writerow([gene, format_number(selected[gene])])


def input_files(input_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in input_dir.iterdir()
        if path.is_file() and path.suffix.lower() in {".tsv", ".txt"}
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prepare(args: argparse.Namespace) -> list[FileSummary]:
    proteins = load_uniprot_fasta(args.fasta)
    files = input_files(args.input_dir)
    if not files:
        raise ValueError(f"{args.input_dir}: no TSV or text files found")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summaries: list[FileSummary] = []
    unresolved: list[list[str]] = []
    for path in files:
        observations, summary, file_unresolved = read_observations(
            path,
            proteins,
            args.target_suffix.upper(),
        )
        raw, nsaf, collapsed_rows = aggregate_observations(
            observations,
            args.aggregation,
            args.nsaf_scale,
        )
        summary.output_genes = min(len(raw), args.top_genes or len(raw))
        summary.collapsed_rows = collapsed_rows
        summaries.append(summary)
        unresolved.extend(file_unresolved)
        if raw:
            output_name = f"{path.stem}.tsv"
            write_query(args.output_dir / "raw" / output_name, raw, args.top_genes)
            write_query(args.output_dir / "nsaf" / output_name, nsaf, args.top_genes)
    with (args.output_dir / "summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SUMMARY_HEADER)
        writer.writerows(summary.as_row() for summary in summaries)
    with (args.output_dir / "unresolved.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["file", "line", "locus", "reason"])
        writer.writerows(unresolved)
    with (args.output_dir / "provenance.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["property", "value"])
        writer.writerows(
            [
                ["input_dir", str(args.input_dir.resolve())],
                ["fasta", str(args.fasta.resolve())],
                ["fasta_sha256", sha256(args.fasta)],
                ["target_suffix", args.target_suffix.upper()],
                ["aggregation", args.aggregation],
                ["nsaf_scale", format_number(args.nsaf_scale)],
                ["top_genes", str(args.top_genes) if args.top_genes is not None else "all"],
            ]
        )
    return summaries


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--target-suffix",
        default="HUMAN",
        help="UniProt entry suffix to retain, default HUMAN",
    )
    parser.add_argument(
        "--aggregation",
        choices=["max", "sum"],
        default="max",
        help="How to combine multiple protein accessions for one gene, default max",
    )
    parser.add_argument(
        "--nsaf-scale",
        type=float,
        default=1_000_000.0,
        help="Scale applied after sample-wise NSAF normalization, default 1000000",
    )
    parser.add_argument(
        "--top-genes",
        type=int,
        help="Optionally retain only the highest-valued genes in each output query",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.top_genes is not None and args.top_genes < 1:
        raise ValueError("--top-genes must be at least 1")
    summaries = prepare(args)
    mapped_files = sum(summary.output_genes > 0 for summary in summaries)
    mapped_rows = sum(summary.mapped_rows for summary in summaries)
    output_genes = sum(summary.output_genes for summary in summaries)
    print(f"prepared {mapped_files}/{len(summaries)} files")
    print(f"mapped {mapped_rows} protein rows to {output_genes} per-file genes")
    print(f"reports written to {args.output_dir}")


if __name__ == "__main__":
    main()
