"""Output generation: TSV files, charts, and ZIP packaging."""

from __future__ import annotations

import logging
import math
import zipfile
from pathlib import Path

from pctsea.models import CellTypeClassification, ScoringSchema

log = logging.getLogger(__name__)

# Column definitions matching CellTypesOutputTableColumns.java
COLUMNS = [
    "cell_type",
    "num_cells_of_type",
    "num_total_cells",
    "num_cells_of_type_corr",
    "num_cells_corr",
    "hyperG_p-value",
    "log2_ratio",
    "ews",
    "norm-ews",
    "supX",
    "norm-supX",
    "empirical_p-value",
    "FDR",
    "2nd_ews",
    "2nd_supX",
    "size_a_type",
    "size_b_others",
    "Dab",
    "KS_p-value",
    "KS_p-value_BH_corrected",
    "KS_significance_level",
    "Num_Genes_Significant",
    "Umap_1",
    "Umap_2",
    "Umap_3",
    "Umap_4",
    "genes",
]


def _fmt(value: float | int | None) -> str:
    """Format a value for TSV output, matching Java's toString behavior."""
    if value is None:
        return ""
    if isinstance(value, float):
        if math.isnan(value):
            return "NaN"
        if math.isinf(value):
            return "Infinity" if value > 0 else "-Infinity"
    return str(value)


def _get_umap(ct: CellTypeClassification, dim: int) -> str:
    if ct.umap_components and dim < len(ct.umap_components):
        val = ct.umap_components[dim]
        if not math.isnan(val):
            return _fmt(val)
    return ""


def _row_values(
    ct: CellTypeClassification,
    num_total_cells: int,
    num_passing: int,
) -> list[str]:
    """Extract column values for a cell type classification row."""
    return [
        ct.name,
        _fmt(ct.num_cells_of_type),
        _fmt(num_total_cells),
        _fmt(ct.num_cells_of_type_passing_threshold),
        _fmt(num_passing),
        _fmt(ct.hypergeometric_pvalue),
        _fmt(ct.casimirs_enrichment_score),
        _fmt(ct.enrichment_score),
        _fmt(ct.get_normalized_enrichment_score()),
        _fmt(ct.supremum_x),
        _fmt(ct.normalized_supremum_x),
        _fmt(ct.enrichment_significance),
        _fmt(ct.enrichment_fdr),
        _fmt(ct.secondary_enrichment_score),
        _fmt(ct.secondary_supremum_x) if ct.secondary_supremum_x is not None else "",
        _fmt(ct.size_a),
        _fmt(ct.size_b),
        _fmt(ct.d_statistic),
        _fmt(ct.ks_pvalue),
        _fmt(ct.ks_corrected_pvalue),
        ct.significance_string,
        _fmt(ct.num_genes_significant),
        _get_umap(ct, 0),
        _get_umap(ct, 1),
        _get_umap(ct, 2),
        _get_umap(ct, 3),
        ct.get_gene_ranking_string(),
    ]


def _sort_cell_types(cell_types: list[CellTypeClassification]) -> list[CellTypeClassification]:
    """Sort cell types: positive NES first (by FDR asc), then negative NES."""
    positive = []
    negative = []

    for ct in cell_types:
        nes = ct.get_normalized_enrichment_score()
        if math.isnan(nes) or nes < 0:
            negative.append(ct)
        else:
            positive.append(ct)

    positive.sort(key=lambda c: (c.enrichment_fdr, c.enrichment_significance))
    negative.sort(
        key=lambda c: (
            c.enrichment_fdr if not math.isnan(c.enrichment_fdr) else float("inf"),
            c.enrichment_significance
            if not math.isnan(c.enrichment_significance)
            else float("inf"),
        )
    )

    return positive + negative


def write_enrichment_results(
    cell_types: list[CellTypeClassification],
    num_total_cells: int,
    num_passing: int,
    output_dir: Path,
    prefix: str,
    scoring_schema: ScoringSchema,
) -> Path:
    """Write the main 27-column cell type enrichment TSV file.

    Returns the path to the written file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    filename = f"{prefix}_{scoring_schema.method.score_name}_cell_types.txt"
    filepath = output_dir / filename

    sorted_types = _sort_cell_types(cell_types)

    with open(filepath, "w") as f:
        # Header
        f.write("\t".join(COLUMNS) + "\n")

        # Data rows
        for ct in sorted_types:
            values = _row_values(ct, num_total_cells, num_passing)
            f.write("\t".join(values) + "\n")

    log.info("Wrote enrichment results to %s (%d cell types)", filepath, len(sorted_types))
    return filepath


def write_parameters_file(
    output_dir: Path,
    prefix: str,
    scoring_schema: ScoringSchema,
    num_input_genes: int,
    num_mapped_genes: int,
    num_total_cells: int,
    num_passing_cells: int,
    num_cell_types: int,
    num_permutations: int,
    datasets: list[str] | None,
) -> Path:
    """Write a parameters summary file."""
    output_dir.mkdir(parents=True, exist_ok=True)
    filename = f"{prefix}_{scoring_schema.method.score_name}_parameters.txt"
    filepath = output_dir / filename

    with open(filepath, "w") as f:
        f.write(f"Output prefix: {prefix}\n")
        f.write(f"Scoring method: {scoring_schema.method.name}\n")
        f.write(f"Score threshold: {scoring_schema.threshold}\n")
        f.write(f"Min genes per cell: {scoring_schema.min_genes_cells}\n")
        f.write(f"Num permutations: {num_permutations}\n")
        f.write(f"Datasets: {', '.join(datasets) if datasets else 'all'}\n")
        f.write(f"Input genes: {num_input_genes}\n")
        f.write(f"Mapped genes: {num_mapped_genes}\n")
        f.write(f"Total single cells: {num_total_cells}\n")
        f.write(f"Cells passing threshold: {num_passing_cells}\n")
        f.write(f"Cell types found: {num_cell_types}\n")

    log.info("Wrote parameters to %s", filepath)
    return filepath


def write_genes_file(
    cell_types: list[CellTypeClassification],
    output_dir: Path,
    prefix: str,
    scoring_schema: ScoringSchema,
) -> Path:
    """Write gene ranking per cell type."""
    output_dir.mkdir(parents=True, exist_ok=True)
    filename = f"{prefix}_{scoring_schema.method.score_name}_genes.txt"
    filepath = output_dir / filename

    with open(filepath, "w") as f:
        f.write("cell_type\tgene\toccurrence\n")
        for ct in cell_types:
            for go in ct.gene_occurrences:
                f.write(f"{ct.name}\t{go.gene}\t{go.occurrence}\n")

    log.info("Wrote genes file to %s", filepath)
    return filepath


def create_zip(output_dir: Path, zip_path: Path) -> Path:
    """Create a ZIP archive of the output directory."""
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for file in output_dir.rglob("*"):
            if file.is_file() and file != zip_path:
                zf.write(file, file.relative_to(output_dir.parent))

    log.info("Created ZIP archive at %s", zip_path)
    return zip_path
