"""CLI entry point for PCTSEA."""

from __future__ import annotations

import logging
import math
import sys
from datetime import datetime
from pathlib import Path

import click

from pctsea.config import (
    MIN_CELLS_PASSING_THRESHOLD,
    MIN_NUM_MAPPED_GENES,
    CellTypeBranch,
    InputDataType,
    ScoringMethod,
)
from pctsea.models import ScoringSchema, SingleCell, SingleCellSet

log = logging.getLogger(__name__)


@click.command()
@click.option(
    "--eef", required=True, type=click.Path(exists=True), help="Experimental expression file"
)
@click.option("--out", required=True, help="Output prefix")
@click.option("--scoring-method", required=True, help="Comma-separated scoring methods")
@click.option("--min-score", required=True, help="Comma-separated min score thresholds")
@click.option("--min-genes-cells", required=True, help="Comma-separated min genes per cell")
@click.option(
    "--input-data-type",
    required=True,
    type=click.Choice([e.name for e in InputDataType], case_sensitive=False),
    help="Type of input proteomics data",
)
@click.option("--perm", default=10, type=int, help="Number of permutations")
@click.option("--datasets", default=None, help="Comma-separated dataset tags")
@click.option(
    "--cell-types-classification",
    default="ORIGINAL",
    type=click.Choice([e.name for e in CellTypeBranch], case_sensitive=False),
    help="Cell type classification level",
)
@click.option(
    "--plot-negative-enriched",
    is_flag=True,
    default=False,
    help="Include negatively enriched cell types",
)
@click.option("--write-scores", is_flag=True, default=False, help="Write per-cell score files")
@click.option("--min-corr", default=-1.0, type=float, help="Minimum correlation threshold")
@click.option("--create-zip/--no-create-zip", default=True, help="Create ZIP archive of results")
@click.option("--db-host", default="localhost", help="MongoDB host")
@click.option("--db-port", default=27017, type=int, help="MongoDB port")
@click.option("-v", "--verbose", is_flag=True, default=False, help="Verbose logging")
def main(
    eef: str,
    out: str,
    scoring_method: str,
    min_score: str,
    min_genes_cells: str,
    input_data_type: str,
    perm: int,
    datasets: str | None,
    cell_types_classification: str,
    plot_negative_enriched: bool,
    write_scores: bool,
    min_corr: float,
    create_zip: bool,
    db_host: str,
    db_port: int,
    verbose: bool,
) -> None:
    """PCTSEA - Proteomics Cell Type Single-cell Enrichment Analysis."""
    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
    )

    # Parse scoring schemas
    methods = [s.strip() for s in scoring_method.split(",")]
    scores = [s.strip() for s in min_score.split(",")]
    genes_cells = [s.strip() for s in min_genes_cells.split(",")]

    if len(methods) != len(scores) or len(methods) != len(genes_cells):
        msg = "Error: --scoring-method, --min-score, and --min-genes-cells "
        msg += "must have equal number of values"
        click.echo(msg, err=True)
        sys.exit(1)

    scoring_schemas: list[ScoringSchema] = []
    for m, s, g in zip(methods, scores, genes_cells, strict=True):
        method = ScoringMethod.from_name(m)
        if method is None:
            click.echo(
                f"Error: Unknown scoring method '{m}'. Options: {ScoringMethod.choices()}",
                err=True,
            )
            sys.exit(1)
        scoring_schemas.append(
            ScoringSchema(method=method, threshold=float(s), min_genes_cells=int(g))
        )

    cell_type_branch = CellTypeBranch.from_name(cell_types_classification)
    if cell_type_branch is None:
        click.echo(
            f"Error: Unknown cell type classification '{cell_types_classification}'", err=True
        )
        sys.exit(1)

    dataset_tags = [d.strip() for d in datasets.split(",")] if datasets else None
    min_corr_val = min_corr if min_corr > -1.0 else None

    # Run the pipeline
    _run_pipeline(
        eef_path=Path(eef),
        prefix=out,
        scoring_schemas=scoring_schemas,
        cell_type_branch=cell_type_branch,
        num_permutations=perm,
        dataset_tags=dataset_tags,
        min_corr=min_corr_val,
        create_zip_file=create_zip,
        write_scores_file=write_scores,
        db_host=db_host,
        db_port=db_port,
    )


def _run_pipeline(
    eef_path: Path,
    prefix: str,
    scoring_schemas: list[ScoringSchema],
    cell_type_branch: CellTypeBranch,
    num_permutations: int,
    dataset_tags: list[str] | None,
    min_corr: float | None,
    create_zip_file: bool,
    write_scores_file: bool,
    db_host: str,
    db_port: int,
) -> None:
    """Execute the PCTSEA analysis pipeline."""
    from pctsea.cell_types import resolve_cell_type
    from pctsea.db import PctseaDB
    from pctsea.enrichment import calculate_enrichment_scores
    from pctsea.hypergeometric import calculate_hypergeometric_stats
    from pctsea.input_parser import read_expression_file
    from pctsea.output import (
        create_zip,
        write_enrichment_results,
        write_genes_file,
        write_parameters_file,
    )
    from pctsea.permutation import calculate_significance_by_permutations
    from pctsea.scoring import filter_by_min_genes, score_single_cells
    from pctsea.umap_clustering import build_gene_occurrences, cluster_cell_types

    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    base_output_dir = eef_path.parent / f"{timestamp}_{prefix}"

    # Step 1: Read input expression file
    log.info("Reading experimental expression file: %s", eef_path)
    experiment_expressions_by_name = read_expression_file(eef_path)
    num_input_genes = len(experiment_expressions_by_name)
    log.info("Read %d genes from input file", num_input_genes)

    # Step 2: Connect to MongoDB
    log.info("Connecting to MongoDB at %s:%d", db_host, db_port)
    db = PctseaDB(host=db_host, port=db_port)

    # Validate datasets
    if dataset_tags:
        for tag in dataset_tags:
            ds = db.find_dataset_by_tag(tag)
            if ds is None:
                log.error("Dataset '%s' not found in database", tag)
                sys.exit(1)

    # Step 3: Load single cell metadata
    log.info("Loading single cell metadata...")
    single_cell_set = SingleCellSet()
    for cell_id_counter, doc in enumerate(db.get_single_cells(dataset_tags)):
        cell = SingleCell(
            id=cell_id_counter,
            name=doc.get("name"),
            dataset_tag=doc.get("datasetTag"),
            biomaterial=doc.get("biomaterial"),
        )
        raw_type = doc.get("type", "")
        if raw_type:
            resolved_name, ct_id = resolve_cell_type(raw_type, cell_type_branch)
            cell.cell_type = resolved_name
            cell.original_cell_type = raw_type
            cell.cell_type_id = ct_id
        single_cell_set.add_single_cell(cell)

    log.info("Loaded %d single cells", single_cell_set.num_cells)

    # Step 4: Load gene expressions from DB
    log.info("Loading gene expressions from database...")
    gene_name_to_id: dict[str, int] = {}
    gene_id_to_name: dict[int, str] = {}
    experiment_expressions: dict[int, float] = {}  # gene_id -> expression

    input_gene_names = list(experiment_expressions_by_name.keys())

    for gene_id, gene_name in enumerate(input_gene_names):
        gene_name_to_id[gene_name] = gene_id
        gene_id_to_name[gene_id] = gene_name
        experiment_expressions[gene_id] = experiment_expressions_by_name[gene_name]

    # Query expressions for all input genes
    mapped_genes = 0
    for doc in db.get_expressions_by_genes(input_gene_names, dataset_tags):
        gene_name = doc.get("gene", "").upper()
        cell_name = doc.get("cellName")
        expression = doc.get("expression", 0.0)

        if gene_name not in gene_name_to_id:
            continue
        if cell_name is None:
            continue

        cell_id = single_cell_set.get_cell_id_by_name(cell_name)
        if cell_id == -1:
            continue

        gene_id = gene_name_to_id[gene_name]
        cell = single_cell_set.get_cell_by_id(cell_id)
        if cell:
            cell.add_gene_expression(gene_id, float(expression))
            mapped_genes += 1

    log.info("Mapped %d gene-cell expression pairs", mapped_genes)

    # Remove cells with no overlapping genes
    all_cells = [
        c
        for c in single_cell_set.cell_list
        if any(c.get_gene_expression(gid) > 0 for gid in experiment_expressions)
    ]
    log.info("Cells with at least one input gene expressed: %d", len(all_cells))

    num_mapped = sum(
        1
        for gid in experiment_expressions
        if any(c.get_gene_expression(gid) > 0 for c in all_cells)
    )
    log.info("Input genes found in database: %d / %d", num_mapped, num_input_genes)

    if num_mapped < MIN_NUM_MAPPED_GENES:
        log.error(
            "Only %d genes mapped (minimum %d required). Aborting.",
            num_mapped,
            MIN_NUM_MAPPED_GENES,
        )
        sys.exit(1)

    # Step 5: For each scoring schema, run the analysis
    for schema in scoring_schemas:
        log.info(
            "=== Running scoring method: %s (threshold=%.3f, min_genes=%d) ===",
            schema.method.name,
            schema.threshold,
            schema.min_genes_cells,
        )

        output_dir = base_output_dir / schema.method.score_name
        output_dir.mkdir(parents=True, exist_ok=True)
        cell_types_dir = output_dir / "cell_types_charts"

        # 5a: Filter cells by min genes
        filtered_cells = filter_by_min_genes(
            all_cells, gene_id_to_name, experiment_expressions, schema.min_genes_cells
        )

        if len(filtered_cells) < MIN_CELLS_PASSING_THRESHOLD:
            log.warning(
                "Only %d cells pass min_genes filter (need %d). Skipping.",
                len(filtered_cells),
                MIN_CELLS_PASSING_THRESHOLD,
            )
            continue

        # 5b: Score cells
        num_passing = score_single_cells(
            filtered_cells,
            gene_id_to_name,
            experiment_expressions,
            schema.method,
            threshold=schema.threshold,
            min_corr=min_corr,
        )
        log.info("Cells passing score threshold: %d / %d", num_passing, len(filtered_cells))

        if num_passing < MIN_CELLS_PASSING_THRESHOLD:
            log.warning(
                "Only %d cells pass score threshold (need %d). Skipping.",
                num_passing,
                MIN_CELLS_PASSING_THRESHOLD,
            )
            continue

        # Sort cells by score descending (NaN goes to end)
        ranked_cells = sorted(
            filtered_cells,
            key=lambda c: c.score if not math.isnan(c.score) else float("-inf"),
            reverse=True,
        )

        # Keep only cells passing threshold for enrichment
        passing_cells = [
            c for c in ranked_cells if not math.isnan(c.score) and c.score >= schema.threshold
        ]

        # 5c: Hypergeometric test
        cell_type_classifications = calculate_hypergeometric_stats(
            passing_cells, schema.threshold, num_passing
        )
        log.info("Cell types identified: %d", len(cell_type_classifications))

        # 5d: Enrichment scores
        calculate_enrichment_scores(
            cell_type_classifications,
            passing_cells,
            is_permutation=False,
            output_dir=cell_types_dir,
        )

        # 5e: Permutation testing
        if num_permutations > 0:
            calculate_significance_by_permutations(
                cell_type_classifications,
                passing_cells,
                num_permutations,
            )

        # 5f: Build gene occurrences
        build_gene_occurrences(cell_type_classifications, schema.threshold)

        # 5g: UMAP clustering
        cluster_cell_types(cell_type_classifications)

        # 5h: Write output files
        write_enrichment_results(
            cell_type_classifications,
            num_total_cells=len(all_cells),
            num_passing=num_passing,
            output_dir=output_dir,
            prefix=prefix,
            scoring_schema=schema,
        )

        write_parameters_file(
            output_dir=output_dir,
            prefix=prefix,
            scoring_schema=schema,
            num_input_genes=num_input_genes,
            num_mapped_genes=num_mapped,
            num_total_cells=len(all_cells),
            num_passing_cells=num_passing,
            num_cell_types=len(cell_type_classifications),
            num_permutations=num_permutations,
            datasets=dataset_tags,
        )

        write_genes_file(cell_type_classifications, output_dir, prefix, schema)

        log.info("Results written to %s", output_dir)

    # ZIP packaging
    if create_zip_file:
        zip_path = base_output_dir.with_suffix(".zip")
        create_zip(base_output_dir, zip_path)

    db.close()
    log.info("PCTSEA analysis complete.")


if __name__ == "__main__":
    main()
