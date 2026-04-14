"""Tests for output generation."""

from __future__ import annotations

from pathlib import Path

from pctsea.config import ScoringMethod
from pctsea.models import CellTypeClassification, ScoringSchema
from pctsea.output import COLUMNS, write_enrichment_results


def test_column_count() -> None:
    """Verify we have exactly 27 columns."""
    assert len(COLUMNS) == 27


def test_write_enrichment_results_creates_file(tmp_path: Path) -> None:
    """Test that the enrichment file is written correctly."""
    ct = CellTypeClassification(name="test_type", cell_type_id=1)
    ct.hypergeometric_pvalue = 0.01
    ct.enrichment_score = 0.5
    ct.enrichment_significance = 0.03
    ct.enrichment_fdr = 0.04
    ct.random_enrichment_scores = [0.1, 0.2, 0.3, 0.4, 0.5]

    schema = ScoringSchema(
        method=ScoringMethod.PEARSONS_CORRELATION,
        threshold=0.2,
        min_genes_cells=100,
    )

    path = write_enrichment_results(
        [ct],
        num_total_cells=1000,
        num_passing=500,
        output_dir=tmp_path,
        prefix="test",
        scoring_schema=schema,
    )

    assert path.exists()
    lines = path.read_text().strip().split("\n")
    assert len(lines) == 2  # header + 1 data row
    header = lines[0].split("\t")
    assert len(header) == 27
    assert header[0] == "cell_type"
