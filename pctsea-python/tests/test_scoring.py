"""Tests for scoring methods."""

from __future__ import annotations

import math

from pctsea.config import ScoringMethod
from pctsea.models import SingleCell
from pctsea.scoring import filter_by_min_genes, score_single_cells


def test_pearson_score_basic(
    sample_cells: list[SingleCell],
    experiment_expressions: dict[int, float],
    gene_id_to_name: dict[int, str],
) -> None:
    """Test that Pearson scoring assigns valid scores."""
    score_single_cells(
        sample_cells,
        gene_id_to_name,
        experiment_expressions,
        ScoringMethod.PEARSONS_CORRELATION,
        threshold=None,
        min_corr=None,
    )
    # All cells with at least 2 nonzero pairs should get a score
    scored = [c for c in sample_cells if not math.isnan(c.score)]
    assert len(scored) > 0
    # Scores should be between -1 and 1
    for c in scored:
        assert -1.0 <= c.score <= 1.0


def test_simple_score_adds_count(
    sample_cells: list[SingleCell],
    experiment_expressions: dict[int, float],
    gene_id_to_name: dict[int, str],
) -> None:
    """Simple score should be larger than raw correlation (correlation + count)."""
    score_single_cells(
        sample_cells,
        gene_id_to_name,
        experiment_expressions,
        ScoringMethod.SIMPLE_SCORE,
        threshold=None,
        min_corr=None,
    )
    scored = [c for c in sample_cells if not math.isnan(c.score)]
    # Simple score = correlation + nonzero count, so it should be > 1 for cells with multiple genes
    assert any(c.score > 1.0 for c in scored)


def test_dot_product_score(
    sample_cells: list[SingleCell],
    experiment_expressions: dict[int, float],
    gene_id_to_name: dict[int, str],
) -> None:
    """Dot product scores should be between 0 and 1 for normalized vectors."""
    score_single_cells(
        sample_cells,
        gene_id_to_name,
        experiment_expressions,
        ScoringMethod.DOT_PRODUCT,
    )
    scored = [c for c in sample_cells if not math.isnan(c.score)]
    assert len(scored) > 0
    for c in scored:
        assert 0.0 <= c.score <= 1.0 + 1e-10


def test_regression_score(
    sample_cells: list[SingleCell],
    experiment_expressions: dict[int, float],
    gene_id_to_name: dict[int, str],
) -> None:
    """R-squared should be between 0 and 1."""
    score_single_cells(
        sample_cells,
        gene_id_to_name,
        experiment_expressions,
        ScoringMethod.REGRESSION,
    )
    scored = [c for c in sample_cells if not math.isnan(c.score)]
    assert len(scored) > 0
    for c in scored:
        assert 0.0 <= c.score <= 1.0 + 1e-10


def test_filter_by_min_genes(
    sample_cells: list[SingleCell],
    experiment_expressions: dict[int, float],
    gene_id_to_name: dict[int, str],
) -> None:
    """Filtering should reduce cell count."""
    filtered = filter_by_min_genes(sample_cells, gene_id_to_name, experiment_expressions, 3)
    # Only cells with gene 2 (every other cell) will have all 3 genes
    assert len(filtered) < len(sample_cells)
    assert len(filtered) > 0
