"""Tests for hypergeometric test."""

from __future__ import annotations

import math

from pctsea.config import ScoringMethod
from pctsea.hypergeometric import calculate_hypergeometric_stats
from pctsea.models import SingleCell
from pctsea.scoring import score_single_cells


def test_hypergeometric_basic(
    sample_cells: list[SingleCell],
    experiment_expressions: dict[int, float],
    gene_id_to_name: dict[int, str],
) -> None:
    """Test that hypergeometric test produces valid p-values."""
    # First score the cells
    score_single_cells(
        sample_cells,
        gene_id_to_name,
        experiment_expressions,
        ScoringMethod.PEARSONS_CORRELATION,
        threshold=0.0,
        min_corr=None,
    )
    num_passing = sum(1 for c in sample_cells if not math.isnan(c.score) and c.score >= 0.0)

    results = calculate_hypergeometric_stats(sample_cells, 0.0, num_passing)

    assert len(results) > 0
    for ct in results:
        assert 0.0 <= ct.hypergeometric_pvalue <= 1.0
        assert ct.num_cells_of_type > 0
