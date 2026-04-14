"""Tests for KS-style enrichment score calculation."""

from __future__ import annotations

import math

from pctsea.enrichment import calculate_enrichment_scores
from pctsea.models import CellTypeClassification, SingleCell


def test_enrichment_basic() -> None:
    """Test enrichment with a small synthetic ranked list."""
    # Create ranked cells: type_a scores high, type_b scores low
    cells: list[SingleCell] = []
    for i in range(20):
        cell = SingleCell(id=i, name=f"cell_{i}")
        cell.score = 20.0 - i  # Descending scores
        if i < 5:
            cell.cell_type = "type_a"
            cell.cell_type_id = 1
        else:
            cell.cell_type = "type_b"
            cell.cell_type_id = 2
        cells.append(cell)

    ct_a = CellTypeClassification(name="type_a", cell_type_id=1)
    ct_b = CellTypeClassification(name="type_b", cell_type_id=2)

    calculate_enrichment_scores([ct_a, ct_b], cells, is_permutation=False)

    # type_a should have positive enrichment (concentrated at top)
    assert not math.isnan(ct_a.enrichment_score)
    assert ct_a.enrichment_score > 0
    assert ct_a.supremum_x > 0

    # type_b should have negative enrichment
    assert not math.isnan(ct_b.enrichment_score)
    assert ct_b.enrichment_score < 0


def test_enrichment_permutation() -> None:
    """Test that permutation mode stores random scores."""
    cells: list[SingleCell] = []
    for i in range(10):
        cell = SingleCell(id=i, name=f"cell_{i}")
        cell.score = 10.0 - i
        cell.cell_type = "type_a" if i < 3 else "type_b"
        cell.cell_type_id = 1 if i < 3 else 2
        cells.append(cell)

    ct = CellTypeClassification(name="type_a", cell_type_id=1)
    calculate_enrichment_scores([ct], cells, is_permutation=True)

    assert len(ct.random_enrichment_scores) == 1
    assert len(ct.random_ks_statistics) == 1
