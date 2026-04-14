"""Hypergeometric test for cell type enrichment."""

from __future__ import annotations

import logging
import math

from scipy.stats import hypergeom

from pctsea.cell_types import get_cell_type_name
from pctsea.models import CellTypeClassification, SingleCell

log = logging.getLogger(__name__)


def calculate_hypergeometric_stats(
    cells: list[SingleCell],
    threshold: float,
    total_passing: int,
) -> list[CellTypeClassification]:
    """Compute hypergeometric enrichment per cell type.

    Args:
        cells: All single cells (scored, with cell_type_id assigned).
        threshold: Score threshold for "passing".
        total_passing: Number of cells that passed the score threshold.

    Returns:
        List of CellTypeClassification objects with hypergeometric p-values.
    """
    # Group cells by cell_type_id (skip unknown = -1)
    cell_type_ids = sorted({c.cell_type_id for c in cells if c.cell_type_id != -1})
    log.info("%d distinct cell types", len(cell_type_ids))

    n_total = len(cells)  # N: population size
    n_passing = total_passing  # K: successes in population

    results: list[CellTypeClassification] = []
    num_significant = 0

    for ct_id in cell_type_ids:
        cells_of_type = [c for c in cells if c.cell_type_id == ct_id]
        n_of_type = len(cells_of_type)  # n: sample size

        # Count cells of this type passing threshold
        n_of_type_passing = sum(
            1 for c in cells_of_type if not math.isnan(c.score) and c.score >= threshold
        )

        # Hypergeometric test (Salva's approach)
        # P(X >= k) = sf(k-1, N, K, n)
        p_value = float(hypergeom.sf(n_of_type_passing - 1, n_total, n_passing, n_of_type))
        p_value = max(p_value, 0.0)

        if p_value < 0.05:
            num_significant += 1

        cell_type_name = get_cell_type_name(ct_id) or str(ct_id)

        ct_class = CellTypeClassification(name=cell_type_name, cell_type_id=ct_id)
        ct_class.hypergeometric_pvalue = p_value
        ct_class.num_cells_of_type = n_of_type
        ct_class.num_cells_of_type_passing_threshold = n_of_type_passing
        ct_class.single_cells_of_type = cells_of_type

        # Casimir's enrichment score: log2((k/K) / (n/N))
        if n_of_type > 0 and n_passing > 0:
            ratio = (n_of_type_passing / n_passing) / (n_of_type / n_total)
            if ratio > 0:
                ct_class.casimirs_enrichment_score = float(math.log2(ratio))
            else:
                ct_class.casimirs_enrichment_score = float("-inf")
        else:
            ct_class.casimirs_enrichment_score = 0.0

        results.append(ct_class)

    log.info(
        "Hypergeometric statistics calculated. %d cell types significant (p<0.05)",
        num_significant,
    )
    return results
