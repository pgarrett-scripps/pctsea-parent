"""Permutation testing for enrichment significance and FDR calculation."""

from __future__ import annotations

import logging
import math
import random
from bisect import bisect_left

from pctsea.enrichment import calculate_enrichment_scores
from pctsea.models import CellTypeClassification, SingleCell

log = logging.getLogger(__name__)


def calculate_significance_by_permutations(
    cell_types: list[CellTypeClassification],
    ranked_cells: list[SingleCell],
    num_permutations: int,
) -> None:
    """Estimate significance via permutation of cell type labels.

    Shuffles cell type labels, recalculates enrichment scores, then computes
    empirical p-values, NES, and FDR. Modifies cell_types in place.
    """
    # Save original cell types
    original_types = [(c.cell_type, c.cell_type_id) for c in ranked_cells]

    log.info("Running %d permutations...", num_permutations)
    for iteration in range(1, num_permutations + 1):
        if iteration % max(1, num_permutations // 10) == 0:
            log.info("Permutation %d/%d", iteration, num_permutations)

        # Shuffle cell type labels
        shuffled = [(ct, cid) for ct, cid in original_types]
        random.shuffle(shuffled)
        for i, cell in enumerate(ranked_cells):
            cell.cell_type = shuffled[i][0]
            cell.cell_type_id = shuffled[i][1]

        # Recalculate enrichment with shuffled labels
        calculate_enrichment_scores(cell_types, ranked_cells, is_permutation=True, parallel=False)

    # Restore original labels
    for i, cell in enumerate(ranked_cells):
        cell.cell_type = original_types[i][0]
        cell.cell_type_id = original_types[i][1]

    log.info("Permutations complete. Calculating significance...")

    # Calculate empirical p-values
    _calculate_empirical_pvalues(cell_types)

    # Calculate FDR
    _calculate_fdr(cell_types)

    # BH correction of KS p-values
    _bh_correct_ks_pvalues(cell_types)


def _calculate_empirical_pvalues(cell_types: list[CellTypeClassification]) -> None:
    """Compute empirical p-value for each cell type.

    p-value = (count(random >= real) + 1) / (num_permutations + 1) for positive scores
    p-value = (count(random <= real) + 1) / (num_permutations + 1) for negative scores
    """
    for ct in cell_types:
        real_score = ct.enrichment_score
        if math.isnan(real_score):
            ct.enrichment_significance = float("nan")
            continue

        valid_randoms = [s for s in ct.random_enrichment_scores if not math.isnan(s)]
        if not valid_randoms:
            ct.enrichment_significance = float("nan")
            continue

        if real_score >= 0.0:
            count = sum(1 for s in valid_randoms if s >= real_score)
        else:
            count = sum(1 for s in valid_randoms if s <= real_score)

        pvalue = (count + 1) / (len(valid_randoms) + 1)
        ct.enrichment_significance = pvalue


def _calculate_fdr(cell_types: list[CellTypeClassification]) -> None:
    """Calculate False Discovery Rate using the GSEA-style approach.

    FDR = (snull/sobs) * (nobs/nnull), where:
    - snull = # random NES scores >= real NES
    - sobs = # real NES scores >= real NES
    - nobs = total # real NES scores
    - nnull = total # random NES scores
    """
    # Collect all positive normalized scores
    real_scores: list[float] = []
    random_scores: list[float] = []

    for ct in cell_types:
        nes = ct.get_normalized_enrichment_score()
        if math.isnan(nes) or nes < 0.0:
            continue
        real_scores.append(nes)
        for rnes in ct.normalized_random_enrichment_scores:
            if rnes >= 0.0:
                random_scores.append(rnes)

    if not real_scores or not random_scores:
        for ct in cell_types:
            ct.enrichment_fdr = float("nan")
        return

    real_sorted = sorted(real_scores)
    random_sorted = sorted(random_scores)
    nobs = len(real_sorted)
    nnull = len(random_sorted)

    for ct in cell_types:
        nes = ct.get_normalized_enrichment_score()
        if math.isnan(nes) or nes < 0.0:
            ct.enrichment_fdr = float("nan")
            continue

        # Count real scores >= nes
        idx_real = bisect_left(real_sorted, nes)
        sobs = nobs - idx_real

        # Count random scores >= nes
        idx_random = bisect_left(random_sorted, nes)
        snull = nnull - idx_random

        if snull == 0:
            ct.enrichment_fdr = 0.0
        elif sobs == 0:
            ct.enrichment_fdr = float("nan")
        else:
            ct.enrichment_fdr = (snull / sobs) * (nobs / nnull)


def _bh_correct_ks_pvalues(cell_types: list[CellTypeClassification]) -> None:
    """Apply Benjamini-Hochberg correction to KS test p-values."""
    pvalues = [ct.ks_pvalue for ct in cell_types]
    corrected = benjamini_hochberg(pvalues)

    for i, ct in enumerate(cell_types):
        ct.ks_corrected_pvalue = corrected[i]
        if ct.ks_corrected_pvalue < 0.001:
            ct.significance_string = "***"
        elif ct.ks_corrected_pvalue < 0.01:
            ct.significance_string = "**"
        elif ct.ks_corrected_pvalue < 0.05:
            ct.significance_string = "*"
        else:
            ct.significance_string = ""


def benjamini_hochberg(pvalues: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction.

    Returns adjusted p-values in the same order as input.
    """
    n = len(pvalues)
    if n == 0:
        return []

    # Sort indices by p-value
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])

    adjusted = [0.0] * n
    cummin = 1.0

    for rank_from_end, (orig_idx, pval) in enumerate(reversed(indexed)):
        rank = n - rank_from_end  # 1-based rank from the sorted order
        bh_val = pval * n / rank
        cummin = min(cummin, bh_val)
        adjusted[orig_idx] = min(cummin, 1.0)

    return adjusted
