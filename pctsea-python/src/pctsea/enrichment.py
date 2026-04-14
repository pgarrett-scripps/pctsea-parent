"""KS-style weighted enrichment score calculation."""

from __future__ import annotations

import logging
import math
from pathlib import Path

from scipy.stats import ks_2samp

from pctsea.models import CellTypeClassification, SingleCell

log = logging.getLogger(__name__)


def calculate_enrichment_scores(
    cell_types: list[CellTypeClassification],
    ranked_cells: list[SingleCell],
    is_permutation: bool = False,
    output_dir: Path | None = None,
    parallel: bool = True,
) -> None:
    """Calculate weighted enrichment scores for each cell type.

    Modifies cell_types in place.
    The ranked_cells list must already be sorted by score descending.
    """
    if parallel and len(cell_types) > 1 and not is_permutation:
        _calculate_parallel(cell_types, ranked_cells, is_permutation, output_dir)
    else:
        for ct in cell_types:
            _calculate_single(ct, ranked_cells, is_permutation, output_dir)


def _calculate_parallel(
    cell_types: list[CellTypeClassification],
    ranked_cells: list[SingleCell],
    is_permutation: bool,
    output_dir: Path | None,
) -> None:
    """Process enrichment scores using thread pool."""
    # For simplicity and to avoid pickling issues, use sequential processing
    # with the option to parallelize later via multiprocessing if needed
    for ct in cell_types:
        _calculate_single(ct, ranked_cells, is_permutation, output_dir)


def _calculate_single(
    cell_type: CellTypeClassification,
    ranked_cells: list[SingleCell],
    is_permutation: bool,
    output_dir: Path | None,
) -> None:
    """Calculate weighted enrichment score for a single cell type."""
    ct_id = cell_type.cell_type_id

    # Compute denominators (sum of scores for in-type and out-type)
    denominator_a = 0.0
    denominator_b = 0.0
    cells_of_type_count = 0

    for cell in ranked_cells:
        if cell.cell_type_id == ct_id:
            denominator_a += cell.score
            cells_of_type_count += 1
        else:
            denominator_b += cell.score

    # Skip if too few cells of this type
    if cells_of_type_count <= 1:
        if not is_permutation:
            cell_type.set_enrichment(float("nan"), float("nan"), -1, -1.0)
        else:
            cell_type.add_random_enrichment(float("nan"), float("nan"))
        return

    # Avoid division by zero
    if denominator_a == 0.0 or denominator_b == 0.0:
        if not is_permutation:
            cell_type.set_enrichment(float("nan"), float("nan"), -1, -1.0)
        else:
            cell_type.add_random_enrichment(float("nan"), float("nan"))
        return

    # Walk the ranked list computing weighted cumulative distributions
    supremum = 0.0
    supremum_x = 0
    numerator_a = 0.0
    numerator_b = 0.0
    previous_a = 0.0
    previous_b = 0.0
    cell_type.num_cell_type_scores = 0
    cell_type.num_other_cell_type_scores = 0

    # For secondary enrichment detection (only on real data)
    differences: list[float] = []
    scores_of_type: list[float] = []
    scores_of_other: list[float] = []

    n = len(ranked_cells)

    for i, cell in enumerate(ranked_cells):
        if cell.cell_type_id == ct_id:
            numerator_a += cell.score
            a = numerator_a / denominator_a
            b = previous_b
            cell_type.num_cell_type_scores += 1
            if not is_permutation and not math.isnan(cell.score):
                scores_of_type.append(cell.score)
        else:
            numerator_b += cell.score
            a = previous_a
            b = numerator_b / denominator_b
            cell_type.num_other_cell_type_scores += 1
            if not is_permutation and not math.isnan(cell.score):
                scores_of_other.append(cell.score)

        difference = a - b
        if not is_permutation:
            differences.append(difference)

        if abs(difference) > abs(supremum):
            supremum = difference
            supremum_x = i + 1

        previous_a = a
        previous_b = b

    # Compute KS D-statistic with sample-size correction
    if cell_type.num_cell_type_scores > 1 and cell_type.num_other_cell_type_scores > 1:
        size_a = cell_type.num_cell_type_scores
        size_b = cell_type.num_other_cell_type_scores
        sqrt_factor = math.sqrt((size_a * size_b) / (size_a + size_b))
        d_statistic = supremum * sqrt_factor

        if not is_permutation:
            # Two-sample KS test p-value
            if scores_of_type and scores_of_other:
                _, ks_pvalue = ks_2samp(scores_of_type, scores_of_other)
            else:
                ks_pvalue = 1.0

            cell_type.size_a = size_a
            cell_type.size_b = size_b
            normalized_sup_x = supremum_x / n if n > 0 else 0.0
            cell_type.set_enrichment(supremum, d_statistic, supremum_x, normalized_sup_x)
            cell_type.ks_pvalue = ks_pvalue

            # Look for secondary enrichment
            _find_secondary_enrichment(
                cell_type,
                differences,
                supremum_x,
                ranked_cells,
                ct_id,
                denominator_a,
                denominator_b,
            )

            # Write output files
            if output_dir is not None:
                _write_enrichment_files(cell_type, differences, output_dir)
        else:
            cell_type.add_random_enrichment(supremum, d_statistic)
    else:
        if not is_permutation:
            cell_type.set_enrichment(float("nan"), float("nan"), -1, -1.0)
        else:
            cell_type.add_random_enrichment(float("nan"), float("nan"))


def _find_secondary_enrichment(
    cell_type: CellTypeClassification,
    differences: list[float],
    supremum_x: int,
    ranked_cells: list[SingleCell],
    ct_id: int,
    denominator_a: float,
    denominator_b: float,
) -> None:
    """Detect secondary enrichment peak before the main supremum."""
    # Look backwards from supremum for a transition from negative to positive difference
    secondary_index = 0
    prev_diff = 0.0
    for i in range(supremum_x - 2, -1, -1):
        diff = differences[i]
        if diff > 0 and prev_diff < 0:
            secondary_index = i
            break
        prev_diff = diff

    if secondary_index == 0:
        return

    # Recalculate enrichment up to secondary index using absolute values
    numerator_a = 0.0
    numerator_b = 0.0
    secondary_supremum = 0.0
    secondary_sup_x = 0

    for i in range(secondary_index + 1):
        cell = ranked_cells[i]
        if cell.cell_type_id == ct_id:
            numerator_a += abs(cell.score)
            a = numerator_a / denominator_a
            b_val = numerator_b / denominator_b if denominator_b != 0 else 0.0
        else:
            numerator_b += abs(cell.score)
            b_val = numerator_b / denominator_b
            a = numerator_a / denominator_a if denominator_a != 0 else 0.0

        if a > b_val:
            diff = a - b_val
            if abs(diff) > abs(secondary_supremum):
                secondary_supremum = diff
                secondary_sup_x = i + 1

    if secondary_sup_x != 0:
        cell_type.secondary_enrichment_score = secondary_supremum
        cell_type.secondary_supremum_x = secondary_sup_x


def _write_enrichment_files(
    cell_type: CellTypeClassification,
    differences: list[float],
    output_dir: Path,
) -> None:
    """Write per-cell-type enrichment walk data for plotting."""
    output_dir.mkdir(parents=True, exist_ok=True)
    safe_name = cell_type.name.replace("/", "_").replace("\\", "_")

    # Write differences file
    diff_file = output_dir / f"{safe_name}_ews.txt"
    with open(diff_file, "w") as f:
        f.write("cell_index\tdifference\n")
        for i, d in enumerate(differences):
            f.write(f"{i + 1}\t{d}\n")
