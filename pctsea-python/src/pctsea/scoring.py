"""Scoring methods for comparing single cells against input protein expressions."""

from __future__ import annotations

import logging
import math
import random

import numpy as np
from scipy.stats import linregress, pearsonr

from pctsea.config import ScoringMethod
from pctsea.models import SingleCell

log = logging.getLogger(__name__)


def filter_by_min_genes(
    cells: list[SingleCell],
    gene_id_to_name: dict[int, str],
    experiment_expressions: dict[int, float],
    min_genes: int,
) -> list[SingleCell]:
    """Keep only cells that express at least min_genes from the input gene list."""
    gene_ids = list(experiment_expressions.keys())
    if min_genes > len(gene_ids):
        msg = f"min_genes_cells ({min_genes}) exceeds the number of input genes ({len(gene_ids)})"
        raise ValueError(msg)

    filtered = []
    for cell in cells:
        count = 0
        for gid in gene_ids:
            cell_expr = cell.get_gene_expression(gid)
            input_expr = experiment_expressions.get(gid, 0.0)
            if cell_expr > 0.0 and input_expr > 0.0:
                count += 1
        if count >= min_genes:
            filtered.append(cell)

    log.info(
        "Filtered cells by min_genes=%d: %d -> %d",
        min_genes,
        len(cells),
        len(filtered),
    )
    return filtered


def score_single_cells(
    cells: list[SingleCell],
    gene_id_to_name: dict[int, str],
    experiment_expressions: dict[int, float],
    method: ScoringMethod,
    threshold: float | None = None,
    min_corr: float | None = None,
) -> int:
    """Score each cell and return the count passing the threshold.

    Modifies cells in place (sets cell.score and cell.genes_used_for_score).
    """
    gene_ids = list(experiment_expressions.keys())
    count_passing = 0

    for cell in cells:
        if method == ScoringMethod.PEARSONS_CORRELATION:
            _pearson_score(cell, gene_ids, gene_id_to_name, experiment_expressions, min_corr)
        elif method == ScoringMethod.SIMPLE_SCORE:
            _simple_score(cell, gene_ids, gene_id_to_name, experiment_expressions, min_corr)
        elif method == ScoringMethod.DOT_PRODUCT:
            _dot_product_score(cell, gene_ids, gene_id_to_name, experiment_expressions)
        elif method == ScoringMethod.REGRESSION:
            _regression_score(cell, gene_ids, gene_id_to_name, experiment_expressions)
        else:
            msg = f"Unsupported scoring method: {method}"
            raise ValueError(msg)

        if not math.isnan(cell.score) and (threshold is None or cell.score >= threshold):
            count_passing += 1

    log.info(
        "Scored %d cells with %s, %d passing threshold",
        len(cells),
        method.score_name,
        count_passing,
    )
    return count_passing


def _add_zero_variance_perturbation(values: list[float]) -> list[float]:
    """Add tiny random perturbation if variance is zero."""
    if len(values) < 2:
        return values
    variance = float(np.var(values))
    if variance == 0.0:
        return [v + random.random() / 1_000_000 for v in values]
    return values


def _collect_nonzero_pairs(
    cell: SingleCell,
    gene_ids: list[int],
    gene_id_to_name: dict[int, str],
    experiment_expressions: dict[int, float],
) -> tuple[list[float], list[float], list[str], int]:
    """Collect expression pairs where both cell and input are > 0.

    Returns (cell_exprs, input_exprs, gene_names, num_nonzero_in_cell).
    """
    cell_exprs: list[float] = []
    input_exprs: list[float] = []
    gene_names: list[str] = []
    num_nonzero = 0

    for gid in gene_ids:
        cell_expr = cell.get_gene_expression(gid)
        input_expr = experiment_expressions.get(gid, 0.0)
        if cell_expr > 0.0 and input_expr > 0.0:
            cell_exprs.append(cell_expr)
            input_exprs.append(input_expr)
            gene_names.append(gene_id_to_name.get(gid, str(gid)))
        if cell_expr != 0.0:
            num_nonzero += 1

    return cell_exprs, input_exprs, gene_names, num_nonzero


def _pearson_score(
    cell: SingleCell,
    gene_ids: list[int],
    gene_id_to_name: dict[int, str],
    experiment_expressions: dict[int, float],
    min_corr: float | None,
) -> None:
    """Pearson correlation between cell and input expressions."""
    cell_exprs, input_exprs, gene_names, _ = _collect_nonzero_pairs(
        cell, gene_ids, gene_id_to_name, experiment_expressions
    )
    cell.genes_used_for_score = gene_names

    if len(cell_exprs) < 2:
        cell.score = float("nan")
        return

    cell_exprs = _add_zero_variance_perturbation(cell_exprs)
    input_exprs = _add_zero_variance_perturbation(input_exprs)

    corr, _ = pearsonr(input_exprs, cell_exprs)
    score = float(corr)

    if min_corr is not None and score < min_corr:
        cell.score = float("nan")
    else:
        cell.score = score


def _simple_score(
    cell: SingleCell,
    gene_ids: list[int],
    gene_id_to_name: dict[int, str],
    experiment_expressions: dict[int, float],
    min_corr: float | None,
) -> None:
    """Simple score = Pearson correlation + count of non-zero pairs."""
    cell_exprs, input_exprs, gene_names, num_nonzero = _collect_nonzero_pairs(
        cell, gene_ids, gene_id_to_name, experiment_expressions
    )
    cell.genes_used_for_score = gene_names

    if len(cell_exprs) < 2:
        cell.score = float("nan")
        return

    cell_exprs = _add_zero_variance_perturbation(cell_exprs)
    input_exprs = _add_zero_variance_perturbation(input_exprs)

    corr, _ = pearsonr(input_exprs, cell_exprs)
    score = float(corr)

    if min_corr is not None and score < min_corr:
        cell.score = float("nan")
    else:
        cell.score = score + num_nonzero


def _dot_product_score(
    cell: SingleCell,
    gene_ids: list[int],
    gene_id_to_name: dict[int, str],
    experiment_expressions: dict[int, float],
) -> None:
    """Normalized dot product of expression vectors."""
    cell_exprs_all: list[float] = []
    input_exprs_all: list[float] = []
    gene_names: list[str] = []

    for gid in gene_ids:
        cell_expr = cell.get_gene_expression(gid)
        input_expr = experiment_expressions.get(gid, 0.0)
        # Include pairs where both > 0
        if cell_expr > 0.0 and input_expr > 0.0:
            cell_exprs_all.append(cell_expr)
            input_exprs_all.append(input_expr)
        if cell_expr > 0.0:
            gene_names.append(gene_id_to_name.get(gid, str(gid)))

    cell.genes_used_for_score = gene_names

    if not cell_exprs_all:
        cell.score = 0.0
        return

    # Normalize each vector and compute dot product
    cell_arr = np.array(cell_exprs_all, dtype=np.float64)
    input_arr = np.array(input_exprs_all, dtype=np.float64)

    cell_norm = np.linalg.norm(cell_arr)
    input_norm = np.linalg.norm(input_arr)

    if cell_norm == 0.0 or input_norm == 0.0:
        cell.score = 0.0
        return

    cell_arr = cell_arr / cell_norm
    input_arr = input_arr / input_norm

    cell.score = float(np.dot(input_arr, cell_arr))


def _regression_score(
    cell: SingleCell,
    gene_ids: list[int],
    gene_id_to_name: dict[int, str],
    experiment_expressions: dict[int, float],
) -> None:
    """R-squared from linear regression."""
    cell_exprs, input_exprs, gene_names, _ = _collect_nonzero_pairs(
        cell, gene_ids, gene_id_to_name, experiment_expressions
    )
    cell.genes_used_for_score = gene_names

    if len(cell_exprs) < 2:
        cell.score = float("nan")
        return

    cell_exprs = _add_zero_variance_perturbation(cell_exprs)
    input_exprs = _add_zero_variance_perturbation(input_exprs)

    result = linregress(input_exprs, cell_exprs)
    r_squared = result.rvalue**2
    cell.score = float(r_squared)
