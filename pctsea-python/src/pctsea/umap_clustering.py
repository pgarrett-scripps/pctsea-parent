"""UMAP dimensionality reduction for cell type clustering."""

from __future__ import annotations

import logging
from collections import Counter

import numpy as np

from pctsea.config import UMAP_DIMENSIONS
from pctsea.models import CellTypeClassification, GeneOccurrence

log = logging.getLogger(__name__)


def cluster_cell_types(
    cell_types: list[CellTypeClassification],
    fdr_threshold: float | None = None,
) -> None:
    """Run UMAP clustering on cell types based on gene occurrence profiles.

    Modifies cell_types in place, setting umap_components.
    """
    # Filter to cell types with positive enrichment scores
    eligible = []
    for ct in cell_types:
        if fdr_threshold is not None and ct.enrichment_fdr > fdr_threshold:
            continue
        if ct.enrichment_score > 0 and ct.gene_occurrences:
            eligible.append(ct)

    if len(eligible) < 3:
        log.warning("Too few cell types (%d) for UMAP clustering, skipping", len(eligible))
        return

    # Build gene set from all eligible cell types
    gene_set: set[str] = set()
    for ct in eligible:
        for go in ct.gene_occurrences:
            gene_set.add(go.gene)

    gene_list = sorted(gene_set)
    gene_index = {g: i for i, g in enumerate(gene_list)}

    # Build feature matrix (cell types x genes)
    matrix = np.zeros((len(eligible), len(gene_list)), dtype=np.float32)
    for i, ct in enumerate(eligible):
        for go in ct.gene_occurrences:
            if go.gene in gene_index:
                matrix[i, gene_index[go.gene]] = float(go.occurrence)

    # Configure UMAP
    n_neighbors = max(3, min(int(len(eligible) * 0.05), 50))
    n_components = min(UMAP_DIMENSIONS, len(eligible) - 1)

    try:
        import umap

        reducer = umap.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=0.1,
            metric="euclidean",
        )
        embedding = reducer.fit_transform(matrix)

        for i, ct in enumerate(eligible):
            components = list(float(x) for x in embedding[i])
            # Pad to UMAP_DIMENSIONS if needed
            while len(components) < UMAP_DIMENSIONS:
                components.append(float("nan"))
            ct.umap_components = components

        log.info("UMAP clustering complete for %d cell types", len(eligible))
    except Exception:
        log.exception("UMAP clustering failed")


def build_gene_occurrences(
    cell_types: list[CellTypeClassification],
    threshold: float,
) -> None:
    """Build gene occurrence rankings for each cell type.

    Counts how many cells of each type used each gene for scoring.
    """
    for ct in cell_types:
        gene_counts: Counter[str] = Counter()
        for cell in ct.single_cells_of_type:
            if cell.score >= threshold:
                for gene in set(cell.genes_used_for_score):
                    gene_counts[gene] += 1

        ct.gene_occurrences = [
            GeneOccurrence(gene=g, occurrence=c) for g, c in gene_counts.most_common()
        ]
