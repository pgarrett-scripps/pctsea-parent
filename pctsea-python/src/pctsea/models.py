"""Data models for PCTSEA."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from pctsea.config import CellTypeBranch, ScoringMethod


@dataclass
class ScoringSchema:
    """A scoring method with its threshold and minimum genes per cell."""

    method: ScoringMethod
    threshold: float
    min_genes_cells: int


@dataclass
class GeneOccurrence:
    """Tracks how many cells use a particular gene in a cell type."""

    gene: str
    occurrence: int = 0

    def increment(self) -> None:
        self.occurrence += 1


@dataclass
class CellTypeBranched:
    """Hierarchical cell type classification with multiple levels."""

    original: str
    type_: str | None = None
    subtype: str | None = None
    characteristic: str | None = None

    def get_branch(self, branch: CellTypeBranch) -> str:
        if branch == CellTypeBranch.ORIGINAL:
            return self.original
        if branch == CellTypeBranch.CHARACTERISTIC:
            return self.characteristic or ""
        if branch == CellTypeBranch.TYPE:
            return self.type_ or ""
        if branch == CellTypeBranch.TYPE_SUBTYPE:
            ret = self.type_ or ""
            if self.subtype:
                ret += "-" + self.subtype
            return ret
        if branch == CellTypeBranch.TYPE_SUBTYPE_CHARACTERISTIC:
            ret = self.type_ or ""
            if self.subtype:
                ret += "-" + self.subtype
            if self.characteristic:
                ret += "-" + self.characteristic
            return ret
        msg = f"{branch} not supported"
        raise ValueError(msg)


@dataclass
class SingleCell:
    """A single cell with sparse gene expression data and scoring results."""

    id: int
    name: str | None = None
    cell_type: str | None = None
    original_cell_type: str | None = None
    cell_type_id: int = -1
    score: float = float("nan")
    expressions: dict[int, float] = field(default_factory=dict)
    genes_used_for_score: list[str] = field(default_factory=list)
    dataset_tag: str | None = None
    biomaterial: str | None = None

    def add_gene_expression(self, gene_id: int, value: float) -> None:
        existing = self.expressions.get(gene_id, 0.0)
        if value > existing:
            self.expressions[gene_id] = value

    def get_gene_expression(self, gene_id: int) -> float:
        return self.expressions.get(gene_id, 0.0)

    def num_expressed_genes(self, gene_ids: list[int]) -> int:
        return sum(1 for gid in gene_ids if self.expressions.get(gid, 0.0) > 0.0)


@dataclass
class Gene:
    """A gene with expression data across single cells."""

    gene_id: int
    name: str
    expressions_by_cell_id: dict[int, float] = field(default_factory=dict)
    num_cells_by_type: dict[str, int] = field(default_factory=dict)

    def add_expression(
        self, cell_id: int, value: float, cell_type_name: str | None = None
    ) -> None:
        existing = self.expressions_by_cell_id.get(cell_id, 0.0)
        if value > existing:
            self.expressions_by_cell_id[cell_id] = value
        if cell_type_name:
            self.num_cells_by_type[cell_type_name] = (
                self.num_cells_by_type.get(cell_type_name, 0) + 1
            )

    def get_num_cells_expressing(self, cell_type_name: str) -> int:
        return self.num_cells_by_type.get(cell_type_name, 0)


@dataclass
class CellTypeClassification:
    """Result for a single cell type with all enrichment metrics."""

    name: str
    cell_type_id: int = -1

    # Hypergeometric test
    hypergeometric_pvalue: float = 1.0
    casimirs_enrichment_score: float = 0.0
    num_cells_of_type: int = 0
    num_cells_of_type_passing_threshold: int = 0

    # Enrichment scores
    enrichment_score: float = float("nan")
    enrichment_unweighted_score: float = 0.0
    secondary_enrichment_score: float | None = None
    secondary_supremum_x: int | None = None

    # KS test statistics
    d_statistic: float = 0.0
    ks_pvalue: float = 1.0
    ks_corrected_pvalue: float = 1.0
    supremum_x: int = 0
    normalized_supremum_x: float = 0.0
    size_a: int = 0
    size_b: int = 0

    # Distribution tracking
    num_cell_type_scores: int = 0
    num_other_cell_type_scores: int = 0

    # Significance
    enrichment_significance: float = 1.0
    enrichment_fdr: float = 1.0
    significance_string: str = ""

    # Normalized scores
    _normalized_enrichment_score: float | None = None
    normalized_random_enrichment_scores: list[float] = field(default_factory=list)

    # Permutation results
    random_ks_statistics: list[float] = field(default_factory=list)
    random_enrichment_scores: list[float] = field(default_factory=list)

    # UMAP
    umap_components: list[float] | None = None

    # Gene contributions
    gene_occurrences: list[GeneOccurrence] = field(default_factory=list)
    num_genes_significant: int = 0

    # Associated cells
    single_cells_of_type: list[SingleCell] = field(default_factory=list)

    def set_enrichment(
        self,
        supremum: float,
        d_statistic: float,
        supremum_x: int,
        normalized_supremum_x: float,
    ) -> None:
        self.enrichment_score = supremum
        self.d_statistic = d_statistic
        self.supremum_x = supremum_x
        self.normalized_supremum_x = normalized_supremum_x

    def add_random_enrichment(self, supremum: float, d_statistic: float) -> None:
        self.random_enrichment_scores.append(supremum)
        self.random_ks_statistics.append(d_statistic)

    def get_normalized_enrichment_score(self) -> float:
        if self._normalized_enrichment_score is not None:
            return self._normalized_enrichment_score

        if math.isnan(self.enrichment_score):
            self._normalized_enrichment_score = float("nan")
            return self._normalized_enrichment_score

        positive_randoms = [s for s in self.random_enrichment_scores if s >= 0.0]
        negative_randoms = [s for s in self.random_enrichment_scores if s < 0.0]

        if self.enrichment_score >= 0.0:
            if not positive_randoms:
                self._normalized_enrichment_score = float("nan")
                return self._normalized_enrichment_score
            expected = abs(sum(positive_randoms) / len(positive_randoms))
        else:
            if not negative_randoms:
                self._normalized_enrichment_score = float("nan")
                return self._normalized_enrichment_score
            expected = abs(sum(negative_randoms) / len(negative_randoms))

        if expected == 0.0:
            self._normalized_enrichment_score = float("nan")
            return self._normalized_enrichment_score

        self._normalized_enrichment_score = self.enrichment_score / expected

        # Also normalize the random scores
        self.normalized_random_enrichment_scores = [
            r / expected for r in self.random_enrichment_scores
        ]

        return self._normalized_enrichment_score

    def set_name(self, name: str) -> None:
        name = name.strip()
        try:
            float(name)
            name = "_" + name
        except ValueError:
            pass
        if not name:
            name = "_"
        self.name = name

    def get_gene_ranking_string(self) -> str:
        parts = []
        for go in self.gene_occurrences:
            parts.append(f"{go.gene}[{go.occurrence}]")
        return ",".join(parts)


@dataclass
class SingleCellSet:
    """Container for single cells with lookup by ID and name."""

    cells_by_id: dict[int, SingleCell] = field(default_factory=dict)
    cell_id_by_name: dict[str, int] = field(default_factory=dict)
    total_num_cells_for_dataset: int = 0

    def add_single_cell(self, cell: SingleCell) -> None:
        self.cells_by_id[cell.id] = cell
        if cell.name:
            self.cell_id_by_name[cell.name] = cell.id

    def get_cell_by_id(self, cell_id: int) -> SingleCell | None:
        return self.cells_by_id.get(cell_id)

    def get_cell_id_by_name(self, name: str) -> int:
        return self.cell_id_by_name.get(name, -1)

    @property
    def cell_list(self) -> list[SingleCell]:
        return list(self.cells_by_id.values())

    @property
    def num_cells(self) -> int:
        return len(self.cells_by_id)

    def clear(self) -> None:
        self.cells_by_id.clear()
        self.cell_id_by_name.clear()
