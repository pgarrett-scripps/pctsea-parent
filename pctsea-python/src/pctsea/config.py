"""Enums and constants for PCTSEA."""

from __future__ import annotations

from enum import Enum


class ScoringMethod(Enum):
    """Scoring method for measuring similarities between input and single cell expressions."""

    SIMPLE_SCORE = ("SimpleScore", "Sum of matching genes and normalized intensity differences")
    PEARSONS_CORRELATION = (
        "correlation",
        "Pearson correlation between input protein abundances and single cell expressions",
    )
    DOT_PRODUCT = (
        "Dot-product",
        "Dot product between input protein expressions and single cell gene expressions",
    )
    QUICK_SCORE = (
        "Quick_score",
        "Product of detection factor (cells with gene / cells in type) per gene",
    )
    REGRESSION = (
        "Linear Regression",
        "R2 value of linear regression between protein abundances and cell gene expressions",
    )

    def __init__(self, score_name: str, description: str) -> None:
        self.score_name = score_name
        self.description = description

    def __str__(self) -> str:
        return self.score_name

    @classmethod
    def from_name(cls, name: str) -> ScoringMethod | None:
        for method in cls:
            if method.score_name == name or method.name == name:
                return method
        return None

    @classmethod
    def choices(cls) -> list[str]:
        return [m.name for m in cls]


class CellTypeBranch(Enum):
    """Cell type classification hierarchy level."""

    ORIGINAL = "ORIGINAL"
    TYPE = "TYPE"
    TYPE_SUBTYPE = "TYPE_SUBTYPE"
    TYPE_SUBTYPE_CHARACTERISTIC = "TYPE_SUBTYPE_CHARACTERISTIC"
    CHARACTERISTIC = "CHARACTERISTIC"

    @classmethod
    def from_name(cls, name: str) -> CellTypeBranch | None:
        for branch in cls:
            if branch.name.upper() == name.upper():
                return branch
        return None


class InputDataType(Enum):
    """Type of input proteomics data."""

    IP = "Immunoprecipitation"
    PROTEOME_OF_CELL_LINE = "Proteome of a cell line"
    PROTEOME_OF_SINGLE_CELL = "Proteome of a single cell"
    PROTEOME_OF_TISSUE = "Proteome of a tissue"
    PARTIAL_PROTEOME_OF_CELL_LINE = "Partial proteome of cell line"
    PARTIAL_PROTEOME_OF_TISSUE = "Partial proteome of a tissue"


# --- Constants ---

FDR_THRESHOLD: float = 0.05
MIN_CELLS_PASSING_THRESHOLD: int = 1000
UMAP_DIMENSIONS: int = 4
MIN_NUM_MAPPED_GENES: int = 20

# Default thresholds per scoring method
DEFAULT_MIN_SCORE_SIMPLE_SCORE: float = 0.0
DEFAULT_MIN_SCORE_PEARSON: float = 0.2
DEFAULT_MIN_GENES_CELLS_SIMPLE_SCORE: int = 10
DEFAULT_MIN_GENES_CELLS_PEARSON: int = 100
