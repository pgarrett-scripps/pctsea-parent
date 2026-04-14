"""Tests for cell type ontology mapping."""

from __future__ import annotations

from pctsea.cell_types import (
    get_branched,
    get_cell_type_id,
    get_cell_type_name,
    get_original_cell_types,
    parse_cell_type_typos,
)


def test_load_cell_types() -> None:
    """Test that the ontology file loads successfully."""
    types = get_original_cell_types()
    assert len(types) > 600  # The file has 668 entries


def test_get_branched() -> None:
    """Test hierarchical lookup."""
    branched = get_branched("activated t cell")
    assert branched is not None
    assert branched.original == "activated t cell"
    assert branched.type_ is not None


def test_cell_type_id_assignment() -> None:
    """Test that IDs are assigned consistently."""
    id1 = get_cell_type_id("test_type_a")
    id2 = get_cell_type_id("test_type_b")
    id1_again = get_cell_type_id("test_type_a")
    assert id1 != id2
    assert id1 == id1_again
    assert get_cell_type_name(id1) == "test_type_a"


def test_parse_typos() -> None:
    """Test typo correction."""
    assert parse_cell_type_typos("activative t cell") == "activated t cell"
    assert parse_cell_type_typos("unknown1") == "unknown"
    assert parse_cell_type_typos("epithelial") == "epithelial cell"
    assert parse_cell_type_typos("kerationcyte") == "keratinocyte"
    assert parse_cell_type_typos("beta cell") == "b cell"
