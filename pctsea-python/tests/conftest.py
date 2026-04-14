"""Shared test fixtures for PCTSEA tests."""

from __future__ import annotations

from pathlib import Path

import pytest

from pctsea.models import SingleCell


@pytest.fixture
def sample_input_path() -> Path:
    return Path(__file__).parent / "fixtures" / "sample_input.txt"


@pytest.fixture
def sample_cells() -> list[SingleCell]:
    """Create a small set of test single cells with known expressions."""
    cells = []
    for i in range(20):
        cell = SingleCell(id=i, name=f"cell_{i}")
        # First 10 cells are "type_a", last 10 are "type_b"
        if i < 10:
            cell.cell_type = "type_a"
            cell.cell_type_id = 1
        else:
            cell.cell_type = "type_b"
            cell.cell_type_id = 2

        # Gene 0: expressed in all cells
        cell.add_gene_expression(0, float(i + 1))
        # Gene 1: expressed more in type_a
        if i < 10:
            cell.add_gene_expression(1, float(10 - i))
        else:
            cell.add_gene_expression(1, float(i - 9) * 0.1)
        # Gene 2: expressed in some cells
        if i % 2 == 0:
            cell.add_gene_expression(2, float(i * 0.5 + 1))

        cells.append(cell)
    return cells


@pytest.fixture
def experiment_expressions() -> dict[int, float]:
    """Sample experiment expression values for 3 genes."""
    return {0: 5.0, 1: 3.0, 2: 2.0}


@pytest.fixture
def gene_id_to_name() -> dict[int, str]:
    return {0: "GENE_A", 1: "GENE_B", 2: "GENE_C"}
