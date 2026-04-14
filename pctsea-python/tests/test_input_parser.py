"""Tests for input file parsing."""

from __future__ import annotations

from pathlib import Path

from pctsea.input_parser import read_expression_file


def test_read_expression_file(sample_input_path: Path) -> None:
    expressions = read_expression_file(sample_input_path)
    assert len(expressions) == 10
    assert expressions["TP53"] == 2.45
    assert expressions["KRAS"] == 3.50
    assert "BRCA1" in expressions


def test_read_expression_file_gene_names_uppercased(tmp_path: Path) -> None:
    f = tmp_path / "test.txt"
    f.write_text("Gene\tExpression\ntp53\t1.5\negfr\t2.0\n")
    expressions = read_expression_file(f)
    assert "TP53" in expressions
    assert "EGFR" in expressions
