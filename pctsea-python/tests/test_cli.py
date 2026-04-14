"""Tests for CLI entry point."""

from __future__ import annotations

from click.testing import CliRunner

from pctsea.cli import main


def test_cli_help() -> None:
    """Test that --help works."""
    runner = CliRunner()
    result = runner.invoke(main, ["--help"])
    assert result.exit_code == 0
    assert "PCTSEA" in result.output
    assert "--eef" in result.output
    assert "--scoring-method" in result.output


def test_cli_missing_required() -> None:
    """Test that missing required args produces an error."""
    runner = CliRunner()
    result = runner.invoke(main, [])
    assert result.exit_code != 0
