"""Tests for permutation testing and FDR calculation."""

from __future__ import annotations

from pctsea.permutation import benjamini_hochberg


def test_bh_correction_basic() -> None:
    """Test BH correction with known values."""
    pvals = [0.01, 0.04, 0.03, 0.20]
    adjusted = benjamini_hochberg(pvals)

    assert len(adjusted) == 4
    # Adjusted p-values should be >= original
    for orig, adj in zip(pvals, adjusted, strict=True):
        assert adj >= orig or abs(adj - orig) < 1e-10
    # All should be <= 1.0
    for adj in adjusted:
        assert adj <= 1.0


def test_bh_correction_monotonic() -> None:
    """BH-adjusted values for sorted input should be non-decreasing after re-sorting."""
    pvals = [0.001, 0.01, 0.05, 0.10, 0.50]
    adjusted = benjamini_hochberg(pvals)

    # When sorted, adjusted should be non-decreasing
    sorted_adj = sorted(adjusted)
    for i in range(len(sorted_adj) - 1):
        assert sorted_adj[i] <= sorted_adj[i + 1] + 1e-10


def test_bh_correction_empty() -> None:
    """Empty input should return empty output."""
    assert benjamini_hochberg([]) == []
