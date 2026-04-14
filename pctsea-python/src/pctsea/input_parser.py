"""Parse experimental expression input files."""

from __future__ import annotations

import logging
from pathlib import Path

log = logging.getLogger(__name__)


def read_expression_file(path: Path) -> dict[str, float]:
    """Read a tab-separated experimental expression file.

    Expected format:
        Header line (skipped)
        gene_name<TAB>expression_value

    Returns a dict mapping uppercased gene names to expression values.
    """
    expressions: dict[str, float] = {}
    with open(path) as f:
        for line_num, line in enumerate(f, start=1):
            line = line.strip()
            if not line:
                continue
            if line_num == 1:
                # Skip header
                continue

            # Try tab first, then whitespace
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 2:
                log.warning("Skipping line %d: fewer than 2 columns", line_num)
                continue

            gene_name = parts[0].strip().upper()
            try:
                value = float(parts[1].strip())
            except ValueError:
                log.warning(
                    "Skipping line %d: cannot parse expression value '%s'",
                    line_num,
                    parts[1],
                )
                continue

            if gene_name in expressions:
                log.debug(
                    "Duplicate gene %s at line %d, keeping higher value", gene_name, line_num
                )
                expressions[gene_name] = max(expressions[gene_name], value)
            else:
                expressions[gene_name] = value

    log.info("Read %d genes from %s", len(expressions), path)
    return expressions
