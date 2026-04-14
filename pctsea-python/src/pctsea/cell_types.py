"""Cell type ontology mapping and hierarchical classification."""

from __future__ import annotations

import importlib.resources
import logging

from pctsea.config import CellTypeBranch
from pctsea.models import CellTypeBranched

log = logging.getLogger(__name__)

# Module-level state
_by_original: dict[str, CellTypeBranched] = {}
_by_type: dict[str, set[str]] = {}
_by_subtype: dict[str, set[str]] = {}
_by_characteristic: dict[str, set[str]] = {}
_original_cell_types: set[str] = set()
_loaded = False

_cell_type_ids: dict[str, int] = {}
_cell_type_names: dict[int, str] = {}
_next_cell_type_id = 1

# Typo correction cache
_typo_cache: dict[str, str] = {}


def _ensure_loaded() -> None:
    global _loaded
    if _loaded:
        return
    _load_hierarchical_cell_types()
    _loaded = True


def _load_hierarchical_cell_types() -> None:
    """Load the spike_cell_types_mapping_CB.txt resource file."""
    data_files = importlib.resources.files("pctsea") / "data"
    mapping_file = data_files / "spike_cell_types_mapping_CB.txt"
    text = mapping_file.read_text(encoding="latin-1")
    for num_line, line in enumerate(text.splitlines(), start=1):
        if num_line == 1:
            continue
        if not line.strip():
            continue
        parts = line.split("\t")
        original = parts[0].strip()
        _original_cell_types.add(original)

        type_ = None
        subtype = None
        characteristic = None

        if len(parts) > 3:
            type_ = parts[3].strip().replace(" ", "_") or None
            if type_:
                _by_type.setdefault(type_, set()).add(original)

        if len(parts) > 4:
            subtype = parts[4].strip() or None
            if subtype:
                _by_subtype.setdefault(subtype, set()).add(original)

        if len(parts) > 5:
            characteristic = parts[5].strip() or None
            if characteristic:
                _by_characteristic.setdefault(characteristic, set()).add(original)

        if original not in _by_original:
            _by_original[original] = CellTypeBranched(
                original=original,
                type_=type_,
                subtype=subtype,
                characteristic=characteristic,
            )


def get_cell_type_id(cell_type: str) -> int:
    """Get or create a numeric ID for a cell type name."""
    global _next_cell_type_id
    if cell_type in _cell_type_ids:
        return _cell_type_ids[cell_type]
    cid = _next_cell_type_id
    _next_cell_type_id += 1
    _cell_type_ids[cell_type] = cid
    _cell_type_names[cid] = cell_type
    return cid


def get_cell_type_name(cell_type_id: int) -> str | None:
    """Reverse lookup: ID -> name."""
    return _cell_type_names.get(cell_type_id)


def get_branched(original_type: str) -> CellTypeBranched | None:
    """Look up the hierarchical classification for an original cell type."""
    _ensure_loaded()
    if original_type is None:
        return None
    return _by_original.get(original_type.strip())


def parse_cell_type_typos(cell_type: str) -> str:
    """Fix common typos and normalize cell type names.

    Port of SingleCell.parseCellTypeTypos() from Java.
    """
    if cell_type in _typo_cache:
        return _typo_cache[cell_type]

    original = cell_type

    # Parse heterogeneity
    if "_" in cell_type:
        cell_type = cell_type.split("_")[0].strip()

    # Fix specific typos
    if cell_type in ("activative t cell", "actived t cell"):
        cell_type = "activated t cell"
    elif cell_type in ("unknown1", "unknown2"):
        cell_type = "unknown"

    cell_type = cell_type.replace("acinar", "acniar")
    cell_type = cell_type.replace("  ", " ")
    cell_type = cell_type.replace("soomth", "smooth")
    cell_type = cell_type.replace("muscel", "muscle")

    if cell_type.startswith("b cell") or " b cell" in cell_type:
        cell_type = "b cell"

    if cell_type == "antigen presenting cell":
        cell_type = "antigen-presenting cell"
    elif cell_type == "astrocyte(bergmann glia)":
        cell_type = "astrocyte"
    elif cell_type == "epithelial":
        cell_type = "epithelial cell"
    elif cell_type == "kerationcyte":
        cell_type = "keratinocyte"
    elif cell_type == "mast":
        cell_type = "mast cell"
    elif cell_type == "megakaryocyte/erythroid progenitor":
        cell_type = "megakaryocyte/erythtoid progenitor cell"
    elif cell_type == "syncytiotrophoblast":
        cell_type = "syncytiotrophoblast cell"
    elif cell_type == "neutriophil":
        cell_type = "neutrophil"
    elif cell_type == "kuppfer cell":
        cell_type = "kupffer cell"
    elif cell_type == "beta cell":
        cell_type = "b cell"

    cell_type = cell_type.strip()
    _typo_cache[original] = cell_type
    return cell_type


def resolve_cell_type(raw_type: str, branch: CellTypeBranch) -> tuple[str, int]:
    """Resolve a raw cell type string to its branch name and numeric ID.

    Returns (resolved_name, cell_type_id).
    """
    _ensure_loaded()
    cleaned = parse_cell_type_typos(raw_type)

    if branch == CellTypeBranch.ORIGINAL:
        cid = get_cell_type_id(cleaned)
        return cleaned, cid

    branched = get_branched(cleaned)
    if branched is None:
        cid = get_cell_type_id(cleaned)
        return cleaned, cid

    resolved = branched.get_branch(branch)
    if not resolved:
        resolved = cleaned
    cid = get_cell_type_id(resolved)
    return resolved, cid


def get_original_cell_types() -> set[str]:
    _ensure_loaded()
    return _original_cell_types.copy()
