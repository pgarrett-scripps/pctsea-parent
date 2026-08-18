# Product direction

The Rust project is a scientifically improved successor to legacy PCTSEA, not a
numerical reimplementation. Legacy results are useful as regression references,
but clearer statistics, better validation, and reproducibility take priority over
matching historical output.

## Version-one scope

- Human data only.
- The Human Cell Landscape (HCL) is the first supported atlas.
- H5AD is the public expression format because it is distributed by the atlas
  authors. The published HCL Figure 1 H5AD uses an older dense, chunked layout
  and stores tissue/batch but not cell type, so the paired official
  `HCL_Fig1_cell_Info.xlsx` supplies cell-type annotations.
- The CLI will be able to acquire a known, versioned HCL release and will also
  accept a user-provided H5AD path.
- The deliverable is a Rust library and CLI. A web or desktop interface is out of
  scope until the scientific workflow is validated.

## Storage policy

Users do not choose an internal storage technology. Expression remains in H5AD
and is read by selected gene column. The CLI automatically builds a compact
binary cache of the reconciled cell metadata because repeatedly parsing the
600,000-row Excel workbook is slow. Parquet is unnecessary for this first atlas;
the prototype `.pctsea` format remains available for portable custom imports.

## Validation policy

Validation will cover input integrity, deterministic execution, known statistical
cases, and several legacy PCTSEA analyses. Differences from legacy behavior are
acceptable when documented and scientifically justified.
