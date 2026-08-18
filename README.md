# pCTSEA

[![CI](https://github.com/pgarrett-scripps/pctsea-parent/actions/workflows/ci.yml/badge.svg)](https://github.com/pgarrett-scripps/pctsea-parent/actions/workflows/ci.yml)

pCTSEA is a Rust command-line tool and library for proteomic cell-type
enrichment analysis. It ranks single cells by similarity to a quantitative
gene or protein query, then tests which cell types are concentrated near the
top of the ranking.

The tool reads H5AD/HDF5 directly. It does not require Java, R, MongoDB, or a
running service.

## Install

Download the archive for your computer from the
[GitHub Releases page](https://github.com/pgarrett-scripps/pctsea-parent/releases):

| Computer | Release target |
| --- | --- |
| Linux x86-64 | `x86_64-unknown-linux-gnu` |
| Windows x86-64 | `x86_64-pc-windows-msvc` |
| Intel Mac | `x86_64-apple-darwin` |
| Apple Silicon Mac | `aarch64-apple-darwin` |

Each archive has a matching `.sha256` checksum file. The executable is named
`pctsea` on Linux and macOS, and `pctsea.exe` on Windows.

To build from source, install stable Rust and run:

```bash
cargo build --release --locked
./target/release/pctsea --version
```

## Quick start

Download the official Human Cell Landscape expression matrix and its paired
cell annotations:

```bash
pctsea atlas download hcl --output ./data/HCL_Fig1_adata.h5ad
```

Run a query:

```bash
pctsea analyze \
  --atlas ./data/HCL_Fig1_adata.h5ad \
  --input examples/t_cell_query.tsv \
  --permutations 1000 \
  --output results.tsv
```

The HCL download is about 811 MiB. If `--output` is omitted during download,
pCTSEA uses `PCTSEA_DATA_DIR`, the platform data directory, or a local
`.pctsea` fallback. Run `pctsea atlas path hcl` to print the resolved paths.

## Query format

Input is a tab-separated gene and quantitative value, with an optional header:

```text
gene    value
CD3D    10
CD3E    9
TRBC1   8
```

Gene names are matched without regard to case. The score can be Pearson
correlation, cosine similarity, or dot product. Run `pctsea help` for all
analysis options.

Protein-level inputs should be converted into this small canonical query
format before analysis. The validation adapter in
[`tools/prepare_proteomics_queries.py`](tools/prepare_proteomics_queries.py)
maps UniProt accessions and can compute length-corrected NSAF values from
spectral counts. See [`tools/README.md`](tools/README.md) for its assumptions.

## Supported atlas formats

- Modern H5AD files
- The older dense Human Cell Landscape H5AD plus its official annotation file
- Native sparse `.pctsea` files built from portable TSV data

The native format is an optional import and cache format. H5AD remains the
normal public atlas format.

## Scientific status

This is a new implementation, not a bit-for-bit port of the historical Java
workflow. It uses weighted running-sum enrichment, deterministic label
permutations, plus-one corrected p-values, and Benjamini-Hochberg FDRs. Results
should be validated on reference datasets before publication.

The scientific scope and validation policy are in
[`DIRECTION.md`](DIRECTION.md).

## Development and releases

Pull requests and pushes to `main` run formatting, Clippy, Rust tests, Python
adapter tests, a release build, and native tests on Windows and both Mac
architectures.

To publish a release:

1. Update the version in `Cargo.toml`.
2. Merge the change into `main`.
3. Create and push the matching tag, such as `v0.1.0`.

GitHub Actions verifies the tag, runs native release tests, builds all four
archives, generates checksums, and creates the GitHub release. More local
development commands are in [`CONTRIBUTING.md`](CONTRIBUTING.md).

## License

Apache License 2.0. See [`LICENSE`](LICENSE).
