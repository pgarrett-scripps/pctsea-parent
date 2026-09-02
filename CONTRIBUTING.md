# Contributing

Install the stable Rust toolchain. The repository toolchain file also installs
Rustfmt and Clippy.

Run the local quality gate before opening a pull request:

```bash
cargo fmt --package pctsea --check
cargo clippy --locked --all-targets -- -D warnings
cargo test --locked --all-targets
python3 -m unittest discover -s tools -p 'test_*.py' -v
cargo build --locked --release
```

The official HCL integration test is ignored by default because the reference
files are large. Run it when those files are available:

```bash
PCTSEA_HCL_PATH=validation/hcl/HCL_Fig1_adata.h5ad \
  cargo test --test hcl_official -- --ignored
```

## Release checklist

1. Update `package.version` in `Cargo.toml`.
2. Update `Cargo.lock` with `cargo check`.
3. Confirm the quality gate passes.
4. Merge the version change into `main`.
5. Tag that commit with the matching `v` prefix and push the tag.

For example, Cargo version `0.2.0` must use tag `v0.2.0`. The release workflow
rejects a mismatch before publishing anything.
