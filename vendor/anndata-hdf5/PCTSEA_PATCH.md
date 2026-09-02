# Local Windows compatibility patch

This directory contains `anndata-hdf5` 0.5.3 under its MIT license.

The published crate enables both the `static` and `threadsafe` features of
`hdf5-metno-sys` on every platform. Native HDF5 rejects that combination when
building a static library with MSVC.

The local patch keeps `static`, `zlib`, and `threadsafe` on non-Windows
platforms. Windows keeps `static` and `zlib` but omits the unsupported native
thread-safety option. No Rust source code has been changed.

Remove this patch when the upstream crate makes its HDF5 features
platform-specific.
