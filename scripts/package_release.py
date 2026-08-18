#!/usr/bin/env python3
"""Create a platform release archive and SHA-256 checksum."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import shutil
import tarfile
import tempfile
import zipfile


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True, type=Path)
    parser.add_argument("--target", required=True)
    parser.add_argument("--version", required=True)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    if not args.binary.is_file():
        raise SystemExit(f"binary does not exist: {args.binary}")

    version = args.version.removeprefix("v")
    bundle_name = f"pctsea-{version}-{args.target}"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="pctsea-release-") as temp_dir:
        bundle = Path(temp_dir) / bundle_name
        bundle.mkdir()
        shutil.copy2(args.binary, bundle / args.binary.name)
        shutil.copy2("README.md", bundle / "README.md")
        shutil.copy2("LICENSE", bundle / "LICENSE")

        if args.target.endswith("windows-msvc"):
            archive = args.output_dir / f"{bundle_name}.zip"
            with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as output:
                for path in sorted(bundle.iterdir()):
                    output.write(path, f"{bundle_name}/{path.name}")
        else:
            archive = args.output_dir / f"{bundle_name}.tar.gz"
            with tarfile.open(archive, "w:gz") as output:
                output.add(bundle, arcname=bundle_name)

    checksum = file_sha256(archive)
    checksum_path = archive.with_name(f"{archive.name}.sha256")
    checksum_path.write_text(f"{checksum}  {archive.name}\n", encoding="utf-8")
    print(archive)
    print(checksum_path)


if __name__ == "__main__":
    main()
