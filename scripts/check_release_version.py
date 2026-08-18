#!/usr/bin/env python3
"""Verify that a release tag matches the Cargo package version."""

from __future__ import annotations

import argparse
from pathlib import Path
import tomllib


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("tag", help="Release tag such as v0.1.0")
    args = parser.parse_args()

    tag_version = args.tag.removeprefix("v")
    with Path("Cargo.toml").open("rb") as cargo_file:
        package_version = tomllib.load(cargo_file)["package"]["version"]

    if args.tag == tag_version or tag_version != package_version:
        raise SystemExit(
            f"release tag {args.tag!r} must be v{package_version} to match Cargo.toml"
        )

    print(f"release version {package_version} verified")


if __name__ == "__main__":
    main()
