#!/usr/bin/env python3
"""Create or verify SHA256SUMS; verification also works without Git."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "SHA256SUMS"


def digest(path):
    with path.open("rb") as handle:
        result = hashlib.sha256()
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true", help="maintainer-only: freeze reviewed files"
    )
    args = parser.parse_args()
    if args.write:
        paths = (
            subprocess.check_output(
                ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
                cwd=ROOT,
            )
            .decode()
            .split("\0")
        )
        names = sorted(
            set(
                p
                for p in paths
                if p and p != "SHA256SUMS" and Path(p).name != ".DS_Store"
            )
        )
        lines = []
        for name in names:
            path = ROOT / name
            if path.is_symlink() or not path.is_file() or "\n" in name or "\r" in name:
                raise ValueError(f"Unsupported release file: {name}")
            lines.append(f"{digest(path)}  {name}\n")
        OUTPUT.write_text("".join(lines))
    lines = OUTPUT.read_text().splitlines()
    seen = set()
    for line in lines:
        expected, name = line.split("  ", 1)
        path = ROOT / name
        if (
            name in seen
            or path.is_symlink()
            or not path.resolve().is_relative_to(ROOT)
            or len(expected) != 64
        ):
            raise ValueError(f"Invalid manifest entry: {name}")
        seen.add(name)
        if digest(path) != expected:
            raise ValueError(f"Checksum mismatch: {name}")
    if not lines:
        raise ValueError("Empty checksum manifest")
    print(f"PASS: {len(lines)} release file SHA-256 checksums")


if __name__ == "__main__":
    main()
