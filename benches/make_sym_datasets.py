#!/usr/bin/env python3
"""
Create symmetric datasets by appending opposite values to each line.

Example:
  "0;23;14;24" -> "0;23;14;24;-23;-14;-24"
"""

from __future__ import annotations

import argparse
from pathlib import Path


def ensure_output_root(output_root: Path, overwrite: bool) -> None:
    if output_root.exists() and any(output_root.iterdir()):
        if not overwrite:
            raise SystemExit(
                f"Output root '{output_root}' exists and is not empty. "
                "Use --overwrite or choose a new folder."
            )


def process_line(line: str) -> str:
    stripped = line.strip()
    if not stripped:
        return ""

    parts = [p.strip() for p in stripped.split(";")]
    values = []
    for p in parts:
        try:
            values.append(int(p))
        except ValueError:
            # If the line isn't purely numeric, keep it unchanged.
            return stripped

    existing = set(values)
    out = list(values)
    for v in values:
        if v == 0:
            continue
        neg = -v
        if neg not in existing:
            out.append(neg)
            existing.add(neg)

    return ";".join(str(v) for v in out)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create symmetric datasets by appending opposite values."
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        default=Path("benches/dataSets/2ADD"),
        help="Input dataset root.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("benches/symDataSets"),
        help="Output dataset root.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow writing into an existing, non-empty output folder.",
    )
    args = parser.parse_args()

    input_root = args.input_root
    output_root = args.output_root
    if not input_root.exists():
        raise SystemExit(f"Input root not found: {input_root}")

    ensure_output_root(output_root, args.overwrite)

    for path in input_root.rglob("*"):
        rel = path.relative_to(input_root)
        out_path = output_root / rel
        if path.is_dir():
            out_path.mkdir(parents=True, exist_ok=True)
            continue
        if path.suffix.lower() != ".txt":
            continue

        out_path.parent.mkdir(parents=True, exist_ok=True)
        with path.open() as fin, out_path.open("w") as fout:
            for line in fin:
                fout.write(process_line(line) + "\n")


if __name__ == "__main__":
    main()
