#!/usr/bin/env python3
"""
Scan a range of integer targets with `rpag` and record those with a desired
number of adders. The script runs `rpag <value>`, parses `no_of_adders` from
stdout, and writes matching values to an output file (one per line).

Usage: ./filter_adders.py <min_value> <max_value> <adders> [-o output_file] 
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path
import re


NO_OF_ADDERS_RE = re.compile(r"no_of_adders\s*=\s*(\d+)", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run `rpag` for each integer in the given range and keep the ones "
            "whose output reports the requested number of adders."
        )
    )
    parser.add_argument("min_value", type=int, help="Inclusive start of the target range.")
    parser.add_argument("max_value", type=int, help="Inclusive end of the target range.")
    parser.add_argument(
        "adders",
        type=int,
        help="Only keep values for which `no_of_adders` equals this count.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("matching_values.txt"),
        help="File to write matching values to (default: matching_values.txt).",
    )
    parser.add_argument(
        "--rpag-cmd",
        default="rpag",
        help="rpag executable to invoke (default: rpag).",
    )
    return parser.parse_args()


def extract_no_of_adders(output: str) -> int | None:
    """Return the parsed `no_of_adders` from rpag output, or None if missing."""
    match = NO_OF_ADDERS_RE.search(output)
    return int(match.group(1)) if match else None


def run_rpag(target: int, rpag_cmd: str) -> tuple[int | None, str]:
    """Run rpag for the given target and return (no_of_adders, combined_output)."""
    try:
        result = subprocess.run(
            [rpag_cmd, str(target)],
            capture_output=True,
            text=True,
            check=False,
        )
    except FileNotFoundError:
        print(f"Error: rpag command '{rpag_cmd}' not found.", file=sys.stderr)
        sys.exit(1)

    combined_output = (result.stdout or "") + (result.stderr or "")
    if result.returncode != 0:
        print(f"[warn] rpag exited with {result.returncode} for target {target}", file=sys.stderr)
    return extract_no_of_adders(combined_output), combined_output


def main() -> None:
    args = parse_args()
    if args.min_value > args.max_value:
        print("Error: min_value must be less than or equal to max_value.", file=sys.stderr)
        sys.exit(1)

    matches: list[int] = []
    for value in range(args.min_value, args.max_value + 1):
        no_of_adders, output = run_rpag(value, args.rpag_cmd)
        if no_of_adders is None:
            print(f"[warn] Could not find no_of_adders in output for {value}", file=sys.stderr)
            continue
        if no_of_adders <= args.adders:
            matches.append(value)

    if matches:
        args.output.write_text("\n".join(str(v) for v in matches) + "\n")
        print(f"Wrote {len(matches)} matching value(s) to {args.output}")
    else:
        # still create/clear the file so output is predictable
        args.output.write_text("")
        print("No values matched; output file left empty.")


if __name__ == "__main__":
    main()
