#!/usr/bin/env python3
"""
Run the solver on snapshot.rscm cases and summarize LUT estimates vs synthesis.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import tempfile
from pathlib import Path
from typing import Optional


def resolve_cases_root(
    cases_root: Optional[Path],
    dataset_root: Optional[Path],
    bits: Optional[str],
    constants: Optional[str],
) -> Path:
    if cases_root is not None:
        return cases_root
    if dataset_root is None or bits is None or constants is None:
        raise SystemExit("Provide --cases-root or --dataset-root + --bits + --constants.")
    return dataset_root / bits / constants


def load_synth_luts(cost_path: Path) -> Optional[float]:
    try:
        data = json.loads(cost_path.read_text())
    except Exception:
        return None
    val = data.get("synthetised_luts")
    if not isinstance(val, (int, float)) or not math.isfinite(val):
        return None
    return float(val)


def run_solver(solver: Path, snapshot: Path, tmp_json: Path) -> Optional[float]:
    snapshot_abs = snapshot.resolve()
    snapshot_arg = f"//{str(snapshot_abs).lstrip('/')}"
    args = [
        str(solver),
        f"--recompute-snapshot={snapshot_arg}",
        f"--json={tmp_json}",
    ]
    result = subprocess.run(args, capture_output=True, text=True)
    if result.returncode != 0:
        return None
    try:
        solution = json.loads(tmp_json.read_text())
    except Exception:
        return None
    est = solution.get("meta", {}).get("costs", {}).get("luts")
    if not isinstance(est, (int, float)) or not math.isfinite(est):
        return None
    return float(est)


def print_table(rows: list[tuple[str, float, float, float, float]]) -> None:
    header = ("case", "estimated", "synthetised", "relative_%", "diff_luts")
    widths = [len(h) for h in header]
    for batch, est, synth, rel, diff in rows:
        widths[0] = max(widths[0], len(batch))
        widths[1] = max(widths[1], len(f"{est:.2f}"))
        widths[2] = max(widths[2], len(f"{synth:.2f}"))
        widths[3] = max(widths[3], len(f"{rel:.2f}"))
        widths[4] = max(widths[4], len(f"{diff:.2f}"))

    fmt = (
        f"{{:{widths[0]}}}  "
        f"{{:>{widths[1]}}}  "
        f"{{:>{widths[2]}}}  "
        f"{{:>{widths[3]}}}  "
        f"{{:>{widths[4]}}}"
    )
    print(fmt.format(*header))
    print(
        fmt.format(
            "-" * widths[0],
            "-" * widths[1],
            "-" * widths[2],
            "-" * widths[3],
            "-" * widths[4],
        )
    )
    for batch, est, synth, rel, diff in rows:
        print(fmt.format(batch, f"{est:.2f}", f"{synth:.2f}", f"{rel:.2f}", f"{diff:.2f}"))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Summarize estimated vs synthesized LUTs in batches."
    )
    parser.add_argument(
        "--solver",
        type=Path,
        default=Path("/home/smith/CLionProjects/rscmSolverCpp/src/cmake-build-release/rscmSolver"),
        help="Path to the solver executable.",
    )
    parser.add_argument(
        "--cases-root",
        type=Path,
        help="Root folder containing case_### subfolders.",
    )
    parser.add_argument(
        "--dataset-root",
        type=Path,
        default=Path("benches/bench_new_luts_recomputed_test"),
        help="Dataset root if --cases-root is not provided.",
    )
    parser.add_argument(
        "--bits",
        help="Bits folder name (e.g. 8bits) if using --dataset-root.",
    )
    parser.add_argument(
        "--constants",
        help="Constants folder name (e.g. 8on8bitsR100) if using --dataset-root.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="Optional maximum number of cases to report.",
    )
    args = parser.parse_args()

    cases_root = resolve_cases_root(args.cases_root, args.dataset_root, args.bits, args.constants)
    if not cases_root.exists():
        raise SystemExit(f"Cases root not found: {cases_root}")

    solver = args.solver
    if not solver.exists():
        raise SystemExit(f"Solver not found: {solver}")

    cases = sorted([p for p in cases_root.glob("case_*") if p.is_dir()])
    if not cases:
        raise SystemExit(f"No case_* folders under {cases_root}")

    rows: list[tuple[str, float, float, float, float]] = []

    for idx, case_dir in enumerate(cases, start=1):
        if args.limit is not None and idx > args.limit:
            break
        snapshot = case_dir / "snapshot.rscm"
        cost_path = case_dir / "cost.json"
        if not snapshot.exists() or not cost_path.exists():
            continue

        synth = load_synth_luts(cost_path)
        if synth is None:
            continue

        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_json = Path(tmp_dir) / "solution.json"
            est = run_solver(solver, snapshot, tmp_json)

        if est is None:
            continue

        if synth != 0:
            rel = (est - synth) / synth * 100.0
        else:
            rel = 0.0
        diff = est - synth
        rows.append((case_dir.name, est, synth, rel, diff))

    if not rows:
        raise SystemExit("No valid cases to summarize.")

    print_table(rows)


if __name__ == "__main__":
    main()
