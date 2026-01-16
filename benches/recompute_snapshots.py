#!/usr/bin/env python3
"""
Recompute costs from snapshot.rscm files using the solver binary.

This mirrors the snapshot folder structure into a fresh output root and writes:
  - solution.json (solver JSON dump)
  - console.txt (stdout/stderr)
  - params.json (run metadata)
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path
import shutil
from typing import List


def collect_snapshots(root: Path, name: str) -> List[Path]:
    return sorted(root.rglob(name))


def ensure_fresh_output(output_root: Path, overwrite: bool) -> None:
    if output_root.exists() and any(output_root.iterdir()):
        if not overwrite:
            raise SystemExit(
                f"Output root '{output_root}' exists and is not empty. "
                "Use --overwrite or choose a new folder."
            )


def run_snapshot(
    solver: Path,
    snapshot: Path,
    out_dir: Path,
    keep_going: bool,
    copy_snapshot: bool,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    console_path = out_dir / "console.txt"
    json_path = out_dir / "solution.json"
    params_path = out_dir / "params.json"
    snapshot_out = out_dir / snapshot.name

    if copy_snapshot:
        shutil.copy2(snapshot, snapshot_out)

    # Use an absolute path and add a leading slash for compatibility with
    # older binaries that trim one extra character from the snapshot flag.
    snapshot_abs = snapshot.resolve()
    snapshot_arg = f"//{str(snapshot_abs).lstrip('/')}"
    args = [
        str(solver),
        f"--recompute-snapshot={snapshot_arg}",
        f"--json={json_path}",
    ]

    start = time.perf_counter()
    result = subprocess.run(args, capture_output=True, text=True)
    duration = time.perf_counter() - start

    console_path.write_text(result.stdout + "\n" + result.stderr)

    params = {
        "snapshot": str(snapshot),
        "cmd": args,
        "returncode": result.returncode,
        "duration_seconds": duration,
    }
    params_path.write_text(json.dumps(params, indent=2))

    if result.returncode != 0 and not keep_going:
        raise SystemExit(f"Solver failed for snapshot: {snapshot}")


def aggregate_results(base_out_dir: Path, baseline_summary: Path | None) -> None:
    cost_sums: dict[str, float] = {}
    cost_counts: dict[str, int] = {}
    cost_names: set[str] = set()
    runtime_sum = 0.0
    runtime_count = 0

    for case_dir in sorted(base_out_dir.glob("case_*")):
        params_path = case_dir / "params.json"
        if params_path.exists():
            try:
                params = json.loads(params_path.read_text())
                runtime_sum += float(params.get("duration_seconds", 0.0))
                runtime_count += 1
            except Exception:
                pass

        solution_path = case_dir / "solution.json"
        if solution_path.exists():
            try:
                solution = json.loads(solution_path.read_text())
                costs = solution.get("meta", {}).get("costs", {})
                for name, val in costs.items():
                    cost_names.add(name)
                    if isinstance(val, (int, float)):
                        cost_sums[name] = cost_sums.get(name, 0.0) + float(val)
                        cost_counts[name] = cost_counts.get(name, 0) + 1
            except Exception:
                pass

    avg_costs: dict[str, object] = {}
    for name in cost_names:
        if cost_counts.get(name):
            avg_costs[name] = cost_sums[name] / cost_counts[name]
        else:
            avg_costs[name] = "not available"

    avg_runtime = (runtime_sum / runtime_count) if runtime_count else None
    if baseline_summary and baseline_summary.exists():
        try:
            baseline = json.loads(baseline_summary.read_text())
            if "avg_runtime_seconds" in baseline:
                avg_runtime = baseline.get("avg_runtime_seconds")
        except Exception:
            pass

    summary = {
        "runs": runtime_count,
        "avg_runtime_seconds": avg_runtime,
        "costs": avg_costs,
    }
    (base_out_dir / "summary.json").write_text(json.dumps(summary, indent=2))


def collect_case_roots(output_root: Path) -> List[Path]:
    parents: set[Path] = set()
    for case_dir in output_root.rglob("case_*"):
        if case_dir.is_dir():
            parents.add(case_dir.parent)
    return sorted(parents)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Recompute costs from snapshot files into a fresh output folder."
    )
    parser.add_argument(
        "--snapshots-root",
        type=Path,
        default=Path("benches/bench_new_luts"),
        help="Root folder containing snapshot.rscm files.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("benches/bench_new_luts_recomputed"),
        help="Fresh output folder to write solution.json and logs.",
    )
    parser.add_argument(
        "--solver",
        type=Path,
        default=Path("build/rscmSolver"),
        help="Path to the solver executable (rscmSolver/rscmBuilder).",
    )
    parser.add_argument(
        "--snapshot-name",
        default="snapshot.rscm",
        help="Snapshot filename to search for.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow writing into an existing, non-empty output folder.",
    )
    parser.add_argument(
        "--keep-going",
        action="store_true",
        help="Continue even if a snapshot fails.",
    )
    parser.add_argument(
        "--no-copy-snapshot",
        action="store_true",
        help="Do not copy snapshot.rscm into the output folder.",
    )
    parser.add_argument(
        "--no-summary",
        action="store_true",
        help="Skip recomputing summary.json files.",
    )
    args = parser.parse_args()

    solver = args.solver
    if not solver.exists():
        raise SystemExit(f"Solver executable not found: {solver}")

    snapshots = collect_snapshots(args.snapshots_root, args.snapshot_name)
    if not snapshots:
        raise SystemExit(f"No snapshots named '{args.snapshot_name}' under {args.snapshots_root}")

    ensure_fresh_output(args.output_root, args.overwrite)

    for idx, snapshot in enumerate(snapshots, start=1):
        rel = snapshot.relative_to(args.snapshots_root)
        out_dir = args.output_root / rel.parent
        run_snapshot(
            solver,
            snapshot,
            out_dir,
            args.keep_going,
            copy_snapshot=not args.no_copy_snapshot,
        )
        if idx % 25 == 0 or idx == len(snapshots):
            print(f"Recomputed {idx}/{len(snapshots)}")

    if not args.no_summary:
        for base_dir in collect_case_roots(args.output_root):
            rel = base_dir.relative_to(args.output_root)
            baseline = args.snapshots_root / rel / "summary.json"
            aggregate_results(base_dir, baseline)


if __name__ == "__main__":
    main()
