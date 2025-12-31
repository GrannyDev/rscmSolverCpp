#!/usr/bin/env python3
"""
Benchmark runner for rmcmSolver over the datasets in benches/dataSets/2ADD/.

It mirrors the dataset hierarchy into a results folder and, for each target
set, creates a folder per test case containing:
  - console.txt : solver stdout/stderr
  - solution.json : solver JSON dump (if produced)
  - params.json : parameters passed to the solver plus timing

Usage examples:
    python3 benches/bench_runner.py --cost=mux_bits
    python3 benches/bench_runner.py --set timeout --max-cases 10 --cost=area
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path
from typing import List, Optional


def parse_int_list(spec: str) -> List[int]:
    return [int(x) for x in spec.split(",") if x]


def parse_cost_model(name: str) -> str:
    allowed = {"area", "area_cost", "mux_count", "mux_bits", "luts", "fpga_delay", "asic_delay"}
    if name not in allowed:
        raise ValueError(f"Unknown cost model '{name}'. Allowed: {', '.join(sorted(allowed))}")
    # keep backward-compatible alias
    return "area_cost" if name == "area" else name


def extract_beta_from_folder(bits_folder: Path) -> int:
    # Folder names look like "10bits"; beta is folder number + 1
    digits = "".join(ch for ch in bits_folder.name if ch.isdigit())
    if not digits:
        raise ValueError(f"Cannot parse beta from folder name '{bits_folder.name}'")
    return int(digits) + 1


def load_targets(file_path: Path, max_cases: Optional[int]) -> List[str]:
    targets: List[str] = []
    with file_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # dataset uses ';' separators, solver expects comma-separated ints
            targets.append(line.replace(";", ","))
            if max_cases is not None and len(targets) >= max_cases:
                break
    return targets


def run_case(
    solver_path: Path,
    out_dir: Path,
    beta: int,
    nb_input_bits: int,
    layout: str,
    targets: str,
    cost: str,
    heuristic: Optional[int],
    timeout_minutes: Optional[int],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    console_path = out_dir / "console.txt"
    json_path = out_dir / "solution.json"
    params_path = out_dir / "params.json"

    args = [
        str(solver_path),
        f"--beta={beta}",
        f"--nb-input-bits={nb_input_bits}",
        f"--layout={layout}",
        f"--targets={targets}",
        f"--cost={cost}",
        f"--json={json_path}",
    ]
    if heuristic is not None:
        args.append(f"--heuristic={heuristic}")
    timeout_seconds = timeout_minutes * 60 if timeout_minutes is not None else None
    if timeout_seconds is not None:
        args.append(f"--timeout={timeout_seconds}")

    start = time.perf_counter()
    # Note: solver timeout is passed via --timeout=<seconds> (see src/main.cpp).
    # We intentionally do not use Python's subprocess timeout so the solver controls timing.
    result = subprocess.run(args, capture_output=True, text=True)
    duration = time.perf_counter() - start

    console_path.write_text(result.stdout + "\n" + result.stderr)

    params = {
        "beta": beta,
        "nb_input_bits": nb_input_bits,
        "layout": layout,
        "targets": targets,
        "cost": cost,
        "heuristic": heuristic,
        "timeout_minutes": timeout_minutes,
        "timeout_seconds": timeout_seconds,
        "cmd": args,
        "returncode": result.returncode,
        "duration_seconds": duration,
    }
    params_path.write_text(json.dumps(params, indent=2))


def aggregate_results(base_out_dir: Path) -> None:
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

    summary = {
        "runs": runtime_count,
        "avg_runtime_seconds": (runtime_sum / runtime_count) if runtime_count else None,
        "costs": avg_costs,
    }
    (base_out_dir / "summary.json").write_text(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description="Run rmcmSolver benchmarks.")
    parser.add_argument("--dataset-root", default="benches/dataSets/2ADD", type=Path)
    parser.add_argument("--results-root", default="results", type=Path)
    parser.add_argument("--solver", default=Path("cmake-build-release/rmcmSolver"), type=Path)
    parser.add_argument("--set", choices=["timeout", "no_timeout", "all"], default="all",
                        help="Choose which dataset subset to run.")
    parser.add_argument("--cost", default="mux_bits", type=parse_cost_model,
                        help="Cost model: area, mux_count, mux_bits, luts, fpga_delay, asic_delay.")
    parser.add_argument("--max-cases", type=int, default=100,
                        help="Maximum number of test cases per file (<=100).")
    parser.add_argument("--heuristic", type=int, default=None,
                        help="Optional heuristic limit passed to solver.")
    parser.add_argument("--timeout-minutes", type=int, default=None,
                        help="Per-target-set timeout in minutes.")
    parser.add_argument("--nb-input-bits", type=int, default=8,
                        help="nb-input-bits passed to solver (default 8).")
    parser.add_argument("--layout", default="1,1",
                        help="Layout to pass to solver (default 1,1).")
    args = parser.parse_args()

    subsets = ["timeout", "no_timeout"] if args.set == "all" else [args.set]
    solver_path = args.solver.resolve()
    if not solver_path.exists():
        raise FileNotFoundError(f"Solver binary not found at {solver_path}")

    for subset in subsets:
        subset_root = args.dataset_root / subset
        if not subset_root.exists():
            print(f"Skipping missing subset {subset_root}")
            continue

        for bits_dir in sorted(subset_root.iterdir()):
            if not bits_dir.is_dir():
                continue
            beta = extract_beta_from_folder(bits_dir)
            for target_file in sorted(bits_dir.glob("*.txt")):
                targets_list = load_targets(target_file, args.max_cases)
                rel_base = target_file.relative_to(args.dataset_root)
                base_out_dir = args.results_root / rel_base.with_suffix("")

                for idx, targets in enumerate(targets_list, start=1):
                    case_dir = base_out_dir / f"case_{idx:03d}"
                    run_case(
                        solver_path=solver_path,
                        out_dir=case_dir,
                        beta=beta,
                        nb_input_bits=args.nb_input_bits,
                        layout=args.layout,
                        targets=targets,
                        cost=args.cost,
                        heuristic=args.heuristic,
                        timeout_minutes=args.timeout_minutes,
                    )
                aggregate_results(base_out_dir)


if __name__ == "__main__":
    main()
