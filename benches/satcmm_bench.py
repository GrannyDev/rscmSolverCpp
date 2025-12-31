#!/usr/bin/env python3
"""
Benchmark runner for satcmm.

Runs satcmm for each target set in a dataset folder (e.g., timeout or
no_timeout), mirrors results alongside the dataset structure, captures
stdout/stderr, the last "Adder graph" block, feeds that block to an analysis
script, and aggregates the "Number of equivalent 2:1 MUXs" as cost.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import List, Optional


def load_target_sets(path: Path, max_cases: Optional[int]) -> List[str]:
    sets: List[str] = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if line:
                sets.append(line.replace(";", "|"))
                if max_cases is not None and len(sets) >= max_cases:
                    break
    return sets


def extract_adder_graph(output: str) -> Optional[str]:
    marker = "Adder graph:"
    idx = output.rfind(marker)
    if idx == -1:
        return None
    return output[idx + len(marker):].strip()


def run_analyzer(analyzer: Path, adder_graph: str, x_bits: int, const_count: int) -> tuple[Optional[int], str, float]:
    # If analyzer is a Python script, invoke via the current interpreter to avoid exec format issues.
    if analyzer.suffix == ".py":
        args = [sys.executable, str(analyzer), adder_graph, str(x_bits), str(const_count)]
    else:
        args = [str(analyzer), adder_graph, str(x_bits), str(const_count)]
    start = time.perf_counter()
    result = subprocess.run(args, capture_output=True, text=True)
    duration = time.perf_counter() - start
    output = result.stdout + "\n" + result.stderr
    mux_cost = None
    m = re.search(r"Number of equivalent 2:1 MUXs\s*=\s*(\d+)", output)
    if m:
        mux_cost = int(m.group(1))
    return mux_cost, output, duration


def extract_bits_from_folder(bits_folder: Path) -> int:
    digits = "".join(ch for ch in bits_folder.name if ch.isdigit())
    if not digits:
        raise ValueError(f"Cannot parse bit-width from folder name '{bits_folder.name}'")
    return int(digits)


def run_case(
    idx: int,
    total: int,
    targets: str,
    satcmm: Path,
    solver_binary: Path,
    analyzer: Path,
    out_root: Path,
    allow_negative: bool,
    min_adders: int,
    max_adders: int,
    timeout_seconds: Optional[int],
    sign_inversion: int,
    x_bits: int,
    const_count: Optional[int],
) -> None:
    case_dir = out_root / f"case_{idx:03d}"
    case_dir.mkdir(parents=True, exist_ok=True)
    console_path = case_dir / "console.txt"
    adder_graph_path = case_dir / "adder_graph.txt"
    analyzer_path = case_dir / "analyzer.txt"
    params_path = case_dir / "params.json"

    sat_args = [
        str(satcmm),
        targets,
        f"allow_negative_coefficients={1 if allow_negative else 0}",
        "solver_name=executable",
        f"executable_binary={solver_binary}",
        f"min_num_adders={min_adders}",
        f"max_num_adders={max_adders}",
        f"allow_coefficient_sign_inversion={sign_inversion}",
    ]
    if timeout_seconds is not None:
        sat_args.append(f"timeout={timeout_seconds}")

    start = time.perf_counter()
    sat_res = subprocess.run(sat_args, capture_output=True, text=True)
    sat_duration = time.perf_counter() - start
    sat_output = sat_res.stdout + "\n" + sat_res.stderr
    console_path.write_text(sat_output)

    adder_graph = extract_adder_graph(sat_output)
    mux_cost = None
    analyzer_output = ""
    analyzer_duration = 0.0
    if adder_graph:
        adder_graph_path.write_text(adder_graph)
        cc = const_count if const_count is not None else len([v for v in targets.split("|") if v])
        mux_cost, analyzer_output, analyzer_duration = run_analyzer(analyzer, adder_graph, x_bits, cc)
        analyzer_path.write_text(analyzer_output)

    params = {
        "targets": targets,
        "satcmm_args": sat_args,
        "returncode": sat_res.returncode,
        "satcmm_duration_seconds": sat_duration,
        "analyzer_duration_seconds": analyzer_duration if adder_graph else None,
        "mux_cost": mux_cost,
        "adder_graph_found": adder_graph is not None,
    }
    params_path.write_text(json.dumps(params, indent=2))

    progress = (idx / total) * 100.0
    print(f"\rProgress: {idx}/{total} ({progress:5.1f}%)", end="", flush=True)


def aggregate_results(out_root: Path) -> None:
    mux_sum = 0
    mux_count = 0
    sat_time = 0.0
    sat_runs = 0
    analyzer_time = 0.0
    analyzer_runs = 0

    for case_dir in sorted(out_root.glob("case_*")):
        params_file = case_dir / "params.json"
        if not params_file.exists():
            continue
        try:
            params = json.loads(params_file.read_text())
        except Exception:
            continue
        if params.get("mux_cost") is not None:
            mux_sum += params["mux_cost"]
            mux_count += 1
        sat_time += float(params.get("satcmm_duration_seconds", 0.0))
        sat_runs += 1
        if params.get("analyzer_duration_seconds") is not None:
            analyzer_time += float(params["analyzer_duration_seconds"])
            analyzer_runs += 1

    summary = {
        "runs": sat_runs,
        "avg_satcmm_runtime_seconds": (sat_time / sat_runs) if sat_runs else None,
        "avg_analyzer_runtime_seconds": (analyzer_time / analyzer_runs) if analyzer_runs else None,
        "avg_mux_cost": (mux_sum / mux_count) if mux_count else None,
    }
    (out_root / "summary.json").write_text(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark satcmm over a dataset folder.")
    parser.add_argument("--satcmm", required=True, type=Path, help="Path to satcmm executable.")
    parser.add_argument("--solver-binary", required=True, type=Path, help="Path to SAT solver binary (executable_binary).")
    parser.add_argument("--analyzer", required=True, type=Path, help="Path to analyze_result.py (or equivalent).")
    parser.add_argument("--dataset-root", default=Path("dataSets/2ADD"), type=Path, help="Dataset root (expects timeout/no_timeout subfolders).")
    parser.add_argument("--set", choices=["timeout", "no_timeout", "all"], default="all",
                        help="Choose which dataset subset to run.")
    parser.add_argument("--results-root", required=True, type=Path, help="Where to write results.")
    parser.add_argument("--max-cases", type=int, default=None, help="Maximum number of target sets per file.")
    parser.add_argument("--allow-negative", action="store_true", help="Set allow_negative_coefficients=1 (default off).")
    parser.add_argument("--min-adders", type=int, default=1)
    parser.add_argument("--max-adders", type=int, default=2)
    parser.add_argument("--x-bits", type=int, default=8, help="Input bit-width passed to analyzer.")
    parser.add_argument("--const-count", type=int, default=None, help="Optional constant count override for analyzer.")
    parser.add_argument("--timeout", type=int, default=None, help="Optional timeout (seconds) passed to satcmm.")
    parser.add_argument(
        "--sign-inversion",
        type=int,
        choices=[0, 1, 2],
        default=0,
        help="allow_coefficient_sign_inversion: 0=normalized, 1=SAT may invert, 2=exact signs.",
    )
    args = parser.parse_args()

    subsets = ["timeout", "no_timeout"] if args.set == "all" else [args.set]

    files_data: list[tuple[Path, int, List[str]]] = []
    total_cases = 0
    for subset in subsets:
        subset_root = args.dataset_root / subset
        if not subset_root.exists():
            print(f"Skipping missing subset {subset_root}")
            continue
        for bits_dir in sorted(subset_root.iterdir()):
            if not bits_dir.is_dir():
                continue
            try:
                x_bits = extract_bits_from_folder(bits_dir)
            except ValueError as e:
                print(f"Skipping folder '{bits_dir}': {e}")
                continue
            for target_file in sorted(bits_dir.glob("*.txt")):
                targets_sets = load_target_sets(target_file, args.max_cases)
                if not targets_sets:
                    continue
                rel_base = target_file.relative_to(args.dataset_root)
                base_out_dir = args.results_root / rel_base.with_suffix("")
                files_data.append((base_out_dir, x_bits, targets_sets))
                total_cases += len(targets_sets)

    if total_cases == 0:
        print("No target sets found.")
        return

    args.results_root.mkdir(parents=True, exist_ok=True)
    case_idx = 0
    for base_out_dir, x_bits, targets_sets in files_data:
        for targets in targets_sets:
            case_idx += 1
            run_case(
                idx=case_idx,
                total=total_cases,
                targets=targets,
                satcmm=args.satcmm.resolve(),
                solver_binary=args.solver_binary.resolve(),
                analyzer=args.analyzer.resolve(),
                out_root=base_out_dir,
            allow_negative=args.allow_negative,
            min_adders=args.min_adders,
            max_adders=args.max_adders,
            timeout_seconds=args.timeout,
            sign_inversion=args.sign_inversion,
            x_bits=x_bits,
            const_count=args.const_count,
        )
        aggregate_results(base_out_dir)
    print()  # newline after progress


if __name__ == "__main__":
    main()
