#!/usr/bin/env python3
"""
Benchmark runner for kmult over the datasets in benches/dataSets/2ADD/.

It mirrors the dataset hierarchy into a results folder and, for each target
set, creates a folder per test case containing:
  - console.txt : kmult stdout/stderr
  - params.json : parameters passed to kmult plus timing and parsed metrics
  - module.dot : dot file emitted by kmult (if generated)

Usage examples:
    python3 benches/pagfuse.py --set timeout --kmult /path/to/kmult
    python3 benches/pagfuse.py --results-root benches/results/pagfuse --kmult-args "-r100000 -f0"
"""

from __future__ import annotations

import argparse
import json
import re
import shlex
import shutil
import subprocess
import time
from pathlib import Path
from typing import List, Optional


def parse_dot_file(dot_path: Path) -> tuple[int, int]:
    node_labels: dict[str, str] = {}
    edges: list[tuple[str, str]] = []
    node_re = re.compile(r'(node\d+)\s*\[.*label\s*=\s*"([^"]+)"')
    edge_re = re.compile(r'(node\d+)\s*->\s*(node\d+)')

    with dot_path.open() as f:
        for line in f:
            line = line.strip()
            node_match = node_re.match(line)
            if node_match:
                node_labels[node_match.group(1)] = node_match.group(2)
            edge_match = edge_re.match(line)
            if edge_match:
                edges.append((edge_match.group(1), edge_match.group(2)))

    adder_count = sum(1 for label in node_labels.values() if label in {"Add", "Sub", "AddSub"})

    mux_entries: dict[str, int] = {}
    for _, dst in edges:
        if node_labels.get(dst) == "Mux":
            mux_entries[dst] = mux_entries.get(dst, 0) + 1
    mux_incoming_sum = sum(entries - 1 for entries in mux_entries.values())

    return adder_count, mux_incoming_sum


def extract_bits_from_folder(bits_folder: Path) -> int:
    digits = "".join(ch for ch in bits_folder.name if ch.isdigit())
    if not digits:
        raise ValueError(f"Cannot parse bit-width from folder name '{bits_folder.name}'")
    return int(digits)


def extract_bits_from_filename(target_file: Path) -> Optional[int]:
    match = re.search(r"on(\d+)bits", target_file.stem)
    if not match:
        return None
    return int(match.group(1))


def load_targets(file_path: Path, max_cases: Optional[int]) -> List[str]:
    targets: List[str] = []
    with file_path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            targets.append(line)
            if max_cases is not None and len(targets) >= max_cases:
                break
    return targets


def parse_targets(line: str) -> List[str]:
    return [s.strip() for s in line.split(";") if s.strip()]


def has_c_arg(kmult_args: List[str]) -> bool:
    for arg in kmult_args:
        if arg == "-c":
            return True
        if arg.startswith("-c") and len(arg) > 2:
            return True
    return False


def run_case(
    kmult_path: Path,
    out_dir: Path,
    dot_file: Path,
    targets_line: str,
    kmult_args: List[str],
    c_bits: int,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    console_path = out_dir / "console.txt"
    params_path = out_dir / "params.json"

    targets = parse_targets(targets_line)
    args = [str(kmult_path)] + kmult_args[:]
    if not has_c_arg(args):
        args.append(f"-c{c_bits}")
    args += targets

    start = time.perf_counter()
    result = subprocess.run(args, capture_output=True, text=True, cwd=out_dir)
    duration = time.perf_counter() - start
    output = (result.stdout or "") + "\n" + (result.stderr or "")
    console_path.write_text(output)

    dot_path = dot_file if dot_file.is_absolute() else out_dir / dot_file
    adder_count = None
    mux_entries_sum = None
    parse_error = None
    if dot_path.exists():
        try:
            adder_count, mux_entries_sum = parse_dot_file(dot_path)
        except Exception as exc:
            parse_error = str(exc)
    else:
        parse_error = f"Missing dot file: {dot_path}"

    params = {
        "targets": targets_line,
        "kmult_args": args,
        "returncode": result.returncode,
        "duration_seconds": duration,
        "adder_count": adder_count,
        "mux_entries_sum": mux_entries_sum,
        "dot_file": str(dot_path) if dot_path.exists() else None,
        "dot_parse_error": parse_error,
    }
    params_path.write_text(json.dumps(params, indent=2))


def aggregate_results(base_out_dir: Path) -> None:
    runs = 0
    runtime_sum = 0.0
    adder_sum = 0
    adder_count = 0
    mux_sum = 0
    mux_count = 0
    parse_errors = 0

    for case_dir in sorted(base_out_dir.glob("case_*")):
        params_path = case_dir / "params.json"
        if not params_path.exists():
            continue
        try:
            params = json.loads(params_path.read_text())
        except Exception:
            continue
        runs += 1
        runtime_sum += float(params.get("duration_seconds", 0.0))
        if params.get("adder_count") is not None:
            adder_sum += int(params["adder_count"])
            adder_count += 1
        if params.get("mux_entries_sum") is not None:
            mux_sum += int(params["mux_entries_sum"])
            mux_count += 1
        if params.get("dot_parse_error"):
            parse_errors += 1

    summary = {
        "runs": runs,
        "avg_runtime_seconds": (runtime_sum / runs) if runs else None,
        "avg_adders": (adder_sum / adder_count) if adder_count else None,
        "avg_mux_entries_sum": (mux_sum / mux_count) if mux_count else None,
        "dot_parse_errors": parse_errors,
    }
    (base_out_dir / "summary.json").write_text(json.dumps(summary, indent=2))


def resolve_executable(path: Path) -> Path:
    if path.exists():
        return path
    resolved = shutil.which(str(path))
    if resolved:
        return Path(resolved)
    raise FileNotFoundError(f"Executable not found at {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run kmult over dataset target sets.")
    parser.add_argument("--dataset-root", default=Path("benches/dataSets/2ADD"), type=Path)
    parser.add_argument("--results-root", default=Path("results"), type=Path)
    parser.add_argument("--kmult", default=Path("/home/smith/Compiled/muxmcmtool/kmult"), type=Path)
    parser.add_argument("--kmult-args", default="-r100000 -f0",
                        help="Extra args passed to kmult (excluding -c and targets).")
    parser.add_argument("--dot-file", default="module.dot",
                        help="Dot file name emitted by kmult (default: module.dot).")
    parser.add_argument("--set", choices=["timeout", "no_timeout", "all"], default="all",
                        help="Choose which dataset subset to run.")
    parser.add_argument("--max-cases", type=int, default=100,
                        help="Maximum number of target sets per file (<=100).")
    parser.add_argument("--c-bits", type=int, default=None,
                        help="Override -c<value> passed to kmult; default derives from file name.")
    args = parser.parse_args()

    subsets = ["timeout", "no_timeout"] if args.set == "all" else [args.set]
    kmult_path = resolve_executable(args.kmult)
    kmult_args = shlex.split(args.kmult_args)
    dot_file = Path(args.dot_file)

    for subset in subsets:
        subset_root = args.dataset_root / subset
        if not subset_root.exists():
            print(f"Skipping missing subset {subset_root}")
            continue

        for bits_dir in sorted(subset_root.iterdir()):
            if not bits_dir.is_dir():
                continue
            try:
                bits = extract_bits_from_folder(bits_dir)
            except ValueError as exc:
                print(f"Skipping folder '{bits_dir}': {exc}")
                continue

            for target_file in sorted(bits_dir.glob("*.txt")):
                file_bits = extract_bits_from_filename(target_file)
                c_bits = args.c_bits if args.c_bits is not None else (file_bits or bits)
                targets_list = load_targets(target_file, args.max_cases)
                if not targets_list:
                    continue
                rel_base = target_file.relative_to(args.dataset_root)
                base_out_dir = args.results_root / rel_base.with_suffix("")
                total_cases = len(targets_list)

                for idx, targets_line in enumerate(targets_list, start=1):
                    case_dir = base_out_dir / f"case_{idx:03d}"
                    run_case(
                        kmult_path=kmult_path,
                        out_dir=case_dir,
                        dot_file=dot_file,
                        targets_line=targets_line,
                        kmult_args=kmult_args,
                        c_bits=c_bits,
                    )
                    progress = (idx / total_cases) * 100 if total_cases else 100.0
                    print(f"\r{target_file.name}: {idx}/{total_cases} ({progress:5.1f}%)", end="", flush=True)
                if total_cases:
                    print()
                aggregate_results(base_out_dir)


if __name__ == "__main__":
    main()
