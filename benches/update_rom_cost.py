#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path


def count_both_ops(solution: dict) -> int:
    count = 0
    for key, block in solution.items():
        if key == "meta":
            continue
        if not isinstance(block, dict):
            continue
        add_sub = block.get("ADD_SUB")
        if isinstance(add_sub, dict) and add_sub.get("type") == "both":
            count += 1
    return count


def compute_rom_cost(mux_bits: float, both_count: int, w_conf: int) -> float:
    return 442 + 1.35 * (mux_bits * (2 ** w_conf))


def update_solution(src_path: Path, dst_path: Path) -> None:
    solution = json.loads(src_path.read_text())
    meta = solution.get("meta", {})
    if not isinstance(meta, dict):
        raise ValueError(f"{src_path} missing meta object")

    costs = meta.get("costs", {})
    if not isinstance(costs, dict):
        raise ValueError(f"{src_path} missing meta.costs object")

    if "mux_bits" not in costs:
        raise ValueError(f"{src_path} missing meta.costs.mux_bits")
    if "area_cost" not in costs:
        raise ValueError(f"{src_path} missing meta.costs.area_cost")
    if "wConf" not in meta:
        raise ValueError(f"{src_path} missing meta.wConf")

    mux_bits = float(costs["mux_bits"])
    area_cost = float(costs["area_cost"])
    w_conf = int(meta["wConf"])

    both_count = count_both_ops(solution)
    rom_cost = compute_rom_cost(mux_bits, both_count, w_conf)

    costs["rom_cost"] = rom_cost
    costs["area_cost"] = area_cost - rom_cost

    dst_path.write_text(json.dumps(solution, indent=2) + "\n")


def mirror_and_update(src_root: Path, dst_root: Path) -> None:
    if not src_root.exists():
        raise FileNotFoundError(f"Input folder not found: {src_root}")
    if dst_root.exists() and any(dst_root.iterdir()):
        raise FileExistsError(f"Output folder already exists and is not empty: {dst_root}")

    for root, dirs, files in os.walk(src_root):
        rel_root = Path(root).relative_to(src_root)
        dst_dir = dst_root / rel_root
        dst_dir.mkdir(parents=True, exist_ok=True)

        for name in files:
            src_path = Path(root) / name
            dst_path = dst_dir / name
            if name == "solution.json":
                update_solution(src_path, dst_path)
            else:
                shutil.copy2(src_path, dst_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Clone a bench results folder and inject rom_cost into solution.json files.",
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("benches/bench_new_luts_recomputed"),
        help="Folder containing bench results (default: benches/bench_new_luts_recomputed).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("benches/bench_new_luts_recomputed_rom_cost"),
        help="Output folder to create (default: benches/bench_new_luts_recomputed_rom_cost).",
    )
    args = parser.parse_args()

    mirror_and_update(args.input, args.output)


if __name__ == "__main__":
    main()
