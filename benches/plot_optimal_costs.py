#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt


def load_costs(folder: Path):
    luts = []
    area = []
    missing = 0
    lowest_luts = None
    lowest_luts_path = None
    for path in sorted(folder.glob("*.json")):
        try:
            with path.open("r", encoding="utf-8") as f:
                data = json.load(f)
        except (OSError, json.JSONDecodeError):
            print(f"Skipping unreadable json: {path}", file=sys.stderr)
            continue

        costs = data.get("meta", {}).get("costs", {})
        luts_val = costs.get("luts")
        area_val = costs.get("area_cost")
        if isinstance(luts_val, int):
            luts.append(luts_val)
            if lowest_luts is None or luts_val < lowest_luts:
                lowest_luts = luts_val
                lowest_luts_path = path
        if isinstance(area_val, int):
            area.append(area_val)
        if not isinstance(luts_val, int) or not isinstance(area_val, int):
            missing += 1
    return luts, area, missing, lowest_luts, lowest_luts_path


def main():
    parser = argparse.ArgumentParser(description="Plot histograms of LUT and area costs for optimal solutions.")
    parser.add_argument(
        "--dir",
        default="all_optimal_solution",
        help="Directory containing JSON solutions (default: all_optimal_solution)",
    )
    parser.add_argument(
        "--out",
        default="optimal_cost_histograms.png",
        help="Output image path (default: optimal_cost_histograms.png)",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Show the plot window instead of only saving",
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=20,
        help="Number of histogram bins (default: 20)",
    )
    args = parser.parse_args()

    folder = Path(args.dir)
    if not folder.exists():
        print(f"Directory not found: {folder}", file=sys.stderr)
        return 1

    luts, area, missing, lowest_luts, lowest_luts_path = load_costs(folder)
    if not luts and not area:
        print("No cost data found in JSON files.", file=sys.stderr)
        return 1

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    if luts:
        axes[0].hist(luts, bins=args.bins, color="#2E86AB", edgecolor="black")
    axes[0].set_title("LUT Cost Histogram")
    axes[0].set_xlabel("Estimated LUTs")
    axes[0].set_ylabel("Count")
    axes[0].set_yscale("log")

    if area:
        axes[1].hist(area, bins=args.bins, color="#E27D60", edgecolor="black")
    axes[1].set_title("Area Cost Histogram")
    axes[1].set_xlabel("Estimated Area Cost")
    axes[1].set_ylabel("Count")
    axes[1].set_yscale("log")

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"Saved plot to {args.out}")
    if missing:
        print(f"Warning: {missing} files missing luts or area_cost values.", file=sys.stderr)
    if lowest_luts_path is not None:
        print(f"Lowest LUTs: {lowest_luts} in {lowest_luts_path.name}")

    if args.show:
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
