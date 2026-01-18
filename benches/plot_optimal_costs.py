#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("pgf")
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
        default=None,
        help="Base output image path (if set, writes <stem>_luts and <stem>_area images)",
    )
    parser.add_argument(
        "--out-luts",
        default="optimal_cost_histogram_luts.png",
        help="Output path for LUTs histogram (default: optimal_cost_histogram_luts.png)",
    )
    parser.add_argument(
        "--out-area",
        default="optimal_cost_histogram_area.png",
        help="Output path for area histogram (default: optimal_cost_histogram_area.png)",
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
    parser.add_argument(
        "--font-size",
        type=float,
        default=10.0,
        help="Base font size in points (default: 10)",
    )
    parser.add_argument(
        "--legend-font-size",
        type=float,
        default=None,
        help="Legend font size in points (defaults to --font-size)",
    )
    args = parser.parse_args()

    legend_font_size = args.legend_font_size if args.legend_font_size is not None else args.font_size
    matplotlib.rcParams.update(
        {
            "text.usetex": True,
            "pgf.rcfonts": False,
            "font.family": "serif",
            "font.size": args.font_size,
            "legend.fontsize": legend_font_size,
        }
    )

    fig_width = 3.2
    fig_height = 1.8

    folder = Path(args.dir)
    if not folder.exists():
        print(f"Directory not found: {folder}", file=sys.stderr)
        return 1

    luts, area, missing, lowest_luts, lowest_luts_path = load_costs(folder)
    if not luts and not area:
        print("No cost data found in JSON files.", file=sys.stderr)
        return 1

    out_luts = Path(args.out_luts).with_suffix(".pgf")
    out_area = Path(args.out_area).with_suffix(".pgf")
    if args.out:
        base = Path(args.out)
        stem = base.stem
        suffix = base.suffix or ".png"
        out_luts = base.with_name(f"{stem}_luts{suffix}").with_suffix(".pgf")
        out_area = base.with_name(f"{stem}_area{suffix}").with_suffix(".pgf")

    figs = []
    if luts:
        fig_luts, ax_luts = plt.subplots(1, 1, figsize=(fig_width, fig_height))
        ax_luts.hist(luts, bins=args.bins, color="#2E86AB", edgecolor="black")
        ax_luts.set_xlabel("LUTs estimate")
        ax_luts.set_ylabel("Count (log scale)")
        ax_luts.set_yscale("log")
        fig_luts.tight_layout()
        fig_luts.savefig(out_luts)
        print(f"Saved LUTs plot to {out_luts}")
        figs.append(fig_luts)
    else:
        print("No LUT cost data found; LUTs plot not created.", file=sys.stderr)

    if area:
        fig_area, ax_area = plt.subplots(1, 1, figsize=(fig_width, fig_height))
        ax_area.hist(area, bins=args.bins, color="#E27D60", edgecolor="black")
        ax_area.set_xlabel("ASIC area estimate")
        ax_area.set_ylabel("Count (log scale)")
        ax_area.set_yscale("log")
        fig_area.tight_layout()
        fig_area.savefig(out_area)
        print(f"Saved area plot to {out_area}")
        figs.append(fig_area)
    else:
        print("No area cost data found; area plot not created.", file=sys.stderr)
    if missing:
        print(f"Warning: {missing} files missing luts or area_cost values.", file=sys.stderr)
    if lowest_luts_path is not None:
        print(f"Lowest LUTs: {lowest_luts} in {lowest_luts_path.name}")
    if luts:
        print(f"LUTs min/max: {min(luts)} / {max(luts)}")
    if area:
        print(f"Area min/max: {min(area)} / {max(area)}")

    if args.show:
        if figs:
            plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
