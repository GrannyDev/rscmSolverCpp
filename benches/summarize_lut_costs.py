#!/usr/bin/env python3
"""
Summarize LUT estimation error from cost.json files.
"""

from __future__ import annotations

import argparse
import json
import math
import statistics
from pathlib import Path
from typing import Iterable, Optional


def iter_cost_files(root: Path) -> Iterable[Path]:
    return root.rglob("cost.json")


def load_estimated_luts(solution_path: Path) -> Optional[float]:
    try:
        solution = json.loads(solution_path.read_text())
    except Exception:
        return None
    costs = solution.get("meta", {}).get("costs", {})
    est = costs.get("luts")
    if not isinstance(est, (int, float)) or not math.isfinite(est):
        return None
    return float(est)


def extract_features(solution: dict) -> dict[str, float]:
    meta = solution.get("meta", {})
    muxes_bw = meta.get("muxes_bw", [])
    if not isinstance(muxes_bw, list):
        muxes_bw = []

    features: dict[str, float] = {
        "w_in": float(meta.get("wIn", 0)),
        "w_conf": float(meta.get("wConf", 0)),
        "w_selec": float(meta.get("wSelec", 0)),
        "w_out": float(meta.get("wOut", 0)),
        "needs_table": 1.0 if meta.get("needs_table") else 0.0,
        "mux_count": float(len(muxes_bw)),
        "mux_bw_sum": float(sum(muxes_bw)) if muxes_bw else 0.0,
        "mux_bw_max": float(max(muxes_bw)) if muxes_bw else 0.0,
    }
    features["mux_bw_avg"] = (
        features["mux_bw_sum"] / features["mux_count"] if features["mux_count"] else 0.0
    )

    adders = [k for k in solution.keys() if k.startswith("ADDER_")]
    features["adder_count"] = float(len(adders))

    def init_type_metrics(prefix: str) -> None:
        features[f"{prefix}_bw_sum"] = 0.0
        features[f"{prefix}_mux_count"] = 0.0
        features[f"{prefix}_mux_bw_sum"] = 0.0
        features[f"{prefix}_options_sum"] = 0.0

    for key in ("left_inputs", "right_inputs", "left_shifts", "right_shifts", "outputs_shifts", "add_sub"):
        init_type_metrics(key)

    sum_max_lr_shift_bw = 0.0
    for adder_key in adders:
        adder = solution.get(adder_key, {})

        def handle_param(param_name: str, feature_prefix: str) -> tuple[float, float]:
            param = adder.get(param_name, {})
            bw = float(param.get("computed_bitwidth", 0))
            features[f"{feature_prefix}_bw_sum"] += bw

            options = 0
            if "inputs" in param and isinstance(param["inputs"], list):
                options = len(param["inputs"])
            elif "shifts" in param and isinstance(param["shifts"], list):
                options = len(param["shifts"])
            elif param_name == "ADD_SUB":
                options = 2 if param.get("type") == "both" else 1
            features[f"{feature_prefix}_options_sum"] += float(options)

            mux_nb = param.get("mux_nb")
            if isinstance(mux_nb, int) and 0 <= mux_nb < len(muxes_bw):
                features[f"{feature_prefix}_mux_count"] += 1.0
                features[f"{feature_prefix}_mux_bw_sum"] += float(muxes_bw[mux_nb])

            return bw, float(options)

        left_bw, _ = handle_param("LEFT_SHIFTS", "left_shifts")
        right_bw, _ = handle_param("RIGHT_SHIFTS", "right_shifts")
        sum_max_lr_shift_bw += max(left_bw, right_bw)

        handle_param("LEFT_INPUTS", "left_inputs")
        handle_param("RIGHT_INPUTS", "right_inputs")
        handle_param("OUTPUTS_SHIFTS", "outputs_shifts")
        handle_param("ADD_SUB", "add_sub")

    features["sum_max_lr_shift_bw"] = sum_max_lr_shift_bw
    return features


def solve_linear_system(a: list[list[float]], b: list[float]) -> list[float]:
    n = len(b)
    for i in range(n):
        pivot = max(range(i, n), key=lambda r: abs(a[r][i]))
        if abs(a[pivot][i]) < 1e-12:
            raise ValueError("Singular matrix in regression")
        if pivot != i:
            a[i], a[pivot] = a[pivot], a[i]
            b[i], b[pivot] = b[pivot], b[i]
        for r in range(i + 1, n):
            factor = a[r][i] / a[i][i]
            if factor == 0:
                continue
            for c in range(i, n):
                a[r][c] -= factor * a[i][c]
            b[r] -= factor * b[i]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        acc = b[i]
        for c in range(i + 1, n):
            acc -= a[i][c] * x[c]
        x[i] = acc / a[i][i]
    return x


def linear_regression(rows: list[dict[str, float]], features: list[str], target_key: str) -> dict[str, object]:
    x_rows = []
    y_vals = []
    for row in rows:
        if target_key not in row:
            continue
        try:
            x = [1.0] + [float(row[f]) for f in features]
        except KeyError:
            continue
        x_rows.append(x)
        y_vals.append(float(row[target_key]))

    if not x_rows:
        raise ValueError("No rows available for regression")

    cols = len(x_rows[0])
    xtx = [[0.0] * cols for _ in range(cols)]
    xty = [0.0] * cols
    for x, y in zip(x_rows, y_vals):
        for i in range(cols):
            xty[i] += x[i] * y
            for j in range(cols):
                xtx[i][j] += x[i] * x[j]

    coeffs = solve_linear_system(xtx, xty)
    preds = [sum(c * v for c, v in zip(coeffs, x)) for x in x_rows]

    mean_y = statistics.mean(y_vals)
    ss_tot = sum((y - mean_y) ** 2 for y in y_vals)
    ss_res = sum((y - yhat) ** 2 for y, yhat in zip(y_vals, preds))
    r2 = 1.0 - (ss_res / ss_tot) if ss_tot else 0.0
    mae = statistics.mean(abs(y - yhat) for y, yhat in zip(y_vals, preds))
    rmse = math.sqrt(statistics.mean((y - yhat) ** 2 for y, yhat in zip(y_vals, preds)))

    return {
        "coeffs": coeffs,
        "features": ["intercept"] + features,
        "r2": r2,
        "mae": mae,
        "rmse": rmse,
        "rows": len(y_vals),
    }


def extract_bits_folder(path: Path) -> Optional[str]:
    for part in path.parts:
        if part.endswith("bits") and part[:-4].isdigit():
            return part
    return None


def compute_stats(errors: list[float], ratios: list[float]) -> dict[str, float]:
    abs_errors = [abs(e) for e in errors]
    stats = {
        "mean_err": statistics.mean(errors),
        "median_err": statistics.median(errors),
        "mean_abs_err": statistics.mean(abs_errors),
        "median_abs_err": statistics.median(abs_errors),
        "min_err": min(errors),
        "max_err": max(errors),
    }
    if ratios:
        stats["mean_ratio"] = statistics.mean(ratios)
        stats["median_ratio"] = statistics.median(ratios)
    return stats


def main() -> None:
    parser = argparse.ArgumentParser(description="Summarize LUT estimation errors.")
    parser.add_argument(
        "--root",
        type=Path,
        default=Path("benches/bench_new_luts_recomputed_test"),
        help="Root folder containing cost.json files.",
    )
    parser.add_argument(
        "--per-bits",
        action="store_true",
        help="Print a per-<N>bits breakdown.",
    )
    parser.add_argument(
        "--emit-csv",
        type=Path,
        help="Write per-case features to a CSV file.",
    )
    parser.add_argument(
        "--emit-json",
        type=Path,
        help="Write per-case features to a JSON file.",
    )
    parser.add_argument(
        "--list-features",
        action="store_true",
        help="List available regression features and exit.",
    )
    parser.add_argument(
        "--regress",
        action="store_true",
        help="Run a linear regression to fit synthesized LUTs from features.",
    )
    parser.add_argument(
        "--regress-features",
        default="mux_bw_sum,sum_max_lr_shift_bw,outputs_shifts_bw_sum,add_sub_bw_sum",
        help="Comma-separated list of feature names for regression.",
    )
    args = parser.parse_args()

    total_files = 0
    valid = 0
    missing_solution = 0
    errors: list[float] = []
    ratios: list[float] = []
    percent_errors: list[float] = []
    by_bits: dict[str, list[float]] = {}
    rows: list[dict[str, float]] = []

    for path in iter_cost_files(args.root):
        total_files += 1
        try:
            data = json.loads(path.read_text())
        except Exception:
            continue
        solution_path = path.with_name("solution.json")
        if not solution_path.exists():
            missing_solution += 1
            continue
        try:
            solution = json.loads(solution_path.read_text())
        except Exception:
            continue
        costs = solution.get("meta", {}).get("costs", {})
        est = costs.get("luts")
        synth = data.get("synthetised_luts")
        if not isinstance(est, (int, float)) or not math.isfinite(est):
            continue
        if not isinstance(synth, (int, float)):
            continue
        if not math.isfinite(synth):
            continue

        features = extract_features(solution)

        valid += 1
        error = float(est) - float(synth)
        errors.append(error)
        if synth != 0:
            ratios.append(float(est) / float(synth))
            percent_errors.append((float(est) - float(synth)) / float(synth) * 100.0)

        row = dict(features)
        row["estimated_luts"] = float(est)
        row["synthetised_luts"] = float(synth)
        row["error_luts"] = float(est) - float(synth)
        row["error_percent"] = (
            (float(est) - float(synth)) / float(synth) * 100.0 if synth != 0 else 0.0
        )
        rows.append(row)

        if args.per_bits:
            bits = extract_bits_folder(path)
            if bits is not None:
                by_bits.setdefault(bits, []).append(error)

    if not errors:
        print(f"No valid entries found under {args.root}")
        return

    if args.list_features:
        if rows:
            print("features:")
            for name in sorted(rows[0].keys()):
                if name in {"estimated_luts", "synthetised_luts", "error_luts", "error_percent"}:
                    continue
                print(f"- {name}")
        else:
            print("No rows to list features from.")
        return

    stats = compute_stats(errors, ratios)
    print(f"root: {args.root}")
    print(f"cost.json files: {total_files}")
    print(f"valid entries: {valid}")
    print(f"missing solution.json: {missing_solution}")
    print(f"mean err (est - synth): {stats['mean_err']:.6g}")
    print(f"median err: {stats['median_err']:.6g}")
    print(f"mean abs err: {stats['mean_abs_err']:.6g}")
    print(f"median abs err: {stats['median_abs_err']:.6g}")
    print(f"min err: {stats['min_err']:.6g}")
    print(f"max err: {stats['max_err']:.6g}")
    if "mean_ratio" in stats:
        print(f"mean ratio: {stats['mean_ratio']:.6g}")
        print(f"median ratio: {stats['median_ratio']:.6g}")
    if percent_errors:
        abs_percent = [abs(p) for p in percent_errors]
        print(f"mean % err: {statistics.mean(percent_errors):.6g}%")
        print(f"median % err: {statistics.median(percent_errors):.6g}%")
        print(f"mean abs % err: {statistics.mean(abs_percent):.6g}%")
        print(f"median abs % err: {statistics.median(abs_percent):.6g}%")

    if args.per_bits:
        for bits in sorted(by_bits):
            per = by_bits[bits]
            if not per:
                continue
            print(
                f"{bits}: count={len(per)} "
                f"mean_err={statistics.mean(per):.6g} "
                f"median_err={statistics.median(per):.6g}"
            )

    if args.emit_csv and rows:
        feature_names = sorted(rows[0].keys())
        header = ",".join(feature_names)
        lines = [header]
        for row in rows:
            values = [str(row.get(name, "")) for name in feature_names]
            lines.append(",".join(values))
        args.emit_csv.write_text("\n".join(lines))
        print(f"Wrote CSV: {args.emit_csv}")

    if args.emit_json and rows:
        args.emit_json.write_text(json.dumps(rows, indent=2))
        print(f"Wrote JSON: {args.emit_json}")

    if args.regress:
        feature_list = [f.strip() for f in args.regress_features.split(",") if f.strip()]
        try:
            reg = linear_regression(rows, feature_list, "synthetised_luts")
        except ValueError as exc:
            print(f"Regression failed: {exc}")
            return
        print("regression:")
        for name, coeff in zip(reg["features"], reg["coeffs"]):
            print(f"  {name} = {coeff:.6g}")
        print(f"  r2 = {reg['r2']:.6g}")
        print(f"  mae = {reg['mae']:.6g}")
        print(f"  rmse = {reg['rmse']:.6g}")
        print(f"  rows = {reg['rows']}")


if __name__ == "__main__":
    main()
