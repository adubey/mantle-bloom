#!/usr/bin/env python3
"""Aggregates bin/debug/run_sweep.py's raw per-(parameter, multiplier, seed, checkpoint) JSONL
output into a per-(parameter, multiplier, checkpoint) mean/std summary across the SWEEP_SEEDS,
folding in the shared BASELINE_PARAM rows as every real parameter's own multiplier=1.0 point
(see sweep_lib.build_jobs / BASELINE_PARAM for why baseline is stored once, not per parameter).
Writes a single JSON file meant to be embedded directly into the dashboard artifact -- small
enough (8 params x 5 multipliers x 4 checkpoints x a handful of stats) that there's no need for
a live data source.

Usage: backend/.venv/bin/python bin/debug/aggregate_sweep.py [--in results/sweep_results.jsonl] [--out results/sweep_summary.json]
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sweep_lib import BASELINE_PARAM, PARAM_SPECS  # noqa: E402

THIS_DIR = Path(__file__).resolve().parent
DEFAULT_IN = THIS_DIR / "results" / "sweep_results.jsonl"
DEFAULT_OUT = THIS_DIR / "results" / "sweep_summary.json"

STAT_KEYS = ("land_fraction", "ice_cap_fraction", "plains_fraction", "land_volume_above_sea_km3")


def _mean_std(values: list[float]) -> tuple[float, float]:
    values = [v for v in values if v is not None]
    if not values:
        return None, None
    mean = statistics.fmean(values)
    std = statistics.pstdev(values) if len(values) > 1 else 0.0
    return mean, std


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--in", dest="inp", type=Path, default=DEFAULT_IN)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    if not args.inp.exists():
        parser.error(f"{args.inp} does not exist -- run run_sweep.py first")

    rows = []
    with args.inp.open() as f:
        for line in f:
            line = line.strip()
            if line:
                rows.append(json.loads(line))

    # Rows are keyed by literal parameter name in the file, but a BASELINE_PARAM row belongs to
    # every real parameter's own multiplier=1.0 point -- fan it out here so grouping below is a
    # simple per-parameter groupby with no special-casing.
    baseline_rows = [r for r in rows if r["parameter"] == BASELINE_PARAM]
    fanned_out = []
    for real_param in PARAM_SPECS:
        for row in baseline_rows:
            fanned_out.append({**row, "parameter": real_param})
    all_rows = [r for r in rows if r["parameter"] != BASELINE_PARAM] + fanned_out

    grouped: dict[tuple[str, float, int], dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    absolute_values: dict[tuple[str, float], set[float]] = defaultdict(set)
    for row in all_rows:
        key = (row["parameter"], row["multiplier"], row["checkpoint_years"])
        for stat in STAT_KEYS:
            grouped[key][stat].append(row[stat])
        if row["value"] is not None:
            absolute_values[(row["parameter"], row["multiplier"])].add(row["value"])

    summary = {"parameters": {}, "checkpoint_years": sorted({r["checkpoint_years"] for r in all_rows})}
    for param, spec in PARAM_SPECS.items():
        points = []
        multipliers = sorted({m for (p, m, _cp) in grouped if p == param})
        for multiplier in multipliers:
            checkpoints = sorted(cp for (p, m, cp) in grouped if p == param and m == multiplier)
            abs_values = absolute_values.get((param, multiplier))
            point = {
                "multiplier": multiplier,
                "value": (sorted(abs_values)[0] if abs_values else spec.baseline * multiplier),
                "n_seeds": len(grouped[(param, multiplier, checkpoints[0])][STAT_KEYS[0]]) if checkpoints else 0,
                "by_checkpoint": {},
            }
            for cp in checkpoints:
                stats_at_cp = grouped[(param, multiplier, cp)]
                point["by_checkpoint"][str(cp)] = {
                    stat: dict(zip(("mean", "std"), _mean_std(stats_at_cp[stat]))) for stat in STAT_KEYS
                }
            points.append(point)
        summary["parameters"][param] = {
            "label": spec.label,
            "unit": spec.unit,
            "baseline": spec.baseline,
            "points": points,
        }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2))
    n_complete = sum(
        1
        for param_data in summary["parameters"].values()
        for point in param_data["points"]
        if len(point["by_checkpoint"]) == len(summary["checkpoint_years"])
    )
    print(f"Wrote {args.out} ({n_complete} fully-covered parameter/multiplier points).")


if __name__ == "__main__":
    main()
