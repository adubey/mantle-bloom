"""Plate speed distribution from `plate_torque_breakdown.py` outputs: by world, crust kind, plate and time (#942).

Usage (repo root):
    backend/.venv/bin/python bin/debug/summarize_speed_census.py ~/mantle-bloom-captures/issue942/census \
        [--out analysis/issue942/speed-census.json]

Speeds are area-weighted means of `|omega x r|` in cm/yr over every plate-sample. Reports, per world and pooled:
crust speed (all, continental, oceanic) and the speed mantle flow alone would give; quantiles of plate speed over
plate area; speed in four time windows; speed of plates with and without a subducting boundary; and the rank
correlation across worlds between solved speed and mantle-flow-only speed. Earth, for comparison: a mean of about
4–5 cm/yr and a range of about 1–10 cm/yr.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

EARTH_MEAN_CM_PER_YR = (4.0, 5.0)
WINDOWS_MYR = ((0, 50), (50, 100), (100, 150), (150, 200))


def _mean(pairs: list[tuple[float, float]]) -> float | None:
    """Weighted mean of (value, weight) pairs."""
    weight = sum(w for _, w in pairs)
    return sum(v * w for v, w in pairs) / weight if weight else None


def _quantiles(pairs: list[tuple[float, float]], qs=(0.1, 0.5, 0.9)) -> dict[str, float]:
    values = np.array([v for v, _ in pairs])
    weights = np.array([w for _, w in pairs])
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order]) / weights.sum()
    return {f"p{round(q * 100)}": float(values[order][np.searchsorted(cumulative, q)]) for q in qs}


def summarize(samples: list[dict]) -> dict:
    rows = [(sample["myr"], row) for sample in samples for row in sample["plates"]]
    plates = [(row["speed_cm_per_yr"]["all"], row["area_km2"]) for _, row in rows]
    kinds = {}
    for kind in ("continental", "oceanic"):
        kinds[kind] = _mean([
            (row["speed_cm_per_yr"][kind], row["continental_area_km2"] if kind == "continental" else row["area_km2"] - row["continental_area_km2"])
            for _, row in rows if row["speed_cm_per_yr"][kind] is not None
        ])
    area = sum(w for _, w in plates)
    subducting = [(row["speed_cm_per_yr"]["all"], row["area_km2"]) for _, row in rows if row.get("slab")]
    quiet = [(row["speed_cm_per_yr"]["all"], row["area_km2"]) for _, row in rows if not row.get("slab")]
    return {
        "plate_samples": len(rows),
        "plates_per_sample": len(rows) / len(samples),
        "speed_cm_per_yr": {"all": _mean(plates), **kinds},
        "mantle_flow_only_cm_per_yr": _mean([(row["counterfactual_speed_cm_per_yr"]["mantle_flow_only"]["all"], row["area_km2"]) for _, row in rows]),
        "plate_speed_quantiles_by_area": _quantiles(plates),
        "area_share_under_1_cm_per_yr": sum(w for v, w in plates if v < 1.0) / area,
        "area_share_over_4_cm_per_yr": sum(w for v, w in plates if v >= 4.0) / area,
        "fastest_plate_cm_per_yr": max(v for v, _ in plates),
        "by_window_myr": {
            f"{lo}-{hi}": _mean([(row["speed_cm_per_yr"]["all"], row["area_km2"]) for myr, row in rows if lo < myr <= hi]) for lo, hi in WINDOWS_MYR
        },
        "with_subducting_boundary": {"area_share": sum(w for _, w in subducting) / area, "speed_cm_per_yr": _mean(subducting)},
        "without_subducting_boundary": {"speed_cm_per_yr": _mean(quiet)},
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("runs", type=Path)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    worlds = {}
    pooled: list[dict] = []
    for path in sorted(args.runs.glob("*.json")):
        run = json.loads(path.read_text())
        worlds[str(run["seed"])] = {
            "solver_converged_all": all(s["solver"]["solver_converged"] for s in run["samples"]),
            "world_knobs": run["world_knobs"],
            **summarize(run["samples"]),
        }
        pooled += run["samples"]
    result = {"earth_mean_cm_per_yr": EARTH_MEAN_CM_PER_YR, "worlds": worlds, "pooled": summarize(pooled)}
    speed = result["pooled"]["speed_cm_per_yr"]["all"]
    result["pooled"]["earth_mean_over_model"] = [e / speed for e in EARTH_MEAN_CM_PER_YR]
    solved = [w["speed_cm_per_yr"]["all"] for w in worlds.values()]
    flow = [w["mantle_flow_only_cm_per_yr"] for w in worlds.values()]
    result["across_worlds"] = {
        "worlds": len(worlds),
        "speed_min_median_max": [float(np.min(solved)), float(np.median(solved)), float(np.max(solved))],
        "mantle_flow_only_min_median_max": [float(np.min(flow)), float(np.median(flow)), float(np.max(flow))],
        "spearman_solved_vs_mantle_flow_only": float(spearmanr(solved, flow).statistic) if len(worlds) > 2 else None,
    }
    text = json.dumps(result, indent=1) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
