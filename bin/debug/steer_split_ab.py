#!/usr/bin/env python3
"""A/B the RFC 0003 split-cut steering (#389) on generated worlds, one seed at a time.

Wraps rebaseline_rfcs._measure (default generated quad worlds, cratons on, 5 Myr steps,
100 Myr checkpoints) and adds read-only taps that work on both the baseline and the branch:
- every successful split's daughters' RMS flow-fit residual against their set Euler poles
  (RFC 0003 acceptance: the mean rises by less than 10%);
- time spent in maybe_split_plate, and per call in the steering (branch only), with how
  often a cut whose band held strength (steer_needed) was steered.
`--set NAME=VALUE` overrides a merge_split constant (e.g. SPLIT_CUT_BAND_MULT=5.0 for the
band-width sensitivity run). compare_cratons puts its own tree's backend first on sys.path,
so run the baseline from a copy of this script inside an origin/main worktree's bin/debug;
the output records which app it imported. Saves are scratch artifacts; only compact JSON belongs in analysis/.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/steer_split_ab.py \
        --out /tmp/branch.json --save-dir /tmp/issue389
    cp bin/debug/steer_split_ab.py ../mantle-bloom-origin-main/bin/debug/
    (cd ../mantle-bloom-origin-main && PYTHONPATH=backend backend/.venv/bin/python \
        bin/debug/steer_split_ab.py --out /tmp/main.json --save-dir /tmp/issue389-main)
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/steer_split_ab.py --seeds 456862899 \
        --years 3e9 --out /tmp/long.json --save-dir /tmp/issue389
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import rebaseline_rfcs
from app import mantle, merge_split


def _taps() -> tuple[dict, list]:
    stats = {"split_attempt_s": 0.0, "split_attempts": 0, "residuals": [],
             "steer_calls": 0, "steer_s": 0.0, "steer_s_max": 0.0, "steer_max_nodes": 0, "steer_needed": 0, "steered": 0}
    originals = [("maybe_split_plate", merge_split.maybe_split_plate)]

    census: dict[str, int] = {}
    pending: list[dict] = []
    stats["gate_census"] = census

    def maybe_split_plate(world, plate):
        pending.clear()
        cuts_before = _cut_count(world)
        started = time.perf_counter()
        result = originals[0][1](world, plate)
        if pending and _cut_count(world) > cuts_before:
            for key, value in pending[0].items():
                census[key] = census.get(key, 0) + int(value)
        stats["split_attempt_s"] += time.perf_counter() - started
        stats["split_attempts"] += 1
        for daughter in result or ():
            points, _ = daughter.all_points_and_elevation()
            velocities = mantle.flow_at(points, world.mantle_centers)
            stats["residuals"].append(merge_split._fit_residual_rms(points, velocities, daughter.omega))
        return result

    merge_split.maybe_split_plate = maybe_split_plate
    steer = getattr(merge_split, "_steer_cut_off_strong_lithosphere", None)
    if steer is not None:
        originals.append(("_steer_cut_off_strong_lithosphere", steer))

        def tapped_steer(points, *args):
            started = time.perf_counter()
            result = steer(points, *args)
            elapsed = time.perf_counter() - started
            areas, strength, cut_normal, band_sin = args[0], args[1], args[4], args[5]
            stats["steer_calls"] += 1
            stats["steer_needed"] += float((strength * areas)[np.abs(points @ cut_normal) < band_sin].sum()) > 0.0
            stats["steer_s"] += elapsed
            if elapsed > stats["steer_s_max"]:
                stats["steer_s_max"], stats["steer_max_nodes"] = elapsed, len(points)
            stats["steered"] += result[1] is not cut_normal
            pending.append(_gate_census(points, *args, result[1], last_craton[0]))
            return result

        merge_split._steer_cut_off_strong_lithosphere = tapped_steer
        last_craton = [None]
        cut_strength = merge_split._cut_strength
        originals.append(("_cut_strength", cut_strength))

        def tapped_cut_strength(craton_crust_m, hm):
            last_craton[0] = np.asarray(craton_crust_m, dtype=float)
            return cut_strength(craton_crust_m, hm)

        merge_split._cut_strength = tapped_cut_strength
    return stats, originals


# RFC 0003's sweep ranges for SPLIT_CUT_MAX_TILT_DEG and SPLIT_CUT_MAX_MISFIT, plus 45 deg.
CENSUS_TILTS_DEG = (10.0, 20.0, 30.0, 45.0)
CENSUS_MISFITS = (0.05, 0.1, 0.2)


def _cut_count(world) -> float:
    ledger = getattr(world, "craton_ledger", {})
    return ledger.get("split_cuts", 0.0) + ledger.get("failed_rift_cuts", 0.0)


def _gate_census(points, areas, strength, velocities, mask_a, cut_normal, band_sin, threshold, min_nodes, chosen, craton) -> dict:
    """For a cut that went ahead: did its flow-derived band hold craton, and did the chosen
    cut miss it. For each census tilt: did a craton-free candidate exist in the cone, which
    gates did those candidates pass, and was one admissible at each census misfit. Recomputes
    the steering's own arrays, read-only, on this run's trajectory."""
    def through(normal):
        return bool(np.any(craton[np.abs(points @ normal) < band_sin] > 0.0))

    row = {"cuts": 1, "n0_through_craton": through(cut_normal), "chosen_through_craton": through(chosen)}
    if not row["n0_through_craton"]:
        return row
    default_tilt = merge_split.SPLIT_CUT_MAX_TILT_DEG
    for tilt in CENSUS_TILTS_DEG:
        merge_split.SPLIT_CUT_MAX_TILT_DEG = tilt
        try:
            normals = merge_split._candidate_cut_normals(cut_normal)
        finally:
            merge_split.SPLIT_CUT_MAX_TILT_DEG = default_tilt
        craton_scores, misfit = merge_split._cut_scores(points, areas, (craton > 0.0).astype(float), mask_a, normals, band_sin)
        fits = merge_split._batched_half_fits(points, velocities, normals)
        free = craton_scores == 0.0
        gates = {
            "valid_and_size": fits["valid"] & (np.minimum(fits["count_a"], fits["count_b"]) >= min_nodes),
            "separation": np.linalg.norm(fits["pole_a"] - fits["pole_b"], axis=1) >= threshold,
            "residual": fits["rms"] <= (1.0 + merge_split.SPLIT_CUT_MAX_RESIDUAL_INCREASE) * fits["rms"][0],
        }
        others = np.logical_and.reduce(list(gates.values()))
        tag = f"tilt{tilt:g}"
        row[f"{tag}_craton_free_candidate"] = bool(free.any())
        for name, ok in gates.items():
            row[f"{tag}_craton_free_passes_{name}"] = bool((free & ok).any())
        for max_misfit in CENSUS_MISFITS:
            row[f"{tag}_misfit{max_misfit:g}_craton_free_admissible"] = bool((free & others & (misfit <= max_misfit)).any())
    return row


def _measure(seed: int, years: float, save_dir: Path) -> dict:
    stats, originals = _taps()
    started = time.perf_counter()
    try:
        run = rebaseline_rfcs._measure(seed, years, save_dir)
    finally:
        for name, original in originals:
            setattr(merge_split, name, original)
    residuals = np.asarray(stats.pop("residuals"))
    rows = run["checkpoints"]
    return {
        "seed": seed,
        "wall_s": time.perf_counter() - started,
        "per_gyr": run["per_gyr"],
        "rift_cuts": rows[-1]["rift_cuts"],
        "events": rows[-1]["events"],
        "mean_plate_count_per_gyr": [
            float(np.mean([row["plate_count"] for row in rows if start <= row["myr"] <= start + 1000]))
            for start in range(0, int(rows[-1]["myr"]), 1000)
        ],
        "craton_balance_error_km3": rows[-1]["craton_balance_error_km3"],
        "post_split_residual_cm_per_yr": {
            "count": len(residuals),
            "mean": float(residuals.mean() / mantle.cm_per_yr_to_rad_per_yr(1.0)) if len(residuals) else None,
        },
        "timing": stats,
        "save_sha256": run["save_sha256"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default="7,11,2024")
    parser.add_argument("--years", type=float, default=1e9)
    parser.add_argument("--set", action="append", default=[], metavar="NAME=VALUE")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--save-dir", type=Path, required=True)
    args = parser.parse_args()
    overrides = {}
    for item in args.set:
        name, value = item.split("=", 1)
        if not hasattr(merge_split, name):
            parser.error(f"merge_split has no {name}")
        overrides[name] = float(value)
        setattr(merge_split, name, float(value))
    app_root = Path(merge_split.__file__).resolve().parents[2]
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=app_root, text=True).strip()
    dirty = bool(subprocess.check_output(["git", "status", "--porcelain", "--", "backend/app"], cwd=app_root, text=True))
    result = {"app_path": str(Path(merge_split.__file__).resolve()), "model_commit": commit, "model_dirty": dirty, "overrides": overrides,
              "step_years": 5e6, "every_years": 100e6, "years": args.years, "runs": []}
    for seed in map(int, args.seeds.split(",")):
        result["runs"].append(_measure(seed, args.years, args.save_dir))
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
    main()
