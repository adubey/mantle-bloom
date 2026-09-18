#!/usr/bin/env python3
"""Per-step land-area-vs-crust-volume ledger for GitHub issue #177 -- follow-up to #171's
parameter sweep, which found that average rotation rate (mantle.MANTLE_FLOW_REFERENCE_RATE) is
the strongest single lever on land-fraction decline, and that unlike every other knob swept, its
effect *compounds* over time rather than settling. #177's own hypothesis: this is a land-
creation/land-destruction imbalance (destruction scaling with collision frequency, which rises
with rotation rate, while creation may not scale the same way), not just a rate to retune.

First attempt at this (see git history) resampled elevation onto a fixed world-frame lat/lon
grid, before vs after each step, and diffed which cells crossed the sea-level threshold. That
measures apparent motion, not tectonic creation/destruction: at this model's plate speeds (the
4x-multiplier run moves a plate ~1600 km in a single 10 My step), a continent's coastline sweeps
across most of the fixed grid every step, so "cells that flipped land<->ocean" is dominated by
translation churn (up to ~20M km^2/step gained+lost, on a ~50M km^2 land base) with the real
tectonic signal buried in it.

This version instead tracks whole-world aggregates that are constructed to be translation-
invariant by definition (each node's own elevation/crust-type/thickness -- what a rigid rotation
of a plate's frame never touches):
  - land_node_count / continental_node_count / total_node_count (see
    plates.collect_all_crust_type_view_codes -- CRUST_TYPE_INHERIT-resolved per-node continental
    flag), each * lithosphere.node_area_m2 for land_area_km2 / continental_area_km2
  - total_continental_crust_volume_km3 (stats.py, unchanged) and mean_continental_hc_m derived
    from it -- if lithosphere_plate.py's row-drop/suture-retreat machinery is conserving *volume*
    (thrusting a dropped row's Hc onto the plate's surviving edge, per the
    LEADING_ROW_CONTESTED_FRACTION / SUTURE_ACCRETION_SPREAD_NODES comments in
    lithosphere_plate.py) while *area* shrinks (fewer nodes = smaller footprint, even if each
    survivor is thicker), continental_area declining much faster than continental_volume -- and
    mean_continental_hc_m rising -- is the direct fingerprint of exactly that: area being lost to
    node attrition at retreating margins even though the crust itself isn't disappearing.

Usage (run from anywhere; the venv with numpy/scipy/etc is backend/.venv):
    backend/.venv/bin/python bin/debug/rotation_ledger.py
    backend/.venv/bin/python bin/debug/rotation_ledger.py --seed 829071382 --multipliers 0.25,4.0
    backend/.venv/bin/python bin/debug/rotation_ledger.py --years 180000000 --step-years 10000000
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import numpy as np  # noqa: E402

from app import elevation_lines, lithosphere, plates, stats  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app.world import World, generate_world, step_world  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sweep_lib import (  # noqa: E402
    BASELINE_AVG_RATE_CM_YR,
    BASELINE_MAX_RATE_CM_YR,
    _apply_rotation_rate_overrides,
    _land_volume_above_sea_km3,
)

DEFAULT_SEED = 829071382  # #171's first "reported-problem" seed
DEFAULT_MULTIPLIERS = (0.25, 4.0)  # the two ends of #171's swept avg-rotation-rate range
DEFAULT_STEP_YEARS = 10_000_000
DEFAULT_TOTAL_YEARS = 120_000_000
DEFAULT_NODE_DENSITY = 1.0


def _effective_is_continental(world: World) -> np.ndarray:
    """Per-node continental flag, world order (matches `plates.collect_all_points`'s own
    per-plate/per-line concatenation order) -- CRUST_TYPE_INHERIT nodes resolve to their own
    plate's nominal `crust_type`, exactly `plates.collect_all_crust_type_view_codes`'s own
    per-plate logic, just returned as a plain bool mask instead of a 4-way view code."""
    chunks = []
    for plate in world.plates:
        n = plate.node_count()
        if n == 0:
            continue
        plate_is_continental = plate.crust_type == "continental"
        chunks.append(
            elevation_lines.effective_is_continental_from_codes(plate.collect("crust_type_code"), plate_is_continental)
        )
    return np.concatenate(chunks, axis=0) if chunks else np.zeros(0, dtype=bool)


def snapshot(world: World) -> dict:
    area_m2 = lithosphere.node_area_m2(line_spacing_rad(world.node_density))
    is_continental = _effective_is_continental(world)
    collected = plates.collect_all_points(world.plates)
    if collected is None:
        elevation = np.zeros(0)
    else:
        _, elevation, _ = collected
    is_land = elevation > world.sea_level_m

    total_node_count = int(is_continental.shape[0])
    continental_node_count = int(np.count_nonzero(is_continental))
    land_node_count = int(np.count_nonzero(is_land))
    continental_land_node_count = int(np.count_nonzero(is_land & is_continental))

    snap = stats.compute_stats(world)
    continental_crust_volume_km3 = snap["total_continental_crust_volume_km3"]
    mean_continental_hc_m = (
        (continental_crust_volume_km3 * 1.0e9) / (continental_node_count * area_m2)
        if continental_node_count > 0
        else None
    )

    return {
        "total_node_count": total_node_count,
        "continental_node_count": continental_node_count,
        "oceanic_node_count": total_node_count - continental_node_count,
        "land_node_count": land_node_count,
        "continental_land_node_count": continental_land_node_count,
        "land_area_km2": land_node_count * area_m2 / 1.0e6,
        "continental_area_km2": continental_node_count * area_m2 / 1.0e6,
        "continental_crust_volume_km3": continental_crust_volume_km3,
        "mean_continental_hc_m": mean_continental_hc_m,
        "land_volume_above_sea_km3": _land_volume_above_sea_km3(world),
        "land_fraction": snap["land_fraction"],
    }


def run_ledger(
    seed: int,
    multiplier: float,
    node_density: float,
    step_years: int,
    total_years: int,
) -> list[dict]:
    _apply_rotation_rate_overrides(BASELINE_AVG_RATE_CM_YR * multiplier, BASELINE_MAX_RATE_CM_YR)
    world = generate_world(seed=seed, node_density=node_density)

    records = []
    prev = snapshot(world)
    years_done = 0.0
    while years_done < total_years:
        step = min(step_years, total_years - years_done)
        step_world(world, years=step)
        years_done += step
        cur = snapshot(world)

        record = {"seed": seed, "multiplier": multiplier, "years_elapsed": years_done}
        for key, value in cur.items():
            record[key] = value
            prior = prev[key]
            if value is not None and prior is not None:
                record[f"delta_{key}"] = value - prior
        records.append(record)
        prev = cur
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--multipliers", type=str, default=",".join(str(m) for m in DEFAULT_MULTIPLIERS))
    parser.add_argument("--node-density", type=float, default=DEFAULT_NODE_DENSITY)
    parser.add_argument("--step-years", type=int, default=DEFAULT_STEP_YEARS)
    parser.add_argument("--years", type=int, default=DEFAULT_TOTAL_YEARS)
    parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "results" / "rotation_ledger.jsonl")
    args = parser.parse_args()

    multipliers = [float(m) for m in args.multipliers.split(",")]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    all_records = []
    for multiplier in multipliers:
        print(f"--- seed={args.seed} multiplier={multiplier}x ---", file=sys.stderr)
        records = run_ledger(args.seed, multiplier, args.node_density, args.step_years, args.years)
        all_records.extend(records)
        for r in records:
            print(
                f"  {r['years_elapsed']/1e6:6.0f} My  "
                f"land_area={r['land_area_km2']/1e6:7.3f}M km2 ({r['delta_land_area_km2']/1e6:+6.3f})  "
                f"cont_area={r['continental_area_km2']/1e6:7.3f}M km2 ({r['delta_continental_area_km2']/1e6:+6.3f})  "
                f"cont_vol={r['continental_crust_volume_km3']/1e6:7.2f}M km3 ({r['delta_continental_crust_volume_km3']/1e6:+6.3f})  "
                f"mean_hc={r['mean_continental_hc_m']:7.0f}m ({r['delta_mean_continental_hc_m']:+6.1f})  "
                f"cont_nodes={r['continental_node_count']} ({r['delta_continental_node_count']:+d})",
                file=sys.stderr,
            )

    with open(args.out, "a") as f:
        for r in all_records:
            f.write(json.dumps(r) + "\n")
    print(f"wrote {len(all_records)} records to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
