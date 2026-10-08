#!/usr/bin/env python3
"""Where does the quad "pile-up" come from? (issue #289, found on #177)

On quad worlds, `avg_rotation_rate` 4x leaves less land than 0.25x but much more land
volume above sea level. The sweep's original land-volume metric weighted every node by the
nominal `node_area_m2` and summed every plate, so part of the rise was measurement: quad cells
aren't equal-area, and overlapping plates count the same ground twice. This probe splits the
metric up per (multiplier, seed, checkpoint):

- `land_vol_nominal_km3`: the pre-#289 sweep metric (`sweep_lib._land_volume_above_sea_nominal_km3`);
- `land_vol_exact_km3`: the same, weighted by `Plate.accounting_areas_m2`;
- `land_vol_exact_dedup_km3`: exact, with each node's area divided by the number of plates on
  it (`stats.overlap_area_weights`) -- the sweep's current `land_volume_above_sea_km3`;
- `represented_sphere_fraction` / `overlap_area_fraction`: total accounting area and the
  overlapped part of it, as fractions of the sphere;
- land area, mean land elevation and mean Hc on land (exact areas);
- the #273 continental-material ledger (`continental_ledger.inventories`): live surface
  volume, every source/sink account and the closure residual.

    backend/.venv/bin/python bin/debug/probe_rotation_pileup.py --seeds 350921662,57685824 \
        --out results/rotation_pileup_probe.jsonl
"""

from __future__ import annotations

import argparse
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import sweep_lib  # noqa: E402  (sets up sys.path for app)

import numpy as np  # noqa: E402

from app import continental_ledger, lithosphere, stats  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402

SPHERE_M2 = stats.SPHERE_AREA_M2


def measure(world) -> dict:
    spacing = line_spacing_rad(world.node_density)
    nominal = lithosphere.node_area_m2(spacing)
    live = [p for p in world.plates if p.node_count() > 0]
    weights = stats.overlap_area_weights(world)
    sea = world.sea_level_m
    out = dict.fromkeys(
        ("land_vol_nominal_km3", "land_vol_exact_km3", "land_vol_exact_dedup_km3", "land_area_exact_km2",
         "land_area_dedup_km2", "represented_m2", "overlap_m2", "land_elev_area_m3", "land_hc_area_m3"),
        0.0,
    )
    for p in live:
        _, elev = p.all_points_and_elevation()
        areas = p.accounting_areas_m2(spacing)
        hc = p.collect("crustal_thickness_m")
        weight = weights[p.plate_id]
        mask = weight < 1.0
        above = np.clip(elev - sea, 0.0, None)
        land = elev > sea
        out["land_vol_nominal_km3"] += float(above.sum()) * nominal / 1e9
        out["land_vol_exact_km3"] += float(np.dot(above, areas)) / 1e9
        out["land_vol_exact_dedup_km3"] += float(np.dot(above, areas * weight)) / 1e9
        out["land_area_exact_km2"] += float(areas[land].sum()) / 1e6
        out["land_area_dedup_km2"] += float((areas * weight)[land].sum()) / 1e6
        out["represented_m2"] += float(areas.sum())
        out["overlap_m2"] += float(areas[mask].sum())
        out["land_elev_area_m3"] += float(np.dot(above[land], areas[land]))
        out["land_hc_area_m3"] += float(np.dot(hc[land], areas[land]))
    land_m2 = out["land_area_exact_km2"] * 1e6
    record = {
        "land_vol_nominal_km3": out["land_vol_nominal_km3"],
        "land_vol_exact_km3": out["land_vol_exact_km3"],
        "land_vol_exact_dedup_km3": out["land_vol_exact_dedup_km3"],
        "land_area_exact_km2": out["land_area_exact_km2"],
        "land_area_dedup_km2": out["land_area_dedup_km2"],
        "represented_sphere_fraction": out["represented_m2"] / SPHERE_M2,
        "overlap_area_fraction": out["overlap_m2"] / SPHERE_M2,
        "mean_land_height_above_sea_m": out["land_elev_area_m3"] / land_m2 if land_m2 else None,
        "mean_land_hc_m": out["land_hc_area_m3"] / land_m2 if land_m2 else None,
        "sea_level_m": sea,
    }
    record.update({f"ledger_{k}": v for k, v in continental_ledger.inventories(world).items()})
    return record


def run(multiplier: float, seed: int, checkpoints: tuple[int, ...]) -> list[dict]:
    sweep_lib._apply_rotation_rate_overrides(
        sweep_lib.BASELINE_AVG_RATE_CM_YR * multiplier, sweep_lib.BASELINE_MAX_RATE_CM_YR
    )
    world = generate_world(seed=seed, node_density=1.0)
    rows, done = [], 0
    for cp in checkpoints:
        while done < cp:
            step_world(world, years=sweep_lib.STEP_YEARS)
            done += sweep_lib.STEP_YEARS
        rows.append({"multiplier": multiplier, "seed": seed, "checkpoint_years": cp, **measure(world)})
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seeds", required=True)
    parser.add_argument("--multipliers", default="0.25,4.0")
    parser.add_argument("--checkpoints", default="30,120", help="Myr, comma-separated")
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    cps = tuple(int(float(c) * 1e6) for c in args.checkpoints.split(","))
    jobs = [(float(m), int(seed), cps) for m in args.multipliers.split(",") for seed in args.seeds.split(",")]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("a") as f, ProcessPoolExecutor(max_workers=args.workers) as pool:
        for rows in pool.map(run, *zip(*jobs)):
            for r in rows:
                f.write(json.dumps(r) + "\n")
            f.flush()
            print(f"done x{rows[0]['multiplier']} seed={rows[0]['seed']}", flush=True)


if __name__ == "__main__":
    main()
