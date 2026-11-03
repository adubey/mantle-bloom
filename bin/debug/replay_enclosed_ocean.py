#!/usr/bin/env python3
"""Issue #758: replay a save and record, every step, how below-sea-level water is classed.

Per step it records sea level and the eustasy closure (ocean volume at the solved level plus
trapped water, against the conserved total), the ocean's node count and area, whole-world
lake water (liquid plus ice), the normalized lake and land ledger residuals, every enclosed
below-sea-level body of at least `--min-body-km2` with its area, oceanic-crust share and class,
and precipitation over land near the save's largest enclosed body (land climate cells within
`--near-km`, and beyond `--exclude-km`, of where that body's cells sat at the start).

A **class change** is a connected group of cells, below sea level both before and after a
step, whose `is_ocean` flipped in that step, covering at least `--min-body-km2`. Shoreline
moves (cells crossing sea level) aren't class changes. Last step's cells are matched to this
step's by nearest position. The summary counts them, split by direction.

Run from `backend/` so the worktree's environment is active:

    .venv/bin/python ../bin/debug/replay_enclosed_ocean.py <save> --myr 5 \\
        --step-years 100000 --out <scratch>/branch.json
    .venv/bin/python ../bin/debug/replay_enclosed_ocean.py <save> --myr 5 \\
        --step-years 100000 --backend <baseline>/backend --out <scratch>/main.json
"""

from __future__ import annotations

import json
import resource
import time

import numpy as np
from replay_lib import load_replay, parse_replay_args
from scipy.spatial import cKDTree


def _arguments(ap):
    ap.add_argument("--min-body-km2", type=float, default=200_000.0)
    ap.add_argument("--near-km", type=float, default=1_000.0)
    ap.add_argument("--exclude-km", type=float, default=100.0)


args = parse_replay_args(__doc__, _arguments)
world, result = load_replay(args, issue=758)

from app import plates  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.elevation_lines import PLANET_RADIUS_KM, line_spacing_rad  # noqa: E402
from app.hydroclimate import eustasy, hydrology  # noqa: E402

world.ice_age_period_years = 0.0
spacing = line_spacing_rad(world.node_density)


def node_state():
    hydro = world.hydrology_cache
    order = hydro.plates_in_order
    elevation = plates.collect_all_elevation(order)
    areas = plates.collect_all_accounting_areas_m2(order, spacing)
    oceanic = ~np.concatenate([plates.node_is_continental(p) for p in order])
    below = elevation <= world.sea_level_m
    labels = hydrology._below_sea_components(below, hydro.neighbor_idx)
    return hydro, elevation, areas, oceanic, below, labels


def enclosed_bodies(hydro, areas, oceanic, below, labels):
    lab = labels[below]
    area = np.bincount(lab, weights=areas[below])
    oceanic_area = np.bincount(lab, weights=(areas * oceanic)[below])
    ocean_area = np.bincount(lab, weights=(areas * hydro.is_ocean)[below])
    largest = int(np.bincount(lab).argmax())
    keep = np.flatnonzero(area >= args.min_body_km2 * 1e6)
    return [
        {
            "label": int(k),
            "area_km2": float(area[k] / 1e6),
            "oceanic_share": float(oceanic_area[k] / area[k]),
            "is_ocean": bool(ocean_area[k] >= 0.5 * area[k]),
        }
        for k in sorted(keep, key=lambda k: -area[k]) if k != largest
    ]


def class_changes(prev, hydro, areas, below):
    """Groups of cells below sea level before and after this step whose class flipped."""
    prev_points, prev_ocean, prev_below = prev
    _, idx = cKDTree(prev_points).query(hydro.points, workers=-1)
    flipped = below & prev_below[idx] & (prev_ocean[idx] != hydro.is_ocean)
    if not flipped.any():
        return []
    groups = hydrology._below_sea_components(flipped, hydro.neighbor_idx)
    lab = groups[flipped]
    area = np.bincount(lab, weights=areas[flipped])
    to_ocean = np.bincount(lab, weights=(areas * hydro.is_ocean)[flipped])
    return [
        {"area_km2": float(area[k] / 1e6), "to": "ocean" if to_ocean[k] >= 0.5 * area[k] else "lake"}
        for k in np.flatnonzero(area >= args.min_body_km2 * 1e6)
    ]


def near_body_cells(hydro, below, labels, bodies):
    """Climate cells within `--near-km` of the largest enclosed body's cells at the start, but
    more than `--exclude-km` from them, so the body's own floor (land on one arm, water on the
    other) isn't in the sample."""
    fields = world.climate_cache
    if not bodies or fields is None:
        return None
    body_points = hydro.points[below & (labels == bodies[0]["label"])]
    chord = 2.0 * np.sin(args.near_km / PLANET_RADIUS_KM / 2.0)
    inner = 2.0 * np.sin(args.exclude_km / PLANET_RADIUS_KM / 2.0)
    dist, _ = cKDTree(body_points).query(fields.world_xyz.reshape(-1, 3), distance_upper_bound=chord, workers=-1)
    return ((dist > inner) & np.isfinite(dist)).reshape(fields.precipitation_mm.shape)


def near_land_precipitation(near):
    fields = world.climate_cache
    if near is None or fields is None or fields.precipitation_mm.shape != near.shape:
        return None
    land = near & ~fields.is_ocean & (fields.lake_depth_m < hydrology.LAKE_MIN_VISIBLE_DEPTH_M)
    weight = np.cos(np.radians(fields.lat_deg))[:, None] * land
    return {
        "land_cells": int(land.sum()),
        "mean_mm": float((fields.precipitation_mm * weight).sum() / weight.sum()) if weight.sum() > 0 else None,
    }


def snapshot(step, seconds, changes, near):
    hydro, elevation, areas, oceanic, below, labels = node_state()
    lake_liquid = plates.collect_all_lake_depth(hydro.plates_in_order)
    lake_ice = plates.collect_all_lake_ice(hydro.plates_in_order)
    trapped = eustasy.trapped_water_volume_m3(world)
    ocean_volume = eustasy.water_volume_for_sea_level(world, world.sea_level_m)
    row = {
        "step": step,
        "elapsed_myr": world.elapsed_years / 1e6,
        "sea_level_m": world.sea_level_m,
        "eustasy_closure_m3": ocean_volume + trapped - world.ocean_water_volume_m3,
        "eustasy_tolerance_m3": eustasy._SOLVE_TOLERANCE_M * float(areas.mean()),
        "ocean_nodes": int(hydro.is_ocean.sum()),
        "ocean_area_km2": float(areas[hydro.is_ocean].sum() / 1e6),
        "lake_liquid_km3": float(np.dot(areas, np.clip(lake_liquid, 0.0, None)) / 1e9),
        "lake_ice_km3": float(np.dot(areas, np.clip(lake_ice, 0.0, None)) / 1e9),
        "enclosed_bodies": enclosed_bodies(hydro, areas, oceanic, below, labels),
        "class_changes": changes,
        "near_land_precipitation": near_land_precipitation(near),
        "step_seconds": seconds,
    }
    lake, land = hydro.lake_water_ledger, getattr(hydro, "land_water_ledger", None)
    if lake is not None and step > 0:
        scale = max(abs(lake.external_supply_m3), abs(lake.overflow_out_m3), 1.0)
        row["lake_residual_normalized"] = float(lake.residual_m3 / scale)
    if land is not None and step > 0:
        row["land_residual_normalized"] = float(land.residual_m3 / max(land.rain_m3 + land.meltwater_m3, 1.0))
    return row, (hydro.points.copy(), hydro.is_ocean.copy(), below.copy())


steps = int(round(args.myr * 1e6 / args.step_years))
result.update(min_body_km2=args.min_body_km2, near_km=args.near_km, exclude_km=args.exclude_km, steps=[])
args.out.parent.mkdir(parents=True, exist_ok=True)
hydro0, _, areas0, oceanic0, below0, labels0 = node_state()
near = near_body_cells(hydro0, below0, labels0, enclosed_bodies(hydro0, areas0, oceanic0, below0, labels0))
row, prev = snapshot(0, 0.0, [], near)
result["steps"].append(row)
print(json.dumps({k: row[k] for k in ("elapsed_myr", "sea_level_m", "ocean_nodes")}), flush=True)
for i in range(1, steps + 1):
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    seconds = time.perf_counter() - started
    hydro, _, areas, _, below, _ = node_state()
    row, prev = snapshot(i, seconds, class_changes(prev, hydro, areas, below), near)
    result["steps"].append(row)
    print(json.dumps({k: row[k] for k in ("elapsed_myr", "sea_level_m", "ocean_nodes", "class_changes")}), flush=True)
    args.out.write_text(json.dumps(result, indent=1))

changes = [c for r in result["steps"] for c in r["class_changes"]]
result["summary"] = {
    "class_changes": len(changes),
    "to_ocean": sum(c["to"] == "ocean" for c in changes),
    "to_lake": sum(c["to"] == "lake" for c in changes),
    "steps_with_change": sum(bool(r["class_changes"]) for r in result["steps"]),
    "max_abs_eustasy_closure_over_tolerance": max(
        abs(r["eustasy_closure_m3"]) / r["eustasy_tolerance_m3"] for r in result["steps"]
    ),
    "max_abs_lake_residual_normalized": max(abs(r.get("lake_residual_normalized", 0.0)) for r in result["steps"]),
    "max_abs_land_residual_normalized": max(abs(r.get("land_residual_normalized", 0.0)) for r in result["steps"]),
    "median_step_seconds": float(np.median([r["step_seconds"] for r in result["steps"][1:]])) if steps else None,
    "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6,
}
args.out.write_text(json.dumps(result, indent=1))
print(json.dumps(result["summary"]))
print(f"wrote {args.out}")
