#!/usr/bin/env python3
"""Anchor lake volume after one step against the lake solver's slice length (#853).

`solve_lakes_prognostic` integrates a step in `substeps` equal slices, and `hydrology` calls it
with the default of one. This script reloads the save for every configuration, takes one
`step_world(step_years)`, and reports the anchor's lake component after it. `--slice-years`
forces the solver's slice length by wrapping `lakes.solve_lakes_prognostic`; everything else in
the step (climate, hydrology inputs, terrain) is computed before the solve, so it is identical
across slice lengths: this is the "climate held fixed" arm.

The "climate refreshed" arm is `--refresh`: the same total span run as `span / slice` ordinary
steps of `slice` years each, so climate and hydrology are recomputed every slice. It reports
the anchor after the last step (and the per-step series with `--series`).

    backend/.venv/bin/python bin/debug/measure_lake_slice_length.py <save> \\
        --anchors analysis/rfc0016/target-anchor-window-d.json \\
        --step-years 5000 --slice-years 5000 2500 1000 200 --out /tmp/out.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import resource
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
MERGED_AREA_KM2 = 3.0e6  # window D/E merged inland sea spans 6-9M km2, the small basin under 1M


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("save", type=Path)
    ap.add_argument("--anchors", type=Path, required=True)
    ap.add_argument("--backend", type=Path, default=ROOT / "backend")
    ap.add_argument("--step-years", type=float, nargs="+", required=True)
    ap.add_argument("--slice-years", type=float, nargs="+", required=True)
    ap.add_argument("--refresh", action="store_true", help="run span/slice ordinary steps instead of substeps")
    ap.add_argument("--series", action="store_true", help="with --refresh, keep every step's anchor volume")
    ap.add_argument("--out", type=Path, required=True)
    return ap.parse_args()


args = parse_args()
sys.path.insert(0, str(args.backend.resolve()))
from scipy.spatial import cKDTree  # noqa: E402

from app import persistence  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402

data = args.save.read_bytes()
anchor_points = np.asarray(json.loads(args.anchors.read_text())["anchor_points"], dtype=float)
visible = hydrology.LAKE_MIN_VISIBLE_DEPTH_M
real_solve = lakes.solve_lakes_prognostic
forced = {"substeps": None, "last": None}


def solve_with_slices(forest, elevation, prev, years, climate, hypsometry, **kw):
    """Same solve, with the slice count set from the forced slice length."""
    n = forced["substeps"]
    kw["substeps"] = max(1, math.ceil(years / forced["slice"] - 1e-9)) if n is None else n
    forced["last"] = kw["substeps"]
    return real_solve(forest, elevation, prev, years, climate, hypsometry, **kw)


lakes.solve_lakes_prognostic = solve_with_slices


def load():
    world = persistence.load_world_bytes(data)
    world.ice_age_period_years = 0.0
    return world


def anchor_state(world) -> dict:
    fields = world.hydrology_cache
    points = np.asarray(fields.points)
    distance, nearest = cKDTree(points).query(anchor_points, k=1)
    wet = np.asarray(fields.lake_depth) > visible
    components = hydrology.lake_components(wet, fields.neighbor_idx)
    hyps = getattr(fields, "hypsometry", None) or lakes.CellHypsometry.whole_cells(fields.elevation)
    area_m2 = hyps.area_m2
    label = np.full(len(wet), -1, dtype=np.int64)
    for i, nodes in enumerate(components):
        label[np.asarray(nodes, dtype=np.int64)] = i
    anchor_label = int(label[nearest[0]])
    anchor_wet = anchor_label >= 0
    if not anchor_wet and (label >= 0).any():
        # Same rule as replay_rfc0016_river_water.target_lake: a dry anchor reports the nearest lake.
        wet_nodes = np.flatnonzero(label >= 0)
        anchor_label = int(label[wet_nodes[cKDTree(points[wet_nodes]).query(points[nearest[0]], k=1)[1]]])
    out = {
        "anchor_distance_km": float(distance[0] * 6371.0),
        "anchor_wet": anchor_wet,
        "world_lake_km3": float(np.dot(fields.lake_depth, area_m2) / 1e9),
    }
    if anchor_label >= 0:
        nodes = np.asarray(components[anchor_label], dtype=np.int64)
        out["anchor_lake_km3"] = float(np.dot(fields.lake_depth[nodes], area_m2[nodes]) / 1e9)
        out["anchor_lake_area_km2"] = float(area_m2[nodes].sum() / 1e6)
    else:
        out["anchor_lake_km3"] = 0.0
        out["anchor_lake_area_km2"] = 0.0
    out["merged"] = out["anchor_lake_area_km2"] >= MERGED_AREA_KM2
    ledger = getattr(fields, "lake_water_ledger", None)
    if ledger is not None:
        out["ledger_residual_norm"] = abs(float(ledger.residual_m3)) / max(float(ledger.external_supply_m3), 1.0)
    return out


records = []
for step_years in args.step_years:
    for slice_years in args.slice_years:
        if slice_years > step_years + 1e-9:
            continue
        world = load()
        started = time.perf_counter()
        rec = {"step_years": step_years, "slice_years": slice_years,
               "arm": "refreshed" if args.refresh else "fixed_climate"}
        if args.refresh:
            forced["substeps"] = 1
            count = int(round(step_years / slice_years))
            series = []
            for _ in range(count):
                world_mod.step_world(world, slice_years)
                if args.series:
                    series.append(anchor_state(world)["anchor_lake_km3"])
            rec["steps_run"] = count
            if args.series:
                rec["series_anchor_km3"] = series
        else:
            forced["substeps"] = None
            forced["slice"] = slice_years
            world_mod.step_world(world, step_years)
            rec["substeps"] = forced["last"]
        rec.update(anchor_state(world))
        rec["runtime_s"] = time.perf_counter() - started
        rec["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
        records.append(rec)
        print(json.dumps({k: v for k, v in rec.items() if k != "series_anchor_km3"}), flush=True)
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps({
            "issue": 853, "save": str(args.save), "save_sha256": hashlib.sha256(data).hexdigest(),
            "anchors": str(args.anchors), "merged_area_threshold_km2": MERGED_AREA_KM2,
            "records": records}, indent=1))
print(f"wrote {args.out}")
