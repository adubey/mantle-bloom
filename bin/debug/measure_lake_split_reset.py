#!/usr/bin/env python3
"""Issue #521: replay a saved world and measure the lake split reset, on `main` or on the branch.

Steps the save `--steps` times by `--step-years` and records, per step:

- global visible lake volume (km^3) and area (km^2), and the count of `merge` and `split`
  lake events;
- each tracked body's wet extent (nodes), volume, and the smallest lake in the tree that holds
  it: that lake's level, saddle (`min_depth`), floor and rim, and whether the level rests on
  the saddle (within 1e-6 m) with evaporation at that level against its inflow;
- discharge (m^3/s) arriving at land sinks (`flow_target == -1`), split into frozen (annual
  mean air below `hydrology.FREEZE_POINT_C`) and unfrozen, and reaching the ocean;
- the step's liquid precipitation and Budyko runoff over land (m^3/s), the water routing starts
  from, read off `hydrology`'s `breaching.runoff_mm` call;
- step runtime (s) and peak RSS.

Tracked bodies are connected components of visible lake water (`LAKE_MIN_VISIBLE_DEPTH_M`) in
the *saved* state, numbered in order of their lowest node index (the ids #519 used, e.g. sea 26
and lake 13 of `seed308161711-121.5Myr`). Each step re-locates a body as the wet component
that most of its last wet nodes fall in, or by its saved extent once it has dried. Without
`--track`, the three largest bodies by saved volume are tracked.

`--backend` points at the backend to import, so the same script measures both trees:

    cd backend
    .venv/bin/python -I ../bin/debug/measure_lake_split_reset.py <save> --steps 10 \\
        --track 26 13 --out <scratch>/branch.json
    .venv/bin/python -I ../bin/debug/measure_lake_split_reset.py <save> --steps 10 \\
        --track 26 13 --backend <baseline-worktree>/backend --out <scratch>/main.json
"""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("save", type=Path)
ap.add_argument("--steps", type=int, default=10)
ap.add_argument("--step-years", type=float, default=100_000.0)
ap.add_argument("--track", type=int, nargs="*")
ap.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
ap.add_argument("--out", type=Path, required=True)
args = ap.parse_args()
sys.path.insert(0, str(args.backend.resolve()))

from scipy.spatial import cKDTree  # noqa: E402

from app import persistence  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402
from app.hydroclimate.erosion import climate_grid_indices  # noqa: E402

VIS = hydrology.LAKE_MIN_VISIBLE_DEPTH_M
SECONDS_PER_YEAR = hydrology.SECONDS_PER_YEAR
RSS_PER_MB = 2**20 if sys.platform == "darwin" else 2**10  # ru_maxrss: bytes on macOS, KiB on Linux

# The resolve inputs, kept so each tracked basin's inflow and evaporation can be read back.
captured: dict = {}
_orig_resolve_lakes = lakes.resolve_lakes


def _capture(*args):
    # `compute_hydrology` passes climate and hypsometry last, both before #539 (when an
    # `is_frozen` array came before them) and after, so either backend can be measured.
    captured.update(climate=args[-2], hyps=args[-1])
    return _orig_resolve_lakes(*args)


lakes.resolve_lakes = _capture
hydrology.lakes.resolve_lakes = _capture
_orig_runoff_mm = hydrology.breaching.runoff_mm


def _capture_runoff(precipitation_mm, potential_evaporation_mm):
    runoff = _orig_runoff_mm(precipitation_mm, potential_evaporation_mm)
    if sys._getframe(1).f_globals.get("__name__") != hydrology.__name__:
        return runoff  # breaching's own breach sizing calls it too
    captured.update(precip_mm=np.asarray(precipitation_mm).copy(), runoff_mm=np.asarray(runoff).copy())
    return runoff


hydrology.breaching.runoff_mm = _capture_runoff


def components(f):
    comps = hydrology.lake_components(f.lake_depth > VIS, f.neighbor_idx)
    comps.sort(key=lambda m: int(m.min()))
    return comps


def smallest_containing(forest, nodes):
    """The deepest lake in the tree whose members include every node in `nodes`."""
    want = set(nodes.tolist())
    for root in forest:
        if not want.issubset(root.members.tolist()):
            continue
        lake = root
        while True:
            inner = next((c for c in lake.children if want.issubset(c.members.tolist())), None)
            if inner is None:
                return lake
            lake = inner
    return None


def frozen_sink_budget(world, f):
    c = world.climate_cache
    h, w = c.air_temperature_c.shape
    row, col = climate_grid_indices(f.points, h, w)
    frozen = c.air_temperature_c[row, col] < hydrology.FREEZE_POINT_C
    inflow = np.zeros(len(f.flow_target))
    routed = f.flow_target >= 0
    np.add.at(inflow, f.flow_target[routed], f.discharge_m3_s[routed])
    sink = ~routed & ~f.is_ocean
    return {
        "frozen_sink_inflow_m3_s": float(inflow[sink & frozen].sum()),
        "unfrozen_sink_inflow_m3_s": float(inflow[sink & ~frozen].sum()),
        "ocean_discharge_m3_s": float(f.discharge_m3_s[f.is_ocean].sum()),
    }


def land_source(f):
    land = ~f.is_ocean
    area = f.hypsometry.area_m2
    to_m3_s = area[land] / 1000.0 / SECONDS_PER_YEAR
    return {
        "land_precip_m3_s": float((captured["precip_mm"][land] * to_m3_s).sum()),
        "land_runoff_m3_s": float((captured["runoff_mm"][land] * to_m3_s).sum()),
    }


def body(f, wet, footprint):
    area = f.hypsometry.area_m2
    out = {"wet_nodes": int(len(wet)), "volume_km3": float((f.lake_depth[wet] * area[wet]).sum() / 1e9)}
    lake = smallest_containing(f.lake_forest, wet if len(wet) else footprint)
    if lake is None:
        out["basin"] = None
        return out
    hyps, clim = captured["hyps"], captured["climate"]
    m = lake.members
    level = lake.current_water_elevation
    inflow = float(clim.inflow_m3_per_yr[m].sum())
    evap = lakes.lake_evaporation_m3_per_yr(level, m, hyps, clim)
    out["basin"] = {
        "members": int(len(m)), "children": len(lake.children),
        "level_m": level, "saddle_m": lake.min_depth, "floor_m": lake.floor_elevation, "rim_m": lake.max_depth,
        "on_saddle": bool(lake.children) and abs(level - lake.min_depth) <= 1e-6,
        "evap_over_inflow": evap / inflow if inflow > 0 else None,
    }
    return out


world = persistence.load_world_bytes(args.save.read_bytes())
f = world.hydrology_cache
comps = components(f)
saved_volume = [float((f.lake_depth[c] * f.hypsometry.area_m2[c]).sum() / 1e9) for c in comps]
track_ids = args.track or sorted(range(len(comps)), key=lambda i: -saved_volume[i])[:3]
last_wet = {tid: f.points[comps[tid]] for tid in track_ids}
footprint_xyz = {tid: f.points[comps[tid]] for tid in track_ids}
result = {
    "elapsed_years": world.elapsed_years, "nodes": len(f.elevation), "steps": [],
    "saved_tracks": {tid: {"wet_nodes": int(len(comps[tid])), "volume_km3": saved_volume[tid]} for tid in track_ids},
}
print(json.dumps({k: v for k, v in result.items() if k != "steps"}), flush=True)

for step in range(1, args.steps + 1):
    start = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    runtime_s = time.perf_counter() - start
    f = world.hydrology_cache
    area = f.hypsometry.area_m2
    wet_mask = f.lake_depth > VIS
    comps = components(f)
    comp_of = {int(n): i for i, c in enumerate(comps) for n in c.tolist()}
    tree = cKDTree(f.points)
    events = f.lake_events or []
    rec = {
        "step": step, "elapsed_years": world.elapsed_years, "runtime_s": runtime_s,
        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / RSS_PER_MB,
        "lake_volume_km3": float((f.lake_depth * area).sum() / 1e9),
        "lake_area_km2": float(area[wet_mask].sum() / 1e6),
        "merge_events": sum(e.kind == "merge" for e in events),
        "split_events": sum(e.kind == "split" for e in events),
        **frozen_sink_budget(world, f),
        **land_source(f),
        "tracks": {},
    }
    for tid in track_ids:
        _, near = tree.query(last_wet[tid])
        _, foot = tree.query(footprint_xyz[tid])
        hits = [comp_of[int(i)] for i in near if int(i) in comp_of]
        wet = comps[max(set(hits), key=hits.count)] if hits else np.zeros(0, dtype=np.int64)
        rec["tracks"][tid] = body(f, wet, np.unique(foot))
        if len(wet):
            last_wet[tid] = f.points[wet]
    result["steps"].append(rec)
    print(json.dumps(rec), flush=True)
    args.out.write_text(json.dumps(result, indent=1))
