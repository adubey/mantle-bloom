#!/usr/bin/env python3
"""Issue #539 (RFC 0014 D2(b)): replay a saved world and measure liquid water at cold nodes.

Steps the save `--steps` times by `--step-years` and records, per step:

- liquid water deleted by freezing (m^3/s): the ice ledger's frozen river ice and frozen lake
  water, converted back to water at `hydrology.WATER_TO_ICE` (0 once #539 lands);
- inflow (m^3/s) to land sinks (`flow_target == -1`) whose annual-mean air temperature is below
  `hydrology.FREEZE_POINT_C`, and their count; before #539 every cold land node is such a sink;
- the five largest cold sinks: their inflow, and how much of it arrived through a spill (a
  `should_spill` node routed onward, counted at the spill nearest the sink so nested spills
  aren't counted twice; transit losses below the spill are ignored, so it's an upper bound),
  and how much of that came from basins that flood and then empty: the spilling sink held
  visible water this step and holds none, and doesn't spill, next step;
- cold sinks that spill (`should_spill` on a cold land node);
- discharge reaching the ocean (m^3/s);
- visible lake volume (km^3) and area (km^2), all lakes and cold lakes (cold air at the cell);
- glacier volume (km^3);
- breached pits (nodes with a breach notch cut this step);
- floor-draining lake leaves (`Lake.drains_at_floor`), and how many switched on or off since
  the previous step, matched by position within `--match-km`;
- step runtime (s) and peak RSS.

Nodes move with their plates, so a node is matched across steps by position (nearest node
within `--match-km`), never by index.

`--backend` points at the backend to import, so the same script measures both trees:

    cd backend
    .venv/bin/python -I ../bin/debug/measure_cold_sink_water.py <save> --steps 10 \\
        --out <scratch>/branch.json
    .venv/bin/python -I ../bin/debug/measure_cold_sink_water.py <save> --steps 10 \\
        --backend <baseline-worktree>/backend --out <scratch>/main.json
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
ap.add_argument("--match-km", type=float, default=20.0)
ap.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
ap.add_argument("--out", type=Path, required=True)
args = ap.parse_args()
sys.path.insert(0, str(args.backend.resolve()))

from scipy.spatial import cKDTree  # noqa: E402

from app import persistence  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.plates import PLANET_RADIUS_KM  # noqa: E402
from app.hydroclimate import hydrology, lakes  # noqa: E402
from app.hydroclimate.erosion import climate_grid_indices  # noqa: E402

VIS = hydrology.LAKE_MIN_VISIBLE_DEPTH_M
SECONDS_PER_YEAR = hydrology.SECONDS_PER_YEAR
RSS_PER_MB = 2**20 if sys.platform == "darwin" else 2**10  # ru_maxrss: bytes on macOS, KiB on Linux
MATCH = args.match_km / PLANET_RADIUS_KM  # chord on the unit sphere, near enough for small distances


def cold_mask(world, f):
    c = world.climate_cache
    h, w = c.air_temperature_c.shape
    row, col = climate_grid_indices(f.points, h, w)
    return c.air_temperature_c[row, col] < hydrology.FREEZE_POINT_C


def deleted_m3_s(f, years):
    ledger = f.ice.ledger
    ice_m3 = ledger.frozen_river_ice_m3 + ledger.frozen_lake_m3
    return ice_m3 / hydrology.WATER_TO_ICE / years / SECONDS_PER_YEAR


def spill_delivery(f):
    """For every terminal node, the discharge delivered through the spills nearest it, and
    those spills' node indices."""
    spilling = np.flatnonzero(f.should_spill & (f.flow_target >= 0))
    is_spill = np.zeros(len(f.flow_target), dtype=bool)
    is_spill[spilling] = True
    delivered: dict[int, list[int]] = {}
    for start in spilling.tolist():
        node, hops = int(f.flow_target[start]), 0
        while node >= 0 and not is_spill[node] and hops <= len(f.flow_target):
            nxt = int(f.flow_target[node])
            if nxt < 0:
                break
            node, hops = nxt, hops + 1
        if node >= 0 and not is_spill[node]:
            delivered.setdefault(node, []).append(start)
    return delivered


def leaves(forest):
    for lake in lakes.iter_all_lakes(forest):
        if not lake.children:
            yield lake


def matched(tree, xyz):
    if len(xyz) == 0:
        return np.zeros(0, dtype=np.int64), np.zeros(0, dtype=bool)
    dist, idx = tree.query(xyz)
    return idx, dist <= MATCH


world = persistence.load_world_bytes(args.save.read_bytes())
result = {"elapsed_years": world.elapsed_years, "nodes": len(world.hydrology_cache.elevation), "steps": []}
print(json.dumps({k: v for k, v in result.items() if k != "steps"}), flush=True)
prev_sources: list[dict] = []  # last step's top-sink spill sources, judged once this step is known
prev_drain_xyz = None

for step in range(1, args.steps + 1):
    start = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    runtime_s = time.perf_counter() - start
    f = world.hydrology_cache
    area = f.hypsometry.area_m2
    cold = cold_mask(world, f)
    land = ~f.is_ocean
    wet = f.lake_depth > VIS
    tree = cKDTree(f.points)

    inflow = np.zeros(len(f.flow_target))
    routed = f.flow_target >= 0
    np.add.at(inflow, f.flow_target[routed], f.discharge_m3_s[routed])
    cold_sink = ~routed & land & cold

    # Judge last step's spill sources: did each basin empty?
    if prev_sources:
        idx, ok = matched(tree, np.array([s["xyz"] for s in prev_sources]))
        emptied = ok & ~wet[idx] & ~f.should_spill[idx]
        for s, e in zip(prev_sources, emptied.tolist()):
            s["sink"]["flood_then_empty_m3_s"] += s["through_m3_s"] if (s["was_wet"] and e) else 0.0
        result["steps"][-1]["top_cold_sinks_flood_then_empty_m3_s"] = sum(
            t["flood_then_empty_m3_s"] for t in result["steps"][-1]["top_cold_sinks"]
        )

    delivered = spill_delivery(f)
    top = np.argsort(-np.where(cold_sink, inflow, -1.0))[:5]
    top_sinks, prev_sources = [], []
    for t in top.tolist():
        if not cold_sink[t]:
            continue
        sources = delivered.get(t, [])
        sink = {
            "inflow_m3_s": float(inflow[t]),
            "spill_delivered_m3_s": float(f.discharge_m3_s[sources].sum()) if sources else 0.0,
            "flood_then_empty_m3_s": 0.0,
            "lake_depth_m": float(f.lake_depth[t]),
        }
        top_sinks.append(sink)
        prev_sources += [
            {"sink": sink, "xyz": f.points[s].tolist(), "through_m3_s": float(f.discharge_m3_s[s]), "was_wet": bool(wet[s])}
            for s in sources
        ]

    drains = [lake.sink_node_idx for lake in leaves(f.lake_forest or []) if lake.drains_at_floor]
    drain_xyz = f.points[drains] if drains else np.zeros((0, 3))
    toggles = None
    if prev_drain_xyz is not None:
        on_now = np.zeros(len(f.flow_target), dtype=bool)
        on_now[drains] = True
        idx, ok = matched(tree, prev_drain_xyz)
        switched_off = int(np.count_nonzero(~ok | ~on_now[idx]))
        prev_tree = cKDTree(prev_drain_xyz) if len(prev_drain_xyz) else None
        if prev_tree is None:
            switched_on = len(drains)
        else:
            dist, _ = prev_tree.query(drain_xyz) if len(drains) else (np.zeros(0), None)
            switched_on = int(np.count_nonzero(dist > MATCH))
        toggles = switched_off + switched_on
    prev_drain_xyz = drain_xyz

    notch = f.breach_notch_m if len(f.breach_notch_m) == len(f.flow_target) else np.zeros(len(f.flow_target))
    rec = {
        "step": step, "elapsed_years": world.elapsed_years, "runtime_s": runtime_s,
        "peak_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / RSS_PER_MB,
        "deleted_by_freezing_m3_s": deleted_m3_s(f, args.step_years),
        "cold_sink_count": int(np.count_nonzero(cold_sink)),
        "cold_sink_inflow_m3_s": float(inflow[cold_sink].sum()),
        "cold_spilling_sinks": int(np.count_nonzero(f.should_spill & cold & land)),
        "ocean_discharge_m3_s": float(f.discharge_m3_s[f.is_ocean].sum()),
        "lake_volume_km3": float((f.lake_depth * area).sum() / 1e9),
        "lake_area_km2": float(area[wet].sum() / 1e6),
        "cold_lake_volume_km3": float((f.lake_depth[cold] * area[cold]).sum() / 1e9),
        "cold_lake_area_km2": float(area[wet & cold].sum() / 1e6),
        "glacier_volume_km3": float((f.glacier_depth * area).sum() / 1e9),
        "breached_pits": int(np.count_nonzero(notch > 0.0)),
        "floor_drain_leaves": len(drains),
        "floor_drain_toggles": toggles,
        "top_cold_sinks": top_sinks,
    }
    result["steps"].append(rec)
    print(json.dumps({k: v for k, v in rec.items() if k != "top_cold_sinks"}), flush=True)
    args.out.write_text(json.dumps(result, indent=1))
