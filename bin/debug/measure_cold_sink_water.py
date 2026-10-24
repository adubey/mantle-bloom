#!/usr/bin/env python3
"""Issue #539 (RFC 0014 D2(b)): replay a saved world and measure liquid water at cold nodes.

Steps the save for `--myr` in `--step-years` steps and records, per step:

- liquid water deleted by freezing (m^3/s): the ice ledger's frozen river ice and frozen lake
  water, converted back to water at `hydrology.WATER_TO_ICE`. From #539 on the ledger never
  books any, so this is 0 by construction there; the next measure is the one that can show a
  loss on either backend;
- unbooked lake inflow (m^3/s): for every lake body resolved at its own balance this step (each
  `lakes._water_balance` call whose level was written out, not one discarded by a split or
  overwritten by a merge), the inflow that its evaporation at that level doesn't remove, summed
  over bodies that neither spill, drain at their floor nor sit at a depth cap. A lake at its
  balance has none; a lake dried by a frozen member (before #539) leaves all its inflow
  unbooked. Bodies held at a depth cap are summed separately (`held_at_depth_cap_m3_s`): the
  cap holds that water back, which is RFC 0014 D4's question, not this issue's. Before #539
  water settling at a frozen node is already in the ledger term, so on such a backend the
  inflow at frozen nodes (annual mean below 0 C, as hydrology reads it) is left out here;
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

The command line and save loading are `replay_lib`'s, so `--backend` measures `origin/main`
from a detached worktree. Ten 100 kyr steps, branch then baseline:

    cd backend
    .venv/bin/python ../bin/debug/measure_cold_sink_water.py <save> --myr 1 --out <scratch>/branch.json
    .venv/bin/python ../bin/debug/measure_cold_sink_water.py <save> --myr 1 \\
        --backend <baseline-worktree>/backend --out <scratch>/main.json
"""

from __future__ import annotations

import json
import resource
import sys
import time

import numpy as np
from replay_lib import load_replay, parse_replay_args

args = parse_replay_args(__doc__, lambda ap: ap.add_argument("--match-km", type=float, default=20.0))
world, result = load_replay(args, issue=539)

from scipy.spatial import cKDTree  # noqa: E402

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


import inspect  # noqa: E402

# Before #539 `_compute_flow_direction` took the freeze override, and frozen nodes deleted the
# water settling there through the ledger.
FREEZES_AT_COLD_NODES = "apply_freeze" in inspect.signature(hydrology._compute_flow_direction).parameters
hydrology_temperature: dict = {}
_orig_compute_hydrology = hydrology.compute_hydrology


def _record_temperature(world_, precipitation, temperature, *a, **k):
    hydrology_temperature["c"] = np.asarray(temperature)
    return _orig_compute_hydrology(world_, precipitation, temperature, *a, **k)


hydrology.compute_hydrology = _record_temperature

# Every `_water_balance` call this step: (lake, level, tier cap depth, climate, hypsometry).
# `_resolve` passes them positionally, with climate and hypsometry last on either backend.
balance_calls: list[tuple] = []
_orig_water_balance = lakes._water_balance


def _record_balance(*a):
    level = _orig_water_balance(*a)
    balance_calls.append((a[0], level, a[-5], a[-2], a[-1]))
    return level


lakes._water_balance = _record_balance


def lake_balance(f):
    """(unbooked inflow, inflow held at a depth cap), m^3/s, over this step's final bodies."""
    unbooked = held = 0.0
    frozen = hydrology_temperature["c"] < hydrology.FREEZE_POINT_C if FREEZES_AT_COLD_NODES else None
    for lake, level, tier_depth, climate, hyps in balance_calls:
        m = lake.members
        if lake.children and level < lake.min_depth:
            continue  # a split: discarded, each child resolved on its own
        if not np.allclose(f.lake_depth[m], hyps.mean_depth(level, m), rtol=1e-9, atol=1e-9):
            continue  # overwritten by a merge or a parent's resolution
        if lake.drains_at_floor or (lake.max_depth is not None and level >= lake.max_depth - 1e-6):
            continue  # its water passes on through the outlet
        inflow = climate.inflow_m3_per_yr[m] if frozen is None else np.where(frozen[m], 0.0, climate.inflow_m3_per_yr[m])
        surplus = float(inflow.sum()) - lakes.lake_evaporation_m3_per_yr(level, m, hyps, climate)
        rim_depth = None if lake.max_depth is None else lake.max_depth - lake.floor_elevation
        at_cap = (rim_depth is None or tier_depth < rim_depth) and level >= lake.floor_elevation + tier_depth - 1e-6
        if at_cap:
            held += max(surplus, 0.0)
        else:
            unbooked += max(surplus, 0.0)
    return unbooked / SECONDS_PER_YEAR, held / SECONDS_PER_YEAR


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


result.update(elapsed_years=world.elapsed_years, nodes=len(world.hydrology_cache.elevation), steps=[])
del result["checkpoints"]  # this script records every step instead
print(json.dumps({k: v for k, v in result.items() if k != "steps"}), flush=True)
prev_sources: list[dict] = []  # last step's top-sink spill sources, judged once this step is known
prev_drain_xyz = None

for step in range(1, int(round(args.myr * 1e6 / args.step_years)) + 1):
    start = time.perf_counter()
    balance_calls.clear()
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
        **dict(zip(("unbooked_lake_inflow_m3_s", "held_at_depth_cap_m3_s"), lake_balance(f))),
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
