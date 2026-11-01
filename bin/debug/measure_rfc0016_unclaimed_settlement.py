#!/usr/bin/env python3
"""RFC 0016 G0: classify where routed river water settles on land in one hydrology step.

The lake solver only books supply at the members of a lake body it advances, so water settling
anywhere else leaves the rain-only land ledger as its residual. This steps a saved world once
and sorts each reservoir pass's land settlement by the kind of node it landed on.

Example:

    backend/.venv/bin/python bin/debug/measure_rfc0016_unclaimed_settlement.py <checkpoint.mbworld> \\
        --myr 0.1 --out <scratch>/rfc0016-settlement.json
"""

from __future__ import annotations

import json
from collections import defaultdict

import numpy as np

from replay_lib import load_replay, parse_replay_args


def node_classes(forest, n: int) -> tuple[np.ndarray, list[str]]:
    from app.hydroclimate import lakes

    names = [
        "no_lake_body", "body_member", "leaf_sink", "floor_drain_member", "floor_drain_sink",
    ]
    classes = np.zeros(n, dtype=np.int8)
    bodies = list(lakes.iter_all_lakes(forest))
    for body in bodies:
        if body.drains_at_floor:
            continue
        classes[np.asarray(body.members)] = 1
    for body in bodies:
        if body.children:
            continue
        if body.drains_at_floor:
            classes[np.asarray(body.members)] = 3
            classes[body.sink_node_idx] = 4
        else:
            classes[body.sink_node_idx] = 2
    return classes, names


def main() -> None:
    args = parse_replay_args(__doc__)
    world, result = load_replay(args, issue=725)
    world.ice_age_period_years = 0.0
    from app import world as world_mod  # noqa: E402
    from app.hydroclimate import hydrology, lakes  # noqa: E402

    routes = []
    captured = {}
    original_route = hydrology.route_river_reservoirs
    original_solve = lakes.solve_lakes_prognostic

    def route(*a, **k):
        routes.append((a[3], original_route(*a, **k), a[2]))
        return routes[-1][1]

    def solve(forest, *a, **k):
        captured["forest"] = forest
        return original_solve(forest, *a, **k)

    hydrology.route_river_reservoirs = route
    lakes.solve_lakes_prognostic = solve
    world_mod.step_world(world, args.step_years)

    n = len(routes[0][0])
    classes, names = node_classes(captured["forest"], n)
    for name, (targets, routed, is_ocean) in zip(("runoff", "overflow"), routes):
        settled = np.where(is_ocean, 0.0, routed.settled_m3)
        totals = defaultdict(lambda: {"m3": 0.0, "nodes": 0})
        for index, class_name in enumerate(names):
            mask = (classes == index) & (settled > 0.0)
            totals[class_name] = {"m3": float(settled[mask].sum()), "nodes": int(mask.sum())}
        sinks = (classes == 0) & (settled > 0.0)
        totals["no_lake_body_with_target"] = {
            "m3": float(settled[sinks & (targets >= 0)].sum()), "nodes": int((sinks & (targets >= 0)).sum()),
        }
        result[name] = dict(totals)
    hydro = world.hydrology_cache
    # Where the largest stranded overflow settles, and the path that ends there.
    targets, routed, is_ocean = routes[1]
    settled = np.where(is_ocean, 0.0, routed.settled_m3)
    bodies = list(lakes.iter_all_lakes(captured["forest"]))
    stranded = []
    for node in np.argsort(-settled)[:5].tolist():
        if settled[node] <= 0.0:
            continue
        holders = [body for body in bodies if node in set(np.asarray(body.members).tolist())]
        hop, seen = node, [node]
        for _ in range(6):
            hop = int(targets[hop])
            if hop < 0 or hop in seen:
                seen.append(hop)
                break
            seen.append(hop)
        stranded.append({
            "node": node, "settled_m3": float(settled[node]), "target": int(targets[node]),
            "path": seen, "class": names[int(classes[node])],
            "elevation_m": float(hydro.elevation[node]), "lake_depth_m": float(hydro.lake_depth[node]),
            "should_spill": bool(hydro.should_spill[node]) if hasattr(hydro, "should_spill") else None,
            "bodies": [
                {
                    "members": len(body.members), "sink": int(body.sink_node_idx),
                    "is_spilling": bool(body.is_spilling), "drains_at_floor": bool(body.drains_at_floor),
                    "outlet_target": int(body.outlet_target_idx), "children": len(body.children),
                }
                for body in holders if not body.children
            ],
        })
    result["stranded_overflow"] = stranded
    overflow = routes[1][1].ledger
    result["overflow_input_m3"] = overflow.input_m3
    result["lake_ledger_dry_runoff_m3"] = hydro.lake_water_ledger.dry_runoff_m3
    result["lake_ledger_river_storage_m3"] = hydro.lake_water_ledger.river_storage_m3
    result["lake_ledger_meltwater_m3"] = hydro.lake_water_ledger.meltwater_m3
    result["land_residual_m3"] = hydro.land_water_ledger.residual_m3
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=1))
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
