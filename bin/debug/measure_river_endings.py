#!/usr/bin/env python3
"""Issue #503: replay a saved world and measure where its rivers end, on `main` or on the branch.

Each step, from that step's hydrology, it measures:

- water reaching no outlet (m^3/s): inflow into warm (>= 0C), ice-free land sinks
  (`flow_target == -1`) that hold no visible lake, plus inflow into `flow_target` cycles from
  outside them. Both are rivers that reach neither the sea nor a lake;
- flow cycles: their count;
- discharge into the ocean (m^3/s);
- drawn river lines that end off the water on the elevation map (`_rivers_to_draw`): drawn
  segments whose downstream node is neither drawn on, ocean, nor a visible lake;
- length-visible river networks by mouth type (ocean / lake / other);
- warm single-node ice plugs: isolated nodes with visible ice at >= 0C;
- the ice ledger's frozen river ice and snow accumulation (km^3), every other ice-ledger
  term (`ice_<term>_km3`), and visible lake nodes.

Per-step values are averaged over each checkpoint window, beside the end-of-window glacier
volume, step runtime and peak RSS. `--backend` points at the backend to import, so the same
script measures both trees:

    cd backend
    .venv/bin/python ../bin/debug/measure_river_endings.py <local-save> --myr 5 --out <scratch>/branch.json
    .venv/bin/python ../bin/debug/measure_river_endings.py <local-save> --myr 5 \\
        --backend <baseline-worktree>/backend --out <scratch>/main.json
"""

from __future__ import annotations

import numpy as np
from replay_lib import load_replay, parse_replay_args, run_replay


def flow_cycle_mask(flow_target: np.ndarray) -> np.ndarray:
    """Nodes on a `flow_target` cycle: whatever doesn't reach -1 after `n` hops lands on one;
    walking each landing node's cycle marks the rest."""
    n = len(flow_target)
    on_cycle = np.zeros(n, dtype=bool)
    jump = np.append(np.where(flow_target >= 0, flow_target, n), n)
    hops = 1
    while hops < n:
        jump = jump[jump]
        hops *= 2
    for start in np.unique(jump[:n][jump[:n] < n]).tolist():
        node = start
        while not on_cycle[node]:
            on_cycle[node] = True
            node = int(flow_target[node])
    return on_cycle


def main() -> None:
    args = parse_replay_args(__doc__)
    world, result = load_replay(args, issue=503)
    from app import render_image  # noqa: E402
    from app.hydroclimate import hydrology  # noqa: E402
    from app.hydroclimate.erosion import climate_grid_indices  # noqa: E402

    def measure() -> dict:
        hydro = world.hydrology_cache
        flow_target = hydro.flow_target
        discharge = hydro.discharge_m3_s
        climate = world.climate_cache
        height, width = climate.precipitation_mm.shape
        row, col = climate_grid_indices(hydro.points, height, width)
        temperature = climate.air_temperature_c[row, col]
        inflow = np.zeros(len(flow_target))
        routed = flow_target >= 0
        np.add.at(inflow, flow_target[routed], discharge[routed])
        ice = hydro.glacier_depth >= hydrology.GLACIER_VISIBLE_DEPTH_M
        lake = hydro.lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
        invisible_sink = ~routed & ~hydro.is_ocean & ~ice & (temperature >= 0.0) & ~lake
        on_cycle = flow_cycle_mask(flow_target)
        into_cycle = routed & ~on_cycle
        into_cycle[into_cycle] = on_cycle[flow_target[into_cycle]]
        sink_inflow = float(inflow[invisible_sink].sum())
        cycle_inflow = float(discharge[into_cycle].sum())

        src, _, _ = render_image._rivers_to_draw(world, hide_minor_networks=True)
        drawn = np.zeros(len(flow_target), dtype=bool)
        drawn[src] = True
        ends = flow_target[src]
        off_water = ~drawn[ends] & ~hydro.is_ocean[ends] & ~lake[ends]

        visible = [r for r in hydrology.group_rivers(hydro) if render_image._river_network_is_visible(r)]
        isolated = ice & ~np.any(ice[hydro.neighbor_idx], axis=1)
        ledger = hydro.ice.ledger
        return {
            "water_to_nowhere_m3_s": sink_inflow + cycle_inflow,
            "invisible_sink_inflow_m3_s": sink_inflow,
            "cycle_inflow_m3_s": cycle_inflow,
            "flow_cycles": sum(1 for _ in _cycle_starts(flow_target, on_cycle)),
            "ocean_discharge_m3_s": float(discharge[hydro.is_ocean].sum()),
            "drawn_ends_off_water": int(off_water.sum()),
            "drawn_ends_total": int((~drawn[ends]).sum()),
            "visible_networks_ocean": sum(r.mouth_type == "ocean" for r in visible),
            "visible_networks_lake": sum(r.mouth_type == "lake" for r in visible),
            "visible_networks_other": sum(r.mouth_type == "other" for r in visible),
            "warm_ice_plugs": int((isolated & (temperature >= 0.0)).sum()),
            "frozen_river_ice_km3": ledger.frozen_river_ice_m3 / 1e9,
            "snow_accumulation_km3": ledger.accumulation_m3 / 1e9,
            "visible_lake_nodes": int(lake.sum()),
            **{f"ice_{key}_km3": value / 1e9 for key, value in vars(ledger).items() if key.endswith("_m3")},
        }

    def snapshot() -> dict:
        hydro = world.hydrology_cache
        out = {"elapsed_myr": world.elapsed_years / 1e6}
        if hydro is not None:
            out["glacier_volume_km3"] = hydro.ice.ledger.closing_m3 / 1e9 if hydro.ice is not None else None
            out["at_checkpoint"] = measure()
        return out

    window: dict = {"seconds": 0.0, "steps": []}

    def after_step(seconds: float) -> None:
        window["seconds"] += seconds
        window["steps"].append(measure())

    def close_window(steps: int) -> dict:
        rows = window["steps"]
        closed = {
            "window_mean": {key: float(np.mean([r[key] for r in rows])) for key in rows[0]},
            "window_std": {key: float(np.std([r[key] for r in rows])) for key in rows[0]},
            "step_seconds_mean": window["seconds"] / max(len(rows), 1),
        }
        window["seconds"] = 0.0
        window["steps"] = []
        return closed

    run_replay(world, args, result, snapshot, after_step, close_window)


def _cycle_starts(flow_target: np.ndarray, on_cycle: np.ndarray):
    """One node per cycle in `on_cycle`."""
    seen = np.zeros(len(flow_target), dtype=bool)
    for start in np.flatnonzero(on_cycle).tolist():
        if seen[start]:
            continue
        node = start
        while not seen[node]:
            seen[node] = True
            node = int(flow_target[node])
        yield start


if __name__ == "__main__":
    main()
