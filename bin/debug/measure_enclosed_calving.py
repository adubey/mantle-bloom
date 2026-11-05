#!/usr/bin/env python3
"""Issue #552: how much calved ice leaves from nodes that are not connected ocean?

Replays a save and records, every step: total calved ice (m^3 of ice), the part calved at nodes
where `is_ocean` is False (an enclosed basin below sea level calves its thin ice straight to
the ocean on `origin/main`; the fix grounds it), that part's share, the glacier ice volume,
sea level and the worst ice-ledger residual. The same script measures both trees:

    backend/.venv/bin/python bin/debug/measure_enclosed_calving.py <save>.mbworld \\
        --myr 2 --step-years 100000 --checkpoint-myr 2 --out <scratch>/branch.json
    backend/.venv/bin/python bin/debug/measure_enclosed_calving.py <save>.mbworld \\
        --myr 2 --step-years 100000 --checkpoint-myr 2 --backend <baseline>/backend --out <scratch>/main.json
"""

from __future__ import annotations

import numpy as np

from replay_lib import load_replay, parse_replay_args, run_replay


def main() -> None:
    args = parse_replay_args(__doc__)
    world, result = load_replay(args, issue=552)
    window: list[dict] = []

    def snapshot() -> dict:
        # A freshly loaded world has no ice step yet, so its ice volume is unknown until step 1.
        ice = world.hydrology_cache.ice if world.hydrology_cache is not None else None
        return {"elapsed_myr": world.elapsed_years / 1e6, "sea_level_m": float(world.sea_level_m),
                "ice_volume_m3": float(ice.ledger.closing_m3) if ice is not None else None}

    def after_step(_seconds: float) -> None:
        hydro = world.hydrology_cache
        calved = np.asarray(hydro.ice.calved_m3, dtype=float)
        off = calved[~hydro.is_ocean]
        window.append({"calved_m3": float(calved.sum()), "calved_off_ocean_m3": float(off.sum()),
                       "off_ocean_nodes": int(np.count_nonzero(off)), "residual_m3": float(hydro.ice.ledger.residual_m3)})

    def close_window(_steps: int) -> dict:
        calved = sum(w["calved_m3"] for w in window)
        off = sum(w["calved_off_ocean_m3"] for w in window)
        out = {"calved_m3": calved, "calved_off_ocean_m3": off, "calved_off_ocean_share": off / calved if calved else 0.0,
               "max_off_ocean_nodes": max((w["off_ocean_nodes"] for w in window), default=0),
               "max_abs_ice_residual_m3": max((abs(w["residual_m3"]) for w in window), default=0.0),
               "per_step_off_ocean_m3": [w["calved_off_ocean_m3"] for w in window]}
        window.clear()
        return out

    run_replay(world, args, result, snapshot, after_step, close_window)


if __name__ == "__main__":
    main()
