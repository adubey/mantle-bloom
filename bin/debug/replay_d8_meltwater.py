#!/usr/bin/env python3
"""Replay RFC 0015 G6 meltwater diagnostics on a saved world.

Run with a detached baseline or the D8 backend, one arm at a time:

    backend/.venv/bin/python bin/debug/replay_d8_meltwater.py <save> \
        --backend <backend> --myr 10 --step-years 250000 \
        --checkpoint-myr 1 --out <output.json>

Use --ice-age-period-years 8000000 for the report-only ice-age arm. The
baseline records null for D8 diagnostics it does not provide.
"""

from __future__ import annotations

from dataclasses import asdict

import numpy as np

from replay_lib import load_replay, parse_replay_args, run_replay


def _arguments(parser) -> None:
    parser.add_argument("--ice-age-period-years", type=float, default=0.0)


def _glacier_rivers(hydro, hydrology) -> list[float]:
    """Mouth discharges of networks with a headwater next to visible glacier ice."""
    ice = np.asarray(hydro.glacier_depth) >= hydrology.GLACIER_VISIBLE_DEPTH_M
    near_ice = ice | np.any(ice[hydro.neighbor_idx], axis=1)
    discharges = []
    for river in hydrology.group_rivers(hydro):
        members = river.member_idx
        downstream = hydro.flow_target[members]
        internal = downstream[np.isin(downstream, members)]
        headwaters = members[~np.isin(members, internal)]
        if np.any(near_ice[headwaters]):
            discharges.append(river.flow_rate)
    return sorted(discharges, reverse=True)


def _ledger(value):
    if value is None:
        return None
    return {**asdict(value), "residual_m3": value.residual_m3}


def main() -> None:
    args = parse_replay_args(__doc__, _arguments)
    world, result = load_replay(args, issue=1019)
    from app.hydroclimate import hydrology  # noqa: E402

    world.ice_age_period_years = args.ice_age_period_years
    result["ice_age_period_years"] = world.ice_age_period_years
    world.simulate_glacier_march = True
    result["simulate_glacier_march"] = world.simulate_glacier_march
    result["glacier_stop_mode"] = getattr(world, "glacier_stop_mode", "exhaust")
    result["glacier_step_years"] = getattr(world, "glacier_step_years", None)
    result["steps"] = []

    def snapshot() -> dict:
        return {"elapsed_myr": world.elapsed_years / 1e6}

    def after_step(seconds: float) -> None:
        hydro = world.hydrology_cache
        rivers = _glacier_rivers(hydro, hydrology)
        drainage = getattr(hydro, "glacier_drainage", None)
        result["steps"].append({
            "elapsed_myr": world.elapsed_years / 1e6,
            "wall_s": seconds,
            "glacier_source": _ledger(getattr(hydro, "glacier_meltwater_ledger", None)),
            "bare_snow_source": _ledger(getattr(hydro, "bare_snowmelt_ledger", None)),
            "drainage": asdict(drainage) if drainage is not None else None,
            "warm_thick_surface_melt_m3": getattr(hydro, "warm_thick_surface_melt_m3", None),
            "ice_margin_networks": len(rivers),
            "largest_ice_margin_discharges_m3_s": rivers[:10],
        })

    run_replay(world, args, result, snapshot, after_step, lambda _steps: {})


if __name__ == "__main__":
    main()
