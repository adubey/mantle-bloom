#!/usr/bin/env python3
"""Report RFC 0015's global ice scale, eustatic equivalent and sliding speeds.

Run the same physical-rate 10 Myr replay used by RFC 0015's diagnosis/held-out
work, then measure global ice and ocean area plus area-weighted sliding-speed
quantiles for interior and margin ice at each 1 Myr checkpoint:

    backend/.venv/bin/python bin/debug/measure_rfc0015_scale.py <save> \
        --myr 10 --step-years 1000000 --checkpoint-myr 1 --out <output.json>

The flat-ocean sea-level equivalent is grounded ice volume converted to water
equivalent and divided by the connected ocean area at that checkpoint. It is a
scale comparison with the LGM reference, not the world's solved sea level.
Velocities are the model's calibrated sliding-speed estimate from mass-balance
flux, not velocities from a time-integrated ice-flow solver.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from replay_lib import load_replay, parse_replay_args, run_replay


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantiles: tuple[float, ...]) -> list[float | None]:
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(valid):
        return [None] * len(quantiles)
    values = values[valid]
    weights = weights[valid]
    order = np.argsort(values)
    values = values[order]
    cumulative = np.cumsum(weights[order])
    targets = np.asarray(quantiles) * cumulative[-1]
    return [float(values[min(np.searchsorted(cumulative, target, side="left"), len(values) - 1)]) for target in targets]


def speed_summary(mask: np.ndarray, speed: np.ndarray, areas: np.ndarray) -> dict:
    quantiles = weighted_quantile(speed[mask], areas[mask], (0.5, 0.9, 0.95))
    values = speed[mask]
    values = values[np.isfinite(values)]
    return {
        "node_count": int(np.count_nonzero(mask)),
        "area_km2": float(areas[mask].sum() / 1e6),
        "area_weighted_speed_m_per_yr": {
            "p50": quantiles[0],
            "p90": quantiles[1],
            "p95": quantiles[2],
            "max": float(values.max()) if values.size else None,
        },
    }


def snapshot(world, plates, ice_sheet, hydrology) -> dict:
    hydro = world.hydrology_cache
    if hydro is None:
        return {"elapsed_myr": world.elapsed_years / 1e6, "state_available": False}

    from app.elevation_lines import line_spacing_rad  # noqa: E402

    areas = plates.collect_all_accounting_areas_m2(hydro.plates_in_order, line_spacing_rad(world.node_density))
    depth = np.asarray(hydro.glacier_depth, dtype=float)
    ocean = np.asarray(hydro.is_ocean, dtype=bool)
    ice = hydro.ice
    speed = np.asarray(ice.sliding_speed_m_per_yr if ice is not None else np.zeros_like(depth), dtype=float)
    deflection = plates.collect_all_ice_load_deflection(hydro.plates_in_order)
    unloaded_bed = np.asarray(hydro.elevation, dtype=float) - np.asarray(deflection, dtype=float)
    grounded = ice_sheet.is_grounded(depth, unloaded_bed, world.sea_level_m, ocean)

    ice_volume_m3 = float(np.dot(areas, depth))
    grounded_ice_volume_m3 = float(np.dot(areas, np.where(grounded, depth, 0.0)))
    ocean_area_m2 = float(areas[ocean].sum())
    ice_bearing = depth >= hydrology.GLACIER_VISIBLE_DEPTH_M
    neighbors = np.asarray(hydro.neighbor_idx, dtype=np.int64)
    touches_ice_free = np.any(~ice_bearing[neighbors], axis=1)
    margin = ice_bearing & touches_ice_free
    interior = ice_bearing & ~touches_ice_free
    rho_ice_water = float(ice_sheet.RHO_ICE / ice_sheet.RHO_WATER)
    eustatic_equivalent_m = (grounded_ice_volume_m3 * rho_ice_water / ocean_area_m2) if ocean_area_m2 > 0 else None
    lgm_volume_fraction = grounded_ice_volume_m3 / (52e6 * 1e9)

    return {
        "elapsed_myr": world.elapsed_years / 1e6,
        "state_available": True,
        "sea_level_m": float(world.sea_level_m),
        "ocean_area_km2": ocean_area_m2 / 1e6,
        "ice_volume_km3": ice_volume_m3 / 1e9,
        "grounded_ice_volume_km3": grounded_ice_volume_m3 / 1e9,
        "grounded_ice_as_share_of_lgm_52m_km3": lgm_volume_fraction,
        "grounded_ice_flat_ocean_sea_level_equivalent_m": eustatic_equivalent_m,
        "sea_level_equivalent_as_share_of_lgm_134m": (
            eustatic_equivalent_m / 134.0 if eustatic_equivalent_m is not None else None
        ),
        "visible_ice_area_km2": float(areas[ice_bearing].sum() / 1e6),
        "visible_ice_margin_area_share": (
            float(areas[margin].sum() / areas[ice_bearing].sum()) if np.any(ice_bearing) else None
        ),
        "sliding_speed": {
            "all_visible_ice": speed_summary(ice_bearing, speed, areas),
            "interior": speed_summary(interior, speed, areas),
            "margin": speed_summary(margin, speed, areas),
        },
    }


def main() -> None:
    args = parse_replay_args(__doc__)
    world, result = load_replay(args, issue=544)
    from app import plates  # noqa: E402
    from app.hydroclimate import hydrology, ice_sheet  # noqa: E402

    # Keep machine-specific paths out of checked-in compact reports.
    result["backend"] = "backend"
    result["save"] = Path(result["save"]).name
    world.ice_age_period_years = 0.0
    world.simulate_glacier_march = True
    result.update({
        "purpose": "RFC 0015 scale check",
        "glacier_solver_backend": "rust" if ice_sheet._glacier_rust.RUST_AVAILABLE else "python",
        "ice_age_period_years": world.ice_age_period_years,
        "simulate_glacier_march": world.simulate_glacier_march,
        "glacier_stop_mode": getattr(world, "glacier_stop_mode", "exhaust"),
        "glacier_step_years": getattr(world, "glacier_step_years", None),
        "lgm_grounded_ice_reference_km3": 52e6,
        "lgm_sea_level_reference_m": 134.0,
        "checkpoints": [],
    })

    def read_snapshot() -> dict:
        return snapshot(world, plates, ice_sheet, hydrology)

    run_replay(world, args, result, read_snapshot, lambda _seconds: None, lambda _steps: {})


if __name__ == "__main__":
    main()
