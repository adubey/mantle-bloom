#!/usr/bin/env python3
"""Issue #522: replay a saved world and measure its climate fields, on `main` or on the branch.

Replays `--myr` from the save in `--step-years` steps. At every checkpoint it computes a fresh
climate snapshot of the world (`climate.compute_climate` at the world's own climate grid) and
records, area-weighted by cos(latitude):

- mean sea-level-equivalent air temperature over land and ocean (C);
- mean humidity over land and ocean, in g/kg. A `main` backend's unitless index is converted at
  its reference point (1.0 == air over 20 C open water == 11.56 g/kg), so the two arms compare
  directly; the raw index is kept as `humidity_index` too;
- mean precipitation over land and ocean (mm/yr) and land's share of all precipitation (%);
- the same temperature, humidity and precipitation in 15-degree latitude bands, both
  hemispheres folded together (`bands`);
- Köppen group shares of land area (A tropical, B arid, C temperate, D continental, E polar);
- land fraction and the ice cap -- land under visible glacier ice (>= GLACIER_VISIBLE_DEPTH_M),
  from the hydrology cache -- and ice volume (km^3);
- step runtime and peak RSS.

`--backend` points at the backend to import, so the same script measures both trees:

    cd backend
    .venv/bin/python ../bin/debug/measure_climate_fields.py <save>.mbworld --myr 5 --out <scratch>/branch.json
    .venv/bin/python ../bin/debug/measure_climate_fields.py <save>.mbworld --myr 5 \\
        --backend ../../mantle-bloom-origin-main/backend --out <scratch>/main.json
"""

from __future__ import annotations

import numpy as np
from replay_lib import load_replay, parse_replay_args, run_replay

# The old unitless humidity index's 1.0, in g/kg (air over 20 C open water, #522).
LEGACY_HUMIDITY_INDEX_G_KG = 11.56
BAND_EDGES_DEG = np.arange(0.0, 90.1, 15.0)


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float | None:
    total = float(weights.sum())
    return None if total <= 0.0 else round(float((values * weights).sum()) / total, 3)


def main() -> None:
    args = parse_replay_args(__doc__)
    world, result = load_replay(args, issue=522)
    from app import plates as plates_mod  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.hydroclimate import biomes, climate, hydrology  # noqa: E402

    land_classes = len(biomes.KOPPEN_CODES)
    group_of_biome = np.array([code[0] for code in biomes.KOPPEN_CODES] + [""] * len(biomes.PELAGIC_NAMES))

    def humidity_g_kg(fields) -> tuple[np.ndarray, np.ndarray | None]:
        if hasattr(fields, "humidity_g_kg"):
            return fields.humidity_g_kg, None
        return fields.humidity * LEGACY_HUMIDITY_INDEX_G_KG, fields.humidity

    def climate_snapshot() -> dict:
        fields = climate.compute_climate(world, *climate.grid_dimensions(world.climate_density))
        weight = np.repeat(np.cos(np.radians(fields.lat_deg))[:, None], len(fields.lon_deg), axis=1)
        ocean = fields.is_ocean
        land = ~ocean
        humidity, humidity_index = humidity_g_kg(fields)
        precip = fields.precipitation_mm
        land_precip = float((precip * weight)[land].sum())
        out = {
            "air_temp_land_c": _weighted_mean(fields.air_temperature_c[land], weight[land]),
            "air_temp_ocean_c": _weighted_mean(fields.air_temperature_c[ocean], weight[ocean]),
            "humidity_land_g_kg": _weighted_mean(humidity[land], weight[land]),
            "humidity_ocean_g_kg": _weighted_mean(humidity[ocean], weight[ocean]),
            "precip_land_mm": _weighted_mean(precip[land], weight[land]),
            "precip_ocean_mm": _weighted_mean(precip[ocean], weight[ocean]),
            "land_precip_share_pct": round(100.0 * land_precip / max(float((precip * weight).sum()), 1e-9), 2),
        }
        if humidity_index is not None:
            out["humidity_index_ocean"] = _weighted_mean(humidity_index[ocean], weight[ocean])
        abs_lat = np.abs(np.repeat(fields.lat_deg[:, None], len(fields.lon_deg), axis=1))
        bands = []
        for lo, hi in zip(BAND_EDGES_DEG[:-1], BAND_EDGES_DEG[1:]):
            in_band = (abs_lat >= lo) & (abs_lat < hi if hi < 90.0 else abs_lat <= hi)
            bl, bo = in_band & land, in_band & ocean
            bands.append({
                "band": f"{int(lo)}-{int(hi)}",
                "air_temp_land_c": _weighted_mean(fields.air_temperature_c[bl], weight[bl]),
                "sst_c": _weighted_mean(fields.ocean_temperature_c[bo], weight[bo]),
                "humidity_ocean_g_kg": _weighted_mean(humidity[bo], weight[bo]),
                "precip_land_mm": _weighted_mean(precip[bl], weight[bl]),
                "precip_ocean_mm": _weighted_mean(precip[bo], weight[bo]),
            })
        out["bands"] = bands
        land_biomes = fields.biome_ids[land]
        land_weight = weight[land]
        land_total = max(float(land_weight.sum()), 1e-9)
        groups = group_of_biome[np.clip(land_biomes, 0, len(group_of_biome) - 1)]
        out["koppen_land_pct"] = {
            g: round(100.0 * float(land_weight[(groups == g) & (land_biomes < land_classes)].sum()) / land_total, 2)
            for g in "ABCDE"
        }
        return out

    def snapshot() -> dict:
        out = {"elapsed_myr": world.elapsed_years / 1e6}
        hydro = world.hydrology_cache
        if hydro is not None:
            areas = plates_mod.collect_all_accounting_areas_m2(hydro.plates_in_order, line_spacing_rad(world.node_density))
            land = ~hydro.is_ocean
            on_land_ice = land & (hydro.glacier_depth >= hydrology.GLACIER_VISIBLE_DEPTH_M)
            out["land_pct"] = round(100.0 * float(areas[land].sum() / areas.sum()), 3)
            out["ice_cap_pct_of_land"] = round(100.0 * float(areas[on_land_ice].sum()) / max(float(areas[land].sum()), 1.0), 3)
            out["ice_volume_km3"] = round(float(np.dot(areas, hydro.glacier_depth)) / 1e9, 1)
        out.update(climate_snapshot())
        return out

    window = {"seconds": 0.0}

    def after_step(seconds: float) -> None:
        window["seconds"] += seconds

    def close_window(steps: int) -> dict:
        closed = {"step_seconds_mean": round(window["seconds"] / steps, 2)}
        window["seconds"] = 0.0
        return closed

    run_replay(world, args, result, snapshot, after_step, close_window)


if __name__ == "__main__":
    main()
