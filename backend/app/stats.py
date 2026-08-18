"""Aggregate world statistics for the frontend's Stats panel.

Fully stateless, same philosophy as climate.py: every field here is recomputed from scratch
from the *current* world state on every call, via climate.compute_climate's fixed
equirectangular grid (so land/ocean, temperature, and precipitation stats all agree with
what the climate map views actually display) -- no history, no caching. History across time
is a frontend concern (see App.tsx), matching plate-sim's own precedent: the backend has no
analogous per-step storage, only a stateless snapshot endpoint.

Land vs ocean, and land/ocean fractions, use climate.py's own `is_ocean` mask (elevation <=
sea level, sea level = 0.0 -- see plates.py, there's no separate named SEA_LEVEL constant)
rather than crust_type, for the same reason climate.py itself does: a submerged continental
shelf is physically ocean. Fractions are a plain count over grid cells, not cos(lat)-weighted
-- the climate grid is a plain equirectangular lattice, not an equal-area projection, so this
is an approximation, but it's the same one plate-sim's own stats endpoint makes.
"""

from __future__ import annotations

import numpy as np

from . import climate
from .world import World


def _min_max_mean(values: np.ndarray) -> tuple[float | None, float | None, float | None]:
    if values.size == 0:
        return None, None, None
    return float(values.min()), float(values.max()), float(values.mean())


def compute_stats(world: World) -> dict:
    fields = climate.compute_climate(world)
    is_ocean = fields.is_ocean
    is_land = ~is_ocean
    total = is_ocean.size

    elevation_min, elevation_max, elevation_mean = _min_max_mean(fields.elevation_m)
    land_temp_min, land_temp_max, land_temp_mean = _min_max_mean(fields.land_temperature_c[is_land])
    air_temp_min, air_temp_max, air_temp_mean = _min_max_mean(fields.air_temperature_c[is_land])
    ocean_temp_min, ocean_temp_max, ocean_temp_mean = _min_max_mean(fields.ocean_temperature_c[is_ocean])
    precip_min, precip_max, precip_mean = _min_max_mean(fields.precipitation_mm)

    return {
        "elapsed_years": world.elapsed_years,
        "land_fraction": float(is_land.sum()) / total,
        "ocean_fraction": float(is_ocean.sum()) / total,
        "elevation_min_m": elevation_min,
        "elevation_max_m": elevation_max,
        "elevation_mean_m": elevation_mean,
        "land_temperature_min_c": land_temp_min,
        "land_temperature_max_c": land_temp_max,
        "land_temperature_mean_c": land_temp_mean,
        "air_temperature_min_c": air_temp_min,
        "air_temperature_max_c": air_temp_max,
        "air_temperature_mean_c": air_temp_mean,
        "ocean_temperature_min_c": ocean_temp_min,
        "ocean_temperature_max_c": ocean_temp_max,
        "ocean_temperature_mean_c": ocean_temp_mean,
        "precipitation_min_mm": precip_min,
        "precipitation_max_mm": precip_max,
        "precipitation_mean_mm": precip_mean,
    }
