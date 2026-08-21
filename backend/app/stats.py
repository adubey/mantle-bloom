"""Aggregate world statistics for the frontend's Stats panel.

Same philosophy as climate.py: every field here is derived from the *current* world state
on every call, via climate.py's fixed equirectangular grid (so land/ocean, temperature, and
precipitation stats all agree with what the climate map views actually display). Uses
`climate.compute_climate_cached` rather than `compute_climate` directly, reusing whatever
erosion.py already computed this step instead of triggering a second recomputation -- see
that function's own docstring for what "cached" means here (same-turn reuse, up to one step
stale, not a correctness mechanism). History across time is a frontend concern (see
App.tsx) -- the backend has no per-step storage of its own, only a snapshot endpoint.

Land vs ocean, and land/ocean fractions, use climate.py's own `is_ocean` mask (elevation <=
world.sea_level_m, live-adjustable -- see World.sea_level_m) rather than crust_type, for the
same reason climate.py itself does: a submerged continental shelf is physically ocean.
Fractions are a plain count over grid cells, not cos(lat)-weighted -- the climate grid is a
plain equirectangular lattice, not an equal-area projection, so this is an approximation.

`elevation_*_m` covers land cells only (height above sea level); `ocean_depth_*_m` is the
mirror for ocean cells (positive depth below sea level, i.e. `sea_level_m - elevation`) --
kept as two separate stats rather than one combined min/max/mean the way it used to be,
since lumping a -11000m trench and a 9000m peak into the same range made neither number
very informative on its own.

`biome_land_fraction` reuses biomes.classify_biomes the same way the "biome" map view itself
does (see render_image.py's `_render_biome_view`: land cells use air_temperature_c, ocean
cells are forced to Ocean regardless), though against climate.py's native (coarser) grid
here rather than that view's own finer render grid -- an aggregate land fraction doesn't need
the extra resolution the way a rendered map's coastlines visibly do --
the denominator is land cells only (Ocean is always 0% by construction, so it's omitted from
the dict entirely rather than reported as a permanent zero).
"""

from __future__ import annotations

import numpy as np

from . import biomes, climate
from .render_image import grid_slope
from .world import World


def _min_max_mean(values: np.ndarray) -> tuple[float | None, float | None, float | None]:
    if values.size == 0:
        return None, None, None
    return float(values.min()), float(values.max()), float(values.mean())


def compute_stats(world: World) -> dict:
    fields = climate.compute_climate_cached(world)
    is_ocean = fields.is_ocean
    is_land = ~is_ocean
    total = is_ocean.size

    elevation_min, elevation_max, elevation_mean = _min_max_mean(fields.elevation_m[is_land])
    ocean_depth = world.sea_level_m - fields.elevation_m[is_ocean]
    ocean_depth_min, ocean_depth_max, ocean_depth_mean = _min_max_mean(ocean_depth)
    land_temp_min, land_temp_max, land_temp_mean = _min_max_mean(fields.land_temperature_c[is_land])
    air_temp_min, air_temp_max, air_temp_mean = _min_max_mean(fields.air_temperature_c[is_land])
    ocean_temp_min, ocean_temp_max, ocean_temp_mean = _min_max_mean(fields.ocean_temperature_c[is_ocean])
    precip_min, precip_max, precip_mean = _min_max_mean(fields.precipitation_mm)

    display_temp = np.where(fields.is_ocean, fields.ocean_temperature_c, fields.air_temperature_c)
    slope = grid_slope(fields.elevation_m, fields.lat_deg)
    biome_ids = biomes.classify_biomes(display_temp, fields.precipitation_mm, fields.elevation_m, slope, fields.is_ocean, world.sea_level_m)
    land_biome_ids = biome_ids[is_land]
    n_land = int(is_land.sum())
    biome_land_fraction = {
        name: float(np.count_nonzero(land_biome_ids == i)) / n_land
        for i, name in enumerate(biomes.BIOME_NAMES)
        if i != biomes.OCEAN and n_land > 0
    }

    return {
        "elapsed_years": world.elapsed_years,
        "land_fraction": float(is_land.sum()) / total,
        "ocean_fraction": float(is_ocean.sum()) / total,
        "elevation_min_m": elevation_min,
        "elevation_max_m": elevation_max,
        "elevation_mean_m": elevation_mean,
        "ocean_depth_min_m": ocean_depth_min,
        "ocean_depth_max_m": ocean_depth_max,
        "ocean_depth_mean_m": ocean_depth_mean,
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
        "biome_land_fraction": biome_land_fraction,
    }
