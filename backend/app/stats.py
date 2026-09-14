"""Aggregate world statistics for the frontend's Stats panel.

Same philosophy as climate.py: every field here is derived from the *current* world state
on every call, via climate.py's fixed equirectangular grid (so land/ocean, temperature, and
precipitation stats all agree with what the climate map views actually display). Uses
`climate.compute_climate_cached` rather than `compute_climate` directly, reusing whatever
erosion.py already computed this step instead of triggering a second recomputation -- see
that function's own docstring for what "cached" means here (same-turn reuse, up to one step
stale, not a correctness mechanism). `compute_stats` itself is a stateless, single-instant
snapshot -- `GET /world/stats`; the recorded time series that snapshot feeds is `World.
stats_history` (see world.py's `World.record_stats`/`GET /world/stats_history`), which is
what actually persists across a save/load, not this module.

Land vs ocean, and land/ocean fractions, use climate.py's own `is_ocean` mask (below
world.sea_level_m, live-adjustable -- see World.sea_level_m -- and connected to the world
ocean, see hydrology.connected_ocean_mask) rather than crust_type, for the same reason
climate.py itself does: a submerged continental shelf is physically ocean, an enclosed
interior depression is not.

A landlocked lake or sea isn't part of that connected ocean (that's exactly what makes it a
lake, not a bay), but it's still standing water, not land -- `land_fraction`/`ocean_fraction`
(labelled "Land"/"Water" in the frontend's Stats panel, see StatsModal.tsx) fold in
`fields.lake_depth_m` (the same resampled, persisted per-node depth hydrology.py's own
`is_lake` check uses, cleared past `hydrology.LAKE_MIN_VISIBLE_DEPTH_M` the same way) so a
basin covered by a real, visible lake counts toward "Water" the instant it floods, at any
depth -- not just once it happens to sink deep enough to trip `_reconcile_land_ocean`'s own
stale-cache fallback below. Every *other* stat here (`elevation_*`, `ocean_depth_*`, the
temperature stats, `biome_land_fraction`/`biome_ocean_fraction`) deliberately keeps reading
the narrower `is_land`/`is_ocean` split unchanged -- a lake sitting over dry-land elevation
isn't ocean bathymetry, and biome classification has no separate "lake" class to route it to.
Fractions are a plain count over grid cells, not cos(lat)-weighted -- the climate grid is a
plain equirectangular lattice, not an equal-area projection, so this is an approximation.

`elevation_*_m` covers land cells only (height above sea level); `ocean_depth_*_m` is the
mirror for ocean cells (positive depth below sea level, i.e. `sea_level_m - elevation`) --
kept as two separate stats rather than one combined min/max/mean the way it used to be,
since lumping a -11000m trench and a 9000m peak into the same range made neither number
very informative on its own.

`elevation_*_m`/`ocean_depth_*_m` reconcile `fields.is_ocean` against the *current* elevation
before splitting land from ocean, rather than trusting it outright: `is_ocean` is resampled
from last step's hydrology cache (see hydrology.sample_is_ocean's docstring), so mid-step
tectonics (a collision uplifting a former seabed, a rift dropping crust out from under an
existing coastline) can leave it briefly disagreeing with where elevation actually sits versus
world.sea_level_m -- a mountain that just rose still reads "ocean" until hydrology's next
recompute, and the reverse for a newly-subsided trench. Left alone this shows up as a
physically nonsensical negative "ocean depth" (a peak above sea level tagged ocean) or an
elevation "land" value plunging to MIN_ELEVATION_M (a trench tagged land, far past any real
endorheic basin -- Earth's deepest, the Dead Sea, bottoms out around -430 m, hence
`MAX_ENDORHEIC_BASIN_DEPTH_M`). `_reconcile_land_ocean` folds both mismatches the physically
sane way: an "ocean" cell sitting above sea level is recounted as land (it plainly is, whatever
last step's cache says), and a "land" cell sitting deeper than any plausible desert basin is
recounted as ocean instead. A genuine, shallow endorheic basin -- still tagged land, still
below sea level, but within that plausible band -- is untouched, so elevation_min_m can
legitimately land slightly below 0 rather than being hard-floored there.

`biome_land_fraction` reads `ClimateFields.biome_ids` -- the same stored classification
`compute_climate` computes once (via `biomes.smooth_biome_field`, so the boundary-cleanup
pass is already baked in) and every other biome-consuming caller now shares, at climate.py's
native (coarser) grid rather than the "biome" map view's own finer render grid
(render_image.py's `_render_biome_view`) -- an aggregate land fraction doesn't need the extra
resolution the way a rendered map's coastlines visibly do. The denominator is land cells only
(Ocean is always 0% by construction, so it's omitted from the dict entirely rather than
reported as a permanent zero).

`biome_ocean_fraction` is the exact mirror for the pelagic (ocean) classes -- the same
`ClimateFields.biome_ids` field, but counted over ocean cells with `biomes.OCEAN_IDS` as the
denominator, so the Stats panel's Biome tab can chart the ocean provinces over time the same
way it charts the Köppen land classes. Land Köppen classes are omitted from it (0% of ocean),
and it's `{}` when there are no ocean cells at all.

`plate_count`/`elevation_point_count` are two of the exceptions to "every stat here is a
spatial min/max/mean snapshot of the current world": each is a single running total (plate
count, and the sum of every plate's own `node_count()`), with no per-call distribution to
take a min/max/mean of. The frontend's Simulation tab is what turns a run of these single
numbers into a min/max/mean/std-dev over time, the same "backend snapshot, frontend
accumulates history" split every other stat here already uses.

`total_land_area_km2`/`total_continental_crust_volume_km3` are the same kind of running
total, added to give the Simulation tab a mass/volume-conservation check independent of
land_fraction -- read straight off `world.plates`' own Hc columns and current elevation,
*not* the climate grid `is_land`/`is_ocean` above (so, unlike every stat above, immune to
the hydrology-cache staleness `_reconcile_land_ocean`'s own docstring describes, and not an
approximation of a non-equal-area grid either -- every lattice node really is the same
physical footprint, `lithosphere.node_area_m2`). The distinction they're for: land area can
swing a lot from tectonics moving existing crust above/below sea level (isostasy, sea-level
change, redistribution within a colliding pair) while the underlying crustal *volume* barely
moves -- vs. a genuine mass-conservation bug (a topology change that drops a column's volume
instead of preserving it, e.g. an unaccounted-for retreat) actually shrinking the volume
number itself. Continental-crust-only (oceanic crust is routinely created/destroyed by
spreading/subduction by design, so summing it in wouldn't isolate a conservation bug the way
continental-only does); see docs/TODO.md's collision/land-fraction investigation notes for
the mechanisms this was added to help tell apart.
"""

from __future__ import annotations

import numpy as np

from . import biomes, climate, hydrology, lithosphere
from .elevation_lines import line_spacing_rad
from .world import World


def _min_max_mean_std(values: np.ndarray) -> tuple[float | None, float | None, float | None, float | None]:
    if values.size == 0:
        return None, None, None, None
    return float(values.min()), float(values.max()), float(values.mean()), float(values.std())


# How far below sea level an ordinary endorheic desert basin can plausibly sit -- Earth's most
# extreme example, the Dead Sea depression, bottoms out around -430 m, so 500 m leaves a little
# headroom. Below this, a "land" cell is almost certainly a stale-cache mismatch (see
# `_reconcile_land_ocean`), not a real closed basin.
MAX_ENDORHEIC_BASIN_DEPTH_M = 500.0


def _reconcile_land_ocean(fields: "climate.ClimateFields", sea_level_m: float) -> tuple[np.ndarray, np.ndarray]:
    """(is_ocean, is_land) for stats purposes, correcting `fields.is_ocean` (last-step-cached,
    see module docstring) against the *current* elevation: a cell it calls "ocean" despite
    sitting above sea level is recounted as land, and a cell it calls "land" despite sitting
    deeper than any plausible endorheic basin (`MAX_ENDORHEIC_BASIN_DEPTH_M`) is recounted as
    ocean. Everything else -- including a genuine, shallow below-sea-level closed basin --
    keeps its original classification."""
    elevation = fields.elevation_m
    is_ocean = fields.is_ocean & (elevation <= sea_level_m)
    is_ocean |= (~fields.is_ocean) & (elevation < sea_level_m - MAX_ENDORHEIC_BASIN_DEPTH_M)
    return is_ocean, ~is_ocean


def _total_land_area_and_continental_volume(world: World) -> tuple[float, float]:
    """(total land area m^2, total continental crustal volume m^3) straight off
    `world.plates` -- see this module's own docstring for why these are computed here rather
    than off the climate grid `compute_stats` otherwise uses throughout. `node_area_m2` is
    (almost exactly) constant across latitude by construction, so both are a plain node-count
    sum, not an integral -- one pass per plate, no grid resample."""
    area_m2 = lithosphere.node_area_m2(line_spacing_rad(world.node_density))
    land_area = 0.0
    continental_volume = 0.0
    for plate in world.plates:
        _, elevation = plate.all_points_and_elevation()
        if len(elevation) == 0:
            continue
        land_area += float(np.count_nonzero(elevation > world.sea_level_m)) * area_m2
        if plate.crust_type == "continental":
            continental_volume += float(np.sum(plate.collect("crustal_thickness_m"))) * area_m2
    return land_area, continental_volume


def _is_water(fields: "climate.ClimateFields", is_ocean: np.ndarray) -> np.ndarray:
    """`is_ocean` plus every cell covered by a real, visible lake or sea (see this module's
    own docstring) -- the mask `land_fraction`/`ocean_fraction` actually want, as opposed to
    the narrower `is_ocean` every other stat below still uses unchanged."""
    is_lake_or_sea = fields.lake_depth_m > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    return is_ocean | is_lake_or_sea


def compute_stats(world: World) -> dict:
    fields = climate.compute_climate_cached(world)
    is_ocean, is_land = _reconcile_land_ocean(fields, world.sea_level_m)
    is_water = _is_water(fields, is_ocean)
    total = is_ocean.size

    elevation_min, elevation_max, elevation_mean, elevation_std = _min_max_mean_std(fields.elevation_m[is_land])
    ocean_depth = world.sea_level_m - fields.elevation_m[is_ocean]
    ocean_depth_min, ocean_depth_max, ocean_depth_mean, ocean_depth_std = _min_max_mean_std(ocean_depth)
    land_temp_min, land_temp_max, land_temp_mean, land_temp_std = _min_max_mean_std(fields.land_temperature_c[is_land])
    air_temp_min, air_temp_max, air_temp_mean, air_temp_std = _min_max_mean_std(fields.air_temperature_c[is_land])
    ocean_temp_min, ocean_temp_max, ocean_temp_mean, ocean_temp_std = _min_max_mean_std(fields.ocean_temperature_c[is_ocean])
    precip_min, precip_max, precip_mean, precip_std = _min_max_mean_std(fields.precipitation_mm)

    land_area_m2, continental_crust_volume_m3 = _total_land_area_and_continental_volume(world)

    land_biome_ids = fields.biome_ids[is_land]
    n_land = int(is_land.sum())
    biome_land_fraction = {
        name: float(np.count_nonzero(land_biome_ids == i)) / n_land
        for i, name in enumerate(biomes.BIOME_NAMES)
        if i not in biomes.OCEAN_IDS and n_land > 0
    }
    ocean_biome_ids = fields.biome_ids[is_ocean]
    n_ocean = int(is_ocean.sum())
    biome_ocean_fraction = {
        name: float(np.count_nonzero(ocean_biome_ids == i)) / n_ocean
        for i, name in enumerate(biomes.BIOME_NAMES)
        if i in biomes.OCEAN_IDS and n_ocean > 0
    }

    return {
        "elapsed_years": world.elapsed_years,
        "plate_count": len(world.plates),
        "elevation_point_count": sum(p.node_count() for p in world.plates),
        "sea_level_m": world.sea_level_m,
        "total_land_area_km2": land_area_m2 / 1.0e6,
        "total_continental_crust_volume_km3": continental_crust_volume_m3 / 1.0e9,
        "land_fraction": float((~is_water).sum()) / total,
        "ocean_fraction": float(is_water.sum()) / total,
        "elevation_min_m": elevation_min,
        "elevation_max_m": elevation_max,
        "elevation_mean_m": elevation_mean,
        "elevation_std_m": elevation_std,
        "ocean_depth_min_m": ocean_depth_min,
        "ocean_depth_max_m": ocean_depth_max,
        "ocean_depth_mean_m": ocean_depth_mean,
        "ocean_depth_std_m": ocean_depth_std,
        "land_temperature_min_c": land_temp_min,
        "land_temperature_max_c": land_temp_max,
        "land_temperature_mean_c": land_temp_mean,
        "land_temperature_std_c": land_temp_std,
        "air_temperature_min_c": air_temp_min,
        "air_temperature_max_c": air_temp_max,
        "air_temperature_mean_c": air_temp_mean,
        "air_temperature_std_c": air_temp_std,
        "ocean_temperature_min_c": ocean_temp_min,
        "ocean_temperature_max_c": ocean_temp_max,
        "ocean_temperature_mean_c": ocean_temp_mean,
        "ocean_temperature_std_c": ocean_temp_std,
        "precipitation_min_mm": precip_min,
        "precipitation_max_mm": precip_max,
        "precipitation_mean_mm": precip_mean,
        "precipitation_std_mm": precip_std,
        "biome_land_fraction": biome_land_fraction,
        "biome_ocean_fraction": biome_ocean_fraction,
    }
