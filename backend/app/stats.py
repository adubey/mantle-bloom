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
All spatial fractions, means, and standard deviations derived from that equirectangular
climate grid use exact spherical-strip area weights. Extrema remain ordinary extrema.

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
approximation of a non-equal-area grid either -- each node is weighted by its accounting
area, `Plate.accounting_areas_m2`). The distinction they're for: land area can
swing a lot from tectonics moving existing crust above/below sea level (isostasy, sea-level
change, redistribution within a colliding pair) while the underlying crustal *volume* barely
moves -- vs. a genuine mass-conservation bug (a topology change that drops a column's volume
instead of preserving it, e.g. an unaccounted-for retreat) actually shrinking the volume
number itself. Continental-crust-only (oceanic crust is routinely created/destroyed by
spreading/subduction by design, so summing it in wouldn't isolate a conservation bug the way
continental-only does); see GitHub issues #119 and #120's collision/land-fraction investigation notes for
the mechanisms this was added to help tell apart.

`land_fraction_node`/`land_fraction_stale` address a narrower gap in the same family: `land_
fraction` above is resampled from `world.hydrology_cache` (see `hydrology.sample_is_ocean`),
which is populated only inside `erosion.apply_erosion` -- so with `simulate_climate_biomes`
off, it freezes at whatever the world looked like when erosion last ran and silently stops
tracking tectonics/eustasy moving the coastline underneath it. `land_fraction_node` is a raw,
always-fresh `elevation > sea_level_m` node count (like `total_land_area_km2`, immune to the
staleness, but a *fraction* comparable to `land_fraction` rather than an area) and `land_
fraction_stale` flags the frozen-cache case explicitly rather than leaving the Stats panel to
guess. See GitHub issue #121.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import cratons, crust_transfer, lithosphere, orogeny
from .hydroclimate import biomes, climate, hydrology
from .elevation_lines import PLANET_RADIUS_KM, line_spacing_rad
from .world import World


@dataclass(frozen=True)
class Stat4:
    """A min/max/mean/std-dev quadruple for one quantity's spatial distribution over the
    current world (elevation, ocean depth, land/air/ocean temperature, precipitation, Hc/Hm)
    -- the same grouping the frontend's own `StatGroup` (see StatsModal.tsx) plots as one
    avg/min/max graph plus a std-dev band, rather than four separately-named numbers. Replaces
    what used to be four parallel tuple-unpacked locals per quantity here, and `to_dict` below
    replaces the four individually-spelled keys each quantity contributed to `compute_stats`'
    returned dict -- the dict's actual key names (and thus the frontend/save-file contract)
    are unchanged."""

    min: float | None
    max: float | None
    mean: float | None
    std: float | None

    @classmethod
    def of(cls, values: np.ndarray, weights: np.ndarray | None = None) -> "Stat4":
        if values.size == 0:
            return cls(None, None, None, None)
        if weights is None:
            mean = float(values.mean())
            variance = float(values.var())
        else:
            mean = float(np.average(values, weights=weights))
            variance = float(np.average((values - mean) ** 2, weights=weights))
        return cls(float(values.min()), float(values.max()), mean, float(np.sqrt(variance)))

    def to_dict(self, prefix: str, suffix: str = "") -> dict:
        """`{prefix}_min{suffix}`/`_max{suffix}`/`_mean{suffix}`/`_std{suffix}` -- matches the
        naming every quantity below already used, so this is a drop-in replacement via `**`."""
        return {
            f"{prefix}_min{suffix}": self.min,
            f"{prefix}_max{suffix}": self.max,
            f"{prefix}_mean{suffix}": self.mean,
            f"{prefix}_std{suffix}": self.std,
        }
def _spherical_area_weights(lat_deg: np.ndarray, width: int) -> np.ndarray:
    """Relative cell areas for a full-globe equirectangular grid.

    Latitude values are row centers, ordered north-to-south as in ``climate._build_grid``.
    Boundaries lie halfway between adjacent centers, with the outer boundaries at the
    poles. The area of a spherical strip is proportional to ``sin(north) - sin(south)``;
    equal-width longitude cells therefore share that row weight. The common radius and
    longitude-width factors cancel from every normalized statistic.
    """
    lat_deg = np.asarray(lat_deg, dtype=float)
    if lat_deg.ndim != 1 or lat_deg.size == 0:
        raise ValueError("lat_deg must be a non-empty one-dimensional array")
    if width <= 0:
        raise ValueError("width must be positive")
    boundaries_deg = np.empty(lat_deg.size + 1, dtype=float)
    boundaries_deg[0] = 90.0
    boundaries_deg[-1] = -90.0
    boundaries_deg[1:-1] = 0.5 * (lat_deg[:-1] + lat_deg[1:])
    row_weights = np.sin(np.radians(boundaries_deg[:-1])) - np.sin(np.radians(boundaries_deg[1:]))
    if np.any(row_weights <= 0.0):
        raise ValueError("lat_deg must be strictly ordered north-to-south")
    return np.broadcast_to(row_weights[:, None], (lat_deg.size, width))


def _weighted_fraction(mask: np.ndarray, weights: np.ndarray) -> float:
    if np.all(mask):
        return 1.0
    if not np.any(mask):
        return 0.0
    return float(weights[mask].sum() / weights.sum())


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


SPHERE_AREA_M2 = 4.0 * np.pi * (PLANET_RADIUS_KM * 1000.0) ** 2


def _node_overlap(world: World) -> dict[int, dict]:
    from .plates import compute_node_overlap

    return compute_node_overlap(world.plates)


def overlap_area_weights(world: World, overlap: dict[int, dict] | None = None) -> dict[int, np.ndarray]:
    """Per plate (keyed by plate_id, plates with nodes only), each node's share of its own
    ground: 1 / (number of plates on that spot), from `plates.compute_node_overlap`'s
    `cover_count`. Multiplying `Plate.accounting_areas_m2` by this counts overlapped ground
    once in a whole-world sum over every plate (issue #289). Pass `overlap` to reuse an
    existing `compute_node_overlap` result."""
    overlap = _node_overlap(world) if overlap is None else overlap
    return {plate_id: 1.0 / (1.0 + info["cover_count"]) for plate_id, info in overlap.items()}


def _continental_overlap_weights(overlap: dict[int, dict]) -> dict[int, np.ndarray]:
    """`overlap_area_weights` for a sum over continental plates only: the denominator counts
    only continental plates on that spot, so a continental node under an oceanic plate keeps
    its full area. Plate-level `crust_type`, matching which plates the continental sums here
    include."""
    return {plate_id: 1.0 / (1.0 + info["continental_cover_count"]) for plate_id, info in overlap.items()}


def _overlap_once_totals(world: World) -> dict[str, float]:
    """The overlap-aware sums described in this module's own docstring."""
    spacing_rad = line_spacing_rad(world.node_density)
    overlap = _node_overlap(world)
    weights = overlap_area_weights(world, overlap)
    continental_weights = _continental_overlap_weights(overlap)
    land_area = continental_volume = represented = overlapped = 0.0
    for plate in world.plates:
        if plate.plate_id not in weights:
            continue
        _, elevation = plate.all_points_and_elevation()
        area_m2 = plate.accounting_areas_m2(spacing_rad)
        once_m2 = area_m2 * weights[plate.plate_id]
        represented += float(area_m2.sum())
        overlapped += float(area_m2[weights[plate.plate_id] < 1.0].sum())
        land_area += float(once_m2[elevation > world.sea_level_m].sum())
        if plate.crust_type == "continental":
            continental_once_m2 = area_m2 * continental_weights[plate.plate_id]
            continental_volume += float(np.dot(plate.collect("crustal_thickness_m"), continental_once_m2))
    return {
        "total_land_area_dedup_km2": land_area / 1.0e6,
        "total_continental_crust_volume_dedup_km3": continental_volume / 1.0e9,
        "plate_overlap_area_fraction": overlapped / SPHERE_AREA_M2,
        "represented_area_fraction": represented / SPHERE_AREA_M2,
    }


def _total_land_area_and_continental_volume(world: World) -> tuple[float, float, int]:
    """(total land area m^2, total continental crustal volume m^3, total land node count)
    straight off `world.plates` -- see this module's own docstring for why these are computed
    here rather than off the climate grid `compute_stats` otherwise uses throughout.
    Each node is weighted by its accounting area (`Plate.accounting_areas_m2`: exact cell
    areas, since cells are not equal-area) -- one pass per plate, no
    grid resample. The node count is also `land_fraction_node`'s numerator (see `compute_stats`) -- a raw
    `elevation > sea_level_m` count, immune to the same hydrology-cache staleness as
    `land_area`/`continental_volume`, for the same reason."""
    spacing_rad = line_spacing_rad(world.node_density)
    land_nodes = 0
    land_area = 0.0
    continental_volume = 0.0
    for plate in world.plates:
        _, elevation = plate.all_points_and_elevation()
        if len(elevation) == 0:
            continue
        area_m2 = plate.accounting_areas_m2(spacing_rad)
        is_land = elevation > world.sea_level_m
        land_nodes += int(np.count_nonzero(is_land))
        land_area += float(np.sum(area_m2[is_land]))
        if plate.crust_type == "continental":
            continental_volume += float(np.dot(plate.collect("crustal_thickness_m"), area_m2))
    return land_area, continental_volume, land_nodes


def _is_water(fields: "climate.ClimateFields", is_ocean: np.ndarray) -> np.ndarray:
    """`is_ocean` plus every cell covered by a real, visible lake or sea (see this module's
    own docstring) -- the mask `land_fraction`/`ocean_fraction` actually want, as opposed to
    the narrower `is_ocean` every other stat below still uses unchanged."""
    is_lake_or_sea = fields.lake_depth_m > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
    return is_ocean | is_lake_or_sea


def compute_stats(world: World) -> dict:
    from . import continental_ledger

    hc = np.concatenate([p.collect("crustal_thickness_m") for p in world.plates]) if world.plates else np.empty(0)
    hm = np.concatenate([p.collect("mantle_lithosphere_thickness_m") for p in world.plates]) if world.plates else np.empty(0)
    hc_stats = Stat4.of(hc)
    hm_stats = Stat4.of(hm)
    fields = climate.compute_climate_cached(world)
    area_weights = _spherical_area_weights(fields.lat_deg, fields.elevation_m.shape[1])
    is_ocean, is_land = _reconcile_land_ocean(fields, world.sea_level_m)
    is_water = _is_water(fields, is_ocean)

    land_weights = area_weights[is_land]
    ocean_weights = area_weights[is_ocean]
    elevation_stats = Stat4.of(fields.elevation_m[is_land], land_weights)
    # Use the same reconciled land mask as the elevation stats. Measure the 5% band
    # relative to sea level, including its lower boundary; abs also handles worlds
    # whose highest land is a below-sea-level endorheic basin.
    land_height = fields.elevation_m[is_land] - world.sea_level_m
    land_near_max_elevation_fraction = None
    if land_height.size:
        peak = float(land_height.max())
        near_peak = land_height >= peak - 0.05 * abs(peak)
        land_near_max_elevation_fraction = float(land_weights[near_peak].sum() / land_weights.sum())

    ocean_depth = world.sea_level_m - fields.elevation_m[is_ocean]
    ocean_depth_stats = Stat4.of(ocean_depth, ocean_weights)
    land_temp_stats = Stat4.of(fields.land_temperature_c[is_land], land_weights)
    air_temp_stats = Stat4.of(fields.air_temperature_c[is_land], land_weights)
    ocean_temp_stats = Stat4.of(fields.ocean_temperature_c[is_ocean], ocean_weights)
    precip_stats = Stat4.of(fields.precipitation_mm.ravel(), area_weights.ravel())

    land_area_m2, continental_crust_volume_m3, land_node_count = _total_land_area_and_continental_volume(world)
    elevation_point_count = sum(p.node_count() for p in world.plates)

    craton = cratons.diagnostics(world)

    land_biome_ids = fields.biome_ids[is_land]
    land_weight = float(land_weights.sum())
    biome_land_fraction = {
        name: float(land_weights[land_biome_ids == i].sum()) / land_weight
        for i, name in enumerate(biomes.BIOME_NAMES)
        if i not in biomes.OCEAN_IDS and land_weight > 0.0
    }
    ocean_biome_ids = fields.biome_ids[is_ocean]
    ocean_weight = float(ocean_weights.sum())
    biome_ocean_fraction = {
        name: float(ocean_weights[ocean_biome_ids == i].sum()) / ocean_weight
        for i, name in enumerate(biomes.BIOME_NAMES)
        if i in biomes.OCEAN_IDS and ocean_weight > 0.0
    }

    return {
        "continental_material_ledger": continental_ledger.inventories(world),
        "orogenic_relief_budget": dict(orogeny.ensure_budget(world)),
        "suture_crust_transfer": crust_transfer.summary(world),
        "crust_state_area_m2": orogeny.crust_state_areas_m2(world),
        "hc_at_max_fraction": float(np.mean(hc >= lithosphere.MAX_CRUSTAL_THICKNESS_M - 1e-6)) if hc.size else None,
        **hc_stats.to_dict("hc", "_m"),
        "hm_at_max_fraction": float(np.mean(hm >= lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M - 1e-6)) if hm.size else None,
        **hm_stats.to_dict("hm", "_m"),
        "elapsed_years": world.elapsed_years,
        "plate_count": len(world.plates),
        "elevation_point_count": elevation_point_count,
        "sea_level_m": world.sea_level_m,
        "total_land_area_km2": land_area_m2 / 1.0e6,
        "total_continental_crust_volume_km3": continental_crust_volume_m3 / 1.0e9,
        **_overlap_once_totals(world),
        "land_fraction": _weighted_fraction(~is_water, area_weights),
        "ocean_fraction": _weighted_fraction(is_water, area_weights),
        # A raw `elevation > sea_level_m` node count (see `_total_land_area_and_continental_
        # volume`'s docstring) -- unlike `land_fraction` above, not connectivity-aware (an
        # enclosed sub-sea-level pit with no lake fill yet still counts as land here) and not
        # cos(lat)-weighted, but always reflects the *current* world, never `hydrology_cache`.
        # See GitHub issue #121.
        "land_fraction_node": float(land_node_count) / elevation_point_count if elevation_point_count > 0 else None,
        # True when `land_fraction`/`ocean_fraction` above were resampled from a
        # `world.hydrology_cache` that wasn't (re)built this step -- i.e. `simulate_climate_
        # biomes` has been off for at least one step since the cache was last refreshed, so
        # tectonics/eustasy may have moved the coastline out from under it. Never true while
        # `hydrology_cache` is `None`: that case computes `is_ocean` fresh (elevation-only
        # fallback, see `hydrology.sample_is_ocean`), not from a frozen cache. See issue #121.
        "land_fraction_stale": (
            world.hydrology_cache is not None
            and world.hydrology_cache_step != world.steps_taken
        ),
        "land_near_max_elevation_fraction": land_near_max_elevation_fraction,
        **elevation_stats.to_dict("elevation", "_m"),
        **ocean_depth_stats.to_dict("ocean_depth", "_m"),
        **land_temp_stats.to_dict("land_temperature", "_c"),
        **air_temp_stats.to_dict("air_temperature", "_c"),
        **ocean_temp_stats.to_dict("ocean_temperature", "_c"),
        **precip_stats.to_dict("precipitation", "_mm"),
        "biome_land_fraction": biome_land_fraction,
        "biome_ocean_fraction": biome_ocean_fraction,
        # Cratons (cratons.py): live extent, and the persisted ledger's cumulative formation
        # and destruction by mechanism -- the per-snapshot history makes destruction events
        # attributable over a long run.
        "craton_area_km2": craton["craton_area_m2"] / 1.0e6,
        "craton_continental_fraction": craton["craton_area_fraction_of_continental"],
        "craton_volume_km3": craton["craton_volume_m3"] / 1.0e9,
        "craton_destroyed_km3": sum(craton[f"craton_{key}"] for key in cratons.CRATON_SINK_ACCOUNTS) / 1.0e9,
        "craton_ledger_km3": {key: craton[f"craton_{key}"] / 1.0e9 for key in cratons.CRATON_ACCOUNTS},
    }
