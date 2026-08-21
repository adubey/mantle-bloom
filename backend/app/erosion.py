"""Weather-driven erosion: rain/sheet erosion, river-channelized erosion, and weathering,
all directly reducing elevation every step; slow/big rivers deposit part of what they carry
downstream instead of losing it all to the coast. The other direction of the weather<->
geology coupling -- terrain influencing weather (lapse-rate cooling, mountain wind
deflection, orographic rain shadow) -- already exists in climate.py, as part of its climate
pipeline. This module is the new half: climate.py's fixed grid, and hydrology.py's flow
routing over the geology node cloud, feeding back onto elevation.

**The mapping problem, and why it's easier in this direction.** climate.py already solves
node-cloud -> grid (`_sample_elevation_and_crust`'s cKDTree nearest-neighbor resample) to
build the grid in the first place. This module needs the reverse, grid -> node-cloud, which
turns out to be simpler: the climate grid is a plain regular lat/lon lattice, so a node's
own world position converts straight to (lat, lon) (geometry.xyz_to_latlon) and then to a
grid (row, col) by direct arithmetic -- no tree, no resampling, just array indexing (see
`climate_grid_indices`, whose row/col convention mirrors climate._build_grid exactly -- public,
not private, since geology.py also needs it, see that module).

**Slope is the one genuinely new piece of math climate.py's grid can't hand over for free**
(it gets slope from neighbor-index differences; an irregular node cloud has no such
structure) -- this reuses reassign.py's whole-world cKDTree pattern (build once, query k
nearest neighbors) instead: for each node, the elevation drop to the *lowest* of its nearest
neighbors (the "slope to lowest neighbor" definition used throughout this module), divided
by the real great-circle distance to that neighbor. This is a genuine dimensionless rise/run
-- not elevation-drop-per-grid-step, which would conflate slope with grid resolution -- so
RAIN_EROSION_COEFFICIENT below is tuned for this rise/run scale (order 0.001-0.1), not for
meters-of-drop-per-grid-cell (order 10-100s of meters); a coefficient tuned against that
coarser, resolution-dependent scale would make rain erosion negligible here.
WEATHERING_COEFFICIENT needs no such rescaling, since it depends only on wind speed, whose
scale (MERIDIONAL_BASE_SPEED = 6.0) is fixed by this module's own climate pipeline rather
than by any grid convention. River erosion needs the same reasoning as rain erosion, one
level further removed: it depends on `water_accum` (see hydrology.py), which is itself
downstream-accumulated *precipitation* (not a raw grid-cell count, which would implicitly
scale with resolution) -- RIVER_EROSION_COEFFICIENT is derived the same order-of-magnitude
way RAIN_EROSION_COEFFICIENT is.

**Erosion sources implemented here, and what's still out of scope.** Ocean/coastal erosion is
out of scope (a distinct source never implemented here). Weathering's vegetation boost is out
of scope (no vegetation field). River-channelized erosion, downstream deposition, and glacier
erosion -- previously out of scope for the same reason ("needs flow routing over a rotating,
irregular per-plate lattice, a separate, harder problem") -- are now implemented, see
hydrology.py for how that flow-routing graph is built and how glacier_depth itself is
grown/melted/flowed. Glacier-driven **flattening** (broad terrain smoothing under an ice
sheet, distinct from the directional erosion term) is a mantle-bloom-original addition -- see
`_flatten` below.

**No lag.** This module always erodes against the *current* step's climate and flow routing,
never a stale previous-step snapshot: climate.py here is fully stateless and cheap enough to
call fresh every step. Flow routing, though, is comparatively expensive (no JIT), so rather
than computing it twice, this module computes it once (via hydrology.compute_hydrology) and
reuses the result for both erosion and the world's cached river/lake fields
(World.hydrology_cache).
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import climate, geometry, hydrology
from .boundary import MAX_ELEVATION_M, MIN_ELEVATION_M
from .plates import PLANET_RADIUS_KM, Plate

if TYPE_CHECKING:
    from .world import World

SLOPE_NEIGHBOR_COUNT = 4

# Starting points, not final -- see module docstring. Tuned by rough order-of-magnitude
# reasoning against boundary.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR (800 m/Myr): at a
# moderately steep slope (~0.05) and moderate precipitation (~1000 mm/yr), this coefficient
# gives rain erosion the same order of magnitude as mountain-building uplift -- the design
# goal being that erosion should roughly balance typical uplift rates, not swamp them or be
# swamped by them (the raw number can't be picked from a grid-cell-based slope convention;
# see above for why).
RAIN_EROSION_COEFFICIENT = 6000.0
# Needs no rescaling the way the rain coefficient does, since wind speed's scale
# (MERIDIONAL_BASE_SPEED = 6.0) is fixed by this module's own climate pipeline, not tied to a
# grid-resolution-dependent slope convention.
WEATHERING_COEFFICIENT = 3.0
# Humidity level at which the weathering-humidity factor saturates to 1.0.
HUMIDITY_REFERENCE = 1.0
# Unlike rain/river erosion (both already multiply by `slope`), weathering had no relief
# dependence at all until this was added -- confirmed as the main reason flat, low-lying
# coastal land was drowning within a step or two of being created: wind/humidity-driven
# weathering ate a real coastal plain (whose slope is close to 0) exactly as fast as a
# mountain flank, even though its whole elevation buffer is often just a few meters. Real
# denudation is strongly relief-correlated -- flat lowlands erode far slower than steep
# terrain -- so `relief_factor` reintroduces that, same normalize-and-saturate idiom as
# `humidity_norm` above. WEATHERING_RELIEF_REFERENCE_SLOPE is the slope at which weathering
# reaches full strength; picked against a real generated world's own slope distribution
# (median land slope ~2.7e-4, far below this -- so typical flat terrain gets a small
# fraction of full weathering -- while land already at the ~90th percentile of steepness,
# ~3e-3, reaches full strength).
WEATHERING_RELIEF_REFERENCE_SLOPE = 0.002

# Stream-power river erosion: coefficient*channel_boost*water_accum^FLOW_EXPONENT*
# slope^SLOPE_EXPONENT. Coefficient re-derived by the same order-of-magnitude reasoning as
# RAIN_EROSION_COEFFICIENT (see module docstring); the two exponents need no such
# re-derivation, being dimensionless and not tied to any grid/unit convention.
RIVER_EROSION_COEFFICIENT = 100.0
RIVER_FLOW_EXPONENT = 0.5
RIVER_SLOPE_EXPONENT = 1.0
# A river preferentially re-carves its own established channel (real rivers meander within,
# not across, their valley). Channel depth is a real length (meters), so these values are
# plain physical quantities, not tied to any grid/unit convention that would need
# re-derivation.
CHANNEL_EROSION_BOOST = 0.6
CHANNEL_BOOST_REFERENCE_M = 200.0
MAX_CHANNEL_DEPTH_M = 2000.0

# Channel width is a mantle-bloom-original addition. "Larger flows make a channel larger":
# standard hydraulic-geometry scaling (width ~ discharge^b, b commonly observed near 0.5 for
# real rivers) -- same discharge exponent as RIVER_FLOW_EXPONENT above, but width is purely a
# function of how much water
# passes through, not of slope (a wide, lazy lowland river and a narrow, steep mountain
# torrent can carry comparable discharge, but width is driven by the water, not the
# gradient). Grows the same way channel_depth does -- persistent, monotonically
# non-decreasing, capped -- rather than a stateless function of this step's flow alone, so a
# river that temporarily dries up doesn't instantly narrow either.
WIDTH_GROWTH_COEFFICIENT = 50.0
WIDTH_FLOW_EXPONENT = 0.5
MAX_CHANNEL_WIDTH_M = 5000.0

# A big, slow river drops part of its sediment load locally (floodplain/delta) instead of
# carrying all of it to the coast. river_speed here is a stylized, unitless quantity (see
# hydrology.compute_river_speed) rather than a literal speed -- so its threshold isn't a
# physical speed either, just derived against this codebase's own river_speed scale.
# DEPOSITION_MIN_FLOW_M and DEPOSITION_FRACTION need no such derivation (both are already
# dimensionless/fractional, not tied to any particular scale convention).
DEPOSITION_SPEED_THRESHOLD = 2.0
DEPOSITION_MIN_FLOW_M = 0.05
DEPOSITION_FRACTION = 0.15

# Glacier erosion: scales with slope and actual accumulated ice depth (hydrology.py's
# glacier_depth, a real persistent field, not a stateless cold proxy) -- depth*slope
# approximates basal shear stress, a standard real glacial-erosion proxy: a flat-bottomed
# accumulation bowl still correctly erodes near zero regardless of ice depth.
# Temperature/precipitation-driven ice depth is a real length in consistent units (unlike the
# slope-based rain/river coefficients above, which needed re-derivation for a rise/run
# convention), so this coefficient doesn't need the same treatment.
GLACIER_EROSION_COEFFICIENT = 0.05
GLACIER_EROSION_REFERENCE_DEPTH_M = 100.0
GLACIER_EROSION_MAX_FACTOR = 2.0

# Glacier flattening is a mantle-bloom-original addition modeling how real continental ice
# sheets grind down local relief over broad areas (e.g. the Canadian Shield/Fennoscandia read
# as glacially smoothed bedrock, not just eroded-lower). Implemented as a relaxation of each
# node's elevation toward the mean of its hydrology.py flow-graph neighbors (a genuine local
# blur, not a directional erosion/deposition), scaled by the same ice_factor glacier erosion
# uses -- reuses hydrology's own k=FLOW_NEIGHBOR_COUNT neighbor graph rather than a separate
# query. GLACIER_FLATTEN_RATE_PER_MYR has no established real-world value to draw on; picked
# as a starting point, checked against a live run the same way every other from-scratch rate
# in this codebase was.
GLACIER_FLATTEN_RATE_PER_MYR = 0.2


@dataclass
class ErosionResult:
    """This step's per-node erosion breakdown, returned by apply_erosion for geology.py's own
    soil/coal/oil-gas formation (see that module) to reuse directly rather than re-deriving --
    the same "don't recompute, reuse what a prior step already produced" precedent World.
    climate_cache/World.hydrology_cache already set, except this one is threaded as a plain
    function return/argument rather than cached on World, since (unlike climate/hydrology) it
    has exactly one consumer, called from the same step_world function that produced it -- no
    later same-turn caller (rendering, /world/stats) ever needs it. Index-aligned with
    World.hydrology_cache's own arrays (points/line_refs/is_ocean in particular -- this
    dataclass deliberately doesn't repeat those, geology.py reads them off
    World.hydrology_cache directly, including its own sea_level_m-aware is_ocean rather than
    this module's own internal elevation<=0.0 shorthand): both come from an identical
    per-plate/per-line gather over the same (unreordered) world.plates within the same
    step_world call."""

    points: np.ndarray
    elevation: np.ndarray  # this step's *pre*-erosion elevation (same array hydro.elevation holds)
    slope: np.ndarray
    rain: np.ndarray
    river: np.ndarray
    weathering: np.ndarray
    sediment_deposited: np.ndarray
    temperature_c: np.ndarray
    precipitation_mm: np.ndarray


def _gather_nodes(
    world: "World",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[tuple[Plate, int, int, int]]]:
    """Every node's world position, elevation, and prior channel_depth/channel_width/
    glacier_depth, concatenated, alongside (plate, line_index, start, end) references --
    unlike reassign.py's own _gather_nodes (which needs per-node plate/line identity, since
    nodes there can move between lines), this only needs enough to slice the erosion result
    straight back onto each line's own contiguous elevation range, since erosion never moves
    a node or changes line topology."""
    points_list, elev_list, channel_list, width_list, glacier_list = [], [], [], [], []
    line_refs: list[tuple[Plate, int, int, int]] = []
    offset = 0
    for plate in world.plates:
        for line_index, line in enumerate(plate.lines):
            n = len(line.theta)
            if n == 0:
                continue
            points_list.append(line.world_xyz(plate.frame))
            elev_list.append(line.elevation)
            channel_list.append(line.channel_depth)
            width_list.append(line.channel_width)
            glacier_list.append(line.glacier_depth)
            line_refs.append((plate, line_index, offset, offset + n))
            offset += n
    if not points_list:
        return np.zeros((0, 3)), np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0), []
    return (
        np.concatenate(points_list, axis=0),
        np.concatenate(elev_list, axis=0),
        np.concatenate(channel_list, axis=0),
        np.concatenate(width_list, axis=0),
        np.concatenate(glacier_list, axis=0),
        line_refs,
    )


def compute_slope(points: np.ndarray, elevation: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-node dimensionless rise/run -- elevation drop to the *lowest* of each node's
    SLOPE_NEIGHBOR_COUNT nearest neighbors (0 if this node is already a local minimum -- the
    "slope to lowest neighbor" definition used throughout this module), divided by the real
    great-circle distance to that specific neighbor -- alongside the raw drop in meters,
    unnormalized. Returns (slope, drop_m): `slope` drives the erosion rate formula, `drop_m`
    caps it (see apply_erosion) so a single step can't erode a node *past* its lowest
    neighbor's own elevation, which would carve a new pit lower than the valley it drains
    into. The two are tracked as separate return values because `slope` here is normalized
    (dimensionless rise/run) while the cap needs the raw, unnormalized drop -- capping against
    the normalized value would bound elevation change in the wrong units entirely."""
    n = len(points)
    if n <= SLOPE_NEIGHBOR_COUNT:
        return np.zeros(n), np.zeros(n)

    tree = cKDTree(points)
    _, neighbor_idx = tree.query(points, k=SLOPE_NEIGHBOR_COUNT + 1)
    neighbor_idx = neighbor_idx[:, 1:]  # column 0 is always the point itself, at distance 0

    neighbor_elevation = elevation[neighbor_idx]
    rows = np.arange(n)
    lowest_col = np.argmin(neighbor_elevation, axis=1)
    lowest_elevation = neighbor_elevation[rows, lowest_col]
    lowest_idx = neighbor_idx[rows, lowest_col]

    drop_m = np.clip(elevation - lowest_elevation, 0.0, None)
    run_m = geometry.angular_distance(points, points[lowest_idx]) * PLANET_RADIUS_KM * 1000.0
    run_m = np.maximum(run_m, 1.0)  # avoid a divide-by-zero for (near-)coincident points
    return drop_m / run_m, drop_m


def climate_grid_indices(world_xyz: np.ndarray, height: int, width: int) -> tuple[np.ndarray, np.ndarray]:
    """Nearest climate-grid (row, col) for each node's world position -- direct array
    indexing, not a tree lookup, since (unlike the geology side) the climate grid is
    already a plain regular lat/lon lattice. Mirrors climate._build_grid's own convention
    exactly: row 0 = north pole, row increases southward; column increases eastward,
    wrapping at the antimeridian."""
    lat, lon = geometry.xyz_to_latlon(world_xyz)
    lat_deg = np.degrees(lat)
    lon_deg = np.degrees(lon)
    row = np.clip(np.floor((90.0 - lat_deg) / (180.0 / height)).astype(int), 0, height - 1)
    col = np.floor((lon_deg + 180.0) / (360.0 / width)).astype(int) % width
    return row, col


def _flatten(hydro: "hydrology.HydrologyFields", ice_factor: np.ndarray, years: float) -> np.ndarray:
    """Glacier flattening (mantle-bloom-original, see module docstring): relaxes each node's
    elevation toward the mean of its hydrology.py flow-graph neighbors, scaled by
    GLACIER_FLATTEN_RATE_PER_MYR and the same ice_factor glacier erosion uses -- a genuine
    local blur (can raise a valley or lower a peak), not a directional erosion/deposition
    term, so it's returned as a signed delta rather than folded into erosion_amount."""
    years_myr = years / 1_000_000.0
    local_mean = hydro.elevation[hydro.neighbor_idx].mean(axis=1)
    relax = 1.0 - np.exp(-GLACIER_FLATTEN_RATE_PER_MYR * ice_factor * years_myr)
    return (local_mean - hydro.elevation) * relax


def apply_erosion(world: "World", years: float) -> ErosionResult | None:
    """Erodes every plate's elevation nodes based on the world's current climate and flow
    routing -- rain/sheet erosion (precipitation x slope), river-channelized erosion
    (accumulated flow x slope, boosted by the node's own established channel), weathering
    (wind speed x humidity), and glacier erosion (accumulated ice depth x slope) -- then
    routes the combined eroded material downstream, redepositing part of it wherever a big,
    slow river drops its load (a floodplain/delta) instead of losing everything to the coast.
    Separately relaxes elevation under thick ice toward its local neighborhood mean (glacial
    flattening, see `_flatten`). Also grows channel_depth (from this step's river-erosion
    term) and channel_width (from discharge alone -- larger flows carve a wider channel);
    lake_depth/glacier_depth/silt_depth are hydrology.py's own state transitions, read directly from
    World.hydrology_cache. All persistent, see plates.ElevationLine. Mutates world.plates'
    line elevations in place; never touches node positions or line topology, so this can't
    interact with line regularization or point reassignment at all (both of those are
    purely about node density/position/ownership).

    Always computes climate fresh (never reuses World.climate_cache itself -- this runs
    right after this step's own tectonic/topology changes, so a cache from a previous step
    would already be stale for erosion's own purposes) and stores the result back onto
    World.climate_cache/World.hydrology_cache, so /world/stats and a map render don't each
    also trigger their own recomputation this same turn -- see
    climate.compute_climate_cached.

    Returns this step's ErosionResult (or None for an empty world, mirroring the
    World.hydrology_cache = None branch below) so geology.py's own soil/coal/oil-gas formation
    (called right after this, from world.step_world) can reuse these same per-node terms
    instead of re-deriving them -- see ErosionResult's own docstring."""
    fields = climate.compute_climate(world, *climate.grid_dimensions(world.climate_density))
    world.climate_cache = fields

    points, elevation, prior_channel_depth, prior_channel_width, prior_glacier_depth, line_refs = _gather_nodes(world)
    n = len(points)
    if n == 0:
        world.hydrology_cache = None
        return None

    height, width = fields.precipitation_mm.shape
    row, col = climate_grid_indices(points, height, width)
    precipitation_mm = fields.precipitation_mm[row, col]
    wind_speed = np.hypot(fields.wind_u, fields.wind_v)[row, col]
    humidity = fields.humidity[row, col]
    is_ocean_node = elevation <= 0.0
    # The same real temperature a node actually experiences that render_image.py's own
    # temperature view displays -- ocean surface over water, moderated air over land.
    temperature = np.where(is_ocean_node, fields.ocean_temperature_c[row, col], fields.air_temperature_c[row, col])

    slope, drop_to_lowest_neighbor_m = compute_slope(points, elevation)
    dt_myr = years / 1_000_000.0

    hydro = hydrology.compute_hydrology(world, precipitation_mm, temperature, years)
    world.hydrology_cache = hydro
    for message in hydro.lake_events:
        world.log_event(message)
    water_accum_m = hydro.flow_accum / 1000.0

    rain = RAIN_EROSION_COEFFICIENT * slope * (precipitation_mm / 1000.0) * dt_myr
    channel_boost = 1.0 + CHANNEL_EROSION_BOOST * np.clip(prior_channel_depth / CHANNEL_BOOST_REFERENCE_M, 0.0, 1.0)
    river = (
        RIVER_EROSION_COEFFICIENT
        * channel_boost
        * np.power(np.clip(water_accum_m, 0.0, None), RIVER_FLOW_EXPONENT)
        * np.power(slope, RIVER_SLOPE_EXPONENT)
        * dt_myr
    )
    humidity_norm = np.clip(humidity / HUMIDITY_REFERENCE, 0.0, 1.0)
    relief_factor = np.clip(slope / WEATHERING_RELIEF_REFERENCE_SLOPE, 0.0, 1.0)
    weathering = WEATHERING_COEFFICIENT * wind_speed * humidity_norm * relief_factor * dt_myr
    ice_factor = np.clip(prior_glacier_depth / GLACIER_EROSION_REFERENCE_DEPTH_M, 0.0, GLACIER_EROSION_MAX_FACTOR)
    glacier = GLACIER_EROSION_COEFFICIENT * slope * ice_factor * dt_myr
    # Capped at the drop to the lowest neighbor so a single step can't erode a node below the
    # valley floor it drains into. Zeroed over ocean nodes (elevation <= sea level, the same
    # convention climate.py/plates.py use everywhere else): every source here is a subaerial
    # process -- coastal/ocean erosion is a separate source that would touch the seafloor, and
    # that's the one source this module still doesn't implement (see module docstring).
    erosion_amount = np.where(is_ocean_node, 0.0, np.clip(rain + river + weathering + glacier, 0.0, None))
    erosion_amount = np.minimum(erosion_amount, drop_to_lowest_neighbor_m)

    # Deposition: wherever a big (water_accum_m > DEPOSITION_MIN_FLOW_M), slow
    # (river_speed < DEPOSITION_SPEED_THRESHOLD) river passes through, DEPOSITION_FRACTION
    # of the material passing through settles right there instead of continuing downstream
    # -- route_downstream still conserves the total exactly either way.
    river_speed = hydrology.compute_river_speed(slope, hydro.flow_accum)
    is_depositing = (river_speed < DEPOSITION_SPEED_THRESHOLD) & (water_accum_m > DEPOSITION_MIN_FLOW_M)
    retain_fraction = np.where(is_depositing, DEPOSITION_FRACTION, 0.0)
    _, sediment_deposited = hydrology.route_downstream(elevation, is_ocean_node, hydro.flow_target, erosion_amount, retain_fraction=retain_fraction)

    flatten_delta = _flatten(hydro, ice_factor, years)

    new_elevation = np.clip(elevation - erosion_amount + sediment_deposited + flatten_delta, MIN_ELEVATION_M, MAX_ELEVATION_M)
    new_channel_depth = np.where(is_ocean_node, 0.0, np.clip(prior_channel_depth + river, 0.0, MAX_CHANNEL_DEPTH_M))
    # Width grows with discharge alone (no slope/channel_boost term -- see module constants'
    # own comment for why), same persistent/monotonic/capped shape as depth.
    width_growth = WIDTH_GROWTH_COEFFICIENT * np.power(np.clip(water_accum_m, 0.0, None), WIDTH_FLOW_EXPONENT) * dt_myr
    new_channel_width = np.where(is_ocean_node, 0.0, np.clip(prior_channel_width + width_growth, 0.0, MAX_CHANNEL_WIDTH_M))

    # theta (and therefore every other parallel array's shape) is never touched here -- see
    # plates.ElevationLine's own docstring for why dataclasses.replace, not an explicit
    # field-by-field reconstruction, is the pattern that stays correct as more persistent
    # fields get added (this site used to silently wipe is_volcano to False every step).
    for plate, line_index, start, end in line_refs:
        line = plate.lines[line_index]
        plate.lines[line_index] = dataclasses.replace(
            line,
            elevation=new_elevation[start:end],
            channel_depth=new_channel_depth[start:end],
            channel_width=new_channel_width[start:end],
            lake_depth=hydro.lake_depth[start:end],
            glacier_depth=hydro.glacier_depth[start:end],
            silt_depth=hydro.silt_depth[start:end],
        )

    return ErosionResult(
        points=points,
        elevation=elevation,
        slope=slope,
        rain=rain,
        river=river,
        weathering=weathering,
        sediment_deposited=sediment_deposited,
        temperature_c=temperature,
        precipitation_mm=precipitation_mm,
    )
