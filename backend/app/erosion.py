"""Weather-driven erosion: rain/sheet erosion, river-channelized erosion, and weathering,
all directly reducing elevation every step; slow/big rivers deposit part of what they carry
downstream instead of losing it all to the coast. The other direction of the weather<->
geology coupling -- terrain influencing weather (lapse-rate cooling, mountain wind
deflection, orographic rain shadow) -- already exists in climate.py, ported from plate-sim
alongside the rest of its climate pipeline. This module is the new half: climate.py's fixed
grid, and hydrology.py's flow routing over the geology node cloud, feeding back onto
elevation.

**The mapping problem, and why it's easier in this direction.** climate.py already solves
node-cloud -> grid (`_sample_elevation_and_crust`'s cKDTree nearest-neighbor resample) to
build the grid in the first place. This module needs the reverse, grid -> node-cloud, which
turns out to be simpler: the climate grid is a plain regular lat/lon lattice, so a node's
own world position converts straight to (lat, lon) (geometry.xyz_to_latlon) and then to a
grid (row, col) by direct arithmetic -- no tree, no resampling, just array indexing (see
`_climate_grid_indices`, whose row/col convention mirrors climate._build_grid exactly).

**Slope is the one genuinely new piece of math climate.py's grid can't hand over for free**
(it gets slope from neighbor-index differences; an irregular node cloud has no such
structure) -- this reuses reassign.py's whole-world cKDTree pattern (build once, query k
nearest neighbors) instead: for each node, the elevation drop to the *lowest* of its nearest
neighbors (matching plate-sim's own "slope to lowest neighbor" definition), divided by the
real great-circle distance to that neighbor. This is a genuine dimensionless rise/run, unlike
plate-sim's own slope -- documented there as a known simplification, "elevation drop per grid
step, not drop-over-real-distance" -- so RAIN_EROSION_COEFFICIENT below is *not* plate-sim's
own 0.15: that value was tuned against meters-of-drop-per-grid-cell (order 10-100s of
meters), not a true rise/run (order 0.001-0.1), and porting it verbatim would make rain
erosion negligible. WEATHERING_COEFFICIENT ports more directly, since wind speed uses the
exact same scale in both codebases (MERIDIONAL_BASE_SPEED = 6.0 in each). River erosion
needs the same reasoning as rain erosion, one level further removed: it depends on
`water_accum` (see hydrology.py), which is itself downstream-accumulated *precipitation*
(not a raw grid-cell count the way plate-sim's own accumulation implicitly scales with
resolution) -- RIVER_EROSION_COEFFICIENT is re-derived the same way RAIN_EROSION_COEFFICIENT
was, not ported verbatim.

**Scope cut from plate-sim's five erosion sources.** Ocean/coastal erosion is still dropped
(a distinct source never implemented here). Weathering's vegetation boost is still dropped
(no vegetation field). River-channelized erosion, downstream deposition, and glacier erosion
-- previously dropped for the same reason ("needs flow routing over a rotating, irregular
per-plate lattice, a separate, harder problem") -- are now implemented, see hydrology.py for
how that flow-routing graph is built and how glacier_depth itself is grown/melted/flowed.
Glacier-driven **flattening** (broad terrain smoothing under an ice sheet, distinct from the
directional erosion term) is a mantle-bloom-original addition, not a plate-sim port -- see
`_flatten` below.

**No lag, unlike plate-sim.** plate-sim's erosion reads the *previous* step's climate and
recomputes its own flow routing independently of hydrology.py's own pass (an accepted
redundancy there, cheap under numba JIT). climate.py here is fully stateless and cheap
enough to call fresh; flow routing here is comparatively expensive (no JIT), so this module
computes it once (via hydrology.compute_hydrology) and reuses the result for both erosion
and the world's cached river/lake fields (World.hydrology_cache), rather than paying for it
twice.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import climate, geometry, hydrology
from .boundary import MAX_ELEVATION_M, MIN_ELEVATION_M
from .plates import PLANET_RADIUS_KM, ElevationLine, Plate

if TYPE_CHECKING:
    from .world import World

SLOPE_NEIGHBOR_COUNT = 4

# Starting points, not final -- see module docstring. Tuned by rough order-of-magnitude
# reasoning against boundary.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR (800 m/Myr): at a
# moderately steep slope (~0.05) and moderate precipitation (~1000 mm/yr), this coefficient
# gives rain erosion the same order of magnitude as mountain-building uplift, matching
# plate-sim's own stated design goal ("tuned so it roughly balances typical uplift rates")
# without being able to reuse its actual number (different slope units -- see above).
RAIN_EROSION_COEFFICIENT = 6000.0
# Ported close to plate-sim's own 3.0 -- wind speed uses the same scale in both codebases
# (MERIDIONAL_BASE_SPEED = 6.0), so this needs less re-derivation than the rain coefficient.
WEATHERING_COEFFICIENT = 3.0
# Humidity level at which the weathering-humidity factor saturates to 1.0 -- same value and
# meaning as plate-sim's own HUMIDITY_REFERENCE.
HUMIDITY_REFERENCE = 1.0

# Stream-power river erosion: coefficient*channel_boost*water_accum^FLOW_EXPONENT*
# slope^SLOPE_EXPONENT. Coefficient re-derived (not ported, see module docstring) by the
# same order-of-magnitude reasoning as RAIN_EROSION_COEFFICIENT; the two exponents port
# directly from plate-sim (dimensionless, not tied to any grid/unit convention).
RIVER_EROSION_COEFFICIENT = 100.0
RIVER_FLOW_EXPONENT = 0.5
RIVER_SLOPE_EXPONENT = 1.0
# A river preferentially re-carves its own established channel (real rivers meander within,
# not across, their valley) -- ported values/meaning directly from plate-sim; channel depth
# is a real length (meters) in both codebases, so these don't need re-derivation.
CHANNEL_EROSION_BOOST = 0.6
CHANNEL_BOOST_REFERENCE_M = 200.0
MAX_CHANNEL_DEPTH_M = 2000.0

# A big, slow river drops part of its sediment load locally (floodplain/delta) instead of
# carrying all of it to the coast. river_speed here is a stylized, unitless quantity (see
# hydrology.compute_river_speed) rather than a literal speed, same as plate-sim's own -- so
# its threshold isn't a physical speed either, just re-derived (not ported) against this
# codebase's own river_speed scale. DEPOSITION_MIN_FLOW_M and DEPOSITION_FRACTION port
# directly (both already dimensionless/fractional, not grid-unit-dependent).
DEPOSITION_SPEED_THRESHOLD = 2.0
DEPOSITION_MIN_FLOW_M = 0.05
DEPOSITION_FRACTION = 0.15

# Glacier erosion: scales with slope and actual accumulated ice depth (hydrology.py's
# glacier_depth, a real persistent field, not a stateless cold proxy) -- depth*slope
# approximates basal shear stress, a standard real glacial-erosion proxy: a flat-bottomed
# accumulation bowl still correctly erodes near zero regardless of ice depth. Ported
# directly from plate-sim -- temperature/precipitation-driven ice depth uses the same units
# in both codebases (unlike the slope-based rain/river coefficients above), so this doesn't
# need the same re-derivation.
GLACIER_EROSION_COEFFICIENT = 0.05
GLACIER_EROSION_REFERENCE_DEPTH_M = 100.0
GLACIER_EROSION_MAX_FACTOR = 2.0

# Glacier flattening -- NOT a plate-sim mechanic (confirmed: no flatten/smooth/scour effect
# tied to ice cover exists there at all), a mantle-bloom-original addition modeling how real
# continental ice sheets grind down local relief over broad areas (e.g. the Canadian
# Shield/Fennoscandia read as glacially smoothed bedrock, not just eroded-lower). Implemented
# as a relaxation of each node's elevation toward the mean of its hydrology.py flow-graph
# neighbors (a genuine local blur, not a directional erosion/deposition), scaled by the same
# ice_factor glacier erosion uses -- reuses hydrology's own k=FLOW_NEIGHBOR_COUNT neighbor
# graph rather than a separate query. GLACIER_FLATTEN_RATE_PER_MYR has no real-world number
# to port; picked as a starting point, checked against a live run the same way every other
# from-scratch rate in this codebase was.
GLACIER_FLATTEN_RATE_PER_MYR = 0.2


def _gather_nodes(
    world: "World",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[tuple[Plate, int, int, int]]]:
    """Every node's world position, elevation, and prior channel_depth/glacier_depth,
    concatenated, alongside (plate, line_index, start, end) references -- unlike
    reassign.py's own _gather_nodes (which needs per-node plate/line identity, since nodes
    there can move between lines), this only needs enough to slice the erosion result
    straight back onto each line's own contiguous elevation range, since erosion never moves
    a node or changes line topology."""
    points_list, elev_list, channel_list, glacier_list = [], [], [], []
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
            glacier_list.append(line.glacier_depth)
            line_refs.append((plate, line_index, offset, offset + n))
            offset += n
    if not points_list:
        return np.zeros((0, 3)), np.zeros(0), np.zeros(0), np.zeros(0), []
    return (
        np.concatenate(points_list, axis=0),
        np.concatenate(elev_list, axis=0),
        np.concatenate(channel_list, axis=0),
        np.concatenate(glacier_list, axis=0),
        line_refs,
    )


def _compute_slope(points: np.ndarray, elevation: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-node dimensionless rise/run -- elevation drop to the *lowest* of each node's
    SLOPE_NEIGHBOR_COUNT nearest neighbors (0 if this node is already a local minimum,
    matching plate-sim's own "slope to lowest neighbor" definition), divided by the real
    great-circle distance to that specific neighbor -- alongside the raw drop in meters,
    unnormalized. Returns (slope, drop_m): `slope` drives the erosion rate formula, `drop_m`
    caps it (see apply_erosion) so a single step can't erode a node *past* its lowest
    neighbor's own elevation, which would carve a new pit lower than the valley it drains
    into -- the same cap plate-sim's own erosion applies (there, "slope" -- its own raw
    per-grid-cell drop, not a rise/run -- doubles directly as that cap; here the two need to
    be tracked separately, since this module's `slope` is normalized and plate-sim's isn't)."""
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


def _climate_grid_indices(world_xyz: np.ndarray, height: int, width: int) -> tuple[np.ndarray, np.ndarray]:
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


def apply_erosion(world: "World", years: float) -> None:
    """Erodes every plate's elevation nodes based on the world's current climate and flow
    routing -- rain/sheet erosion (precipitation x slope), river-channelized erosion
    (accumulated flow x slope, boosted by the node's own established channel), weathering
    (wind speed x humidity), and glacier erosion (accumulated ice depth x slope) -- then
    routes the combined eroded material downstream, redepositing part of it wherever a big,
    slow river drops its load (a floodplain/delta) instead of losing everything to the coast.
    Separately relaxes elevation under thick ice toward its local neighborhood mean (glacial
    flattening, see `_flatten`). Also grows channel_depth (from this step's river-erosion
    term); lake_depth/glacier_depth are hydrology.py's own state transitions, read directly
    from World.hydrology_cache. All three persistent, see plates.ElevationLine. Mutates
    world.plates' line elevations in place; never touches node positions or line topology,
    so this can't interact with line regularization or point reassignment at all (both of
    those are purely about node density/position/ownership).

    Always computes climate fresh (never reuses World.climate_cache itself -- this runs
    right after this step's own tectonic/topology changes, so a cache from a previous step
    would already be stale for erosion's own purposes) and stores the result back onto
    World.climate_cache/World.hydrology_cache, so /world/stats and a map render don't each
    also trigger their own recomputation this same turn -- see
    climate.compute_climate_cached."""
    fields = climate.compute_climate(world)
    world.climate_cache = fields

    points, elevation, prior_channel_depth, prior_glacier_depth, line_refs = _gather_nodes(world)
    n = len(points)
    if n == 0:
        world.hydrology_cache = None
        return

    height, width = fields.precipitation_mm.shape
    row, col = _climate_grid_indices(points, height, width)
    precipitation_mm = fields.precipitation_mm[row, col]
    wind_speed = np.hypot(fields.wind_u, fields.wind_v)[row, col]
    humidity = fields.humidity[row, col]
    is_ocean_node = elevation <= 0.0
    # The same real temperature a node actually experiences that render_image.py's own
    # temperature view displays -- ocean surface over water, moderated air over land.
    temperature = np.where(is_ocean_node, fields.ocean_temperature_c[row, col], fields.air_temperature_c[row, col])

    slope, drop_to_lowest_neighbor_m = _compute_slope(points, elevation)
    dt_myr = years / 1_000_000.0

    hydro = hydrology.compute_hydrology(world, precipitation_mm, temperature, years)
    world.hydrology_cache = hydro
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
    weathering = WEATHERING_COEFFICIENT * wind_speed * humidity_norm * dt_myr
    ice_factor = np.clip(prior_glacier_depth / GLACIER_EROSION_REFERENCE_DEPTH_M, 0.0, GLACIER_EROSION_MAX_FACTOR)
    glacier = GLACIER_EROSION_COEFFICIENT * slope * ice_factor * dt_myr
    # Capped at the drop to the lowest neighbor -- same reason plate-sim caps its own
    # erosion at "slope" -- so a single step can't erode a node below the valley floor it
    # drains into. Zeroed over ocean nodes (elevation <= sea level, the same convention
    # climate.py/plates.py use everywhere else): every source here is a subaerial process,
    # and plate-sim itself excludes ocean cells from these same sources (its separate
    # coastal-current erosion source is the only one that touches the seafloor, and that's
    # the one source this module still doesn't implement -- see module docstring).
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

    for plate, line_index, start, end in line_refs:
        line = plate.lines[line_index]
        plate.lines[line_index] = ElevationLine(
            phi=line.phi,
            theta=line.theta,
            elevation=new_elevation[start:end],
            channel_depth=new_channel_depth[start:end],
            lake_depth=hydro.lake_depth[start:end],
            glacier_depth=hydro.glacier_depth[start:end],
        )
