"""Weather-driven erosion: rain/sheet erosion and weathering, both directly reducing
elevation every step. The other direction of the weather<->geology coupling -- terrain
influencing weather (lapse-rate cooling, mountain wind deflection, orographic rain shadow)
-- already exists in climate.py, ported from plate-sim alongside the rest of its climate
pipeline. This module is the new half: climate.py's fixed grid feeding back onto geology's
irregular per-plate node cloud.

**The mapping problem, and why it's easier in this direction.** climate.py already solves
node-cloud -> grid (`_sample_elevation_and_crust`'s cKDTree nearest-neighbor resample) to
build the grid in the first place. This module needs the reverse, grid -> node-cloud, which
turns out to be simpler: the climate grid is a plain regular lat/lon lattice, so a node's
own world position converts straight to (lat, lon) (geometry.xyz_to_latlon) and then to a
grid (row, col) by direct arithmetic -- no tree, no resampling, just array indexing (see
`_climate_grid_indices`, whose row/col convention mirrors climate._build_grid exactly).

**Slope is the one genuinely new piece of math.** climate.py's grid gets slope for free
from neighbor-index differences; an irregular node cloud has no such structure. This reuses
reassign.py's whole-world cKDTree pattern (build once, query k nearest neighbors) instead:
for each node, the elevation drop to the *lowest* of its nearest neighbors (matching
plate-sim's own "slope to lowest neighbor" definition), divided by the real great-circle
distance to that neighbor. This is a genuine dimensionless rise/run, unlike plate-sim's own
slope -- documented there as a known simplification, "elevation drop per grid step, not
drop-over-real-distance" -- so RAIN_EROSION_COEFFICIENT below is *not* plate-sim's own
0.15: that value was tuned against meters-of-drop-per-grid-cell (order 10-100s of meters),
not a true rise/run (order 0.001-0.1), and porting it verbatim would make rain erosion
negligible. WEATHERING_COEFFICIENT ports more directly, since wind speed uses the exact
same scale in both codebases (MERIDIONAL_BASE_SPEED = 6.0 in each).

**Scope cut from plate-sim's five erosion sources, for two reasons -- see the module
docstring precedent (climate.py's own "why not ported" section) for the same style of
call.** River/channelized erosion and deposition both need downhill flow-routing
(plate-sim's `flow_target`/D8 accumulation), which assumes persistent per-cell state on a
fixed grid -- building that over a rotating, irregular per-plate lattice is a separate, much
harder problem, not attempted here. Glacier erosion is dropped because mantle-bloom has no
glacier field at all. Vegetation's boost to weathering is dropped because mantle-bloom has
no vegetation field (deliberately -- see climate.py's own "deliberately not ported" note).
What's left, rain/sheet erosion and weathering, are the two sources that don't need any of
that missing infrastructure -- and erosion here is one-way: material is removed, never
redeposited elsewhere.

**No lag, unlike plate-sim.** plate-sim's erosion reads the *previous* step's climate,
because climate there is only computed once per step, late in a long per-step pipeline
(`_finalize`), after erosion has already run. climate.py here is already fully stateless and
cheap to call on demand, so this module just calls `climate.compute_climate(world)` fresh,
immediately before using it -- no staleness to reason about.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import climate, geometry
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


def _gather_nodes(world: "World") -> tuple[np.ndarray, np.ndarray, list[tuple[Plate, int, int, int]]]:
    """Every node's world position and elevation, concatenated, alongside
    (plate, line_index, start, end) references -- unlike reassign.py's own _gather_nodes
    (which needs per-node plate/line identity, since nodes there can move between lines),
    this only needs enough to slice the erosion result straight back onto each line's own
    contiguous elevation range, since erosion never moves a node or changes line topology."""
    points_list, elev_list = [], []
    line_refs: list[tuple[Plate, int, int, int]] = []
    offset = 0
    for plate in world.plates:
        for line_index, line in enumerate(plate.lines):
            n = len(line.theta)
            if n == 0:
                continue
            points_list.append(line.world_xyz(plate.frame))
            elev_list.append(line.elevation)
            line_refs.append((plate, line_index, offset, offset + n))
            offset += n
    if not points_list:
        return np.zeros((0, 3)), np.zeros(0), []
    return np.concatenate(points_list, axis=0), np.concatenate(elev_list, axis=0), line_refs


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


def apply_erosion(world: "World", years: float) -> None:
    """Erodes every plate's elevation nodes based on the world's current climate --
    rain/sheet erosion (precipitation x slope) plus weathering (wind speed x humidity).
    Mutates world.plates' line elevations in place; never touches node positions or line
    topology, so this can't interact with line regularization or point reassignment at all
    (both of those are purely about node density/position/ownership).

    Always computes climate fresh (never reuses World.climate_cache itself -- this runs
    right after this step's own tectonic/topology changes, so a cache from a previous step
    would already be stale for erosion's own purposes) and stores the result back onto
    World.climate_cache, so /world/stats and a climate map render don't each also trigger
    their own recomputation this same turn -- see climate.compute_climate_cached."""
    fields = climate.compute_climate(world)
    world.climate_cache = fields

    points, elevation, line_refs = _gather_nodes(world)
    n = len(points)
    if n == 0:
        return

    height, width = fields.precipitation_mm.shape
    row, col = _climate_grid_indices(points, height, width)
    precipitation_mm = fields.precipitation_mm[row, col]
    wind_speed = np.hypot(fields.wind_u, fields.wind_v)[row, col]
    humidity = fields.humidity[row, col]

    slope, drop_to_lowest_neighbor_m = _compute_slope(points, elevation)
    dt_myr = years / 1_000_000.0

    rain = RAIN_EROSION_COEFFICIENT * slope * (precipitation_mm / 1000.0) * dt_myr
    humidity_norm = np.clip(humidity / HUMIDITY_REFERENCE, 0.0, 1.0)
    weathering = WEATHERING_COEFFICIENT * wind_speed * humidity_norm * dt_myr
    # Capped at the drop to the lowest neighbor -- same reason plate-sim caps its own
    # erosion at "slope" -- so a single step can't carve a node below the valley floor it
    # drains into. Zeroed over ocean nodes (elevation <= sea level, the same convention
    # climate.py/plates.py use everywhere else): rain/weathering erosion is a subaerial
    # process, and plate-sim itself excludes ocean cells from these same two sources (its
    # separate coastal-current erosion source is the only one that touches the seafloor,
    # and that's one of the sources this module deliberately doesn't implement -- see
    # module docstring).
    is_ocean_node = elevation <= 0.0
    erosion_amount = np.where(is_ocean_node, 0.0, np.clip(rain + weathering, 0.0, None))
    erosion_amount = np.minimum(erosion_amount, drop_to_lowest_neighbor_m)

    for plate, line_index, start, end in line_refs:
        line = plate.lines[line_index]
        new_elevation = np.clip(line.elevation - erosion_amount[start:end], MIN_ELEVATION_M, MAX_ELEVATION_M)
        plate.lines[line_index] = ElevationLine(phi=line.phi, theta=line.theta, elevation=new_elevation)
