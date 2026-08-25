"""Plates as spherical polygons, each carrying its own set of `ElevationLine`s.

Each plate owns a rotation matrix (`frame`) mapping its local (phi, theta) spherical
coordinates to world unit vectors (see `geometry.plate_frame_from_seed`). Rotating a plate
rigidly only ever updates `frame` -- the (phi, theta) node coordinates themselves never
change, so rotation never needs resampling. See docs/simulation-model.md for the full design
writeup, and elevation_lines.py for the node representation itself (`ElevationLine`, node
density/spacing, and periodic line regularization) -- this module is about the plates that
carry it: identity, motion, territory, and generation.
"""

from __future__ import annotations

import abc
from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterator, Protocol

import numpy as np
from scipy.spatial import ConvexHull, QhullError, cKDTree

from . import ellipse, geometry, mantle
from .elevation_lines import (
    DEFAULT_NODE_DENSITY,
    ERUPTION_ELEVATION_M,
    MAX_ELEVATION_M,
    MIN_ELEVATION_M,
    NODE_DENSITY_CHOICES,
    PLANET_RADIUS_KM,
    TARGET_LINE_SPACING_KM,
    TARGET_LINE_SPACING_RAD,
    VOLCANO_ACTIVE_MAX_YEARS,
    VOLCANO_ACTIVE_MIN_YEARS,
    ElevationLine,
    ElevationPoint,
    build_lines_from_lattice,
    install_point_field_accessors,
    iter_local_lattice,
    line_spacing_rad,
    needs_regularizing,
    regularize_line,
)
from .lat_long_grid import LatLongGrid
from .noise import SphereNoise
from .rtree_index import RTree

if TYPE_CHECKING:
    from .world import World

CONTINENTAL_FRACTION = 0.4
BASE_CONTINENTAL_M = 200.0
BASE_OCEANIC_M = -3800.0
CONTINENTAL_NOISE_AMPLITUDE_M = 1200.0
OCEANIC_NOISE_AMPLITUDE_M = 500.0

# Plate count is chosen automatically (see generate_plates) rather than asked of the user --
# an inclusive range of plausible Earth-like plate counts.
MIN_AUTO_PLATES = 8
MAX_AUTO_PLATES = 20
# Both user-facing (UI sliders) -- see generate_plates' continental_fraction/land_fraction.
DEFAULT_CONTINENTAL_FRACTION = 0.70
DEFAULT_LAND_FRACTION = 0.29
# However high the requested continental fraction, still leave room for real ocean floor.
MIN_OCEANIC_PLATES = 3
# Resolution for the one-off whole-sphere sweep generate_plates uses to translate a
# requested land_fraction into a concrete noise threshold (see _land_noise_threshold) --
# coarser than the simulation/render grids since this only needs to be a statistically
# representative sample, not something visually smooth or physically carried.
LAND_FRACTION_SAMPLE_SPACING_KM = 150.0
LAND_FRACTION_SAMPLE_SPACING_RAD = LAND_FRACTION_SAMPLE_SPACING_KM / PLANET_RADIUS_KM


class SpherePolygon(Protocol):
    """A region of the sphere's surface that can answer point-membership queries.
    `Plate` (below) is the only implementer today, but this is kept representation-agnostic
    -- structural, not a base class -- so anything shaped like a polygon on a sphere can
    satisfy it without inheriting from `Plate`."""

    def contains(self, lat: float, lon: float) -> bool:
        """True if the geographic point (lat, lon, radians) falls inside this polygon."""
        ...


# Two plates count as neighbours once the closest points of their two outlines come within
# this many multiples of the default line spacing -- generous enough to still catch plates
# separated by a boundary's own elevation-transition zone (see boundary.py's
# FAR_THRESHOLD_RAD, a similar multiple of spacing_rad), while not matching plates that
# merely share the same hemisphere.
NEIGHBOUR_DISTANCE_RAD = 6.0 * TARGET_LINE_SPACING_RAD


def _plates_within(plate: "Plate", all_plates: list["Plate"], threshold_rad: float) -> list["Plate"]:
    """Every other plate in `all_plates` whose outline (`Plate.get_bounding_polygon`) comes
    within `threshold_rad` of `plate`'s own -- shared by both PlateWithLines.get_neighbours
    and PlateWithRTree.get_neighbours, which differ only in what outline_world() itself does
    for that representation. Uses the cached get_bounding_polygon() rather than outline_world()
    directly since this runs once per plate per call, each time re-reading every other
    plate's own outline -- an O(n) set of calls across all_plates that would otherwise
    recompute the same unchanged outlines from scratch every time (see world.py's callers,
    e.g. volcanism.merge_close_volcanic_fields' own per-field loop). Same bounding-sphere
    prefilter as boundary.py's step_boundaries (cheap enough to compute per plate, and enough
    to rule out most pairs before a real nearest-point query)."""
    own_points = plate.get_bounding_polygon()
    if len(own_points) == 0:
        return []
    own_centroid, own_radius = geometry.bounding_sphere(own_points)

    neighbours = []
    for other in all_plates:
        if other.plate_id == plate.plate_id:
            continue
        other_points = other.get_bounding_polygon()
        if len(other_points) == 0:
            continue
        other_centroid, other_radius = geometry.bounding_sphere(other_points)
        centroid_dist = float(geometry.angular_distance(own_centroid, other_centroid))
        if centroid_dist - own_radius - other_radius > threshold_rad:
            continue
        other_tree = other.get_bounding_polygon_tree()
        if other_tree is None:
            continue
        # Only the <= threshold_rad decision below is ever used, not the actual distance --
        # distance_upper_bound lets cKDTree stop descending into a branch as soon as it can
        # prove that branch can't beat the bound, rather than finding every point's true
        # global nearest neighbour only to immediately compare it against the same bound.
        # Points with no neighbour within threshold_rad come back as +inf, which still
        # correctly fails the comparison below.
        closest_dist = float(other_tree.query(own_points, distance_upper_bound=threshold_rad)[0].min())
        if closest_dist <= threshold_rad:
            neighbours.append(other)
    return neighbours


def _contested_by_any(points_xyz: np.ndarray, neighbours: list["Plate"]) -> np.ndarray:
    """`geometry.points_in_any_spherical_polygon`, OR-ed across every neighbour's own
    `contains_batch` instead of a shared polygon-list winding test -- lets a `PlateWithLines`
    neighbour answer via its own O(log rows) fast path (see `PlateWithLines.contains_batch`)
    rather than every neighbour paying the full winding-number cost regardless of
    representation. Same semantics otherwise: stops early once every point is already
    contested by some earlier neighbour, all-`False` if either input is empty."""
    n = len(points_xyz)
    contested = np.zeros(n, dtype=bool)
    if n == 0 or not neighbours:
        return contested
    for neighbour in neighbours:
        contested |= neighbour.contains_batch(points_xyz)
        if np.all(contested):
            break
    return contested


def _update_to_lat_long_grid(plate: "Plate", grid: LatLongGrid) -> None:
    """The plate -> grid half of the `LatLongGrid` round trip (see lat_long_grid.py's own
    module docstring): resamples this plate's current node elevations onto whichever grid
    cell each node's world position falls nearest to. Shared by PlateWithLines.
    update_to_lat_long_grid and PlateWithRTree.update_to_lat_long_grid, which differ only in
    what all_points_and_elevation() itself gathers -- the same "representation differs, the
    traversal against the abstract Plate interface doesn't" precedent _plates_within already
    sets for get_neighbours, above."""
    points, elevation = plate.all_points_and_elevation()
    if len(points) == 0:
        return
    rows, cols = grid.row_col_for_world_xyz(points)
    grid.set_elevation(rows, cols, elevation)


def _update_deltas_from_lat_long_grid(plate: "Plate", grid: LatLongGrid) -> None:
    """The grid -> plate half of the round trip: only the *net* elevation change grid-space
    code made since update_to_lat_long_grid populated the grid (see
    LatLongGrid.change_elevation) is applied back, added onto each node's own current
    elevation -- never the grid's absolute value, which is only a coarse nearest-cell resample
    and would otherwise flatten every node sharing a cell down to the exact same elevation.
    Clipped the same way every other elevation-modifying pass in this codebase is (boundary.py,
    erosion.py, volcanism.py), since a node's own elevation plus a grid delta sized for a
    different, coarser node can otherwise drift past the world's physical bounds."""
    points, _ = plate.all_points_and_elevation()
    if len(points) == 0:
        return
    rows, cols = grid.row_col_for_world_xyz(points)
    deltas = grid.delta_at(rows, cols)
    for (point, _), delta in zip(plate.map_world_points(), deltas):
        if delta != 0.0:
            point.set_elevation(float(np.clip(point.get_elevation() + delta, MIN_ELEVATION_M, MAX_ELEVATION_M)))


def _lines_from_resample(
    frame: np.ndarray,
    points: np.ndarray,
    elevation: np.ndarray,
    coverage_radius_rad: float,
    spacing_rad: float,
    exclude_tree: cKDTree | None = None,
) -> list[ElevationLine]:
    """Sweep `frame`'s own local lattice and keep whichever candidate nodes fall within
    `coverage_radius_rad` of `points` -- and, if `exclude_tree` is given, strictly closer to
    `points` than to whatever `exclude_tree` indexes (every *other* live plate's own current
    nodes, see `PlateWithLines._merge_nodes_with`) -- each claimed node's elevation coming
    from its own nearest point in `points`. Shared by `PlateWithLines._merge_nodes_with`
    (which needs the exclusivity check, since a merging pair's old points can't be trusted to
    sit clear of a third plate's own territory) and `PlateWithLines.grow_into` (which doesn't:
    gaps.py's own gap points are already pre-filtered, before this is ever called, to sit far
    from every existing plate)."""
    tree = cKDTree(points)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        own_dist, _ = tree.query(world_pts)
        if exclude_tree is None:
            return own_dist < coverage_radius_rad
        other_dist, _ = exclude_tree.query(world_pts)
        return (own_dist < coverage_radius_rad) & (own_dist < other_dist)

    def elevation_at(world_pts: np.ndarray) -> np.ndarray:
        _, idx = tree.query(world_pts)
        return elevation[idx]

    return build_lines_from_lattice(frame, is_owned, elevation_at, spacing_rad=spacing_rad)


# --- Deformation constants (PlateWithLines.deform) ---
#
# `shift()`/`deform()` replace the old boundary.py `step_boundaries` pipeline: instead of
# classifying a boundary node as convergent/divergent/transform from the *velocity*
# decomposition `closing_rate`, deform() classifies it from *geometry* -- did this plate's
# rotated territory end up overlapping a neighbour's polygon (contested -> convergent), or
# is it in open space nobody else claims (-> divergent/rift), or is it merely close to a
# neighbour without overlapping (-> transform)? The elevation-delta rates/reaches below are
# carried over unchanged from the old model -- only the classification predicate and the
# grow/shrink node-count cap (now `D`, this step's actual max node displacement, rather than
# `closing_rate * years`) changed. See docs/simulation-model.md's "Boundary evolution"
# section for the physical reasoning behind each rate/reach.

CONVERGENT_MOUNTAIN_RATE_M_PER_MYR = 800.0
CONVERGENT_TRENCH_RATE_M_PER_MYR = 700.0
DIVERGENT_RIDGE_TARGET_M = -1500.0  # new oceanic crust at a mid-ocean ridge
DIVERGENT_RIFT_TARGET_M = -200.0  # new continental crust in a rift valley
DIVERGENT_RELAX_RATE_PER_MYR = 0.5

# Continental rifting stretches and thins the crust over a much wider zone than oceanic
# ridge spreading (which keeps FAR_THRESHOLD_RAD's narrower reach, below).
RIFT_RANGE_KM = 300.0
RIFT_RANGE_RAD = RIFT_RANGE_KM / PLANET_RADIUS_KM

# Continent-continent collision crumples a much broader belt than a plain trench/mountain
# boundary (e.g. the Himalaya/Tibetan Plateau deformation zone).
COLLISION_RANGE_KM = 400.0
COLLISION_RANGE_RAD = COLLISION_RANGE_KM / PLANET_RADIUS_KM

# Collision's second, much weaker and much farther-reaching band -- zero out to
# FAR_FIELD_COLLISION_INNER_RAD, then ramping down to zero by FAR_FIELD_COLLISION_OUTER_RAD.
FAR_FIELD_COLLISION_INNER_KM = 1000.0
FAR_FIELD_COLLISION_OUTER_KM = 3000.0
FAR_FIELD_COLLISION_INNER_RAD = FAR_FIELD_COLLISION_INNER_KM / PLANET_RADIUS_KM
FAR_FIELD_COLLISION_OUTER_RAD = FAR_FIELD_COLLISION_OUTER_KM / PLANET_RADIUS_KM
FAR_FIELD_MOUNTAIN_RATE_M_PER_MYR = 60.0

# Oceanic-under-continental subduction: the volcanic arc forms inland of the trench, not at
# it -- a band (see _band_intensity), zero at the boundary, peaking at the band's midpoint,
# zero again past the outer edge.
SUBDUCTION_ARC_INNER_KM = 100.0
SUBDUCTION_ARC_OUTER_KM = 300.0
SUBDUCTION_ARC_INNER_RAD = SUBDUCTION_ARC_INNER_KM / PLANET_RADIUS_KM
SUBDUCTION_ARC_OUTER_RAD = SUBDUCTION_ARC_OUTER_KM / PLANET_RADIUS_KM

# Transform (strike-slip) boundaries: narrower reach and gentler rate than either
# convergent case -- real motion here produces at most local pressure-ridge relief.
TRANSFORM_RANGE_KM = 50.0
TRANSFORM_RANGE_RAD = TRANSFORM_RANGE_KM / PLANET_RADIUS_KM
TRANSFORM_UPLIFT_RATE_M_PER_MYR = 200.0

# Reference (World.node_density == 1.0) values for the density-scaled thresholds below.
FAR_THRESHOLD_RAD = 1.6 * TARGET_LINE_SPACING_RAD
EXTEND_THRESHOLD_RAD = 1.3 * TARGET_LINE_SPACING_RAD
MAX_BOUNDARY_EFFECT_RAD = max(
    FAR_THRESHOLD_RAD,
    COLLISION_RANGE_RAD,
    FAR_FIELD_COLLISION_OUTER_RAD,
    SUBDUCTION_ARC_OUTER_RAD,
    TRANSFORM_RANGE_RAD,
    RIFT_RANGE_RAD,
)
# Hard safety ceiling only, not the primary cap any more -- D (this step's actual max node
# displacement, see Plate.shift) is deform()'s real physical bound on how much a line's end
# can grow/shrink in one call.
MAX_EXTEND_NODES_PER_STEP = 400

# Growth at a line end -- ordinarily plain ridge/rift fill -- instead comes back as a fresh
# volcano (guaranteed one immediate eruption, then the ordinary per-step eruption roll in
# volcanism.py takes over) with this probability per growth *event*, representing "the plate
# has been stretched too thin to keep filling with plain crust." Deliberately probabilistic
# rather than a hard threshold on any per-call quantity -- two threshold-based designs were
# tried and rejected during development:
#   - "the line's own existing gap already exceeds target spacing": dead code, confirmed
#     directly -- elevation_lines.regularize_line runs at the end of every deform() call and
#     resamples every line back to (within tolerance of) exact target spacing, so by the time
#     the *next* call's growth check runs, any such gap has already been smoothed away by the
#     *previous* call's own regularize pass.
#   - "this call is inserting at least N new nodes at once": also confirmed empirically
#     unreachable at realistic step sizes/plate rates -- sampled 1392 real growth events
#     across a running simulation and 100% of them inserted exactly 1 node, since ordinary
#     per-step divergence rarely outruns a single spacing unit's worth of growth in one call
#     regardless of how the threshold was tuned.
# A small per-event probability sidesteps needing any persistent "how long has this been
# thinning" state to track (which line/end bookkeeping would have to survive regularize,
# split, and merge) while still producing "occasionally, not constantly" volcanic crust at
# active rifts over the course of a real run -- the same shape volcanism.py's own eruption
# roll already uses for "occasional" events elsewhere in this codebase.
STRETCH_VOLCANO_PROBABILITY = 0.02

def _divergent_target(crust_type: str) -> float:
    return DIVERGENT_RIDGE_TARGET_M if crust_type == "oceanic" else DIVERGENT_RIFT_TARGET_M


def _band_intensity(dist: np.ndarray, inner: float, outer: float) -> np.ndarray:
    """Triangular profile: 0 at and outside [inner, outer], peaking at 1.0 at the band's
    midpoint -- for the subduction volcanic arc, strongest *offset* from the boundary."""
    mid = (inner + outer) / 2.0
    half_width = (outer - inner) / 2.0
    return np.clip(1.0 - np.abs(dist - mid) / half_width, 0.0, 1.0)


def _far_field_intensity(dist: np.ndarray, inner: float, outer: float) -> np.ndarray:
    """One-sided ramp: 0 below `inner`, 1.0 right at `inner`, decaying linearly to 0 by
    `outer` -- the collision far-field band, offset inland rather than continuous with the
    boundary itself."""
    ramp = np.clip(1.0 - (dist - inner) / (outer - inner), 0.0, 1.0)
    return np.where(dist < inner, 0.0, ramp)


def _far_threshold_rad(spacing_rad: float) -> float:
    return 1.6 * spacing_rad


def _extend_threshold_rad(spacing_rad: float) -> float:
    return 1.3 * spacing_rad


def _max_boundary_effect_rad(spacing_rad: float) -> float:
    return max(
        _far_threshold_rad(spacing_rad),
        COLLISION_RANGE_RAD,
        FAR_FIELD_COLLISION_OUTER_RAD,
        SUBDUCTION_ARC_OUTER_RAD,
        TRANSFORM_RANGE_RAD,
        RIFT_RANGE_RAD,
    )


def _max_extend_nodes_per_step(node_density: float) -> int:
    # A 1D count, not an area -- scales by sqrt(node_density), same reasoning as
    # MAX_EXTEND_NODES_PER_STEP's own comment.
    return max(1, round(MAX_EXTEND_NODES_PER_STEP * np.sqrt(node_density)))


class Plate(abc.ABC):
    """A plate's shared identity/motion state plus an abstract interface over however it
    represents its own terrain nodes -- `PlateWithLines` (parallel `ElevationLine`s, see
    elevation_lines.py) and `PlateWithRTree` (an R-tree-indexed point cloud, see below).
    Every method below that doesn't depend on node representation (motion, identity) is
    implemented once here; `node_count`/`all_points_and_elevation`/`outline_world`/`collect`
    are representation-specific and left abstract."""

    def __init__(
        self,
        plate_id: int,
        frame: np.ndarray,
        crust_type: str,
        omega: np.ndarray | None = None,
        age_steps: int = 0,
    ) -> None:
        self._plate_id = plate_id
        self._frame = frame
        self._crust_type = crust_type
        self._omega = omega if omega is not None else np.zeros(3)
        self._age_steps = age_steps
        # Lazily (re)computed by get_bounding_polygon() below -- None means "stale, recompute
        # on next call," not "empty polygon" (an empty plate's real outline is a valid
        # np.zeros((0, 3)), which must stay distinguishable from "not computed yet").
        # Invalidated by rotate() here and by whichever of set_lines/replace_line
        # (PlateWithLines) or set_nodes (PlateWithRTree) actually changes node positions --
        # elevation-only mutations (erosion, uplift, ...) don't touch outline_world's inputs
        # (each line's/point cloud's own theta/phi), so they leave the cache untouched.
        self._bounding_polygon_cache: np.ndarray | None = None
        # A cKDTree over that same cached outline, built lazily by get_bounding_polygon_tree()
        # below and invalidated in lockstep with it (same _invalidate_bounding_polygon call) --
        # _plates_within queries this plate's outline once per *other* plate whose own
        # get_neighbours pass considers it a candidate neighbour, so without this a busy
        # plate's outline got re-treed from scratch on every one of those incoming checks
        # this step, not just once.
        self._bounding_polygon_tree_cache: cKDTree | None = None

    @property
    def plate_id(self) -> int:
        return self._plate_id

    @property
    def frame(self) -> np.ndarray:
        """3x3 rotation matrix, local -> world. Only ever changes via `rotate`."""
        return self._frame

    @property
    def crust_type(self) -> str:
        """\"continental\" or \"oceanic\"."""
        return self._crust_type

    @property
    def omega(self) -> np.ndarray:
        """Angular velocity, world frame. Only ever changes via `set_omega`."""
        return self._omega

    @property
    def age_steps(self) -> int:
        """Steps since this plate was created (by generation, merge, or split). Gates split
        eligibility in merge_split.py so a plate can't fragment repeatedly in quick
        succession -- see the note there on why that runaway is a real failure mode."""
        return self._age_steps

    @property
    def seed_world(self) -> np.ndarray:
        """World position of this plate's local (phi=0, theta=0) reference point."""
        return self._frame[:, 0]

    def set_omega(self, omega: np.ndarray) -> None:
        self._omega = omega

    def rotate(self, increment: np.ndarray) -> None:
        """Apply an incremental rotation matrix to this plate's frame -- the one place a
        plate's rigid motion actually advances `frame` each step (see world.py)."""
        self._frame = increment @ self._frame
        self._invalidate_bounding_polygon()

    def age_one_step(self) -> None:
        self._age_steps += 1

    def reset_age(self) -> None:
        self._age_steps = 0

    @abc.abstractmethod
    def node_count(self) -> int: ...

    @abc.abstractmethod
    def all_points_and_elevation(self) -> tuple[np.ndarray, np.ndarray]:
        """Every node's world position and elevation, concatenated."""
        ...

    @abc.abstractmethod
    def outline_world(self) -> np.ndarray:
        """A live approximation of this plate's current territory outline."""
        ...

    def get_bounding_polygon(self) -> np.ndarray:
        """`outline_world()`, cached -- for any caller that doesn't need a guaranteed-fresh
        recompute (most don't: nothing about a plate's outline changes except by rotate() or
        a representation's own node-set mutation, both of which invalidate this cache
        themselves). Prefer this over calling `outline_world()` directly wherever the same
        plate's outline might reasonably be asked for more than once before its geometry
        next changes -- e.g. `get_neighbours`, run once per plate per pass, otherwise
        recomputing every other plate's outline from scratch each time (see
        `_plates_within`)."""
        if self._bounding_polygon_cache is None:
            self._bounding_polygon_cache = self.outline_world()
        return self._bounding_polygon_cache

    def get_bounding_polygon_tree(self) -> cKDTree | None:
        """A `cKDTree` over `get_bounding_polygon()`, cached the same way -- `None` if this
        plate currently has no outline (mirrors `get_bounding_polygon()`'s own empty-array
        case; a cKDTree can't be built over zero points). See `_plates_within`, the one
        caller: without this, every *other* plate's `get_neighbours` pass that considers this
        plate a candidate neighbour re-treed the same unchanged outline from scratch."""
        if self._bounding_polygon_tree_cache is None:
            polygon = self.get_bounding_polygon()
            if len(polygon) == 0:
                return None
            self._bounding_polygon_tree_cache = cKDTree(polygon)
        return self._bounding_polygon_tree_cache

    def _invalidate_bounding_polygon(self) -> None:
        self._bounding_polygon_cache = None
        self._bounding_polygon_tree_cache = None

    def contains_batch(self, points_xyz: np.ndarray) -> np.ndarray:
        """True for every point in `points_xyz` (world unit vectors) currently inside this
        plate's territory -- the batched form of `contains`, and what `deform`'s own
        contested/open classification actually calls (see `_contested_by_any`). Default
        implementation: the generic winding-number test against `get_bounding_polygon()`.
        `PlateWithLines` overrides this with a much faster row-lookup path (see there);
        `PlateWithRTree` has no equivalent row structure to exploit, so it's left on this
        default."""
        return geometry.points_in_spherical_polygon(points_xyz, self.get_bounding_polygon())

    @abc.abstractmethod
    def collect(self, field_name: str) -> np.ndarray:
        """Every node's current `field_name` value (elevation or any ElevationLine
        OPTIONAL_FIELDS name), concatenated in this plate's own node order. Empty
        (`np.zeros(0)`, or `dtype=bool` for "is_volcano") if this plate has no nodes."""
        ...

    @abc.abstractmethod
    def contains(self, lat: float, lon: float) -> bool:
        """True if the geographic point (lat, lon, radians) falls within this plate's
        current territory -- see `SpherePolygon`, which this satisfies."""
        ...

    @abc.abstractmethod
    def get_neighbours(self, all_plates: list["Plate"], threshold_rad: float = NEIGHBOUR_DISTANCE_RAD) -> list["Plate"]:
        """Every other plate in `all_plates` (this plate need not be excluded by the caller
        -- it's excluded here) whose outline comes within `threshold_rad` of this plate's
        own -- defaults to NEIGHBOUR_DISTANCE_RAD, but callers with their own notion of
        "close enough" (e.g. volcanism.py's merge/isolation distances) can pass their own."""
        ...

    @abc.abstractmethod
    def __iter__(self) -> Iterator[ElevationPoint]:
        """Every node this plate owns, as `ElevationPoint`s, in whatever order this plate's
        own representation stores them -- `PlateWithLines` line-by-line, `PlateWithRTree` in
        its flat array's own order. Lets code that just wants "every node, read or write"
        (not a bulk array op) work the same way against either representation instead of
        reaching into `PlateWithLines.lines`/`ElevationLine`'s arrays or `PlateWithRTree`'s
        own flat `_theta`/`_phi`/`_elevation`/`_fields`."""
        ...

    @abc.abstractmethod
    def map_world_points(self) -> Iterator[tuple[ElevationPoint, np.ndarray]]:
        """Every node this plate owns, paired with its own world xyz position -- the same
        nodes/order as `__iter__`, just with each `ElevationPoint` accompanied by the world
        coordinate a caller would otherwise have had to derive itself (e.g. via
        `ElevationLine.world_xyz(plate.frame)`). Each `ElevationPoint` is a live view, so a
        value computed as a function of world position (noise, distance, sampled field) can be
        written straight back with the point's own `set_*` -- in place, no
        `replace`/`replace_line`/`set_nodes` round-trip needed."""
        ...

    @abc.abstractmethod
    def map_world_points_on_plate(self) -> Iterator[tuple[ElevationPoint, np.ndarray, float]]:
        """Same as `map_world_points`, with each pair additionally carrying how far across
        this plate the node sits, normalized to [0, 1] -- 0 and 1 at the plate's own
        boundary, values in between toward the interior. `PlateWithLines` measures this
        along the node's own row: 0/1 at its `ElevationLine`'s own low/high theta endpoints,
        the same two points `outline_world()` already traces as that row's edge.
        `PlateWithRTree` has no row structure to measure along, so it approximates the same
        thing against the plate's own overall theta range instead -- cheaper than a true
        per-node distance-to-outline query, at the cost of not accounting for phi."""
        ...

    @abc.abstractmethod
    def set_fields_on_plate(self, **fields: np.ndarray) -> None:
        """Bulk in-place write for `elevation` and/or any `ElevationLine.OPTIONAL_FIELDS`
        name: each keyword's array must be exactly this plate's own node count, in the same
        order `map_world_points_on_plate`/`collect` already read/traverse it in. The
        vectorized counterpart to looping `map_world_points_on_plate` and calling each
        point's own `set_*` -- for a caller that already has a full per-node array computed
        (erosion/bathymetry/geology's per-step recompute), this writes it back without
        constructing a `Plate.__iter__`-style point object per node."""
        ...

    @abc.abstractmethod
    def update_to_lat_long_grid(self, grid: LatLongGrid) -> None:
        """Write this plate's current node elevations into `grid`'s nearest cells -- the
        plate -> grid half of the bidirectional round trip `update_deltas_from_lat_long_grid`
        (below) completes. See lat_long_grid.py's own module docstring for why the grid keeps
        this as each cell's *original* elevation rather than folding it into whatever that
        cell's own value already was."""
        ...

    @abc.abstractmethod
    def update_deltas_from_lat_long_grid(self, grid: LatLongGrid) -> None:
        """Apply `grid`'s accumulated elevation delta (see LatLongGrid.change_elevation) back
        onto this plate's own nodes -- the grid -> plate half of the round trip
        `update_to_lat_long_grid` starts. Only the *net change* since that call populated the
        grid is applied, never the grid's own absolute elevation (see lat_long_grid.py)."""
        ...

    def merge_with(
        self,
        other: "Plate",
        spacing_rad: float,
        coverage_radius_rad: float,
        other_points_xyz: np.ndarray,
    ) -> None:
        """Fuse `other` into this plate in place: keep this plate's own identity/frame, absorb
        `other`'s territory, and drop `other` (the caller still owns removing `other` from
        whatever list it lives in -- this only mutates `self`). Blending `omega` and
        resetting `age_steps` is the same regardless of representation, so it's done once
        here; folding `other`'s actual nodes into this plate's own node set is representation-
        specific (a `PlateWithLines` has to resample onto a fresh lattice to keep a fixed-row
        structure, a `PlateWithRTree` can just union the two point clouds) and left to
        `_merge_nodes_with` below.

        `spacing_rad` is this world's own lattice density (see elevation_lines.line_spacing_rad),
        `coverage_radius_rad` how far a candidate node may sit from the merging pair's old
        points and still be claimed (merge_split.MERGE_COVERAGE_RADIUS_RAD, scaled the same
        way). `other_points_xyz` is every currently-live node belonging to any *other* plate
        (not self, not `other`), world-frame xyz -- `np.zeros((0, 3))` if there are none --
        so a representation that resamples can stay clear of a still-living neighbor's own
        territory the same way plates.generate_plates' nearest-seed tiling does (see the
        merge overlap bug `merge_split.merge_plates` used to have without this)."""
        self._merge_nodes_with(other, spacing_rad, coverage_radius_rad, other_points_xyz)
        self.set_omega(mantle.clamp_rate((self.omega + other.omega) / 2.0))
        self.reset_age()

    @abc.abstractmethod
    def _merge_nodes_with(
        self,
        other: "Plate",
        spacing_rad: float,
        coverage_radius_rad: float,
        other_points_xyz: np.ndarray,
    ) -> None:
        """The representation-specific half of `merge_with`: fold `other`'s current nodes
        into this plate's own node set in place. See `merge_with` for what each parameter
        means."""
        ...

    @abc.abstractmethod
    def split(self, new_id: int, cut_normal: np.ndarray, min_nodes: int) -> tuple["Plate", "Plate"] | None:
        """Partition this plate's own nodes by which side of the great circle with normal
        `cut_normal` (world-frame unit vector) each falls on, into two new plates of this same
        representation/frame/crust_type: the `cut_normal`-positive half keeps this plate's own
        `plate_id`, the negative half gets `new_id`. `omega`/`age_steps` are left at their
        defaults on both halves -- the caller (merge_split.maybe_split_plate) fits and sets
        its own Euler pole per half afterward, from mantle flow computed over this (undivided)
        plate's own points, not from the geometric partition itself. Returns `None`, changing
        nothing, if either half would end up with fewer than `min_nodes` nodes."""
        ...

    @abc.abstractmethod
    def grow_into(
        self,
        new_points_xyz: np.ndarray,
        new_elevation: np.ndarray,
        coverage_radius_rad: float,
        spacing_rad: float,
    ) -> None:
        """Claim `new_points_xyz` (world-frame unit vectors, each with its own
        `new_elevation`) as new territory, in place -- the same "fold new nodes into an
        existing node set" shape `merge_with` uses for fusing two whole plates, just for a
        raw batch of new points instead of another `Plate` (deform()'s own "claim adjacent
        territory" sub-step is the only caller today). `coverage_radius_rad`/`spacing_rad`
        only matter to a representation that has to resample onto a fresh lattice
        (`PlateWithLines`); one that doesn't (`PlateWithRTree`) can ignore them and just
        append the new points directly. Every `ElevationLine.OPTIONAL_FIELDS` value for a
        newly-claimed node starts at its default (zero/False) -- no history to carry, same
        convention `ElevationLine.with_new_nodes` already uses."""
        ...

    @abc.abstractmethod
    def shift(self, world: "World", years: float) -> float:
        """Refit this plate's Euler pole from the mantle-flow field (damped toward the new
        target -- same fit-then-clamp as generation-time, just exponentially smoothed via
        `mantle.VELOCITY_DAMPING`), then rotate rigidly by `years` at that rate -- exact for
        every node this plate carries, no resampling. Returns `D`, the greatest angular
        distance (radians) any of this plate's own nodes actually moved this step -- the
        physical bound `deform()` must not stretch, crush, or subduct past."""
        ...

    @abc.abstractmethod
    def deform(self, world: "World", other_plates: list["Plate"], years: float, max_distance: float) -> None:
        """Reconcile this plate's actual post-`shift()` footprint against the footprint it's
        entitled to occupy -- the sphere minus every other currently-live plate's own
        bounding polygon (`other_plates`, each plate's `get_bounding_polygon()`). A node now
        geometrically inside some other plate's polygon is contested (collision/subduction:
        apply uplift/trench elevation, then shrink the affected end, consuming no more than
        `max_distance` worth of nodes -- see `shift()`). Territory nobody else claims is open
        (a rift: grow into it, or -- if already stretched thin -- spawn a fresh volcano
        instead). Everything else close to a neighbour but not actually overlapping it is
        transform. `years` sizes the elevation-delta magnitudes (uplift/trench/relax rates
        are per-Myr); `max_distance` (`D`, from `shift()`) caps how many nodes any single
        grow/shrink/claim may touch this call."""
        ...


@dataclass
class _RowLookup:
    """`PlateWithLines.contains_batch`'s cached fast-path data -- see that method's own
    docstring for the algorithm. `phis` is every line's own `phi`, sorted ascending;
    `low_thetas[i]`/`high_thetas[i]` are that same row's own theta[0]/theta[-1] (`outline_
    world`'s own "low"/"high" labels, kept index-aligned with `phis` the same way). The rest
    are all derived once here rather than per query:

    `margin_rad` -- see `_row_lookup_bulge_margin_rad`'s own docstring: a query point whose
    nearest-row interval test says "outside" isn't necessarily outside -- this bounds how
    far the *true* boundary can lie beyond the idealized per-row interval.
    `phi_min_pad`/`phi_max_pad` -- `phis[0]`/`phis[-1]`, padded outward by `margin_rad`.
    `padded_low`/`padded_high` -- `low_thetas`/`high_thetas`, each widened by `margin_rad`
    *and* by whichever of its own immediate row-neighbours (index - 1, index + 1) reaches
    further -- covers a query point landing just past a shelf-step boundary, whose relevant
    bulge belongs jointly to the two rows either side of that step, not to its own nearest
    row alone."""

    phis: np.ndarray
    low_thetas: np.ndarray
    high_thetas: np.ndarray
    margin_rad: float
    phi_min_pad: float
    phi_max_pad: float
    padded_low: np.ndarray
    padded_high: np.ndarray


def _row_lookup_bulge_margin_rad(phis: np.ndarray, low_thetas: np.ndarray, high_thetas: np.ndarray) -> float:
    """How far a `PlateWithLines` outline's real boundary can lie beyond the idealized
    "nearest row, check its own theta interval" model `contains_batch` uses as its fast path.

    `outline_world` connects consecutive vertices with straight 3D chords, not the curves an
    idealized fixed-phi/fixed-theta staircase would trace. A chord between two vertices that
    share a theta (the "vertical" step edges, phi changing) is an exact meridian -- no
    deviation. A chord between two vertices that share a *phi* instead (a row's own closing
    edge at the plate's extreme phi, or a shelf step's own "horizontal" jump between two
    adjacent rows' theta) is not: on a sphere, two equal-latitude points' connecting chord
    bulges *toward the nearer pole*, exactly the way a great-circle flight path between two
    equal-latitude cities bulges poleward. For two points at shared latitude phi separated by
    a theta half-span `d`, the chord's own peak latitude phi_mid satisfies
    `sin(phi_mid) = sin(phi) / sqrt(1 - cos(phi)^2 * sin(d)^2)` (derived by normalizing the
    two points' own vector sum, which -- for equal-magnitude equal-latitude inputs -- gives
    exactly that arc's midpoint). `phi_mid - phi` is that edge's own bulge; this returns the
    max over every such edge in this plate's outline (every shelf step between adjacent rows,
    plus the two extreme closing edges) -- a global, not per-edge, bound, so
    `contains_batch`'s per-row padding stays a simple lookup rather than tracking which
    specific edges bound which specific phi range.

    Purely additive, not proportional to plate size: verified by direct measurement (real
    plates from stepped worlds) to only ever be a small fraction of a plate's own line
    spacing, however large the plate -- see `_FAR_FIELD_PAD_RAD` in geometry.py for the same
    "additive, not multiplicative" reasoning applied to a related but distinct problem (how
    far the winding-number test itself stays reliable)."""

    def bulge(phi_edge: np.ndarray, half_span: np.ndarray) -> np.ndarray:
        denom = np.sqrt(np.clip(1.0 - np.cos(phi_edge) ** 2 * np.sin(half_span) ** 2, 1e-12, None))
        sin_phi_mid = np.clip(np.sin(phi_edge) / denom, -1.0, 1.0)
        return np.abs(np.arcsin(sin_phi_mid)) - np.abs(phi_edge)

    margins = [0.0]
    if len(phis) >= 2:
        boundary_phi = (phis[:-1] + phis[1:]) / 2.0
        margins.append(float(np.max(bulge(boundary_phi, np.abs(high_thetas[1:] - high_thetas[:-1]) / 2.0))))
        margins.append(float(np.max(bulge(boundary_phi, np.abs(low_thetas[1:] - low_thetas[:-1]) / 2.0))))
    margins.append(float(bulge(np.asarray(phis[0]), np.asarray((high_thetas[0] - low_thetas[0]) / 2.0))))
    margins.append(float(bulge(np.asarray(phis[-1]), np.asarray((high_thetas[-1] - low_thetas[-1]) / 2.0))))
    return max(margins)


class PlateWithLines(Plate):
    """A plate whose terrain is a set of parallel `ElevationLine`s at fixed plate-local
    latitudes -- see the module docstring for why this representation makes rigid rotation
    exact and resampling-free."""

    def __init__(
        self,
        plate_id: int,
        frame: np.ndarray,
        crust_type: str,
        lines: list[ElevationLine] | None = None,
        omega: np.ndarray | None = None,
        age_steps: int = 0,
    ) -> None:
        super().__init__(plate_id, frame, crust_type, omega=omega, age_steps=age_steps)
        self._lines: list[ElevationLine] = list(lines) if lines is not None else []
        # Lazily (re)built by _get_row_lookup() below, invalidated in lockstep with the
        # bounding-polygon cache (same rotate()/set_lines()/replace_line() call sites) --
        # see contains_batch's own docstring for what this backs.
        self._row_lookup_cache: _RowLookup | None = None

    @property
    def lines(self) -> tuple[ElevationLine, ...]:
        """Read-only -- use `set_lines`/`replace_line` to change this plate's lines."""
        return tuple(self._lines)

    def set_lines(self, new_lines: list[ElevationLine]) -> None:
        self._lines = list(new_lines)
        self._invalidate_bounding_polygon()

    def replace_line(self, index: int, new_line: ElevationLine) -> None:
        self._lines[index] = new_line
        self._invalidate_bounding_polygon()

    def _invalidate_bounding_polygon(self) -> None:
        super()._invalidate_bounding_polygon()
        self._row_lookup_cache = None

    def outline_world(self) -> np.ndarray:
        """Derived directly from each line's current two endpoints -- the actual edge
        deform() maintains -- rather than a separately-tracked polygon that could drift out
        of sync with the real data. Traces a *staircase*, not a smooth scanline: the
        high-theta edge across lines in ascending phi, then the low-theta edge back down,
        stepping at the midpoint phi between each pair of adjacent rows so the loop only
        ever claims territory out to whichever row's own actual extent applies on its own
        side of that midpoint.

        A straight diagonal between two rows whose theta extents differ a lot (an ordinary
        outcome once a plate's shape is no longer convex-ish -- deform() growing/shrinking
        each row's ends independently means adjacent rows routinely end up at quite
        different theta bounds) cuts across whatever concave notch sits between them,
        silently claiming sphere area this plate doesn't actually cover. That's fatal once
        `deform()` uses this same outline for its own contested/open classification --
        `PlateWithLines.deform`'s own docstring, and the invariant test in
        `unit_tests/test_plates.py`/`stress_tests/test_world_stepping.py` checking that no
        plate's own node ever falls inside a different plate's `get_bounding_polygon()`,
        both depend on this being a tight fit, not just a "reasonable envelope." Confirmed
        directly: the smooth-diagonal version this replaced put real nodes 50-120km inside a
        neighbor's polygon -- well past that neighbor's own actual nearest node -- purely
        from this concave-notch-filling effect, not genuine territory overlap."""
        lines_with_nodes = [line for line in self._lines if len(line) > 0]
        if not lines_with_nodes:
            return np.zeros((0, 3))
        ordered = sorted(lines_with_nodes, key=lambda line: line.phi)
        phis = [line.phi for line in ordered]
        high_thetas = [line[-1].get_theta() for line in ordered]
        low_thetas = [line[0].get_theta() for line in ordered]
        n = len(ordered)

        high_side: list[tuple[float, float]] = []
        for i in range(n):
            high_side.append((phis[i], high_thetas[i]))
            if i + 1 < n and high_thetas[i] != high_thetas[i + 1]:
                boundary_phi = (phis[i] + phis[i + 1]) / 2.0
                high_side.append((boundary_phi, high_thetas[i]))
                high_side.append((boundary_phi, high_thetas[i + 1]))

        low_side: list[tuple[float, float]] = []
        for i in range(n - 1, -1, -1):
            low_side.append((phis[i], low_thetas[i]))
            if i - 1 >= 0 and low_thetas[i] != low_thetas[i - 1]:
                boundary_phi = (phis[i] + phis[i - 1]) / 2.0
                low_side.append((boundary_phi, low_thetas[i]))
                low_side.append((boundary_phi, low_thetas[i - 1]))

        loop = high_side + low_side
        phi_arr = np.array([p for p, _ in loop])
        theta_arr = np.array([t for _, t in loop])
        loop_local = geometry.local_xyz(phi_arr, theta_arr)
        return geometry.to_world(self._frame, loop_local)

    def _get_row_lookup(self) -> _RowLookup | None:
        """`_RowLookup`, cached and invalidated the same way `get_bounding_polygon` is (see
        `_invalidate_bounding_polygon`) -- `None` if this plate currently has no lines."""
        if self._row_lookup_cache is not None:
            return self._row_lookup_cache
        lines_with_nodes = [line for line in self._lines if len(line) > 0]
        if not lines_with_nodes:
            return None
        ordered = sorted(lines_with_nodes, key=lambda line: line.phi)
        phis = np.array([line.phi for line in ordered])
        low_thetas = np.array([line.theta[0] for line in ordered])
        high_thetas = np.array([line.theta[-1] for line in ordered])

        margin_rad = _row_lookup_bulge_margin_rad(phis, low_thetas, high_thetas)
        if len(phis) >= 2:
            low_prev = np.concatenate([low_thetas[:1], low_thetas[:-1]])
            low_next = np.concatenate([low_thetas[1:], low_thetas[-1:]])
            high_prev = np.concatenate([high_thetas[:1], high_thetas[:-1]])
            high_next = np.concatenate([high_thetas[1:], high_thetas[-1:]])
        else:
            low_prev = low_next = low_thetas
            high_prev = high_next = high_thetas
        padded_low = np.minimum(np.minimum(low_prev, low_thetas), low_next) - margin_rad
        padded_high = np.maximum(np.maximum(high_prev, high_thetas), high_next) + margin_rad

        self._row_lookup_cache = _RowLookup(
            phis=phis,
            low_thetas=low_thetas,
            high_thetas=high_thetas,
            margin_rad=margin_rad,
            phi_min_pad=float(phis[0] - margin_rad),
            phi_max_pad=float(phis[-1] + margin_rad),
            padded_low=padded_low,
            padded_high=padded_high,
        )
        return self._row_lookup_cache

    def contains_batch(self, points_xyz: np.ndarray) -> np.ndarray:
        """Overrides `Plate.contains_batch` with an O(log rows) fast path exploiting this
        representation's own structure, exact-fallback for the rest: `outline_world`'s
        staircase is sorted by phi, so which row governs a given query phi is a `searchsorted`
        away, not a full winding-number test over every polygon vertex.

        A query point (converted to this plate's own local phi/theta) is:
        - definitely inside if its phi falls within this plate's own row range *and* its
          theta falls within its nearest row's own [low_theta, high_theta] -- verified exact
          (zero false positives across 96k+ points sampled from real captured production
          calls, plus 70k+ adversarial synthetic cases): the idealized per-row interval can
          only *underclaim* territory relative to the true outline (see
          `_row_lookup_bulge_margin_rad`'s own docstring for why), never overclaim it.
        - definitely outside if it's not even within `_RowLookup`'s own padded margin of the
          idealized region -- that margin already bounds the maximum the true boundary can
          deviate from the idealized one, so anything beyond it truly cannot be inside.
        - otherwise (idealized says outside, but within the padded margin -- a thin band that
          only matters near a plate's own boundary) exactly resolved via the same winding-
          number test `Plate.contains_batch`'s default uses, for just those points.

        Bit-exact against that same winding-number test on every real call captured from a
        multi-world, multi-step run (0 mismatches / 123k+ points) once geometry.py's own
        far-field guard (`_plausibly_near`) was fixed -- see that function's docstring for
        the unrelated pre-existing bug that surfaced during this validation."""
        n = len(points_xyz)
        if n == 0:
            return np.zeros(0, dtype=bool)
        row_lookup = self._get_row_lookup()
        if row_lookup is None:
            return np.zeros(n, dtype=bool)
        phis, low_thetas, high_thetas = row_lookup.phis, row_lookup.low_thetas, row_lookup.high_thetas

        local_xyz = geometry.to_local(self._frame, points_xyz)
        phi_q, theta_q = geometry.xyz_to_latlon(local_xyz)
        idx = np.searchsorted(phis, phi_q)
        idx_lo = np.clip(idx - 1, 0, len(phis) - 1)
        idx_hi = np.clip(idx, 0, len(phis) - 1)
        nearer_to_lo = np.abs(phi_q - phis[idx_lo]) <= np.abs(phis[idx_hi] - phi_q)
        nearest = np.where(nearer_to_lo, idx_lo, idx_hi)

        idealized_inside = (
            (phi_q >= phis[0])
            & (phi_q <= phis[-1])
            & (theta_q >= low_thetas[nearest])
            & (theta_q <= high_thetas[nearest])
        )
        maybe_boundary = (
            ~idealized_inside
            & (phi_q >= row_lookup.phi_min_pad)
            & (phi_q <= row_lookup.phi_max_pad)
            & (theta_q >= row_lookup.padded_low[nearest])
            & (theta_q <= row_lookup.padded_high[nearest])
        )
        result = idealized_inside.copy()
        if np.any(maybe_boundary):
            result[maybe_boundary] = geometry.points_in_spherical_polygon(
                points_xyz[maybe_boundary], self.get_bounding_polygon()
            )
        return result

    def node_count(self) -> int:
        return sum(len(line) for line in self._lines)

    def all_points_and_elevation(self) -> tuple[np.ndarray, np.ndarray]:
        """Every elevation-line node's world position and elevation, concatenated."""
        if not self._lines:
            return np.zeros((0, 3)), np.zeros(0)
        points = np.concatenate([line.world_xyz(self._frame) for line in self._lines], axis=0)
        elevation = np.concatenate([line.elevation for line in self._lines], axis=0)
        return points, elevation

    def collect(self, field_name: str) -> np.ndarray:
        chunks = [getattr(line, field_name) for line in self._lines if len(line) > 0]
        if not chunks:
            return np.zeros(0, dtype=bool) if field_name == "is_volcano" else np.zeros(0)
        return np.concatenate(chunks, axis=0)

    def contains(self, lat: float, lon: float) -> bool:
        point_xyz = geometry.latlon_to_xyz(np.asarray(lat), np.asarray(lon))
        return geometry.point_in_spherical_polygon(point_xyz, self.get_bounding_polygon())

    def get_neighbours(self, all_plates: list["Plate"], threshold_rad: float = NEIGHBOUR_DISTANCE_RAD) -> list["Plate"]:
        return _plates_within(self, all_plates, threshold_rad)

    def __iter__(self) -> Iterator[ElevationPoint]:
        for line in self._lines:
            yield from line

    def map_world_points(self) -> Iterator[tuple[ElevationPoint, np.ndarray]]:
        for line in self._lines:
            if len(line) == 0:
                continue
            world_pts = line.world_xyz(self._frame)
            for point, world_xyz in zip(line, world_pts):
                yield point, world_xyz

    def map_world_points_on_plate(self) -> Iterator[tuple[ElevationPoint, np.ndarray, float]]:
        for line in self._lines:
            if len(line) == 0:
                continue
            world_pts = line.world_xyz(self._frame)
            low_theta = line.theta[0]
            span = line.theta[-1] - low_theta
            for point, world_xyz in zip(line, world_pts):
                fraction = 0.5 if span == 0 else float((point.get_theta() - low_theta) / span)
                yield point, world_xyz, fraction

    def set_fields_on_plate(self, **fields: np.ndarray) -> None:
        offset = 0
        for line in self._lines:
            n = len(line)
            if n == 0:
                continue
            line.set_fields(**{name: values[offset : offset + n] for name, values in fields.items()})
            offset += n

    def update_to_lat_long_grid(self, grid: LatLongGrid) -> None:
        _update_to_lat_long_grid(self, grid)

    def update_deltas_from_lat_long_grid(self, grid: LatLongGrid) -> None:
        _update_deltas_from_lat_long_grid(self, grid)

    def _merge_nodes_with(
        self,
        other: "Plate",
        spacing_rad: float,
        coverage_radius_rad: float,
        other_points_xyz: np.ndarray,
    ) -> None:
        """A one-time resample onto a fresh local lattice, elevation carried over by nearest-
        neighbor lookup into the pre-merge combined point cloud -- see the module docstring
        for why that's an acceptable cost for a merge (rare, discrete) versus routine per-step
        motion. Only `elevation` survives the resample; every `ElevationLine.OPTIONAL_FIELDS`
        value (rivers, lakes, soil, ...) resets to its default on the newly-built lines, the
        same tradeoff this had before being made representation-generic.

        Exclusivity against every other *live* plate (`other_points_xyz`) is enforced the same
        way initial generation guarantees it (see `generate_plates`' nearest-seed Voronoi
        query): a candidate lattice point is only kept if it's both within
        `coverage_radius_rad` of the merged pair's own old points *and* strictly closer to
        those old points than to any other plate's current nearest node -- without the second
        half of that check, a candidate near one of the merged pair's own (possibly scattered,
        if either parent had already been through an earlier merge) old points would be
        claimed even inside a completely unrelated, still-living plate's own territory."""
        keep_pts, keep_elev = self.all_points_and_elevation()
        absorb_pts, absorb_elev = other.all_points_and_elevation()
        old_points = np.concatenate([keep_pts, absorb_pts], axis=0)
        old_elevation = np.concatenate([keep_elev, absorb_elev], axis=0)
        exclude_tree = cKDTree(other_points_xyz) if len(other_points_xyz) else None
        self.set_lines(
            _lines_from_resample(self.frame, old_points, old_elevation, coverage_radius_rad, spacing_rad, exclude_tree)
        )

    def grow_into(
        self,
        new_points_xyz: np.ndarray,
        new_elevation: np.ndarray,
        coverage_radius_rad: float,
        spacing_rad: float,
    ) -> None:
        """A one-time resample onto a fresh local lattice, same shape (and same
        elevation-only-survives tradeoff) as `_merge_nodes_with`, just folding in a raw batch
        of new points instead of another plate's -- no exclusivity check against other plates
        needed here (unlike a merge), since the caller (gaps.py) has already restricted
        `new_points_xyz` to territory no live plate currently covers."""
        old_points, old_elevation = self.all_points_and_elevation()
        combined_points = np.concatenate([old_points, new_points_xyz], axis=0)
        combined_elevation = np.concatenate([old_elevation, new_elevation], axis=0)
        self.set_lines(_lines_from_resample(self.frame, combined_points, combined_elevation, coverage_radius_rad, spacing_rad))

    def split(self, new_id: int, cut_normal: np.ndarray, min_nodes: int) -> tuple["Plate", "Plate"] | None:
        """Cuts along existing node data rather than resampling -- unlike `merge_with`, a
        split doesn't need a fresh lattice, just a partition of each line's own nodes by which
        side of `cut_normal` they fall on (`ElevationLine.masked`), so every
        `OPTIONAL_FIELDS` value survives exactly, not just `elevation`."""
        lines_a: list[ElevationLine] = []
        lines_b: list[ElevationLine] = []
        for line in self._lines:
            world_pts = line.world_xyz(self._frame)
            side = np.sum(world_pts * cut_normal, axis=-1) > 0
            if np.any(side):
                lines_a.append(line.masked(side))
            if np.any(~side):
                lines_b.append(line.masked(~side))

        if sum(len(l) for l in lines_a) < min_nodes or sum(len(l) for l in lines_b) < min_nodes:
            return None

        plate_a = PlateWithLines(plate_id=self.plate_id, frame=self._frame.copy(), crust_type=self._crust_type, lines=lines_a)
        plate_b = PlateWithLines(plate_id=new_id, frame=self._frame.copy(), crust_type=self._crust_type, lines=lines_b)
        return plate_a, plate_b

    def shift(self, world: "World", years: float) -> float:
        old_points, _ = self.all_points_and_elevation()

        if len(old_points) > 0:
            velocities = mantle.flow_at(old_points, world.mantle_centers)
            target_omega = mantle.fit_euler_pole(old_points, velocities)
            new_omega = self.omega + mantle.VELOCITY_DAMPING * (target_omega - self.omega)
            self.set_omega(mantle.clamp_rate(new_omega))

        increment = geometry.rotation_matrix_from_omega(self.omega, years)
        self.rotate(increment)

        if len(old_points) == 0:
            return 0.0
        new_points, _ = self.all_points_and_elevation()
        return float(geometry.angular_distance(old_points, new_points).max())

    def deform(self, world: "World", other_plates: list["Plate"], years: float, max_distance: float) -> None:
        spacing_rad = line_spacing_rad(world.node_density)
        far_threshold_rad = _far_threshold_rad(spacing_rad)
        extend_threshold_rad = _extend_threshold_rad(spacing_rad)
        max_boundary_effect_rad = _max_boundary_effect_rad(spacing_rad)
        max_extend_nodes = _max_extend_nodes_per_step(world.node_density)

        own_points, _ = self.all_points_and_elevation()
        if not self._lines or len(own_points) == 0:
            return

        neighbours = self.get_neighbours(other_plates, threshold_rad=max_boundary_effect_rad)

        if neighbours:
            pieces = [p.all_points_and_elevation()[0] for p in neighbours]
            owners = [np.full(len(pts), p.plate_id) for p, pts in zip(neighbours, pieces)]
            neighbour_points = np.concatenate(pieces, axis=0)
            neighbour_owner = np.concatenate(owners, axis=0)
        else:
            neighbour_points = np.zeros((0, 3))
            neighbour_owner = np.zeros(0, dtype=int)

        if len(neighbour_points) > 0:
            tree = cKDTree(neighbour_points, balanced_tree=False, compact_nodes=False)
            dist_all, idx_all = tree.query(own_points, workers=query_workers(len(own_points)))
            neighbor_owner_all = neighbour_owner[idx_all]
        else:
            dist_all = np.full(len(own_points), np.inf)
            neighbor_owner_all = np.zeros(len(own_points), dtype=int)

        neighbour_by_id = {p.plate_id: p for p in neighbours}

        # Polygon-containment classification -- the one place this replaces the old
        # velocity-based closing_rate. Cheap prefilter first (near_mask, from the k-d tree
        # distance already computed above): a node far from every neighbour can never be
        # contested, so the more expensive per-point polygon test only runs where it could
        # actually matter. Deliberately far_threshold_rad here, not max_boundary_effect_rad
        # (which reaches out to FAR_FIELD_COLLISION_OUTER_RAD, 3000km, for the far-field
        # mountain-uplift *intensity* curve) -- genuine polygon overlap is a boundary-local
        # phenomenon; a node thousands of km from its own nearest cross-plate node is never
        # actually going to land inside that neighbour's real territory, so testing it would
        # only spend the expensive per-point polygon check on nodes that can never come back
        # contested. Confirmed directly: this was the dominant per-step cost at realistic
        # node counts (a 10-plate, default-density step_world call went from ~40s to well
        # under a second after narrowing this).
        contested_all = np.zeros(len(own_points), dtype=bool)
        near_mask = dist_all < far_threshold_rad
        near_points = own_points[near_mask]
        if len(near_points) > 0:
            near_contested = _contested_by_any(near_points, neighbours)
            contested_all[near_mask] = near_contested

        default_intensity_all = np.clip(1.0 - dist_all / far_threshold_rad, 0.0, 1.0)
        collision_intensity_all = np.clip(1.0 - dist_all / COLLISION_RANGE_RAD, 0.0, 1.0)
        far_field_intensity_all = _far_field_intensity(dist_all, FAR_FIELD_COLLISION_INNER_RAD, FAR_FIELD_COLLISION_OUTER_RAD)
        arc_intensity_all = _band_intensity(dist_all, SUBDUCTION_ARC_INNER_RAD, SUBDUCTION_ARC_OUTER_RAD)
        transform_intensity_all = np.clip(1.0 - dist_all / TRANSFORM_RANGE_RAD, 0.0, 1.0)
        rift_intensity_all = np.clip(1.0 - dist_all / RIFT_RANGE_RAD, 0.0, 1.0)

        neighbor_is_oceanic_all = np.array(
            [neighbour_by_id[pid].crust_type == "oceanic" if pid in neighbour_by_id else False for pid in neighbor_owner_all]
        )

        convergent_all = contested_all
        # Distance-only partition of the remaining (uncontested) near-boundary nodes: within
        # TRANSFORM_RANGE_RAD's narrow reach counts as transform, the wider band beyond it
        # (out to RIFT_RANGE_RAD/FAR_THRESHOLD_RAD) as divergent relaxation. The old model
        # kept these mutually exclusive via closing_rate's sign; without a continuous
        # velocity signal, distance is the next-best proxy -- a deliberate simplification,
        # not a physical claim that transform never grades into divergence.
        transform_all = ~contested_all & (dist_all < TRANSFORM_RANGE_RAD)
        wide_reach = RIFT_RANGE_RAD if self.crust_type == "continental" else FAR_THRESHOLD_RAD
        divergent_all = ~contested_all & ~transform_all & (dist_all < wide_reach)
        subduction_all = convergent_all & neighbor_is_oceanic_all
        collision_all = convergent_all & ~neighbor_is_oceanic_all

        target = _divergent_target(self.crust_type)
        years_myr = years / 1_000_000.0
        relax_factor = 1.0 - np.exp(-DIVERGENT_RELAX_RATE_PER_MYR * years_myr)

        new_lines: list[ElevationLine] = []
        offset = 0
        for line_index, line in enumerate(self._lines):
            n = len(line)
            dist = dist_all[offset : offset + n]
            default_intensity = default_intensity_all[offset : offset + n]
            collision_intensity = collision_intensity_all[offset : offset + n]
            far_field_intensity = far_field_intensity_all[offset : offset + n]
            arc_intensity = arc_intensity_all[offset : offset + n]
            transform_intensity = transform_intensity_all[offset : offset + n]
            rift_intensity = rift_intensity_all[offset : offset + n]
            contested = contested_all[offset : offset + n]
            transform = transform_all[offset : offset + n]
            divergent = divergent_all[offset : offset + n]
            subduction = subduction_all[offset : offset + n]
            collision = collision_all[offset : offset + n]
            offset += n

            elevation = line.elevation.copy()
            if self.crust_type == "continental":
                elevation[subduction] += CONVERGENT_MOUNTAIN_RATE_M_PER_MYR * years_myr * arc_intensity[subduction]
                elevation[collision] += CONVERGENT_MOUNTAIN_RATE_M_PER_MYR * years_myr * collision_intensity[collision]
                elevation[collision] += FAR_FIELD_MOUNTAIN_RATE_M_PER_MYR * years_myr * far_field_intensity[collision]
            else:
                elevation[contested] -= CONVERGENT_TRENCH_RATE_M_PER_MYR * years_myr * default_intensity[contested]

            elevation[transform] += TRANSFORM_UPLIFT_RATE_M_PER_MYR * years_myr * transform_intensity[transform]

            divergent_intensity = rift_intensity if self.crust_type == "continental" else default_intensity
            elevation[divergent] += (target - elevation[divergent]) * relax_factor * divergent_intensity[divergent]

            elevation = np.clip(elevation, MIN_ELEVATION_M, MAX_ELEVATION_M)
            updated_line = line.replace(elevation=elevation)
            grown_line = self._grow_or_shrink_line_for_deform(
                updated_line,
                dist,
                contested,
                spacing_rad,
                extend_threshold_rad,
                max_extend_nodes,
                max_distance,
                world,
                line_index,
                neighbours,
            )
            if len(grown_line) > 0:
                new_lines.append(grown_line)

        self.set_lines(new_lines)
        self._claim_adjacent_territory(world, neighbours, spacing_rad)

        for line_index, line in enumerate(self._lines):
            if needs_regularizing(line, spacing_rad):
                self.replace_line(line_index, regularize_line(line, spacing_rad))

    def _count_open_prefix(self, theta_candidates: np.ndarray, phi: float, neighbours: list["Plate"]) -> int:
        """How many of `theta_candidates` (in order, starting closest to the existing edge)
        are NOT contested by any neighbour -- growth stops at the first candidate that
        would land inside a neighbour's own current territory, rather than blindly
        inserting every node a plain distance estimate suggested (the nearest neighbour
        *node* can be much closer or farther than the actual polygon boundary in the
        specific direction growth is extending)."""
        if len(theta_candidates) == 0 or not neighbours:
            return len(theta_candidates)
        world_pts = geometry.to_world(self.frame, geometry.local_xyz(np.full_like(theta_candidates, phi), theta_candidates))
        contested = _contested_by_any(world_pts, neighbours)
        first_contested = np.argmax(contested) if np.any(contested) else len(contested)
        return int(first_contested)

    def _grow_or_shrink_line_for_deform(
        self,
        line: ElevationLine,
        dist: np.ndarray,
        contested: np.ndarray,
        spacing_rad: float,
        extend_threshold_rad: float,
        max_extend_nodes: int,
        max_distance: float,
        world: "World",
        line_index: int,
        neighbours: list["Plate"],
    ) -> ElevationLine:
        """Shrink `line`'s two ends by however many *consecutive* contested nodes sit
        there, then grow whichever end is left both uncontested and far from any
        neighbour -- the `deform()` counterpart to the old `boundary._grow_or_shrink_line`.

        Deliberately end-only, not "remove any contested node anywhere in the line": every
        other piece of this codebase that touches an `ElevationLine` (outline_world's own
        polygon trace, `elevation_lines.regularize_line`'s endpoint-preserving resample, the
        old `boundary._grow_or_shrink_line` this replaces) assumes a line's node set is one
        *contiguous* span of territory at its own phi -- removing a node stranded in the
        middle would puncture a hole outline_world() has no way to represent (it only reads
        each line's own first/last theta), so the hole would keep reading as still-claimed
        territory. Confirmed directly as a real bug during development: removing interior
        contested nodes made the "no plate's node sits inside a neighbour's polygon"
        invariant *worse*, not better, since the resulting holes over-claimed exactly where
        they'd just been hollowed out. An interior-only contested patch (rare -- it needs a
        neighbour's own growth to have reached past this row's two ends without yet
        registering at either) is left alone for this call; the neighbour's own continued
        growth in subsequent turns reaches this row's nearer end before long, at which point
        the ordinary end-shrink below picks it up.

        Node-count caps: `max_distance` (`D`, this step's actual max node displacement -- see
        `Plate.shift`) and the hard safety ceiling `max_extend_nodes`."""
        theta = line.theta.copy()
        elevation = line.elevation.copy()
        contested = contested.copy()
        dist = dist.copy()
        persistent_fields = {name: getattr(line, name).copy() for name in ElevationLine.OPTIONAL_FIELDS}
        if len(theta) == 0:
            return ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)

        dtheta = spacing_rad / max(np.cos(line.phi), 1e-3)
        target = _divergent_target(self.crust_type)
        n_distance_cap = max(1, int(max_distance / spacing_rad))

        def contested_run_from_end(mask: np.ndarray, from_high: bool) -> int:
            ordered = mask[::-1] if from_high else mask
            run = 0
            for value in ordered:
                if not value:
                    break
                run += 1
            return run

        # High end first so the low-end index (0) is unaffected by any change made here.
        if contested[-1]:
            n_remove = min(contested_run_from_end(contested, from_high=True), n_distance_cap, max_extend_nodes, len(theta) - 1)
            if n_remove > 0:
                theta = theta[:-n_remove]
                elevation = elevation[:-n_remove]
                contested = contested[:-n_remove]
                dist = dist[:-n_remove]
                persistent_fields = {name: values[:-n_remove] for name, values in persistent_fields.items()}

        if len(theta) == 0:
            return ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)

        if contested[0]:
            n_remove = min(contested_run_from_end(contested, from_high=False), n_distance_cap, max_extend_nodes, len(theta) - 1)
            if n_remove > 0:
                theta = theta[n_remove:]
                elevation = elevation[n_remove:]
                contested = contested[n_remove:]
                dist = dist[n_remove:]
                persistent_fields = {name: values[n_remove:] for name, values in persistent_fields.items()}

        if len(theta) == 0:
            return ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)

        def grow_end(n_new: int, end_tag: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            rng = np.random.default_rng((world.seed, round(world.elapsed_years), self.plate_id, line_index, end_tag))
            overstretched = rng.random() < STRETCH_VOLCANO_PROBABILITY
            if not overstretched:
                return np.full(n_new, target), np.zeros(n_new, dtype=bool), np.zeros(n_new)
            elev = np.full(n_new, target) + ERUPTION_ELEVATION_M  # guaranteed immediate first eruption
            is_volcano = np.ones(n_new, dtype=bool)
            remaining = rng.uniform(VOLCANO_ACTIVE_MIN_YEARS, VOLCANO_ACTIVE_MAX_YEARS, size=n_new)
            return elev, is_volcano, remaining

        # High end first so the low-end index (0) is unaffected by any change made here.
        # dist can be +inf (no neighbour anywhere -- e.g. a genuinely isolated plate/field),
        # so the gap-derived candidate count is computed against a finite stand-in distance
        # before dividing; the real cap either way is n_distance_cap/max_extend_nodes.
        if not contested[-1] and dist[-1] > extend_threshold_rad:
            gap_estimate = min(dist[-1], (n_distance_cap + 1) * spacing_rad)
            n_candidates = min(max(int(gap_estimate / spacing_rad), 1), n_distance_cap, max_extend_nodes)
            candidate_theta = theta[-1] + dtheta * np.arange(1, n_candidates + 1)
            n_new = self._count_open_prefix(candidate_theta, line.phi, neighbours)
            if n_new > 0:
                new_theta = candidate_theta[:n_new]
                new_elevation, new_is_volcano, new_remaining = grow_end(n_new, 0)
                theta = np.append(theta, new_theta)
                elevation = np.append(elevation, new_elevation)
                for name, values in persistent_fields.items():
                    if name == "is_volcano":
                        fill = new_is_volcano
                    elif name == "volcano_active_years_remaining":
                        fill = new_remaining
                    else:
                        fill = np.zeros(n_new, dtype=values.dtype)
                    persistent_fields[name] = np.append(values, fill)

        if not contested[0] and dist[0] > extend_threshold_rad:
            gap_estimate = min(dist[0], (n_distance_cap + 1) * spacing_rad)
            n_candidates = min(max(int(gap_estimate / spacing_rad), 1), n_distance_cap, max_extend_nodes)
            candidate_theta = theta[0] - dtheta * np.arange(1, n_candidates + 1)
            n_new = self._count_open_prefix(candidate_theta, line.phi, neighbours)
            if n_new > 0:
                new_theta = candidate_theta[:n_new][::-1]
                new_elevation, new_is_volcano, new_remaining = grow_end(n_new, 1)
                theta = np.insert(theta, 0, new_theta)
                elevation = np.insert(elevation, 0, new_elevation)
                for name, values in persistent_fields.items():
                    if name == "is_volcano":
                        fill = new_is_volcano
                    elif name == "volcano_active_years_remaining":
                        fill = new_remaining
                    else:
                        fill = np.zeros(n_new, dtype=values.dtype)
                    persistent_fields[name] = np.insert(values, 0, fill)

        return ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)

    def _claim_adjacent_territory(self, world: "World", neighbours: list["Plate"], spacing_rad: float) -> None:
        """Claim a whole new phi row just beyond this plate's current phi extremes, where
        that row is open (unclaimed) territory -- the one case ordinary per-line end-growth
        structurally can't reach, since growth only ever extends an *existing* line's own
        theta range, never adds a brand new line. This is what actually lets a plate grow
        toward its own pole (the old gaps.py's other role -- reclaiming ground a
        subducted neighbour vacated *within* an existing row's own theta range -- doesn't
        need a separate mechanism at all: the very next time that row's end-growth check
        runs, the vacated neighbour is simply gone from `dist`/`contested`, and ordinary
        end-growth already extends into it).

        Deliberately NOT `Plate.grow_into` (a full lattice resample of the *entire* plate,
        representation-generic but rare-event-priced -- merge/absorb-scale, not "every plate,
        every turn"): confirmed directly during development that calling it here every turn
        made a plate balloon by several times its own size in one call, since the resample's
        own coverage radius around a handful of newly-claimed ring points reconstructs far
        more lattice area than just those points. Building the new row directly, the same
        way `_grow_or_shrink_line_for_deform`'s own growth builds new nodes, keeps this the
        same gradual, bounded-per-turn shape as every other change in this model."""
        lines_with_nodes = [line for line in self._lines if len(line) > 0]
        if not lines_with_nodes:
            return
        ordered = sorted(lines_with_nodes, key=lambda line: line.phi)
        max_phi_limit = np.pi / 2 - spacing_rad / 2
        base = base_elevation(self.crust_type)
        amp = noise_amplitude(self.crust_type)
        new_lines: list[ElevationLine] = []

        for reference, direction in ((ordered[0], -1), (ordered[-1], 1)):
            new_phi = reference.phi + direction * spacing_rad
            if abs(new_phi) > max_phi_limit:
                continue
            dtheta = spacing_rad / max(np.cos(new_phi), 1e-3)
            span = reference.theta[-1] - reference.theta[0]
            n_cols = max(int(round(span / dtheta)) + 1, 1)
            theta_candidates = reference.theta[0] + dtheta * np.arange(n_cols)
            world_pts = geometry.to_world(self.frame, geometry.local_xyz(np.full(n_cols, new_phi), theta_candidates))

            contested = _contested_by_any(world_pts, neighbours)
            open_mask = ~contested
            if not np.any(open_mask):
                continue

            direction_tag = 0 if direction < 0 else 1
            rng = np.random.default_rng((world.seed, round(world.elapsed_years), self.plate_id, direction_tag))
            noise = SphereNoise(rng, octaves=3, base_freq=2.5)
            theta_open = theta_candidates[open_mask]
            elevation_open = base + amp * noise.sample(world_pts[open_mask])
            new_lines.append(ElevationLine(phi=new_phi, theta=theta_open, elevation=elevation_open))

        if new_lines:
            self.set_lines(list(self._lines) + new_lines)


# outline_world's boundary-detection pass (see PlateWithRTree below) needs at least this many
# nodes for "boundary node" to be a meaningful distinct subset of "every node" -- below it,
# every node is returned as-is rather than running (and likely degenerating) the density
# check.
OUTLINE_MIN_NODES_FOR_HULL = 4
# Bounds the cost of _typical_spacing's own nearest-neighbor sampling for a very large plate
# -- a plain median over a few hundred samples is already a stable estimate of "how far
# apart nodes normally sit" without needing to query every single node.
OUTLINE_SPACING_SAMPLE_SIZE = 200
# Half-width of each node's own neighborhood box, in multiples of the plate's typical
# spacing -- wide enough that an interior node's box reliably contains several rings of
# neighbors (not just its single nearest one, which would be too noisy a density signal).
OUTLINE_NEIGHBORHOOD_SPACING_FACTOR = 3.0
# A node counts as "boundary" once its own neighborhood box holds no more than this fraction
# of the densest node's own count -- relative to the plate's own densest node (not an
# absolute count) so this works the same whether the plate as a whole is sparse or dense.
OUTLINE_BOUNDARY_DENSITY_FRACTION = 0.6


@install_point_field_accessors
class ElevationPointInCloud:
    """One node of a `PlateWithRTree`'s flat point cloud: a pointer to the plate plus its
    index into that plate's own flat `_theta`/`_phi`/`_elevation`/`_fields` arrays -- the
    `PlateWithRTree` counterpart to `ElevationPointOnLine`, sharing the same field-list-driven
    `get_*`/`set_*` machinery (see `install_point_field_accessors`) so both `ElevationPoint`
    implementations expose an identical surface despite backing storage that doesn't otherwise
    look alike (fixed-phi rows vs. an unstructured cloud). Like `ElevationPointOnLine`, a live
    view, not a snapshot -- stale once `set_nodes` rebuilds the plate's arrays out from under
    it."""

    def __init__(self, plate: "PlateWithRTree", index: int) -> None:
        n = plate.node_count()
        if not -n <= index < n:
            raise IndexError(f"PlateWithRTree point index {index} out of range for length {n}")
        self._plate = plate
        self._index = index % n

    @property
    def plate(self) -> "PlateWithRTree":
        return self._plate

    @property
    def index(self) -> int:
        """This point's (always non-negative) index into `plate`'s own flat arrays."""
        return self._index

    @property
    def phi(self) -> float:
        return float(self._plate.phi[self._index])

    def _field_array(self, name: str) -> np.ndarray:
        return self._plate.field_array(name)


class PlateWithRTree(Plate):
    """A plate whose terrain is a flat, unstructured cloud of nodes at plate-local (phi,
    theta) coordinates, spatially indexed by an `RTree` (see rtree_index.py) rather than
    organized into `PlateWithLines`' fixed-phi rows. Where `PlateWithLines` gets exact,
    resampling-free rotation from a grid that never needs reordering, this representation
    instead gets O(log n) box/nearest-neighbor spatial queries over an arbitrary point set --
    useful for anything that wants to ask "what's near this point" without caring what row,
    if any, it's on. The index is a static, bulk-loaded structure (see RTree.build) rebuilt
    from scratch by `set_nodes` whenever the node set changes, rather than incrementally
    maintained -- see rtree_index.py's own module docstring for why that's the right
    tradeoff here.

    Carries the same per-node fields as `ElevationLine` (elevation plus every
    `ElevationLine.OPTIONAL_FIELDS` name) so it satisfies the same `Plate.collect` contract
    -- just as flat parallel arrays over every node at once, rather than one line at a time."""

    def __init__(
        self,
        plate_id: int,
        frame: np.ndarray,
        crust_type: str,
        theta: np.ndarray | None = None,
        phi: np.ndarray | None = None,
        elevation: np.ndarray | None = None,
        omega: np.ndarray | None = None,
        age_steps: int = 0,
        **fields: np.ndarray,
    ) -> None:
        super().__init__(plate_id, frame, crust_type, omega=omega, age_steps=age_steps)
        self.set_nodes(
            theta if theta is not None else np.zeros(0),
            phi if phi is not None else np.zeros(0),
            elevation if elevation is not None else np.zeros(0),
            **fields,
        )

    def set_nodes(self, theta: np.ndarray, phi: np.ndarray, elevation: np.ndarray, **fields: np.ndarray) -> None:
        """Replace every node at once and rebuild the spatial index -- there's no fixed row
        structure to preserve here (unlike PlateWithLines.replace_line), so a full swap is
        the natural granularity for this representation. `fields` accepts any of
        ElevationLine.OPTIONAL_FIELDS by name; anything not passed defaults the same way
        ElevationLine itself does (zeros, or False for is_volcano)."""
        self._theta = theta
        self._phi = phi
        self._elevation = elevation
        self._fields = {
            name: fields[name] if name in fields else self._default_field(name, theta)
            for name in ElevationLine.OPTIONAL_FIELDS
        }
        local_xy = np.stack([theta, phi], axis=1) if len(theta) else np.zeros((0, 2))
        self._rtree = RTree.build(local_xy)
        self._invalidate_bounding_polygon()

    @staticmethod
    def _default_field(name: str, theta: np.ndarray) -> np.ndarray:
        return np.zeros_like(theta, dtype=bool) if name == "is_volcano" else np.zeros_like(theta)

    @property
    def theta(self) -> np.ndarray:
        return self._theta

    @property
    def phi(self) -> np.ndarray:
        return self._phi

    @property
    def elevation(self) -> np.ndarray:
        return self._elevation

    @property
    def rtree(self) -> RTree:
        """The current spatial index over this plate's own (theta, phi) nodes -- exposed for
        callers that want their own box/nearest-neighbor queries, beyond just this class's
        own outline_world use of it below."""
        return self._rtree

    def field_array(self, name: str) -> np.ndarray:
        """This plate's own flat array backing `name` ("theta", "elevation", or any
        `ElevationLine.OPTIONAL_FIELDS` name) -- what `ElevationPointInCloud` indexes into,
        the flat-array counterpart to `ElevationLine`'s own `_<name>` attributes."""
        if name == "theta":
            return self._theta
        if name == "elevation":
            return self._elevation
        return self._fields[name]

    def node_count(self) -> int:
        return len(self._theta)

    def all_points_and_elevation(self) -> tuple[np.ndarray, np.ndarray]:
        if len(self._theta) == 0:
            return np.zeros((0, 3)), np.zeros(0)
        local = geometry.local_xyz(self._phi, self._theta)
        return geometry.to_world(self._frame, local), self._elevation

    def collect(self, field_name: str) -> np.ndarray:
        if len(self._theta) == 0:
            return np.zeros(0, dtype=bool) if field_name == "is_volcano" else np.zeros(0)
        return self._elevation if field_name == "elevation" else self._fields[field_name]

    def contains(self, lat: float, lon: float) -> bool:
        point_xyz = geometry.latlon_to_xyz(np.asarray(lat), np.asarray(lon))
        return geometry.point_in_spherical_polygon(point_xyz, self.get_bounding_polygon())

    def get_neighbours(self, all_plates: list["Plate"], threshold_rad: float = NEIGHBOUR_DISTANCE_RAD) -> list["Plate"]:
        return _plates_within(self, all_plates, threshold_rad)

    def __iter__(self) -> Iterator[ElevationPoint]:
        for i in range(len(self._theta)):
            yield ElevationPointInCloud(self, i)

    def map_world_points(self) -> Iterator[tuple[ElevationPoint, np.ndarray]]:
        if len(self._theta) == 0:
            return
        world_pts = geometry.to_world(self._frame, geometry.local_xyz(self._phi, self._theta))
        for i, world_xyz in enumerate(world_pts):
            yield ElevationPointInCloud(self, i), world_xyz

    def map_world_points_on_plate(self) -> Iterator[tuple[ElevationPoint, np.ndarray, float]]:
        """No row structure to measure a node's own edge-to-edge position along (unlike
        PlateWithLines), so this approximates it against the plate's overall theta range --
        an O(n) min/max reduction, not a per-node spatial query -- rather than something more
        faithful like distance-to-outline_world(), which would cost a hull/point-in-polygon
        query per node."""
        if len(self._theta) == 0:
            return
        world_pts = geometry.to_world(self._frame, geometry.local_xyz(self._phi, self._theta))
        low_theta = float(self._theta.min())
        span = float(self._theta.max()) - low_theta
        for i, world_xyz in enumerate(world_pts):
            fraction = 0.5 if span == 0 else (float(self._theta[i]) - low_theta) / span
            yield ElevationPointInCloud(self, i), world_xyz, fraction

    def set_fields_on_plate(self, **fields: np.ndarray) -> None:
        for name, values in fields.items():
            self.field_array(name)[:] = values

    def _typical_spacing(self) -> float:
        """Median nearest-neighbor distance over a bounded sample of this plate's own nodes
        -- "how far apart nodes normally sit," the same role TARGET_LINE_SPACING_RAD plays
        for PlateWithLines, just measured directly from the actual node cloud rather than
        assumed from a fixed generation-time constant (nothing here guarantees this
        representation's own nodes came from that same lattice)."""
        n = len(self._theta)
        local_xy = np.stack([self._theta, self._phi], axis=1)
        sample_idx = (
            np.arange(n)
            if n <= OUTLINE_SPACING_SAMPLE_SIZE
            else np.linspace(0, n - 1, OUTLINE_SPACING_SAMPLE_SIZE).astype(int)
        )
        dists = [
            result[1]
            for i in sample_idx
            if (result := self._rtree.nearest_one(local_xy[i], exclude_index=int(i))) is not None
        ]
        return float(np.median(dists)) if dists else 1.0

    def outline_world(self) -> np.ndarray:
        """A live approximation of this plate's territory outline, built from the R-tree
        rather than PlateWithLines' fixed-row scanline: first finds *boundary* nodes (ones
        whose own local neighborhood -- an axis-aligned box queried via the R-tree -- holds
        fewer other nodes than an interior node's would, since an edge node's box is partly
        empty space outside the plate while an interior one's isn't), then returns the 2D
        convex hull of just those nodes' local (theta, phi) coordinates, mapped to world
        space. Restricting the hull to boundary nodes rather than every node is the same
        "hull of a full point set == hull of its own boundary" shortcut
        plate_bounding_ellipse's own docstring relies on, just applied here to build the
        outline itself rather than to fit an ellipse around it."""
        n = len(self._theta)
        if n == 0:
            return np.zeros((0, 3))
        local_xy = np.stack([self._theta, self._phi], axis=1)
        if n < OUTLINE_MIN_NODES_FOR_HULL:
            return geometry.to_world(self._frame, geometry.local_xyz(self._phi, self._theta))

        half_width = OUTLINE_NEIGHBORHOOD_SPACING_FACTOR * self._typical_spacing()
        counts = np.array(
            [self._rtree.count_in_box(local_xy[i] - half_width, local_xy[i] + half_width) for i in range(n)]
        )
        boundary_mask = counts <= OUTLINE_BOUNDARY_DENSITY_FRACTION * counts.max()
        boundary_xy = local_xy[boundary_mask]
        if len(boundary_xy) < 3:
            boundary_xy = local_xy  # not enough boundary nodes for a hull -- fall back to every node

        try:
            hull = ConvexHull(boundary_xy)
        except QhullError:
            # A degenerate point cloud (collinear, or too few distinct points) has no real
            # 2D hull -- every node is its own outline, same fallback as the too-few-nodes
            # case above.
            return geometry.to_world(self._frame, geometry.local_xyz(self._phi, self._theta))
        hull_xy = boundary_xy[hull.vertices]
        return geometry.to_world(self._frame, geometry.local_xyz(hull_xy[:, 1], hull_xy[:, 0]))

    def update_to_lat_long_grid(self, grid: LatLongGrid) -> None:
        _update_to_lat_long_grid(self, grid)

    def update_deltas_from_lat_long_grid(self, grid: LatLongGrid) -> None:
        _update_deltas_from_lat_long_grid(self, grid)

    def _merge_nodes_with(
        self,
        other: "Plate",
        spacing_rad: float,
        coverage_radius_rad: float,
        other_points_xyz: np.ndarray,
    ) -> None:
        """No fixed row structure to preserve (unlike `PlateWithLines`), so merging is just
        the union of the two plates' own node clouds -- `other`'s nodes converted into this
        plate's own local (theta, phi) frame, no lattice resample and no elevation-nearest-
        neighbor lookup needed. `spacing_rad`/`coverage_radius_rad`/`other_points_xyz` only
        matter to a representation that has to rebuild a fresh lattice (see
        `PlateWithLines._merge_nodes_with`); a plain union can't overlap a third plate's
        territory any more than `self`'s and `other`'s own pre-merge node clouds already
        didn't, so there's nothing here for them to do."""
        other_world_pts, other_elevation = other.all_points_and_elevation()
        other_local = geometry.to_local(self._frame, other_world_pts)
        other_phi, other_theta = geometry.xyz_to_latlon(other_local)
        combined_fields = {
            name: np.concatenate([self._fields[name], other.collect(name)]) for name in ElevationLine.OPTIONAL_FIELDS
        }
        self.set_nodes(
            np.concatenate([self._theta, other_theta]),
            np.concatenate([self._phi, other_phi]),
            np.concatenate([self._elevation, other_elevation]),
            **combined_fields,
        )

    def split(self, new_id: int, cut_normal: np.ndarray, min_nodes: int) -> tuple["Plate", "Plate"] | None:
        """Same partition-by-side-of-`cut_normal` idea as `PlateWithLines.split`, just over
        the flat node cloud instead of per-line -- a boolean mask on every one of this plate's
        own flat arrays (theta/phi/elevation/fields), each half becoming a fresh
        `PlateWithRTree` (rebuilding its own R-tree via `set_nodes`, see `__init__`)."""
        if len(self._theta) == 0:
            return None
        world_pts = geometry.to_world(self._frame, geometry.local_xyz(self._phi, self._theta))
        side = np.sum(world_pts * cut_normal, axis=-1) > 0
        if np.count_nonzero(side) < min_nodes or np.count_nonzero(~side) < min_nodes:
            return None

        def half(plate_id: int, mask: np.ndarray) -> "PlateWithRTree":
            fields = {name: self._fields[name][mask] for name in ElevationLine.OPTIONAL_FIELDS}
            return PlateWithRTree(
                plate_id=plate_id,
                frame=self._frame.copy(),
                crust_type=self._crust_type,
                theta=self._theta[mask],
                phi=self._phi[mask],
                elevation=self._elevation[mask],
                **fields,
            )

        return half(self.plate_id, side), half(new_id, ~side)

    def grow_into(
        self,
        new_points_xyz: np.ndarray,
        new_elevation: np.ndarray,
        coverage_radius_rad: float,
        spacing_rad: float,
    ) -> None:
        """No lattice to resample (unlike `PlateWithLines`), so this just appends
        `new_points_xyz` (converted to this plate's own local frame) directly to the node
        cloud -- same "no resample needed, just union the points" reasoning as
        `_merge_nodes_with`. `coverage_radius_rad`/`spacing_rad` are unused here for the same
        reason `_merge_nodes_with` doesn't need them."""
        new_local = geometry.to_local(self._frame, new_points_xyz)
        new_phi, new_theta = geometry.xyz_to_latlon(new_local)
        new_fields = {
            name: self._default_field(name, new_theta) for name in ElevationLine.OPTIONAL_FIELDS
        }
        self.set_nodes(
            np.concatenate([self._theta, new_theta]),
            np.concatenate([self._phi, new_phi]),
            np.concatenate([self._elevation, new_elevation]),
            **{name: np.concatenate([self._fields[name], values]) for name, values in new_fields.items()},
        )

    def shift(self, world: "World", years: float) -> float:
        # TODO(PlateWithRTree.step): identical to PlateWithLines.shift -- refit omega from
        # mantle.flow_at/fit_euler_pole (damped via mantle.VELOCITY_DAMPING, clamped via
        # mantle.clamp_rate), rotate via geometry.rotation_matrix_from_omega + self.rotate,
        # and measure D as geometry.angular_distance(old_points, new_points).max() over
        # all_points_and_elevation() before/after. Nothing here depends on the flat-cloud
        # representation, so this could plausibly be pulled up into Plate itself and shared
        # by both subclasses rather than duplicated -- left as its own override for now so
        # PlateWithRTree's whole step lands as one self-contained TODO.
        raise NotImplementedError("PlateWithRTree.shift is not yet implemented")

    def deform(self, world: "World", other_plates: list["Plate"], years: float, max_distance: float) -> None:
        # TODO(PlateWithRTree.step): PlateWithLines.deform grows/shrinks at each *line's* two
        # theta-ends; this representation has no row structure, so growth/shrink instead
        # means inserting/removing individual points directly into/out of the flat
        # theta/phi/elevation/fields arrays (see set_nodes), then rebuilding the R-tree.
        # Concretely:
        #   1. Classification is representation-agnostic and reusable as-is: for each of
        #      this plate's own nodes, `contested = any(geometry.point_in_spherical_polygon(
        #      point, n.get_bounding_polygon()) for n in neighbours)`, prefiltered by a
        #      cheap cKDTree nearest-neighbour-distance query exactly like PlateWithLines
        #      does (dist_all/near_mask).
        #   2. Elevation deltas (mountain/trench/transform/rift-relax) apply the same six
        #      intensity curves over the whole flat node array at once -- no per-line loop
        #      needed, just boolean-mask indexing into self._elevation directly.
        #   3. Grow/shrink has no "line end" concept -- the R-tree analogue of "the boundary
        #      nodes of this plate" is exactly what outline_world()'s own boundary-node
        #      detection already identifies (the low-density nodes it hands to ConvexHull).
        #      Shrinking removes contested boundary nodes (and, if the contested region goes
        #      deeper, whichever of *their* own nearest neighbours are also contested,
        #      propagating inward via the R-tree rather than via a fixed row) capped by
        #      max_distance/spacing_rad. Growing inserts new points spaced spacing_rad
        #      outward from an uncontested, far-from-any-neighbour boundary node, along the
        #      local outward normal (e.g. away from this plate's own centroid through that
        #      node) -- there's no single well-defined "theta direction" to extend along the
        #      way a line has, so the outward direction has to be estimated locally instead
        #      (e.g. from the boundary node's own two hull neighbours).
        #   4. The "claim adjacent territory" sub-step becomes: sweep a local lattice (same
        #      iter_local_lattice usage PlateWithLines uses) restricted to a ring just
        #      outside this plate's own R-tree box/nearest-neighbour reach, same
        #      not-contested filter, then Plate.grow_into (already representation-generic --
        #      PlateWithRTree.grow_into above just appends points, no resample needed).
        #   5. Regularizing has no analogue here -- there's no fixed-row spacing to drift out
        #      of alignment, so this step is likely just a no-op for this representation.
        raise NotImplementedError("PlateWithRTree.deform is not yet implemented")


ELLIPSE_OUTLINE_POINTS = 72


@dataclass
class BoundingEllipse:
    center_xyz: np.ndarray  # unit vector, true/un-rotated world frame
    diameter_a_km: float  # major
    diameter_b_km: float  # minor
    outline_xyz: np.ndarray  # (ELLIPSE_OUTLINE_POINTS, 3) unit vectors, true world frame


def plate_bounding_ellipse(points_world_xyz: np.ndarray) -> BoundingEllipse | None:
    """The minimum-area ellipse enclosing `points_world_xyz` (real diameters in km, "rotated
    to fit as closely as possible" per the map-view feature this backs) -- fit in a local
    azimuthal-equidistant projection centered on the point cloud's own `bounding_sphere`
    centroid (exact true-km radial distance from that center; not exact between two
    arbitrary points -- see geometry.azimuthal_equidistant_forward -- but fitting is always
    done relative to that one shared center, so this doesn't matter here).

    Fit against *every* node point, not just `outline_world()`: the minimum enclosing
    ellipse of a full point set is identical to that of just its convex hull, and
    `outline_world()`'s own docstring admits it isn't a guaranteed hull for a concave plate
    ("exact for convex-ish plates, a reasonable envelope otherwise") -- using it risks
    silently missing an interior extremal point. Cost is negligible either way (Khachiyan is
    O(N) per iteration, no O(N^2) anywhere, and a plate's node count is a few thousand at
    most).

    `None` for an empty plate (`node_count() == 0`, e.g. one fully consumed by subduction but
    not yet pruned)."""
    if len(points_world_xyz) == 0:
        return None
    centroid, _ = geometry.bounding_sphere(points_world_xyz)
    east, north = geometry.local_tangent_basis(centroid)
    xy_km = geometry.azimuthal_equidistant_forward(centroid, east, north, points_world_xyz) * PLANET_RADIUS_KM

    fit = ellipse.min_enclosing_ellipse(xy_km)

    t = np.linspace(0.0, 2.0 * np.pi, ELLIPSE_OUTLINE_POINTS, endpoint=False)
    local = np.stack([fit.semi_major * np.cos(t), fit.semi_minor * np.sin(t)], axis=-1)
    cos_a, sin_a = np.cos(fit.angle_rad), np.sin(fit.angle_rad)
    rotate = np.array([[cos_a, -sin_a], [sin_a, cos_a]])
    boundary_km = fit.center + local @ rotate.T
    outline_xyz = geometry.azimuthal_equidistant_inverse(centroid, east, north, boundary_km / PLANET_RADIUS_KM)
    center_xyz = geometry.azimuthal_equidistant_inverse(centroid, east, north, (fit.center / PLANET_RADIUS_KM)[None, :])[0]

    return BoundingEllipse(
        center_xyz=center_xyz,
        diameter_a_km=2.0 * fit.semi_major,
        diameter_b_km=2.0 * fit.semi_minor,
        outline_xyz=outline_xyz,
    )


# Below this many query points, spinning up cKDTree.query's workers=-1 thread pool costs
# more than it saves -- benchmarked against an ~80k-point tree: workers=-1 was ~15x *slower*
# than workers=1 (the default, serial) at 5 query points, still slightly slower at 500, and
# only became a clear win (~2x faster) at 5000+. boundary.py's/hydrology.py's/volcanism.py's
# own per-plate cKDTree queries range from a handful of points (a small or freshly-spawned
# plate) to tens of thousands (an established one at this simulation's default node_density),
# so query_workers below decides per call rather than hardcoding one choice for every plate.
PARALLEL_QUERY_MIN_POINTS = 2000


def query_workers(n: int) -> int:
    """-1 (parallel across every core) if `n` query points is large enough for that to pay
    for itself, else 1 (serial) -- see PARALLEL_QUERY_MIN_POINTS's own comment for the
    benchmark behind the cutoff."""
    return -1 if n >= PARALLEL_QUERY_MIN_POINTS else 1


def gather_node_positions(plate_list: list[Plate]) -> tuple[np.ndarray, list[Plate]]:
    """Every elevation-node's current world position, concatenated, alongside the ordered
    list of contributing plates (every plate in `plate_list` with `node_count() > 0`, in the
    order their nodes appear in the returned array) -- the position-only half of the
    near-identical per-step `_gather_nodes` helpers in erosion.py/hydrology.py/bathymetry.py,
    and climate.py's own `_sample_elevation_and_crust`. Factored out here so a single
    step_world call can compute every node's current world position once and pass the same
    (points, plates_in_order) into all of them, rather than each independently re-deriving
    identical world positions from plate-local data that hasn't moved since the last rotation
    (see docs/architecture.md's World.climate_cache/hydrology_cache notes for the same
    "compute once this step, reuse" precedent).

    `plates_in_order` -- not, as an earlier version of this function returned, (plate,
    line_index, start, end) references into `PlateWithLines`' own `.lines` -- is what makes
    this representation-agnostic: any bulk per-field gather (`collect_all_elevation` and
    friends, below) or per-plate write-back loop (`Plate.map_world_points_on_plate`) already
    visits nodes in this same plate-major order, so a caller never needs to reach into any one
    representation's own storage just to stay aligned with `points`. Each caller still gathers
    its own elevation/other per-node fields fresh (via those bulk collectors) -- only the
    position/plate-order gather itself is shared, since some fields (elevation in particular)
    do change mid-step between callers."""
    plates_in_order = [p for p in plate_list if p.node_count() > 0]
    if not plates_in_order:
        return np.zeros((0, 3)), []
    points = np.concatenate([p.all_points_and_elevation()[0] for p in plates_in_order], axis=0)
    return points, plates_in_order


def collect_all_points(plate_list: list[Plate]) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    """Every plate's current elevation-node positions, elevations, and owning plate_id,
    concatenated -- shared by the render grid's nearest-node resample
    (render_image._render_grid_arrays) and nearest_plate_id's click hit-test below. Takes a
    plain plate list (not World) to avoid a plates.py -> world.py import cycle."""
    points_list, elevation_list, owner_list = [], [], []
    for plate in plate_list:
        pts, elev = plate.all_points_and_elevation()
        if len(pts) == 0:
            continue
        points_list.append(pts)
        elevation_list.append(elev)
        owner_list.append(np.full(len(pts), plate.plate_id))
    if not points_list:
        return None
    return (
        np.concatenate(points_list, axis=0),
        np.concatenate(elevation_list, axis=0),
        np.concatenate(owner_list, axis=0),
    )


def _collect_all(plate_list: list[Plate], field_name: str) -> np.ndarray:
    """Every plate's current `field_name` (elevation or an ElevationLine OPTIONAL_FIELDS
    name), concatenated in the exact same per-plate/per-node order collect_all_points uses --
    so results from two different `_collect_all` calls can still be indexed together with the
    same nearest-neighbor result (see render_image._render_grid_arrays). Delegates to each
    plate's own `collect` -- representation-independent, works against the abstract `Plate`
    interface, not just PlateWithLines."""
    chunks = [p.collect(field_name) for p in plate_list if p.node_count() > 0]
    if not chunks:
        return np.zeros(0, dtype=bool) if field_name == "is_volcano" else np.zeros(0)
    return np.concatenate(chunks, axis=0)


def collect_all_lake_depth(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "lake_depth")


def collect_all_glacier_depth(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "glacier_depth")


def collect_all_silt_depth(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "silt_depth")


def collect_all_channel_depth(plate_list: list[Plate]) -> np.ndarray:
    """Used by climate.py to size a river's own evaporative surface for its moisture-
    recycling humidity source (see that module)."""
    return _collect_all(plate_list, "channel_depth")


def collect_all_channel_width(plate_list: list[Plate]) -> np.ndarray:
    """Used by render_image.py to draw a wide river thicker than a narrow one."""
    return _collect_all(plate_list, "channel_width")


def collect_all_is_volcano(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "is_volcano")


def collect_all_elevation(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "elevation")


def collect_all_soil_depth(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "soil_depth")


def collect_all_soil_mineral_content(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "soil_mineral_content")


def collect_all_soil_organic_content(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "soil_organic_content")


def collect_all_coal_deposit(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "coal_deposit_m")


def collect_all_oil_gas_deposit(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "oil_gas_deposit_m")


def collect_all_mineral_deposit(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "mineral_deposit_m")


def nearest_plate_id(plate_list: list[Plate], query_xyz: np.ndarray) -> int | None:
    """Which plate owns the node nearest `query_xyz` -- the Plate Inspector's click
    hit-test. `None` if every plate is empty (shouldn't happen via the API, but a
    freshly-constructed empty World has no plates at all)."""
    collected = collect_all_points(plate_list)
    if collected is None:
        return None
    points, _, owner = collected
    _, idx = cKDTree(points).query(query_xyz)
    return int(owner[idx])


def base_elevation(crust_type: str) -> float:
    return BASE_CONTINENTAL_M if crust_type == "continental" else BASE_OCEANIC_M


def noise_amplitude(crust_type: str) -> float:
    return CONTINENTAL_NOISE_AMPLITUDE_M if crust_type == "continental" else OCEANIC_NOISE_AMPLITUDE_M


def _build_lines_for_plate(
    plate_index: int,
    frame: np.ndarray,
    crust_type: str,
    owner_tree: cKDTree,
    noise: SphereNoise,
    land_threshold: float | None = None,
    spacing_rad: float = TARGET_LINE_SPACING_RAD,
) -> list[ElevationLine]:
    """Keep only lattice nodes whose nearest seed is this plate's own seed (i.e. nodes
    actually inside this plate's spherical Voronoi cell), and assign each a base elevation
    plus noise texture. `land_threshold` (continental crust only, see
    _land_noise_threshold) overrides the usual fixed BASE_CONTINENTAL_M floor with one
    derived from the requested land_fraction, so elevation = amp * (noise - threshold) is
    positive for exactly the fraction of continental crust needed to hit that target."""
    amp = noise_amplitude(crust_type)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        _, nearest_idx = owner_tree.query(world_pts)
        return nearest_idx == plate_index

    if crust_type == "continental" and land_threshold is not None:

        def elevation_at(world_pts: np.ndarray) -> np.ndarray:
            return amp * (noise.sample(world_pts) - land_threshold)
    else:
        base = base_elevation(crust_type)

        def elevation_at(world_pts: np.ndarray) -> np.ndarray:
            return base + amp * noise.sample(world_pts)

    return build_lines_from_lattice(frame, is_owned, elevation_at, spacing_rad=spacing_rad)


def _land_noise_threshold(
    owner_tree: cKDTree, crust_types: list[str], noise: SphereNoise, land_fraction: float
) -> float | None:
    """Translate a requested whole-sphere land_fraction into a concrete noise threshold for
    continental crust's elevation formula (see _build_lines_for_plate).

    A one-off whole-sphere sweep (independent of any plate's own lattice, at the coarser
    LAND_FRACTION_SAMPLE_SPACING_RAD -- this only needs to be a statistically representative
    sample) measures both which crust_type each sample point would land in (nearest-seed,
    the same rule that decides real plate territory) and that point's noise value. The
    measured continental *area* fraction -- not just the continental *plate count* fraction
    passed in as continental_fraction, which can differ meaningfully since Voronoi cells
    from random seed points aren't equal-area -- sets how much of that continental area
    needs to end up above sea level to hit the requested whole-sphere land_fraction: e.g. if
    continental crust only covers 40% of the sphere but 29% land was requested, ~72% of
    continental crust needs to be land. Returns None if there's no continental crust at all
    to place land on."""
    sample_pts = np.concatenate(
        [
            world_pts
            for _, _, world_pts in iter_local_lattice(np.eye(3), spacing_rad=LAND_FRACTION_SAMPLE_SPACING_RAD)
        ],
        axis=0,
    )
    _, nearest_idx = owner_tree.query(sample_pts)
    is_continental = np.array([crust_types[i] == "continental" for i in nearest_idx])
    continental_area_fraction = float(np.mean(is_continental))
    if continental_area_fraction <= 0.0:
        return None

    target_sub_fraction = min(land_fraction / continental_area_fraction, 1.0)
    continental_noise = noise.sample(sample_pts[is_continental])
    return float(np.quantile(continental_noise, 1.0 - target_sub_fraction))


def generate_plates(
    seed: int,
    num_plates: int | None = None,
    continental_fraction: float | None = None,
    land_fraction: float | None = None,
    node_density: float = DEFAULT_NODE_DENSITY,
) -> list[Plate]:
    """Tile the whole sphere into plates. `num_plates` is optional -- when omitted, a
    plausible Earth-like count is drawn from the seed's own RNG stream (so it's still fully
    determined by `seed`, just not something the caller has to pick). `continental_fraction`
    is also optional -- when given (the UI's "continental plates" slider, 0 to 1), that
    fraction of plates (rounded, `num_plates` bumped up if needed so there's still room for
    at least MIN_OCEANIC_PLATES of real ocean floor) are made continental, instead of the
    usual independent CONTINENTAL_FRACTION coin flip per plate. `land_fraction` (the UI's
    "initial land" slider, also 0 to 1) similarly overrides how much of the *whole sphere*
    -- not just of continental crust -- starts above sea level; see
    _land_noise_threshold for how that target is actually hit. `node_density` (the UI's
    "point density" choice, see NODE_DENSITY_CHOICES) scales how many elevation-line nodes
    each plate starts with -- see line_spacing_rad.

    Every plate's territory comes from the same nearest-seed test (`owner_tree.query`
    below): each lattice node is claimed by exactly one plate, so the tiling has no gaps
    and no overlaps by construction -- there's no separate polygon-boundary step that could
    fall out of sync with it (see Plate.outline_world for the live, rendering-only outline
    derived from this same data after the world has evolved)."""
    rng = np.random.default_rng(seed)
    if num_plates is None:
        num_plates = int(rng.integers(MIN_AUTO_PLATES, MAX_AUTO_PLATES + 1))

    num_continents: int | None = None
    if continental_fraction is not None:
        continental_fraction = max(0.0, min(continental_fraction, 1.0))
        num_continents = round(continental_fraction * num_plates)
        num_plates = max(num_plates, num_continents + MIN_OCEANIC_PLATES)

    seed_xyz = rng.normal(size=(num_plates, 3))
    seed_xyz /= np.linalg.norm(seed_xyz, axis=-1, keepdims=True)

    if num_continents is None:
        crust_types = [
            "continental" if rng.random() < CONTINENTAL_FRACTION else "oceanic" for _ in range(num_plates)
        ]
    else:
        continental_indices = set(rng.choice(num_plates, size=num_continents, replace=False).tolist())
        crust_types = ["continental" if i in continental_indices else "oceanic" for i in range(num_plates)]

    owner_tree = cKDTree(seed_xyz)
    noise = SphereNoise(rng, octaves=4, base_freq=2.5)

    land_threshold = None
    if land_fraction is not None:
        land_fraction = max(0.0, min(land_fraction, 1.0))
        land_threshold = _land_noise_threshold(owner_tree, crust_types, noise, land_fraction)

    spacing_rad = line_spacing_rad(node_density)
    plates: list[Plate] = []
    for i in range(num_plates):
        frame = geometry.plate_frame_from_seed(seed_xyz[i])
        lines = _build_lines_for_plate(i, frame, crust_types[i], owner_tree, noise, land_threshold, spacing_rad=spacing_rad)
        plates.append(PlateWithLines(plate_id=i, frame=frame, crust_type=crust_types[i], lines=lines))
    return plates
