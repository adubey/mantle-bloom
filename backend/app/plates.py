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
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull, QhullError, cKDTree

from . import ellipse, geometry
from . import elevation_lines
from .elevation_lines import (
    DEFAULT_NODE_DENSITY,
    NODE_DENSITY_CHOICES,
    PLANET_RADIUS_KM,
    TARGET_LINE_SPACING_KM,
    TARGET_LINE_SPACING_RAD,
    ElevationLine,
    ElevationPoint,
    install_point_field_accessors,
    iter_local_lattice,
    line_spacing_rad,
    split_into_contiguous_runs,
)

if TYPE_CHECKING:
    from . import terrain_noise

CONTINENTAL_FRACTION = 0.4
BASE_CONTINENTAL_M = 200.0
BASE_OCEANIC_M = -3800.0
# Widened from the original 1200/500 so a freshly generated world already shows real relief --
# rolling hills and real ocean-basin variation, not a flat plain/seafloor waiting for tectonics
# to draw the first contours. Still well inside MIN/MAX_ELEVATION_M (-11000/9000) either way,
# and continental crust's own land/sea split near sea level is unaffected (BASE_CONTINENTAL_M,
# or the land_fraction-derived threshold when one is given, is untouched -- only how far the
# noise texture swings around whichever baseline is already used).
CONTINENTAL_NOISE_AMPLITUDE_M = 2000.0
OCEANIC_NOISE_AMPLITUDE_M = 900.0

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
    recompute the same unchanged outlines from scratch every time (see get_neighbours' own
    callers, e.g. torque.py's per-plate neighbour torque and lithosphere_plate.py's own
    boundary-reach query). Same bounding-sphere
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


# --- Deformation constants (LithospherePlate.deform) ---
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
DIVERGENT_RELAX_RATE_PER_MYR = 0.15  # was 0.5 -- see the comment below

# The "divergent" classification is purely geometric (uncontested and within reach of a
# neighbour, see deform()'s own comment) -- it can't by itself tell a genuinely active,
# still-subsiding rift apart from a long-settled passive margin that simply happens to sit
# near a neighbour (most of a real continent's own coastline, in reality: the Atlantic
# seaboard hasn't been actively rifting for tens of millions of years, yet still reads as
# "uncontested and close to an oceanic neighbour" under this same geometric test every
# turn). Confirmed directly as a real, previously-unnoticed drain -- and confirmed as two
# separate contributors, needing both constants below to actually fix: at the old rate
# (0.5, ~78% of the remaining gap closed per 3 Myr step), already-elevated, long-stable
# coastline got yanked most of the way to DIVERGENT_RIFT_TARGET_M/DIVERGENT_RIDGE_TARGET_M
# within a step or two of ever qualifying as divergent -- so fast that DIVERGENT_YOUNG_AGE_
# MYR's own age gate (below) mostly closes the barn door after the horse has left: a plate's
# own ongoing rotation continuously sweeps *fresh* coastline into the divergent band, so most
# of the loss was already-done first-time hits, not repeat hits on the same settled land the
# age gate alone can prevent. Slowing the rate itself (0.5 -> 0.15) shrinks *every* hit,
# first time included -- the tradeoff being that genuinely active rifting/ridge spreading
# also settles more gradually now, not just the stale-coastline case, so a fresh rift no
# longer reaches its target in a step or two the way it used to.
#
# DIVERGENT_YOUNG_AGE_MYR still matters on top of the slower rate: deform() only relaxes a
# node while ElevationLine.divergent_age_myr (Myr spent *continuously* divergent, reset to 0
# the instant it isn't) stays under this threshold, so a margin that's stayed divergent this
# long is treated as mature and left alone from then on, no matter how much longer it keeps
# testing as geometrically divergent -- still needed even at the slower rate, since given
# enough consecutive divergent steps a node would otherwise keep creeping toward the target
# forever rather than ever actually settling.
DIVERGENT_YOUNG_AGE_MYR = 10.0

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

# Reverse faults: real shortening in a collision belt isn't smooth vertical uplift spread
# evenly across the whole zone -- fold-thrust belts partition it into discrete thrust sheets
# (fast-rising ridges) separated by footwall synclines/intermontane basins that keep rising far
# more slowly, a real, well-documented process (the north-south rift valleys cutting straight
# across the Tibetan Plateau's own overall convergent uplift; Basin-and-Range-style extension
# nested inside the Anatolian collision zone). Modeled as a smooth, deterministic noise field
# sampled in the plate's own *local* frame (geometry.local_xyz(line.phi, line.theta), not world
# xyz -- see LithospherePlate.deform's own use) so a given downthrown block stays attached to the
# same crust as the plate rotates, the same "attached to the crust, not the world" property
# every other persistent field in this codebase already has (see docs/simulation-model.md's
# "Why not a grid"). Seeded from (world.seed, plate_id) only, not elapsed_years, so the fault
# pattern is a fixed geological feature of this plate rather than reshuffling every step.
REVERSE_FAULT_SEED_TAG = 9001  # arbitrary, distinguishes this RNG stream from any other keyed by (world.seed, plate_id, ...)
REVERSE_FAULT_NOISE_FREQ = 9.0
REVERSE_FAULT_VALLEY_THRESHOLD = -0.15  # noise below this reads as a downthrown fault block
REVERSE_FAULT_VALLEY_UPLIFT_FACTOR = 0.15  # a valley block still rises, just far slower than a thrust ridge

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

# How many target spacings of margin `_claim_adjacent_territory` keeps between a plate's
# poleward-most row and its own local pole (+-pi/2). A row is a circle of local latitude, so
# its theta step (spacing_rad / cos(phi)) blows up as cos(phi) -> 0: without a margin a plate
# that grows to encircle its own pole ends up with a handful of degenerate near-pole rings,
# and -- since nothing here treats theta as periodic -- ordinary end-growth just keeps winding
# those rings past a full revolution, covering the same ground many times over (the
# concentric-circle / moire "holes" artifacts in the Plate Inspector, plus a real unbounded
# contribution to plate overlap and node count on long runs). This margin keeps near-pole
# rows at a sane circumference; `_grow_or_shrink_line_for_deform` separately refuses to extend
# any row past 2*pi (`_ROW_FULL_REVOLUTION_SLACK`). Generation's own lattice sweep
# (`iter_local_lattice`) still fills to the pole -- a plate that legitimately owns the pole at
# generation keeps its small Voronoi-clipped rings; only *growth* toward the pole is capped.
POLE_CAP_MARGIN_MULT = 4.0
# A row whose theta span is within this many target-spacing steps of a full 2*pi revolution
# is treated as a closed ring and never grown further (nor, in regularize_line, resampled to
# more than one revolution).
_ROW_FULL_REVOLUTION_SLACK = 1.0

# Minimum length (nodes) of an interior `shrinkable` run for `_grow_or_shrink_line_for_deform`
# to carve it out and split the row into separate contiguous `ElevationLine`s -- see that
# method. Shorter transient interior contests are left for when they reach an end: a 1-2 node
# gap wouldn't clear `elevation_lines.CONTIGUOUS_RUN_GAP_MULT` to survive as a real split, and
# would just be refilled by the next regularize pass. Set to CONTIGUOUS_RUN_GAP_MULT so the
# post-removal gap is always wide enough for `split_into_contiguous_runs` to actually break.
_INTERIOR_SUBDUCTION_MIN_RUN = int(elevation_lines.CONTIGUOUS_RUN_GAP_MULT)

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


def _row_median_step(line: ElevationLine) -> float | None:
    """The typical theta step of `line` before it's masked by a partition -- passed to
    `largest_contiguous_run` as its reference spacing so it can still recognise a
    partition-stranded two-node row (too short to estimate a spacing from what survives). A
    pre-partition row is contiguous and regularized, so its median step is a clean estimate;
    `None` for a row with fewer than two nodes."""
    if len(line) < 2:
        return None
    return float(np.median(np.diff(line.theta)))


def node_components(points_xyz: np.ndarray, connect_radius_rad: float) -> np.ndarray:
    """Label each of `points_xyz` (world unit vectors) with a connected-component id, where
    two nodes are connected if they sit within `connect_radius_rad` of each other. Used by
    `Plate.defragment` to tell a plate that's been physically severed into two landmasses
    (still carried as one `Plate`) from one that's merely shed a few stranded nodes -- see
    that method. Component ids are contiguous from 0 but otherwise arbitrary (not size-
    ordered)."""
    n = len(points_xyz)
    if n == 0:
        return np.zeros(0, dtype=int)
    pairs = cKDTree(points_xyz).query_pairs(connect_radius_rad, output_type="ndarray")
    if len(pairs) == 0:
        return np.arange(n)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return labels


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
        internal_stress: float = 0.0,
    ) -> None:
        self._plate_id = plate_id
        self._frame = frame
        self._crust_type = crust_type
        self._omega = omega if omega is not None else np.zeros(3)
        self._age_steps = age_steps
        # Accumulated breakup pressure -- see merge_split.accumulate_plate_stress (background,
        # size-driven accumulation plus an overlap-driven top-up biased toward whichever plate
        # in an overlapping pair is the larger one) and maybe_split_plate (folds this into the
        # existing size-based split-gate relaxation). Reset to 0 on a successful split, halved
        # on a failed rift -- see reset_age's own precedent for "a topology event releases
        # accumulated pressure."
        self._internal_stress = internal_stress
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
        # A cKDTree over this plate's full node cloud (all_points_and_elevation()[0]), built
        # lazily by get_node_kdtree() and invalidated in lockstep with the two caches above --
        # node *positions*, like the outline, only change via rotate() or a node-set mutation,
        # never an elevation-only edit. torque.gather_boundary_force_inputs queries every
        # neighbour's tree once per plate per pass; sharing one cached tree per plate turns
        # ~24 fresh per-call tree builds a step into one build per plate.
        self._node_kdtree_cache: cKDTree | None = None

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

    @property
    def internal_stress(self) -> float:
        """Accumulated breakup pressure -- see merge_split.accumulate_plate_stress /
        maybe_split_plate. Dimensionless (a relaxation-time-scaled accumulator, not a literal
        Pa stress -- see that module's own comment); 0.0 for a quiet, small, non-overlapped
        plate, and for every plate on a save written before this field existed (see
        __getattr__)."""
        return self._internal_stress

    def set_internal_stress(self, value: float) -> None:
        self._internal_stress = value

    def __getattr__(self, name: str):
        """A `Plate` unpickled from a save written before `internal_stress` existed has no
        `_internal_stress` in its restored `__dict__` (pickle bypasses `__init__` entirely) --
        default it to 0.0, the same "quiet, unstressed plate" reading a fresh Plate.__init__
        gives, rather than raising. Mirrors ElevationLine.__getattr__'s own precedent for the
        same class of backward-compatibility gap."""
        if name == "_internal_stress":
            object.__setattr__(self, name, 0.0)
            return 0.0
        raise AttributeError(name)

    @abc.abstractmethod
    def node_count(self) -> int: ...

    def has_negligible_territory(self) -> bool:
        """True once this plate has been whittled down to no real remaining territory --
        the "no land left" half of merge_split.remove_defunct_plates/apply_topology_changes'
        own pruning (the other half, `node_count() == 0`, is checked separately by both
        callers). Default: fewer nodes than can form a real 2D hull (see
        `OUTLINE_MIN_NODES_FOR_HULL`) -- a representation-agnostic floor any subclass without
        its own more specific notion can fall back to. `PlateWithLines` overrides this with
        its own exact, pre-existing "at most one line left" definition, so that behavior is
        unchanged for it; `PlateWithMesh`/`PlateWithRTree` use this default."""
        return self.node_count() < OUTLINE_MIN_NODES_FOR_HULL

    def defragment(
        self, next_id: int, connect_radius_rad: float, min_fragment_nodes: int
    ) -> tuple[list["Plate"], int] | None:
        """Reconcile "one `Plate` object" with "one contiguous patch of crust."

        Ordinary per-step `deform()` only ever grows/shrinks a line's *ends*, and its shrink
        rule deliberately never deletes a line's last node (see
        `_grow_or_shrink_line_for_deform`) -- so subduction/transform can carve a plate's
        node cloud into two disconnected landmasses, or strand a comb of one-node rows far
        from the plate body, and nothing notices: `maybe_split_plate` only cuts on mantle-
        flow *disagreement*, not geometry, and two co-moving lobes never trip it.

        This finds those cases directly. Connected components of this plate's nodes at
        `connect_radius_rad` (see `node_components`); each component with at least
        `min_fragment_nodes` nodes becomes its own plate (the largest keeps this plate's own
        id/frame/omega/age -- see `_plates_from_node_masks`), everything smaller is dropped
        as stranded crust.

        Returns `None` -- nothing to do -- when the plate is already a single contiguous
        patch (the overwhelmingly common case), when it has exactly one component big
        enough to anchor a plate and no stranded nodes to shed, or when it's debris with no
        component large enough to anchor a plate at all (left for `has_negligible_territory`
        / `remove_defunct_plates` to prune). Otherwise returns
        `(replacement_plates, n_new_ids_consumed)`, where `replacement_plates[0]` reuses
        this plate's own id and `next_id, next_id + 1, ...` are consumed for the rest, in
        descending component-size order. `next_id` is `World.next_plate_id`."""
        points, _ = self.all_points_and_elevation()
        if len(points) < 2:
            return None

        labels = node_components(points, connect_radius_rad)
        component_ids, counts = np.unique(labels, return_counts=True)
        if len(component_ids) == 1:
            return None

        # Largest component first, so it's the one that keeps this plate's identity.
        order = np.argsort(counts)[::-1]
        kept = [int(component_ids[i]) for i in order if counts[i] >= min_fragment_nodes]
        dropped_nodes = len(points) - int(counts[np.isin(component_ids, kept)].sum())
        # No component big enough to anchor a plate -- the whole thing is debris. Leave it
        # for merge_split.remove_defunct_plates / has_negligible_territory to prune; defrag
        # never deletes a whole plate itself (that path is fragile against small synthetic
        # plates and adds nothing the negligible-territory check doesn't already do).
        if not kept:
            return None
        if len(kept) == 1 and dropped_nodes == 0:
            return None

        n_new_ids = len(kept) - 1
        masks = [labels == cid for cid in kept]
        ids = [self.plate_id, *range(next_id, next_id + n_new_ids)]
        return self._plates_from_node_masks(masks, ids), n_new_ids

    def _plates_from_node_masks(self, masks: list[np.ndarray], ids: list[int]) -> list["Plate"]:
        """Build one plate per mask in `masks` (each a boolean array over this plate's nodes
        in `all_points_and_elevation` order), assigning `ids[k]` to `masks[k]`'s plate --
        `ids[0]` is always this plate's own id, so the first mask should be the one that
        keeps this plate's identity. Representation-specific; only `PlateWithLines`
        implements it (the one representation `defragment` is wired up for)."""
        raise NotImplementedError

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

    def get_node_kdtree(self) -> cKDTree | None:
        """A `cKDTree` over `all_points_and_elevation()[0]`, cached and invalidated the same
        way `get_bounding_polygon_tree()` is -- `None` if this plate currently has no nodes.
        Shared across `torque.gather_boundary_force_inputs`' per-neighbour nearest-node
        queries (one plate is a neighbour of several others, and is queried in both the shift
        and deform pass) so its node cloud is treed once per step, not once per query."""
        if self._node_kdtree_cache is None:
            points = self.all_points_and_elevation()[0]
            if len(points) == 0:
                return None
            self._node_kdtree_cache = cKDTree(points, balanced_tree=False, compact_nodes=False)
        return self._node_kdtree_cache

    def _invalidate_bounding_polygon(self) -> None:
        self._bounding_polygon_cache = None
        self._bounding_polygon_tree_cache = None
        self._node_kdtree_cache = None

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
        "close enough" (e.g. a boundary-effect or force-reach radius) can pass their own."""
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


@dataclass
class _RowLookup:
    """`PlateWithLines.contains_batch`'s cached fast-path data -- see that method's own
    docstring for the algorithm. `phis` is every *distinct* line phi, sorted ascending (a
    row split into arcs by interior subduction has several lines at one phi -- see
    `_grow_or_shrink_line_for_deform`). `low_thetas[i]`/`high_thetas[i]` are that row's
    overall theta envelope (min arc low / max arc high), index-aligned with `phis`;
    `interval_lo`/`interval_hi` are `(n_rows, max_arcs)` with each row's actual arc
    intervals, absent slots padded (+inf / -inf) so they never match. The rest are derived
    once here rather than per query:

    `margin_rad` -- see `_row_lookup_bulge_margin_rad`'s own docstring: a query point whose
    nearest-row interval test says "outside" isn't necessarily outside -- this bounds how
    far the *true* boundary can lie beyond the idealized per-row interval.
    `phi_min_pad`/`phi_max_pad` -- `phis[0]`/`phis[-1]`, padded outward by `margin_rad`.
    `padded_low`/`padded_high` -- `low_thetas`/`high_thetas`, each widened by `margin_rad`
    *and* by whichever of its own immediate row-neighbours (index - 1, index + 1) reaches
    further -- covers a query point landing just past a shelf-step boundary, whose relevant
    bulge belongs jointly to the two rows either side of that step, not to its own nearest
    row alone. (An interior hole's rim within `margin_rad` falls back to the winding test
    against the keyholed `get_bounding_polygon()`, same as any other near-boundary point.)"""

    phis: np.ndarray
    low_thetas: np.ndarray
    high_thetas: np.ndarray
    interval_lo: np.ndarray
    interval_hi: np.ndarray
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


def _row_intervals(lines: list[ElevationLine]) -> list[tuple[float, list[tuple[float, float]]]]:
    """`(phi, [(theta_lo, theta_hi), ...])` per *distinct* phi across `lines`, phi ascending
    and each row's intervals sorted, non-overlapping. One `ElevationLine` is one interval;
    a row split by interior subduction (see `_grow_or_shrink_line_for_deform`) contributes
    two or more. Touching/overlapping intervals at one phi (not expected -- arcs are carved
    with a real gap between them) are merged so downstream gap detection stays clean."""
    by_phi: dict[float, list[tuple[float, float]]] = {}
    for line in lines:
        if len(line) == 0:
            continue
        by_phi.setdefault(line.phi, []).append((float(line.theta[0]), float(line.theta[-1])))
    rows: list[tuple[float, list[tuple[float, float]]]] = []
    for phi in sorted(by_phi):
        merged: list[list[float]] = []
        for lo, hi in sorted(by_phi[phi]):
            if merged and lo <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], hi)
            else:
                merged.append([lo, hi])
        rows.append((phi, [(lo, hi) for lo, hi in merged]))
    return rows


def _interval_complement(lo: float, hi: float, cover: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """The sub-intervals of `[lo, hi]` left uncovered by any interval in `cover`."""
    segs = [(lo, hi)]
    for c_lo, c_hi in cover:
        nxt: list[tuple[float, float]] = []
        for s_lo, s_hi in segs:
            if c_hi <= s_lo or c_lo >= s_hi:
                nxt.append((s_lo, s_hi))
                continue
            if s_lo < c_lo:
                nxt.append((s_lo, c_lo))
            if c_hi < s_hi:
                nxt.append((c_hi, s_hi))
        segs = nxt
    return segs


def _plate_outline_loops(rows: list[tuple[float, list[tuple[float, float]]]]) -> list[list[tuple[float, float]]]:
    """Every closed boundary loop of a plate's territory, as the exact boundary of the union
    of its rows' theta-intervals (each row's band runs from the midpoint phi with the row
    below to the midpoint with the row above; the two extreme rows' outer edges sit at their
    own phi, matching the old single-interval staircase). One outer loop for a simple plate;
    additional inner loops, wound opposite, for every hole a split row leaves between its
    arcs (interior subduction -- see `_grow_or_shrink_line_for_deform`) or a notch a
    partial override cuts. `_stitch_loops` joins them into the single vertex array
    `get_bounding_polygon()` returns.

    Handles the general case directly (holes, one-sided notches, disjoint pieces) rather
    than the old "outer staircase + separately detected enclosed holes" split, which
    mishandled a hole that stays open where an adjacent row was end-eroded instead of
    interior-split."""
    n = len(rows)
    if n == 0:
        return []
    phis = [phi for phi, _ in rows]
    band_lo = [phis[0], *[(phis[i - 1] + phis[i]) / 2.0 for i in range(1, n)]]
    band_hi = [*[(phis[i] + phis[i + 1]) / 2.0 for i in range(n - 1)], phis[-1]]

    # Axis-aligned boundary segments in (phi, theta): every interval's two vertical edges
    # (full band height), plus the parts of its top/bottom edges not shared with the
    # neighbouring row's coverage.
    adjacency: dict[tuple[float, float], list[tuple[float, float]]] = {}

    def add_segment(a: tuple[float, float], b: tuple[float, float]) -> None:
        if a == b:
            return
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)

    for i, (_phi, ivs) in enumerate(rows):
        above = rows[i + 1][1] if i + 1 < n else []
        below = rows[i - 1][1] if i - 1 >= 0 else []
        lo_b, hi_b = band_lo[i], band_hi[i]
        for lo, hi in ivs:
            add_segment((lo_b, lo), (hi_b, lo))
            add_segment((lo_b, hi), (hi_b, hi))
            for s_lo, s_hi in _interval_complement(lo, hi, above):
                add_segment((hi_b, s_lo), (hi_b, s_hi))
            for s_lo, s_hi in _interval_complement(lo, hi, below):
                add_segment((lo_b, s_lo), (lo_b, s_hi))

    # Walk the segment graph into closed loops, keeping the covered region on the left at
    # every vertex (turn as far counterclockwise as the available edges allow). That orients
    # the outer boundary CCW and every hole CW -- opposite windings, so a keyhole stitch
    # cancels to zero inside the holes.
    used: set[tuple[tuple[float, float], tuple[float, float]]] = set()

    def turn_key(incoming: tuple[float, float], outgoing: tuple[float, float]) -> int:
        # incoming/outgoing are unit cardinal directions (dx, dy). Rank: left(0) < straight(1)
        # < right(2) < back(3), by the cross product and dot of the two headings.
        ix, iy = incoming
        ox, oy = outgoing
        cross = ix * oy - iy * ox
        dot = ix * ox + iy * oy
        if cross > 0:
            return 0
        if cross == 0 and dot > 0:
            return 1
        if cross < 0:
            return 2
        return 3

    def direction(a: tuple[float, float], b: tuple[float, float]) -> tuple[float, float]:
        dx, dy = b[1] - a[1], b[0] - a[0]  # (theta, phi) as (x, y)
        if abs(dx) >= abs(dy):
            return (1.0 if dx > 0 else -1.0, 0.0)
        return (0.0, 1.0 if dy > 0 else -1.0)

    def covered(phi_t: float, theta_t: float) -> bool:
        for i in range(n):
            if band_lo[i] <= phi_t <= band_hi[i]:
                if any(lo <= theta_t <= hi for lo, hi in rows[i][1]):
                    return True
        return False

    def left_is_covered(loop: list[tuple[float, float]]) -> bool:
        # Test a point a hair to the left of the loop's first edge -- keep the loop only if
        # covered territory sits on its left (drops the unbounded face and each hole's own
        # interior face, keeps the outer boundary and every real hole rim).
        (ay, ax), (by, bx) = loop[0], loop[1]
        dx, dy = bx - ax, by - ay
        norm = (dx * dx + dy * dy) ** 0.5
        if norm == 0:
            return False
        eps = 1e-6
        theta_t = (ax + bx) / 2.0 + eps * (-dy / norm)
        phi_t = (ay + by) / 2.0 + eps * (dx / norm)
        return covered(phi_t, theta_t)

    loops: list[list[tuple[float, float]]] = []
    for start, neighbours in adjacency.items():
        for first in neighbours:
            if (start, first) in used:
                continue
            loop = [start]
            prev, cur = start, first
            while True:
                used.add((prev, cur))
                loop.append(cur)
                incoming = direction(prev, cur)
                candidates = [nb for nb in adjacency[cur] if (cur, nb) not in used and nb != prev]
                if not candidates and (cur, prev) not in used and prev != cur:
                    candidates = [prev]  # dead end -- retrace the spur
                if not candidates:
                    break
                nxt = min(candidates, key=lambda nb: turn_key(incoming, direction(cur, nb)))
                prev, cur = cur, nxt
                if cur == start:
                    used.add((prev, cur))
                    break
                if (prev, cur) in used:
                    break
            if len(loop) >= 4 and left_is_covered(loop):
                loops.append(loop)
    return loops


def _stitch_loops(loops: list[list[tuple[float, float]]]) -> list[tuple[float, float]]:
    """Join several boundary loops into one vertex loop via zero-width keyhole seams -- a
    degenerate out-and-back seam edge pair contributes zero winding everywhere off itself,
    so the winding-number test still reads covered territory as inside and holes as outside,
    and every `get_bounding_polygon()` consumer keeps working on one plain array. The seam
    for each loop runs between its nearest vertex and the nearest vertex of the loop stitched
    so far."""
    if not loops:
        return []
    combined = list(loops[0])
    for extra in loops[1:]:
        best = None
        for i, (pi, ti) in enumerate(combined):
            for j, (pj, tj) in enumerate(extra):
                d = (pi - pj) ** 2 + (ti - tj) ** 2
                if best is None or d < best[0]:
                    best = (d, i, j)
        _, i, j = best
        rotated = extra[j:] + extra[:j]
        combined = combined[: i + 1] + rotated + [extra[j], combined[i]] + combined[i + 1 :]
    return combined


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
        internal_stress: float = 0.0,
    ) -> None:
        super().__init__(plate_id, frame, crust_type, omega=omega, age_steps=age_steps, internal_stress=internal_stress)
        self._lines: list[ElevationLine] = list(lines) if lines is not None else []
        # Lazily (re)built by _get_row_lookup() below, invalidated in lockstep with the
        # bounding-polygon cache (same rotate()/set_lines()/replace_line() call sites) --
        # see contains_batch's own docstring for what this backs.
        self._row_lookup_cache: _RowLookup | None = None
        # Every non-empty line's node world-xyz, concatenated in line order -- a pure
        # function of each line's plate-local (phi, theta) and this plate's frame, so it
        # only changes on rotate() or a node-set mutation, never an elevation-only edit.
        # Built in one vectorized pass by _get_world_points() (a single local_xyz + frame
        # rotation over the whole plate, not a small pair of numpy calls per line) and
        # invalidated in lockstep with the bounding-polygon caches. Backs
        # all_points_and_elevation, itself called dozens of times per step for the same
        # unchanged geometry (torque's per-neighbour shift/deform passes, erosion,
        # merge/defrag checks). Read-only for callers, same as get_bounding_polygon().
        self._world_points_cache: np.ndarray | None = None

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

    def has_negligible_territory(self) -> bool:
        """"No real remaining territory" for this representation. Two cases:

        - At most one non-empty line (the original definition -- a sliver along one
          latitude).
        - A *comb of stubs*: many lines but barely more than one node each on average
          (< 2). Ordinary `deform()` shrinks a line only from its ends and never deletes a
          line's last node, so a heavily-subducted oceanic plate decays into 100+ rows of
          one stranded node apiece -- a high line count masking that there's no 2D patch
          left. The original `len(self.lines) <= 1` never caught this; the ratio test is
          scale-free (same at any node_density) and sits far below any legitimate plate
          (whose rows carry tens of nodes -- `maybe_split_plate`/defrag both floor a real
          plate well above this). The base-class node-count floor still covers
          `PlateWithMesh`/`PlateWithRTree`."""
        nonempty = [line for line in self._lines if len(line) > 0]
        if len(nonempty) <= 1:
            return True
        return sum(len(line) for line in nonempty) < 2.0 * len(nonempty)

    def _invalidate_bounding_polygon(self) -> None:
        super()._invalidate_bounding_polygon()
        self._row_lookup_cache = None
        self._world_points_cache = None

    def outline_world(self) -> np.ndarray:
        """Derived directly from each line's current endpoints -- the actual edge deform()
        maintains -- rather than a separately-tracked polygon that could drift out of sync
        with the real data. The exact boundary of the union of every row's theta-interval(s),
        each row's band spanning the midpoint phi to the row below and above (a straight
        diagonal between two rows with very different theta bounds would cut across the
        concave notch between them, silently claiming sphere area this plate doesn't cover --
        fatal once `deform()` uses this same outline for its own contested/open
        classification; see `LithospherePlate.deform` and the no-node-inside-a-neighbour's-
        polygon invariant test in `unit_tests/test_plates.py` / `stress_tests/
        test_world_stepping.py`).

        A row split into two arcs by interior subduction (see
        `_grow_or_shrink_line_for_deform`), or by a split/defragment partition, contributes
        both intervals; the gap between them is a genuine hole in this plate's territory
        (a neighbour's lobe that punched through the middle of it, or a stranded sibling),
        traced as its own loop wound opposite to the outer boundary and stitched in via a
        zero-width keyhole seam (`_plate_outline_loops` / `_stitch_loops`) so the winding-
        number test reads it as outside the plate. Still one plain `(n, 3)` array, so every
        `get_bounding_polygon()` consumer is unchanged."""
        lines_with_nodes = [line for line in self._lines if len(line) > 0]
        if not lines_with_nodes:
            return np.zeros((0, 3))
        loop = _stitch_loops(_plate_outline_loops(_row_intervals(lines_with_nodes)))
        if len(loop) < 3:
            return np.zeros((0, 3))
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
        rows = _row_intervals(lines_with_nodes)
        phis = np.array([phi for phi, _ in rows])
        low_thetas = np.array([ivs[0][0] for _, ivs in rows])
        high_thetas = np.array([ivs[-1][1] for _, ivs in rows])
        max_arcs = max(len(ivs) for _, ivs in rows)
        interval_lo = np.full((len(rows), max_arcs), np.inf)
        interval_hi = np.full((len(rows), max_arcs), -np.inf)
        for i, (_, ivs) in enumerate(rows):
            for k, (lo, hi) in enumerate(ivs):
                interval_lo[i, k] = lo
                interval_hi[i, k] = hi

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
            interval_lo=interval_lo,
            interval_hi=interval_hi,
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
        phis = row_lookup.phis

        local_xyz = geometry.to_local(self._frame, points_xyz)
        phi_q, theta_q = geometry.xyz_to_latlon(local_xyz)
        idx = np.searchsorted(phis, phi_q)
        idx_lo = np.clip(idx - 1, 0, len(phis) - 1)
        idx_hi = np.clip(idx, 0, len(phis) - 1)
        nearer_to_lo = np.abs(phi_q - phis[idx_lo]) <= np.abs(phis[idx_hi] - phi_q)
        nearest = np.where(nearer_to_lo, idx_lo, idx_hi)

        # theta_q inside *any* of the nearest row's arc intervals (usually one; two+ when a
        # row was split by interior subduction). Absent arc slots are +inf/-inf padded, so
        # they never match -- a point in the gap between two arcs reads as outside here and,
        # if beyond the bulge margin, as definitely outside below.
        lo_sel = row_lookup.interval_lo[nearest]
        hi_sel = row_lookup.interval_hi[nearest]
        in_any_interval = np.any((theta_q[:, None] >= lo_sel) & (theta_q[:, None] <= hi_sel), axis=1)
        idealized_inside = (phi_q >= phis[0]) & (phi_q <= phis[-1]) & in_any_interval
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

    def _get_world_points(self) -> np.ndarray:
        """Every non-empty line's node world-xyz `(n, 3)`, concatenated in line order --
        cached (see `_world_points_cache`). Rebuilt with a single `local_xyz` +
        frame-rotation over the whole plate's `(phi, theta)` rather than a per-line pair of
        small numpy calls, which is what drove the profiled `world_xyz`/`latlon_to_xyz`
        call counts (see docs/profiling.md). Read-only for callers."""
        if self._world_points_cache is None:
            lines = [line for line in self._lines if len(line) > 0]
            if not lines:
                self._world_points_cache = np.zeros((0, 3))
            else:
                theta = np.concatenate([line.theta for line in lines])
                phi = np.repeat(
                    np.array([line.phi for line in lines], dtype=float),
                    [len(line) for line in lines],
                )
                self._world_points_cache = geometry.to_world(self._frame, geometry.local_xyz(phi, theta))
        return self._world_points_cache

    def all_points_and_elevation(self) -> tuple[np.ndarray, np.ndarray]:
        """Every elevation-line node's world position and elevation, concatenated. The
        positions come from `_get_world_points()`'s cache (read-only, like
        `get_bounding_polygon()`); elevation is gathered fresh every call since it changes
        without a node-set mutation."""
        return self._get_world_points(), self.collect("elevation")

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

    def _plates_from_node_masks(self, masks: list[np.ndarray], ids: list[int]) -> list["Plate"]:
        """Partition this plate's lines by per-node membership (same `ElevationLine.masked`
        machinery as `split`, so every `OPTIONAL_FIELDS` value survives exactly -- no
        resample). `ids[0]` keeps this plate's own id, omega, and age; the rest are fresh
        fragments carrying a copy of this plate's omega (they were co-moving with it, which
        is exactly why `maybe_split_plate` never separated them) and age 0. `type(self)` so a
        `LithospherePlate` stays a `LithospherePlate`."""
        plates: list["Plate"] = []
        for k, (mask, pid) in enumerate(zip(masks, ids)):
            offset = 0
            lines: list[ElevationLine] = []
            for line in self._lines:
                n = len(line)
                sub = mask[offset : offset + n]
                offset += n
                if np.any(sub):
                    # A single connected component can still wrap a row into two arcs (a
                    # U-shape closed through other rows) -- carry each arc as its own
                    # contiguous `ElevationLine` so the fragment's envelope keyholes the gap
                    # out rather than claiming it (see `split_into_contiguous_runs`).
                    lines.extend(split_into_contiguous_runs(line.masked(sub), _row_median_step(line)))
            if not lines:
                continue
            # A fragment's own crust_type is the majority of what its nodes actually are, not
            # a blind copy of the parent's -- see elevation_lines.majority_crust_type. A no-op
            # for every plate that's never had a magma-typing event (rift decompression
            # melting or gap-fill), which is every fragment before that feature existed.
            fragment_crust_type = elevation_lines.majority_crust_type(lines, self._crust_type)
            plates.append(
                type(self)(
                    plate_id=pid,
                    frame=self._frame.copy(),
                    crust_type=fragment_crust_type,
                    lines=lines,
                    omega=self._omega.copy(),
                    age_steps=self._age_steps if k == 0 else 0,
                    internal_stress=self._internal_stress if k == 0 else 0.0,
                )
            )
        return plates


# outline_world's boundary-detection pass needs at least this many nodes for "boundary node"
# to be a meaningful distinct subset of "every node" -- below it, every node is returned as-is
# rather than running (and likely degenerating) the density check.
OUTLINE_MIN_NODES_FOR_HULL = 4


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


# Two plates' node clouds are "co-located" -- overlapping the same patch of sphere rather
# than merely adjacent -- when nodes land within this multiple of a target spacing of each
# other. Ordinary shared boundaries sit ~one full spacing apart, so half a spacing only fires
# on genuine territory overlap (a stalled collision, a bad split partition, a plate drifting
# over a neighbour it can't merge with). Shared by main._plate_overlaps (the Plate Inspector /
# diagnostics view) and merge_split.update_overlap_tracking (the per-node onset stamp).
OVERLAP_TOLERANCE_MULT = 0.5


def compute_node_overlap(plate_list: list[Plate], tol_rad: float) -> dict[int, dict]:
    """Per plate (keyed by plate_id, only plates with nodes), a genuine node-cloud overlap
    read against every *other* plate:

    - `overlap_mask`: bool array aligned to this plate's own node order
      (`all_points_and_elevation()` / `collect` order) -- True where this node sits within
      `tol_rad` of some other plate's node.
    - `by_partner`: {other_plate_id: count of this plate's own unique nodes on top of it},
      sorted-desc when iterated is up to the caller.

    One global `cKDTree.query_pairs` over every node, so O(N log N) once rather than a
    per-pair envelope test -- the same construction main._plate_overlaps used inline before
    this was factored out so the API view and merge_split's onset tracker can't drift."""
    active = [p for p in plate_list if p.node_count() > 0]
    result: dict[int, dict] = {
        p.plate_id: {"overlap_mask": np.zeros(p.node_count(), dtype=bool), "by_partner": {}} for p in active
    }
    if len(active) < 2:
        return result

    clouds = [p.all_points_and_elevation()[0] for p in active]
    counts = [len(c) for c in clouds]
    offsets = np.cumsum([0, *counts])
    owner = np.concatenate([np.full(n, i) for i, n in enumerate(counts)])
    pairs = cKDTree(np.concatenate(clouds)).query_pairs(tol_rad, output_type="ndarray")
    if len(pairs) == 0:
        return result

    owners_lo, owners_hi = owner[pairs[:, 0]], owner[pairs[:, 1]]
    cross = owners_lo != owners_hi
    pairs, owners_lo, owners_hi = pairs[cross], owners_lo[cross], owners_hi[cross]

    for glob, src, dst in ((pairs[:, 0], owners_lo, owners_hi), (pairs[:, 1], owners_hi, owners_lo)):
        for i, src_plate in enumerate(active):
            here = src == i
            if not here.any():
                continue
            local = glob[here] - offsets[i]
            result[src_plate.plate_id]["overlap_mask"][local] = True
            dst_here = dst[here]
            for j, dst_plate in enumerate(active):
                on_j = dst_here == j
                if not on_j.any():
                    continue
                result[src_plate.plate_id]["by_partner"][dst_plate.plate_id] = int(len(np.unique(local[on_j])))
    return result


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


def collect_all_volcano_active_years_remaining(plate_list: list[Plate]) -> np.ndarray:
    """Per-node countdown to dormancy for volcano nodes (0 elsewhere) -- see volcanism.py.
    Used by GET /world/volcanoes to flag which volcanoes are still erupting-capable."""
    return _collect_all(plate_list, "volcano_active_years_remaining")


def collect_all_overlap_onset_years(plate_list: list[Plate]) -> np.ndarray:
    """Used by render_image.py's `overlapAge` debug view -- see
    merge_split.update_overlap_tracking / ElevationLine.overlap_onset_years."""
    return _collect_all(plate_list, "overlap_onset_years")


def collect_all_elevation(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "elevation")


def collect_all_crustal_thickness(plate_list: list[Plate]) -> np.ndarray:
    """Every node's Hc (v2 LithospherePlate only -- all-zero for v1 PlateWithLines). Read by
    erosion.py so a step's net rock-column change lands on Hc and Airy isostasy can rebound
    it, not just on the bare `elevation` cache."""
    return _collect_all(plate_list, "crustal_thickness_m")


def collect_all_mantle_lithosphere_thickness(plate_list: list[Plate]) -> np.ndarray:
    """Every node's Hm (v2 LithospherePlate only -- all-zero for v1). Paired with
    `collect_all_crustal_thickness` for erosion.py's isostatic-rebound bookkeeping."""
    return _collect_all(plate_list, "mantle_lithosphere_thickness_m")


def collect_all_elev_change_reason(plate_list: list[Plate]) -> np.ndarray:
    """Every node's elevation-change provenance code (see elevation_lines.ELEV_CHANGE_*) --
    read by erosion.py to preserve a quiescent node's older provenance, and by
    render_image.py's "elevReason" debug view."""
    return _collect_all(plate_list, "elev_change_reason")


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


# Categories for render_image.py's "crustType" debug view -- resolved per-plate (unlike the
# plain _collect_all helpers above, a raw crust_type_code is meaningless without its owning
# plate's own crust_type to resolve CRUST_TYPE_INHERIT against, see
# elevation_lines.effective_is_continental_from_codes), so this can't just be another
# `_collect_all(plate_list, "crust_type_code")` one-liner.
CRUST_TYPE_VIEW_OCEANIC = 0  # oceanic, matches its plate's own nominal type (the common case)
CRUST_TYPE_VIEW_CONTINENTAL = 1  # continental, matches its plate's own nominal type
CRUST_TYPE_VIEW_OCEANIC_ANOMALY = 2  # oceanic node on a nominally continental plate
CRUST_TYPE_VIEW_CONTINENTAL_ANOMALY = 3  # continental node on a nominally oceanic plate


def collect_all_crust_type_view_codes(plate_list: list[Plate]) -> np.ndarray:
    """Every node's CRUST_TYPE_VIEW_* category -- see the constants above. The two "anomaly"
    categories are exactly the nodes a magma-typing event (rift decompression melting,
    gap-fill) stamped a different composition than the plate they sit on, e.g. a drowned
    continental margin that finally melted through to real oceanic crust, or a volcanic island
    breaching the surface on an oceanic plate -- everywhere else this is a no-op reading of
    the plate's own nominal crust_type."""
    chunks = []
    for p in plate_list:
        n = p.node_count()
        if n == 0:
            continue
        plate_is_continental = p.crust_type == "continental"
        is_continental = elevation_lines.effective_is_continental_from_codes(p.collect("crust_type_code"), plate_is_continental)
        anomaly = is_continental != plate_is_continental
        codes = np.where(
            anomaly,
            CRUST_TYPE_VIEW_OCEANIC_ANOMALY if plate_is_continental else CRUST_TYPE_VIEW_CONTINENTAL_ANOMALY,
            CRUST_TYPE_VIEW_CONTINENTAL if plate_is_continental else CRUST_TYPE_VIEW_OCEANIC,
        )
        chunks.append(codes)
    return np.concatenate(chunks, axis=0) if chunks else np.zeros(0, dtype=np.int8)


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


def _land_noise_threshold(
    owner_tree: cKDTree,
    crust_types: list[str],
    noise: "terrain_noise.ReliefField",
    land_fraction: float,
    sealevel_noise_offset: float = 0.0,
) -> float | None:
    """Translate a requested whole-sphere land_fraction into a concrete noise threshold for
    continental crust's elevation formula (each caller applies this threshold in its own
    per-node elevation/thickness formula).

    `noise` is anything with a `sample(xyz) -> array` method -- a bare `SphereNoise` or
    `terrain_noise.ContinentalRelief` (whose `sample()` is the land/sea-deciding component,
    deliberately kept at the same low-frequency character so this coarse quantile stays a
    good estimator; its orogenic `uplift()` is added elsewhere and never crosses sea level).

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
    to place land on.

    `sealevel_noise_offset` corrects for the reference continental column not sitting
    exactly at sea level: a node is land when its noise value exceeds `threshold +
    offset`, not `threshold` (the caller knows `offset` -- it is `(Hc_at_sealevel - Hc0) /
    amplitude` through the isostasy formula, a small negative number). Subtracting it here
    means `quantile(1 - target)` lands on the actual land/sea crossing, so the measured
    land fraction tracks the request instead of overshooting it. Default 0.0 keeps the
    old behaviour for a caller that adds the noise straight onto elevation."""
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
    return float(np.quantile(continental_noise, 1.0 - target_sub_fraction)) - sealevel_noise_offset


