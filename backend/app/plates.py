"""Plates: identity, motion, territory, and the `PlateSurface` contract.

Each plate owns a rotation matrix (`frame`) mapping its plate-local coordinates to world unit
vectors (see `geometry.plate_frame_from_seed`). Rotating a plate rigidly only ever updates
`frame` -- the plate-local node positions themselves never change, so rotation never needs
resampling. See docs/simulation-model.md for the full design writeup, and
sparse_quad_patch.py for the surface representation itself (`PlateWithSparseQuadPatch`).
"""

from __future__ import annotations

import abc
from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterator, Mapping, Protocol

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from . import ellipse, geometry, healpix_grid
from . import elevation_lines
from .elevation_lines import (
    PLANET_RADIUS_KM,
    TARGET_LINE_SPACING_RAD,
    ElevationPoint,
    iter_local_lattice,
)
from .surface_fields import SURFACE_FIELDS, SurfaceField

if TYPE_CHECKING:
    from . import terrain_noise
    from .world import World

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


@dataclass(frozen=True)
class SurfaceNodes:
    """Representation-neutral bulk view of a plate surface's live nodes.

    All arrays share one stable ordering for the lifetime of ``topology_revision``. The
    position arrays are ``(n, 3)`` unit vectors; ``node_ids`` is an ``(n, 2)`` opaque uint64
    identity; ``area_m2`` and every requested field are ``(n,)`` arrays. ``area_is_exact``
    is true for cell-backed areas (every current surface). Callers
    must treat this container as read-only and use ``set_fields_on_plate`` for write-back.
    """

    local_xyz: np.ndarray
    world_xyz: np.ndarray
    node_ids: np.ndarray
    area_m2: np.ndarray
    area_is_exact: bool
    fields: Mapping[str, np.ndarray]


@dataclass(frozen=True)
class SurfaceAdjacency:
    """CSR adjacency over the node ordering used by :class:`SurfaceNodes`."""

    offsets: np.ndarray
    neighbours: np.ndarray


class PlateSurface(abc.ABC):
    """Storage-neutral compatibility boundary for terrain carried by a plate.

    Phase 1 deliberately describes capabilities, not rows, vertices, or quads. Topology
    operations and conservative remapping remain implementation-specific until later phases.
    """

    @property
    @abc.abstractmethod
    def topology_revision(self) -> int: ...

    @property
    @abc.abstractmethod
    def geometry_revision(self) -> int: ...

    @abc.abstractmethod
    def surface_nodes(self, *field_names: str) -> SurfaceNodes: ...

    @abc.abstractmethod
    def set_fields_on_plate(self, **fields: np.ndarray) -> None: ...

    @abc.abstractmethod
    def adjacency(self) -> SurfaceAdjacency: ...

    @abc.abstractmethod
    def boundary_loops_world(self) -> tuple[np.ndarray, ...]: ...

    @abc.abstractmethod
    def contains_batch(self, points_xyz: np.ndarray) -> np.ndarray: ...

    def field_metadata(self) -> Mapping[str, SurfaceField]:
        """The complete persistent-field registry shared by every representation."""
        return SURFACE_FIELDS

    def node_index_for_id(self, node_id: np.ndarray | tuple[int, int]) -> int | None:
        """Resolve an opaque surface node ID in the current topology revision."""
        target = np.asarray(node_id, dtype=np.uint64)
        if target.shape != (2,):
            raise ValueError("surface node ID must contain exactly two uint64 words")
        matches = np.flatnonzero(np.all(self.surface_nodes().node_ids == target, axis=1))
        return None if len(matches) == 0 else int(matches[0])


# Two plates count as neighbours once the closest points of their two outlines come within
# this many multiples of the default line spacing -- generous enough to still catch plates
# separated by a boundary's own elevation-transition zone (see boundary.py's
# FAR_THRESHOLD_RAD, a similar multiple of spacing_rad), while not matching plates that
# merely share the same hemisphere.
NEIGHBOUR_DISTANCE_RAD = 6.0 * TARGET_LINE_SPACING_RAD


def _plates_within(plate: "Plate", all_plates: list["Plate"], threshold_rad: float) -> list["Plate"]:
    """Every other plate in `all_plates` whose outline (`Plate.get_bounding_polygon`) comes
    within `threshold_rad` of `plate`'s own -- backs `get_neighbours`. Uses the cached get_bounding_polygon() rather than outline_world()
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
    `contains_batch` instead of a shared polygon-list winding test -- each neighbour answers
    with its own exact cell lookup rather than paying the full winding-number cost. Same
    semantics otherwise: stops early once every point is already
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


class Plate(PlateSurface, abc.ABC):
    """A plate's shared identity/motion state plus an abstract interface over however it
    represents its own terrain nodes (`PlateWithSparseQuadPatch`, see sparse_quad_patch.py).
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
        self._topology_revision = 0
        self._geometry_revision = 0
        # Lazily (re)computed by get_bounding_polygon() below -- None means "stale, recompute
        # on next call," not "empty polygon" (an empty plate's real outline is a valid
        # np.zeros((0, 3)), which must stay distinguishable from "not computed yet").
        # Invalidated by rotate() here and by the subclass's own topology changes --
        # elevation-only mutations (erosion, uplift, ...) don't touch outline_world's inputs,
        # so they leave the cache untouched.
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
        self._geometry_revision += 1
        self._invalidate_bounding_polygon()

    def shift(self, world: "World", years: float) -> float:  # noqa: F821
        """Integrate this plate's rigid motion through the shared torque model.

        Torque inputs consume only the representation-neutral surface API, so rigid motion
        belongs here rather than in the tectonic engine. Topology remains fixed by this
        operation: only the frame and geometry revision move.
        """
        # Local import avoids plates <-> torque's module-level dependency cycle.
        from . import torque

        # Debug-world-only escape hatch (see World.pinned_omegas): scripted scenarios move
        # exactly as configured rather than through the real torque balance.
        pinned = world.pinned_omegas.get(self.plate_id)
        if pinned is not None:
            old_points, _ = self.all_points_and_elevation()
            if len(old_points) == 0:
                return 0.0
            return torque.apply_omega_and_rotate(
                self, old_points, np.asarray(pinned, dtype=float), years
            )
        other_plates = [plate for plate in world.plates if plate.plate_id != self.plate_id]
        return torque.shift_plate(self, world, other_plates, years)

    @property
    def topology_revision(self) -> int:
        return self._topology_revision

    @property
    def geometry_revision(self) -> int:
        return self._geometry_revision

    def age_one_step(self) -> None:
        self._age_steps += 1

    def reset_age(self) -> None:
        self._age_steps = 0

    def set_age_steps(self, age_steps: int) -> None:
        self._age_steps = age_steps

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
        gives, rather than raising."""
        defaults = {"_internal_stress": 0.0, "_topology_revision": 0, "_geometry_revision": 0}
        if name in defaults:
            object.__setattr__(self, name, defaults[name])
            return defaults[name]
        raise AttributeError(name)

    @abc.abstractmethod
    def node_count(self) -> int: ...

    def has_negligible_territory(self) -> bool:
        """True once this plate has been whittled down to no real remaining territory --
        the "no land left" half of merge_split.remove_defunct_plates/apply_topology_changes'
        own pruning (the other half, `node_count() == 0`, is checked separately by both
        callers). Default: fewer nodes than can form a real 2D hull (see
        `OUTLINE_MIN_NODES_FOR_HULL`)."""
        return self.node_count() < OUTLINE_MIN_NODES_FOR_HULL

    def defragment(
        self, next_id: int, connect_radius_rad: float, min_fragment_nodes: int, world: "World"  # noqa: F821
    ) -> tuple[list["Plate"], int, list["Plate"]] | None:
        """Reconcile "one `Plate` object" with "one contiguous patch of crust."

        Per-step `deform()` adds and removes cells only at the boundary, so
        subduction/transform can carve a plate's cells into two disconnected landmasses, or
        strand a few cells far from the plate body, and nothing notices: `maybe_split_plate` only cuts on mantle-
        flow *disagreement*, not geometry, and two co-moving lobes never trip it.

        This finds those cases directly. Connected components of this plate's nodes at
        `connect_radius_rad` (see `node_components`); each component with at least
        `min_fragment_nodes` nodes becomes its own plate (the largest keeps this plate's own
        id/frame/omega/age -- see `_plates_from_node_masks`), everything smaller is cut off
        as stranded crust, one fragment per component.

        Returns `None` -- nothing to do -- when the plate is already a single contiguous
        patch (the overwhelmingly common case), when it has exactly one component big
        enough to anchor a plate and no stranded nodes to shed, or when it's debris with no
        component large enough to anchor a plate at all (left for `has_negligible_territory`
        / `remove_defunct_plates` to prune). Otherwise returns
        `(replacement_plates, n_new_ids_consumed, stranded_fragments)`, where
        `replacement_plates[0]` reuses this plate's own id and `next_id, next_id + 1, ...` are
        consumed for the rest, in descending component-size order. `stranded_fragments` are
        the cut-off components as plates that aren't live: they carry this plate's id, omega
        and frame and age 0, for the caller to place or drop (see
        `merge_split.defragment_plates`). `next_id` is `World.next_plate_id`. `world` is used
        only to record any stranded/dropped nodes into `World.removed_points_log` (see the
        "Added/Removed Points" debug view) -- never mutated otherwise."""
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
        dropped_mask = ~np.isin(labels, kept)
        dropped_nodes = int(dropped_mask.sum())
        # No component big enough to anchor a plate -- the whole thing is debris. Leave it
        # for merge_split.remove_defunct_plates / has_negligible_territory to prune; defrag
        # never deletes a whole plate itself (that path is fragile against small synthetic
        # plates and adds nothing the negligible-territory check doesn't already do).
        if not kept:
            return None
        if len(kept) == 1 and dropped_nodes == 0:
            return None

        if dropped_nodes > 0:
            world.record_removed_points(points[dropped_mask], self.plate_id)

        n_new_ids = len(kept) - 1
        stranded = [int(cid) for cid in component_ids if int(cid) not in kept]
        masks = [labels == cid for cid in [*kept, *stranded]]
        ids = [self.plate_id, *range(next_id, next_id + n_new_ids), *[self.plate_id] * len(stranded)]
        plates = self._plates_from_node_masks(masks, ids)
        return plates[: len(kept)], n_new_ids, plates[len(kept) :]

    @abc.abstractmethod
    def _plates_from_node_masks(self, masks: list[np.ndarray], ids: list[int]) -> list["Plate"]:
        """Build one plate per mask in `masks` (each a boolean array over this plate's nodes
        in `all_points_and_elevation` order), assigning `ids[k]` to `masks[k]`'s plate --
        `ids[0]` is always this plate's own id, so the first mask should be the one that
        keeps this plate's identity."""
        ...

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

    @abc.abstractmethod
    def collect(self, field_name: str) -> np.ndarray:
        """Every node's current `field_name` value (any `surface_fields.SURFACE_FIELDS` name),
        concatenated in this plate's own node order. Empty
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
        """Every node this plate owns, as `ElevationPoint`s, in this plate's own node order.
        Lets code that just wants "every node, read or write" (not a bulk array op) work
        without reaching into the representation's own storage."""
        ...

    @abc.abstractmethod
    def map_world_points(self) -> Iterator[tuple[ElevationPoint, np.ndarray]]:
        """Every node this plate owns, paired with its own world xyz position -- the same
        nodes/order as `__iter__`, just with each `ElevationPoint` accompanied by the world
        coordinate a caller would otherwise have had to derive itself. Each `ElevationPoint` is
        a live view, so a value computed as a function of world position (noise, distance,
        sampled field) can be written straight back with the point's own `set_*` -- in place,
        no topology round-trip needed."""
        ...

    @abc.abstractmethod
    def set_fields_on_plate(self, **fields: np.ndarray) -> None:
        """Bulk in-place write for any `surface_fields.SURFACE_FIELDS` name: each keyword's array must be exactly this plate's own node count, in the same
        order `map_world_points`/`collect` already read/traverse it in. The
        vectorized counterpart to looping `map_world_points` and calling each
        point's own `set_*` -- for a caller that already has a full per-node array computed
        (erosion/bathymetry/geology's per-step recompute), this writes it back without
        constructing a `Plate.__iter__`-style point object per node."""
        ...


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

    `plates_in_order` is what makes this representation-agnostic: any bulk per-field gather (`collect_all_elevation` and
    friends, below) or per-plate write-back loop (`Plate.set_fields_on_plate`) already
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


def cached_node_position_tree(world: "World | None", points: np.ndarray) -> cKDTree:  # noqa: F821
    """A `cKDTree` over `points` -- this step's full node cloud, as returned by
    `gather_node_positions` -- shared via `World.node_position_tree_cache` across every
    caller that runs after this step's shift/deform/topology changes have already settled
    (see that field's own docstring on `World` for exactly which window that is -- climate.py's
    `_sample_elevation_and_crust` already established the pattern this generalizes). Every
    such caller is handed the *same* `points` array this step (the one `world.py`'s
    `step_world` gathers once via `gather_node_positions` and threads through
    `climate.compute_climate`/`hydrology.compute_hydrology`/`erosion.apply_erosion`), so
    whichever of them runs first this step pays the build and every later one reuses it.
    `world=None` -- a direct unit-test call against a bare points array, not a real step --
    always builds fresh instead of risking a stale hit against some other test's leftover
    cache."""
    if world is None:
        return cKDTree(points)
    cached = world.node_position_tree_cache
    if cached is not None:
        return cached[1]
    tree = cKDTree(points)
    world.node_position_tree_cache = (points, tree)
    return tree


def cached_node_healpix_index(world: "World | None", points: np.ndarray) -> healpix_grid.NodePixelIndex:  # noqa: F821
    """The `"healpix"` `World.node_cloud_resample_mode` analogue of `cached_node_position_tree`
    immediately above -- shared via `World.node_healpix_grid_cache`/`node_healpix_index_cache`
    across every this-step-node-cloud caller under that mode
    (`render_image._node_cloud_and_tree`, `climate._sample_elevation_and_crust` -- issue #133
    phase 2), the same split phase 1 established: the `HealpixGrid` itself
    (`node_healpix_grid_cache`, keyed by `nside`, a pure function of node *count*) is reused for
    a world's whole life unless that count crosses a bracket, while the per-step scatter+fill
    (`node_healpix_index_cache`) is rebuilt whenever `step_world` resets it (same invalidation
    event as `node_position_tree_cache`). `world=None` (direct unit-test call) always builds
    fresh, matching `cached_node_position_tree`'s own testing convention."""
    if world is None:
        nside = healpix_grid.nside_for_node_count(points.shape[0])
        return healpix_grid.build_node_pixel_index(healpix_grid.build(nside), points)
    cached_index = world.node_healpix_index_cache
    if cached_index is not None:
        return cached_index
    nside = healpix_grid.nside_for_node_count(points.shape[0])
    cached_grid = world.node_healpix_grid_cache
    if cached_grid is not None and cached_grid[0] == nside:
        grid = cached_grid[1]
    else:
        grid = healpix_grid.build(nside)
        world.node_healpix_grid_cache = (nside, grid)
    index = healpix_grid.build_node_pixel_index(grid, points)
    world.node_healpix_index_cache = index
    return index


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
# on genuine territory overlap. Used by faults.py's point-overlap fault spawning.
OVERLAP_TOLERANCE_MULT = 0.5


def compute_node_overlap(plate_list: list[Plate]) -> dict[int, dict]:
    """Per plate (keyed by plate_id, only plates with nodes), a genuine territory overlap
    read against every *other* plate:

    - `overlap_mask`: bool array aligned to this plate's own node order
      (`all_points_and_elevation()` / `collect` order) -- True where this node's centre lies
      inside some other plate's territory.
    - `by_partner`: {other_plate_id: count of this plate's own unique nodes on top of it},
      sorted-desc when iterated is up to the caller.
    - `cover_count`: int array aligned like `overlap_mask` -- how many *other* plates this
      node sits on (`overlap_mask` is `cover_count > 0`). Whole-world area sums divide by
      `1 + cover_count` to count overlapped ground once (issue #289).
    - `continental_cover_count`: the same, counting only other plates whose `crust_type` is
      continental -- the denominator for a continental-plates-only sum, where an oceanic
      plate on top duplicates nothing.

    Containment, not node proximity: the same test deform's contested classification uses.
    Two cell lattices in different frames never line up node for node, so a proximity
    tolerance misses about a third of the nodes that really sit on another plate (issue #228
    Phase 4). Shared by main._plate_overlaps (the Plate Inspector / diagnostics view) and
    merge_split.update_overlap_tracking (the per-node onset stamp), so the two can't drift."""
    active = [p for p in plate_list if p.node_count() > 0]
    result: dict[int, dict] = {
        p.plate_id: {
            "overlap_mask": np.zeros(p.node_count(), dtype=bool),
            "by_partner": {},
            "cover_count": np.zeros(p.node_count(), dtype=np.int64),
            "continental_cover_count": np.zeros(p.node_count(), dtype=np.int64),
        }
        for p in active
    }
    if len(active) >= 2:
        _contained_node_overlap(active, result)
    return result


def _contained_node_overlap(active: list[Plate], result: dict[int, dict]) -> None:
    """`compute_node_overlap`'s containment test, filling `result` in place. Candidate
    pairs come from bounding caps rather than outline proximity (`_plates_within`), so a plate
    buried wholly inside another -- the superimposed case the forced merge exists for -- is
    still tested."""
    clouds = [p.all_points_and_elevation()[0] for p in active]
    caps = [geometry.bounding_sphere(c) for c in clouds]
    # Cell centres sit up to a cell's half-diagonal inside the territory's edge, so each cap is
    # padded by one of its plate's largest cells before a pair is ruled out.
    pads = [float(np.sqrt(np.max(p.surface_nodes().area_m2))) / (elevation_lines.PLANET_RADIUS_KM * 1000.0) for p in active]
    for i, plate in enumerate(active):
        centre_i, radius_i = caps[i]
        for j, other in enumerate(active):
            if i == j:
                continue
            centre_j, radius_j = caps[j]
            if float(geometry.angular_distance(centre_i, centre_j)) > radius_i + radius_j + pads[j]:
                continue
            inside = other.contains_batch(clouds[i])
            if np.any(inside):
                result[plate.plate_id]["overlap_mask"] |= inside
                result[plate.plate_id]["cover_count"] += inside
                if other.crust_type == "continental":
                    result[plate.plate_id]["continental_cover_count"] += inside
                result[plate.plate_id]["by_partner"][other.plate_id] = int(np.count_nonzero(inside))


def _collect_all(plate_list: list[Plate], field_name: str) -> np.ndarray:
    """Every plate's current `field_name` (any `surface_fields.SURFACE_FIELDS` name),
    concatenated in the exact same per-plate/per-node order collect_all_points uses --
    so results from two different `_collect_all` calls can still be indexed together with the
    same nearest-neighbor result (see render_image._render_grid_arrays). Delegates to each
    plate's own `collect`, through the abstract `Plate` interface."""
    chunks = [p.collect(field_name) for p in plate_list if p.node_count() > 0]
    if not chunks:
        return np.zeros(0, dtype=bool) if field_name == "is_volcano" else np.zeros(0)
    return np.concatenate(chunks, axis=0)


def collect_all_lake_depth(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "lake_depth")


def collect_all_glacier_depth(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "glacier_depth")


def collect_all_ice_load_deflection(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "ice_load_deflection_m")


def collect_all_breach_notch_depth(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "breach_notch_depth_m")


def collect_all_channel_reference_elevation(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "channel_reference_elevation_m")


def collect_all_silt_depth(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "silt_depth")


def collect_all_channel_depth(plate_list: list[Plate]) -> np.ndarray:
    """Used by climate.py to size a river's own evaporative surface for its moisture-
    recycling humidity source (see that module), and by render_image.py's
    `_hillshade_for_world` to carve a wide-enough channel's own incision into the elevation
    hillshade lights."""
    return _collect_all(plate_list, "channel_depth")


def collect_all_channel_width(plate_list: list[Plate]) -> np.ndarray:
    """Used by render_image.py to draw a wide river's line thicker and more strongly tinted
    toward river-blue than a narrow one (`_rivers_to_draw`/`_draw_rivers`), and to gate how much
    of a channel's depth actually incises the hillshade relief (`_hillshade_for_world`) -- a
    channel too narrow for either purpose reads as a barely-there creek instead."""
    return _collect_all(plate_list, "channel_width")


def collect_all_is_volcano(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "is_volcano")


def collect_all_volcano_active_years_remaining(plate_list: list[Plate]) -> np.ndarray:
    """Per-node countdown to dormancy for volcano nodes (0 elsewhere) -- see volcanism.py.
    Used by GET /world/volcanoes to flag which volcanoes are still erupting-capable."""
    return _collect_all(plate_list, "volcano_active_years_remaining")


def collect_all_overlap_onset_years(plate_list: list[Plate]) -> np.ndarray:
    """Used by render_image.py's `overlapAge` debug view -- see
    merge_split.update_overlap_tracking / the `overlap_onset_years` field."""
    return _collect_all(plate_list, "overlap_onset_years")


def collect_all_node_created_years(plate_list: list[Plate]) -> np.ndarray:
    """Every live node's `world.elapsed_years` at creation (-1.0 = predates tracking / a
    genesis node from initial world generation) -- used by render_image.py's `nodeAge` debug
    view. See the `node_created_years` field."""
    return _collect_all(plate_list, "node_created_years")


def collect_all_accounting_areas_m2(plate_list: list[Plate], spacing_rad: float) -> np.ndarray:
    """Every node's `Plate.accounting_areas_m2`, in the same order as the other
    `collect_all_*` gathers -- the per-node weight for any whole-world area or volume sum."""
    chunks = [p.accounting_areas_m2(spacing_rad) for p in plate_list if p.node_count() > 0]
    return np.concatenate(chunks) if chunks else np.zeros(0)


def collect_all_elevation(plate_list: list[Plate]) -> np.ndarray:
    return _collect_all(plate_list, "elevation")


def collect_all_crustal_thickness(plate_list: list[Plate]) -> np.ndarray:
    """Every node's Hc. Read by
    erosion.py so a step's net rock-column change lands on Hc and Airy isostasy can rebound
    it, not just on the bare `elevation` cache."""
    return _collect_all(plate_list, "crustal_thickness_m")


def collect_all_mantle_lithosphere_thickness(plate_list: list[Plate]) -> np.ndarray:
    """Every node's Hm. Paired with
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


def nearest_node_index(plate_list: list[Plate], query_xyz: np.ndarray) -> int | None:
    """The index of the single node nearest `query_xyz`, in the same concatenated per-plate/
    per-node order every `collect_all_*` function uses -- so this index can be reused directly
    against `collect_all_points`/`collect_all_elevation`/`collect_all_node_created_years`/etc
    without a second lookup. `GET /world/node_at`'s click-to-inspect hit-test (see
    `docs/debugging.md`'s "Added/Removed Points" view) -- `nearest_plate_id` above answers
    "which plate" alone; this answers "which node exactly," so a caller wanting both per-node
    fields *and* plate ownership needn't repeat the same k-d tree build twice. `None` if every
    plate is empty."""
    collected = collect_all_points(plate_list)
    if collected is None:
        return None
    points, _, _owner = collected
    _, idx = cKDTree(points).query(query_xyz)
    return int(idx)


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
