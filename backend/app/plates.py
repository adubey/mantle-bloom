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
from typing import Iterator, Protocol

import numpy as np
from scipy.spatial import ConvexHull, QhullError, cKDTree

from . import ellipse, geometry
from .elevation_lines import (
    DEFAULT_NODE_DENSITY,
    NODE_DENSITY_CHOICES,
    PLANET_RADIUS_KM,
    TARGET_LINE_SPACING_KM,
    TARGET_LINE_SPACING_RAD,
    ElevationLine,
    ElevationPoint,
    build_lines_from_lattice,
    install_point_field_accessors,
    iter_local_lattice,
    line_spacing_rad,
)
from .noise import SphereNoise
from .rtree_index import RTree

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
    """Every other plate in `all_plates` whose outline (`Plate.outline_world`) comes within
    `threshold_rad` of `plate`'s own -- shared by both PlateWithLines.get_neighbours and
    PlateWithRTree.get_neighbours, which differ only in what outline_world() itself does for
    that representation. Same bounding-sphere prefilter as boundary.py's step_boundaries
    (cheap enough to compute per plate, and enough to rule out most pairs before a real
    nearest-point query)."""
    own_points = plate.outline_world()
    if len(own_points) == 0:
        return []
    own_centroid, own_radius = geometry.bounding_sphere(own_points)

    neighbours = []
    for other in all_plates:
        if other.plate_id == plate.plate_id:
            continue
        other_points = other.outline_world()
        if len(other_points) == 0:
            continue
        other_centroid, other_radius = geometry.bounding_sphere(other_points)
        centroid_dist = float(geometry.angular_distance(own_centroid, other_centroid))
        if centroid_dist - own_radius - other_radius > threshold_rad:
            continue
        closest_dist = float(cKDTree(other_points).query(own_points)[0].min())
        if closest_dist <= threshold_rad:
            neighbours.append(other)
    return neighbours


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

    @property
    def lines(self) -> tuple[ElevationLine, ...]:
        """Read-only -- use `set_lines`/`replace_line` to change this plate's lines."""
        return tuple(self._lines)

    def set_lines(self, new_lines: list[ElevationLine]) -> None:
        self._lines = list(new_lines)

    def replace_line(self, index: int, new_line: ElevationLine) -> None:
        self._lines[index] = new_line

    def outline_world(self) -> np.ndarray:
        """Derived directly from each line's current two endpoints -- the actual edge
        boundary evolution maintains (see boundary.py) -- rather than a separately-tracked
        polygon that could drift out of sync with the real data. Traces the high-theta edge
        across lines in ascending phi, then the low-theta edge back down: a standard
        scanline-to-polygon conversion, exact for convex-ish plates and a reasonable
        envelope otherwise. Always non-overlapping with a live-computed neighbor's outline in
        the same sense the underlying elevation data is (see plates.iter_local_lattice /
        boundary.step_boundaries), since it's read from that same data, not duplicated
        state."""
        lines_with_nodes = [line for line in self._lines if len(line) > 0]
        if not lines_with_nodes:
            return np.zeros((0, 3))
        ordered = sorted(lines_with_nodes, key=lambda line: line.phi)
        high_phi = np.array([line.phi for line in ordered])
        high_theta = np.array([line[-1].get_theta() for line in ordered])
        low_theta = np.array([line[0].get_theta() for line in ordered])
        loop_local = np.concatenate(
            [
                geometry.local_xyz(high_phi, high_theta),
                geometry.local_xyz(high_phi[::-1], low_theta[::-1]),
            ],
            axis=0,
        )
        return geometry.to_world(self._frame, loop_local)

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
        return geometry.point_in_spherical_polygon(point_xyz, self.outline_world())

    def get_neighbours(self, all_plates: list["Plate"], threshold_rad: float = NEIGHBOUR_DISTANCE_RAD) -> list["Plate"]:
        return _plates_within(self, all_plates, threshold_rad)

    def __iter__(self) -> Iterator[ElevationPoint]:
        for line in self._lines:
            yield from line


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
        return geometry.point_in_spherical_polygon(point_xyz, self.outline_world())

    def get_neighbours(self, all_plates: list["Plate"], threshold_rad: float = NEIGHBOUR_DISTANCE_RAD) -> list["Plate"]:
        return _plates_within(self, all_plates, threshold_rad)

    def __iter__(self) -> Iterator[ElevationPoint]:
        for i in range(len(self._theta)):
            yield ElevationPointInCloud(self, i)

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


def gather_node_positions(
    plate_list: list[PlateWithLines],
) -> tuple[np.ndarray, list[tuple[PlateWithLines, int, int, int]]]:
    """Every elevation-line node's current world position, concatenated, alongside
    (plate, line_index, start, end) references into the flat array -- the position-only half
    of the near-identical per-step `_gather_nodes` helpers in erosion.py/hydrology.py/
    bathymetry.py, and climate.py's own `_sample_elevation_and_crust`. Factored out here so a
    single step_world call can compute each node's `line.world_xyz(plate.frame)` rotation
    once and pass the same (points, line_refs) into all of them, rather than each
    independently re-deriving identical world positions from plate-local data that hasn't
    moved since the last rotation (see docs/architecture.md's World.climate_cache/
    hydrology_cache notes for the same "compute once this step, reuse" precedent). Each
    caller still gathers its own elevation/other per-node fields fresh off `line_refs` --
    only the rotation itself is shared, since some of those fields (elevation in particular)
    do change mid-step between callers."""
    points_list = []
    line_refs: list[tuple[PlateWithLines, int, int, int]] = []
    offset = 0
    for plate in plate_list:
        for line_index, line in enumerate(plate.lines):
            n = len(line)
            if n == 0:
                continue
            points_list.append(line.world_xyz(plate.frame))
            line_refs.append((plate, line_index, offset, offset + n))
            offset += n
    if not points_list:
        return np.zeros((0, 3)), []
    return np.concatenate(points_list, axis=0), line_refs


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
