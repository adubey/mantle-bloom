"""Plates as spherical polygons carrying their own parallel elevation lines.

Each plate owns a rotation matrix (`frame`) mapping its local (phi, theta) spherical
coordinates to world unit vectors (see `geometry.plate_frame_from_seed`), and a set of
`ElevationLine`s at fixed plate-local latitudes `phi`. Rotating a plate rigidly only ever
updates `frame` -- the (phi, theta) node coordinates themselves never change, so rotation
never needs resampling. See docs/simulation-model.md for the full design writeup.
"""

from __future__ import annotations

import abc
from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from . import ellipse, geometry
from .noise import SphereNoise

PLANET_RADIUS_KM = 6371.0
CONTINENTAL_FRACTION = 0.4
BASE_CONTINENTAL_M = 200.0
BASE_OCEANIC_M = -3800.0
CONTINENTAL_NOISE_AMPLITUDE_M = 1200.0
OCEANIC_NOISE_AMPLITUDE_M = 500.0

# Halving this doubles resolution in each dimension (phi rows and theta samples per row),
# i.e. ~4x the nodes per plate. Several other modules define *absolute node-count*
# thresholds (not distances, which already scale automatically as multiples of
# TARGET_LINE_SPACING_RAD) that represent a physical area or distance in terms of the *old*
# density -- those were rescaled alongside this (merge_split.SPLIT_MIN_NODES,
# gaps.MIN_GAP_POINTS/MAX_ABSORB_NODES_PER_PLATE_PER_CALL by ~4x for area,
# boundary.MAX_EXTEND_NODES_PER_STEP by ~2x for a 1D distance) -- see each for the reasoning.
# This is the reference value for the default node_density=1.0 -- see line_spacing_rad below
# for how a world's own chosen density (World.node_density, set once at generation and read
# by every module in this same list for the rest of that world's life) scales it at runtime,
# now that density is a per-world user choice rather than a hardcoded, one-off code change.
TARGET_LINE_SPACING_KM = 125.0
TARGET_LINE_SPACING_RAD = TARGET_LINE_SPACING_KM / PLANET_RADIUS_KM

# UI-facing choices for World.node_density -- a discrete set (not a free-form slider) since
# there's no natural continuous unit for "how many points," only "how many times as many."
# 2.0 (half the default multiplier) is a lower-resolution middle ground -- fewer nodes than the
# default, so plate-movement-only stepping (World.simulate_plate_movement,
# World.simulate_climate_biomes off) runs faster, without dropping all the way to 1.0's much
# coarser boundary geometry. 0.5 (an eighth of the default) is coarser still -- the fastest,
# lowest-fidelity option, useful where even 1.0's geometry is more than a given step needs.
NODE_DENSITY_CHOICES = (0.5, 1.0, 2.0, 4.0)
DEFAULT_NODE_DENSITY = 4.0


def line_spacing_rad(node_density: float) -> float:
    """The line spacing (radians) that gives a plate ~node_density times as many nodes as
    the default TARGET_LINE_SPACING_RAD would. Node count for a fixed physical area scales
    with the *square* of resolution (see TARGET_LINE_SPACING_KM's own comment -- halving
    spacing quadruples node count), so this divides by sqrt(node_density), not
    node_density itself. Every module that derives a distance threshold or an absolute
    node-count cap from TARGET_LINE_SPACING_RAD calls this (with the world's own
    node_density) instead of reading the bare module constant directly, so that a world
    generated at a non-default density stays self-consistent for its entire life -- not just
    at generation, but through every later regularize/gap-fill/merge/split/volcanism pass
    too (each of those modules' own docstrings/comments explain why its own particular
    thresholds need this)."""
    return TARGET_LINE_SPACING_RAD / np.sqrt(node_density)

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


class ElevationLine:
    """A fixed plate-local latitude `phi` holding elevation (and other persistent, land-only
    or lake/volcano/soil/resource) samples at plate-local longitude nodes `theta`.

    All fields are read-only -- an existing line is never mutated in place. Instead:
    - `replace(...)` swaps in a subset of fields, keeping `theta`'s shape/order untouched
      (the common case -- most steps only ever change `elevation` or one or two persistent
      fields for the same set of nodes).
    - `masked(mask)` filters and/or reorders every field together by a boolean mask or
      fancy index (plate split, node removal, node reassignment reordering).
    - `with_new_nodes(theta, elevation)` appends brand-new nodes (zero/False for every
      persistent field -- no history to carry) to the end, unsorted.

    Threading every persistent field through a single generic method here (keyed off
    `OPTIONAL_FIELDS`, one list) is what earlier avoided a real, previously confirmed bug:
    erosion.py's and bathymetry.py's own hand-written reconstruction sites were both written
    before is_volcano/volcano_active_years_remaining existed, so neither passed them through
    -- silently wiping every node's volcanic status to False every single step. A call site
    that constructs a new field-by-field `ElevationLine(...)` directly (bypassing these
    methods) reintroduces exactly that risk."""

    # All persistent, land-only, meters (unless noted), same shape as theta -- see
    # hydrology.py/lakes.py/volcanism.py/geology.py. Because the grid is plate-local and
    # rotates with a plate's `frame` rather than sitting fixed in world space, these ride
    # along for free just by being an ordinary parallel array on this same line, exactly
    # like elevation itself -- no explicit semi-Lagrangian advection needed every step:
    # rotating a plate only ever touches `frame`, never these arrays.
    OPTIONAL_FIELDS = (
        "channel_depth",  # river channel incision, self-reinforcing
        "channel_width",  # river channel width, grows with flow -- see erosion.py
        "lake_depth",  # standing lake water depth
        "glacier_depth",  # accumulated ice, meters ice-equivalent
        # Sediment settled on a lake's own bed, monotonically increasing (never erodes back
        # away, same self-reinforcing character as channel_depth) -- raises the *effective*
        # floor a lake's own depth is measured against without touching real terrain
        # `elevation` itself, see lakes.py's own module docstring for why. Always 0 outside
        # an active lake.
        "silt_depth",
        # Two more of the same "rides along for free" persistent fields, see volcanism.py.
        # is_volcano never reverts to False once set (permanent provenance -- a dormant
        # volcano is still excluded from being redetected as a fresh rift gap);
        # volcano_active_years_remaining is a countdown, 0 once dormant (whether or not
        # is_volcano is set).
        "is_volcano",  # bool
        "volcano_active_years_remaining",  # years
        # Soil, land-only -- see geology.py. Unlike every other field here, these three can
        # both rise *and* fall (real soil forms and erodes), not just accumulate.
        "soil_depth",  # meters, regolith/soil thickness
        "soil_mineral_content",  # [0, 1], weathered/hydrothermal richness
        "soil_organic_content",  # [0, 1], accumulated organic matter
        # Resource deposits -- see geology.py/volcanism.py. All monotonically non-decreasing,
        # the same self-reinforcing "once formed, never erodes back away" convention
        # silt_depth already uses (buried peat/hydrocarbons/ore aren't un-buried by a later
        # climate shift).
        "coal_deposit_m",  # land-only
        "oil_gas_deposit_m",  # ocean-only
        "mineral_deposit_m",  # either -- grown by volcanism.py's own eruptions
    )

    def __init__(
        self,
        phi: float,
        theta: np.ndarray,
        elevation: np.ndarray,
        channel_depth: np.ndarray | None = None,
        channel_width: np.ndarray | None = None,
        lake_depth: np.ndarray | None = None,
        glacier_depth: np.ndarray | None = None,
        silt_depth: np.ndarray | None = None,
        is_volcano: np.ndarray | None = None,
        volcano_active_years_remaining: np.ndarray | None = None,
        soil_depth: np.ndarray | None = None,
        soil_mineral_content: np.ndarray | None = None,
        soil_organic_content: np.ndarray | None = None,
        coal_deposit_m: np.ndarray | None = None,
        oil_gas_deposit_m: np.ndarray | None = None,
        mineral_deposit_m: np.ndarray | None = None,
    ) -> None:
        self._phi = phi
        self._theta = theta
        self._elevation = elevation
        self._channel_depth = channel_depth if channel_depth is not None else np.zeros_like(theta)
        self._channel_width = channel_width if channel_width is not None else np.zeros_like(theta)
        self._lake_depth = lake_depth if lake_depth is not None else np.zeros_like(theta)
        self._glacier_depth = glacier_depth if glacier_depth is not None else np.zeros_like(theta)
        self._silt_depth = silt_depth if silt_depth is not None else np.zeros_like(theta)
        self._is_volcano = is_volcano if is_volcano is not None else np.zeros_like(theta, dtype=bool)
        self._volcano_active_years_remaining = (
            volcano_active_years_remaining if volcano_active_years_remaining is not None else np.zeros_like(theta)
        )
        self._soil_depth = soil_depth if soil_depth is not None else np.zeros_like(theta)
        self._soil_mineral_content = soil_mineral_content if soil_mineral_content is not None else np.zeros_like(theta)
        self._soil_organic_content = soil_organic_content if soil_organic_content is not None else np.zeros_like(theta)
        self._coal_deposit_m = coal_deposit_m if coal_deposit_m is not None else np.zeros_like(theta)
        self._oil_gas_deposit_m = oil_gas_deposit_m if oil_gas_deposit_m is not None else np.zeros_like(theta)
        self._mineral_deposit_m = mineral_deposit_m if mineral_deposit_m is not None else np.zeros_like(theta)

    @property
    def phi(self) -> float:
        return self._phi

    @property
    def theta(self) -> np.ndarray:
        return self._theta

    @property
    def elevation(self) -> np.ndarray:
        return self._elevation

    @property
    def channel_depth(self) -> np.ndarray:
        return self._channel_depth

    @property
    def channel_width(self) -> np.ndarray:
        return self._channel_width

    @property
    def lake_depth(self) -> np.ndarray:
        return self._lake_depth

    @property
    def glacier_depth(self) -> np.ndarray:
        return self._glacier_depth

    @property
    def silt_depth(self) -> np.ndarray:
        return self._silt_depth

    @property
    def is_volcano(self) -> np.ndarray:
        return self._is_volcano

    @property
    def volcano_active_years_remaining(self) -> np.ndarray:
        return self._volcano_active_years_remaining

    @property
    def soil_depth(self) -> np.ndarray:
        return self._soil_depth

    @property
    def soil_mineral_content(self) -> np.ndarray:
        return self._soil_mineral_content

    @property
    def soil_organic_content(self) -> np.ndarray:
        return self._soil_organic_content

    @property
    def coal_deposit_m(self) -> np.ndarray:
        return self._coal_deposit_m

    @property
    def oil_gas_deposit_m(self) -> np.ndarray:
        return self._oil_gas_deposit_m

    @property
    def mineral_deposit_m(self) -> np.ndarray:
        return self._mineral_deposit_m

    def world_xyz(self, frame: np.ndarray) -> np.ndarray:
        phi_arr = np.full_like(self.theta, self.phi)
        local = geometry.local_xyz(phi_arr, self.theta)
        return geometry.to_world(frame, local)

    def replace(self, **overrides: np.ndarray) -> "ElevationLine":
        """A new line with the given fields (elevation and/or any of OPTIONAL_FIELDS)
        swapped in and every other field copied from this one unchanged -- `theta`/`phi`
        are never touched here, so only use this when the node set itself isn't changing."""
        kwargs: dict[str, np.ndarray] = {name: getattr(self, name) for name in self.OPTIONAL_FIELDS}
        kwargs["elevation"] = self.elevation
        kwargs.update(overrides)
        return ElevationLine(phi=self.phi, theta=self.theta, **kwargs)

    def masked(self, mask) -> "ElevationLine":
        """A new line with `theta`, `elevation`, and every OPTIONAL_FIELDS array filtered
        and/or reordered together by a boolean mask or fancy index -- for removing nodes
        (plate split, node reassignment) or reordering them (after concatenating in new
        nodes at the end)."""
        kwargs = {name: getattr(self, name)[mask] for name in self.OPTIONAL_FIELDS}
        return ElevationLine(phi=self.phi, theta=self.theta[mask], elevation=self.elevation[mask], **kwargs)

    def with_new_nodes(self, theta: np.ndarray, elevation: np.ndarray) -> "ElevationLine":
        """A new line with `theta`/`elevation` nodes appended at the end -- every
        OPTIONAL_FIELDS value for the new nodes starts at zero/False, no history to carry.
        The result is unsorted by theta; follow with `.masked(np.argsort(new_line.theta))`
        if ascending order matters to the caller."""
        n = len(theta)
        kwargs = {
            name: np.concatenate([getattr(self, name), np.zeros(n, dtype=getattr(self, name).dtype)])
            for name in self.OPTIONAL_FIELDS
        }
        return ElevationLine(
            phi=self.phi,
            theta=np.concatenate([self.theta, theta]),
            elevation=np.concatenate([self.elevation, elevation]),
            **kwargs,
        )


class Plate(abc.ABC):
    """A plate's shared identity/motion state plus an abstract interface over however it
    represents its own terrain nodes -- currently only `PlateWithLines` (parallel
    `ElevationLine`s), with an R-tree-backed representation planned as a future second
    subclass. Every method below that doesn't depend on node representation (motion,
    identity) is implemented once here; `node_count`/`all_points_and_elevation`/
    `outline_world`/`collect` are representation-specific and left abstract."""

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
        lines_with_nodes = [line for line in self._lines if len(line.theta) > 0]
        if not lines_with_nodes:
            return np.zeros((0, 3))
        ordered = sorted(lines_with_nodes, key=lambda line: line.phi)
        high_phi = np.array([line.phi for line in ordered])
        high_theta = np.array([line.theta[-1] for line in ordered])
        low_theta = np.array([line.theta[0] for line in ordered])
        loop_local = np.concatenate(
            [
                geometry.local_xyz(high_phi, high_theta),
                geometry.local_xyz(high_phi[::-1], low_theta[::-1]),
            ],
            axis=0,
        )
        return geometry.to_world(self._frame, loop_local)

    def node_count(self) -> int:
        return sum(len(line.theta) for line in self._lines)

    def all_points_and_elevation(self) -> tuple[np.ndarray, np.ndarray]:
        """Every elevation-line node's world position and elevation, concatenated."""
        if not self._lines:
            return np.zeros((0, 3)), np.zeros(0)
        points = np.concatenate([line.world_xyz(self._frame) for line in self._lines], axis=0)
        elevation = np.concatenate([line.elevation for line in self._lines], axis=0)
        return points, elevation

    def collect(self, field_name: str) -> np.ndarray:
        chunks = [getattr(line, field_name) for line in self._lines if len(line.theta) > 0]
        if not chunks:
            return np.zeros(0, dtype=bool) if field_name == "is_volcano" else np.zeros(0)
        return np.concatenate(chunks, axis=0)


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
            n = len(line.theta)
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


def iter_local_lattice(frame: np.ndarray, spacing_rad: float = TARGET_LINE_SPACING_RAD):
    """Sweep a full plate-local (phi, theta) lattice at `spacing_rad` resolution, yielding
    (phi, theta_candidates, world_pts) per row. Shared by initial generation and by
    plate-merge resampling (see merge_split.py), and, at a resolution independent of the
    physical line spacing, by the render-grid sweep (see render_image.py's
    _render_grid_arrays) that gives the rendered map full coverage regardless of how sparse
    the underlying physical data is once projected."""
    max_abs_phi = np.pi / 2 - spacing_rad / 2
    phi_values = np.arange(-max_abs_phi, max_abs_phi, spacing_rad)
    for phi in phi_values:
        dtheta = spacing_rad / max(np.cos(phi), 1e-3)
        n_theta = max(int(np.round(2 * np.pi / dtheta)), 1)
        theta_candidates = np.linspace(-np.pi, np.pi, n_theta, endpoint=False)

        local_pts = geometry.local_xyz(np.full_like(theta_candidates, phi), theta_candidates)
        world_pts = geometry.to_world(frame, local_pts)
        yield float(phi), theta_candidates, world_pts


def build_lines_from_lattice(frame: np.ndarray, is_owned, elevation_at, spacing_rad: float = TARGET_LINE_SPACING_RAD) -> list[ElevationLine]:
    """Build a plate's elevation lines by sweeping its local lattice and keeping whichever
    nodes `is_owned(world_pts) -> bool array` selects, with elevation from
    `elevation_at(owned_world_pts) -> array`. `spacing_rad` defaults to the reference
    density (1.0) -- every caller that has a `World` in hand should instead pass
    `line_spacing_rad(world.node_density)`, so newly-built lines (initial generation, gap
    absorption/spawning, plate merges, volcanic fields) match whatever density that world was
    actually generated at, not silently fall back to the default."""
    lines: list[ElevationLine] = []
    for phi, theta_candidates, world_pts in iter_local_lattice(frame, spacing_rad=spacing_rad):
        owned = is_owned(world_pts)
        if not np.any(owned):
            continue
        theta_owned = theta_candidates[owned]
        elevation = elevation_at(world_pts[owned])
        lines.append(ElevationLine(phi=phi, theta=theta_owned, elevation=elevation))
    return lines


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
