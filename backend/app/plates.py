"""Plates as spherical polygons carrying their own parallel elevation lines.

Each plate owns a rotation matrix (`frame`) mapping its local (phi, theta) spherical
coordinates to world unit vectors (see `geometry.plate_frame_from_seed`), and a set of
`ElevationLine`s at fixed plate-local latitudes `phi`. Rotating a plate rigidly only ever
updates `frame` -- the (phi, theta) node coordinates themselves never change, so rotation
never needs resampling. See docs/simulation-model.md for the full design writeup.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from . import geometry
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
TARGET_LINE_SPACING_KM = 125.0
TARGET_LINE_SPACING_RAD = TARGET_LINE_SPACING_KM / PLANET_RADIUS_KM

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


@dataclass
class ElevationLine:
    phi: float  # plate-local latitude, radians, constant along the line
    theta: np.ndarray  # plate-local longitudes of nodes, radians, ascending
    elevation: np.ndarray  # meters, same shape as theta

    def world_xyz(self, frame: np.ndarray) -> np.ndarray:
        phi_arr = np.full_like(self.theta, self.phi)
        local = geometry.local_xyz(phi_arr, self.theta)
        return geometry.to_world(frame, local)


@dataclass
class Plate:
    plate_id: int
    frame: np.ndarray  # 3x3 rotation matrix, local -> world
    crust_type: str  # "continental" or "oceanic"
    omega: np.ndarray = field(default_factory=lambda: np.zeros(3))  # angular velocity, world frame
    lines: list[ElevationLine] = field(default_factory=list)
    # Steps since this plate was created (by generation, merge, or split). Gates split
    # eligibility in merge_split.py so a plate can't fragment repeatedly in quick
    # succession -- see the note there on why that runaway is a real failure mode.
    age_steps: int = 0

    @property
    def seed_world(self) -> np.ndarray:
        """World position of this plate's local (phi=0, theta=0) reference point."""
        return self.frame[:, 0]

    def outline_world(self) -> np.ndarray:
        """A live approximation of this plate's current territory outline, derived
        directly from each line's current two endpoints -- the actual edge boundary
        evolution maintains (see boundary.py) -- rather than a separately-tracked polygon
        that could drift out of sync with the real data. Traces the high-theta edge across
        lines in ascending phi, then the low-theta edge back down: a standard scanline-to-
        polygon conversion, exact for convex-ish plates and a reasonable envelope
        otherwise. Always non-overlapping with a live-computed neighbor's outline in the
        same sense the underlying elevation data is (see plates.iter_local_lattice /
        boundary.step_boundaries), since it's read from that same data, not duplicated
        state."""
        lines_with_nodes = [line for line in self.lines if len(line.theta) > 0]
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
        return geometry.to_world(self.frame, loop_local)

    def node_count(self) -> int:
        return sum(len(line.theta) for line in self.lines)

    def all_points_and_elevation(self) -> tuple[np.ndarray, np.ndarray]:
        """Every elevation-line node's world position and elevation, concatenated."""
        if not self.lines:
            return np.zeros((0, 3)), np.zeros(0)
        points = np.concatenate([line.world_xyz(self.frame) for line in self.lines], axis=0)
        elevation = np.concatenate([line.elevation for line in self.lines], axis=0)
        return points, elevation


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


def build_lines_from_lattice(frame: np.ndarray, is_owned, elevation_at) -> list[ElevationLine]:
    """Build a plate's elevation lines by sweeping its local lattice and keeping whichever
    nodes `is_owned(world_pts) -> bool array` selects, with elevation from
    `elevation_at(owned_world_pts) -> array`."""
    lines: list[ElevationLine] = []
    for phi, theta_candidates, world_pts in iter_local_lattice(frame):
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

    return build_lines_from_lattice(frame, is_owned, elevation_at)


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
    _land_noise_threshold for how that target is actually hit.

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

    plates: list[Plate] = []
    for i in range(num_plates):
        frame = geometry.plate_frame_from_seed(seed_xyz[i])
        lines = _build_lines_for_plate(i, frame, crust_types[i], owner_tree, noise, land_threshold)
        plates.append(Plate(plate_id=i, frame=frame, crust_type=crust_types[i], lines=lines))
    return plates
