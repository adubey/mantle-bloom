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
from scipy.spatial import SphericalVoronoi, cKDTree

from . import geometry
from .noise import SphereNoise

PLANET_RADIUS_KM = 6371.0
CONTINENTAL_FRACTION = 0.4
BASE_CONTINENTAL_M = 200.0
BASE_OCEANIC_M = -3800.0
CONTINENTAL_NOISE_AMPLITUDE_M = 1200.0
OCEANIC_NOISE_AMPLITUDE_M = 500.0

TARGET_LINE_SPACING_KM = 250.0
TARGET_LINE_SPACING_RAD = TARGET_LINE_SPACING_KM / PLANET_RADIUS_KM

# How close to a plate-local pole (phi = +-pi/2) a line can get before its circumference
# is too small to bother sampling.
_MAX_ABS_PHI = np.pi / 2 - TARGET_LINE_SPACING_RAD / 2


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
    boundary_local: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))  # local xyz
    lines: list[ElevationLine] = field(default_factory=list)
    # Steps since this plate was created (by generation, merge, or split). Gates split
    # eligibility in merge_split.py so a plate can't fragment repeatedly in quick
    # succession -- see the note there on why that runaway is a real failure mode.
    age_steps: int = 0

    @property
    def seed_world(self) -> np.ndarray:
        """World position of this plate's local (phi=0, theta=0) reference point."""
        return self.frame[:, 0]

    def boundary_world(self) -> np.ndarray:
        return geometry.to_world(self.frame, self.boundary_local)

    def node_count(self) -> int:
        return sum(len(line.theta) for line in self.lines)

    def all_points_and_elevation(self) -> tuple[np.ndarray, np.ndarray]:
        """Every elevation-line node's world position and elevation, concatenated."""
        if not self.lines:
            return np.zeros((0, 3)), np.zeros(0)
        points = np.concatenate([line.world_xyz(self.frame) for line in self.lines], axis=0)
        elevation = np.concatenate([line.elevation for line in self.lines], axis=0)
        return points, elevation


def _base_elevation(crust_type: str) -> float:
    return BASE_CONTINENTAL_M if crust_type == "continental" else BASE_OCEANIC_M


def _noise_amplitude(crust_type: str) -> float:
    return CONTINENTAL_NOISE_AMPLITUDE_M if crust_type == "continental" else OCEANIC_NOISE_AMPLITUDE_M


def iter_local_lattice(frame: np.ndarray):
    """Sweep a full plate-local (phi, theta) lattice at TARGET_LINE_SPACING_RAD resolution,
    yielding (phi, theta_candidates, world_pts) per row. Shared by initial generation and
    by plate-merge resampling (see merge_split.py) -- the two places a full-footprint sweep
    is needed."""
    phi_values = np.arange(-_MAX_ABS_PHI, _MAX_ABS_PHI, TARGET_LINE_SPACING_RAD)
    for phi in phi_values:
        dtheta = TARGET_LINE_SPACING_RAD / max(np.cos(phi), 1e-3)
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
) -> list[ElevationLine]:
    """Keep only lattice nodes whose nearest seed is this plate's own seed (i.e. nodes
    actually inside this plate's spherical Voronoi cell), and assign each a base elevation
    plus noise texture."""
    base = _base_elevation(crust_type)
    amp = _noise_amplitude(crust_type)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        _, nearest_idx = owner_tree.query(world_pts)
        return nearest_idx == plate_index

    def elevation_at(world_pts: np.ndarray) -> np.ndarray:
        return base + amp * noise.sample(world_pts)

    return build_lines_from_lattice(frame, is_owned, elevation_at)


def generate_plates(seed: int, num_plates: int = 12) -> list[Plate]:
    rng = np.random.default_rng(seed)

    seed_xyz = rng.normal(size=(num_plates, 3))
    seed_xyz /= np.linalg.norm(seed_xyz, axis=-1, keepdims=True)

    crust_types = [
        "continental" if rng.random() < CONTINENTAL_FRACTION else "oceanic" for _ in range(num_plates)
    ]

    sv = SphericalVoronoi(seed_xyz, radius=1.0)
    sv.sort_vertices_of_regions()

    owner_tree = cKDTree(seed_xyz)
    noise = SphereNoise(rng, octaves=4, base_freq=2.5)

    plates: list[Plate] = []
    for i in range(num_plates):
        frame = geometry.plate_frame_from_seed(seed_xyz[i])
        boundary_world = sv.vertices[sv.regions[i]]
        boundary_local = geometry.to_local(frame, boundary_world)
        lines = _build_lines_for_plate(i, frame, crust_types[i], owner_tree, noise)
        plates.append(
            Plate(
                plate_id=i,
                frame=frame,
                crust_type=crust_types[i],
                boundary_local=boundary_local,
                lines=lines,
            )
        )
    return plates
