"""A simple mantle-convection flow model, plus fitting each plate's rigid rotation (Euler
pole + rate) to that flow field.

Convection cell centers are placed via a cubed-sphere mapping (gnomonic projection of a
cube face, normalized to the unit sphere) purely for even coverage with no pole clustering
-- "flow in a cube, projected to the sphere". Each center pushes (upwelling) or pulls
(downwelling) points near it tangentially, with distance falloff; plate velocity is never
set directly, only ever fit to this field.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PLANET_RADIUS_KM = 6371.0
KM_PER_CM = 1e-5


def cm_per_yr_to_rad_per_yr(cm_per_yr: float) -> float:
    return (cm_per_yr * KM_PER_CM) / PLANET_RADIUS_KM


def rad_per_yr_to_cm_per_yr(rad_per_yr: float) -> float:
    """Inverse of `cm_per_yr_to_rad_per_yr` -- a plate's `|omega|` as a surface speed at the
    planet's radius, for human-readable diagnostics (the Plate Inspector, /world/plates)."""
    return (rad_per_yr * PLANET_RADIUS_KM) / KM_PER_CM


MANTLE_FLOW_REFERENCE_RATE = cm_per_yr_to_rad_per_yr(4.0)
MIN_PLATE_RATE = cm_per_yr_to_rad_per_yr(0.5)
MAX_PLATE_RATE = cm_per_yr_to_rad_per_yr(15.0)
VELOCITY_DAMPING = 0.3  # fraction of the way from current omega toward the freshly-fit target

_CUBE_FACES = ["+x", "-x", "+y", "-y", "+z", "-z"]


def cube_to_sphere(face: str, u: float, v: float) -> np.ndarray:
    """Gnomonic cubed-sphere mapping: a point on a cube face (u, v in [-1, 1]) -> a unit
    vector on the sphere."""
    if face == "+x":
        p = np.array([1.0, u, v])
    elif face == "-x":
        p = np.array([-1.0, u, v])
    elif face == "+y":
        p = np.array([u, 1.0, v])
    elif face == "-y":
        p = np.array([u, -1.0, v])
    elif face == "+z":
        p = np.array([u, v, 1.0])
    elif face == "-z":
        p = np.array([u, v, -1.0])
    else:
        raise ValueError(f"unknown cube face {face!r}")
    return p / np.linalg.norm(p)


@dataclass
class ConvectionCenter:
    position: np.ndarray  # unit vector, world frame
    strength: float  # signed; positive = upwelling (pushes away), negative = downwelling (pulls in)
    falloff: float  # radians; controls how quickly influence decays with angular distance


def generate_convection_centers(
    rng: np.random.Generator,
    n_centers: int = 8,
    reference_rate: float = MANTLE_FLOW_REFERENCE_RATE,
) -> list[ConvectionCenter]:
    centers = []
    for _ in range(n_centers):
        face = _CUBE_FACES[rng.integers(0, len(_CUBE_FACES))]
        u, v = rng.uniform(-1.0, 1.0, size=2)
        position = cube_to_sphere(face, u, v)
        sign = 1.0 if rng.random() < 0.5 else -1.0
        strength = sign * reference_rate * rng.uniform(0.5, 1.5)
        falloff = rng.uniform(0.6, 1.4)  # radians
        centers.append(ConvectionCenter(position=position, strength=strength, falloff=falloff))
    return centers


def flow_at(xyz: np.ndarray, centers: list[ConvectionCenter]) -> np.ndarray:
    """Tangential mantle-flow velocity at each point in `xyz` (..., 3), summed over all
    convection centers. Units are radians/year-equivalent (an angular rate, since these
    points are unit vectors) so they can be fit directly to an Euler pole/rate."""
    xyz = np.asarray(xyz, dtype=float)
    total = np.zeros_like(xyz)
    for c in centers:
        dot = np.clip(np.sum(xyz * c.position, axis=-1, keepdims=True), -1.0, 1.0)
        angular_dist = np.arccos(dot)

        # Tangential direction at each point, pointing away from the center.
        away = dot * xyz - c.position
        away_norm = np.linalg.norm(away, axis=-1, keepdims=True)
        safe_norm = np.where(away_norm < 1e-9, 1.0, away_norm)
        direction = away / safe_norm

        falloff = np.exp(-((angular_dist / c.falloff) ** 2))
        magnitude = c.strength * falloff
        contribution = direction * magnitude
        contribution = np.where(away_norm < 1e-9, 0.0, contribution)
        total = total + contribution
    return total


def _skew_batch_neg(points: np.ndarray) -> np.ndarray:
    """For each point p (N, 3), build -skew(p) (N, 3, 3) so that (-skew(p)) @ omega ==
    omega x p -- i.e. the linear operator mapping omega to the rigid-rotation velocity
    it induces at p."""
    n = points.shape[0]
    a = np.zeros((n, 3, 3))
    px, py, pz = points[:, 0], points[:, 1], points[:, 2]
    a[:, 0, 1] = pz
    a[:, 0, 2] = -py
    a[:, 1, 0] = -pz
    a[:, 1, 2] = px
    a[:, 2, 0] = py
    a[:, 2, 1] = -px
    return a


def fit_euler_pole(points: np.ndarray, velocities: np.ndarray) -> np.ndarray:
    """Least-squares best-fit angular velocity omega such that omega x p_i ~= v_i for each
    sample point/velocity pair. Returns the omega vector (axis = pole, magnitude = rate)."""
    points = np.asarray(points, dtype=float)
    velocities = np.asarray(velocities, dtype=float)
    a = _skew_batch_neg(points).reshape(-1, 3)
    b = velocities.reshape(-1)
    omega, *_ = np.linalg.lstsq(a, b, rcond=None)
    return omega


def clamp_rate(omega: np.ndarray, min_rate: float = MIN_PLATE_RATE, max_rate: float = MAX_PLATE_RATE) -> np.ndarray:
    rate = np.linalg.norm(omega)
    if rate < 1e-15:
        return omega
    clamped_rate = np.clip(rate, min_rate, max_rate)
    return omega * (clamped_rate / rate)
