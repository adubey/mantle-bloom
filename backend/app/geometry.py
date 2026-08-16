"""Sphere vector math: unit-vector <-> lat/lon conversion, axis-angle rotation
(Rodrigues' formula), and per-plate local coordinate frames.

A plate-local frame is a 3x3 rotation matrix `R_plate` whose columns are the world-space
images of the plate-local basis vectors. Plate-local spherical coordinates (phi, theta) use
the standard convention x=cos(phi)cos(theta), y=cos(phi)sin(theta), z=sin(phi), so
local (1,0,0) is the plate's (phi=0, theta=0) reference point and local (0,0,1) is the
plate's local "pole". World position of any plate-local point is `R_plate @ local_xyz(phi, theta)`.
"""

from __future__ import annotations

import numpy as np

_EPS = 1e-12


def latlon_to_xyz(lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
    """Convert geographic lat/lon (radians) to unit vectors. Returns shape (..., 3)."""
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    x = np.cos(lat) * np.cos(lon)
    y = np.cos(lat) * np.sin(lon)
    z = np.sin(lat)
    return np.stack([x, y, z], axis=-1)


def xyz_to_latlon(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert unit vectors (..., 3) to (lat, lon) in radians."""
    xyz = np.asarray(xyz, dtype=float)
    x, y, z = xyz[..., 0], xyz[..., 1], xyz[..., 2]
    lat = np.arcsin(np.clip(z, -1.0, 1.0))
    lon = np.arctan2(y, x)
    return lat, lon


def local_xyz(phi: np.ndarray, theta: np.ndarray) -> np.ndarray:
    """Plate-local spherical coordinates -> plate-local unit vectors, shape (..., 3)."""
    return latlon_to_xyz(phi, theta)


def normalize(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    norm = np.linalg.norm(v, axis=-1, keepdims=True)
    norm = np.where(norm < _EPS, 1.0, norm)
    return v / norm


def skew(axis: np.ndarray) -> np.ndarray:
    """3x3 skew-symmetric cross-product matrix for a single 3-vector."""
    ax, ay, az = axis
    return np.array(
        [
            [0.0, -az, ay],
            [az, 0.0, -ax],
            [-ay, ax, 0.0],
        ]
    )


def rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rodrigues' rotation formula: 3x3 matrix rotating by `angle` radians about `axis`."""
    axis = np.asarray(axis, dtype=float)
    n = np.linalg.norm(axis)
    if n < _EPS or abs(angle) < _EPS:
        return np.eye(3)
    axis = axis / n
    k = skew(axis)
    return np.eye(3) + np.sin(angle) * k + (1.0 - np.cos(angle)) * (k @ k)


def rotation_matrix_from_omega(omega: np.ndarray, dt: float) -> np.ndarray:
    """Angular velocity vector (axis direction = pole, magnitude = rate) -> rotation matrix
    for advancing by time `dt` under that constant angular velocity."""
    omega = np.asarray(omega, dtype=float)
    rate = np.linalg.norm(omega)
    if rate < _EPS or abs(dt) < _EPS:
        return np.eye(3)
    axis = omega / rate
    return rotation_matrix(axis, rate * dt)


def rotate_vectors(vectors: np.ndarray, axis: np.ndarray, angle: float) -> np.ndarray:
    """Rotate an array of vectors (..., 3) by `angle` radians about `axis`."""
    r = rotation_matrix(axis, angle)
    return vectors @ r.T


def plate_frame_from_seed(seed_xyz: np.ndarray, up_hint: np.ndarray | None = None) -> np.ndarray:
    """Build a plate-local frame (3x3 rotation matrix, world = R @ local) whose local
    (phi=0, theta=0) point is `seed_xyz`.

    `up_hint` picks the local pole direction (projected into the tangent plane at the seed);
    defaults to the global +z axis, falling back to +x if the seed is too close to a global pole
    for that to be well-defined.
    """
    s = normalize(np.asarray(seed_xyz, dtype=float))
    if up_hint is None:
        up_hint = np.array([0.0, 0.0, 1.0])
    up_hint = np.asarray(up_hint, dtype=float)

    n = up_hint - np.dot(up_hint, s) * s
    if np.linalg.norm(n) < 1e-6:
        up_hint = np.array([1.0, 0.0, 0.0])
        n = up_hint - np.dot(up_hint, s) * s
    n = normalize(n)

    e2 = normalize(np.cross(n, s))
    return np.stack([s, e2, n], axis=-1)


def angular_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Great-circle angular distance (radians) between unit vectors a and b."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    dot = np.clip(np.sum(a * b, axis=-1), -1.0, 1.0)
    return np.arccos(dot)


def to_local(frame: np.ndarray, world_xyz: np.ndarray) -> np.ndarray:
    """World unit vectors (..., 3) -> plate-local unit vectors, given local frame `frame`."""
    return world_xyz @ frame  # frame.T @ v for each v, vectorized as v @ frame


def to_world(frame: np.ndarray, local_xyz_vec: np.ndarray) -> np.ndarray:
    """Plate-local unit vectors (..., 3) -> world unit vectors, given local frame `frame`."""
    return local_xyz_vec @ frame.T
