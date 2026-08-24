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


def bounding_sphere(points: np.ndarray) -> tuple[np.ndarray, float]:
    """A point cloud's centroid direction and the max angular distance from it to any of
    its own points. Cheap enough to compute per plate per step, and -- via the triangle
    inequality on the sphere -- enough to lower-bound how close two point clouds' *closest*
    points could possibly be without comparing every point in either cloud: if
    `angular_distance(centroid_a, centroid_b) - radius_a - radius_b` exceeds some distance
    threshold, no point in cloud a can be within that threshold of any point in cloud b."""
    centroid = normalize(points.mean(axis=0))
    return centroid, float(angular_distance(points, centroid).max())


def local_tangent_basis(pole: np.ndarray, up_hint: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """(east, north) unit tangent vectors at `pole`. Same up-hint-projected-into-tangent-
    plane-with-near-pole-fallback construction as `plate_frame_from_seed`, but `pole` here
    *is* the frame's actual pole (unlike `plate_frame_from_seed`, where the seed is the
    equatorial phi=0/theta=0 reference point) -- what
    `azimuthal_equidistant_forward`/`_inverse` below need."""
    pole = normalize(np.asarray(pole, dtype=float))
    up = up_hint if up_hint is not None else np.array([0.0, 0.0, 1.0])
    up = np.asarray(up, dtype=float)
    n = up - np.dot(up, pole) * pole
    if np.linalg.norm(n) < 1e-6:
        up = np.array([1.0, 0.0, 0.0])
        n = up - np.dot(up, pole) * pole
    north = normalize(n)
    east = normalize(np.cross(north, pole))
    return east, north


def azimuthal_equidistant_forward(pole: np.ndarray, east: np.ndarray, north: np.ndarray, points: np.ndarray) -> np.ndarray:
    """Unit vectors (..., 3) -> planar (x, y) in radians at unit-sphere radius (multiply by
    a physical radius, e.g. plates.PLANET_RADIUS_KM, for real distances). Distance from the
    origin is *exact* great-circle distance from `pole` -- this is not distance-preserving
    between two arbitrary non-pole points, only from the shared center, which is all
    `plates.plate_bounding_ellipse` needs (everything is fit relative to one center).

    Numerically singular for a point antipodal (or extremely close to antipodal) to `pole`
    (the tangent direction has ~zero length there) -- such a point silently maps to the
    origin instead of its true (theta=pi, arbitrary bearing). Not a concern for the compact,
    Voronoi-seeded plates this is actually used on, but a real limitation for an artificially
    non-convex or hemisphere-spanning point cloud; not solved here (see
    plates.plate_bounding_ellipse's docstring)."""
    pole = normalize(pole)
    cos_c = np.clip(np.sum(points * pole, axis=-1), -1.0, 1.0)
    theta = np.arccos(cos_c)
    tangent = points - cos_c[..., None] * pole
    norm = np.linalg.norm(tangent, axis=-1, keepdims=True)
    safe = np.where(norm < 1e-12, 1.0, norm)
    tangent_dir = tangent / safe
    # These are already cos(bearing)/sin(bearing) -- components of a unit vector along
    # east/north -- so no separate atan2/trig round-trip is needed to get planar coordinates.
    east_comp = np.sum(tangent_dir * east, axis=-1)
    north_comp = np.sum(tangent_dir * north, axis=-1)
    return np.stack([theta * east_comp, theta * north_comp], axis=-1)


def azimuthal_equidistant_inverse(pole: np.ndarray, east: np.ndarray, north: np.ndarray, xy: np.ndarray) -> np.ndarray:
    """Inverse of `azimuthal_equidistant_forward`: planar (x, y) radians -> unit vectors."""
    pole = normalize(pole)
    x, y = xy[..., 0], xy[..., 1]
    theta = np.hypot(x, y)
    safe = np.where(theta < 1e-12, 1.0, theta)
    tangent_dir = (x / safe)[..., None] * east + (y / safe)[..., None] * north
    return np.cos(theta)[..., None] * pole + np.sin(theta)[..., None] * tangent_dir


def point_in_spherical_polygon(point_xyz: np.ndarray, polygon_xyz: np.ndarray) -> bool:
    """True if unit vector `point_xyz` lies inside the simple spherical polygon whose
    ordered (CW or CCW, either works) vertices are `polygon_xyz` (unit vectors) -- the
    spherical analogue of the planar winding-number point-in-polygon test. Projects each
    vertex into `point_xyz`'s own local tangent plane (see `local_tangent_basis`) and sums
    the signed bearing change edge to edge: a point enclosed by a simple polygon
    accumulates a full +-2*pi turn as the vertices sweep around it, one outside accumulates
    ~0. Degenerate (numerically singular) for a polygon vertex antipodal to `point_xyz`,
    not a concern for the compact plate outlines this is actually used on."""
    polygon_xyz = np.asarray(polygon_xyz, dtype=float)
    if len(polygon_xyz) < 3:
        return False
    east, north = local_tangent_basis(point_xyz)
    tangent = polygon_xyz - np.outer(polygon_xyz @ point_xyz, point_xyz)
    bearings = np.arctan2(tangent @ north, tangent @ east)
    diffs = np.diff(np.concatenate([bearings, bearings[:1]]))
    diffs = (diffs + np.pi) % (2.0 * np.pi) - np.pi
    return abs(float(np.sum(diffs))) > np.pi


def points_in_spherical_polygon(points_xyz: np.ndarray, polygon_xyz: np.ndarray) -> np.ndarray:
    """Vectorized `point_in_spherical_polygon`: same winding-number algorithm, batched over
    every point in `points_xyz` at once instead of a Python-level loop calling the scalar
    version per point. Needed once `PlateWithLines.deform` started calling the containment
    test for every one of a plate's own near-boundary nodes, every turn -- profiled directly
    as the dominant per-step cost at realistic node counts (a single step_world call on a
    10-plate, default-density world went from ~46s to well under a second after switching
    to this). Returns a bool array, one entry per point in `points_xyz` (empty if either
    input is empty or the polygon has fewer than 3 vertices)."""
    points_xyz = np.asarray(points_xyz, dtype=float)
    polygon_xyz = np.asarray(polygon_xyz, dtype=float)
    n = len(points_xyz)
    if len(polygon_xyz) < 3 or n == 0:
        return np.zeros(n, dtype=bool)

    up = np.array([0.0, 0.0, 1.0])
    north = up[None, :] - (points_xyz @ up)[:, None] * points_xyz
    norms = np.linalg.norm(north, axis=-1)
    degenerate = norms < 1e-6
    if np.any(degenerate):
        alt_up = np.array([1.0, 0.0, 0.0])
        north[degenerate] = alt_up[None, :] - (points_xyz[degenerate] @ alt_up)[:, None] * points_xyz[degenerate]
        norms = np.linalg.norm(north, axis=-1)
    north = north / norms[:, None]
    east = np.cross(north, points_xyz)
    east = east / np.linalg.norm(east, axis=-1, keepdims=True)

    # tangent: (n_points, n_vertices, 3) -- every polygon vertex projected into each query
    # point's own local tangent plane, all at once.
    dot_vp = polygon_xyz @ points_xyz.T  # (n_vertices, n_points)
    tangent = polygon_xyz[None, :, :] - dot_vp.T[:, :, None] * points_xyz[:, None, :]
    east_comp = np.einsum("nmc,nc->nm", tangent, east)
    north_comp = np.einsum("nmc,nc->nm", tangent, north)
    bearings = np.arctan2(north_comp, east_comp)  # (n_points, n_vertices)
    closed = np.concatenate([bearings, bearings[:, :1]], axis=1)
    diffs = np.diff(closed, axis=1)
    diffs = (diffs + np.pi) % (2.0 * np.pi) - np.pi
    return np.abs(np.sum(diffs, axis=1)) > np.pi


def to_local(frame: np.ndarray, world_xyz: np.ndarray) -> np.ndarray:
    """World unit vectors (..., 3) -> plate-local unit vectors, given local frame `frame`."""
    return world_xyz @ frame  # frame.T @ v for each v, vectorized as v @ frame


def to_world(frame: np.ndarray, local_xyz_vec: np.ndarray) -> np.ndarray:
    """Plate-local unit vectors (..., 3) -> world unit vectors, given local frame `frame`."""
    return local_xyz_vec @ frame.T
