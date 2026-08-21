import numpy as np

from app import geometry


def test_latlon_xyz_roundtrip():
    rng = np.random.default_rng(0)
    lat = rng.uniform(-np.pi / 2, np.pi / 2, size=50)
    lon = rng.uniform(-np.pi, np.pi, size=50)
    xyz = geometry.latlon_to_xyz(lat, lon)
    assert np.allclose(np.linalg.norm(xyz, axis=-1), 1.0)
    lat2, lon2 = geometry.xyz_to_latlon(xyz)
    assert np.allclose(lat, lat2, atol=1e-9)
    assert np.allclose(lon, lon2, atol=1e-9)


def test_rotation_matrix_is_orthonormal():
    r = geometry.rotation_matrix(np.array([0.0, 0.0, 1.0]), 0.7)
    assert np.allclose(r @ r.T, np.eye(3), atol=1e-9)
    assert np.isclose(np.linalg.det(r), 1.0)


def test_rotation_about_z_matches_longitude_shift():
    v = geometry.latlon_to_xyz(np.radians(10.0), np.radians(20.0))
    r = geometry.rotation_matrix(np.array([0.0, 0.0, 1.0]), np.radians(30.0))
    rotated = r @ v
    lat, lon = geometry.xyz_to_latlon(rotated)
    assert np.isclose(lat, np.radians(10.0), atol=1e-9)
    assert np.isclose(lon, np.radians(50.0), atol=1e-9)


def test_rotation_preserves_angular_distances():
    rng = np.random.default_rng(1)
    lat = rng.uniform(-np.pi / 2, np.pi / 2, size=20)
    lon = rng.uniform(-np.pi, np.pi, size=20)
    pts = geometry.latlon_to_xyz(lat, lon)
    axis = geometry.normalize(rng.normal(size=3))
    angle = 1.234
    rotated = geometry.rotate_vectors(pts, axis, angle)

    before = geometry.angular_distance(pts[0], pts[1:])
    after = geometry.angular_distance(rotated[0], rotated[1:])
    assert np.allclose(before, after, atol=1e-9)


def test_rotation_matrix_from_omega_matches_rate_times_dt():
    omega = np.array([0.0, 0.0, 0.5])  # rate 0.5 rad per unit time, about +z
    r = geometry.rotation_matrix_from_omega(omega, dt=2.0)
    v = np.array([1.0, 0.0, 0.0])
    rotated = r @ v
    expected = np.array([np.cos(1.0), np.sin(1.0), 0.0])  # angle = 0.5 * 2.0 = 1.0 rad
    assert np.allclose(rotated, expected, atol=1e-9)


def test_plate_frame_maps_local_origin_to_seed():
    seed = geometry.normalize(np.array([1.0, 2.0, 3.0]))
    frame = geometry.plate_frame_from_seed(seed)
    assert np.allclose(frame @ np.array([1.0, 0.0, 0.0]), seed, atol=1e-9)
    # Frame must be a proper rotation (orthonormal, right-handed).
    assert np.allclose(frame @ frame.T, np.eye(3), atol=1e-9)
    assert np.isclose(np.linalg.det(frame), 1.0)


def test_plate_frame_handles_seed_near_global_pole():
    seed = geometry.normalize(np.array([1e-8, 1e-8, 1.0]))
    frame = geometry.plate_frame_from_seed(seed)
    assert np.allclose(frame @ frame.T, np.eye(3), atol=1e-9)
    assert np.allclose(frame @ np.array([1.0, 0.0, 0.0]), seed, atol=1e-6)


def test_local_world_roundtrip_via_frame():
    rng = np.random.default_rng(2)
    seed = geometry.normalize(rng.normal(size=3))
    frame = geometry.plate_frame_from_seed(seed)

    phi = rng.uniform(-np.pi / 2, np.pi / 2, size=30)
    theta = rng.uniform(-np.pi, np.pi, size=30)
    local_pts = geometry.local_xyz(phi, theta)

    world_pts = geometry.to_world(frame, local_pts)
    assert np.allclose(np.linalg.norm(world_pts, axis=-1), 1.0, atol=1e-9)

    back_to_local = geometry.to_local(frame, world_pts)
    assert np.allclose(back_to_local, local_pts, atol=1e-9)


def test_angular_distance_known_values():
    a = np.array([1.0, 0.0, 0.0])
    b = np.array([0.0, 1.0, 0.0])
    assert np.isclose(geometry.angular_distance(a, b), np.pi / 2)
    assert np.isclose(geometry.angular_distance(a, a), 0.0, atol=1e-9)
    assert np.isclose(geometry.angular_distance(a, -a), np.pi)


def test_local_tangent_basis_is_orthonormal_to_pole():
    pole = geometry.normalize(np.array([1.0, 2.0, 3.0]))
    east, north = geometry.local_tangent_basis(pole)
    assert np.isclose(np.dot(east, north), 0.0, atol=1e-9)
    assert np.isclose(np.dot(east, pole), 0.0, atol=1e-9)
    assert np.isclose(np.dot(north, pole), 0.0, atol=1e-9)
    assert np.isclose(np.linalg.norm(east), 1.0, atol=1e-9)
    assert np.isclose(np.linalg.norm(north), 1.0, atol=1e-9)


def test_local_tangent_basis_handles_pole_at_global_pole():
    pole = np.array([0.0, 0.0, 1.0])
    east, north = geometry.local_tangent_basis(pole)
    assert np.allclose(east @ east, 1.0, atol=1e-9)
    assert np.allclose(north @ north, 1.0, atol=1e-9)
    assert np.isclose(np.dot(east, north), 0.0, atol=1e-9)


def test_azimuthal_equidistant_roundtrips_and_preserves_radius_from_pole():
    rng = np.random.default_rng(3)
    pole = geometry.normalize(rng.normal(size=3))
    east, north = geometry.local_tangent_basis(pole)
    # Points within a modest angular range of the pole -- azimuthal equidistant's own
    # documented safe regime (see its docstring re: antipodal singularity).
    points = geometry.normalize(pole[None, :] * 0.7 + rng.normal(size=(40, 3)) * 0.3)

    xy = geometry.azimuthal_equidistant_forward(pole, east, north, points)
    back = geometry.azimuthal_equidistant_inverse(pole, east, north, xy)
    assert np.allclose(back, points, atol=1e-9)

    true_radius = geometry.angular_distance(points, pole)
    xy_radius = np.hypot(xy[:, 0], xy[:, 1])
    assert np.allclose(true_radius, xy_radius, atol=1e-9)
