import numpy as np

from app import mantle


def test_cube_to_sphere_is_unit_length():
    rng = np.random.default_rng(0)
    for _ in range(50):
        face = ["+x", "-x", "+y", "-y", "+z", "-z"][rng.integers(0, 6)]
        u, v = rng.uniform(-1, 1, size=2)
        p = mantle.cube_to_sphere(face, u, v)
        assert np.isclose(np.linalg.norm(p), 1.0)


def test_cube_to_sphere_face_centers_are_axis_aligned():
    assert np.allclose(mantle.cube_to_sphere("+x", 0, 0), [1, 0, 0])
    assert np.allclose(mantle.cube_to_sphere("-z", 0, 0), [0, 0, -1])


def test_fit_euler_pole_recovers_known_omega_exactly():
    rng = np.random.default_rng(1)
    true_omega = np.array([0.1, -0.05, 0.2])
    points = rng.normal(size=(30, 3))
    points /= np.linalg.norm(points, axis=-1, keepdims=True)
    velocities = np.cross(true_omega, points)

    fitted = mantle.fit_euler_pole(points, velocities)
    assert np.allclose(fitted, true_omega, atol=1e-9)


def test_fit_euler_pole_is_robust_to_noise():
    rng = np.random.default_rng(2)
    true_omega = np.array([0.0, 0.0, 0.3])
    points = rng.normal(size=(200, 3))
    points /= np.linalg.norm(points, axis=-1, keepdims=True)
    velocities = np.cross(true_omega, points) + rng.normal(scale=1e-4, size=(200, 3))

    fitted = mantle.fit_euler_pole(points, velocities)
    assert np.allclose(fitted, true_omega, atol=1e-2)


def test_flow_at_is_tangential():
    rng = np.random.default_rng(3)
    centers = mantle.generate_convection_centers(rng, n_centers=6)
    points = rng.normal(size=(40, 3))
    points /= np.linalg.norm(points, axis=-1, keepdims=True)

    v = mantle.flow_at(points, centers)
    radial_component = np.sum(v * points, axis=-1)
    assert np.allclose(radial_component, 0.0, atol=1e-9)


def test_flow_at_center_itself_is_zero():
    center = mantle.ConvectionCenter(position=np.array([0.0, 0.0, 1.0]), strength=1.0, falloff=1.0)
    v = mantle.flow_at(np.array([[0.0, 0.0, 1.0]]), [center])
    assert np.allclose(v, 0.0)


def test_clamp_rate_leaves_zero_vector_unchanged():
    assert np.allclose(mantle.clamp_rate(np.zeros(3)), 0.0)


def test_clamp_rate_bounds_magnitude():
    below_min = np.array([1e-3 * mantle.MIN_PLATE_RATE, 0, 0])
    clamped = mantle.clamp_rate(below_min)
    assert np.isclose(np.linalg.norm(clamped), mantle.MIN_PLATE_RATE)

    huge = np.array([1.0, 0, 0])
    clamped = mantle.clamp_rate(huge)
    assert np.isclose(np.linalg.norm(clamped), mantle.MAX_PLATE_RATE)
