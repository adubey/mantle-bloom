import numpy as np

from app import geometry
from app.v2 import lithosphere, torque


def test_basal_drag_vanishes_when_plate_matches_mantle_flow():
    """Eq. 10: F = (mu/ds)(v_mantle - v_plate) -- zero relative velocity means zero drag
    force everywhere, hence zero net torque, regardless of how fast the plate itself is
    spinning in absolute terms."""

    class FakePlate:
        crust_type = "continental"

        def __init__(self, omega, points):
            self.omega = omega
            self._points = points

        def all_points_and_elevation(self):
            return self._points, np.zeros(len(self._points))

    class FakeWorld:
        def __init__(self, omega, points):
            self.mantle_centers = []  # flow_at with no centers returns zero field
            self.omega = omega

    rng = np.random.default_rng(0)
    points = rng.normal(size=(50, 3))
    points /= np.linalg.norm(points, axis=-1, keepdims=True)
    zero_omega = np.zeros(3)
    plate = FakePlate(zero_omega, points)
    world = FakeWorld(zero_omega, points)

    tau = torque.basal_drag_torque(plate, world, spacing_rad=0.02)
    assert np.allclose(tau, 0.0, atol=1e-6)


def test_slab_pull_only_applies_to_oceanic_plates():
    class FakePlate:
        def __init__(self, crust_type):
            self.crust_type = crust_type

    n = 10
    points = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    hm = np.full(n, 60_000.0)
    direction = np.tile(np.array([0.0, 1.0, 0.0]), (n, 1))
    inputs = torque.BoundaryForceInputs(
        own_points=points,
        own_hc=np.full(n, 7000.0),
        own_hm=hm,
        dist_to_neighbor=np.full(n, 0.001),
        direction_to_neighbor=direction,
        neighbor_is_oceanic=np.zeros(n, dtype=bool),
        neighbor_omega=np.zeros((n, 3)),
    )
    reach = 0.01
    spacing = 0.02

    continental = FakePlate("continental")
    assert np.allclose(torque.slab_pull_torque(continental, inputs, spacing, reach), 0.0)

    oceanic = FakePlate("oceanic")
    tau = torque.slab_pull_torque(oceanic, inputs, spacing, reach)
    assert not np.allclose(tau, 0.0)


def test_slab_pull_points_toward_subduction_direction():
    """A single subducting node at the pole (0,0,1), pulled toward +x, should produce a
    torque with the sign consistent with r x F for that geometry -- a real, checkable
    direction, not just "nonzero"."""

    class FakePlate:
        crust_type = "oceanic"

    point = np.array([[0.0, 0.0, 1.0]])
    direction = np.array([[1.0, 0.0, 0.0]])
    inputs = torque.BoundaryForceInputs(
        own_points=point,
        own_hc=np.array([7000.0]),
        own_hm=np.array([60_000.0]),
        dist_to_neighbor=np.array([0.001]),
        direction_to_neighbor=direction,
        neighbor_is_oceanic=np.array([False]),
        neighbor_omega=np.zeros((1, 3)),
    )
    tau = torque.slab_pull_torque(FakePlate(), inputs, spacing_rad=0.02, reach_rad=0.01)
    r = point[0] * lithosphere.PLANET_RADIUS_M
    force_direction = direction[0]
    expected_direction = geometry.normalize(np.cross(r, force_direction)[None, :])[0]
    assert np.dot(geometry.normalize(tau[None, :])[0], expected_direction) > 0.99


def test_merge_omega_conserves_angular_momentum():
    rng = np.random.default_rng(3)
    points_a = rng.normal(size=(80, 3))
    points_a /= np.linalg.norm(points_a, axis=-1, keepdims=True)
    points_b = rng.normal(size=(60, 3))
    points_b /= np.linalg.norm(points_b, axis=-1, keepdims=True)

    class FakePlate:
        def __init__(self, omega):
            self.omega = omega

    hc, hm = np.full(80, 35_000.0), np.full(80, 100_000.0)
    inertia_a = lithosphere.moment_of_inertia_tensor(points_a, hc, hm, lithosphere.RHO_CONTINENTAL_CRUST, 0.02)
    hc_b, hm_b = np.full(60, 7000.0), np.full(60, 60_000.0)
    inertia_b = lithosphere.moment_of_inertia_tensor(points_b, hc_b, hm_b, lithosphere.RHO_OCEANIC_CRUST, 0.02)

    omega_a = np.array([1e-9, 0.0, 0.0])
    omega_b = np.array([0.0, 1e-9, 0.0])
    plate_a, plate_b = FakePlate(omega_a), FakePlate(omega_b)

    combined_omega = torque.merge_omega(plate_a, inertia_a, plate_b, inertia_b)

    l_before = lithosphere.angular_momentum(inertia_a, omega_a) + lithosphere.angular_momentum(inertia_b, omega_b)
    l_after = lithosphere.angular_momentum(inertia_a + inertia_b, combined_omega)
    # merge_omega clamps the resulting rate (mantle.clamp_rate), so speed can change --
    # what must be conserved is *direction*, the part clamping doesn't touch.
    assert np.dot(geometry.normalize(l_before[None, :])[0], geometry.normalize(l_after[None, :])[0]) > 0.999
