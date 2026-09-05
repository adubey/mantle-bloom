import numpy as np
import pytest

from app import geometry
from app import lithosphere, mantle, torque


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


def test_subducting_mask_only_flags_a_converging_oceanic_boundary():
    """`subducting_boundary_mask` gates slab pull: a continental plate never subducts, and
    even an oceanic plate only feels pull where it is actually converging on its neighbour
    (a ridge or transform stretch, however close, does not pull a slab down -- the bug that
    railed every oceanic plate at MAX_PLATE_RATE)."""

    class FakePlate:
        def __init__(self, crust_type, omega):
            self.crust_type = crust_type
            self.omega = omega

    n = 10
    points = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    direction = np.tile(np.array([0.0, 1.0, 0.0]), (n, 1))  # neighbour is toward +y
    # A plate spinning about +z carries +x-material toward +y -- i.e. toward the neighbour.
    converging_omega = np.array([0.0, 0.0, 2.0 * mantle.MAX_PLATE_RATE])
    inputs = torque.BoundaryForceInputs(
        own_points=points,
        own_hc=np.full(n, 7000.0),
        own_hm=np.full(n, 60_000.0),
        own_crust_type_codes=np.zeros(n, dtype=np.int8),
        dist_to_neighbor=np.full(n, 0.001),
        direction_to_neighbor=direction,
        neighbor_is_oceanic=np.zeros(n, dtype=bool),
        neighbor_omega=np.zeros((n, 3)),
    )
    reach = 0.01

    assert not torque.subducting_boundary_mask(FakePlate("continental", converging_omega), inputs, reach).any()
    assert torque.subducting_boundary_mask(FakePlate("oceanic", converging_omega), inputs, reach).all()
    # Same geometry, spinning the other way -> the boundary is diverging, no slab pull.
    assert not torque.subducting_boundary_mask(FakePlate("oceanic", -converging_omega), inputs, reach).any()


def test_slab_pull_is_zero_without_a_subducting_node():
    n = 10
    inputs = torque.BoundaryForceInputs(
        own_points=np.tile(np.array([1.0, 0.0, 0.0]), (n, 1)),
        own_hc=np.full(n, 7000.0),
        own_hm=np.full(n, 60_000.0),
        own_crust_type_codes=np.zeros(n, dtype=np.int8),
        dist_to_neighbor=np.full(n, 0.001),
        direction_to_neighbor=np.tile(np.array([0.0, 1.0, 0.0]), (n, 1)),
        neighbor_is_oceanic=np.zeros(n, dtype=bool),
        neighbor_omega=np.zeros((n, 3)),
    )
    assert np.allclose(torque.slab_pull_torque(None, inputs, np.zeros(n, dtype=bool), spacing_rad=0.02), 0.0)
    assert not np.allclose(torque.slab_pull_torque(None, inputs, np.ones(n, dtype=bool), spacing_rad=0.02), 0.0)


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
        own_crust_type_codes=np.zeros(1, dtype=np.int8),
        dist_to_neighbor=np.array([0.001]),
        direction_to_neighbor=direction,
        neighbor_is_oceanic=np.array([False]),
        neighbor_omega=np.zeros((1, 3)),
    )
    tau = torque.slab_pull_torque(FakePlate(), inputs, np.array([True]), spacing_rad=0.02)
    r = point[0] * lithosphere.PLANET_RADIUS_M
    force_direction = direction[0]
    expected_direction = geometry.normalize(np.cross(r, force_direction)[None, :])[0]
    assert np.dot(geometry.normalize(tau[None, :])[0], expected_direction) > 0.99


def test_classify_boundary_nodes_flags_deep_interior_overlap():
    """A node far (> reach_rad) from any neighbour node but geometrically *inside* a
    neighbour's territory -- a plate that has slid deep over another -- must still classify
    `contested`, so rheology thickens/subducts it instead of leaving the overlap frozen.
    Before the bounding-sphere widening, only the reach_rad-local band was ever polygon-
    tested."""
    center = geometry.normalize(np.array([1.0, 0.0, 0.0]))
    east, north = geometry.local_tangent_basis(center)
    # Neighbour node cloud: a ring ~0.3 rad around `center` (centroid ~= center, radius ~0.3).
    ang = np.linspace(0, 2 * np.pi, 24, endpoint=False)
    ring = geometry.normalize(
        np.cos(0.3) * center + np.sin(0.3) * (np.cos(ang)[:, None] * east + np.sin(ang)[:, None] * north)
    )

    class FakeNeighbour:
        crust_type = "continental"
        plate_id = 99

        def all_points_and_elevation(self):
            return ring, np.zeros(len(ring))

        def contains_batch(self, pts):
            return geometry.angular_distance(np.asarray(pts), center) < 0.15

    # Self plate: one node right at `center` (deep inside the ring, ~0.3 rad from any ring
    # node -> not "near"), one node out at the antipode (neither near nor inside).
    own = np.array([center, [-1.0, 0.0, 0.0]])
    reach_rad = 0.06
    inputs = torque.BoundaryForceInputs(
        own_points=own,
        own_hc=np.full(2, 35_000.0),
        own_hm=np.full(2, 100_000.0),
        own_crust_type_codes=np.zeros(2, dtype=np.int8),
        dist_to_neighbor=np.array([np.inf, np.inf]),
        direction_to_neighbor=np.tile(np.array([0.0, 1.0, 0.0]), (2, 1)),
        neighbor_is_oceanic=np.zeros(2, dtype=bool),
        neighbor_omega=np.zeros((2, 3)),
    )

    class SelfPlate:
        crust_type = "continental"
        omega = np.zeros(3)

    convergent, divergent, transform, contested = torque.classify_boundary_nodes(
        SelfPlate(), [FakeNeighbour()], inputs, reach_rad
    )
    assert bool(contested[0]) and not bool(contested[1])
    assert bool(convergent[0])  # a deep-interior overlap folds into convergent (thickens / subducts)
    assert not divergent.any()  # ... and is never divergent
    assert not transform.any()  # ... nor transform -- it is outside the near-boundary band entirely


def test_classify_boundary_nodes_is_motion_based_not_geometric():
    """A node in the near-boundary band that is not overlapping any polygon is classified by
    the *sign* of its closing rate: converging -> convergent (builds an orogen before any
    overlap accumulates), opening -> divergent, near-zero -> transform."""
    # Two nodes on the equator, a neighbour node just east of each (within reach_rad).
    own = geometry.normalize(np.array([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]]))
    reach_rad = 0.06
    nb = geometry.normalize(np.array([[np.cos(0.03), np.sin(0.03), 0.0]] * 2))
    direction = geometry.normalize(nb - own)
    inputs = torque.BoundaryForceInputs(
        own_points=own,
        own_hc=np.full(2, 35_000.0),
        own_hm=np.full(2, 100_000.0),
        own_crust_type_codes=np.zeros(2, dtype=np.int8),
        dist_to_neighbor=np.array([0.03, 0.03]),
        direction_to_neighbor=direction,
        neighbor_is_oceanic=np.zeros(2, dtype=bool),
        neighbor_omega=np.zeros((2, 3)),
    )

    class Neighbour:
        crust_type = "continental"
        plate_id = 1

        def all_points_and_elevation(self):
            return nb, np.zeros(2)

        def contains_batch(self, pts):
            return np.zeros(len(pts), dtype=bool)  # never overlapping

    # Self spinning about +z so its equatorial material moves +y, toward the neighbour -> both
    # nodes converge.
    fast = 5.0 * torque.BOUNDARY_TRANSFORM_RATE_M_PER_S * torque.SECONDS_PER_YEAR / lithosphere.PLANET_RADIUS_M

    class SelfPlate:
        crust_type = "continental"
        omega = np.array([0.0, 0.0, fast])

    convergent, divergent, transform, contested = torque.classify_boundary_nodes(
        SelfPlate(), [Neighbour()], inputs, reach_rad
    )
    assert convergent.all() and not contested.any()  # converging, though nothing overlaps
    assert not divergent.any() and not transform.any()

    class SelfPlateOpening(SelfPlate):
        omega = np.array([0.0, 0.0, -fast])  # equatorial material moves -y, away

    convergent, divergent, transform, _ = torque.classify_boundary_nodes(
        SelfPlateOpening(), [Neighbour()], inputs, reach_rad
    )
    assert divergent.all() and not convergent.any() and not transform.any()

    class SelfPlateStill(SelfPlate):
        omega = np.zeros(3)  # no relative motion -> transform

    convergent, divergent, transform, _ = torque.classify_boundary_nodes(
        SelfPlateStill(), [Neighbour()], inputs, reach_rad
    )
    assert transform.all() and not convergent.any() and not divergent.any()


def test_basal_drag_coefficients_reproduce_the_plain_torque():
    """`basal_drag_torque(omega) == b - K @ omega` exactly, for arbitrary omega -- the affine
    split `integrate_omega` integrates implicitly must be the same physics as the direct
    evaluation, not an approximation."""

    class FakePlate:
        crust_type = "oceanic"

        def __init__(self, omega, points):
            self.omega = omega
            self._points = points

        def all_points_and_elevation(self):
            return self._points, np.zeros(len(self._points))

    class FakeWorld:
        mantle_centers = mantle.generate_convection_centers(np.random.default_rng(1), n_centers=8)

    rng = np.random.default_rng(0)
    points = geometry.normalize(rng.normal(size=(200, 3)))
    world = FakeWorld()
    b, k = torque.basal_drag_coefficients(FakePlate(np.zeros(3), points), world, spacing_rad=0.02)

    for omega in [np.zeros(3), np.array([2e-8, -1e-8, 3e-9]), rng.normal(size=3) * 1e-8]:
        direct = torque.basal_drag_torque(FakePlate(omega, points), world, spacing_rad=0.02)
        assert np.allclose(direct, b - k @ omega, rtol=1e-9, atol=1e-6 * np.linalg.norm(b))


def test_integrate_omega_relaxes_to_the_mantle_rate_instead_of_railing():
    """The historical bug: an explicit Euler step on the (very stiff) basal-drag term
    overshoots by many orders of magnitude, so `clamp_rate` pins every plate at
    `MAX_PLATE_RATE` on the first step and forever after. With drag integrated implicitly a
    drag-only plate settles at its local mantle-flow rate -- well below the clamp -- and
    stays put across many steps regardless of step size."""

    class FakePlate:
        crust_type = "oceanic"

        def __init__(self, omega, points):
            self.omega = omega
            self._points = points

        def all_points_and_elevation(self):
            return self._points, np.zeros(len(self._points))

    class FakeWorld:
        mantle_centers = mantle.generate_convection_centers(np.random.default_rng(936513024), n_centers=8)

    rng = np.random.default_rng(0)
    # A compact cap, not a whole-sphere scatter, so the mantle flow across it is coherent
    # and the fitted rate is a meaningful single number.
    center = geometry.normalize(np.array([1.0, 0.3, 0.2]))
    tangent = geometry.normalize(np.cross(center, np.array([0.0, 0.0, 1.0])))
    bitangent = np.cross(center, tangent)
    ang = rng.uniform(0.0, 0.5, size=400) ** 0.5 * 0.5
    azi = rng.uniform(0.0, 2 * np.pi, size=400)
    points = geometry.normalize(
        np.cos(ang)[:, None] * center
        + (np.sin(ang) * np.cos(azi))[:, None] * tangent
        + (np.sin(ang) * np.sin(azi))[:, None] * bitangent
    )
    world = FakeWorld()

    hc = np.full(len(points), lithosphere.REFERENCE_HC_OCEANIC_M)
    hm = np.full(len(points), lithosphere.REFERENCE_HM_OCEANIC_M)
    inertia = lithosphere.moment_of_inertia_tensor(points, hc, hm, lithosphere.RHO_OCEANIC_CRUST, 0.02)

    fitted_rate = np.linalg.norm(mantle.fit_euler_pole(points, mantle.flow_at(points, world.mantle_centers)))
    assert fitted_rate < mantle.MAX_PLATE_RATE  # the field itself does not demand a railed plate

    plate = FakePlate(np.zeros(3), points)
    for years in (10_000.0, 100_000.0, 3_000_000.0):
        for _ in range(50):
            b, k = torque.basal_drag_coefficients(plate, world, spacing_rad=0.02)
            plate.omega = torque.integrate_omega(plate, np.zeros(3), b, k, inertia, years)
        assert np.linalg.norm(plate.omega) == pytest.approx(fitted_rate, rel=1e-3)
        assert np.linalg.norm(plate.omega) < 0.5 * mantle.MAX_PLATE_RATE


def test_slab_drag_keeps_a_subducting_plate_off_the_clamp():
    """Slab pull with no resistance drives a subducting oceanic plate straight past
    `MAX_PLATE_RATE`, so `integrate_omega`'s `clamp_rate` pins it exactly at the clamp (the
    'every oceanic plate railed at MAX' reading). Folding `slab_drag_coefficient_matrix` into
    the implicit `K` -- the sunk slab's viscous coupling to the deep mantle -- lets the speed
    self-regulate below the clamp instead, and the drag is unconditionally stable at any step
    size the way basal drag is."""
    rng = np.random.default_rng(495717634)
    center = geometry.normalize(np.array([1.0, 0.2, -0.3]))
    tangent = geometry.normalize(np.cross(center, np.array([0.0, 0.0, 1.0])))
    bitangent = np.cross(center, tangent)
    ang = rng.uniform(0.0, 1.0, size=300) ** 0.5 * 0.4
    azi = rng.uniform(0.0, 2 * np.pi, size=300)
    points = geometry.normalize(
        np.cos(ang)[:, None] * center
        + (np.sin(ang) * np.cos(azi))[:, None] * tangent
        + (np.sin(ang) * np.sin(azi))[:, None] * bitangent
    )
    n = len(points)
    # Every node subducts toward +tangent -- a whole-perimeter trench, the worst case.
    inputs = torque.BoundaryForceInputs(
        own_points=points,
        own_hc=np.full(n, lithosphere.REFERENCE_HC_OCEANIC_M),
        own_hm=np.full(n, lithosphere.REFERENCE_HM_OCEANIC_M),
        own_crust_type_codes=np.zeros(n, dtype=np.int8),
        dist_to_neighbor=np.full(n, 0.001),
        direction_to_neighbor=np.tile(tangent, (n, 1)),
        neighbor_is_oceanic=np.ones(n, dtype=bool),
        neighbor_omega=np.zeros((n, 3)),
    )
    subducting = np.ones(n, dtype=bool)
    inertia = lithosphere.moment_of_inertia_tensor(
        points, inputs.own_hc, inputs.own_hm, lithosphere.RHO_OCEANIC_CRUST, 0.02
    )
    pull = torque.slab_pull_torque(None, inputs, subducting, spacing_rad=0.02)

    class FakePlate:
        omega = np.zeros(3)

    without_drag = torque.integrate_omega(FakePlate(), pull, np.zeros(3), np.zeros((3, 3)), inertia, years=100_000.0)
    assert np.linalg.norm(without_drag) == pytest.approx(mantle.MAX_PLATE_RATE, rel=1e-9)

    drag_k = torque.slab_drag_coefficient_matrix(inputs, subducting, spacing_rad=0.02)
    plate = FakePlate()
    for years in (10_000.0, 100_000.0, 5_000_000.0):
        settled = plate.omega
        for _ in range(60):
            settled = torque.integrate_omega(plate, pull, np.zeros(3), drag_k, inertia, years)
            plate = type("P", (), {"omega": settled})()
        assert np.linalg.norm(settled) < mantle.MAX_PLATE_RATE


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


def test_collision_friction_torque_brakes_harder_for_a_more_severe_overlap():
    """`overlap_severity` (see OVERLAP_FRICTION_SEVERITY_GAIN) is a deliberate *extra* brake on
    top of what more contested nodes already sum to -- confirmed here by holding the
    contested-node geometry/velocity fixed and only varying the severity argument."""
    n = 5
    own_points = np.tile(np.array([1.0, 0.0, 0.0]), (n, 1))
    collision_mask = np.ones(n, dtype=bool)
    # A neighbour spinning the opposite way -> real relative motion at the contested nodes,
    # so the resistive torque is nonzero to begin with.
    inputs = torque.BoundaryForceInputs(
        own_points=own_points,
        own_hc=np.full(n, 35_000.0),
        own_hm=np.full(n, 100_000.0),
        own_crust_type_codes=np.zeros(n, dtype=np.int8),
        dist_to_neighbor=np.full(n, 0.001),
        direction_to_neighbor=np.tile(np.array([0.0, 1.0, 0.0]), (n, 1)),
        neighbor_is_oceanic=np.zeros(n, dtype=bool),
        neighbor_omega=np.tile(np.array([0.0, 0.0, -0.02]), (n, 1)),
    )

    class FakePlate:
        crust_type = "continental"
        omega = np.array([0.0, 0.0, 0.02])

    plate = FakePlate()
    light = torque.collision_friction_torque(plate, inputs, collision_mask, spacing_rad=0.02, overlap_severity=0.0)
    severe = torque.collision_friction_torque(plate, inputs, collision_mask, spacing_rad=0.02, overlap_severity=1.0)

    assert np.linalg.norm(severe) > np.linalg.norm(light)
    # Same geometry/velocities throughout -- the resistive force still opposes relative
    # motion, it's just scaled up, not redirected.
    assert np.dot(geometry.normalize(severe[None, :])[0], geometry.normalize(light[None, :])[0]) > 0.999
