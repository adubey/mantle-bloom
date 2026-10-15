from types import SimpleNamespace

import numpy as np
import pytest

from app import geometry
from app import cratons, lithosphere, mantle, torque
from app.elevation_lines import CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC, line_spacing_rad


class _OceanicPlate:
    """Only slab pull's crust-code resolution reads the plate: inherited codes are oceanic."""

    crust_type = "oceanic"


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

        def accounting_areas_m2(self, spacing_rad):
            return np.full(len(self._points), lithosphere.node_area_m2(spacing_rad))

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
    assert np.allclose(torque.slab_pull_torque(_OceanicPlate(), inputs, np.zeros(n, dtype=bool), spacing_rad=0.02), 0.0)
    assert not np.allclose(torque.slab_pull_torque(_OceanicPlate(), inputs, np.ones(n, dtype=bool), spacing_rad=0.02), 0.0)


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

        def accounting_areas_m2(self, spacing_rad):
            return np.full(len(self._points), lithosphere.node_area_m2(spacing_rad))

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

        def accounting_areas_m2(self, spacing_rad):
            return np.full(len(self._points), lithosphere.node_area_m2(spacing_rad))

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
    pull = torque.slab_pull_torque(_OceanicPlate(), inputs, subducting, spacing_rad=0.02)

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


def _contact_band(spacing_rad: float, rows: int, hc: float = 70_000.0, hm: float = 200_000.0, omega=(0.0, 0.0, 0.0)):
    """A straight north-south contact on the equator: `rows` nodes along it, one per spacing,
    and BOUNDARY_FORCE_REACH_MULTIPLIER nodes across the band. The neighbour is still, so the
    plate's omega about z is the relative motion, normal to the contact. Hc/Hm default to 2x
    reference, past COLLISION_RELIEF_FRICTION_SCALE_M of relief, so g = 1."""
    across = int(torque.BOUNDARY_FORCE_REACH_MULTIPLIER)
    lat = (np.arange(rows) - (rows - 1) / 2) * spacing_rad
    lon = np.arange(across) * spacing_rad
    lat_grid, lon_grid = (g.ravel() for g in np.meshgrid(lat, lon, indexing="ij"))
    points = np.column_stack([np.cos(lat_grid) * np.cos(lon_grid), np.cos(lat_grid) * np.sin(lon_grid), np.sin(lat_grid)])
    n = len(points)
    inputs = torque.BoundaryForceInputs(
        own_points=points,
        own_hc=np.full(n, hc),
        own_hm=np.full(n, hm),
        own_crust_type_codes=np.zeros(n, dtype=np.int8),
        dist_to_neighbor=np.full(n, 0.001),
        direction_to_neighbor=np.tile(np.array([0.0, 1.0, 0.0]), (n, 1)),
        neighbor_is_oceanic=np.zeros(n, dtype=bool),
        neighbor_omega=np.zeros((n, 3)),
    )

    class FakePlate:
        crust_type = "continental"

        def __init__(self):
            self.omega = np.array(omega, dtype=float)

        def accounting_areas_m2(self, spacing):
            return np.full(n, lithosphere.node_area_m2(spacing))

    return FakePlate(), inputs, np.ones(n, dtype=bool)


@pytest.mark.parametrize("density", [1.0, 4.0])
def test_fully_braked_contact_resists_with_the_himalaya_tibet_force_per_trench_length(density):
    """RFC 0008 (#399): a fully braked contact at the reference 4 cm/yr resists with
    CONTACT_RESISTANCE_N_PER_M per metre of trench (Copley et al. 2010), whatever the node
    density, because the coefficient is spread over the band's physical width."""
    spacing = line_spacing_rad(density)
    rows = 5
    rate = mantle.cm_per_yr_to_rad_per_yr(torque.CONTACT_RESISTANCE_REF_CM_PER_YR)
    plate, inputs, mask = _contact_band(spacing, rows, omega=(0.0, 0.0, rate))

    tau = torque.collision_friction_torque(plate, inputs, mask, spacing)

    trench_m = rows * spacing * lithosphere.PLANET_RADIUS_M
    # Torque about z over the lever arm R: the summed resisting force, opposing the motion.
    force_per_m = -tau[2] / lithosphere.PLANET_RADIUS_M / trench_m
    assert force_per_m == pytest.approx(torque.CONTACT_RESISTANCE_N_PER_M, rel=0.01)


def test_contact_drag_factor_is_the_larger_of_relief_and_craton_and_stays_bounded():
    """One bounded factor g = max(relief, craton strength) in [0, 1] replaces the multiplicative
    gains (RFC 0008, #399): flat crust without craton doesn't brake, and no combination of relief
    and craton brakes harder than a full Tibet-height or cratonic contact."""
    spacing = line_spacing_rad(1.0)
    full = torque.CONTACT_RESISTANCE_N_PER_M / (
        torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing * lithosphere.PLANET_RADIUS_M * torque.CONTACT_RESISTANCE_REF_CM_PER_YR * 0.01 / torque.SECONDS_PER_YEAR
    )
    for hc, hm in ((35_000.0, 100_000.0), (50_000.0, 140_000.0), (70_000.0, 200_000.0), (84_000.0, 240_000.0)):
        for craton in (0.0, 10_000.0, 35_000.0, 80_000.0):
            plate, inputs, mask = _contact_band(spacing, 2, hc=hc, hm=hm)
            craton_m = np.full(len(mask), craton)
            g = torque.contact_drag_pa_s_per_m(plate, inputs, mask, spacing, craton_m) / full
            relief_only = torque.contact_drag_pa_s_per_m(plate, inputs, mask, spacing) / full
            assert np.all((g >= 0.0) & (g <= 1.0 + 1e-12))
            np.testing.assert_allclose(g, np.maximum(relief_only, cratons.strength(craton_m)))

    flat, inputs, mask = _contact_band(spacing, 2, hc=35_000.0, hm=100_000.0)
    assert np.all(torque.contact_drag_pa_s_per_m(flat, inputs, mask, spacing, np.zeros(len(mask))) == 0.0)
    assert np.all(torque.contact_drag_pa_s_per_m(flat, inputs, mask, spacing, np.full(len(mask), 80_000.0)) == pytest.approx(full))


def test_contact_drag_slows_but_never_reverses_a_plate_in_one_step():
    """The drag is implicit and linear in the relative velocity, so even a long step against a
    fully braked contact leaves the plate moving the same way, only slower (RFC 0008, #399)."""
    spacing = line_spacing_rad(1.0)
    rate = mantle.cm_per_yr_to_rad_per_yr(4.0)
    plate, inputs, mask = _contact_band(spacing, 5, omega=(0.0, 0.0, rate))
    # A tiny plate (the band alone), so the contact drag dominates its inertia.
    rho = np.full(len(mask), lithosphere.RHO_CONTINENTAL_CRUST)
    inertia = lithosphere.moment_of_inertia_tensor(inputs.own_points, inputs.own_hc, inputs.own_hm, rho, spacing)
    b, k = torque.collision_drag_coefficients(plate, inputs, mask, spacing)
    for years in (1e3, 1e6, 5e6, 1e8):
        new_omega = torque.integrate_omega(plate, np.zeros(3), b, k, inertia, years)
        assert 0.0 <= new_omega[2] < rate
        assert np.dot(new_omega, plate.omega) >= 0.0


def test_contact_craton_is_the_weaker_side_so_both_plates_brake_alike():
    # #340: collision resistance is one force across a contact, set by the weaker side.
    class FakePlate:
        def __init__(self, plate_id, craton):
            self.plate_id, self.craton = plate_id, np.asarray(craton, dtype=float)

        def collect(self, name):
            assert name == "craton_crust_m"
            return self.craton

    a = FakePlate(1, [30_000.0, 30_000.0, 0.0, 30_000.0])
    b = FakePlate(2, [20_000.0, 0.0])
    c = FakePlate(3, [5_000.0])
    inputs = torque.BoundaryForceInputs(
        own_points=np.zeros((4, 3)), own_hc=np.zeros(4), own_hm=np.zeros(4), own_crust_type_codes=np.zeros(4, dtype=np.int8),
        dist_to_neighbor=np.zeros(4), direction_to_neighbor=np.zeros((4, 3)), neighbor_is_oceanic=np.zeros(4, dtype=bool),
        neighbor_omega=np.zeros((4, 3)), neighbor_plate_id=np.array([2, 2, 3, -1]), neighbor_node_index=np.array([0, 1, 0, -1]),
    )

    np.testing.assert_array_equal(torque.contact_craton_m(a, [b, c], inputs), [20_000.0, 0.0, 0.0, 0.0])


def _jam_band(spacing_rad: float, rows: int, hc: float, hm: float = 100_000.0, rate_cm_per_yr: float = 0.0, craton: float = 0.0):
    """`_contact_band`, but on an oceanic plate whose middle `rows` rows are a continental-coded
    terrane (RFC 0008, #400), with a row of its ocean floor at each end of the same contact: the
    slab that drags the terrane down (#447). The terrane is so forced down, and in the jam mask
    while the contact isn't diverging. Returns the terrane mask."""
    rate = mantle.cm_per_yr_to_rad_per_yr(rate_cm_per_yr)
    plate, inputs, _ = _contact_band(spacing_rad, rows + 2, hc=hc, hm=hm, omega=(0.0, 0.0, rate))
    n = len(inputs.own_points)
    row = np.arange(n) // int(torque.BOUNDARY_FORCE_REACH_MULTIPLIER)
    terrane = (row >= 1) & (row <= rows)
    plate.crust_type = "oceanic"
    inputs.own_crust_type_codes = np.where(terrane, CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC).astype(np.int8)
    inputs.own_craton_m = np.where(terrane, craton, 0.0)
    inputs.neighbor_plate_id = np.full(n, 2)
    return plate, inputs, terrane


def _pull_per_node(plate, inputs, spacing):
    return torque.slab_pull_force_n(plate, inputs, np.ones(len(inputs.own_points), dtype=bool), spacing)


def test_ocean_floor_pulls_exactly_as_before_buoyancy():
    # RFC 0008 (#400): with no continental crust and no craton, B is the mantle term alone.
    spacing = line_spacing_rad(1.0)
    plate, inputs, _ = _contact_band(spacing, 2, hc=lithosphere.REFERENCE_HC_OCEANIC_M, hm=lithosphere.REFERENCE_HM_OCEANIC_M)
    plate.crust_type = "oceanic"
    expected = (
        (lithosphere.RHO_LITHOSPHERE_MANTLE - lithosphere.RHO_ASTHENOSPHERE)
        * lithosphere.GRAVITY_M_S2
        * lithosphere.REFERENCE_HM_OCEANIC_M
        * torque.SUBDUCTION_LSINK_M
        * spacing
        * lithosphere.PLANET_RADIUS_M
    )
    np.testing.assert_allclose(_pull_per_node(plate, inputs, spacing), expected, rtol=1e-12)


def test_thick_continental_crust_pulls_nothing_and_brakes_while_thin_crust_pulls_without_braking():
    """Cloos (1993): continental crust past ~15 km makes the column buoyant enough to jam the
    trench. A 35 km terrane on 100 km of mantle lithosphere pulls nothing and jams fully; a 10 km
    one on a thick 150 km root still sinks, more weakly than bare ocean floor, and doesn't jam."""
    spacing = line_spacing_rad(1.0)
    reach = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing
    thick, inputs, terrane = _jam_band(spacing, 2, hc=35_000.0, rate_cm_per_yr=3.0)
    assert np.all(_pull_per_node(thick, inputs, spacing)[terrane] == 0.0)
    np.testing.assert_allclose(torque.contact_jam_factor(thick, None, [], inputs, reach)[terrane], 1.0)

    thin, inputs, _ = _jam_band(spacing, 2, hc=10_000.0, hm=150_000.0, rate_cm_per_yr=3.0)
    bare, bare_inputs, _ = _jam_band(spacing, 2, hc=0.0, hm=150_000.0, rate_cm_per_yr=3.0)
    pull = _pull_per_node(thin, inputs, spacing)[terrane]
    assert np.all((pull > 0.0) & (pull < _pull_per_node(bare, bare_inputs, spacing)[terrane]))
    assert np.all(torque.contact_jam_factor(thin, None, [], inputs, reach) == 0.0)


def test_full_strength_cratonic_keel_adds_no_pull():
    # Jordan (1978): a depleted keel is neutrally buoyant, so a thicker one doesn't pull harder.
    hm = np.array([100_000.0, 250_000.0])
    full = np.full(2, cratons.CRATON_FULL_STRENGTH_HC_M)
    keel = torque.column_buoyancy_pa(np.zeros(2), hm, full)
    np.testing.assert_allclose(keel, 0.0)
    half = torque.column_buoyancy_pa(np.zeros(2), hm, full / 2)
    plain = torque.column_buoyancy_pa(np.zeros(2), hm, np.zeros(2))
    np.testing.assert_allclose(half, plain / 2)


def test_jam_and_contact_drag_factors_stay_bounded():
    spacing = line_spacing_rad(1.0)
    hc = np.linspace(-5_000.0, 90_000.0, 200)
    j = torque.jam_factor(hc)
    assert np.all((j >= 0.0) & (j <= 1.0))
    assert np.all(j[hc <= torque.JAM_MIN_CRUST_M] == 0.0) and np.all(j[hc >= torque.JAM_FULL_CRUST_M] == 1.0)
    full = torque.contact_drag_pa_s_per_m(*_contact_band(spacing, 1)[:3], spacing)[0]
    rng = np.random.default_rng(400)
    for hc_band, hm_band in ((35_000.0, 100_000.0), (70_000.0, 200_000.0)):
        plate, inputs, mask = _contact_band(spacing, 4, hc=hc_band, hm=hm_band)
        mask = rng.random(len(mask)) < 0.5
        g = torque.contact_drag_pa_s_per_m(plate, inputs, mask, spacing, rng.uniform(0, 80_000.0, len(mask)), rng.uniform(0, 1, len(mask))) / full
        assert np.all((g >= 0.0) & (g <= 1.0 + 1e-12))


@pytest.mark.parametrize("density", [1.0, 4.0])
def test_fully_jammed_contact_resists_with_the_himalaya_tibet_force_per_trench_length(density):
    # RFC 0008 (#400): j = 1 brakes as hard as solution 0's full brake, with no relief or craton.
    spacing = line_spacing_rad(density)
    rows = 5
    plate, inputs, _ = _jam_band(spacing, rows, hc=35_000.0, rate_cm_per_yr=torque.CONTACT_RESISTANCE_REF_CM_PER_YR)
    reach = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing
    jam = torque.contact_jam_factor(plate, None, [], inputs, reach)
    tau = torque.collision_friction_torque(plate, inputs, np.zeros(len(jam), dtype=bool), spacing, jam=jam)
    force_per_m = -tau[2] / lithosphere.PLANET_RADIUS_M / (rows * spacing * lithosphere.PLANET_RADIUS_M)
    assert force_per_m == pytest.approx(torque.CONTACT_RESISTANCE_N_PER_M, rel=0.01)


def test_jam_brake_slows_but_never_reverses_a_plate_in_one_step():
    spacing = line_spacing_rad(1.0)
    plate, inputs, _ = _jam_band(spacing, 5, hc=35_000.0, rate_cm_per_yr=4.0)
    jam = torque.contact_jam_factor(plate, None, [], inputs, torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing)
    rho = np.full(len(jam), lithosphere.RHO_CONTINENTAL_CRUST)
    inertia = lithosphere.moment_of_inertia_tensor(inputs.own_points, inputs.own_hc, inputs.own_hm, rho, spacing)
    b, k = torque.collision_drag_coefficients(plate, inputs, np.zeros(len(jam), dtype=bool), spacing, jam=jam)
    for years in (1e3, 1e6, 5e6, 1e8):
        new_omega = torque.integrate_omega(plate, np.zeros(3), b, k, inertia, years)
        assert 0.0 <= new_omega[2] < plate.omega[2]


def test_jam_holds_through_a_stalled_contact_and_releases_only_when_it_diverges():
    """The jam mask runs through zero closing (#400): slab pull's 1 cm/yr activation threshold
    would drop the brake as soon as it slowed the contact, and the contact would speed back up.
    Across that threshold the jam stays on and the resisting torque is continuous."""
    spacing = line_spacing_rad(1.0)
    reach = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing
    for rate, jammed in ((3.0, True), (1.01, True), (0.99, True), (0.0, True), (-0.99, True), (-1.01, False), (-3.0, False)):
        plate, inputs, terrane = _jam_band(spacing, 2, hc=35_000.0, rate_cm_per_yr=rate)
        assert np.all(torque.contact_jam_factor(plate, None, [], inputs, reach) == np.where(terrane & jammed, 1.0, 0.0))
    no_collision = np.zeros(len(inputs.own_points), dtype=bool)

    def resisting(rate):
        plate, inputs, _ = _jam_band(spacing, 2, hc=35_000.0, rate_cm_per_yr=rate)
        jam = torque.contact_jam_factor(plate, None, [], inputs, reach)
        return -torque.collision_friction_torque(plate, inputs, no_collision, spacing, jam=jam)[2], torque.subducting_boundary_mask(plate, inputs, reach)

    (above, subducting_above), (below, subducting_below) = resisting(1.01), resisting(0.99)
    assert subducting_above.all() and not subducting_below.any()
    assert below == pytest.approx(above * 0.99 / 1.01, rel=1e-9)

    # A plate driven past 1 cm/yr settles on the speed where driving and the jam balance,
    # below the threshold, and stays there rather than speeding back up.
    plate, inputs, _ = _jam_band(spacing, 5, hc=35_000.0, rate_cm_per_yr=3.0)
    rho = np.full(len(inputs.own_points), lithosphere.RHO_CONTINENTAL_CRUST)
    inertia = lithosphere.moment_of_inertia_tensor(inputs.own_points, inputs.own_hc, inputs.own_hm, rho, spacing)
    jam = torque.contact_jam_factor(plate, None, [], inputs, reach)
    _, k = torque.collision_drag_coefficients(plate, inputs, no_collision, spacing, jam=jam)
    target = mantle.cm_per_yr_to_rad_per_yr(0.5)
    driving = np.array([0.0, 0.0, k[2, 2] * target])
    rates = []
    for _ in range(20):
        jam = torque.contact_jam_factor(plate, None, [], inputs, reach)
        b, k = torque.collision_drag_coefficients(plate, inputs, no_collision, spacing, jam=jam)
        plate.omega = torque.integrate_omega(plate, driving, b, k, inertia, 1e6)
        rates.append(plate.omega[2])
    assert all(later <= earlier + 1e-30 for earlier, later in zip(rates, rates[1:]))
    assert rates[-1] == pytest.approx(target, rel=1e-6)


class _JamPlate:
    def __init__(self, plate_id, crust_type, hc, codes, cell_keys=None):
        self.plate_id, self.crust_type = plate_id, crust_type
        self.fields = {"crustal_thickness_m": np.asarray(hc, dtype=float), "crust_type_code": np.asarray(codes, dtype=np.int8)}
        self.cell_keys = np.asarray(cell_keys if cell_keys is not None else np.arange(len(hc)), dtype=np.int64)

    def collect(self, name):
        return self.fields[name]


def _frame(**lower_keys):
    """A polarity frame whose plate `p<id>` has lower-plate cells `lower_keys["p<id>"]`."""
    masks = {}
    for name, keys in lower_keys.items():
        keys = np.asarray(keys, dtype=np.int64)
        masks[int(name[1:])] = SimpleNamespace(node_keys=np.concatenate([keys, keys + 1000]), lower=np.arange(2 * len(keys)) < len(keys))
    return SimpleNamespace(collision_polarity_frame=SimpleNamespace(masks=masks))


def test_lower_plate_side_is_the_last_frozen_lower_plate_cell():
    # Plates shift before the prepass, so a plate's lower side is the last frame's, matched by
    # cell key: a lost cell drops out, a new one isn't lower until the next prepass.
    land = _JamPlate(2, "continental", [35_000.0] * 4, [0] * 4, cell_keys=[11, 12, 1010, 10])
    np.testing.assert_array_equal(torque.lower_plate_mask(land, _frame(p2=[10, 11, 99]), 4), [True, False, False, True])
    np.testing.assert_array_equal(torque.lower_plate_mask(land, None, 4), [False] * 4)
    np.testing.assert_array_equal(torque.lower_plate_mask(land, _frame(p3=[10]), 4), [False] * 4)


def _mixed_contact(plate_type: str, codes, neighbour_ids, rate_cm_per_yr: float = 2.0):
    """One node per entry on a still neighbour's contact, closing at `rate_cm_per_yr`: a plate
    of nominal `plate_type` whose nodes carry `codes`, each facing plate `neighbour_ids[i]`."""
    n = len(codes)
    plate, inputs, _ = _contact_band(line_spacing_rad(1.0), n, omega=(0.0, 0.0, mantle.cm_per_yr_to_rad_per_yr(rate_cm_per_yr)))
    keep = np.arange(n) * int(torque.BOUNDARY_FORCE_REACH_MULTIPLIER)  # the node on the contact, per row
    for name in ("own_points", "own_hc", "own_hm", "dist_to_neighbor", "direction_to_neighbor", "neighbor_is_oceanic", "neighbor_omega"):
        setattr(inputs, name, getattr(inputs, name)[keep])
    inputs.own_crust_type_codes = np.asarray(codes, dtype=np.int8)
    inputs.own_craton_m = np.zeros(n)
    inputs.neighbor_plate_id = np.asarray(neighbour_ids, dtype=np.int64)
    inputs.neighbor_node_index = np.zeros(n, dtype=np.int64)
    plate.crust_type = plate_type
    return plate, inputs


def test_the_down_going_side_follows_each_nodes_crust_and_its_plates_slab_not_the_majority_type():
    """#447: oceanic nodes go down on any plate; a continental node goes down only behind its
    own plate's slab at that contact (an oceanic node facing the same neighbour). Before, an
    oceanic-majority plate went down everywhere and a continental-majority one nowhere."""
    reach = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * line_spacing_rad(1.0)
    c, o, inherit = CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC, 0
    # A continent with its own ocean floor facing plate 2 and a stretch of plain margin facing
    # plate 3: the ocean floor and the continent behind it at plate 2 go down, not the margin.
    plate, inputs = _mixed_contact("continental", [o, inherit, inherit], [2, 2, 3])
    np.testing.assert_array_equal(torque.down_going_side(plate, inputs, reach), [True, True, False])
    np.testing.assert_array_equal(torque.subducting_boundary_mask(plate, inputs, reach), [True, True, False])
    # A terrane on ocean floor goes down where the plate subducts under that neighbour, not
    # where the terrane alone faces it.
    plate, inputs = _mixed_contact("oceanic", [inherit, c, c], [2, 2, 3])
    np.testing.assert_array_equal(torque.down_going_side(plate, inputs, reach), [True, True, False])
    # Single-type plates: as before, all or nothing.
    plate, inputs = _mixed_contact("oceanic", [inherit] * 3, [2, 2, 3])
    assert torque.down_going_side(plate, inputs, reach).all()
    plate, inputs = _mixed_contact("continental", [inherit] * 3, [2, 2, 3])
    assert not torque.down_going_side(plate, inputs, reach).any()
    assert not torque.subducting_boundary_mask(plate, inputs, reach).any()
    # A diverging contact holds no slab, so the continent behind it isn't dragged down.
    plate, inputs = _mixed_contact("continental", [o, inherit, inherit], [2, 2, 3], rate_cm_per_yr=-3.0)
    np.testing.assert_array_equal(torque.down_going_side(plate, inputs, reach), [True, False, False])


def test_a_continent_behind_its_own_subducting_ocean_floor_jams_without_a_polarized_front():
    # India (#447): before, only a polarized front could put a continental-majority plate's
    # continent on the down side, so it didn't jam until the prepass had seen the front.
    reach = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * line_spacing_rad(1.0)
    plate, inputs = _mixed_contact("continental", [CRUST_TYPE_OCEANIC, 0, 0], [2, 2, 3])
    inputs.own_hc = np.full(3, 35_000.0)
    plate.plate_id, plate.cell_keys = 1, np.arange(3)
    np.testing.assert_array_equal(torque.contact_jam_factor(plate, None, [], inputs, reach), [0.0, 1.0, 0.0])


def test_contact_drag_acts_on_continent_continent_nodes_whatever_either_plates_majority_type():
    plate, inputs = _mixed_contact("oceanic", [0, CRUST_TYPE_CONTINENTAL, CRUST_TYPE_CONTINENTAL], [2, 2, 2])
    inputs.neighbor_is_oceanic = np.array([False, False, True])
    contested = np.ones(3, dtype=bool)
    np.testing.assert_array_equal(torque.continental_contact_mask(plate, inputs, contested), [False, True, False])
    plate.crust_type = "continental"
    inputs.own_crust_type_codes = np.zeros(3, dtype=np.int8)
    np.testing.assert_array_equal(torque.continental_contact_mask(plate, inputs, contested), [True, True, False])


def test_the_neighbour_side_of_a_contact_is_its_nearest_nodes_own_crust():
    """#447: an ocean margin on a continental-majority neighbour reads as oceanic, and a
    continent on an oceanic-majority neighbour as continental."""
    from scipy.spatial import cKDTree

    class Neighbour:
        def __init__(self, plate_id, crust_type, points, codes):
            self.plate_id, self.crust_type, self.omega = plate_id, crust_type, np.zeros(3)
            self._tree, self._codes = cKDTree(points), np.asarray(codes, dtype=np.int8)

        def get_node_kdtree(self):
            return self._tree

        def collect(self, name):
            assert name == "crust_type_code"
            return self._codes

    class Own:
        crust_type, omega = "continental", np.zeros(3)

        def __init__(self, points):
            self._points = points

        def all_points_and_elevation(self):
            return self._points, np.zeros(len(self._points))

        def collect(self, name):
            return np.zeros(len(self._points), dtype=np.int8 if name == "crust_type_code" else float)

    spacing = line_spacing_rad(1.0)
    lat = np.array([0.0, 10.0, 20.0, 30.0]) * spacing
    own = np.column_stack([np.cos(lat), np.zeros(4), np.sin(lat)])
    across = np.column_stack([np.cos(lat) * np.cos(spacing), np.cos(lat) * np.sin(spacing), np.sin(lat)])
    land = Neighbour(2, "continental", across[:2], [CRUST_TYPE_OCEANIC, 0])
    sea = Neighbour(3, "oceanic", across[2:], [CRUST_TYPE_CONTINENTAL, 0])
    inputs = torque.gather_boundary_force_inputs(Own(own), [land, sea], spacing, 3 * spacing)
    np.testing.assert_array_equal(inputs.neighbor_plate_id, [2, 2, 3, 3])
    np.testing.assert_array_equal(inputs.neighbor_is_oceanic, [True, False, False, True])


def test_both_sides_of_a_jammed_contact_feel_the_lower_sides_jam():
    """One resistance across the contact, equal and opposite (#400): the upper plate is never in
    the jam mask, but reads j at the lower plate's nearest node; the lower plate reads its own."""
    spacing = line_spacing_rad(1.0)
    reach = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing
    # Upper: a still continental plate at the contact, its nodes facing a continental front.
    upper, inputs, _ = _contact_band(spacing, 1)
    upper.plate_id = 1
    n = len(inputs.own_points)
    inputs.neighbor_plate_id = np.full(n, 2)
    inputs.neighbor_node_index = np.arange(n)
    inputs.neighbor_omega = np.tile([0.0, 0.0, -mantle.cm_per_yr_to_rad_per_yr(2.0)], (n, 1))  # neighbour closing on it
    lower = _JamPlate(2, "continental", [35_000.0, 22_500.0, 35_000.0], [CRUST_TYPE_CONTINENTAL] * 3)
    inputs.neighbor_is_oceanic = np.zeros(n, dtype=bool)

    np.testing.assert_allclose(torque.contact_jam_factor(upper, _frame(p2=[0, 1, 2]), [lower], inputs, reach), [1.0, 0.5, 1.0])
    # No polarized front and no slab of the neighbour's at the contact, so nothing is forced
    # down: a plain continental contact doesn't jam.
    assert np.all(torque.contact_jam_factor(upper, None, [lower], inputs, reach) == 0.0)
    # The neighbour's ocean floor reaching the contact is a slab: its continent behind it is
    # forced down and jams, front or not (#447); the ocean floor itself doesn't.
    lower.fields["crust_type_code"][2] = CRUST_TYPE_OCEANIC
    inputs.neighbor_is_oceanic[2] = True
    np.testing.assert_allclose(torque.contact_jam_factor(upper, None, [lower], inputs, reach), [1.0, 0.5, 0.0])
