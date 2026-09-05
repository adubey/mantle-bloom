"""Angle-aware rift stretching: closing a divergent gap by stretching existing crust
(`rheology.stretch_components`/`apply_stretch_thinning`) rather than only ever growing
brand-new nodes -- see `LithospherePlate._grow_or_shrink_line_for_deform`'s end-stretch (the
theta/within-row share of a gap) and `_claim_adjacent_territory`'s mass-conserving new-row
claim (the phi/between-row share). `LithospherePlate._separation_components` is what decides,
per node, how much of a gap falls into each share -- from the nearest active fault's own
tangent in "fault" mode (`faults.fault_tangent_components`), or from the already-computed
`direction_to_neighbor` projected into the node's own local tangent basis
(`geometry.local_separation_components`) otherwise.
"""

import types

import numpy as np
import pytest

from app import faults, geometry, rheology
from app.elevation_lines import ElevationLine, line_spacing_rad, regularize_line
from app.lithosphere import MIN_CRUSTAL_THICKNESS_M, MIN_MANTLE_LITHOSPHERE_THICKNESS_M
from app.lithosphere_plate import K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION, new_plate

ANGLES_DEG = [0, 15, 30, 45, 60, 75, 90, 120, 150, 180, 210, 240, 270, 300, 330, 350]


# --------------------------------------------------------------------------------------------
# rheology.stretch_components -- pure vector decomposition, the primitive both call sites share.
# --------------------------------------------------------------------------------------------


def test_stretch_components_pure_theta_and_pure_phi():
    theta_gap, phi_gap = rheology.stretch_components(sep_theta=1.0, sep_phi=0.0, gap_rad=10.0)
    assert theta_gap == pytest.approx(10.0)
    assert phi_gap == pytest.approx(0.0, abs=1e-9)

    theta_gap, phi_gap = rheology.stretch_components(sep_theta=0.0, sep_phi=1.0, gap_rad=10.0)
    assert theta_gap == pytest.approx(0.0, abs=1e-9)
    assert phi_gap == pytest.approx(10.0)


@pytest.mark.parametrize("angle_deg", ANGLES_DEG)
def test_stretch_components_at_many_angles(angle_deg):
    """A `gap_rad` of stretch, in a direction `angle_deg` from the theta axis, decomposes into
    orthogonal (theta, phi) components -- a plain vector decomposition, so the two components'
    squares must always sum back to the total gap squared, regardless of angle."""
    angle = np.radians(angle_deg)
    sep_theta, sep_phi = np.cos(angle), np.sin(angle)
    gap = 10.0

    theta_gap, phi_gap = rheology.stretch_components(sep_theta, sep_phi, gap)

    assert theta_gap >= 0.0
    assert phi_gap >= 0.0
    assert theta_gap == pytest.approx(gap * abs(np.cos(angle)), abs=1e-9)
    assert phi_gap == pytest.approx(gap * abs(np.sin(angle)), abs=1e-9)
    assert theta_gap**2 + phi_gap**2 == pytest.approx(gap**2)


def test_stretch_components_normalizes_a_non_unit_separation_vector():
    theta_gap, phi_gap = rheology.stretch_components(sep_theta=3.0, sep_phi=4.0, gap_rad=5.0)
    assert theta_gap == pytest.approx(3.0)
    assert phi_gap == pytest.approx(4.0)


def test_stretch_components_handles_a_zero_separation_vector():
    theta_gap, phi_gap = rheology.stretch_components(sep_theta=0.0, sep_phi=0.0, gap_rad=5.0)
    assert np.isfinite(theta_gap) and np.isfinite(phi_gap)


# --------------------------------------------------------------------------------------------
# rheology.apply_stretch_thinning -- mass-conserving thinning + the melting/volcano trigger.
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("stretch_fraction", [0.0, 0.1, 0.5, 1.0, 2.0, 5.0])
def test_apply_stretch_thinning_conserves_mass(stretch_fraction):
    hc = np.array([20_000.0])
    hm = np.array([50_000.0])
    old_spacing = np.array([0.01])
    extra_gap = np.array([0.01 * stretch_fraction])

    new_hc, new_hm, _ = rheology.apply_stretch_thinning(hc, hm, old_spacing, extra_gap)

    new_spacing = old_spacing[0] + extra_gap[0]
    # Thickness x footprint is unchanged -- the wider gap between this node and its inward
    # neighbour *is* the lower Hc, not an approximation of it.
    assert new_hc[0] * new_spacing == pytest.approx(hc[0] * old_spacing[0])
    assert new_hm[0] * new_spacing == pytest.approx(hm[0] * old_spacing[0])


def test_apply_stretch_thinning_no_stretch_is_a_no_op():
    hc, hm = np.array([20_000.0]), np.array([50_000.0])
    new_hc, new_hm, melting = rheology.apply_stretch_thinning(hc, hm, np.array([0.01]), np.array([0.0]))
    assert new_hc[0] == pytest.approx(hc[0])
    assert new_hm[0] == pytest.approx(hm[0])
    assert not melting[0]


def test_apply_stretch_thinning_floors_at_minimum_thickness():
    hc, hm = np.array([600.0]), np.array([3_000.0])
    new_hc, new_hm, _ = rheology.apply_stretch_thinning(hc, hm, np.array([0.01]), np.array([100.0]))
    assert new_hc[0] == pytest.approx(MIN_CRUSTAL_THICKNESS_M)
    assert new_hm[0] == pytest.approx(MIN_MANTLE_LITHOSPHERE_THICKNESS_M)


def test_apply_stretch_thinning_erupts_exactly_when_crossing_the_rift_threshold():
    threshold = rheology.RIFT_CRITICAL_THICKNESS_M
    hc = np.array([threshold * 1.05])
    hm = np.array([40_000.0])
    old_spacing = np.array([0.01])

    # A tiny stretch keeps Hc above the threshold -- no eruption.
    _, _, melting_above = rheology.apply_stretch_thinning(hc, hm, old_spacing, np.array([old_spacing[0] * 0.01]))
    assert not melting_above[0]

    # A stretch large enough to push Hc below the threshold -- eruption fires.
    _, _, melting_below = rheology.apply_stretch_thinning(hc, hm, old_spacing, np.array([old_spacing[0] * 10.0]))
    assert melting_below[0]


# --------------------------------------------------------------------------------------------
# geometry.local_separation_components -- projecting a world-frame direction (e.g.
# direction_to_neighbor) into a node's own local (theta, phi) tangent basis.
# --------------------------------------------------------------------------------------------


def _tangent_basis(phi: float, theta: float) -> tuple[np.ndarray, np.ndarray]:
    sin_p, cos_p = np.sin(phi), np.cos(phi)
    sin_t, cos_t = np.sin(theta), np.cos(theta)
    theta_hat = np.array([-sin_t, cos_t, 0.0])
    phi_hat = np.array([-sin_p * cos_t, -sin_p * sin_t, cos_p])
    return theta_hat, phi_hat


@pytest.mark.parametrize("phi_deg,theta_deg", [(0, 0), (10, 20), (-30, 50), (60, -100), (45, 170)])
def test_local_separation_components_recovers_pure_axes(phi_deg, theta_deg):
    phi, theta = np.radians(phi_deg), np.radians(theta_deg)
    frame = np.eye(3)
    theta_hat, phi_hat = _tangent_basis(phi, theta)

    sep_theta, sep_phi = geometry.local_separation_components(frame, phi, theta, theta_hat)
    assert sep_theta == pytest.approx(1.0)
    assert sep_phi == pytest.approx(0.0, abs=1e-9)

    sep_theta, sep_phi = geometry.local_separation_components(frame, phi, theta, phi_hat)
    assert sep_theta == pytest.approx(0.0, abs=1e-9)
    assert sep_phi == pytest.approx(1.0)


@pytest.mark.parametrize("angle_deg", ANGLES_DEG)
def test_local_separation_components_at_many_angles(angle_deg):
    phi, theta = 0.3, 0.7
    frame = np.eye(3)
    theta_hat, phi_hat = _tangent_basis(phi, theta)
    angle = np.radians(angle_deg)
    direction = np.cos(angle) * theta_hat + np.sin(angle) * phi_hat

    sep_theta, sep_phi = geometry.local_separation_components(frame, phi, theta, direction)

    assert sep_theta == pytest.approx(np.cos(angle), abs=1e-9)
    assert sep_phi == pytest.approx(np.sin(angle), abs=1e-9)


def test_local_separation_components_transforms_out_of_a_rotated_frame():
    """The projection must undo the plate's own rotation first -- a world-frame direction
    means nothing in local (theta, phi) terms until it is."""
    phi, theta = 0.2, -0.5
    rng = np.random.default_rng(0)
    axis = geometry.normalize(rng.normal(size=3))
    frame = geometry.rotation_matrix(axis, 0.9)
    theta_hat, _ = _tangent_basis(phi, theta)
    world_direction = geometry.to_world(frame, theta_hat)

    sep_theta, sep_phi = geometry.local_separation_components(frame, phi, theta, world_direction)

    assert sep_theta == pytest.approx(1.0)
    assert sep_phi == pytest.approx(0.0, abs=1e-9)


# --------------------------------------------------------------------------------------------
# faults.fault_tangent_components -- a fault trace's own tangent, swapped (separation runs
# *across* a fault's own strike, not along it).
# --------------------------------------------------------------------------------------------


def _straight_fault(plate_id: int, phi0: float, theta0: float, angle_rad: float) -> "faults.Fault":
    """A short, straight fault trace through (phi0, theta0) at `angle_rad` from the theta
    axis -- physically corrected by cos(phi0) so the trace's own tangent, in real arc-length
    terms, actually points at `angle_rad` (matching `fault_tangent_components`'s own
    cos(phi)-scaling of local_theta)."""
    t = np.linspace(-0.01, 0.01, 5)
    cos_p = max(np.cos(phi0), 1e-3)
    local_theta = theta0 + t * np.cos(angle_rad) / cos_p
    local_phi = phi0 + t * np.sin(angle_rad)
    return faults.Fault(
        fault_id=1,
        plate_id=plate_id,
        kind="normal",
        local_phi=local_phi,
        local_theta=local_theta,
        slip_rate_m_per_myr=1.0,
        dip_deg=60.0,
        strike_sense=1,
        dip_dir_local=np.zeros(3),
        lifespan_myr=5.0,
        birth_years=0.0,
        birth_distance_from_boundary_km=0.0,
    )


@pytest.mark.parametrize("angle_deg", ANGLES_DEG)
def test_fault_tangent_components_swaps_the_faults_own_strike(angle_deg):
    angle = np.radians(angle_deg)
    phi0, theta0 = 0.1, 0.2
    fault = _straight_fault(plate_id=7, phi0=phi0, theta0=theta0, angle_rad=angle)
    plate = types.SimpleNamespace(plate_id=7)
    world = types.SimpleNamespace(faults=[fault], boundary_faults=[])

    result = faults.fault_tangent_components(world, plate, phi0, theta0)

    assert result is not None
    sep_theta, sep_phi = result
    norm = np.hypot(sep_theta, sep_phi)
    assert norm > 0.0
    sep_theta, sep_phi = abs(sep_theta) / norm, abs(sep_phi) / norm
    # Swapped: a fault whose own tangent runs along theta (angle=0) separates along phi, and
    # vice versa.
    assert sep_theta == pytest.approx(abs(np.sin(angle)), abs=1e-6)
    assert sep_phi == pytest.approx(abs(np.cos(angle)), abs=1e-6)


def test_fault_tangent_components_none_when_plate_has_no_active_fault():
    plate = types.SimpleNamespace(plate_id=7)
    world = types.SimpleNamespace(faults=[], boundary_faults=[])
    assert faults.fault_tangent_components(world, plate, 0.1, 0.2) is None


def test_fault_tangent_components_none_when_only_other_plates_have_faults():
    fault = _straight_fault(plate_id=99, phi0=0.1, theta0=0.2, angle_rad=0.0)
    plate = types.SimpleNamespace(plate_id=7)
    world = types.SimpleNamespace(faults=[fault], boundary_faults=[])
    assert faults.fault_tangent_components(world, plate, 0.1, 0.2) is None


# --------------------------------------------------------------------------------------------
# Regularize carries the thinning through: a widened end gap, resampled back to target
# density, must keep the thinned Hc profile rather than reverting to a flat reference value.
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("angle_deg", [0, 30, 45, 60, 90])
def test_regularize_line_carries_stretch_thinned_hc_through(angle_deg):
    node_density = 1.0
    spacing_rad = line_spacing_rad(node_density)
    phi = 0.2
    dtheta = spacing_rad / max(np.cos(phi), 1e-3)
    n = 20
    theta = phi * 0 + dtheta * np.arange(n)  # 0, dtheta, 2*dtheta, ...
    reference_hc = 11_000.0
    hc = np.full(n, reference_hc)
    hm = np.full(n, 40_000.0)

    # Stretch the last node outward by this angle's own theta-attributable share of one full
    # gap, thinning its own Hc/Hm to conserve mass -- exactly what
    # LithospherePlate._grow_or_shrink_line_for_deform's end-stretch leaves behind.
    angle = np.radians(angle_deg)
    gap = spacing_rad * 3.0
    theta_gap, _ = rheology.stretch_components(np.cos(angle), np.sin(angle), gap)
    new_hc_end, new_hm_end, _ = rheology.apply_stretch_thinning(
        hc[-1:], hm[-1:], np.array([dtheta]), np.array([theta_gap])
    )
    theta[-1] = theta[-1] + theta_gap
    hc[-1] = new_hc_end[0]
    hm[-1] = new_hm_end[0]
    elevation = np.zeros(n)

    line = ElevationLine(phi=phi, theta=theta, elevation=elevation, crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm)

    resampled = regularize_line(line, spacing_rad=spacing_rad)

    # The resampled line's Hc, read back at the stretched end's own theta position, must still
    # reflect the thinning -- not have snapped back to the flat reference value every other
    # node started at.
    hc_at_end = float(np.interp(theta[-1], resampled.theta, resampled.crustal_thickness_m))
    if theta_gap > 1e-9:
        assert hc_at_end < reference_hc * 0.999
        assert hc_at_end == pytest.approx(hc[-1], rel=0.15)
    else:
        assert hc_at_end == pytest.approx(reference_hc, rel=1e-6)


# --------------------------------------------------------------------------------------------
# LithospherePlate._separation_components -- the mode-aware entry point both
# _grow_or_shrink_line_for_deform and _claim_adjacent_territory call.
# --------------------------------------------------------------------------------------------


def _small_oceanic_disk(node_density=1.0, radius_rad=0.4, seed=1):
    frame = np.eye(3)
    spacing_rad = line_spacing_rad(node_density)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        return geometry.angular_distance(world_pts, frame[:, 2]) < radius_rad

    return new_plate(0, frame, "oceanic", spacing_rad, seed, is_owned=is_owned), spacing_rad


@pytest.mark.parametrize("angle_deg", ANGLES_DEG)
def test_separation_components_boundary_mode_matches_geometry_projection(angle_deg):
    plate, _ = _small_oceanic_disk()
    world = types.SimpleNamespace(fault_deformation_mode="boundary")
    phi, theta = 0.1, 0.2
    theta_hat, phi_hat = _tangent_basis(phi, theta)
    angle = np.radians(angle_deg)
    direction_local = np.cos(angle) * theta_hat + np.sin(angle) * phi_hat
    direction_world = geometry.to_world(plate.frame, direction_local)

    sep_theta, sep_phi = plate._separation_components(world, phi, theta, direction_world)

    assert sep_theta == pytest.approx(np.cos(angle), abs=1e-9)
    assert sep_phi == pytest.approx(np.sin(angle), abs=1e-9)


def test_separation_components_falls_back_to_pure_theta_with_no_direction_at_all():
    """A genuinely isolated end (torque.gather_boundary_force_inputs leaves
    direction_to_neighbor at zero when there is no neighbour within reach at all) must not
    silently stop growing -- it falls back to this row's own theta direction, matching ordinary
    end growth's long-standing behaviour before rift-stretch existed."""
    plate, _ = _small_oceanic_disk()
    world = types.SimpleNamespace(fault_deformation_mode="boundary")
    sep_theta, sep_phi = plate._separation_components(world, 0.1, 0.2, np.zeros(3))
    assert (sep_theta, sep_phi) == (1.0, 0.0)


def test_separation_components_fault_mode_prefers_the_fault_tangent_over_direction():
    plate, _ = _small_oceanic_disk()
    phi0, theta0 = 0.1, 0.2
    # A fault running purely along phi (angle=90deg): separation should be pure theta.
    fault = _straight_fault(plate_id=plate.plate_id, phi0=phi0, theta0=theta0, angle_rad=np.pi / 2)
    world = types.SimpleNamespace(fault_deformation_mode="fault", faults=[fault], boundary_faults=[])
    # A direction that, if it were used instead, would say the opposite (pure phi).
    _, misleading_phi_hat = _tangent_basis(phi0, theta0)
    misleading_direction = geometry.to_world(plate.frame, misleading_phi_hat)

    sep_theta, sep_phi = plate._separation_components(world, phi0, theta0, misleading_direction)

    norm = np.hypot(sep_theta, sep_phi)
    assert abs(sep_theta) / norm == pytest.approx(1.0, abs=1e-6)
    assert abs(sep_phi) / norm == pytest.approx(0.0, abs=1e-6)


# --------------------------------------------------------------------------------------------
# LithospherePlate._claim_adjacent_territory -- the phi/between-row share: gated on genuine
# phi-direction separation, and mass-conserving (new row + K existing neighbour rows) when it
# fires.
# --------------------------------------------------------------------------------------------


@pytest.mark.parametrize("angle_deg", [0, 20, 45, 70, 90])
def test_claim_adjacent_territory_gates_and_conserves_mass_by_angle(angle_deg, monkeypatch):
    plate, spacing_rad = _small_oceanic_disk(radius_rad=0.5)
    n_lines_before = len(plate.lines)
    ordered_before = sorted([ln for ln in plate.lines if len(ln) > 0], key=lambda ln: ln.phi)
    outer_rows_before = ordered_before[-K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION:]
    mean_hc_before = float(np.mean(np.concatenate([row.crustal_thickness_m for row in outer_rows_before])))

    angle = np.radians(angle_deg)
    monkeypatch.setattr(plate, "_separation_components", lambda world, phi, theta, direction: (np.cos(angle), np.sin(angle)))
    world = types.SimpleNamespace(fault_deformation_mode="boundary", seed=1, elapsed_years=0.0)

    plate._claim_adjacent_territory(world, neighbours=[], spacing_rad=spacing_rad)

    ordered_after = sorted([ln for ln in plate.lines if len(ln) > 0], key=lambda ln: ln.phi)

    if angle_deg == 0:
        # Pure theta-aligned separation: no genuine phi-direction gap, so this row-claim
        # mechanism should not fire at all (the theta-stretch end-growth owns this share).
        assert len(plate.lines) == n_lines_before
        return

    # Some angle with a real phi component: a new outermost row appears...
    assert len(ordered_after) == len(ordered_before) + (2 if angle_deg == 90 else 1) or len(ordered_after) > len(ordered_before)
    new_outer_rows = sorted([ln for ln in plate.lines if len(ln) > 0], key=lambda ln: ln.phi)[-K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION - 1 :]
    mean_hc_after = float(np.mean(np.concatenate([row.crustal_thickness_m for row in new_outer_rows])))
    # ...seeded thin, not at the full reference column "for free" -- the K existing rows plus
    # the new row share the reference volume between them (mass conservation), so the mean
    # across the K+1 rows involved must drop, not stay flat.
    assert mean_hc_after < mean_hc_before
