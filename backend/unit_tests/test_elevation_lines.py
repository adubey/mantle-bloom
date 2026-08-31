import numpy as np
from app.elevation_lines import (
    TARGET_LINE_SPACING_RAD,
    ElevationLine,
    _crumple_elevation,
    iter_local_lattice,
    line_spacing_rad,
    needs_regularizing,
    regularize_line,
)


def test_iter_local_lattice_default_spacing_matches_target():
    rows = list(iter_local_lattice(np.eye(3)))
    assert len(rows) > 0
    phis = [phi for phi, _, _ in rows]
    assert np.allclose(np.diff(sorted(phis)), TARGET_LINE_SPACING_RAD, atol=1e-9)


def test_iter_local_lattice_custom_spacing_changes_row_count():
    coarse = list(iter_local_lattice(np.eye(3), spacing_rad=TARGET_LINE_SPACING_RAD * 4))
    fine = list(iter_local_lattice(np.eye(3), spacing_rad=TARGET_LINE_SPACING_RAD / 2))
    assert len(fine) > len(coarse)
    total_fine = sum(len(theta) for _, theta, _ in fine)
    total_coarse = sum(len(theta) for _, theta, _ in coarse)
    assert total_fine > total_coarse


def test_line_spacing_rad_matches_default_at_density_one():
    assert line_spacing_rad(1.0) == TARGET_LINE_SPACING_RAD


def test_line_spacing_rad_halves_at_4x_density():
    # Node count scales with the square of resolution, so 4x the nodes needs the spacing
    # *halved*, not quartered.
    assert np.isclose(line_spacing_rad(4.0), TARGET_LINE_SPACING_RAD / 2)


def test_elevation_line_soil_and_resource_fields_default_to_zero():
    theta = np.array([0.0, 0.1, 0.2])
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.zeros(3))
    for field in ("soil_depth", "soil_mineral_content", "soil_organic_content", "coal_deposit_m", "oil_gas_deposit_m", "mineral_deposit_m"):
        values = getattr(line, field)
        assert values.shape == theta.shape
        assert np.all(values == 0.0)


def test_needs_regularizing_false_for_evenly_spaced_line():
    dtheta = TARGET_LINE_SPACING_RAD
    theta = np.arange(0.0, 1.0, dtheta)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.zeros_like(theta))
    assert not needs_regularizing(line)


def test_needs_regularizing_true_for_uneven_line():
    theta = np.array([0.0, 0.001, 0.5, 0.9])  # wildly irregular gaps
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.zeros_like(theta))
    assert needs_regularizing(line)


def test_regularize_line_preserves_endpoints():
    theta = np.array([0.0, 0.001, 0.5, 0.9])
    elevation = np.array([100.0, 110.0, -200.0, 50.0])
    line = ElevationLine(phi=0.0, theta=theta, elevation=elevation)

    regularized = regularize_line(line)
    assert np.isclose(regularized.theta[0], theta[0])
    assert np.isclose(regularized.theta[-1], theta[-1])
    assert np.all(np.diff(regularized.theta) > 0)


def test_regularize_line_interpolates_elevation_reasonably():
    theta = np.array([0.0, 0.2, 0.4])
    elevation = np.array([0.0, 100.0, 0.0])
    line = ElevationLine(phi=0.0, theta=theta, elevation=elevation)

    regularized = regularize_line(line)
    # New nodes' elevations must stay within the original data's min/max (no overshoot
    # from linear interpolation).
    assert regularized.elevation.min() >= elevation.min() - 1e-9
    assert regularized.elevation.max() <= elevation.max() + 1e-9


def test_regularize_line_short_line_is_a_no_op():
    theta = np.array([0.0, 0.3])
    elevation = np.array([1.0, 2.0])
    line = ElevationLine(phi=0.0, theta=theta, elevation=elevation)
    result = regularize_line(line)
    assert result is line


def test_regularize_line_too_dense_crumples_instead_of_thinning_out():
    # 30 points packed into a span target spacing would only need ~6 for -- points are too
    # close, the "vulcanism increased density" case, not the "gaps too wide" one.
    theta = np.linspace(0.0, 0.03, 30)
    elevation = 50.0 * np.sin(theta * 300)
    line = ElevationLine(phi=0.0, theta=theta, elevation=elevation)

    regularized = regularize_line(line, spacing_rad=0.005)
    assert len(regularized.theta) < len(theta)
    assert np.isclose(regularized.theta[0], theta[0])
    assert np.isclose(regularized.theta[-1], theta[-1])
    assert regularized.elevation[0] == elevation[0]
    assert regularized.elevation[-1] == elevation[-1]


def test_crumple_elevation_amplifies_peaks_and_valleys():
    n = 40
    x = np.linspace(0.0, 1.0, n)
    elevation = 100.0 * np.sin(x * 6 * np.pi)

    crumpled = _crumple_elevation(elevation, 8)
    assert len(crumpled) == 8
    # Squashing the same shape into fewer points should exaggerate its vertical range, not
    # just resample it -- a plain linear thin-out would stay within the original min/max.
    assert crumpled.max() > elevation.max()
    assert crumpled.min() < elevation.min()


def test_crumple_elevation_preserves_endpoints_exactly():
    elevation = np.array([10.0, 40.0, -20.0, 5.0, 30.0, 0.0, -15.0, 25.0])
    crumpled = _crumple_elevation(elevation, 4)
    assert crumpled[0] == elevation[0]
    assert crumpled[-1] == elevation[-1]


def _wound_ring(phi: float, revolutions: float, spacing_rad: float) -> ElevationLine:
    """A single over-wound row: theta marching well past a full 2*pi (as a plate that grew
    around its own local pole did before the wrap guard existed / on worlds saved then)."""
    dtheta = spacing_rad / max(np.cos(phi), 1e-3)
    theta = np.arange(0.0, revolutions * 2.0 * np.pi, dtheta)
    return ElevationLine(phi=phi, theta=theta, elevation=np.linspace(-1000.0, 1000.0, len(theta)))


def test_needs_regularizing_flags_an_over_wound_ring():
    spacing = line_spacing_rad(1.0)
    assert needs_regularizing(_wound_ring(1.45, revolutions=3.4, spacing_rad=spacing), spacing)
    # An ordinary near-full but sub-revolution ring is fine.
    assert not needs_regularizing(_wound_ring(1.45, revolutions=0.95, spacing_rad=spacing), spacing)


def test_regularize_line_unwinds_an_over_wound_ring_to_one_revolution():
    spacing = line_spacing_rad(1.0)
    wound = _wound_ring(1.45, revolutions=4.0, spacing_rad=spacing)
    fixed = regularize_line(wound, spacing)

    span = fixed.theta[-1] - fixed.theta[0]
    assert span <= 2.0 * np.pi + 1e-9
    assert span > 2.0 * np.pi - 3 * (spacing / np.cos(1.45))  # kept a full revolution, not a sliver
    assert np.all(np.diff(fixed.theta) > 0)
    assert len(fixed) < len(wound)  # the inner windings are gone
    # It kept the *outermost* revolution (theta near the wound row's own high end).
    assert fixed.theta[-1] == wound.theta[-1]
