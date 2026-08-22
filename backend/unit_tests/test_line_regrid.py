import numpy as np
from app import line_regrid
from app.plates import TARGET_LINE_SPACING_RAD, ElevationLine
from app.world import generate_world, step_world


def test_needs_regularizing_false_for_evenly_spaced_line():
    dtheta = TARGET_LINE_SPACING_RAD
    theta = np.arange(0.0, 1.0, dtheta)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.zeros_like(theta))
    assert not line_regrid.needs_regularizing(line)


def test_needs_regularizing_true_for_uneven_line():
    theta = np.array([0.0, 0.001, 0.5, 0.9])  # wildly irregular gaps
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.zeros_like(theta))
    assert line_regrid.needs_regularizing(line)


def test_regularize_line_preserves_endpoints():
    theta = np.array([0.0, 0.001, 0.5, 0.9])
    elevation = np.array([100.0, 110.0, -200.0, 50.0])
    line = ElevationLine(phi=0.0, theta=theta, elevation=elevation)

    regularized = line_regrid.regularize_line(line)
    assert np.isclose(regularized.theta[0], theta[0])
    assert np.isclose(regularized.theta[-1], theta[-1])
    assert np.all(np.diff(regularized.theta) > 0)


def test_regularize_line_interpolates_elevation_reasonably():
    theta = np.array([0.0, 0.2, 0.4])
    elevation = np.array([0.0, 100.0, 0.0])
    line = ElevationLine(phi=0.0, theta=theta, elevation=elevation)

    regularized = line_regrid.regularize_line(line)
    # New nodes' elevations must stay within the original data's min/max (no overshoot
    # from linear interpolation).
    assert regularized.elevation.min() >= elevation.min() - 1e-9
    assert regularized.elevation.max() <= elevation.max() + 1e-9


def test_regularize_line_short_line_is_a_no_op():
    theta = np.array([0.0, 0.3])
    elevation = np.array([1.0, 2.0])
    line = ElevationLine(phi=0.0, theta=theta, elevation=elevation)
    result = line_regrid.regularize_line(line)
    assert result is line


def test_regularize_line_too_dense_crumples_instead_of_thinning_out():
    # 30 points packed into a span target spacing would only need ~6 for -- points are too
    # close, the "vulcanism increased density" case, not the "gaps too wide" one.
    theta = np.linspace(0.0, 0.03, 30)
    elevation = 50.0 * np.sin(theta * 300)
    line = ElevationLine(phi=0.0, theta=theta, elevation=elevation)

    regularized = line_regrid.regularize_line(line, spacing_rad=0.005)
    assert len(regularized.theta) < len(theta)
    assert np.isclose(regularized.theta[0], theta[0])
    assert np.isclose(regularized.theta[-1], theta[-1])
    assert regularized.elevation[0] == elevation[0]
    assert regularized.elevation[-1] == elevation[-1]


def test_crumple_elevation_amplifies_peaks_and_valleys():
    n = 40
    x = np.linspace(0.0, 1.0, n)
    elevation = 100.0 * np.sin(x * 6 * np.pi)

    crumpled = line_regrid._crumple_elevation(elevation, 8)
    assert len(crumpled) == 8
    # Squashing the same shape into fewer points should exaggerate its vertical range, not
    # just resample it -- a plain linear thin-out would stay within the original min/max.
    assert crumpled.max() > elevation.max()
    assert crumpled.min() < elevation.min()


def test_crumple_elevation_preserves_endpoints_exactly():
    elevation = np.array([10.0, 40.0, -20.0, 5.0, 30.0, 0.0, -15.0, 25.0])
    crumpled = line_regrid._crumple_elevation(elevation, 4)
    assert crumpled[0] == elevation[0]
    assert crumpled[-1] == elevation[-1]
