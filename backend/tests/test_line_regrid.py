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


def test_garbage_collection_runs_periodically_during_stepping():
    world = generate_world(seed=30, num_plates=8)
    assert world.steps_since_gc == 0
    for i in range(line_regrid.GC_INTERVAL_STEPS - 1):
        step_world(world, years=1_000_000)
        assert world.steps_since_gc == i + 1
    step_world(world, years=1_000_000)
    assert world.steps_since_gc == 0  # just ran GC and reset


def test_lines_stay_well_formed_after_many_steps_with_gc():
    world = generate_world(seed=31, num_plates=8)
    for _ in range(25):
        step_world(world, years=2_000_000)
    for plate in world.plates:
        for line in plate.lines:
            assert np.all(np.diff(line.theta) > 0)
            assert len(line.theta) == len(line.elevation)
