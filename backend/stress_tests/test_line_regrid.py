import numpy as np
from app import line_regrid
from app.plates import TARGET_LINE_SPACING_RAD, ElevationLine
from app.world import generate_world, step_world


def test_regularization_runs_periodically_during_stepping():
    world = generate_world(seed=30, num_plates=8)
    assert world.steps_since_regularize == 0
    for i in range(line_regrid.REGULARIZE_INTERVAL_STEPS - 1):
        step_world(world, years=1_000_000)
        assert world.steps_since_regularize == i + 1
    step_world(world, years=1_000_000)
    assert world.steps_since_regularize == 0  # just ran regularization and reset


def test_lines_stay_well_formed_after_many_steps_with_regularization():
    world = generate_world(seed=31, num_plates=8)
    for _ in range(25):
        step_world(world, years=2_000_000)
    for plate in world.plates:
        for line in plate.lines:
            assert np.all(np.diff(line.theta) > 0)
            assert len(line.theta) == len(line.elevation)
