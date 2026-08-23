import numpy as np
from app import line_regrid
from app.plates import TARGET_LINE_SPACING_RAD, ElevationLine
from app.world import generate_world, step_world


def test_regularization_runs_periodically_during_stepping():
    # node_density=0.5 (the coarsest choice, an eighth of the default 4.0, see
    # plates.NODE_DENSITY_CHOICES) -- this test only checks the regularize cadence counter,
    # not node positions.
    world = generate_world(seed=30, num_plates=8, node_density=0.5)
    # Regularization cadence only depends on plate/node geometry, never on climate/erosion/
    # hydrology output (see World.simulate_climate_biomes) -- skipping that per-step
    # computation speeds this up without changing what it exercises.
    world.simulate_climate_biomes = False
    assert world.steps_since_regularize == 0
    for i in range(line_regrid.REGULARIZE_INTERVAL_STEPS - 1):
        step_world(world, years=1_000_000)
        assert world.steps_since_regularize == i + 1
    step_world(world, years=1_000_000)
    assert world.steps_since_regularize == 0  # just ran regularization and reset


def test_lines_stay_well_formed_after_many_steps_with_regularization():
    world = generate_world(seed=31, num_plates=8, node_density=0.5)  # see the test above for why this is safe
    world.simulate_climate_biomes = False  # see the test above for why this is safe here
    for _ in range(25):
        step_world(world, years=2_000_000)
    for plate in world.plates:
        for line in plate.lines:
            assert np.all(np.diff(line.theta) > 0)
            assert len(line.theta) == len(line.elevation)
