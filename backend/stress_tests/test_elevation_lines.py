import numpy as np
from app import elevation_lines
from app.world import generate_world, step_world


def test_regularization_runs_every_step_not_periodically():
    # node_density=0.5 (the coarsest choice, an eighth of the default 4.0, see
    # plates.NODE_DENSITY_CHOICES) -- this test only checks line spacing, not node positions.
    # Unlike the old periodic (every REGULARIZE_INTERVAL_STEPS) cadence, LithospherePlate.deform
    # now calls elevation_lines.regularize_line itself at the end of every single call -- so
    # no line should ever need regularizing right after a step, not just every 5th one.
    world = generate_world(seed=30, num_plates=8, node_density=0.5)
    world.simulate_climate_biomes = False
    spacing_rad = elevation_lines.line_spacing_rad(world.node_density)
    for _ in range(3):
        step_world(world, years=1_000_000)
        for plate in world.plates:
            for line in plate.lines:
                assert not elevation_lines.needs_regularizing(line, spacing_rad)


def test_lines_stay_well_formed_after_many_steps_with_regularization():
    world = generate_world(seed=31, num_plates=8, node_density=0.5)  # see the test above for why this is safe
    world.simulate_climate_biomes = False  # see the test above for why this is safe here
    for _ in range(25):
        step_world(world, years=2_000_000)
    for plate in world.plates:
        for line in plate.lines:
            assert np.all(np.diff(line.theta) > 0)
            assert len(line.theta) == len(line.elevation)
