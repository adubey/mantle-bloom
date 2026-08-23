import numpy as np
from app import geometry, reassign
from app.plates import ElevationLine, Plate
from app.world import World, generate_world, step_world


def _plate_with_filler(plate_id, seed_xyz, lines):
    """A plate with the given lines plus a placeholder line, purely so the plate always has
    more than one line -- matching test_merge_split.py's own pattern (a single-line plate is
    pruned as "no land left" by merge_split.remove_defunct_plates before a test's own logic
    would ever run). The filler is 6 tightly-packed nodes (0.01 spacing) rather than 2 or 4
    spread-out ones: every node's 3 nearest same-line siblings then sit far closer (~0.03)
    than anything else in a test's geometry could possibly be, so the filler is always
    "obviously its own plate" regardless of exactly where the other plate's own lines land --
    no need to fine-tune its position relative to a specific test's boundary geometry (an
    earlier version placed the filler by copying its anchor line's own theta values, which
    for the misplaced-point test below accidentally recreated a *second* copy of the exact
    boundary condition under test, since that test's own two lines are constructed to be
    directly adjacent)."""
    frame = geometry.plate_frame_from_seed(seed_xyz)
    anchor = lines[0]
    base_theta = anchor.theta[0] - 1.0
    filler_theta = base_theta + 0.01 * np.arange(6)
    filler = ElevationLine(phi=anchor.phi + 0.15, theta=filler_theta, elevation=np.zeros(6))
    return Plate(plate_id=plate_id, frame=frame, crust_type="oceanic", lines=[*lines, filler])


def test_reassign_runs_periodically_and_never_on_a_regularize_step():
    from app import line_regrid

    # node_density=0.5 (the coarsest choice, an eighth of the default 4.0, see
    # plates.NODE_DENSITY_CHOICES) -- this test only checks the two cadence counters, not node
    # positions.
    world = generate_world(seed=40, num_plates=8, node_density=0.5)
    # This test only checks the two cadence counters, never climate/erosion/hydrology output
    # (see World.simulate_climate_biomes) -- skipping that per-step computation speeds this up
    # without changing what it exercises.
    world.simulate_climate_biomes = False
    assert line_regrid.REGULARIZE_INTERVAL_STEPS == reassign.REASSIGN_INTERVAL_STEPS  # same period by default
    for _ in range(30):
        step_world(world, years=1_000_000)
        # The two cadence counters must never both be reset (i.e. both passes run) on the
        # same step -- the whole point of staggering them (see world.step_world).
        assert not (world.steps_since_regularize == 0 and world.steps_since_reassign == 0)
