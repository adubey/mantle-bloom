import numpy as np

from app import geometry, reassign
from app.plates import ElevationLine, Plate
from app.world import World, generate_world, step_world


def _plate_with_filler(plate_id, seed_xyz, lines):
    """A plate with the given lines plus a placeholder line, purely so the plate always has
    more than one line -- matching test_merge_split.py's own pattern (a single-line plate is
    pruned as "no land left" by merge_split.remove_defunct_plates before a test's own logic
    would ever run). The filler shares its `lines[0]` anchor's theta values, offset only in
    phi, so it's always geometrically right next to that plate's own real data (not the
    exact same point, which -- when two plates share a frame, as these tests' do -- would be
    an exact-distance tie that confuses nearest-neighbor self-exclusion; and not off in some
    arbitrary fixed direction either, which risked landing closer to the *other* plate's own
    cluster than to this one's, since a shared frame gives both plates the same coordinate
    system to be "close" or "far" in)."""
    frame = geometry.plate_frame_from_seed(seed_xyz)
    anchor = lines[0]
    filler = ElevationLine(phi=anchor.phi + 0.15, theta=anchor.theta.copy(), elevation=np.zeros(len(anchor.theta)))
    return Plate(plate_id=plate_id, frame=frame, crust_type="oceanic", lines=[*lines, filler])


def test_reassign_moves_a_point_surrounded_by_a_foreign_plate():
    # Plate 0 owns a line at phi=0 spanning theta in [0, 0.4]; plate 1 (identical frame, so
    # local coordinates line up 1:1) owns a parallel line at phi=0 whose theta continues
    # right where plate 0's leaves off. theta=0.41 is nominally on plate 0's line but sits
    # closer to plate 1's cluster of nodes than to any of its own -- exactly the "3 of 4
    # nearest neighbors belong to another plate" case.
    misplaced_theta = 0.41
    own_line = ElevationLine(
        phi=0.0, theta=np.array([0.0, 0.1, 0.2, misplaced_theta]), elevation=np.array([10.0, 10.0, 10.0, 999.0])
    )
    foreign_line = ElevationLine(
        phi=0.0, theta=np.array([0.42, 0.44, 0.46, 0.48]), elevation=np.array([50.0, 50.0, 50.0, 50.0])
    )
    plate_0 = _plate_with_filler(0, [1.0, 0.0, 0.0], [own_line])
    plate_1 = _plate_with_filler(1, [1.0, 0.0, 0.0], [foreign_line])
    world = World(seed=0, plates=[plate_0, plate_1])

    events = reassign.reassign_misplaced_points(world)

    assert any("reassigned from plate 0 to plate 1" in e for e in events)
    updated_0 = next(p for p in world.plates if p.plate_id == 0)
    updated_1 = next(p for p in world.plates if p.plate_id == 1)
    own_result = next(l for l in updated_0.lines if np.isclose(l.phi, 0.0))
    foreign_result = next(l for l in updated_1.lines if np.isclose(l.phi, 0.0))

    assert misplaced_theta not in own_result.theta
    assert len(own_result.theta) == 3
    assert len(foreign_result.theta) == 5
    assert np.all(np.diff(foreign_result.theta) > 0)  # stays theta-sorted after insertion

    # Elevation interpolated from the point's own prior value (999.0) and its nearest
    # neighbor on the destination line (50.0, at theta=0.42) -- a straight average.
    moved_index = int(np.searchsorted(foreign_result.theta, misplaced_theta))
    assert np.isclose(foreign_result.elevation[moved_index], (999.0 + 50.0) / 2.0)


def test_reassign_leaves_well_placed_points_alone():
    line_a = ElevationLine(phi=0.0, theta=np.array([0.0, 0.1, 0.2, 0.3]), elevation=np.array([1.0, 2.0, 3.0, 4.0]))
    line_b = ElevationLine(phi=0.0, theta=np.array([2.0, 2.1, 2.2, 2.3]), elevation=np.array([5.0, 6.0, 7.0, 8.0]))
    plate_0 = _plate_with_filler(0, [1.0, 0.0, 0.0], [line_a])
    plate_1 = _plate_with_filler(1, [1.0, 0.0, 0.0], [line_b])
    world = World(seed=0, plates=[plate_0, plate_1])

    events = reassign.reassign_misplaced_points(world)

    assert events == []
    assert np.array_equal(next(p for p in world.plates if p.plate_id == 0).lines[0].theta, line_a.theta)
    assert np.array_equal(next(p for p in world.plates if p.plate_id == 1).lines[0].theta, line_b.theta)


def test_reassign_noop_with_too_few_nodes():
    tiny = Plate(
        plate_id=0,
        frame=np.eye(3),
        crust_type="oceanic",
        lines=[ElevationLine(phi=0.0, theta=np.array([0.0, 0.1]), elevation=np.array([1.0, 2.0]))],
    )
    world = World(seed=0, plates=[tiny])
    assert reassign.reassign_misplaced_points(world) == []


def test_reassign_runs_periodically_and_never_on_a_regularize_step():
    world = generate_world(seed=40, num_plates=8)
    from app import line_regrid

    assert line_regrid.REGULARIZE_INTERVAL_STEPS == reassign.REASSIGN_INTERVAL_STEPS  # same period by default
    for _ in range(30):
        step_world(world, years=1_000_000)
        # The two cadence counters must never both be reset (i.e. both passes run) on the
        # same step -- the whole point of staggering them (see world.step_world).
        assert not (world.steps_since_regularize == 0 and world.steps_since_reassign == 0)
