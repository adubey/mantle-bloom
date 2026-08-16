import numpy as np

from app import boundary
from app.plates import ElevationLine
from app.world import generate_world, step_world


def test_closing_rate_positive_when_approaching():
    # Two points near the equator, close together in longitude.
    p_self = np.array([np.cos(0.0), np.sin(0.0), 0.0])
    p_neighbor = np.array([np.cos(0.01), np.sin(0.01), 0.0])

    # self plate spins +z (eastward motion at the equator) -> moves toward neighbor (east of it)
    approaching_omega = np.array([0.0, 0.0, 1.0])
    still_omega = np.zeros(3)
    closing = boundary._closing_rate(
        p_self[None, :], approaching_omega, still_omega, p_neighbor[None, :]
    )
    assert closing[0] > 0

    # self plate spins -z (westward motion) -> moves away from neighbor
    receding_omega = np.array([0.0, 0.0, -1.0])
    closing = boundary._closing_rate(
        p_self[None, :], receding_omega, still_omega, p_neighbor[None, :]
    )
    assert closing[0] < 0


def test_grow_or_shrink_line_extends_when_divergent_and_far():
    line = ElevationLine(phi=0.0, theta=np.array([0.0, 0.1, 0.2]), elevation=np.array([0.0, 0.0, 0.0]))
    far = boundary.EXTEND_THRESHOLD_RAD * 2
    dist = np.array([far, far, far])
    closing = np.array(
        [-boundary.TRANSFORM_RATE_THRESHOLD * 2] * 3
    )  # strongly divergent everywhere

    grown = boundary._grow_or_shrink_line(line, dist, closing, "oceanic")
    assert len(grown.theta) == 5  # one node added at each end
    assert grown.theta[0] < line.theta[0]
    assert grown.theta[-1] > line.theta[-1]
    assert grown.elevation[0] == boundary.DIVERGENT_RIDGE_TARGET_M
    assert grown.elevation[-1] == boundary.DIVERGENT_RIDGE_TARGET_M
    assert np.all(np.diff(grown.theta) > 0)


def test_grow_or_shrink_line_shrinks_when_convergent_and_close():
    line = ElevationLine(
        phi=0.0, theta=np.array([0.0, 0.1, 0.2, 0.3]), elevation=np.array([1.0, 2.0, 3.0, 4.0])
    )
    close = boundary.MERGE_THRESHOLD_RAD / 2
    dist = np.array([close, close, close, close])
    closing = np.array([boundary.TRANSFORM_RATE_THRESHOLD * 2] * 4)  # strongly convergent

    shrunk = boundary._grow_or_shrink_line(line, dist, closing, "continental")
    assert len(shrunk.theta) == 2  # one node dropped from each end
    assert np.array_equal(shrunk.theta, line.theta[1:-1])


def test_grow_or_shrink_line_never_deletes_the_last_node():
    line = ElevationLine(phi=0.0, theta=np.array([0.0]), elevation=np.array([1.0]))
    close = boundary.MERGE_THRESHOLD_RAD / 2
    dist = np.array([close])
    closing = np.array([boundary.TRANSFORM_RATE_THRESHOLD * 2])

    result = boundary._grow_or_shrink_line(line, dist, closing, "continental")
    assert len(result.theta) == 1


def test_grow_or_shrink_line_untouched_when_transform_or_far_middle():
    line = ElevationLine(phi=0.0, theta=np.array([0.0, 0.1, 0.2]), elevation=np.array([1.0, 2.0, 3.0]))
    dist = np.array([boundary.FAR_THRESHOLD_RAD * 10] * 3)  # nowhere near any boundary
    closing = np.array([0.0, 0.0, 0.0])

    result = boundary._grow_or_shrink_line(line, dist, closing, "continental")
    assert np.array_equal(result.theta, line.theta)
    assert np.array_equal(result.elevation, line.elevation)


def test_stepping_with_boundary_evolution_keeps_lines_sorted_and_elevation_bounded():
    world = generate_world(seed=21, num_plates=8)
    for _ in range(15):
        step_world(world, years=3_000_000)

    for plate in world.plates:
        for line in plate.lines:
            assert np.all(np.diff(line.theta) > 0), "line thetas must stay strictly ascending"
            assert np.all(line.elevation >= boundary.MIN_ELEVATION_M)
            assert np.all(line.elevation <= boundary.MAX_ELEVATION_M)


def test_boundary_evolution_changes_node_counts_over_time():
    world = generate_world(seed=22, num_plates=10)
    before = sum(p.node_count() for p in world.plates)
    for _ in range(20):
        step_world(world, years=4_000_000)
    after = sum(p.node_count() for p in world.plates)
    assert after != before
