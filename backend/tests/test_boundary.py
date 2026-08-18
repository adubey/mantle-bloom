import numpy as np

from app import boundary, geometry, mantle
from app.plates import TARGET_LINE_SPACING_RAD, ElevationLine, Plate
from app.world import World, generate_world, step_world


def _plate(plate_id, crust_type, theta, omega, base_elevation):
    """A single-line plate at a fixed seed ([1,0,0], identity frame) so a theta offset
    directly converts to a real distance (theta_rad * PLANET_RADIUS_KM, small-angle exact
    for the scale these tests use) -- lets a test construct nodes at precise, known
    distances from a neighboring plate's own cluster."""
    frame = geometry.plate_frame_from_seed([1.0, 0.0, 0.0])
    theta = np.asarray(theta, dtype=float)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.full(len(theta), base_elevation))
    return Plate(plate_id=plate_id, frame=frame, crust_type=crust_type, omega=np.asarray(omega, dtype=float), lines=[line])


def test_closing_rate_positive_when_approaching():
    # Two points near the equator, close together in longitude.
    p_self = np.array([np.cos(0.0), np.sin(0.0), 0.0])
    p_neighbor = np.array([np.cos(0.01), np.sin(0.01), 0.0])

    # self plate spins +z (eastward motion at the equator) -> moves toward neighbor (east of it)
    approaching_omega = np.array([0.0, 0.0, 1.0])
    still_omega = np.zeros(3)
    closing = boundary.closing_rate(
        p_self[None, :], approaching_omega, still_omega, p_neighbor[None, :]
    )
    assert closing[0] > 0

    # self plate spins -z (westward motion) -> moves away from neighbor
    receding_omega = np.array([0.0, 0.0, -1.0])
    closing = boundary.closing_rate(
        p_self[None, :], receding_omega, still_omega, p_neighbor[None, :]
    )
    assert closing[0] < 0


def test_grow_or_shrink_line_extends_by_one_node_for_a_small_gap():
    line = ElevationLine(phi=0.0, theta=np.array([0.0, 0.1, 0.2]), elevation=np.array([0.0, 0.0, 0.0]))
    just_over = boundary.EXTEND_THRESHOLD_RAD * 1.05  # under 2x target spacing -> one node
    dist = np.array([just_over, just_over, just_over])
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


def test_grow_or_shrink_line_extends_by_many_nodes_for_a_large_gap():
    """A large `years` step can open a gap many spacing units wide in a single step --
    growth must close it fully, not add a fixed one node regardless of gap size (see the
    comment on _grow_or_shrink_line: under-provisioned growth here is what caused gaps.py's
    periodic gap-filling to keep spawning fresh micro-plates at the same busy boundary)."""
    line = ElevationLine(phi=0.0, theta=np.array([0.0, 0.1, 0.2]), elevation=np.array([0.0, 0.0, 0.0]))
    huge_gap = TARGET_LINE_SPACING_RAD * 10.5
    dist = np.array([huge_gap, huge_gap, huge_gap])
    closing = np.array([-boundary.TRANSFORM_RATE_THRESHOLD * 2] * 3)

    grown = boundary._grow_or_shrink_line(line, dist, closing, "oceanic")
    assert len(grown.theta) == 3 + 10 + 10  # 10 new nodes at each end
    assert np.all(np.diff(grown.theta) > 0)
    new_nodes = grown.elevation[(grown.theta > line.theta[-1]) | (grown.theta < line.theta[0])]
    assert np.all(new_nodes == boundary.DIVERGENT_RIDGE_TARGET_M)


def test_grow_or_shrink_line_extend_respects_max_nodes_cap():
    line = ElevationLine(phi=0.0, theta=np.array([0.0, 0.1]), elevation=np.array([0.0, 0.0]))
    absurd_gap = TARGET_LINE_SPACING_RAD * (boundary.MAX_EXTEND_NODES_PER_STEP + 50)
    dist = np.array([absurd_gap, absurd_gap])
    closing = np.array([-boundary.TRANSFORM_RATE_THRESHOLD * 2] * 2)

    grown = boundary._grow_or_shrink_line(line, dist, closing, "oceanic")
    added = len(grown.theta) - 2
    assert added <= boundary.MAX_EXTEND_NODES_PER_STEP * 2


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


def test_band_intensity_zero_outside_band_peaks_at_midpoint():
    dist = np.array([0.0, 0.1, 0.2, 0.3, 0.4])
    result = boundary._band_intensity(dist, inner=0.1, outer=0.3)
    assert result[0] == 0.0  # below inner edge
    assert np.isclose(result[1], 0.0)  # exactly at inner edge
    assert result[2] == 1.0  # at the band's midpoint
    assert np.isclose(result[3], 0.0)  # exactly at outer edge
    assert result[4] == 0.0  # beyond outer edge
    assert np.all((result >= 0.0) & (result <= 1.0))


def test_subduction_creates_a_volcanic_arc_band_offset_from_the_boundary():
    # A continental plate converging on an oceanic one. Points at ~82/162/224/324/424 km
    # from the oceanic cluster -- deliberately all past MERGE_THRESHOLD_RAD (~50km) so
    # _grow_or_shrink_line's own node-deletion at a close, strongly-convergent boundary
    # doesn't remove any of them out from under this test.
    rate = mantle.cm_per_yr_to_rad_per_yr(5.0)
    theta_cont = np.array([0.011, 0.0236, 0.0333, 0.049, 0.0647])
    continental = _plate(0, "continental", theta_cont, omega=[0.0, 0.0, -rate], base_elevation=500.0)
    oceanic = _plate(1, "oceanic", [-0.002, -0.0021, -0.0019], omega=[0.0, 0.0, rate], base_elevation=-3000.0)
    world = World(seed=0, plates=[continental, oceanic])

    boundary.step_boundaries(world, years=1_000_000)

    updated = next(p for p in world.plates if p.plate_id == 0).lines[0].elevation
    delta = updated - 500.0
    # Near the boundary (82km -- inside the arc's inner edge at 100km) and well past its
    # outer edge (324km, 424km): no meaningful volcanic-arc uplift.
    assert delta[0] < 5.0
    assert delta[3] < 5.0
    assert delta[4] < 5.0
    # Inside the 100-300km band (162km, 224km): real uplift, and the point closer to the
    # band's midpoint (200km) rose more than the one farther from it -- confirms the
    # non-monotonic "offset inland" shape, not just a wider plain decay.
    assert delta[1] > 50.0
    assert delta[2] > delta[1]


def test_continental_collision_crumple_zone_reaches_400km():
    # Two continental plates converging. 324km is well past the *old* ~200km reach
    # (FAR_THRESHOLD_RAD) but within COLLISION_RANGE_RAD's 400km.
    rate = mantle.cm_per_yr_to_rad_per_yr(5.0)
    theta_cont = np.array([0.011, 0.0333, 0.049])  # ~82/224/324 km
    plate_a = _plate(0, "continental", theta_cont, omega=[0.0, 0.0, -rate], base_elevation=500.0)
    plate_b = _plate(1, "continental", [-0.002, -0.0021, -0.0019], omega=[0.0, 0.0, rate], base_elevation=500.0)
    world = World(seed=0, plates=[plate_a, plate_b])

    boundary.step_boundaries(world, years=1_000_000)

    updated = next(p for p in world.plates if p.plate_id == 0).lines[0].elevation
    delta = updated - 500.0
    assert np.all(delta > 0.0)  # even the farthest point (324km) shows real uplift
    assert delta[0] > delta[1] > delta[2]  # peaks at the boundary, decays outward


def test_transform_boundary_produces_modest_uplift_within_50km():
    # A near-zero closing rate (well under TRANSFORM_RATE_THRESHOLD) -- neither strongly
    # convergent nor divergent, the definition of a transform boundary.
    rate = mantle.cm_per_yr_to_rad_per_yr(0.05)
    theta_cont = np.array([0.005, 0.010])  # ~30km, ~62km -- inside vs outside the 50km range
    continental = _plate(0, "continental", theta_cont, omega=[0.0, 0.0, -rate], base_elevation=500.0)
    oceanic = _plate(1, "oceanic", [-0.001, -0.0012, -0.0009], omega=[0.0, 0.0, rate], base_elevation=-3000.0)
    world = World(seed=0, plates=[continental, oceanic])

    boundary.step_boundaries(world, years=1_000_000)

    updated = next(p for p in world.plates if p.plate_id == 0).lines[0].elevation
    delta = updated - 500.0
    assert delta[0] > 0.0  # within 50km -- new transform uplift applies
    assert delta[1] < 1e-6  # beyond 50km -- no effect
    # "The rise isn't as big" -- confirm it stays well under real mountain-building rates.
    assert delta[0] < boundary.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR * 0.5
