import numpy as np
from scipy.spatial import cKDTree

from app import geometry, mantle, merge_split
from app.plates import ElevationLine, PlateWithLines, line_spacing_rad
from app.world import World, generate_world, step_world


def _test_plate(plate_id, seed_xyz, crust_type, theta, elevation):
    """A plate with the given line at phi=0 (what each test actually exercises) plus a
    second, far-away placeholder line at a different latitude -- purely so the plate has
    more than one line. Otherwise apply_topology_changes's "no land left" pruning (a plate
    reduced to a single line, see merge_split.remove_defunct_plates) would remove these
    synthetic single-line test plates before the test's own logic ever ran."""
    frame = geometry.plate_frame_from_seed(seed_xyz)
    line = ElevationLine(phi=0.0, theta=np.asarray(theta, dtype=float), elevation=np.asarray(elevation, dtype=float))
    filler = ElevationLine(phi=1.0, theta=np.array([0.0, 0.1]), elevation=np.zeros(2))
    return PlateWithLines(plate_id=plate_id, frame=frame, crust_type=crust_type, lines=[line, filler])


def _converging_omega_pair(rate_cm_per_yr=5.0):
    """omega for (keep, absorb) such that -- given absorb is seeded at a small *positive*
    angle from keep around +z (see the tiny_angle rotation below) -- each plate's material
    moves toward the other's, matching test_boundary.py's closing-rate sign convention."""
    rate = mantle.cm_per_yr_to_rad_per_yr(rate_cm_per_yr)
    return np.array([0.0, 0.0, rate]), np.array([0.0, 0.0, -rate])


def _converging_pair_world(seed=123):
    # See test_find_continental_collision_pairs_detects_close_and_converging_plates for why
    # these offsets are fractions of MERGE_CONTACT_DISTANCE_RAD rather than hardcoded angles.
    d = merge_split.MERGE_CONTACT_DISTANCE_RAD
    keep = _test_plate(0, [1.0, 0.0, 0.0], "continental", np.linspace(-d, 0.0, 6), np.zeros(6))
    seed_absorb = geometry.rotate_vectors(
        np.array([1.0, 0.0, 0.0])[None, :], axis=np.array([0.0, 0.0, 1.0]), angle=d * 0.3
    )[0]
    absorb = _test_plate(1, seed_absorb, "continental", np.linspace(-d * 0.2, d * 0.2, 6), np.full(6, 50.0))
    keep_omega, absorb_omega = _converging_omega_pair()
    keep.set_omega(keep_omega)
    absorb.set_omega(absorb_omega)
    return World(seed=seed, plates=[keep, absorb], next_plate_id=2)


def test_remove_defunct_plates_drops_empty_and_single_line_plates():
    survives = _test_plate(0, [1.0, 0.0, 0.0], "oceanic", [0.0, 0.1], [1.0, 2.0])
    empty = PlateWithLines(plate_id=1, frame=np.eye(3), crust_type="oceanic", lines=[])
    # A single line, however many nodes it carries, counts as "no land left" too -- not just
    # zero lines (see the new condition merge_split.remove_defunct_plates checks for).
    sliver = PlateWithLines(
        plate_id=2,
        frame=np.eye(3),
        crust_type="oceanic",
        lines=[ElevationLine(phi=0.0, theta=np.array([0.0, 0.1, 0.2]), elevation=np.array([1.0, 2.0, 3.0]))],
    )
    world = World(seed=0, plates=[survives, empty, sliver])

    merge_split.remove_defunct_plates(world)
    assert [p.plate_id for p in world.plates] == [0]


def test_find_continental_collision_pairs_detects_close_and_converging_plates():
    # Every offset here is a fraction of MERGE_CONTACT_DISTANCE_RAD, not a hardcoded
    # absolute angle -- stays valid regardless of TARGET_LINE_SPACING_RAD (which that
    # constant scales with). Non-interleaved: keep's world-angle range ([-d, 0]) sits
    # entirely left of absorb's ([0.1d, 0.5d]), so nearest-neighbor pairing is unambiguous
    # and every closing-rate sign comes out consistent (an earlier version of this test used
    # symmetric, overlapping ranges around each plate's own seed, which interleaved the two
    # point sets and made nearest-neighbor pairing -- and therefore the sign of the closing
    # rate -- inconsistent from point to point). Most, but not all, of keep's points end up
    # within d of their nearest absorb point -- comfortably above MERGE_MIN_CONTACT_NODES.
    d = merge_split.MERGE_CONTACT_DISTANCE_RAD
    keep = _test_plate(0, [1.0, 0.0, 0.0], "continental", np.linspace(-d, 0.0, 6), np.zeros(6))
    seed_absorb = geometry.rotate_vectors(
        np.array([1.0, 0.0, 0.0])[None, :], axis=np.array([0.0, 0.0, 1.0]), angle=d * 0.3
    )[0]
    absorb = _test_plate(1, seed_absorb, "continental", np.linspace(-d * 0.2, d * 0.2, 6), np.full(6, 50.0))
    keep_omega, absorb_omega = _converging_omega_pair()
    keep.set_omega(keep_omega)
    absorb.set_omega(absorb_omega)

    world = World(seed=0, plates=[keep, absorb], next_plate_id=2)
    pairs = merge_split.find_continental_collision_pairs(world)
    assert (0, 1) in pairs


def test_find_continental_collision_pairs_ignores_close_but_not_converging_plates():
    """The actual bug this guards against: plates.py's tiling has no gaps, so any two
    neighboring plates are already touching the moment they're generated, regardless of
    whether that boundary is convergent, divergent, or transform. Proximity with zero (or
    non-convergent) relative motion must not be treated as a collision -- confirmed
    directly to fire on the very first simulation step, for any `years`, before this check
    was added."""
    theta = np.linspace(-0.01, 0.01, 5)
    keep = _test_plate(0, [1.0, 0.0, 0.0], "continental", theta, np.zeros(5))

    tiny_angle = merge_split.MERGE_CONTACT_DISTANCE_RAD * 0.3
    seed_absorb = geometry.rotate_vectors(
        np.array([1.0, 0.0, 0.0])[None, :], axis=np.array([0.0, 0.0, 1.0]), angle=tiny_angle
    )[0]
    absorb = _test_plate(1, seed_absorb, "continental", theta, np.full(5, 50.0))
    # Both plates motionless -- close, but nothing is actually colliding.
    world = World(seed=0, plates=[keep, absorb], next_plate_id=2)
    assert merge_split.find_continental_collision_pairs(world) == []


def test_find_continental_collision_pairs_ignores_far_plates():
    theta = np.linspace(-0.01, 0.01, 5)
    keep = _test_plate(0, [1.0, 0.0, 0.0], "continental", theta, np.zeros(5))
    far = _test_plate(1, [0.0, 1.0, 0.0], "continental", theta, np.zeros(5))
    keep_omega, far_omega = _converging_omega_pair()
    keep.set_omega(keep_omega)
    far.set_omega(far_omega)
    world = World(seed=0, plates=[keep, far], next_plate_id=2)
    assert merge_split.find_continental_collision_pairs(world) == []


def test_find_continental_collision_pairs_ignores_oceanic():
    theta = np.linspace(-0.01, 0.01, 5)
    tiny_angle = merge_split.MERGE_CONTACT_DISTANCE_RAD * 0.3
    seed_absorb = geometry.rotate_vectors(
        np.array([1.0, 0.0, 0.0])[None, :], axis=np.array([0.0, 0.0, 1.0]), angle=tiny_angle
    )[0]
    keep = _test_plate(0, [1.0, 0.0, 0.0], "oceanic", theta, np.zeros(5))
    absorb = _test_plate(1, seed_absorb, "continental", theta, np.zeros(5))
    keep_omega, absorb_omega = _converging_omega_pair()
    keep.set_omega(keep_omega)
    absorb.set_omega(absorb_omega)
    world = World(seed=0, plates=[keep, absorb], next_plate_id=2)
    assert merge_split.find_continental_collision_pairs(world) == []


def test_collision_does_not_merge_before_sustained_threshold():
    world = _converging_pair_world()
    threshold = merge_split._collision_threshold_years(world.seed, (0, 1))
    step_years = 5_000_000

    years_run = 0.0
    while years_run + step_years < threshold:
        events = merge_split.apply_topology_changes(world, step_years)
        years_run += step_years
        assert len(world.plates) == 2, f"merged early at {years_run} of {threshold} years"
        assert not any("merged" in e for e in events)


def test_collision_merges_once_sustained_threshold_is_crossed(monkeypatch):
    # This test isolates the sustained-duration mechanism -- neutralize the independent
    # size-based probability gate (see merge_split._merge_probability), since a synthetic
    # 2-plate world makes any pair's combined share of total world nodes trivially large.
    monkeypatch.setattr(merge_split, "MERGE_SIZE_UNLIKELY_FRACTION", 1e9)
    world = _converging_pair_world()
    threshold = merge_split._collision_threshold_years(world.seed, (0, 1))

    # Run right up to (but not past) the threshold first...
    step_years = 5_000_000
    years_run = 0.0
    while years_run + step_years < threshold:
        merge_split.apply_topology_changes(world, step_years)
        years_run += step_years
    assert len(world.plates) == 2

    # ...then one more step should cross it.
    events = merge_split.apply_topology_changes(world, step_years * 2)
    assert len(world.plates) == 1
    assert any("merged" in e and "million years" in e for e in events)


def test_collision_progress_resets_if_convergence_stops():
    world = _converging_pair_world()
    merge_split.apply_topology_changes(world, 10_000_000)
    assert (0, 1) in world.collision_progress

    # Convergence stops (both plates now motionless) -- progress must be dropped, not paused.
    for p in world.plates:
        p.set_omega(np.zeros(3))
    merge_split.apply_topology_changes(world, 10_000_000)
    assert (0, 1) not in world.collision_progress


def test_apply_topology_changes_merges_at_most_one_pair_per_call(monkeypatch):
    # This test isolates the "at most one merge per call" cap -- neutralize the independent
    # size-based probability gate (see merge_split._merge_probability), since a synthetic
    # 3-plate world makes any pair's combined share of total world nodes trivially large.
    monkeypatch.setattr(merge_split, "MERGE_SIZE_UNLIKELY_FRACTION", 1e9)
    # Three mutually close, mutually converging continental plates clustered together --
    # every pairwise combination should register as a collision candidate.
    theta = np.linspace(-0.003, 0.003, 6)
    seeds = [
        [1.0, 0.0, 0.0],
        geometry.rotate_vectors(np.array([1.0, 0.0, 0.0])[None, :], np.array([0.0, 0.0, 1.0]), 0.005)[0],
        geometry.rotate_vectors(np.array([1.0, 0.0, 0.0])[None, :], np.array([0.0, 0.0, 1.0]), 0.010)[0],
    ]
    rate = mantle.cm_per_yr_to_rad_per_yr(5.0)
    omegas = [np.array([0.0, 0.0, rate]), np.zeros(3), np.array([0.0, 0.0, -rate])]
    plates = [
        _test_plate(i, seeds[i], "continental", theta, np.zeros(6)) for i in range(3)
    ]
    for p, om in zip(plates, omegas):
        p.set_omega(om)
    world = World(seed=7, plates=plates, next_plate_id=3)

    pairs = merge_split.find_continental_collision_pairs(world)
    assert len(pairs) >= 2, f"expected at least 2 candidate pairs to test the cap, got {pairs}"

    # Force every candidate pair to already be past its merge threshold.
    for pair in pairs:
        world.collision_progress[pair] = merge_split._collision_threshold_years(world.seed, pair) + 1
    events = merge_split.apply_topology_changes(world, 1_000_000)

    assert len(world.plates) == 2  # exactly one merge happened, not a cascade
    assert sum("merged" in e for e in events) == 1


def test_merge_plates_leaves_one_plate_with_nonzero_nodes():
    theta = np.linspace(-0.01, 0.01, 5)
    keep = _test_plate(0, [1.0, 0.0, 0.0], "continental", theta, np.zeros(5))
    tiny_angle = merge_split.MERGE_CONTACT_DISTANCE_RAD * 0.3
    seed_absorb = geometry.rotate_vectors(
        np.array([1.0, 0.0, 0.0])[None, :], axis=np.array([0.0, 0.0, 1.0]), angle=tiny_angle
    )[0]
    absorb = _test_plate(1, seed_absorb, "continental", theta, np.full(5, 50.0))
    world = World(seed=0, plates=[keep, absorb], next_plate_id=2)

    merge_split.merge_plates(world, id_keep=0, id_absorb=1)

    assert [p.plate_id for p in world.plates] == [0]
    assert world.plates[0].node_count() > 0


def test_merge_plates_does_not_claim_another_plates_territory():
    """Old bug: merge_plates' is_owned only checked distance to the merging pair's own old
    points, so if either parent carried a stray far-flung node (as a plate that's already
    been through an earlier merge can, see merge_plates' own docstring), the resample would
    claim lattice cells near that stray point even where a completely unrelated,
    still-living plate already owns the space -- confirmed directly as the cause of large
    cross-plate node overlap in a real, long-run save file."""
    theta = np.linspace(-0.01, 0.01, 5)
    keep = _test_plate(0, [1.0, 0.0, 0.0], "continental", theta, np.zeros(5))
    tiny_angle = merge_split.MERGE_CONTACT_DISTANCE_RAD * 0.3
    seed_absorb = geometry.rotate_vectors(
        np.array([1.0, 0.0, 0.0])[None, :], axis=np.array([0.0, 0.0, 1.0]), angle=tiny_angle
    )[0]
    absorb = _test_plate(1, seed_absorb, "continental", theta, np.full(5, 50.0))

    # A bystander plate on the opposite side of the sphere, with no relation to keep/absorb's
    # actual collision. Give `keep` one stray node planted right at the bystander's own seed
    # -- standing in for a scattered leftover point a plate that's already been through a
    # prior merge can carry.
    bystander_seed = np.array([-1.0, 0.0, 0.0])
    bystander = _test_plate(2, bystander_seed, "continental", theta, np.zeros(5))
    stray_local = geometry.to_local(keep.frame, bystander_seed)
    stray_phi, stray_theta = geometry.xyz_to_latlon(stray_local)
    stray_line = ElevationLine(phi=float(stray_phi), theta=np.array([float(stray_theta)]), elevation=np.array([0.0]))
    keep.set_lines([*keep.lines, stray_line])

    world = World(seed=0, plates=[keep, absorb, bystander], next_plate_id=3)
    bystander_points_before, _ = bystander.all_points_and_elevation()

    merge_split.merge_plates(world, id_keep=0, id_absorb=1)

    merged = next(p for p in world.plates if p.plate_id == 0)
    merged_points, _ = merged.all_points_and_elevation()
    bystander_after = next(p for p in world.plates if p.plate_id == 2)
    bystander_points_after, _ = bystander_after.all_points_and_elevation()

    # bystander is untouched by the merge...
    assert bystander_points_after.shape == bystander_points_before.shape
    assert np.allclose(np.sort(bystander_points_after, axis=0), np.sort(bystander_points_before, axis=0))

    # ...and the merged plate didn't claim anything near bystander's territory, even though
    # the stray node planted right at bystander's seed would have made it do so under the old
    # (merging-pair-only) exclusivity check.
    spacing_rad = line_spacing_rad(world.node_density)
    coverage_radius_rad = merge_split.MERGE_COVERAGE_RADIUS_RAD * (spacing_rad / merge_split.TARGET_LINE_SPACING_RAD)
    dist, _ = cKDTree(bystander_points_after).query(merged_points)
    assert dist.min() > coverage_radius_rad


def test_merge_probability_decreases_with_combined_size_and_floors():
    theta = np.linspace(-0.01, 0.01, 5)
    small_a = _test_plate(0, [1.0, 0.0, 0.0], "continental", theta, np.zeros(5))
    small_b = _test_plate(1, [1.0, 0.0, 0.0], "continental", theta, np.zeros(5))
    filler = _test_plate(2, [0.0, 1.0, 0.0], "continental", np.linspace(-0.01, 0.01, 500), np.zeros(500))
    world = World(seed=0, plates=[small_a, small_b, filler], next_plate_id=3)

    small_pair_probability = merge_split._merge_probability(world, (0, 1))
    assert small_pair_probability > 0.85  # small combined share of a big filler-padded world

    big_a = _test_plate(3, [1.0, 0.0, 0.0], "continental", np.linspace(-0.01, 0.01, 5000), np.zeros(5000))
    big_b = _test_plate(4, [1.0, 0.0, 0.0], "continental", np.linspace(-0.01, 0.01, 5000), np.zeros(5000))
    world_dominant = World(seed=0, plates=[big_a, big_b, filler], next_plate_id=5)
    big_pair_probability = merge_split._merge_probability(world_dominant, (3, 4))
    # combined share well past the unlikely threshold, so this bottoms out at the floor
    assert np.isclose(big_pair_probability, merge_split.MERGE_PROBABILITY_FLOOR)

    assert big_pair_probability < small_pair_probability


def test_merge_plates_jumps_survivor_to_front_of_check_queue():
    theta = np.linspace(-0.01, 0.01, 5)
    keep = _test_plate(0, [1.0, 0.0, 0.0], "continental", theta, np.zeros(5))
    tiny_angle = merge_split.MERGE_CONTACT_DISTANCE_RAD * 0.3
    seed_absorb = geometry.rotate_vectors(
        np.array([1.0, 0.0, 0.0])[None, :], axis=np.array([0.0, 0.0, 1.0]), angle=tiny_angle
    )[0]
    absorb = _test_plate(1, seed_absorb, "continental", theta, np.full(5, 50.0))
    # id_keep already queued, part-way back -- the merge should move it to the front, not
    # just leave it where reconciliation would otherwise find it.
    world = World(seed=0, plates=[keep, absorb], next_plate_id=2, plate_check_queue=[5, 0, 6])

    merge_split.merge_plates(world, id_keep=0, id_absorb=1)

    assert world.plate_check_queue[0] == 0
    assert world.plate_check_queue.count(0) == 1


def test_maybe_split_plate_returns_none_for_small_plate():
    small = _test_plate(0, [1.0, 0.0, 0.0], "continental", [0.0, 0.1], [0.0, 0.0])
    world = World(seed=0, plates=[small], mantle_centers=[], next_plate_id=1)
    assert merge_split.maybe_split_plate(world, small) is None


def test_maybe_split_plate_splits_under_engineered_flow_divergence():
    seed_xyz = np.array([1.0, 0.0, 0.0])
    frame = geometry.plate_frame_from_seed(seed_xyz)
    # Comfortably more than 2 * SPLIT_MIN_NODES, so each half still clears the threshold.
    theta = np.linspace(-0.5, 0.5, 4 * merge_split.SPLIT_MIN_NODES)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.zeros_like(theta))
    plate = PlateWithLines(plate_id=0, frame=frame, crust_type="continental", lines=[line], age_steps=merge_split.SPLIT_MIN_AGE_STEPS)

    # Two strong, oppositely-signed convection centers straddling the plate so its two
    # halves get pushed in genuinely different directions -- a single rigid rotation can't
    # fit both, which is exactly the split trigger.
    west_pt = geometry.to_world(frame, geometry.local_xyz(np.array([0.0]), np.array([-0.4]))[0])
    east_pt = geometry.to_world(frame, geometry.local_xyz(np.array([0.0]), np.array([0.4]))[0])
    strong_rate = mantle.MANTLE_FLOW_REFERENCE_RATE * 20
    centers = [
        mantle.ConvectionCenter(position=west_pt, strength=strong_rate, falloff=0.3),
        mantle.ConvectionCenter(position=east_pt, strength=-strong_rate, falloff=0.3),
    ]
    # node_density pinned to 1.0 (not DEFAULT_NODE_DENSITY): the plate above is sized to
    # 4 * SPLIT_MIN_NODES total (2x that threshold per resulting half), but maybe_split_plate's
    # own min_nodes scales *up* with node_density (see that function's own comment) -- at
    # DEFAULT_NODE_DENSITY that threshold would be 4x SPLIT_MIN_NODES, well above what either
    # half of this engineered plate has, so the split it's testing for would be rejected.
    world = World(seed=0, plates=[plate], mantle_centers=centers, next_plate_id=1, node_density=1.0)

    points, _ = plate.all_points_and_elevation()
    velocities = mantle.flow_at(points, centers)
    plate.set_omega(mantle.fit_euler_pole(points, velocities))

    result = merge_split.maybe_split_plate(world, plate)
    assert result is not None
    plate_a, plate_b = result
    assert plate_a.plate_id == plate.plate_id
    assert plate_b.plate_id == 1  # drawn from world.next_plate_id
    assert world.next_plate_id == 2

    total_before = sum(len(l.theta) for l in plate.lines)
    total_after = plate_a.node_count() + plate_b.node_count()
    assert total_after == total_before
