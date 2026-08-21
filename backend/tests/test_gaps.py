import numpy as np
from scipy.spatial import cKDTree

from app import gaps, geometry
from app.plates import ElevationLine, Plate
from app.world import World, generate_world, step_world


def test_cluster_groups_nearby_and_splits_far_points():
    pts = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.001, 0.0, 0.0],
            [0.002, 0.0, 0.0],
            [10.0, 0.0, 0.0],
            [10.001, 0.0, 0.0],
        ]
    )
    labels = gaps.cluster_points(pts, radius=0.01)
    assert labels[0] == labels[1] == labels[2]
    assert labels[3] == labels[4]
    assert labels[0] != labels[3]


def test_cluster_empty_input():
    assert len(gaps.cluster_points(np.zeros((0, 3)), radius=0.01)) == 0


def _old_plate(plate_id: int) -> Plate:
    return Plate(plate_id=plate_id, frame=np.eye(3), crust_type="oceanic", age_steps=gaps.YOUNG_PLATE_AGE_STEPS + 1)


def _young_plate(plate_id: int) -> Plate:
    return Plate(plate_id=plate_id, frame=np.eye(3), crust_type="oceanic", age_steps=0)


def test_preferred_border_plate_single_neighbor():
    existing_points = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [-1.0, 0.0, 0.0],
        ]
    )
    existing_owner = np.array([0, 0, 1])
    tree = cKDTree(existing_points)
    plate_by_id = {0: _old_plate(0), 1: _old_plate(1)}

    # A small fraction of BORDER_RADIUS_RAD, not a hardcoded absolute offset -- stays valid
    # regardless of TARGET_LINE_SPACING_RAD (which this radius scales with).
    offset = gaps.BORDER_RADIUS_RAD * 0.5
    near_plate_0_only = geometry.rotate_vectors(
        np.array([1.0, 0.0, 0.0])[None, :], np.array([0.0, 0.0, 1.0]), offset
    )
    assert gaps._preferred_border_plate(near_plate_0_only, tree, existing_owner, plate_by_id) == 0


def test_preferred_border_plate_none_when_evenly_split_between_old_plates():
    existing_points = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    existing_owner = np.array([0, 1])
    tree = cKDTree(existing_points)
    plate_by_id = {0: _old_plate(0), 1: _old_plate(1)}

    evenly_split = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    assert gaps._preferred_border_plate(evenly_split, tree, existing_owner, plate_by_id) is None


def test_preferred_border_plate_wins_when_mostly_one_side():
    existing_points = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    existing_owner = np.array([0, 1])
    tree = cKDTree(existing_points)
    plate_by_id = {0: _old_plate(0), 1: _old_plate(1)}

    # Three points near plate 0, one incidental point near plate 1 -- plate 0 dominates.
    mostly_plate_0 = np.array(
        [[1.0, 0.0, 0.0], [0.99, 0.02, 0.0], [0.98, -0.02, 0.0], [-1.0, 0.0, 0.0]]
    )
    assert gaps._preferred_border_plate(mostly_plate_0, tree, existing_owner, plate_by_id) == 0


def test_preferred_border_plate_none_when_no_plate_nearby():
    existing_points = np.array([[1.0, 0.0, 0.0]])
    existing_owner = np.array([0])
    tree = cKDTree(existing_points)
    plate_by_id = {0: _old_plate(0)}

    far_away = np.array([[0.0, 1.0, 0.0]])
    assert gaps._preferred_border_plate(far_away, tree, existing_owner, plate_by_id) is None


def test_preferred_border_plate_prefers_a_young_plate_over_spawning_new():
    """Neither side dominates (below DOMINANT_BORDER_FRACTION), but plate 1 is young --
    give it the gap rather than returning None (which would spawn yet another new plate
    right next to one that was just created -- the sliver-chain failure this guards
    against, see YOUNG_PLATE_AGE_STEPS)."""
    existing_points = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    existing_owner = np.array([0, 1])
    tree = cKDTree(existing_points)
    plate_by_id = {0: _old_plate(0), 1: _young_plate(1)}

    evenly_split = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    assert gaps._preferred_border_plate(evenly_split, tree, existing_owner, plate_by_id) == 1


def test_preferred_border_plate_ignores_young_plate_with_too_little_border():
    # Three plates: 0 and 2 (old) split most of the border evenly, 1 (young) gets a small
    # sliver -- neither an outright dominant plate nor a young plate with enough presence,
    # so this should still fall through to spawning a new plate (None).
    existing_points = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]])
    existing_owner = np.array([0, 1, 2])
    tree = cKDTree(existing_points)
    plate_by_id = {0: _old_plate(0), 1: _young_plate(1), 2: _old_plate(2)}

    cluster_points = np.array(
        [
            [1.0, 0.0, 0.0], [0.99, 0.02, 0.0], [0.99, -0.02, 0.0], [0.98, 0.03, 0.0],  # plate 0: 4
            [0.0, 1.0, 0.0], [0.02, 0.99, 0.0],  # plate 1 (young): 2
            [-1.0, 0.0, 0.0], [-0.99, 0.02, 0.0], [-0.99, -0.02, 0.0], [-0.98, 0.03, 0.0],  # plate 2: 4
        ]
    )
    assert gaps._preferred_border_plate(cluster_points, tree, existing_owner, plate_by_id) is None


def test_absorb_gap_into_plate_grows_and_preserves_old_elevation():
    seed_xyz = geometry.normalize(np.array([1.0, 0.0, 0.0]))
    frame = geometry.plate_frame_from_seed(seed_xyz)
    theta = np.linspace(-0.05, 0.05, 5)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.full(5, 123.0))
    plate = Plate(plate_id=0, frame=frame, crust_type="continental", lines=[line])
    before_count = plate.node_count()

    # New territory just beyond the line's current theta range, at the same phi.
    gap_theta = np.linspace(0.06, 0.12, 5)
    gap_local = geometry.local_xyz(np.zeros(5), gap_theta)
    gap_points = geometry.to_world(frame, gap_local)

    rng = np.random.default_rng(0)
    claimed = gaps._absorb_gap_into_plate(plate, gap_points, rng, max_new_points=100)

    assert claimed > 0
    assert plate.node_count() > before_count
    all_elevation = np.concatenate([l.elevation for l in plate.lines])
    assert np.any(np.isclose(all_elevation, 123.0))  # old data preserved
    assert np.any(~np.isclose(all_elevation, 123.0))  # new data added, different values


def test_absorb_gap_into_plate_respects_max_new_points():
    seed_xyz = geometry.normalize(np.array([1.0, 0.0, 0.0]))
    frame = geometry.plate_frame_from_seed(seed_xyz)
    theta = np.linspace(-0.05, 0.05, 5)
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.full(5, 0.0))
    plate = Plate(plate_id=0, frame=frame, crust_type="continental", lines=[line])
    before_count = plate.node_count()

    gap_theta = np.linspace(0.06, 0.2, 20)  # plenty of gap points within GROWTH_RING_RAD
    gap_local = geometry.local_xyz(np.zeros(20), gap_theta)
    gap_points = geometry.to_world(frame, gap_local)

    rng = np.random.default_rng(0)
    claimed_capped = gaps._absorb_gap_into_plate(plate, gap_points.copy(), rng, max_new_points=3)
    assert claimed_capped == 3
    grown_capped = plate.node_count()
    assert grown_capped > before_count

    # An uncapped absorb of the same gap should claim (and end up covering) more.
    plate_uncapped = Plate(plate_id=0, frame=frame, crust_type="continental", lines=[line])
    claimed_uncapped = gaps._absorb_gap_into_plate(
        plate_uncapped, gap_points.copy(), rng, max_new_points=len(gap_points)
    )
    assert claimed_uncapped > claimed_capped
    assert plate_uncapped.node_count() > grown_capped


def test_fill_gaps_caps_total_absorption_per_plate_per_call():
    world = generate_world(seed=60, num_plates=8)
    plate = max(world.plates, key=lambda p: p.node_count())
    ordered_lines = sorted(plate.lines, key=lambda l: l.phi)
    mid = len(ordered_lines) // 2
    # Carve several separate, decently large interior holes -- if they're each resolved
    # independently without a shared per-call budget, the plate could grow by far more than
    # MAX_ABSORB_NODES_PER_PLATE_PER_CALL in one fill_gaps call.
    for line in ordered_lines[max(mid - 6, 0) : mid + 6 : 2]:
        if len(line.theta) < 20:
            continue
        center = len(line.theta) // 2
        keep = np.ones(len(line.theta), dtype=bool)
        keep[max(center - 6, 1) : min(center + 7, len(line.theta) - 1)] = False
        line.theta = line.theta[keep]
        line.elevation = line.elevation[keep]

    node_count_before = plate.node_count()
    gaps.fill_gaps(world)

    assert plate.node_count() - node_count_before <= gaps.MAX_ABSORB_NODES_PER_PLATE_PER_CALL


def test_spawn_plate_from_gap_creates_valid_plate():
    gap_theta = np.linspace(-0.05, 0.05, 8)
    gap_local = geometry.local_xyz(np.zeros(8), gap_theta)
    gap_points = geometry.to_world(np.eye(3), gap_local)

    world = World(seed=0, plates=[], mantle_centers=[], next_plate_id=5)
    rng = np.random.default_rng(1)
    new_plate = gaps._spawn_plate_from_gap(world, gap_points, rng)

    assert new_plate.plate_id == 5
    assert world.next_plate_id == 6
    assert new_plate.crust_type == "oceanic"
    assert new_plate.node_count() > 0
    assert np.allclose(new_plate.frame @ new_plate.frame.T, np.eye(3), atol=1e-9)


def test_fill_gaps_noop_when_no_plates():
    world = World(seed=0, plates=[])
    gaps.fill_gaps(world)
    assert world.plates == []


def test_fill_gaps_is_a_noop_immediately_after_generation():
    """Generation already tiles the whole sphere by construction (nearest-seed test), so
    there should be essentially nothing for fill_gaps to do right after generate_world."""
    world = generate_world(seed=8, num_plates=10)
    plate_count_before = len(world.plates)
    node_count_before = sum(p.node_count() for p in world.plates)

    gaps.fill_gaps(world)

    assert len(world.plates) == plate_count_before
    assert sum(p.node_count() for p in world.plates) == node_count_before


def test_fill_gaps_absorbs_an_interior_hole_into_its_one_bordering_plate():
    # node_density pinned to 1.0 (not DEFAULT_NODE_DENSITY): the carved hole below is sized
    # by a fixed line/index-count window, not a fixed physical area, so at a higher density
    # (finer line/theta spacing) the same window covers a physically smaller hole -- while
    # MIN_GAP_POINTS' own threshold scales *up* with density (see _min_gap_points) -- and the
    # carved hole no longer clears it. 1.0 is the density this window was actually tuned
    # against ("comfortably bigger than MIN_GAP_POINTS" below).
    world = generate_world(seed=9, num_plates=8, node_density=1.0)
    plate = max(world.plates, key=lambda p: p.node_count())
    # Carve a sizeable hole out of several consecutive interior lines -- comfortably bigger
    # than MIN_GAP_POINTS once resampled onto the global detection grid, and well inside the
    # plate's own territory, so it can only ever be bordered by this one plate.
    ordered_lines = sorted(plate.lines, key=lambda l: l.phi)
    mid = len(ordered_lines) // 2
    for line in ordered_lines[max(mid - 8, 0) : mid + 8]:
        if len(line.theta) < 40:
            continue
        center = len(line.theta) // 2
        keep = np.ones(len(line.theta), dtype=bool)
        keep[max(center - 14, 1) : min(center + 15, len(line.theta) - 1)] = False
        line.theta = line.theta[keep]
        line.elevation = line.elevation[keep]

    plate_count_before = len(world.plates)
    node_count_before = plate.node_count()

    gaps.fill_gaps(world)

    assert len(world.plates) == plate_count_before  # no new plate -- one dominant plate
    assert plate.node_count() > node_count_before  # hole refilled


def test_fill_gaps_leaves_no_large_gaps_after_a_long_simulation():
    world = generate_world(seed=55, num_plates=10)
    for _ in range(25):
        step_world(world, years=6_000_000)
    # A single fill_gaps call only claims one ring of a gap (see GROWTH_RING_RAD) --
    # gradual healing, not an instant full close -- so give it several passes, matching how
    # step_world actually calls it repeatedly over time.
    for _ in range(10):
        gaps.fill_gaps(world)

    existing_points, _ = gaps._all_existing_points(world)
    tree = cKDTree(existing_points)
    remaining_gap_points = gaps._find_gap_points(tree)
    if len(remaining_gap_points) == 0:
        return
    labels = gaps.cluster_points(remaining_gap_points, gaps.CLUSTER_RADIUS_RAD)
    counts = np.bincount(labels)
    assert counts.max() < gaps.MIN_GAP_POINTS


def test_gap_fill_is_deterministic_for_same_seed_and_step_sequence():
    def run():
        w = generate_world(seed=77, num_plates=10)
        for _ in range(12):
            step_world(w, years=5_000_000)
        return w

    w1, w2 = run(), run()
    p1s = sorted(w1.plates, key=lambda p: p.plate_id)
    p2s = sorted(w2.plates, key=lambda p: p.plate_id)
    assert [p.plate_id for p in p1s] == [p.plate_id for p in p2s]
    for p1, p2 in zip(p1s, p2s):
        assert p1.crust_type == p2.crust_type
        assert np.allclose(p1.frame, p2.frame)
        assert p1.node_count() == p2.node_count()
