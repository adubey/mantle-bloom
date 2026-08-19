import numpy as np

from app import geometry, volcanism
from app.plates import ElevationLine, Plate
from app.world import World, generate_world, step_world


def _line_plate(plate_id: int, theta: np.ndarray) -> Plate:
    line = ElevationLine(phi=0.0, theta=theta, elevation=np.full(len(theta), 200.0))
    return Plate(plate_id=plate_id, frame=np.eye(3), crust_type="continental", lines=[line])


def _two_plate_world(gap: float, d: float = 0.01, n: int = 10) -> World:
    """Two plates, each a single dense line (spacing d) at phi=0 -- plate A's line ends at
    theta=(n-1)*d, plate B's line starts at that same value plus `gap`. Both share the
    identity frame purely so their local theta coordinates line up directly with world-space
    angular distance, for an exact, hand-verifiable test geometry. An *open* chain (plate A's
    own start and plate B's own end are dead ends, near nothing) -- fine for a "does this
    fire" check, but see _two_plate_ring_world for why a real "should NOT fire" control needs
    a closed ring instead."""
    theta_a = d * np.arange(n)
    plate_a = _line_plate(0, theta_a)
    theta_b = theta_a[-1] + gap + d * np.arange(n)
    plate_b = _line_plate(1, theta_b)
    return World(seed=0, plates=[plate_a, plate_b], next_plate_id=2)


def _two_plate_ring_world(gap_multiplier: float, n: int = 10) -> World:
    """Two plates whose lines close into a full ring around the sphere (both junctions --
    A's end to B's start, and B's end wrapping back to A's start -- the same size,
    `gap_multiplier * d`), so neither plate has an open "dead end" boundary point isolated
    from everything. Real generated plates always fully tile the sphere this way (every
    boundary point has a genuine neighbor on some side); a plain open two-plate chain
    (_two_plate_world) doesn't, and its own dead ends read as spuriously "infinitely far from
    any other plate" -- which is real outlier detection working correctly on an unrealistic
    shape, not a bug, but it means a "does this normal boundary NOT fire" control needs this
    closed-ring construction instead to mean what it says."""
    d = np.pi / ((n - 1) + gap_multiplier)
    gap = gap_multiplier * d
    theta_a = d * np.arange(n)
    plate_a = _line_plate(0, theta_a)
    theta_b = theta_a[-1] + gap + d * np.arange(n)
    plate_b = _line_plate(1, theta_b)
    return World(seed=0, plates=[plate_a, plate_b], next_plate_id=2)


def test_whole_world_median_spacing_matches_known_uniform_spacing():
    d = 0.01
    world = _two_plate_world(gap=d, d=d)  # gap == d -- one uniform lattice, no outlier at all
    median = volcanism._whole_world_median_spacing(world)
    assert np.isclose(median, d, rtol=0.05)


def test_whole_world_median_spacing_zero_for_degenerate_world():
    assert volcanism._whole_world_median_spacing(World(seed=0, plates=[])) == 0.0
    single = World(seed=0, plates=[_line_plate(0, np.array([0.0]))])
    assert volcanism._whole_world_median_spacing(single) == 0.0


def test_boundary_points_are_each_lines_two_endpoints():
    world = _two_plate_world(gap=0.01)
    points, owner = volcanism._boundary_points(world)
    # Each plate has exactly one line, so its outline is exactly that line's 2 endpoints.
    assert len(points) == 4
    assert sorted(owner.tolist()) == [0, 0, 1, 1]


def test_nearest_other_plate_boundary_finds_the_true_nearest_cross_plate_point():
    points = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, 0.9], [1.0, 0.0, 0.0], [1.0, 0.1, 0.0]])
    owner = np.array([0, 0, 1, 1])
    dist, idx = volcanism._nearest_other_plate_boundary(points, owner)
    # Point 0 (plate 0) is closer to point 2 than point 3 (both on plate 1).
    assert idx[0] == 2
    # Point 2 (plate 1) is closer to point 1 than point 0 (both on plate 0).
    assert idx[2] == 1
    assert np.all(np.isfinite(dist))


def test_detect_and_spawn_volcanic_fields_fires_on_a_wide_cross_plate_gap():
    d = 0.01
    world = _two_plate_world(gap=50 * d, d=d)  # gap >> GAP_OUTLIER_FACTOR (3x) times the dense-lattice median
    events = volcanism.detect_and_spawn_volcanic_fields(world)
    assert len(events) == 1
    assert "volcanic field" in events[0].lower()

    new_plate = next(p for p in world.plates if p.plate_id not in (0, 1))
    assert new_plate.crust_type == "continental"
    assert new_plate.plate_id in world.volcanic_field_plate_ids
    assert new_plate.node_count() > 0
    for line in new_plate.lines:
        assert np.all(line.is_volcano)
        assert np.all(line.volcano_active_years_remaining >= volcanism.VOLCANO_ACTIVE_MIN_YEARS)
        assert np.all(line.volcano_active_years_remaining <= volcanism.VOLCANO_ACTIVE_MAX_YEARS)


def test_detect_and_spawn_volcanic_fields_does_not_fire_on_a_normal_gap():
    world = _two_plate_ring_world(gap_multiplier=1.2)  # a touching, ordinary boundary -- not an outlier
    events = volcanism.detect_and_spawn_volcanic_fields(world)
    assert events == []
    assert len(world.plates) == 2
    assert world.volcanic_field_plate_ids == set()


def test_detect_and_spawn_volcanic_fields_skips_points_already_on_a_tracked_field():
    d = 0.01
    world = _two_plate_world(gap=50 * d, d=d)
    world.volcanic_field_plate_ids.add(0)  # plate 0 is already an active field
    events = volcanism.detect_and_spawn_volcanic_fields(world)
    assert events == []  # plate 0's boundary points are excluded as source candidates


def test_volcano_fraction():
    line_all_volcano = ElevationLine(phi=0.0, theta=np.zeros(4), elevation=np.zeros(4), is_volcano=np.ones(4, dtype=bool))
    line_none = ElevationLine(phi=0.0, theta=np.zeros(4), elevation=np.zeros(4))
    plate = Plate(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line_all_volcano, line_none])
    assert volcanism._volcano_fraction(plate) == 0.5


def test_apply_volcanic_activity_relabels_a_diluted_field_as_ordinary_continental():
    # 1 volcano node out of 21 total -- well under VOLCANO_FRACTION_DORMANT_THRESHOLD (5%).
    volcano_line = ElevationLine(
        phi=0.0, theta=np.array([0.0]), elevation=np.array([200.0]),
        is_volcano=np.array([True]), volcano_active_years_remaining=np.array([0.0]),
    )
    ordinary_line = ElevationLine(phi=0.1, theta=np.zeros(20), elevation=np.full(20, 200.0))
    plate = Plate(plate_id=5, frame=np.eye(3), crust_type="continental", lines=[volcano_line, ordinary_line])
    world = World(seed=0, plates=[plate], volcanic_field_plate_ids={5})

    events = volcanism.apply_volcanic_activity(world, years=1_000_000)
    assert 5 not in world.volcanic_field_plate_ids
    assert len(events) == 1
    assert "cooled" in events[0].lower()


def test_apply_volcanic_activity_keeps_a_field_tracked_while_still_mostly_volcanic():
    line = ElevationLine(
        phi=0.0, theta=np.zeros(4), elevation=np.full(4, 200.0),
        is_volcano=np.ones(4, dtype=bool), volcano_active_years_remaining=np.full(4, 500_000.0),
    )
    plate = Plate(plate_id=5, frame=np.eye(3), crust_type="continental", lines=[line])
    world = World(seed=0, plates=[plate], volcanic_field_plate_ids={5})

    events = volcanism.apply_volcanic_activity(world, years=1_000_000)
    assert 5 in world.volcanic_field_plate_ids
    assert events == []


def test_apply_volcanic_activity_decrements_and_floors_remaining_active_years():
    line = ElevationLine(
        phi=0.0, theta=np.zeros(4), elevation=np.full(4, 200.0),
        is_volcano=np.ones(4, dtype=bool), volcano_active_years_remaining=np.array([100_000.0, 300_000.0, 1_000.0, 0.0]),
    )
    plate = Plate(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
    world = World(seed=0, plates=[plate])

    volcanism.apply_volcanic_activity(world, years=200_000)
    remaining = world.plates[0].lines[0].volcano_active_years_remaining
    assert remaining.tolist() == [0.0, 100_000.0, 0.0, 0.0]  # clamped at 0, never negative


def test_apply_volcanic_activity_can_erupt_and_add_elevation():
    # A large world of active volcanoes with plenty of remaining life -- over many steps,
    # at ERUPTION_RATE_PER_MYR=3.0/Myr, at least one of them should erupt somewhere.
    n = 200
    line = ElevationLine(
        phi=0.0, theta=np.arange(n) * 0.001, elevation=np.full(n, 200.0),
        is_volcano=np.ones(n, dtype=bool), volcano_active_years_remaining=np.full(n, 1_000_000.0),
    )
    plate = Plate(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
    world = World(seed=0, plates=[plate])

    original_elevation = world.plates[0].lines[0].elevation.copy()
    for _ in range(5):
        volcanism.apply_volcanic_activity(world, years=100_000)
    new_elevation = world.plates[0].lines[0].elevation
    assert np.any(new_elevation > original_elevation)
    assert np.all(new_elevation <= volcanism.MAX_ELEVATION_M)


def test_apply_volcanic_activity_noop_for_empty_world():
    world = World(seed=0, plates=[])
    assert volcanism.apply_volcanic_activity(world, years=1_000_000) == []


def test_volcanism_integration_is_deterministic_for_the_same_seed():
    def run():
        world = generate_world(seed=17, num_plates=10, continental_fraction=0.4)
        for _ in range(10):
            step_world(world, years=2_000_000)
        return sorted(world.volcanic_field_plate_ids), len(world.plates)

    assert run() == run()
