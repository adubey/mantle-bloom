import numpy as np
from app import volcanism
from app.plates import ElevationLine, PlateWithLines
from app.world import World


def test_volcano_fraction():
    line_all_volcano = ElevationLine(phi=0.0, theta=np.zeros(4), elevation=np.zeros(4), is_volcano=np.ones(4, dtype=bool))
    line_none = ElevationLine(phi=0.0, theta=np.zeros(4), elevation=np.zeros(4))
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line_all_volcano, line_none])
    assert volcanism._volcano_fraction(plate) == 0.5


def test_apply_volcanic_activity_relabels_a_diluted_field_as_ordinary_continental():
    # 1 volcano node out of 21 total -- well under VOLCANO_FRACTION_DORMANT_THRESHOLD (5%).
    volcano_line = ElevationLine(
        phi=0.0, theta=np.array([0.0]), elevation=np.array([200.0]),
        is_volcano=np.array([True]), volcano_active_years_remaining=np.array([0.0]),
    )
    ordinary_line = ElevationLine(phi=0.1, theta=np.zeros(20), elevation=np.full(20, 200.0))
    plate = PlateWithLines(plate_id=5, frame=np.eye(3), crust_type="continental", lines=[volcano_line, ordinary_line])
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
    plate = PlateWithLines(plate_id=5, frame=np.eye(3), crust_type="continental", lines=[line])
    world = World(seed=0, plates=[plate], volcanic_field_plate_ids={5})

    events = volcanism.apply_volcanic_activity(world, years=1_000_000)
    assert 5 in world.volcanic_field_plate_ids
    assert events == []


def test_apply_volcanic_activity_decrements_and_floors_remaining_active_years():
    line = ElevationLine(
        phi=0.0, theta=np.zeros(4), elevation=np.full(4, 200.0),
        is_volcano=np.ones(4, dtype=bool), volcano_active_years_remaining=np.array([100_000.0, 300_000.0, 1_000.0, 0.0]),
    )
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
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
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
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


def test_apply_volcanic_activity_erupting_grows_mineral_deposit_monotonically():
    # Same setup as test_apply_volcanic_activity_can_erupt_and_add_elevation -- an eruption
    # should also grow mineral_deposit_m, and never let it fall (monotonic, like silt_depth).
    n = 200
    line = ElevationLine(
        phi=0.0, theta=np.arange(n) * 0.001, elevation=np.full(n, 200.0),
        is_volcano=np.ones(n, dtype=bool), volcano_active_years_remaining=np.full(n, 1_000_000.0),
    )
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
    world = World(seed=0, plates=[plate])

    prior = world.plates[0].lines[0].mineral_deposit_m.copy()
    for _ in range(5):
        volcanism.apply_volcanic_activity(world, years=100_000)
        current = world.plates[0].lines[0].mineral_deposit_m
        assert np.all(current >= prior)  # never decreases
        prior = current.copy()
    assert np.any(prior > 0.0)  # at least one node actually erupted somewhere over 5 steps
    assert np.all(prior <= volcanism.MAX_MINERAL_DEPOSIT_M)
