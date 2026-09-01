import numpy as np
from app import volcanism
from app.plates import ElevationLine, PlateWithLines
from app.world import World


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
    assert volcanism.apply_volcanic_activity(world, years=1_000_000) is None


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
