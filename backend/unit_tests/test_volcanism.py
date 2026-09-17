import numpy as np
from app import lithosphere, volcanism
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


def test_apply_volcanic_activity_backs_erupted_elevation_with_crustal_thickness():
    # Issue #173: an eruption's elevation gain must come with a matching crustal_thickness_m
    # bump when the line actually tracks Hc/Hm, so the added relief is isostatically backed
    # instead of "phantom" (unbacked) relief erosion can tear down for free.
    n = 200
    line = ElevationLine(
        phi=0.0, theta=np.arange(n) * 0.001, elevation=np.full(n, 200.0),
        is_volcano=np.ones(n, dtype=bool), volcano_active_years_remaining=np.full(n, 1_000_000.0),
        crustal_thickness_m=np.full(n, 35_000.0), mantle_lithosphere_thickness_m=np.full(n, 100_000.0),
    )
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
    world = World(seed=0, plates=[plate])

    original_hc = world.plates[0].lines[0].crustal_thickness_m.copy()
    original_elevation = world.plates[0].lines[0].elevation.copy()
    for _ in range(5):
        volcanism.apply_volcanic_activity(world, years=100_000)
    new_line = world.plates[0].lines[0]
    erupted = new_line.elevation > original_elevation
    assert np.any(erupted)
    assert np.all(new_line.crustal_thickness_m[erupted] > original_hc[erupted])
    # An untouched node's crust shouldn't move just because its neighbours erupted.
    untouched = ~erupted & (new_line.elevation == original_elevation)
    if np.any(untouched):
        assert np.allclose(new_line.crustal_thickness_m[untouched], original_hc[untouched])


def test_back_elevation_gain_does_not_launder_pre_existing_unbacked_drift_into_crust():
    # Code-review finding on issue #173's fix: solving crustal_thickness_m from a node's
    # *entire* current elevation (rather than the eruption's own incremental gain) would
    # retroactively bake any pre-existing unbacked drift -- e.g. from faults.py, which mutates
    # elevation directly and never touches Hc -- into real crust the moment that node erupts.
    # The Hc bump one eruption produces must depend only on ERUPTION_ELEVATION_M, not on how
    # much unrelated drift the node happened to be carrying beforehand.
    hc, hm = lithosphere.reference_thickness("continental")
    rho_c = lithosphere.RHO_CONTINENTAL_CRUST
    equilibrium = lithosphere.isostatic_elevation(np.array([hc]), np.array([hm]), rho_c)[0]
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[])

    def hc_gain_from_one_eruption(pre_existing_drift_m: float) -> float:
        line = ElevationLine(
            phi=0.0, theta=np.zeros(1), elevation=np.array([equilibrium + pre_existing_drift_m]),
            crustal_thickness_m=np.array([hc]), mantle_lithosphere_thickness_m=np.array([hm]),
        )
        new_hc, _ = volcanism._back_elevation_gain(
            line, plate, volcanism.ERUPTION_ELEVATION_M, np.array([True])
        )
        return float(new_hc[0] - hc)

    gain_no_drift = hc_gain_from_one_eruption(0.0)
    gain_with_drift = hc_gain_from_one_eruption(2000.0)  # e.g. unrelated fault-driven relief
    assert gain_no_drift > 0.0
    assert abs(gain_with_drift - gain_no_drift) < 1e-6


def test_apply_volcanic_activity_falls_back_to_bare_elevation_without_a_crustal_column():
    # v1 lines with no Hc/Hm tracking (crustal_thickness_m defaults to all-zero) keep the old
    # bare direct-elevation response -- same `has_column` compatibility gate erosion.py uses.
    n = 200
    line = ElevationLine(
        phi=0.0, theta=np.arange(n) * 0.001, elevation=np.full(n, 200.0),
        is_volcano=np.ones(n, dtype=bool), volcano_active_years_remaining=np.full(n, 1_000_000.0),
    )
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
    world = World(seed=0, plates=[plate])

    for _ in range(5):
        volcanism.apply_volcanic_activity(world, years=100_000)
    new_line = world.plates[0].lines[0]
    assert np.all(new_line.crustal_thickness_m == 0.0)
    assert np.any(new_line.elevation > 200.0)


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
