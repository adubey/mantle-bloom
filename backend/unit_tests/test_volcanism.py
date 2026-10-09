import numpy as np
import pytest
from app import continental_ledger, elevation_lines, lithosphere, volcanism
from app.world import World

from .quad_fixtures import quad_plate


def _volcano_field(**fields):
    """200 active, long-lived volcano cells -- enough that some erupt within a few steps."""
    return quad_plate(
        0, "continental", columns=range(20), rows=range(10), elevation=200.0, is_volcano=True,
        volcano_active_years_remaining=1_000_000.0, **fields,
    )


def test_apply_volcanic_activity_decrements_and_floors_remaining_active_years():
    plate = quad_plate(
        0, "continental", columns=range(4), elevation=200.0,
        is_volcano=True, volcano_active_years_remaining=np.array([100_000.0, 300_000.0, 1_000.0, 0.0]),
    )
    world = World(seed=0, plates=[plate])

    volcanism.apply_volcanic_activity(world, years=200_000)
    remaining = world.plates[0].collect("volcano_active_years_remaining")
    assert remaining.tolist() == [0.0, 100_000.0, 0.0, 0.0]  # clamped at 0, never negative


def test_apply_volcanic_activity_can_erupt_and_add_elevation():
    # A large world of active volcanoes with plenty of remaining life -- over many steps,
    # at ERUPTION_RATE_PER_MYR=3.0/Myr, at least one of them should erupt somewhere.
    world = World(seed=0, plates=[_volcano_field()])

    original_elevation = world.plates[0].collect("elevation").copy()
    for _ in range(5):
        volcanism.apply_volcanic_activity(world, years=100_000)
    new_elevation = world.plates[0].collect("elevation")
    assert np.any(new_elevation > original_elevation)
    assert np.all(new_elevation <= elevation_lines.MAX_ELEVATION_M)


def test_apply_volcanic_activity_backs_erupted_elevation_with_crustal_thickness():
    # Issue #173: an eruption's elevation gain must come with a matching crustal_thickness_m
    # bump when the node actually tracks Hc/Hm, so the added relief is isostatically backed
    # instead of "phantom" (unbacked) relief erosion can tear down for free.
    plate = _volcano_field(crustal_thickness_m=35_000.0, mantle_lithosphere_thickness_m=100_000.0)
    world = World(seed=0, plates=[plate])
    continental_ledger.ensure_initialized(world)
    surface_before = continental_ledger.surface_volume_m3(world)

    original_hc = plate.collect("crustal_thickness_m").copy()
    original_elevation = plate.collect("elevation").copy()
    for _ in range(5):
        volcanism.apply_volcanic_activity(world, years=100_000)
    elevation, hc = plate.collect("elevation"), plate.collect("crustal_thickness_m")
    erupted = elevation > original_elevation
    assert np.any(erupted)
    assert np.all(hc[erupted] > original_hc[erupted])
    # An untouched node's crust shouldn't move just because its neighbours erupted.
    untouched = ~erupted & (elevation == original_elevation)
    if np.any(untouched):
        assert np.allclose(hc[untouched], original_hc[untouched])
    surface_gain = continental_ledger.surface_volume_m3(world) - surface_before
    assert world.continental_material_ledger["juvenile_additions_m3"] == pytest.approx(surface_gain)
    continental_ledger.assert_closed(world)


def test_back_elevation_gain_does_not_launder_pre_existing_unbacked_drift_into_crust():
    # Code-review finding on issue #173's fix: solving crustal_thickness_m from a node's
    # *entire* current elevation (rather than the eruption's own incremental gain) would
    # retroactively bake any pre-existing unbacked drift -- e.g. from deform_columns'
    # transform_uplift, kept as a bare elevation delta by design -- into
    # real crust the moment that node erupts. The Hc bump one eruption produces must depend only
    # on ERUPTION_ELEVATION_M, not on how much unrelated drift the node happened to be carrying
    # beforehand.
    hc, hm = lithosphere.reference_thickness("continental")
    rho_c = lithosphere.RHO_CONTINENTAL_CRUST
    equilibrium = lithosphere.isostatic_elevation(np.array([hc]), np.array([hm]), rho_c)[0]

    def hc_gain_from_one_eruption(pre_existing_drift_m: float) -> float:
        new_hc, _ = lithosphere.back_elevation_gain_fields(
            np.array([equilibrium + pre_existing_drift_m]), np.array([hc]), np.array([hm]), np.zeros(1, dtype=np.int8),
            "continental", volcanism.ERUPTION_ELEVATION_M, np.array([True]),
        )
        return float(new_hc[0] - hc)

    gain_no_drift = hc_gain_from_one_eruption(0.0)
    gain_with_drift = hc_gain_from_one_eruption(2000.0)  # e.g. unrelated fault-driven relief
    assert gain_no_drift > 0.0
    assert abs(gain_with_drift - gain_no_drift) < 1e-6


def test_apply_volcanic_activity_falls_back_to_bare_elevation_without_a_crustal_column():
    # Nodes with no Hc/Hm column (crustal_thickness_m zero) keep the bare direct-elevation
    # response -- same `has_column` gate erosion.py uses.
    plate = _volcano_field(crustal_thickness_m=0.0, mantle_lithosphere_thickness_m=0.0)
    world = World(seed=0, plates=[plate])

    for _ in range(5):
        volcanism.apply_volcanic_activity(world, years=100_000)
    assert np.all(plate.collect("crustal_thickness_m") == 0.0)
    assert np.any(plate.collect("elevation") > 200.0)


def test_apply_volcanic_activity_noop_for_empty_world():
    world = World(seed=0, plates=[])
    assert volcanism.apply_volcanic_activity(world, years=1_000_000) is None


def test_apply_volcanic_activity_erupting_grows_mineral_deposit_monotonically():
    # Same setup as test_apply_volcanic_activity_can_erupt_and_add_elevation -- an eruption
    # should also grow mineral_deposit_m, and never let it fall (monotonic, like coal_deposit_m).
    world = World(seed=0, plates=[_volcano_field()])

    prior = world.plates[0].collect("mineral_deposit_m").copy()
    for _ in range(5):
        volcanism.apply_volcanic_activity(world, years=100_000)
        current = world.plates[0].collect("mineral_deposit_m")
        assert np.all(current >= prior)  # never decreases
        prior = current.copy()
    assert np.any(prior > 0.0)  # at least one node actually erupted somewhere over 5 steps
    assert np.all(prior <= volcanism.MAX_MINERAL_DEPOSIT_M)
