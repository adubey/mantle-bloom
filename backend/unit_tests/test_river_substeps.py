"""The #867 river-substep arms: off by default, bit-identical when trivial, and each arm
reaches the model where it says it does. Every test leaves `river_substeps.ARM` off."""

import numpy as np
import pytest

from app import surface_parity
from app.hydroclimate import breaching, hydrology, river_substeps
from app.world import generate_world, step_world


@pytest.fixture(autouse=True)
def _arm_off():
    yield
    river_substeps.set_arm(river_substeps.RiverSubstepArm())


def _small_world():
    world = generate_world(seed=42, num_plates=2, node_density=0.25, climate_density=0.25, fluid_density=0.25)
    world.simulate_glacier_march = False
    return world


def _run(arm: str, steps: int = 2):
    river_substeps.set_arm(river_substeps.parse_arm(arm))
    world = _small_world()
    for _ in range(steps):
        step_world(world, 100_000)
    return world


def test_default_arm_is_off_and_adds_nothing_to_a_step():
    assert river_substeps.ARM == river_substeps.RiverSubstepArm()
    world = _run("off", steps=1)
    assert not river_substeps.LOG
    assert getattr(world.hydrology_cache, "river_source_m3_per_yr", None) is not None  # the march's inputs are kept either way


def test_parse_arm_round_trips_every_spelling():
    for text, mode, trigger, substeps in [
        ("off", "off", "threatened", 1),
        ("outlet-bound", "outlet_bound", "threatened", 1),
        ("outlet-bound:all", "outlet_bound", "all", 1),
        ("uplift_bound", "uplift_bound", "threatened", 1),
        ("march:10", "march", "threatened", 10),
    ]:
        arm = river_substeps.parse_arm(text)
        assert (arm.mode, arm.trigger, arm.substeps) == (mode, trigger, substeps)
        assert river_substeps.parse_arm(arm.label) == arm
    with pytest.raises(ValueError):
        river_substeps.parse_arm("faster-rivers")
    with pytest.raises(ValueError):
        river_substeps.RiverSubstepArm(mode="march", substeps=0)


def test_exported_to_ocean_follows_chains_and_rejects_sinks_and_cycles():
    # 0 -> 1 -> 2 (ocean); 3 -> 4 (land sink); 5 <-> 6 (a cycle); 2 is the ocean itself.
    flow_target = np.array([1, 2, -1, 4, -1, 6, 5])
    is_ocean = np.array([False, False, True, False, False, False, False])
    assert river_substeps.exported_to_ocean(flow_target, is_ocean).tolist() == [True, True, True, False, False, False, False]
    assert river_substeps.exported_to_ocean(np.zeros(0, dtype=np.int64), np.zeros(0, dtype=bool)).tolist() == []


def test_threatened_mask_needs_a_previous_routing():
    world = _small_world()
    # Generation bootstraps a routing, so even step 1 has a previous one; only a cache-less
    # world (a compact save's first step) has none.
    bootstrap = world.hydrology_cache
    world.hydrology_cache = None
    assert river_substeps.threatened_mask(world, np.zeros((3, 3))) is None
    world.hydrology_cache = bootstrap
    step_world(world, 100_000)
    cache = world.hydrology_cache
    mask = river_substeps.threatened_mask(world, cache.points)
    exported = river_substeps.exported_to_ocean(cache.flow_target, cache.is_ocean)
    # Resampled onto the cache's own points, the mask is the cache's own export mask.
    np.testing.assert_array_equal(mask, exported)
    assert mask.any() and not mask.all()


def test_outlet_bound_relief_names_threatened_pits_only():
    world = _small_world()
    step_world(world, 100_000)
    cache = world.hydrology_cache
    river_substeps.set_arm(river_substeps.parse_arm("outlet-bound"))
    relief = river_substeps.outlet_relief(world, cache.points, cache.elevation, cache.plates_in_order, cache.is_ocean, 100_000)
    assert relief is not None
    assert not relief.eligible[cache.is_ocean].any()
    assert relief.extra_cut_m.min() == pytest.approx(100_000.0)  # 1 m/yr x 100 kyr
    river_substeps.set_arm(river_substeps.parse_arm("outlet-bound:all"))
    everywhere = river_substeps.outlet_relief(world, cache.points, cache.elevation, cache.plates_in_order, cache.is_ocean, 100_000)
    np.testing.assert_array_equal(everywhere.eligible, ~cache.is_ocean)
    assert everywhere.eligible.sum() >= relief.eligible.sum()


def test_uplift_bound_relief_is_the_rise_since_the_last_erosion_pass():
    world = _small_world()
    step_world(world, 100_000)
    cache = world.hydrology_cache
    river_substeps.set_arm(river_substeps.parse_arm("uplift-bound"))
    from app import plates
    from app.surface_fields import CHANNEL_REFERENCE_VALID_BELOW_M

    order = cache.plates_in_order
    current = plates.collect_all_elevation(order)
    valid = plates.collect_all_channel_reference_elevation(order) < CHANNEL_REFERENCE_VALID_BELOW_M
    assert valid.any()
    base = river_substeps.uplift_since_reference_m(current, order)
    relief = river_substeps.outlet_relief(world, cache.points, current + 7.0, order, cache.is_ocean, 100_000)
    # The reference is the last erosion pass's elevation, so a uniform 7 m rise adds 7 m of
    # allowance wherever a reference exists and nothing where none does.
    np.testing.assert_allclose((relief.extra_cut_m - base)[valid], 7.0)
    assert not relief.extra_cut_m[~valid].any()
    lowered = river_substeps.outlet_relief(world, cache.points, current - 7.0, order, cache.is_ocean, 100_000)
    np.testing.assert_allclose(lowered.extra_cut_m, np.clip(base - 7.0, 0.0, None))


def test_a_one_slice_march_is_the_baseline_bit_for_bit():
    baseline = surface_parity.state_hash(_run("off"))
    assert surface_parity.state_hash(_run("march:1")) == baseline
    assert not river_substeps.LOG


def test_the_march_runs_logs_and_keeps_the_ledgers_closed():
    world = _run("march:3")
    entry = river_substeps.LOG[-1]
    assert entry["mode"] == "march" and entry["substeps"] == 3
    assert entry["demand_sum_m"] > 0.0 and entry["seconds"] > 0.0
    fields = world.hydrology_cache
    land = fields.land_water_ledger
    residual = hydrology.full_land_water_residual_m3(land, fields.glacier_water_ledger)
    assert abs(residual) <= 1e-9 * (land.rain_m3 + land.meltwater_m3)
    budget = world.erosion_cache.budget
    assert budget["deposited_m3"] + budget["hc_cap_overflow_terminal_m3"] == pytest.approx(budget["removed_m3"], rel=1e-9)
    # The march's cuts went through the same write-back as the baseline's: channels grew.
    channel = np.concatenate([p.collect("channel_depth") for p in world.plates])
    assert channel.max() > 0.0


def test_the_bound_arms_run_through_a_real_step_and_log_their_pits():
    for arm in ("outlet-bound", "uplift-bound"):
        world = _run(arm, steps=2)
        entry = river_substeps.LOG[-1]
        assert entry["mode"] == arm.replace("-", "_")
        assert 0 <= entry["relieved_pits"] <= entry["breached_pits"] <= entry["pits"]
        assert 0 <= entry["eligible_pits"] <= entry["pits"]
        fields = world.hydrology_cache
        land = fields.land_water_ledger
        assert abs(hydrology.full_land_water_residual_m3(land, fields.glacier_water_ledger)) <= 1e-9 * (land.rain_m3 + land.meltwater_m3)


def test_march_inputs_stay_out_of_saves():
    from app import persistence

    world = _run("off", steps=1)
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert not hasattr(loaded.hydrology_cache, "river_source_m3_per_yr")
    assert not hasattr(loaded.hydrology_cache, "breach_water_balance")
    import pickle

    assert not hasattr(pickle.loads(pickle.dumps(world.hydrology_cache)), "river_source_m3_per_yr")
    # A loaded cache can't march on its first step, and says so instead of failing.
    river_substeps.set_arm(river_substeps.parse_arm("march:2"))
    march = river_substeps.march_rivers(
        loaded, loaded.hydrology_cache, loaded.hydrology_cache.points, loaded.hydrology_cache.elevation,
        np.zeros(len(loaded.hydrology_cache.points)), np.ones(len(loaded.hydrology_cache.points)), np.ones(len(loaded.hydrology_cache.points)),
        np.zeros(len(loaded.hydrology_cache.points)), np.ones(len(loaded.hydrology_cache.points)), np.full(len(loaded.hydrology_cache.points), np.inf),
        100_000, 1.0,
    )
    assert march is None
    assert river_substeps.LOG[-1]["skipped"]
    # A real step computes fresh inputs before erosion asks for the march, so the first step
    # after a load marches like any other.
    step_world(loaded, 100_000)
    assert getattr(loaded.hydrology_cache, "river_source_m3_per_yr", None) is not None
    assert river_substeps.LOG[-1]["substeps"] == 2 and "skipped" not in river_substeps.LOG[-1]


def test_relief_without_a_relieved_pit_leaves_breaching_unchanged():
    elevation = np.array([60.0, 0.0, 50.0, 10.0, -50.0])
    is_ocean = np.array([False, False, False, False, True])
    chain = np.array([[1, 1], [0, 2], [1, 3], [2, 4], [3, 3]])
    rate = np.full(5, breaching.BREACH_REFERENCE_CARVE_M_PER_MYR)
    plain = breaching.breach_depressions(elevation, is_ocean, chain, np.zeros(5), rate, years=100_000)
    relief = breaching.OutletRelief(eligible=np.zeros(5, dtype=bool), extra_cut_m=np.full(5, 1e5))
    same = breaching.breach_depressions(elevation, is_ocean, chain, np.zeros(5), rate, years=100_000, relief=relief)
    assert same.breached_pits.tolist() == plain.breached_pits.tolist() == []
    assert same.closed_pits.tolist() == plain.closed_pits.tolist() == [1]
    assert same.relieved_pits.tolist() == []
