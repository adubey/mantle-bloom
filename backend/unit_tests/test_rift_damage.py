import numpy as np
import pytest

from app import cratons, lithosphere, merge_split, persistence, rift_damage
from app import rift_lines as lines
from app.lithosphere_plate import K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION
from app.world import generate_world, step_world

MYR = 1e6
# Easy-to-rift knobs, so a small test world strains within one step.
EAGER = {"k_cap": 0.05, "eps_dot0": 3e-7}


@pytest.fixture(scope="module")
def saved_world() -> bytes:
    world = generate_world(seed=3, num_plates=4, node_density=0.25, climate_density=0.5, fluid_density=0.5)
    world.simulate_climate_biomes = False
    return persistence.save_world_bytes(world)


def _rift_world(saved: bytes):
    world = persistence.load_world_bytes(saved)
    world.split_method, world.rift_line_params = "rift_lines", dict(EAGER)
    return world


def _damage(world, name: str) -> np.ndarray:
    return np.concatenate([plate.collect(name) for plate in world.plates])


def test_the_default_split_method_never_runs_the_rift_lines(saved_world):
    world = persistence.load_world_bytes(saved_world)
    assert world.split_method == "great_circle" and world.rift_line_params == {}
    step_world(world, MYR)
    assert world.rift_runtime is None
    assert not _damage(world, "rift_strain").any() and not _damage(world, "rift_broken").any()


def test_a_damage_pass_strains_a_line_on_each_pulled_plate_and_conserves_every_column(saved_world):
    world = _rift_world(saved_world)
    before = {name: sum(float(np.dot(p.collect(name), p.node_areas_m2())) for p in world.plates)
              for name in ("crustal_thickness_m", "mantle_lithosphere_thickness_m", "continental_material_m")}
    rift_damage.damage_pass(world, MYR)
    state = world.rift_runtime
    assert set(state.steps) == {plate.plate_id for plate in world.plates}  # every plate, not only the one that may split
    assert state.counters["growing_steps"] >= 1 and state.counters["line_steps"] + state.counters["line_none"] == len(world.plates)
    eps, warm = _damage(world, "rift_strain"), _damage(world, "rift_warm_strain")
    assert eps.max() > 0.0 and np.all(eps >= 0.0) and np.all(warm <= eps + 1e-12)
    for step in state.steps.values():
        if step.line is not None:  # strain only on the line's band
            assert step.d_eps[step.line.chain].sum() > 0.0 or not step.d_eps.any()
    for name, total in before.items():
        after = sum(float(np.dot(p.collect(name), p.node_areas_m2())) for p in world.plates)
        assert after == pytest.approx(total, rel=1e-12)
    assert max(state.ledger_gap.values()) < 1e-9


def test_strained_cells_thin_and_subside(saved_world):
    world = _rift_world(saved_world)
    hc0 = {p.plate_id: p.collect("crustal_thickness_m") for p in world.plates}
    z0 = {p.plate_id: p.collect("elevation") for p in world.plates}
    rift_damage.damage_pass(world, MYR)
    thinned = 0
    for plate in world.plates:
        d_eps = world.rift_runtime.steps[plate.plate_id].d_eps
        # Cells far ahead of a rift's tips strain by less than a float can resolve in a thickness; skip those.
        donors = (d_eps > 1e-6) & (hc0[plate.plate_id] > 5_001.0)
        thinned += int(donors.sum())
        assert np.all(plate.collect("crustal_thickness_m")[donors] < hc0[plate.plate_id][donors])
        assert np.all(plate.collect("elevation")[donors] <= z0[plate.plate_id][donors] + 1e-6)
    assert thinned > 0


def test_plume_heat_that_outlasts_the_step_still_counts_when_the_line_is_judged(saved_world):
    # The split test reads the line at the end of the step; heat lasting past it must still load and weaken the line
    # there, and heat that ran out inside the step must not (PR #950 review).
    world = _rift_world(saved_world)
    p = rift_damage.runtime(world).params
    for heat_myr, on_at_end in ((5.0, True), (0.5, False)):
        for plate in world.plates:
            plate.set_fields_on_plate(plume_heat_myr=np.full(plate.node_count(), heat_myr))
        preps = [rift_damage.prepare_plate(world, plate, p, MYR) for plate in world.plates]
        prep = next(prep for prep in preps if prep.line is not None)
        assert np.all(prep.plume_left_yr == heat_myr * MYR)
        prob, _ = rift_damage._line_problem(prep, prep.line, p)
        plume_on = lines.line_state(prob, prep.eps, prep.warm, MYR, p)[2]
        assert bool(plume_on.all()) is on_at_end and bool(plume_on.any()) is on_at_end


def test_a_world_step_with_rift_lines_on_keeps_the_ledgers_and_the_damage_invariants(saved_world):
    world = _rift_world(saved_world)
    step_world(world, MYR)
    assert world.rift_runtime is not None and world.rift_runtime.steps
    eps, warm = _damage(world, "rift_strain"), _damage(world, "rift_warm_strain")
    assert eps.max() > 0.0 and np.all(warm <= eps + 1e-12)  # the remap through the step kept warm strain under strain
    hc = np.concatenate([plate.collect("crustal_thickness_m") for plate in world.plates])
    assert hc.max() <= lithosphere.MAX_CRUSTAL_THICKNESS_M * (1 + 1e-9)
    scale = max(world.craton_ledger["initial_m3"], 1e12)
    assert cratons.balance_error_m3(world) == pytest.approx(0.0, abs=1e-9 * scale)  # rifted_stretch_m3 nets out in the audit
    assert max(world.rift_runtime.ledger_gap.values()) < 1e-9


def _breakable(world):
    """A plate whose rift line gives a valid two-daughter split, with its step record."""
    state = rift_damage.runtime(world)
    for plate in world.plates:
        step = state.steps[plate.plate_id]
        if step.line is not None and lines.split_masks(step.graph, step.areas, step.line.chain, state.params, points=step.points)[0] is not None:
            return plate, step
    raise AssertionError("no plate of the test world has a valid line")


def test_a_broken_line_under_pull_splits_the_plate_in_two(saved_world):
    world = _rift_world(saved_world)
    rift_damage.damage_pass(world, MYR)
    plate, step = _breakable(world)
    nodes, next_id = plate.node_count(), world.next_plate_id
    step.broken_share, step.pull_end, step.strength_end = 1.0, 1.0, 0.0
    daughters = rift_damage.breakthrough(world, plate)
    assert daughters is not None and world.next_plate_id == next_id + 1
    a, b = daughters
    assert {a.plate_id, b.plate_id} == {plate.plate_id, next_id} and a.node_count() + b.node_count() == nodes
    assert min(a.node_count(), b.node_count()) >= 0.1 * nodes  # no microplate (a_min is 15% of area)
    assert a.age_years == 0.0 and b.age_years == 0.0
    event = world.rift_runtime.events[-1]
    assert event["plate_id"] == plate.plate_id and event["small_share"] >= 0.15 and event["path_cells"] == len(step.line.chain)
    assert world.craton_ledger["split_cuts"] == 1.0 and world.craton_ledger["split_band_area_m2"] > 0.0


def test_a_line_splits_only_when_enough_of_it_is_broken_and_the_pull_exceeds_the_strength_left(saved_world):
    world = _rift_world(saved_world)
    rift_damage.damage_pass(world, MYR)
    plate, step = _breakable(world)
    f_break = world.rift_runtime.params.f_break
    step.broken_share, step.pull_end, step.strength_end = f_break - 0.01, 1.0, 0.0
    assert rift_damage.breakthrough(world, plate) is None
    step.broken_share, step.pull_end, step.strength_end = 1.0, 1.0, 2.0
    assert rift_damage.breakthrough(world, plate) is None
    step.pull_end, step.strength_end = 0.0, 0.0
    assert rift_damage.breakthrough(world, plate) is None
    assert world.rift_runtime.events == []


def test_the_rift_line_split_keeps_the_size_and_cooldown_gates(saved_world):
    world = _rift_world(saved_world)
    rift_damage.damage_pass(world, MYR)
    plate, step = _breakable(world)
    step.broken_share, step.pull_end, step.strength_end = 1.0, 1.0, 0.0
    assert plate.age_years < merge_split.SPLIT_MIN_AGE_YEARS
    assert merge_split.maybe_split_plate(world, plate, MYR) is None  # too young, though its line has broken


def test_erupted_cells_lose_their_damage_and_strained_cells_keep_it():
    hc = np.array([lithosphere.REFERENCE_HC_CONTINENTAL_M, lithosphere.REFERENCE_HC_OCEANIC_M, 20_000.0, lithosphere.REFERENCE_HC_CONTINENTAL_M])
    hm = np.array([lithosphere.REFERENCE_HM_CONTINENTAL_M, lithosphere.YOUNG_RIDGE_HM_M, 60_000.0, lithosphere.REFERENCE_HM_CONTINENTAL_M])
    eps, warm, broken = np.array([1.2, 0.4, 0.7, 0.0]), np.array([0.3, 0.4, 0.2, 0.0]), np.array([True, False, True, False])
    rift_damage.reset_erupted(hc, hm, eps, warm, broken)
    assert eps.tolist() == [0.0, 0.0, 0.7, 0.0] and warm.tolist() == [0.0, 0.0, 0.2, 0.0] and broken.tolist() == [False, False, True, False]


def test_a_new_margin_cell_inherits_its_donors_strain_times_its_stretched_share(saved_world):
    world = _rift_world(saved_world)
    plate = max(world.plates, key=lambda p: p.node_count())
    n = plate.node_count()
    rifted = np.array([0, 1])
    assert rift_damage.margin_donor_damage(plate, rifted, rifted, K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION) is None  # no damage yet
    plate.set_fields_on_plate(rift_strain=np.full(n, 0.8), rift_warm_strain=np.full(n, 0.4), rift_broken=np.ones(n, dtype=bool))
    donors = rift_damage.margin_donor_damage(plate, rifted, rifted, K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION)
    assert np.allclose(donors[0], 0.8) and np.allclose(donors[1], 0.4)
    rift_damage.inherit_margin_damage(plate, rifted, np.array([0.5, 0.0]), donors)
    eps, warm, broken = plate.collect("rift_strain"), plate.collect("rift_warm_strain"), plate.collect("rift_broken")
    assert eps[:2].tolist() == pytest.approx([0.4, 0.0]) and warm[:2].tolist() == pytest.approx([0.2, 0.0])
    assert not broken[:2].any() and broken[2:].all() and np.all(eps[2:] == 0.8)


def test_a_mistyped_knob_stops_the_step_instead_of_leaving_the_default(saved_world):
    world = _rift_world(saved_world)
    world.rift_line_params = {"k_cup": 0.1}
    with pytest.raises(ValueError, match="k_cup"):
        rift_damage.damage_pass(world, MYR)


def test_changing_the_knobs_between_steps_takes_effect(saved_world):
    world = _rift_world(saved_world)
    assert rift_damage.runtime(world).params.k_cap == EAGER["k_cap"]
    world.rift_line_params = {"k_cap": 0.5}
    assert rift_damage.runtime(world).params.k_cap == 0.5


def test_the_switch_knobs_and_damage_round_trip_a_save_and_the_scratch_does_not(saved_world):
    world = _rift_world(saved_world)
    rift_damage.damage_pass(world, MYR)
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.split_method == "rift_lines" and loaded.rift_line_params == EAGER and loaded.rift_runtime is None
    for name in ("rift_strain", "rift_warm_strain", "rift_broken"):
        assert np.array_equal(_damage(loaded, name), _damage(world, name))
    assert _damage(loaded, "rift_strain").max() > 0.0


def test_a_save_from_before_rift_lines_loads_on_the_great_circle_method_with_no_damage(saved_world):
    world = persistence.load_world_bytes(saved_world)
    for name in ("split_method", "rift_line_params", "rift_runtime"):
        world.__dict__.pop(name, None)
    for plate in world.plates:
        for name in ("rift_strain", "rift_warm_strain", "rift_broken"):
            plate._fields.pop(name, None)
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.split_method == "great_circle" and loaded.rift_line_params == {} and loaded.rift_runtime is None
    assert not _damage(loaded, "rift_strain").any() and not _damage(loaded, "rift_broken").any()
