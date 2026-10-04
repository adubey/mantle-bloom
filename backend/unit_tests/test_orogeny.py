"""Orogenic relief before delamination (issue #290): the thermal proxy and its lag, crust
states, anatexis, and collapse / ductile flow on standing quad orogens -- see orogeny.py."""

import numpy as np
import pytest

from app import continental_ledger, cratons, lithosphere, orogeny
from app.elevation_lines import CRUST_TYPE_OCEANIC, line_spacing_rad
from app.sparse_quad_patch import PlateWithSparseQuadPatch, cells_per_face_edge, pack_cell_keys, unpack_cell_keys
from app.world import World

DENSITY = 0.5
SPACING = line_spacing_rad(DENSITY)
N = cells_per_face_edge(SPACING)
REF_HC = lithosphere.REFERENCE_HC_CONTINENTAL_M
REF_HM = lithosphere.REFERENCE_HM_CONTINENTAL_M


def _strip(length: int, hc: np.ndarray | None = None, hm: float = REF_HM, refine=None) -> PlateWithSparseQuadPatch:
    ii = np.arange(10, 10 + length)
    keys = pack_cell_keys(np.zeros(length, dtype=int), ii, np.full(length, 20))
    plate = PlateWithSparseQuadPatch(
        1, np.eye(3), "continental", N, keys,
        fields={"crustal_thickness_m": np.full(length, REF_HC), "mantle_lithosphere_thickness_m": np.full(length, hm)},
    )
    if refine is not None:
        plate.refine_cells(plate.cell_keys[refine])
    if hc is not None:
        plate.set_fields_on_plate(crustal_thickness_m=hc(plate))
    plate.set_fields_on_plate(continental_material_m=plate.collect("crustal_thickness_m"))
    lithosphere.sync_plate_elevation(plate)
    return plate


def _world(plate) -> World:
    world = World(seed=0, plates=[plate], next_plate_id=2, node_density=DENSITY)
    continental_ledger.ensure_initialized(world)
    return world


def _peak(thickness: float, width: int = 2):
    def build(plate):
        hc = np.full(plate.node_count(), REF_HC)
        _, _, i, _ = unpack_cell_keys(plate.cell_keys)
        hc[np.argsort(i, kind="stable")[:width]] = thickness
        return hc
    return build


def _volume(plate, name="crustal_thickness_m") -> float:
    return float(plate.collect(name) @ plate.node_areas_m2())


# --- Thermal proxy and crust states -------------------------------------------------------


def test_moho_heats_with_crustal_thickness_and_cools_under_a_thick_lid():
    reference = orogeny.moho_temperature_c(np.array([REF_HC]), np.array([REF_HM]))[0]
    assert 430.0 < reference < 530.0  # stable continent
    thick = orogeny.moho_temperature_c(np.array([2 * REF_HC]), np.array([REF_HM]))[0]
    thick_lid = orogeny.moho_temperature_c(np.array([2 * REF_HC]), np.array([2 * REF_HM]))[0]
    assert thick > thick_lid > reference


def test_crust_states_separate_cold_hot_melting_and_delaminable_columns():
    hc = np.array([REF_HC, 52_000.0, 57_000.0, 70_000.0, 60_000.0, 70_000.0])
    hm = np.array([REF_HM, REF_HM, REF_HM, REF_HM, 240_000.0, REF_HM])
    continental = np.array([True, True, True, True, True, False])
    state = orogeny.crust_state(hc, hm, continental)
    assert list(state) == [
        orogeny.CRUST_STATE_COLD_STRONG,
        orogeny.CRUST_STATE_HOT_WEAK,
        orogeny.CRUST_STATE_MELT_ELIGIBLE,  # hot, but its root is too thin to founder
        orogeny.CRUST_STATE_DELAMINATION_ELIGIBLE,
        orogeny.CRUST_STATE_COLD_STRONG,  # thick crust on a cold, thick lid
        orogeny.CRUST_STATE_NOT_CONTINENTAL,
    ]


def test_delamination_capacity_is_rate_limited_and_zero_without_elapsed_time():
    hc, hm = np.array([80_000.0]), np.array([REF_HM])
    root = orogeny.delaminable_root_m(hc, hm)[0]
    assert root == pytest.approx(80_000.0 - orogeny.ECLOGITE_DEPTH_M)
    assert orogeny.delamination_capacity_m(hc, hm, 0.0)[0] == 0.0
    one = orogeny.delamination_capacity_m(hc, hm, 1.0)[0]
    assert one == pytest.approx(root * (1.0 - np.exp(-orogeny.DELAMINATION_RATE_PER_MYR)))
    assert orogeny.delamination_capacity_m(hc, hm, 1e6)[0] == pytest.approx(root)


# --- Thermal lag --------------------------------------------------------------------------


def test_thickening_buries_the_moho_with_its_old_temperature_and_the_lag_then_decays():
    hc0, hm0 = np.array([REF_HC]), np.array([REF_HM])
    hc1, hm1 = 2 * hc0, 1.5 * hm0
    before = orogeny.moho_temperature_c(hc0, hm0)[0]
    lag = orogeny.bury_moho(np.zeros(1), hc0, hm0, hc1, hm1)
    assert orogeny.moho_temperature_c(hc1, hm1, lag)[0] == pytest.approx(before)
    assert orogeny.moho_temperature_c(hc1, hm1)[0] > before + 200.0

    decayed = orogeny.relax_thermal_lag(lag, orogeny.THERMAL_RELAXATION_MYR)
    assert decayed[0] == pytest.approx(lag[0] / np.e)
    # Exact decay: the same span in many steps lands in the same place.
    stepped = lag
    for _ in range(40):
        stepped = orogeny.relax_thermal_lag(stepped, 1.0)
    assert stepped[0] == pytest.approx(orogeny.relax_thermal_lag(lag, 40.0)[0], rel=1e-12)


def test_thinning_leaves_the_moho_briefly_hotter_than_its_steady_state():
    hc0, hm0 = np.array([2 * REF_HC]), np.array([REF_HM])
    lag = orogeny.bury_moho(np.zeros(1), hc0, hm0, 0.8 * hc0, hm0)
    assert lag[0] < 0.0
    assert orogeny.moho_temperature_c(0.8 * hc0, hm0, lag)[0] == pytest.approx(orogeny.moho_temperature_c(hc0, hm0)[0])


def test_freshly_thickened_crust_only_melts_and_delaminates_once_it_has_heated():
    hc0, hm0 = np.array([REF_HC]), np.array([REF_HM])
    hc1 = np.array([75_000.0])
    lag = orogeny.bury_moho(np.zeros(1), hc0, hm0, hc1, hm0)
    continental = np.array([True])

    assert orogeny.crust_state(hc1, hm0, continental, lag)[0] == orogeny.CRUST_STATE_COLD_STRONG
    assert orogeny.delaminable_root_m(hc1, hm0, lag)[0] == 0.0
    assert orogeny.anatexis_yield_m(hc1, hm0, lag, np.zeros(1), 1.0)[0][0] == 0.0

    warm = orogeny.relax_thermal_lag(lag, 4 * orogeny.THERMAL_RELAXATION_MYR)
    assert orogeny.crust_state(hc1, hm0, continental, warm)[0] == orogeny.CRUST_STATE_DELAMINATION_ELIGIBLE
    assert orogeny.anatexis_yield_m(hc1, hm0, warm, np.zeros(1), 1.0)[0][0] > 0.0


# --- Anatexis and restite -----------------------------------------------------------------


def test_anatexis_needs_a_hot_moho_and_is_rate_limited():
    hm = np.array([REF_HM, REF_HM])
    hc = np.array([REF_HC, 70_000.0])
    melt, residue = orogeny.anatexis_yield_m(hc, hm, np.zeros(2), np.zeros(2), 1.0)
    assert melt[0] == residue[0] == 0.0
    assert melt[1] > 0.0 and residue[1] > melt[1]
    # Never more than the rate allows of the zone above the solidus, and never more than
    # `MAX_MELT_FRACTION` of what it processes.
    moho = orogeny.moho_temperature_c(hc[1:], hm[1:])[0]
    zone = hc[1] * (moho - orogeny.MELT_ONSET_MOHO_C) / (moho - orogeny.SURFACE_TEMPERATURE_C)
    processed = melt[1] + residue[1]
    assert processed == pytest.approx(zone * (1.0 - np.exp(-orogeny.ANATEXIS_RATE_PER_MYR)))
    assert melt[1] <= orogeny.MAX_MELT_FRACTION * processed
    assert orogeny.anatexis_yield_m(hc, hm, np.zeros(2), np.zeros(2), 0.0)[0][1] == 0.0


def test_restite_cannot_melt_again():
    hc, hm = np.array([70_000.0]), np.array([REF_HM])
    fresh, _ = orogeny.anatexis_yield_m(hc, hm, np.zeros(1), np.zeros(1), 1.0)
    depleted, _ = orogeny.anatexis_yield_m(hc, hm, np.zeros(1), np.array([5_000.0]), 1.0)
    spent, _ = orogeny.anatexis_yield_m(hc, hm, np.zeros(1), hc.copy(), 1.0)
    assert fresh[0] > depleted[0] > 0.0
    assert spent[0] == 0.0


def test_restite_adds_a_dense_root_above_the_eclogite_transition():
    # 56 km of hot crust has only 6 km below 50 km -- too thin to founder -- until melting
    # leaves 12 km of restite at its base.
    hc, hm = np.array([56_000.0]), np.array([80_000.0])
    lag = np.zeros(1)
    assert orogeny.moho_temperature_c(hc, hm)[0] >= orogeny.DELAMINATION_MIN_MOHO_C
    assert orogeny.delaminable_root_m(hc, hm, lag)[0] == 0.0
    assert orogeny.delaminable_root_m(hc, hm, lag, np.array([12_000.0]))[0] == pytest.approx(12_000.0)
    assert orogeny.delaminable_root_m(hc, hm, lag, np.array([3_000.0]))[0] == 0.0


def test_anatexis_books_melt_and_residue_and_conserves_crust_and_provenance():
    plate = _strip(12, _peak(72_000.0, width=3))
    plate.set_fields_on_plate(continental_material_m=0.5 * plate.collect("crustal_thickness_m"))
    world = _world(plate)
    hc_before = plate.collect("crustal_thickness_m")
    volume, material = _volume(plate), _volume(plate, "continental_material_m")

    orogeny.anatexis(plate, world, 1_000_000.0)

    budget = world.orogenic_relief_budget
    assert budget["anatexis_melt_extracted_m3"] > 0.0
    assert budget["anatexis_melt_relaminated_m3"] > 0.0
    assert budget["anatexis_melt_emplaced_m3"] + budget["anatexis_melt_relaminated_m3"] == pytest.approx(
        budget["anatexis_melt_extracted_m3"], rel=1e-12
    )
    assert _volume(plate, "restite_m") == pytest.approx(budget["anatexis_residue_m3"], rel=1e-12)
    assert _volume(plate) == pytest.approx(volume, rel=1e-12)
    assert _volume(plate, "continental_material_m") == pytest.approx(material, rel=1e-12)
    hc = plate.collect("crustal_thickness_m")
    # Only the plateau's edge column has a thinner neighbour to relaminate into, and only it
    # sends `ANATEXIS_RELAMINATION_FRACTION` of its melt there; the interior keeps its melt.
    melt, _ = orogeny.anatexis_yield_m(hc_before, plate.collect("mantle_lithosphere_thickness_m"), np.zeros(12), np.zeros(12), 1.0)
    edge_melt = float(melt[2] * plate.node_areas_m2()[2])
    assert budget["anatexis_melt_relaminated_m3"] == pytest.approx(orogeny.ANATEXIS_RELAMINATION_FRACTION * edge_melt, rel=1e-9)
    assert np.array_equal(hc[:2], hc_before[:2])
    assert hc[2] < hc_before[2] and hc[3] > hc_before[3]
    assert np.array_equal(hc[4:], hc_before[4:])
    continental_ledger.assert_closed(world)


def test_relaminated_melt_never_enters_oceanic_cells():
    def build(plate):
        hc = np.full(plate.node_count(), 72_000.0)
        hc[4] = 7_000.0
        return hc

    plate = _strip(5, build)
    codes = plate.collect("crust_type_code")
    codes[4] = CRUST_TYPE_OCEANIC
    plate.set_fields_on_plate(crust_type_code=codes)
    world = _world(plate)

    orogeny.anatexis(plate, world, 1_000_000.0)

    assert plate.collect("crustal_thickness_m")[4] == 7_000.0
    assert plate.collect("restite_m")[4] == 0.0
    budget = world.orogenic_relief_budget
    assert budget["anatexis_melt_extracted_m3"] > 0.0
    assert budget["anatexis_melt_relaminated_m3"] == 0.0


def test_evolving_standing_orogens_relaxes_the_lag():
    plate = _strip(4)
    world = _world(plate)
    plate.set_fields_on_plate(moho_thermal_lag_c=np.full(4, 200.0))

    orogeny.evolve_standing_orogens(plate, world, 10_000_000.0)

    np.testing.assert_allclose(plate.collect("moho_thermal_lag_c"), 200.0 * np.exp(-10.0 / orogeny.THERMAL_RELAXATION_MYR))


# --- Collapse and ductile flow ------------------------------------------------------------


def test_a_saturated_core_collapses_into_its_neighbours_conservatively():
    plate = _strip(12, _peak(lithosphere.MAX_CRUSTAL_THICKNESS_M))
    world = _world(plate)
    hc_before = plate.collect("crustal_thickness_m")
    volume, material = _volume(plate), _volume(plate, "continental_material_m")
    elevation_before = plate.collect("elevation")

    orogeny.relax_orogens(plate, world, 1_000_000.0)

    hc = plate.collect("crustal_thickness_m")
    assert hc[:2].max() < hc_before[:2].max()
    assert hc[2] > hc_before[2]
    assert _volume(plate) == pytest.approx(volume, rel=1e-12)
    assert _volume(plate, "continental_material_m") == pytest.approx(material, rel=1e-12)
    assert np.all(hc[:2] >= orogeny.COLLAPSE_ONSET_HC_M - 1e-6)  # only the excess is mobile
    assert plate.collect("elevation")[0] < elevation_before[0]
    budget = world.orogenic_relief_budget
    moved = budget["collapse_transferred_m3"] + budget["ductile_flow_transferred_m3"]
    assert moved > 0.0
    continental_ledger.assert_closed(world)


def test_crust_below_the_collapse_onset_does_not_flow():
    plate = _strip(8, _peak(orogeny.COLLAPSE_ONSET_HC_M))
    world = _world(plate)
    before = plate.collect("crustal_thickness_m")

    orogeny.relax_orogens(plate, world, 10_000_000.0)

    assert np.array_equal(plate.collect("crustal_thickness_m"), before)


def test_a_lagging_moho_keeps_thickened_crust_from_flowing_ductilely():
    def ductile(lag: float) -> float:
        plate = _strip(12, _peak(70_000.0))
        plate.set_fields_on_plate(moho_thermal_lag_c=np.full(plate.node_count(), lag))
        world = _world(plate)
        orogeny.relax_orogens(plate, world, 1_000_000.0)
        return world.orogenic_relief_budget["ductile_flow_transferred_m3"]

    assert ductile(0.0) > 0.0
    assert ductile(400.0) == 0.0


def test_hot_crust_flows_faster_than_cold_crust_of_the_same_thickness():
    def moved(hm: float) -> tuple[float, float]:
        plate = _strip(12, _peak(60_000.0), hm=hm)
        world = _world(plate)
        orogeny.relax_orogens(plate, world, 1_000_000.0)
        budget = world.orogenic_relief_budget
        return budget["collapse_transferred_m3"], budget["ductile_flow_transferred_m3"]

    cold_collapse, cold_ductile = moved(lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M)
    hot_collapse, hot_ductile = moved(60_000.0)
    assert cold_ductile == 0.0 and cold_collapse > 0.0  # cold: collapse and faulting only
    assert hot_ductile > 0.0
    assert hot_collapse + hot_ductile > 1.5 * cold_collapse


def test_relaxation_does_not_depend_on_the_step_size():
    def run(steps: int) -> np.ndarray:
        plate = _strip(16, _peak(lithosphere.MAX_CRUSTAL_THICKNESS_M, 3))
        world = _world(plate)
        for _ in range(steps):
            orogeny.relax_orogens(plate, world, 4_000_000.0 / steps)
        return plate.collect("crustal_thickness_m")

    one, four, sixteen = run(1), run(4), run(16)
    assert np.max(np.abs(one - sixteen)) < 0.02 * (lithosphere.MAX_CRUSTAL_THICKNESS_M - REF_HC)
    assert np.max(np.abs(four - sixteen)) < 0.01 * (lithosphere.MAX_CRUSTAL_THICKNESS_M - REF_HC)


def test_relaxation_conserves_volume_across_unequal_area_cells():
    plate = _strip(12, _peak(lithosphere.MAX_CRUSTAL_THICKNESS_M, 3), refine=[2, 3])
    world = _world(plate)
    areas = plate.node_areas_m2()
    assert np.ptp(areas) > 0.5 * areas.max()
    volume, material = _volume(plate), _volume(plate, "continental_material_m")

    orogeny.relax_orogens(plate, world, 2_000_000.0)

    assert _volume(plate) == pytest.approx(volume, rel=1e-12)
    assert _volume(plate, "continental_material_m") == pytest.approx(material, rel=1e-12)
    continental_ledger.assert_closed(world)


def test_crust_never_flows_into_an_oceanic_column():
    plate = _strip(6, _peak(lithosphere.MAX_CRUSTAL_THICKNESS_M, 1))
    codes = plate.collect("crust_type_code")
    codes[1] = CRUST_TYPE_OCEANIC
    hc = plate.collect("crustal_thickness_m")
    hc[1] = lithosphere.REFERENCE_HC_OCEANIC_M
    plate.set_fields_on_plate(crust_type_code=codes, crustal_thickness_m=hc, continental_material_m=np.where(codes == CRUST_TYPE_OCEANIC, 0.0, hc))
    world = _world(plate)
    before = plate.collect("crustal_thickness_m")

    orogeny.relax_orogens(plate, world, 5_000_000.0)

    assert np.array_equal(plate.collect("crustal_thickness_m"), before)


def test_cratonic_crust_that_flows_is_booked_as_collision_reworked():
    plate = _strip(8, _peak(70_000.0, 1))
    craton = np.zeros(plate.node_count())
    craton[0] = 10_000.0  # weak enough (40% strength) to flow somewhat
    plate.set_fields_on_plate(craton_crust_m=craton)
    world = _world(plate)
    cratons.ensure_ledger(world)
    areas = plate.node_areas_m2()
    craton_before = float(craton @ areas)

    orogeny.relax_orogens(plate, world, 2_000_000.0)

    craton_after = float(plate.collect("craton_crust_m") @ areas)
    assert craton_after < craton_before
    assert world.craton_ledger["collision_reworked_m3"] == pytest.approx(craton_before - craton_after)


# --- Competition with suture accretion ----------------------------------------------------


def _full_belt_strip():
    from app import quad_tectonics
    from app.lithosphere_plate import SUTURE_ACCRETION_MAX_HC_M

    def build(plate):
        hc = np.full(plate.node_count(), REF_HC)
        hc[1 : quad_tectonics.SUTURE_ACCRETION_MAX_HOPS + 1] = SUTURE_ACCRETION_MAX_HC_M - 500.0
        hc[0] = 30_000.0
        return hc

    plate = _strip(40, build)
    donors = np.zeros(plate.node_count(), dtype=bool)
    donors[0] = True
    return plate, donors


def test_collapse_frees_belt_room_so_less_suture_crust_delaminates():
    from app import quad_tectonics

    def accrete(relax_first_myr: float) -> dict[str, float]:
        plate, donors = _full_belt_strip()
        world = _world(plate)
        if relax_first_myr:
            orogeny.relax_orogens(plate, world, relax_first_myr * 1e6)
            world.orogenic_relief_budget.clear()
        quad_tectonics._accrete_onto_survivors(plate, donors, ~donors, world, years=1_000_000.0)
        return orogeny.ensure_budget(world)

    direct, relaxed = accrete(0.0), accrete(20.0)
    assert direct["delamination_completed_m3"] > 0.0
    assert relaxed["suture_belt_placed_m3"] > direct["suture_belt_placed_m3"]
    assert relaxed["delamination_completed_m3"] < direct["delamination_completed_m3"]


def test_two_fronts_sharing_a_belt_share_its_root_capacity():
    from app import quad_tectonics
    from app.lithosphere_plate import SUTURE_ACCRETION_MAX_HC_M

    def build(plate):
        hc = np.full(plate.node_count(), SUTURE_ACCRETION_MAX_HC_M)
        hc[0] = hc[-1] = 80_000.0
        return hc

    plate = _strip(10, build)
    world = _world(plate)
    donors = np.zeros(plate.node_count(), dtype=bool)
    donors[[0, -1]] = True
    hc = plate.collect("crustal_thickness_m")
    areas = plate.node_areas_m2()
    capacity = float(orogeny.delamination_capacity_m(hc, plate.collect("mantle_lithosphere_thickness_m"), 1.0)[~donors] @ areas[~donors])

    quad_tectonics._accrete_onto_survivors(plate, donors, ~donors, world, years=1_000_000.0)

    budget = world.orogenic_relief_budget
    assert budget["delamination_completed_m3"] == pytest.approx(capacity, rel=1e-9)
    donated = float(hc[donors] @ areas[donors])
    assert budget["no_outlet_delaminated_m3"] == pytest.approx(donated - capacity, rel=1e-9)


def test_stepping_a_quad_world_books_finite_relief_budgets():
    from app.world import generate_world, step_world

    world = generate_world(seed=3, node_density=0.25)
    for _ in range(3):
        step_world(world, 1_000_000.0)
    budget = world.orogenic_relief_budget
    assert set(orogeny.BUDGET_ACCOUNTS) <= set(budget)
    assert all(np.isfinite(v) and v >= 0.0 for v in budget.values())
    assert budget["suture_donated_m3"] >= budget["suture_belt_placed_m3"] + budget["escape_placed_m3"]
