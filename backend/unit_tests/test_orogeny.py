"""Orogenic relief before delamination (issue #290): the thermal proxy, crust states, and
collapse / ductile flow on standing quad orogens -- see orogeny.py."""

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
