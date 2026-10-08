import numpy as np
import pytest

from app import elevation_lines, hm_ledger, lithosphere, phase_budget


class _World:
    debug_diagnostics = True
    node_density = 1.0

    def __init__(self):
        self.phase_budget = {}
        hm_ledger.reset(self)


class _Plate:
    crust_type = "continental"


INHERIT = elevation_lines.CRUST_TYPE_INHERIT


def _record(world, phase, before, after, area_before, area_after):
    phase_budget.record(
        world,
        _Plate(),
        phase,
        np.full(len(before), 35_000.0),
        np.asarray(before, dtype=float),
        np.full(len(before), INHERIT),
        np.full(len(after), 35_000.0),
        np.asarray(after, dtype=float),
        np.full(len(after), INHERIT),
        area_before_m2=np.asarray(area_before, dtype=float),
        area_after_m2=np.asarray(area_after, dtype=float),
    )


def test_positive_source_uses_each_quad_cells_unequal_area():
    world = _World()
    _record(world, "convergent_deformation", [10.0, 10.0], [12.0, 13.0], [2.0, 5.0], [2.0, 5.0])
    scope = world.hm_source_sink_ledger["convergent_strain"]["scopes"]["all"]
    assert scope["source_m3"] == pytest.approx(2.0 * 2.0 + 3.0 * 5.0)
    assert scope["sink_m3"] == 0.0


def test_sink_uses_each_quad_cells_unequal_area():
    world = _World()
    _record(world, "oceanic_cooling_relaxation", [20.0, 30.0], [18.0, 25.0], [3.0, 7.0], [3.0, 7.0])
    scope = world.hm_source_sink_ledger["relaxation"]["scopes"]["all"]
    assert scope["source_m3"] == 0.0
    assert scope["sink_m3"] == pytest.approx(2.0 * 3.0 + 5.0 * 7.0)


def test_cap_clamp_books_sink_and_actual_new_cap_entry_not_net_inference():
    world = _World()
    cap = lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M
    _record(world, "column_cap_clamp", [cap - 2.0, cap + 10.0], [cap, cap], [11.0, 13.0], [11.0, 13.0])
    scope = world.hm_source_sink_ledger["floor_and_cap_clamps"]["scopes"]["all"]
    assert scope["source_m3"] == pytest.approx(2.0 * 11.0)
    assert scope["sink_m3"] == pytest.approx(10.0 * 13.0)
    cap_scope = world.phase_budget["column_cap_clamp"]["hm_cap_transitions"]["scopes"]["all"]
    assert cap_scope["newly_capped_area_m2"] == pytest.approx(11.0)
    assert cap_scope["persistently_capped_area_m2"] == pytest.approx(13.0)


def test_topology_area_change_books_exact_inventory_difference():
    world = _World()
    # The thickness grows, but the represented area shrinks enough that live volume falls:
    # 100 * 2 = 200 before, 150 * 1 = 150 after.
    _record(world, "plate_merge", [100.0], [150.0], [2.0], [1.0])
    scope = world.hm_source_sink_ledger["topology_and_regridding"]["scopes"]["all"]
    assert scope["source_m3"] == 0.0
    assert scope["sink_m3"] == pytest.approx(50.0)


def test_aligned_node_type_change_is_explicit_reclassification():
    world = _World()
    hm_ledger.record_change(
        world,
        _Plate(),
        "convergent_deformation",
        np.array([40_000.0]),
        np.array([elevation_lines.CRUST_TYPE_OCEANIC]),
        np.array([40_000.0]),
        np.array([elevation_lines.CRUST_TYPE_CONTINENTAL]),
        np.array([12.0]),
        np.array([12.0]),
        node_ids_before=np.array([7]),
        node_ids_after=np.array([7]),
    )
    accounts = hm_ledger.cumulative_scopes(world)
    volume = 40_000.0 * 12.0
    assert accounts["all"]["source_m3"] == 0.0
    assert accounts["all"]["sink_m3"] == 0.0
    assert accounts["continental_node"]["reclassification_m3"] == pytest.approx(volume)
    assert accounts["oceanic_node"]["reclassification_m3"] == pytest.approx(-volume)


def test_inherited_node_type_uses_the_before_and_after_plate_defaults():
    world = _World()
    phase_budget.record(
        world,
        _Plate(),
        "plate_merge",
        np.array([35_000.0]),
        np.array([40_000.0]),
        np.array([INHERIT]),
        np.array([35_000.0]),
        np.array([40_000.0]),
        np.array([INHERIT]),
        area_before_m2=np.array([12.0]),
        area_after_m2=np.array([12.0]),
        node_ids_before=np.array([7]),
        node_ids_after=np.array([7]),
        plate_is_continental_before=np.array([False]),
        plate_is_continental_after=np.array([True]),
    )
    accounts = hm_ledger.cumulative_scopes(world)
    volume = 40_000.0 * 12.0
    assert accounts["continental_node"]["reclassification_m3"] == pytest.approx(volume)
    assert accounts["oceanic_node"]["reclassification_m3"] == pytest.approx(-volume)


def test_unaligned_topology_snapshot_does_not_pair_reused_ids_for_cap_transitions():
    world = _World()
    cap = lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M
    phase_budget.record(
        world,
        _Plate(),
        "plate_merge",
        np.array([35_000.0]),
        np.array([cap]),
        np.array([INHERIT]),
        np.array([35_000.0]),
        np.array([cap]),
        np.array([INHERIT]),
        area_before_m2=np.array([3.0]),
        area_after_m2=np.array([5.0]),
        node_ids_before=np.array([[9, 1]], dtype=np.uint64),
        node_ids_after=np.array([[9, 1]], dtype=np.uint64),
    )

    transitions = world.phase_budget["plate_merge"]["hm_cap_transitions"]
    assert transitions["aligned_calls"] == 0
    assert transitions["unaligned_calls"] == 1
    assert transitions["scopes"]["all"]["newly_capped_area_m2"] == pytest.approx(5.0)
    assert transitions["scopes"]["all"]["left_cap_area_m2"] == pytest.approx(3.0)
    assert transitions["scopes"]["all"]["persistently_capped_area_m2"] == 0.0


def test_suture_placement_books_cross_type_volume_as_reclassification():
    world = _World()
    hm_ledger.record_suture_front(
        world,
        donor_plate_id=1,
        neighbour_plate_ids=[2],
        donor_hm_m3=100.0,
        placed_hm_m3=80.0,
        donor_is_continental=False,
        placed_hm_continental_m3=25.0,
    )

    accounts = hm_ledger.cumulative_scopes(world)
    assert accounts["all"]["sink_m3"] == pytest.approx(20.0)
    assert accounts["oceanic_node"]["sink_m3"] == pytest.approx(20.0)
    assert accounts["continental_node"]["reclassification_m3"] == pytest.approx(25.0)
    assert accounts["oceanic_node"]["reclassification_m3"] == pytest.approx(-25.0)


def test_direct_retreat_reclassification_books_existing_hm_without_phase_double_count():
    world = _World()
    hm_ledger.record_reclassification(
        world,
        "subduction_and_suture_transfer",
        6.0 * 14.0,
        from_continental=False,
        to_continental=True,
    )

    accounts = hm_ledger.cumulative_scopes(world)
    assert accounts["all"]["source_m3"] == 0.0
    assert accounts["all"]["sink_m3"] == 0.0
    assert accounts["continental_node"]["reclassification_m3"] == pytest.approx(84.0)
    assert accounts["oceanic_node"]["reclassification_m3"] == pytest.approx(-84.0)
