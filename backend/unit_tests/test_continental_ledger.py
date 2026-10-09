import numpy as np
import pytest

from app import continental_ledger, persistence
from app.elevation_lines import CRUST_TYPE_OCEANIC
from app.sparse_quad_patch import PlateWithSparseQuadPatch, pack_cell_keys
from app.world import World


def _quad_world() -> World:
    keys = pack_cell_keys(
        np.array([0, 0]), np.array([0, 1]), np.array([0, 0])
    )
    plate = PlateWithSparseQuadPatch(
        0,
        np.eye(3),
        "continental",
        4,
        keys,
        fields={
            "crustal_thickness_m": np.array([30_000.0, 20_000.0]),
            "mantle_lithosphere_thickness_m": np.array([100_000.0, 100_000.0]),
        },
    )
    return World(seed=1, plates=[plate], node_density=1.0)


def test_initializes_legacy_world_when_ledger_attribute_is_absent():
    world = _quad_world()
    del world.continental_material_ledger

    continental_ledger.ensure_initialized(world)

    assert np.array_equal(
        world.plates[0].collect("continental_material_m"),
        world.plates[0].collect("crustal_thickness_m"),
    )
    assert world.continental_material_ledger["initial_continental_m3"] > 0.0


def test_assert_closed_rejects_material_above_host_crust():
    world = _quad_world()
    continental_ledger.ensure_initialized(world)
    plate = world.plates[0]
    plate.set_fields_on_plate(
        continental_material_m=plate.collect("crustal_thickness_m") + 1.0
    )

    with pytest.raises(AssertionError, match="within \\[0, Hc\\]"):
        continental_ledger.assert_closed(world)


def test_record_rejects_invalid_or_unknown_entries():
    world = _quad_world()
    continental_ledger.record(world, "juvenile_additions_m3", 12.5)
    assert world.continental_material_ledger["juvenile_additions_m3"] == 12.5
    with pytest.raises(ValueError):
        continental_ledger.record(world, "juvenile_additions_m3", -1.0)
    with pytest.raises(KeyError):
        continental_ledger.record(world, "not_an_account", 1.0)  # type: ignore[arg-type]


def test_ledger_and_node_provenance_round_trip_through_save_load():
    world = _quad_world()
    continental_ledger.ensure_initialized(world)
    continental_ledger.record(world, "deeply_subducted_m3", 42.0)
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.continental_material_ledger == world.continental_material_ledger
    assert np.array_equal(
        loaded.plates[0].collect("continental_material_m"),
        world.plates[0].collect("continental_material_m"),
    )


def test_quad_inventory_uses_exact_cell_areas_and_survives_retyping():
    world = _quad_world()
    continental_ledger.ensure_initialized(world)
    plate = world.plates[0]
    material = plate.collect("continental_material_m")
    expected = float(np.sum(material * plate.node_areas_m2()))
    assert continental_ledger.surface_volume_m3(world) == pytest.approx(expected)

    plate.set_fields_on_plate(
        crust_type_code=np.full(plate.node_count(), CRUST_TYPE_OCEANIC, dtype=np.int8)
    )
    inventory = continental_ledger.inventories(world)
    assert inventory["surface_continental_derived_m3"] == pytest.approx(expected)
    assert inventory["continental_sediment_on_oceanic_hosts_m3"] == pytest.approx(expected)


def test_quad_refinement_conserves_continental_material_volume():
    world = _quad_world()
    continental_ledger.ensure_initialized(world)
    before = continental_ledger.surface_volume_m3(world)
    plate = world.plates[0]
    plate.refine_cells(plate.cell_keys[:1])
    assert continental_ledger.surface_volume_m3(world) == pytest.approx(before)


def test_add_material_thickness_updates_tracer_account_and_balance_together():
    world = _quad_world()
    continental_ledger.ensure_initialized(world)
    plate = world.plates[0]
    delta = np.array([125.0, 0.0])
    expected = float(np.dot(delta, plate.node_areas_m2()))

    booked = continental_ledger.add_material_thickness(
        world, plate, delta, "juvenile_additions_m3"
    )

    assert booked == pytest.approx(expected)
    assert world.continental_material_ledger["juvenile_additions_m3"] == pytest.approx(expected)
    # Relative to the inventory: at ~2.5e17 m^3 one float64 step is 32 m^3, so a fixed 1 m^3
    # tolerance passes or fails with the platform's rounding.
    inventory = world.continental_material_ledger["initial_continental_m3"]
    assert continental_ledger.balance_error_m3(world) == pytest.approx(0.0, abs=1e-10 * inventory)


def test_balance_treats_relaminated_returns_as_a_surface_source():
    world = _quad_world()
    continental_ledger.ensure_initialized(world)
    plate = world.plates[0]
    delta = np.array([0.0, 10.0])
    plate.set_fields_on_plate(
        crustal_thickness_m=plate.collect("crustal_thickness_m") + delta
    )
    continental_ledger.add_material_thickness(
        world, plate, delta, "remelted_relaminated_returns_m3"
    )
    continental_ledger.assert_closed(world)
