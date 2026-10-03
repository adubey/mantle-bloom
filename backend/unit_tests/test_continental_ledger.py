import numpy as np
import pytest

from app import continental_ledger, persistence
from app.elevation_lines import CRUST_TYPE_OCEANIC, ElevationLine
from app.plates import PlateWithLines
from app.sparse_quad_patch import PlateWithSparseQuadPatch, pack_cell_keys
from app.world import World


def _world() -> World:
    line = ElevationLine(
        phi=0.0,
        theta=np.array([0.0, 0.1]),
        elevation=np.array([100.0, 100.0]),
        crustal_thickness_m=np.array([30_000.0, 20_000.0]),
        mantle_lithosphere_thickness_m=np.array([100_000.0, 100_000.0]),
    )
    plate = PlateWithLines(plate_id=0, frame=np.eye(3), crust_type="continental", lines=[line])
    return World(seed=1, plates=[plate], node_density=1.0)


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


def test_initializes_area_weighted_provenance_and_survives_retyping():
    world = _world()
    continental_ledger.ensure_initialized(world)
    line = world.plates[0].lines[0]
    assert np.array_equal(line.continental_material_m, line.crustal_thickness_m)

    initial = continental_ledger.surface_volume_m3(world)
    line._crust_type_code[:] = CRUST_TYPE_OCEANIC
    inventory = continental_ledger.inventories(world)
    assert inventory["surface_continental_derived_m3"] == initial
    assert inventory["continental_sediment_on_oceanic_hosts_m3"] == initial


def test_record_rejects_invalid_or_unknown_entries():
    world = _world()
    continental_ledger.record(world, "juvenile_additions_m3", 12.5)
    assert world.continental_material_ledger["juvenile_additions_m3"] == 12.5
    with pytest.raises(ValueError):
        continental_ledger.record(world, "juvenile_additions_m3", -1.0)
    with pytest.raises(KeyError):
        continental_ledger.record(world, "not_an_account", 1.0)  # type: ignore[arg-type]


def test_ledger_and_node_provenance_round_trip_through_save_load():
    world = _world()
    continental_ledger.ensure_initialized(world)
    continental_ledger.record(world, "deeply_subducted_m3", 42.0)
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    assert loaded.continental_material_ledger == world.continental_material_ledger
    assert np.array_equal(
        loaded.plates[0].lines[0].continental_material_m,
        world.plates[0].lines[0].continental_material_m,
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
    assert continental_ledger.balance_error_m3(world) == pytest.approx(0.0, abs=1.0)


def test_balance_treats_relaminated_returns_as_a_surface_source():
    world = _world()
    continental_ledger.ensure_initialized(world)
    plate = world.plates[0]
    delta = np.array([0.0, 10.0])
    continental_ledger.add_material_thickness(
        world, plate, delta, "remelted_relaminated_returns_m3"
    )
    continental_ledger.assert_closed(world)
