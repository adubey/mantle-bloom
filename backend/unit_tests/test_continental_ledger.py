import numpy as np
import pytest

from app import continental_ledger, persistence
from app.elevation_lines import CRUST_TYPE_OCEANIC, ElevationLine, line_spacing_rad
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


def test_initializes_legacy_world_when_ledger_attribute_is_absent():
    world = _world()
    del world.continental_material_ledger

    continental_ledger.ensure_initialized(world)

    assert np.array_equal(
        world.plates[0].collect("continental_material_m"),
        world.plates[0].collect("crustal_thickness_m"),
    )
    assert world.continental_material_ledger["initial_continental_m3"] > 0.0


def test_assert_closed_rejects_material_above_host_crust():
    world = _world()
    continental_ledger.ensure_initialized(world)
    plate = world.plates[0]
    plate.set_fields_on_plate(
        continental_material_m=plate.collect("crustal_thickness_m") + 1.0
    )

    with pytest.raises(AssertionError, match="within \\[0, Hc\\]"):
        continental_ledger.assert_closed(world)


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
    plate.set_fields_on_plate(
        crustal_thickness_m=plate.collect("crustal_thickness_m") + delta
    )
    continental_ledger.add_material_thickness(
        world, plate, delta, "remelted_relaminated_returns_m3"
    )
    continental_ledger.assert_closed(world)


def test_suture_accretion_returns_material_over_the_cap():
    from app.lithosphere import crust_density, isostatic_elevation
    from app.lithosphere_plate import SUTURE_ACCRETION_MAX_HC_M, _redistribute_accreted_column

    rho_c = crust_density("continental")
    hc = np.full(6, 60_000.0)
    fields = {
        "crustal_thickness_m": hc,
        "mantle_lithosphere_thickness_m": np.full(6, 100_000.0),
        "continental_material_m": hc.copy(),
    }
    elevation = isostatic_elevation(hc, fields["mantle_lithosphere_thickness_m"], rho_c).copy()
    material_before = fields["continental_material_m"].sum()
    removed = np.full(5, 30_000.0)

    unplaced = _redistribute_accreted_column(
        fields, elevation, rho_c, removed, removed, removed, np.ones(5, dtype=bool), from_high=True,
    )

    assert np.all(fields["continental_material_m"] <= fields["crustal_thickness_m"] + 1e-9)
    assert np.all(fields["crustal_thickness_m"] <= SUTURE_ACCRETION_MAX_HC_M + 1e-6)
    assert unplaced > 0.0
    assert fields["continental_material_m"].sum() + unplaced == pytest.approx(material_before + removed.sum())


def test_regularize_line_never_lifts_material_above_host_crust():
    """Halving a line's node count leaves less Hc than material; the remainder must stay out
    of the one column with headroom rather than pile onto it (the caller books it)."""
    from app.elevation_lines import regularize_line

    spacing = line_spacing_rad(1.0)
    material = np.full(9, 40_000.0)
    material[4] = 0.0
    line = ElevationLine(
        phi=0.0,
        theta=np.linspace(0.0, 4 * spacing, 9),
        elevation=np.zeros(9),
        crustal_thickness_m=np.full(9, 40_000.0),
        mantle_lithosphere_thickness_m=np.full(9, 100_000.0),
        continental_material_m=material,
    )

    regularized = regularize_line(line, spacing)

    assert len(regularized) < len(line)
    assert np.all(regularized.continental_material_m <= regularized.crustal_thickness_m + 1e-6)
    assert regularized.continental_material_m.sum() < line.continental_material_m.sum()


def test_leading_row_drop_books_material_it_cannot_place():
    from app.lithosphere import MAX_CRUSTAL_THICKNESS_M, node_area_m2
    from app.lithosphere_plate import LithospherePlate

    spacing = line_spacing_rad(1.0)
    hc = MAX_CRUSTAL_THICKNESS_M - 1_000.0
    lines = [
        ElevationLine(
            phi=k * spacing,
            theta=np.linspace(-0.1, 0.1, 20),
            elevation=np.zeros(20),
            crustal_thickness_m=np.full(20, hc),
            mantle_lithosphere_thickness_m=np.full(20, 100_000.0),
            continental_material_m=np.full(20, hc),
        )
        for k in range(4)
    ]
    plate = LithospherePlate(plate_id=0, frame=np.eye(3), crust_type="continental", lines=lines)
    world = World(seed=0, plates=[plate], mantle_centers=[], node_density=1.0)
    continental_ledger.ensure_initialized(world)
    material_before = sum(float(ln.continental_material_m.sum()) for ln in lines)

    kept = plate._accrete_dropped_row_volume(world, lines, {round(3 * spacing, 6): "hi"})

    material_after = sum(float(ln.continental_material_m.sum()) for ln in kept)
    ledger = world.continental_material_ledger
    assert ledger["delaminated_lower_crust_m3"] > 0.0
    assert (material_before - material_after) * node_area_m2(spacing) == pytest.approx(
        ledger["delaminated_lower_crust_m3"] + ledger["deeply_subducted_m3"]
    )
