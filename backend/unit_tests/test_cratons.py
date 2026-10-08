"""Persistent cratons (issue #274): seeding and formation, the resistances, every destruction
site's ledger booking, conservation through topology changes, and persistence -- see
cratons.py."""

from types import SimpleNamespace

import numpy as np
import pytest

from app import continental_ledger, cratons, lithosphere, persistence, quad_merge, quad_tectonics, rheology
from app import world as world_mod
from app.elevation_lines import CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC, line_spacing_rad
from app.lithosphere_plate import deform_columns
from app.surface_fields import CRATON_UNFORMED_YEARS
from app.sparse_quad_patch import PlateWithSparseQuadPatch, cells_per_face_edge, pack_cell_keys
from app.world import World

DENSITY = 0.5
SPACING = line_spacing_rad(DENSITY)
N = cells_per_face_edge(SPACING)


def _block(i_range, j_range, face: int = 0) -> np.ndarray:
    jj, ii = np.meshgrid(np.arange(*j_range), np.arange(*i_range), indexing="ij")
    return pack_cell_keys(np.full(ii.size, face), ii.ravel(), jj.ravel())


def _plate(plate_id, keys, crust_type="continental", frame=None, **fields) -> PlateWithSparseQuadPatch:
    hc, hm = lithosphere.reference_thickness(crust_type)
    defaults = {
        "crustal_thickness_m": np.full(len(keys), hc),
        "mantle_lithosphere_thickness_m": np.full(len(keys), hm),
    }
    defaults.update(fields)
    plate = PlateWithSparseQuadPatch(plate_id, np.eye(3) if frame is None else frame, crust_type, N, keys, fields=defaults)
    lithosphere.sync_plate_elevation(plate)
    return plate


def _world(*plates) -> World:
    world = World(seed=0, plates=list(plates), next_plate_id=max(p.plate_id for p in plates) + 1, node_density=DENSITY)
    world.fault_deformation_mode = "smooth"
    return world


def _cap(plate_id, frame=None, radius_cells=10, centre=np.array([1.0, 0.0, 0.0])) -> PlateWithSparseQuadPatch:
    radius = radius_cells * SPACING
    plate = PlateWithSparseQuadPatch.from_lattice(
        plate_id, np.eye(3) if frame is None else frame, "continental", SPACING, lambda pts: pts @ centre > np.cos(radius)
    )
    hc, hm = lithosphere.reference_thickness("continental")
    plate.set_fields_on_plate(
        crustal_thickness_m=np.full(plate.node_count(), hc), mantle_lithosphere_thickness_m=np.full(plate.node_count(), hm)
    )
    lithosphere.sync_plate_elevation(plate)
    return plate


def _volume(plate, name) -> float:
    return float(np.dot(plate.collect(name), plate.node_areas_m2()))


def _set_craton(plate, mask, formed=-1.0e9) -> None:
    hc = plate.collect("crustal_thickness_m")
    plate.set_fields_on_plate(
        craton_crust_m=np.where(mask, hc, 0.0), craton_formed_years=np.where(mask, formed, CRATON_UNFORMED_YEARS)
    )


def _unseeded(monkeypatch, *plates) -> World:
    monkeypatch.setattr(cratons, "CRATON_SEED_MARGIN_KM", 1e6)
    world = _world(*plates)
    cratons.ensure_initialized(world)
    return world


# --- Formation ---------------------------------------------------------------------------


def test_seeding_marks_only_deep_interior_cells_and_books_the_initial_volume():
    plate = _cap(1)
    world = _world(plate)

    cratons.ensure_initialized(world)

    craton = plate.collect("craton_crust_m")
    assert np.any(craton > 0.0)
    assert not np.all(craton > 0.0)
    probes = plate._probe_neighbour_indices()
    edge = np.any(np.any(probes < 0, axis=2), axis=1)
    assert np.all(craton[edge] == 0.0)
    seeded = craton > 0.0
    np.testing.assert_allclose(craton[seeded], plate.collect("crustal_thickness_m")[seeded])
    assert np.all(plate.collect("craton_formed_years")[seeded] == -cratons.CRATON_SEED_AGE_YEARS)
    assert world.craton_ledger["initial_m3"] == pytest.approx(_volume(plate, "craton_crust_m"))
    assert cratons.balance_error_m3(world) == pytest.approx(0.0, abs=1.0)


def test_a_seeded_world_is_never_reseeded():
    plate = _cap(1)
    world = _world(plate)
    cratons.ensure_initialized(world)
    plate.set_fields_on_plate(craton_crust_m=np.zeros(plate.node_count()))

    cratons.ensure_initialized(world)

    assert _volume(plate, "craton_crust_m") == 0.0


def test_thin_or_thick_crust_is_never_seeded():
    plate = _cap(1)
    world = _world(plate)
    plate.set_fields_on_plate(crustal_thickness_m=np.full(plate.node_count(), cratons.CRATON_HOST_MAX_HC_M + 5_000.0))

    cratons.ensure_initialized(world)

    assert world.craton_ledger["initial_m3"] == 0.0


def test_a_quiet_interior_cratonises_after_the_formation_time(monkeypatch):
    plate = _cap(1)
    world = _unseeded(monkeypatch, plate)
    world.elapsed_years = 400e6

    cratons.update(world, 0.5 * cratons.CRATON_FORMATION_MYR * 1e6)
    assert _volume(plate, "craton_crust_m") == 0.0
    cratons.update(world, 0.5 * cratons.CRATON_FORMATION_MYR * 1e6)

    craton = plate.collect("craton_crust_m")
    formed = craton > 0.0
    assert np.any(formed) and not np.all(formed)
    assert np.all(plate.collect("craton_formed_years")[formed] == pytest.approx(400e6 - cratons.CRATON_FORMATION_MYR * 1e6))
    assert world.craton_ledger["formed_m3"] == pytest.approx(_volume(plate, "craton_crust_m"))


def test_leaving_the_host_window_resets_the_formation_clock(monkeypatch):
    plate = _cap(1)
    world = _unseeded(monkeypatch, plate)
    cratons.update(world, 0.9 * cratons.CRATON_FORMATION_MYR * 1e6)
    plate.set_fields_on_plate(crustal_thickness_m=np.full(plate.node_count(), cratons.CRATON_HOST_MAX_HC_M + 1.0))
    cratons.update(world, 1e6)

    assert np.all(plate.collect("stable_continental_myr") == 0.0)


# --- Resistance and ledgered destruction --------------------------------------------------


def test_clip_books_only_what_drops_below_the_craton(monkeypatch):
    plate = _plate(1, _block((10, 14), (10, 11)))
    world = _unseeded(monkeypatch, plate)
    _set_craton(plate, np.array([True, True, False, False]))
    craton = plate.collect("craton_crust_m")
    hc = plate.collect("crustal_thickness_m")
    plate.set_fields_on_plate(crustal_thickness_m=hc - 1_000.0)

    lost = cratons.clip_to_column(world, plate, "eroded_m3")

    areas = plate.node_areas_m2()
    assert lost == pytest.approx(1_000.0 * float(areas[:2].sum()))
    assert world.craton_ledger["eroded_m3"] == pytest.approx(lost)
    np.testing.assert_allclose(plate.collect("craton_crust_m")[:2], craton[:2] - 1_000.0)


def test_retyped_columns_lose_their_craton_and_its_date(monkeypatch):
    plate = _plate(1, _block((10, 12), (10, 11)))
    world = _unseeded(monkeypatch, plate)
    _set_craton(plate, np.array([True, True]))
    plate.set_fields_on_plate(crust_type_code=np.array([CRUST_TYPE_OCEANIC, 0], dtype=np.int8))

    cratons.clip_to_column(world, plate, "rifted_m3")

    assert plate.collect("craton_crust_m")[0] == 0.0
    assert plate.collect("craton_formed_years")[0] == CRATON_UNFORMED_YEARS
    assert plate.collect("craton_formed_years")[1] != CRATON_UNFORMED_YEARS


def test_retreat_waits_for_a_sustained_overlap_scaled_by_strength(monkeypatch):
    plate = _plate(1, _block((10, 14), (10, 11)))
    world = _unseeded(monkeypatch, plate)
    _set_craton(plate, np.array([False, True, True, True]))
    world.elapsed_years = 100e6
    onset = np.array([0.0, 0.0, 99e6, 100e6 - cratons.CRATON_RETREAT_DELAY_YEARS])
    plate.set_fields_on_plate(overlap_onset_years=onset)

    np.testing.assert_array_equal(cratons.retreat_allowed(world, plate), [True, False, False, True])


def _divergent_ctx(n, closing_rate):
    zeros = np.zeros(n, dtype=bool)
    return SimpleNamespace(
        convergent=zeros, divergent=~zeros, transform=zeros, closing_rate=np.full(n, closing_rate),
        inputs=SimpleNamespace(neighbor_is_oceanic=zeros), arc_band=zeros, arc_intensity=np.zeros(n),
        fault_influence=np.ones(n), orogen_contested_strength=0.0, orogen_amount=1.0,
        fault_noise=None, own_points=np.zeros((n, 3)),
    )


def test_a_craton_resists_divergent_thinning_and_books_what_it_loses(monkeypatch):
    plate = _plate(1, _block((10, 12), (10, 11)))
    world = _unseeded(monkeypatch, plate)
    _set_craton(plate, np.array([False, True]))
    from app.lithosphere_plate import COLUMN_FIELDS

    fields = {name: plate.collect(name) for name in COLUMN_FIELDS}
    rate = -0.05 / rheology.SECONDS_PER_YEAR
    areas = plate.node_areas_m2()

    out = deform_columns(world, plate, _divergent_ctx(2, rate), fields, lambda: None, areas, 0, 1_000_000)

    start = fields["crustal_thickness_m"]
    thinned = start - out["crustal_thickness_m"]
    assert thinned[0] > 0.0
    assert thinned[1] == pytest.approx((1.0 - cratons.CRATON_RIFT_RESISTANCE) * thinned[0], rel=1e-9)
    np.testing.assert_allclose(out["craton_crust_m"][1], out["crustal_thickness_m"][1])
    assert world.craton_ledger["rifted_m3"] == pytest.approx(thinned[1] * areas[1])


def test_rift_stretch_draws_less_from_a_craton_and_books_its_thinning(monkeypatch):
    def opened(with_craton):
        plate = _plate(1, _block((10, 20), (20, 21)))
        world = _unseeded(monkeypatch, plate)
        if with_craton:
            _set_craton(plate, np.ones(plate.node_count(), dtype=bool))
        inserted = plate.insert_cells(pack_cell_keys(np.array([0]), np.array([20]), np.array([20])))
        assert inserted[0]
        new = plate._index_of_keys(pack_cell_keys(np.array([0]), np.array([20]), np.array([20])))
        before = plate.collect("crustal_thickness_m")
        quad_tectonics._open_rift(plate, world, 1, new, new, np.array([1.0]))
        return plate, world, before, new

    plain, _, plain_before, _ = opened(False)
    craton, world, craton_before, new = opened(True)
    old = np.ones(craton.node_count(), dtype=bool)
    old[new] = False
    plain_loss = float(np.dot(plain_before[old] - plain.collect("crustal_thickness_m")[old], plain.node_areas_m2()[old]))
    craton_loss = float(np.dot(craton_before[old] - craton.collect("crustal_thickness_m")[old], craton.node_areas_m2()[old]))
    # A donor asked for a quarter of the footprint thins by a/(a + D/4): a bit more than a
    # quarter as much, since thinning is not linear in the area asked of it.
    assert 0.0 < craton_loss < 0.5 * plain_loss
    assert world.craton_ledger["rifted_m3"] == pytest.approx(craton_loss, rel=1e-9)
    # What the craton withheld is fresh magmatic crust: the new cell is thinner, not missing.
    assert craton.collect("crustal_thickness_m")[new] < plain.collect("crustal_thickness_m")[new]


def test_subducted_cells_book_their_craton_and_continental_material(monkeypatch):
    plate = _plate(1, _block((10, 16), (10, 11)), "oceanic")
    world = _unseeded(monkeypatch, plate)
    n = plate.node_count()
    plate.set_fields_on_plate(continental_material_m=np.full(n, 500.0), craton_crust_m=np.full(n, 200.0))
    shrinkable = np.zeros(n, dtype=bool)
    shrinkable[-2:] = True
    ctx = SimpleNamespace(
        shrinkable=shrinkable, accrete=np.zeros(n, dtype=bool), spacing_rad=SPACING,
        oceanic_override_retreat_budget_hc=np.zeros(1), own_points=plate.all_points_and_elevation()[0],
        inputs=SimpleNamespace(neighbor_plate_id=np.full(n, -1)), neighbours=[],
    )
    areas = plate.node_areas_m2()

    # Freshly overlapped craton cells hold out...
    assert np.all(quad_tectonics._retreat(plate, world, ctx, 4 * SPACING, 10))
    assert world.craton_ledger["subducted_m3"] == 0.0
    # ...until the overlap has lasted long enough.
    world.elapsed_years = 1.0 + cratons.CRATON_RETREAT_DELAY_YEARS
    plate.set_fields_on_plate(overlap_onset_years=np.where(shrinkable, 1.0, 0.0))
    survivors = quad_tectonics._retreat(plate, world, ctx, 4 * SPACING, 10)

    removed = ~survivors
    assert np.any(removed)
    assert world.craton_ledger["subducted_m3"] == pytest.approx(200.0 * float(areas[removed].sum()))
    assert world.continental_material_ledger["deeply_subducted_m3"] == pytest.approx(500.0 * float(areas[removed].sum()))


def test_suture_accretion_reworks_the_donor_craton_and_delaminates_the_receivers_roots(monkeypatch):
    keys = _block((10, 30), (20, 21))
    near_cap = quad_tectonics.SUTURE_ACCRETION_MAX_HC_M - 1_500.0
    hc = np.full(len(keys), near_cap)
    hc[0] = 30_000.0
    plate = _plate(1, keys, crustal_thickness_m=hc)
    world = _unseeded(monkeypatch, plate)
    continental_ledger.ensure_initialized(world)
    donors = np.zeros(len(keys), dtype=bool)
    donors[0] = True
    # Both the donor and the receiving belt carry cratonic crust, so the test can tell whose
    # crust is reworked and whose sinks.
    craton = np.where(donors, hc, 5_000.0)
    plate.set_fields_on_plate(craton_crust_m=craton)
    areas = plate.node_areas_m2()
    material_before = _volume(plate, "continental_material_m")
    donor_craton = float(craton[0] * areas[0])
    belt_craton_before = float(craton[~donors] @ areas[~donors])

    # The near-cap belts are hot enough to carry eligible dense roots (orogeny.py), so once
    # they fill, part of their own lower crust may founder over this 1 Myr to make room.
    quad_tectonics._accrete_onto_survivors(plate, donors, ~donors, world, years=1_000_000.0)

    ledger = world.craton_ledger
    belt_craton_after = float(plate.collect("craton_crust_m")[~donors] @ areas[~donors])
    assert ledger["collision_reworked_m3"] == pytest.approx(donor_craton)
    assert ledger["delaminated_m3"] > 0.0
    assert ledger["delaminated_m3"] == pytest.approx(belt_craton_before - belt_craton_after)
    material_after = float(np.dot(plate.collect("continental_material_m")[~donors], areas[~donors]))
    delaminated = world.continental_material_ledger["delaminated_lower_crust_m3"]
    assert delaminated > 0.0
    assert material_after + delaminated == pytest.approx(material_before, rel=1e-9)


def test_a_relocated_terrane_keeps_its_craton_and_material_on_its_new_footprint_only(monkeypatch):
    keys = _block((10, 20), (20, 21))
    plate = _plate(1, keys, "oceanic")
    world = _unseeded(monkeypatch, plate)
    codes = plate.collect("crust_type_code")
    codes[0] = CRUST_TYPE_CONTINENTAL
    hc = plate.collect("crustal_thickness_m")
    hc[0] = 35_000.0
    material = np.full(len(keys), 100.0)  # a little continental sediment on the ocean floor
    material[0] = 35_000.0
    plate.set_fields_on_plate(crust_type_code=codes, crustal_thickness_m=hc, continental_material_m=material)
    donors = np.zeros(len(keys), dtype=bool)
    donors[0] = True
    _set_craton(plate, donors, formed=-2e9)
    areas = plate.node_areas_m2()
    craton_before = _volume(plate, "craton_crust_m")
    material_before = _volume(plate, "continental_material_m")

    quad_tectonics._accrete_onto_survivors(plate, donors, ~donors, world)

    survivors = ~donors
    terrane = (plate.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL) & survivors
    craton = plate.collect("craton_crust_m")
    assert np.any(terrane)
    assert float(np.dot(craton[terrane], areas[terrane])) == pytest.approx(craton_before, rel=1e-12)
    assert np.all(craton[survivors & ~terrane] == 0.0)
    assert np.all(plate.collect("craton_formed_years")[terrane] == -2e9)
    material_after = float(np.dot(plate.collect("continental_material_m")[survivors], areas[survivors]))
    assert material_after == pytest.approx(material_before, rel=1e-9)
    assert sum(world.craton_ledger[key] for key in cratons.CRATON_SINK_ACCOUNTS) == 0.0


# --- Topology, phases and persistence -----------------------------------------------------


def test_partition_and_merge_conserve_craton_volume_and_keep_the_oldest_date():
    keep = _cap(1)
    absorb = _cap(2, centre=np.array([np.cos(0.4), np.sin(0.4), 0.0]))
    for plate, formed in ((keep, -2.0e9), (absorb, -1.0e9)):
        _set_craton(plate, np.ones(plate.node_count(), dtype=bool), formed)
    before = _volume(keep, "craton_crust_m") + _volume(absorb, "craton_crust_m")

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    assert _volume(keep, "craton_crust_m") == pytest.approx(before, rel=1e-12)
    assert keep.collect("craton_formed_years").min() == -2.0e9
    halves = keep._plates_from_node_masks([keep.collect("elevation") >= np.median(keep.collect("elevation")), keep.collect("elevation") < np.median(keep.collect("elevation"))], [5, 6])
    assert sum(_volume(p, "craton_crust_m") for p in halves) == pytest.approx(before, rel=1e-12)


def test_a_craton_formed_at_year_zero_keeps_its_date_through_a_merge():
    keep = _cap(1)
    absorb = _cap(2, centre=np.array([np.cos(0.4), np.sin(0.4), 0.0]))
    _set_craton(keep, np.ones(keep.node_count(), dtype=bool), formed=0.0)
    _set_craton(absorb, np.ones(absorb.node_count(), dtype=bool), formed=5.0e8)

    quad_merge.merge(keep, absorb, np.zeros((0, 3)))

    formed = keep.collect("craton_formed_years")[keep.collect("craton_crust_m") > 0.0]
    assert formed.min() == 0.0
    assert np.all(formed <= 5.0e8)


def test_phase_audit_books_vanished_craton_to_the_residual_account(monkeypatch):
    a = _cap(1)
    b = _cap(2, centre=np.array([-1.0, 0.0, 0.0]))
    world = _unseeded(monkeypatch, a, b)
    _set_craton(b, np.ones(b.node_count(), dtype=bool))
    world.craton_ledger["initial_m3"] = cratons.live_volume_m3(world)
    vanished = _volume(b, "craton_crust_m")

    audit = cratons.PhaseAudit(world)
    world.plates = [a]
    audit.settle("unattributed_m3", "topology_removed_m3")

    assert world.craton_ledger["topology_removed_m3"] == pytest.approx(vanished)
    assert world.craton_ledger["unattributed_m3"] == 0.0
    assert cratons.balance_error_m3(world) == pytest.approx(0.0, abs=1.0)


def test_craton_state_round_trips_and_old_saves_are_seeded_on_their_first_step():
    plate = _cap(1)
    world = _world(plate)
    cratons.ensure_initialized(world)
    cratons.record(world, "rifted_m3", 7.0)

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))

    assert loaded.craton_ledger == world.craton_ledger
    np.testing.assert_array_equal(loaded.plates[0].collect("craton_crust_m"), plate.collect("craton_crust_m"))
    np.testing.assert_array_equal(loaded.plates[0].collect("craton_formed_years"), plate.collect("craton_formed_years"))

    del world.craton_ledger
    plate.set_fields_on_plate(craton_crust_m=np.zeros(plate.node_count()))
    old = persistence.load_world_bytes(persistence.save_world_bytes(world))
    # Loading leaves the plates untouched...
    assert "initial_m3" not in old.craton_ledger
    assert _volume(old.plates[0], "craton_crust_m") == 0.0
    # ...and the first step seeds them.
    old.simulate_climate_biomes = False
    world_mod.step_world(old, 1_000_000)
    assert old.craton_ledger["initial_m3"] > 0.0
    assert cratons.balance_error_m3(old) == pytest.approx(0.0, abs=1e-9 * old.craton_ledger["initial_m3"])


def test_record_rejects_negative_or_unknown_entries():
    world = _world(_cap(1))
    with pytest.raises(ValueError):
        cratons.record(world, "eroded_m3", -1.0)
    with pytest.raises(KeyError):
        cratons.record(world, "not_an_account", 1.0)
    cratons.record(world, "unattributed_m3", -1.0)


def test_stepping_a_world_keeps_the_craton_ledger_closed():
    world = world_mod.generate_world(5, node_density=0.25)
    assert world.craton_ledger["initial_m3"] > 0.0

    for _ in range(3):
        world_mod.step_world(world, 5_000_000)

    error = cratons.balance_error_m3(world)
    assert abs(error) <= 1e-9 * world.craton_ledger["initial_m3"]
    assert abs(world.craton_ledger["unattributed_m3"]) <= 1e-6 * world.craton_ledger["initial_m3"]
    stats = world.stats_history[-1]
    assert stats["craton_volume_km3"] == pytest.approx(cratons.live_volume_m3(world) / 1e9)


def test_a_craton_resists_glacial_flattening_like_every_other_erosion_source(monkeypatch):
    from app import erosion

    def run(send_m):
        world = world_mod.generate_world(21, num_plates=8, node_density=0.25)
        plate = next(p for p in world.plates if p.crust_type == "continental")
        hc = plate.collect("crustal_thickness_m")
        target = int(np.argmax(hc))
        craton = np.zeros(plate.node_count())
        craton[target] = hc[target]
        plate.set_fields_on_plate(craton_crust_m=craton)
        offset = sum(p.node_count() for p in world.plates[: world.plates.index(plate)])

        def flatten(hydro, ice_factor, years, removable_m, multiplier=1.0):
            send = np.zeros(hydro.neighbor_idx.shape)
            send[offset + target, 0] = min(send_m, removable_m[offset + target])
            return send

        monkeypatch.setattr(erosion, "_flatten", flatten)
        erosion.apply_erosion(world, years=1_000_000)
        return plate.collect("crustal_thickness_m")[target]

    sent = run(0.0) - run(40.0)
    assert sent == pytest.approx((1.0 - cratons.CRATON_EROSION_RESISTANCE) * 40.0, rel=1e-6)
