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
    assert cratons.balance_error_m3(world) == pytest.approx(0.0, abs=1e-9 * world.craton_ledger["initial_m3"])


def test_a_seeded_world_is_never_reseeded():
    plate = _cap(1)
    world = _world(plate)
    cratons.ensure_initialized(world)
    plate.set_fields_on_plate(craton_crust_m=np.zeros(plate.node_count()))

    cratons.ensure_initialized(world)

    assert _volume(plate, "craton_crust_m") == 0.0


def test_generated_cores_are_booked_and_the_interior_rule_is_skipped():
    """A cratons-first world's generation seeds its cores (RFC 0009, #462): first-step
    seeding books them as the initial volume and adds no deep interior of its own."""
    plate = _cap(1)
    core = plate.all_points_and_elevation()[0][:, 1] > 0.0
    cratons.seed_cores(plate, core)
    world = _world(plate)

    cratons.ensure_initialized(world)

    craton = plate.collect("craton_crust_m")
    assert np.array_equal(craton > 0.0, core)
    np.testing.assert_array_equal(craton[core], plate.collect("crustal_thickness_m")[core])
    assert np.all(plate.collect("craton_formed_years")[core] == -cratons.CRATON_SEED_AGE_YEARS)
    assert np.all(plate.collect("craton_formed_years")[~core] == CRATON_UNFORMED_YEARS)
    assert world.craton_ledger["initial_m3"] == pytest.approx(_volume(plate, "craton_crust_m"))
    assert cratons.balance_error_m3(world) == pytest.approx(0.0, abs=1e-9 * world.craton_ledger["initial_m3"])


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


def test_retreat_never_consumes_a_craton_however_long_it_is_overlapped(monkeypatch):
    plate = _plate(1, _block((10, 14), (10, 11)))
    world = _unseeded(monkeypatch, plate)
    _set_craton(plate, np.array([False, True, True, False]))
    world.elapsed_years = 1e9
    plate.set_fields_on_plate(overlap_onset_years=np.full(4, 1.0))

    np.testing.assert_array_equal(cratons.retreat_allowed(plate), [True, False, False, True])


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
    assert world.craton_ledger["rifted_divergent_m3"] == world.craton_ledger["rifted_m3"]


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
    assert world.craton_ledger["rifted_stretch_m3"] == world.craton_ledger["rifted_m3"]
    # What the craton withheld is fresh magmatic crust: the new cell is thinner, not missing.
    assert craton.collect("crustal_thickness_m")[new] < plain.collect("crustal_thickness_m")[new]


def _retreat_ctx(plate, shrinkable, accrete):
    n = plate.node_count()
    return SimpleNamespace(
        shrinkable=shrinkable, accrete=accrete, terrane=np.zeros(n, dtype=bool), spacing_rad=SPACING,
        oceanic_override_retreat_budget_hc=np.full(1, np.inf), own_points=plate.all_points_and_elevation()[0],
        inputs=SimpleNamespace(neighbor_plate_id=np.full(n, -1), direction_to_neighbor=np.zeros((n, 3))), neighbours=[],
        suture_upper_plate_id=np.full(n, -1), suture_hm_subduct=np.zeros(n, dtype=bool),
    )


def test_a_craton_at_a_trench_wears_back_by_attrition_before_its_cell_retreats(monkeypatch):
    plate = _plate(1, _block((10, 16), (10, 11)))
    world = _unseeded(monkeypatch, plate)
    continental_ledger.ensure_initialized(world)
    n = plate.node_count()
    plate.set_fields_on_plate(continental_material_m=plate.collect("crustal_thickness_m"), craton_crust_m=np.full(n, 2_000.0))
    shrinkable = np.zeros(n, dtype=bool)
    shrinkable[-2:] = True  # overridden by an oceanic plate: a trench, not a suture
    ctx = _retreat_ctx(plate, shrinkable, np.zeros(n, dtype=bool))
    areas = plate.node_areas_m2()
    hc_before = plate.collect("crustal_thickness_m")

    # Over 1 Myr the trench's exposed cells (in a one-cell strip, both) wear back by the
    # attrition rate along their edge, and nothing retreats while it still has craton.
    assert np.all(quad_tectonics._retreat(plate, world, ctx, 4 * SPACING, 10, years=1e6))
    worn = cratons.CRATON_ATTRITION_KM2_PER_MYR * 1e6 / np.sqrt(areas)
    craton = plate.collect("craton_crust_m")
    np.testing.assert_allclose(craton[-2:], 2_000.0 - worn[-2:])
    assert np.all(craton[:-2] == 2_000.0)
    np.testing.assert_allclose(plate.collect("crustal_thickness_m")[-2:], hc_before[-2:] - worn[-2:])
    assert world.craton_ledger["subducted_m3"] == pytest.approx(worn[-2:] @ areas[-2:])
    assert world.continental_material_ledger["deeply_subducted_m3"] == pytest.approx(worn[-2:] @ areas[-2:])

    # Once their craton is gone, the cells retreat down the trench like any others.
    survivors = quad_tectonics._retreat(plate, world, ctx, 4 * SPACING, 10, years=1e9)
    assert not np.any(survivors[-2:])
    assert world.craton_ledger["subducted_m3"] == pytest.approx(2_000.0 * areas[-2:].sum())


def test_a_craton_at_a_collision_front_jams_and_loses_nothing(monkeypatch):
    plate = _plate(1, _block((10, 16), (10, 11)))
    world = _unseeded(monkeypatch, plate)
    continental_ledger.ensure_initialized(world)
    n = plate.node_count()
    _set_craton(plate, np.arange(n) >= n - 2)
    shrinkable = np.zeros(n, dtype=bool)
    shrinkable[-3:] = True
    ctx = _retreat_ctx(plate, shrinkable, shrinkable.copy())
    craton_before = _volume(plate, "craton_crust_m")

    survivors = quad_tectonics._retreat(plate, world, ctx, 4 * SPACING, 10, years=1e9)

    # The weak crust beside the craton takes the shortening; the craton holds, whole.
    np.testing.assert_array_equal(survivors[-3:], [False, True, True])
    assert _volume(plate, "craton_crust_m") == pytest.approx(craton_before, rel=1e-12)
    assert sum(world.craton_ledger[key] for key in cratons.CRATON_SINK_ACCOUNTS) == 0.0


def test_suture_accretion_reworks_the_donor_craton_and_the_receivers_shed_roots_leave_theirs(monkeypatch):
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
    # crust is reworked and whose stays.
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
    # The roots are the belts' non-cratonic lower crust: their craton stays (#340).
    assert ledger["delaminated_m3"] == 0.0
    assert belt_craton_after == belt_craton_before
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

    quad_tectonics._accrete_onto_survivors(plate, donors, ~donors, world, terranes=donors)

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
    assert cratons.balance_error_m3(world) == pytest.approx(0.0, abs=1e-9 * world.craton_ledger["initial_m3"])


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


def test_a_rift_path_account_also_books_rifted_m3_outside_the_balance():
    # RFC 0003 phase 1: the paths are diagnostics; the balance counts rifted_m3 once.
    world = _world(_cap(1))
    cratons.ensure_initialized(world)
    cratons.record(world, "rifted_m3", 2.0)  # booked before a save tracked the paths
    cratons.record(world, "rifted_failed_m3", 3.0)
    cratons.record(world, "rifted_stretch_m3", 5.0)

    ledger = world.craton_ledger
    assert ledger["rifted_m3"] == 10.0
    assert (ledger["rifted_divergent_m3"], ledger["rifted_stretch_m3"], ledger["rifted_failed_m3"]) == (0.0, 5.0, 3.0)
    assert cratons.booked_sinks_m3(world) == 10.0


def test_an_old_ledger_gains_the_rift_diagnostics_on_load():
    world = _world(_cap(1))
    cratons.ensure_initialized(world)
    for key in cratons.CRATON_RIFT_PATH_ACCOUNTS + cratons.CRATON_RIFT_CUT_COUNTERS:
        del world.craton_ledger[key]

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))

    assert all(loaded.craton_ledger[key] == 0.0 for key in cratons.CRATON_RIFT_PATH_ACCOUNTS + cratons.CRATON_RIFT_CUT_COUNTERS)
    cratons.record(loaded, "rifted_divergent_m3", 1.0)
    assert loaded.craton_ledger["rifted_m3"] == 1.0


def test_a_save_from_before_merge_capped_loads_with_it_at_zero_and_balanced():
    # #441: older saves hold their merge losses in delaminated_m3 and lack the new account.
    world = _world(_cap(1))
    cratons.ensure_initialized(world)
    del world.craton_ledger["merge_capped_m3"]
    world.craton_ledger["delaminated_m3"] = 0.0

    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))

    assert loaded.craton_ledger["merge_capped_m3"] == 0.0
    assert abs(cratons.balance_error_m3(loaded)) <= 1e-9 * max(loaded.craton_ledger["initial_m3"], 1.0)


@pytest.mark.parametrize("resistance", [0.0, 0.75, 1.0, 1.5])
def test_rift_thinning_factors_stay_in_unit_range(monkeypatch, resistance):
    # RFC 0003: resistance is clamped, so no strength or tuned constant can make a failed
    # rift thicken a column or thin it past zero.
    from app.merge_split import FAILED_RIFT_THINNING_FRACTION

    monkeypatch.setattr(cratons, "CRATON_RIFT_RESISTANCE", resistance)
    craton = np.linspace(0.0, 2.0 * cratons.CRATON_FULL_STRENGTH_HC_M, 41)
    resisted = cratons.rift_resistance(craton)
    assert np.all((resisted >= 0.0) & (resisted <= 1.0))
    for taper in np.linspace(0.0, 1.0, 11):
        factor = np.clip(1.0 - FAILED_RIFT_THINNING_FRACTION * taper * (1.0 - resisted), 0.0, 1.0)
        assert np.all((factor >= 1.0 - FAILED_RIFT_THINNING_FRACTION) & (factor <= 1.0))


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


def test_an_audit_clip_books_to_its_own_account_not_a_mechanism(monkeypatch):
    # #340: the deform-phase and end-of-step audits clip to `audit_clipped_m3`; the stepped
    # world is covered in stress_tests/test_cratons.py.
    plate = _plate(1, _block((10, 13), (10, 11)))
    world = _unseeded(monkeypatch, plate)
    _set_craton(plate, np.ones(3, dtype=bool))
    world.craton_ledger["initial_m3"] = cratons.live_volume_m3(world)
    audit = cratons.PhaseAudit(world)
    hc = plate.collect("crustal_thickness_m")
    clipped = 0.5 * hc[0] * plate.node_areas_m2()[0]
    hc[0] *= 0.5  # a writer thins a cratonic column without booking it
    plate.set_fields_on_plate(crustal_thickness_m=hc)

    audit.settle("audit_clipped_m3")

    assert "audit_clipped_m3" in cratons.CRATON_SINK_ACCOUNTS
    assert world.craton_ledger["audit_clipped_m3"] == pytest.approx(clipped)
    assert cratons.balance_error_m3(world) == pytest.approx(0.0, abs=1e-9 * world.craton_ledger["initial_m3"])


def test_a_craton_resists_glacial_flattening_like_every_other_erosion_source(monkeypatch):
    from app.hydroclimate import erosion

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


def test_trench_attrition_on_a_continental_margin_is_capped_by_the_arc_budget(monkeypatch):
    # #177: a continental margin under an oceanic neighbour subducts only as much as the arc
    # magmatism allows, and craton attrition is charged to the same budget (#340).
    def worn(budget):
        plate = _plate(1, _block((10, 16), (10, 11)))
        world = _unseeded(monkeypatch, plate)
        continental_ledger.ensure_initialized(world)
        n = plate.node_count()
        plate.set_fields_on_plate(continental_material_m=plate.collect("crustal_thickness_m"), craton_crust_m=np.full(n, 2_000.0))
        shrinkable = np.arange(n) >= n - 2
        ctx = _retreat_ctx(plate, shrinkable, np.zeros(n, dtype=bool))
        ctx.oceanic_override_retreat_budget_hc = np.array([budget])
        quad_tectonics._retreat(plate, world, ctx, 4 * SPACING, 10, years=1e6)
        return world.craton_ledger["subducted_m3"], ctx.oceanic_override_retreat_budget_hc[0]

    assert worn(0.0) == (0.0, 0.0)
    full, _ = worn(np.inf)
    node_area = lithosphere.node_area_m2(SPACING)
    half_budget = 0.5 * full / node_area
    lost, left = worn(half_budget)
    assert lost == pytest.approx(0.5 * full)
    assert left == pytest.approx(0.0, abs=1e-9 * half_budget)


def test_coarsening_keeps_a_column_with_craton_continental():
    # #340: the area vote would make this parent oceanic and strip the craton it averages in.
    plate = _plate(1, _block((10, 12), (10, 11)))
    parent = plate.cell_keys[0]
    children = plate.refine_cells(np.array([parent]))[int(parent)]
    child = np.isin(plate.cell_keys, children)
    codes = plate.collect("crust_type_code")
    codes[np.flatnonzero(child)[:3]] = CRUST_TYPE_OCEANIC
    craton = np.where(child & (codes != CRUST_TYPE_OCEANIC), plate.collect("crustal_thickness_m"), 0.0)
    plate.set_fields_on_plate(crust_type_code=codes, craton_crust_m=craton)
    volume = _volume(plate, "craton_crust_m")

    plate.coarsen_cells(np.array([parent]))

    at = int(np.flatnonzero(plate.cell_keys == parent)[0])
    assert plate.collect("crust_type_code")[at] != CRUST_TYPE_OCEANIC
    assert _volume(plate, "craton_crust_m") == pytest.approx(volume, rel=1e-12)
