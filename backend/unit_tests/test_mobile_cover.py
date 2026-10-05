"""Mobile cover (issue #297 phase 2) -- see erosion.MOBILE_COVER_*."""

import copy

import numpy as np
import pytest

from app import continental_ledger, erosion, lithosphere, mobile_cover, orogeny, persistence, quad_tectonics, volcanism
from app.elevation_lines import line_spacing_rad
from app.sparse_quad_patch import PlateWithSparseQuadPatch, cells_per_face_edge, pack_cell_keys
from app.world import World, generate_world

K = erosion.MOBILE_COVER_ERODIBILITY_FACTOR
COVER_FIELDS = ("mobile_cover_m", "mobile_cover_continental_m")


def _quad_world():
    world = generate_world(seed=3, num_plates=8, surface="quad")
    continental_ledger.ensure_initialized(world)
    return world


def _field(world, name) -> np.ndarray:
    return np.concatenate([p.collect(name) for p in world.plates])


def _volume(world, name) -> float:
    spacing = line_spacing_rad(world.node_density)
    return sum(float(np.dot(p.collect(name), p.accounting_areas_m2(spacing))) for p in world.plates)


def _assert_nested(world) -> None:
    hc, material = _field(world, "crustal_thickness_m"), _field(world, "continental_material_m")
    cover, cover_material = _field(world, "mobile_cover_m"), _field(world, "mobile_cover_continental_m")
    assert cover.min() >= 0.0 and cover_material.min() >= 0.0
    assert np.all(cover <= hc + 1e-6)
    assert np.all(cover_material <= np.minimum(cover, material) + 1e-6)


def test_cover_is_entrained_faster_than_substrate_and_exposes_the_substrate_law():
    entrained, detached = erosion._strip_mobile_cover(np.array([0.0, 100.0, 1.0]), np.full(3, 2.0))
    # Bare rock erodes by the substrate law, unchanged.
    assert entrained[0] == 0.0 and detached[0] == 2.0
    # Thick cover: the same forcing moves K times as much and shields the bed.
    assert entrained[1] == pytest.approx(K * 2.0) and detached[1] == 0.0
    # Thin cover: stripped, then the rest of the step cuts substrate at the bare rate.
    assert entrained[2] == 1.0 and detached[2] == pytest.approx(2.0 - 1.0 / K)


def test_cover_strip_agrees_across_step_splits():
    cover = np.array([0.0, 3.0, 7.5, 40.0])
    capacity = np.array([2.0, 2.0, 1.0, 0.5])
    entrained, detached = erosion._strip_mobile_cover(cover, capacity)
    split_entrained, split_detached, left = np.zeros(4), np.zeros(4), cover.copy()
    for _ in range(16):
        e, d = erosion._strip_mobile_cover(left, capacity / 16)
        split_entrained += e
        split_detached += d
        left -= e
    np.testing.assert_allclose(split_entrained, entrained)
    np.testing.assert_allclose(split_detached, detached)


def test_settled_sediment_becomes_mobile_cover_and_its_inventory_closes():
    world = _quad_world()
    for step in range(2):
        balance_before = continental_ledger.balance_error_m3(world)
        budget = erosion.apply_erosion(world, years=5_000_000).budget
        if step == 0:
            assert budget["mobile_cover_prior_m3"] == 0.0
            assert budget["mobile_cover_entrained_m3"] == 0.0
        else:
            # Last step's deposits are remobilized first.
            assert budget["mobile_cover_entrained_m3"] > 0.0
        # Every settled pathway, lake silt included, enters the cover.
        settled_silt = budget["lake_silt_m3"] - budget["hc_cap_overflow_lake_silt_terminal_m3"]
        assert budget["mobile_cover_deposited_m3"] == pytest.approx(budget["deposited_m3"] + settled_silt, rel=1e-9)
        assert budget["mobile_cover_deposited_m3"] > 0.0
        # prior + deposited - entrained - consolidated - clip == what is left on the surface.
        assert budget["mobile_cover_m3"] == pytest.approx(
            budget["mobile_cover_prior_m3"]
            + budget["mobile_cover_deposited_m3"]
            - budget["mobile_cover_entrained_m3"]
            - budget["mobile_cover_consolidated_m3"]
            - budget["mobile_cover_clip_m3"],
            rel=1e-9,
        )
        assert budget["mobile_cover_m3"] == pytest.approx(_volume(world, "mobile_cover_m"), rel=1e-9)
        assert budget["mobile_cover_continental_m3"] == pytest.approx(_volume(world, "mobile_cover_continental_m"), rel=1e-9)
        assert abs(budget["mobile_cover_clip_m3"]) <= 1e-9 * budget["mobile_cover_deposited_m3"]
        assert budget["bedrock_detached_m3"] == pytest.approx(budget["removed_m3"] - budget["mobile_cover_entrained_m3"])
        assert budget["mobile_cover_continental_entrained_m3"] <= budget["mobile_cover_entrained_m3"]
        # Cover and its provenance stay nested in the column, and the ledger still closes.
        _assert_nested(world)
        assert continental_ledger.balance_error_m3(world) == pytest.approx(balance_before, abs=1e-6 * _volume(world, "continental_material_m"))
    assert _volume(world, "mobile_cover_continental_m") > 0.0


def test_cover_is_remobilized_more_readily_than_bare_substrate():
    bare = _quad_world()
    covered = copy.deepcopy(bare)
    for plate in covered.plates:
        plate.set_fields_on_plate(mobile_cover_m=np.minimum(plate.collect("crustal_thickness_m"), 50.0))
    bare_budget = erosion.apply_erosion(bare, years=5_000_000).budget
    covered_budget = erosion.apply_erosion(covered, years=5_000_000).budget
    # The same forcing removes more where there is loose cover, and the cover shields the bed.
    assert covered_budget["removed_m3"] > 1.05 * bare_budget["removed_m3"]
    assert covered_budget["bedrock_detached_m3"] < bare_budget["bedrock_detached_m3"]
    # Removing the cover exposes the bare law: with no cover, all of it is substrate.
    assert bare_budget["bedrock_detached_m3"] == pytest.approx(bare_budget["removed_m3"])


def test_burial_consolidation_only_moves_the_cover_substrate_split(monkeypatch):
    world = _quad_world()
    for plate in world.plates:
        hc = plate.collect("crustal_thickness_m")
        cover = np.minimum(hc, 3000.0)
        plate.set_fields_on_plate(mobile_cover_m=cover, mobile_cover_continental_m=np.minimum(cover, plate.collect("continental_material_m")))
    never = copy.deepcopy(world)
    budget = erosion.apply_erosion(world, years=5_000_000).budget
    monkeypatch.setattr(erosion, "MOBILE_COVER_CONSOLIDATION_TIMESCALE_MYR", np.inf)
    never_budget = erosion.apply_erosion(never, years=5_000_000).budget
    assert never_budget["mobile_cover_consolidated_m3"] == 0.0
    assert budget["mobile_cover_consolidated_m3"] > 0.0
    assert 0.0 < budget["mobile_cover_continental_consolidated_m3"] <= budget["mobile_cover_consolidated_m3"]
    # Hc, the surface and the continental tracer are untouched; only the split moved.
    for name in ("crustal_thickness_m", "elevation", "continental_material_m"):
        np.testing.assert_array_equal(_field(world, name), _field(never, name))
    assert _volume(never, "mobile_cover_m") - _volume(world, "mobile_cover_m") == pytest.approx(budget["mobile_cover_consolidated_m3"], rel=1e-9)
    assert _volume(never, "mobile_cover_continental_m") - _volume(world, "mobile_cover_continental_m") == pytest.approx(
        budget["mobile_cover_continental_consolidated_m3"], rel=1e-9
    )
    # Never below the burial depth in one step.
    cover_never = _field(never, "mobile_cover_m")
    deep = cover_never > erosion.MOBILE_COVER_CONSOLIDATION_DEPTH_M
    assert np.all(_field(world, "mobile_cover_m")[deep] >= erosion.MOBILE_COVER_CONSOLIDATION_DEPTH_M)


def test_cover_above_a_thinned_column_is_clipped_and_reported_as_stale():
    world = _quad_world()
    plate = world.plates[0]
    hc = plate.collect("crustal_thickness_m")
    plate.set_fields_on_plate(mobile_cover_m=hc + 10.0)
    spacing = line_spacing_rad(world.node_density)
    expected = 10.0 * float(plate.accounting_areas_m2(spacing).sum())
    budget = erosion.apply_erosion(world, years=5_000_000).budget
    assert budget["stale_mobile_cover_excess_m3"] == pytest.approx(expected)
    _assert_nested(world)


def test_mobile_cover_survives_a_save_round_trip_and_quad_remaps():
    world = _quad_world()
    erosion.apply_erosion(world, years=5_000_000)
    loaded = persistence.load_world_bytes(persistence.save_world_bytes(world))
    for name in COVER_FIELDS:
        np.testing.assert_array_equal(_field(loaded, name), _field(world, name))

    plate = max(loaded.plates, key=lambda p: float(p.collect("mobile_cover_m").sum()))
    covered = np.flatnonzero(plate.collect("mobile_cover_m") > 0.0)[:4]
    before = {name: _volume(loaded, name) for name in COVER_FIELDS}
    lineage = plate.refine_cells(plate.cell_keys[covered])
    for name in COVER_FIELDS:
        assert _volume(loaded, name) == pytest.approx(before[name], rel=1e-12)
    _assert_nested(loaded)
    plate.coarsen_cells(np.array(sorted(lineage), dtype=np.int64))
    for name in COVER_FIELDS:
        assert _volume(loaded, name) == pytest.approx(before[name], rel=1e-12)
    _assert_nested(loaded)


# --- The cover through tectonics (mobile_cover.py) -----------------------------------------

DENSITY = 0.5
SPACING = line_spacing_rad(DENSITY)
N = cells_per_face_edge(SPACING)


def _block(i_range, j_range) -> np.ndarray:
    ii, jj = np.meshgrid(np.arange(*i_range), np.arange(*j_range), indexing="ij")
    return pack_cell_keys(np.zeros(ii.size, dtype=int), ii.ravel(), jj.ravel())


def _plate(plate_id, keys, crust_type="oceanic", cover=0.0, **fields) -> PlateWithSparseQuadPatch:
    hc, hm = lithosphere.reference_thickness(crust_type)
    count = len(keys)
    defaults = {"crustal_thickness_m": np.full(count, hc), "mantle_lithosphere_thickness_m": np.full(count, hm)}
    defaults.update(fields)
    plate = PlateWithSparseQuadPatch(plate_id, np.eye(3), crust_type, N, keys, fields=defaults)
    material = plate.collect("crustal_thickness_m") if crust_type == "continental" else np.zeros(count)
    plate.set_fields_on_plate(
        continental_material_m=material,
        mobile_cover_m=np.full(count, cover),
        mobile_cover_continental_m=np.minimum(np.full(count, cover / 2.0), material),
    )
    lithosphere.sync_plate_elevation(plate)
    return plate


def _tectonic_world(*plates) -> World:
    world = World(seed=0, plates=list(plates), next_plate_id=max(p.plate_id for p in plates) + 1, node_density=DENSITY)
    world.fault_deformation_mode = "smooth"
    continental_ledger.ensure_initialized(world)
    # Seeds initial_m3 with the cover the fixture laid down.
    mobile_cover.ensure_ledger(world)
    return world


def _cover_m3(plate, name="mobile_cover_m") -> float:
    return float(plate.collect(name) @ plate.node_areas_m2())


def _assert_plate_nested(plate) -> None:
    cover, cover_material = plate.collect("mobile_cover_m"), plate.collect("mobile_cover_continental_m")
    assert np.all(cover >= 0.0) and np.all(cover <= plate.collect("crustal_thickness_m") + 1e-6)
    assert np.all(cover_material >= 0.0)
    assert np.all(cover_material <= np.minimum(cover, plate.collect("continental_material_m")) + 1e-6)


def test_the_ledger_seeds_its_initial_inventory_and_only_signed_accounts_go_negative():
    plate = _plate(1, _block((10, 14), (20, 24)), "continental", cover=10.0)
    world = _tectonic_world(plate)
    assert world.mobile_cover_ledger["initial_m3"] == pytest.approx(_cover_m3(plate))
    assert mobile_cover.balance_error_m3(world) == pytest.approx(0.0, abs=1e-3)
    mobile_cover.record(world, "fault_advection_m3", -5.0)
    with pytest.raises(ValueError):
        mobile_cover.record(world, "subducted_m3", -5.0)
    with pytest.raises(KeyError):
        mobile_cover.record(world, "initial_m3", 5.0)


def test_orogenic_collapse_carries_the_cover_with_the_crust_it_moves():
    keys = pack_cell_keys(np.zeros(12, dtype=int), np.arange(10, 22), np.full(12, 20))
    hc = np.full(12, lithosphere.REFERENCE_HC_CONTINENTAL_M)
    hc[:2] = lithosphere.MAX_CRUSTAL_THICKNESS_M
    plate = _plate(1, keys, "continental", cover=100.0, crustal_thickness_m=hc)
    world = _tectonic_world(plate)
    before, before_material = _cover_m3(plate), _cover_m3(plate, "mobile_cover_continental_m")

    orogeny.relax_orogens(plate, world, 1_000_000.0)

    cover = plate.collect("mobile_cover_m")
    assert cover[:2].max() < 100.0 and cover[2] > 100.0  # thinned donors, thickened receiver
    assert _cover_m3(plate) == pytest.approx(before, rel=1e-12)
    assert _cover_m3(plate, "mobile_cover_continental_m") == pytest.approx(before_material, rel=1e-12)
    assert mobile_cover.balance_error_m3(world) == pytest.approx(0.0, abs=1e-6 * before)
    _assert_plate_nested(plate)


def test_rift_stretching_spreads_the_cover_onto_the_rifted_cells():
    plate = _plate(1, _block((10, 20), (20, 30)), "continental", cover=50.0)
    world = _tectonic_world(plate)
    before = _cover_m3(plate)
    old = set(map(int, plate.cell_keys))
    n = plate.node_count()

    assert quad_tectonics.grow_frontier(plate, world, np.ones(n, dtype=bool), np.zeros(n, dtype=bool), [], SPACING, 2, 10_000) > 0

    new = np.array([k not in old for k in map(int, plate.cell_keys)])
    assert np.all(plate.collect("mobile_cover_m")[new] > 0.0)
    ledger = world.mobile_cover_ledger
    assert _cover_m3(plate) + ledger["rift_reset_m3"] == pytest.approx(before, rel=1e-9)
    assert mobile_cover.balance_error_m3(world) == pytest.approx(0.0, abs=1e-9 * before)
    _assert_plate_nested(plate)


def test_trench_retreat_books_the_cover_it_subducts():
    a = _plate(1, _block((10, 24), (20, 30)), cover=10.0)
    b = _plate(2, _block((20, 34), (20, 30)))
    world = _tectonic_world(a, b)
    before = _cover_m3(a)

    a.deform(world, [b], 1_000_000, 10 * SPACING)

    assert world.mobile_cover_ledger["subducted_m3"] > 0.0
    assert mobile_cover.balance_error_m3(world) == pytest.approx(0.0, abs=1e-9 * before)
    _assert_plate_nested(a)


def test_suture_donors_metamorphose_their_cover_into_the_accreted_crust():
    keys = _block((10, 24), (20, 30))
    a = _plate(1, keys, "continental", cover=10.0, crustal_thickness_m=np.full(len(keys), 20_000.0))
    b = _plate(2, _block((20, 34), (20, 30)), "continental")
    world = _tectonic_world(a, b)
    before = _cover_m3(a)

    a.deform(world, [b], 1_000_000, 1.5 * SPACING)

    assert world.mobile_cover_ledger["accreted_m3"] > 0.0
    assert world.mobile_cover_ledger["subducted_m3"] == 0.0
    assert mobile_cover.balance_error_m3(world) == pytest.approx(0.0, abs=1e-9 * before)
    _assert_plate_nested(a)


def test_erupted_crust_buries_the_cover_under_it():
    plate = _plate(1, _block((10, 16), (20, 26)), "continental", cover=20.0)
    volcano = np.zeros(plate.node_count(), dtype=bool)
    volcano[[7, 21]] = True
    plate.set_fields_on_plate(is_volcano=volcano, volcano_active_years_remaining=np.where(volcano, 1.0e7, 0.0))
    world = _tectonic_world(plate)
    world.volcanism_multiplier = 1.0e3  # every active vent erupts
    hc_before = plate.collect("crustal_thickness_m")

    volcanism.apply_volcanic_activity(world, 1_000_000)

    gained = plate.collect("crustal_thickness_m") > hc_before
    assert gained[[7, 21]].all()
    cover = plate.collect("mobile_cover_m")
    assert np.all(cover[gained] == 0.0) and np.all(cover[~gained] == 20.0)
    assert world.mobile_cover_ledger["volcanic_buried_m3"] > 0.0
    assert mobile_cover.balance_error_m3(world) == pytest.approx(0.0, abs=1e-6)
