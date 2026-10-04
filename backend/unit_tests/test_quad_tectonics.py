"""Per-step deformation on sparse quad plates (issue #228 Phase 4): the cell-graph topology
primitives, 2D boundary retreat/advance, and the step passes that had to learn about quad
plates (volcanism, gap filling, overlap tracking, merge eligibility, relattice) -- see
quad_tectonics.py."""

import numpy as np
import pytest

from app import continental_ledger, gaps, geometry, lithosphere, merge_split, orogeny, plates, quad_tectonics, volcanism
from app.elevation_lines import CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC, line_spacing_rad
from app.lithosphere_plate import (
    CONTINENTAL_CONTESTED_RETREAT_MIN_RUN,
    EXTEND_THRESHOLD_MULTIPLIER,
    SUTURE_ACCRETION_MAX_HC_M,
    SUTURE_ACCRETION_SPREAD_NODES,
    boundary_context,
    new_plate,
)
from app.sparse_quad_patch import PlateWithSparseQuadPatch, cells_per_face_edge, pack_cell_keys, unpack_cell_keys
from app.world import World

DENSITY = 0.5
SPACING = line_spacing_rad(DENSITY)
N = cells_per_face_edge(SPACING)


def _block(i_range, j_range, face: int = 0) -> np.ndarray:
    jj, ii = np.meshgrid(np.arange(*j_range), np.arange(*i_range), indexing="ij")
    return pack_cell_keys(np.full(ii.size, face), ii.ravel(), jj.ravel())


def _plate(plate_id, keys, crust_type="oceanic", frame=None, **fields) -> PlateWithSparseQuadPatch:
    count = len(keys)
    hc, hm = lithosphere.reference_thickness(crust_type)
    defaults = {
        "crustal_thickness_m": np.full(count, hc),
        "mantle_lithosphere_thickness_m": np.full(count, hm),
    }
    defaults.update(fields)
    plate = PlateWithSparseQuadPatch(plate_id, np.eye(3) if frame is None else frame, crust_type, N, keys, fields=defaults)
    if "elevation" not in fields:
        lithosphere.sync_plate_elevation(plate)
    return plate


def _world(*plates) -> World:
    world = World(seed=0, plates=list(plates), next_plate_id=max(p.plate_id for p in plates) + 1, node_density=DENSITY)
    world.fault_deformation_mode = "smooth"
    return world


def _columns(plate) -> tuple[np.ndarray, np.ndarray]:
    _, _, i, j = unpack_cell_keys(plate.cell_keys)
    return i, j


# --- Topology primitives ------------------------------------------------------------------


def test_insert_cells_keeps_existing_ids_and_skips_active_or_duplicate_keys():
    plate = _plate(1, _block((10, 12), (10, 12)), elevation=np.arange(4.0))
    before = dict(zip(map(int, plate.cell_keys), plate.collect("elevation")))
    fresh = int(pack_cell_keys(0, 12, 10))
    keys = np.array([fresh, fresh, int(plate.cell_keys[0])])

    inserted = plate.insert_cells(keys, {"elevation": np.array([50.0, 60.0, 70.0])})

    np.testing.assert_array_equal(inserted, [True, False, False])
    after = dict(zip(map(int, plate.cell_keys), plate.collect("elevation")))
    assert after.pop(fresh) == 50.0
    assert after == before
    assert plate.topology_revision == 1


def test_insert_cells_rejects_cells_that_would_break_two_to_one_balance():
    plate = _plate(1, _block((10, 12), (10, 12)))
    target = int(pack_cell_keys(0, 11, 11))
    grandchild = plate.refine_cells(np.array([target]))[target][1]
    plate.refine_cells(np.array([grandchild]))
    # A level-0 cell beside the level-2 leaves would touch a leaf two levels finer.
    beside = int(pack_cell_keys(0, 12, 11))

    inserted = plate.insert_cells(np.array([beside]))

    assert not inserted[0]
    assert beside not in set(map(int, plate.cell_keys))
    plate._validate_leaf_topology()


def test_insert_cells_rejects_a_cell_overlapping_a_finer_active_leaf():
    plate = _plate(1, _block((10, 11), (10, 11)))
    parent = int(plate.cell_keys[0])
    plate.refine_cells(np.array([parent]))

    assert not plate.insert_cells(np.array([parent]))[0]


def test_empty_neighbour_keys_cross_a_cube_face_seam():
    corner = pack_cell_keys(np.array([0]), np.array([N - 1]), np.array([N // 2]))
    plate = _plate(1, corner)

    sources, candidates = plate.empty_neighbour_keys(np.array([0]))

    assert len(candidates) == 4
    assert np.all(sources == 0)
    faces = unpack_cell_keys(candidates)[0]
    assert set(faces.tolist()) == {0, 1}
    assert np.all(plate.insert_cells(candidates))
    graph = plate.adjacency()
    centre = int(plate._index_of_keys(corner)[0])
    assert graph.offsets[centre + 1] - graph.offsets[centre] == 4


def test_remove_cells_keeps_survivor_values():
    plate = _plate(1, _block((10, 13), (10, 11)), elevation=np.array([1.0, 2.0, 3.0]))

    plate.remove_cells(np.array([False, True, False]))

    np.testing.assert_array_equal(plate.collect("elevation"), [1.0, 3.0])


def test_hop_distance_and_components_walk_the_cell_graph():
    plate = _plate(1, _block((10, 16), (10, 11)))
    mask = np.zeros(6, dtype=bool)
    mask[0] = True

    np.testing.assert_array_equal(quad_tectonics.hop_distance(plate, mask, 3), [0, 1, 2, 3, 4, 4])

    split = np.array([True, True, False, True, True, True])
    np.testing.assert_array_equal(
        quad_tectonics.components_of_at_least(plate, split, 3), [False, False, False, True, True, True]
    )


# --- Advance ------------------------------------------------------------------------------


def test_boundary_advances_into_open_ground_but_stops_short_of_a_neighbour():
    a = _plate(1, _block((10, 20), (20, 30)))
    b = _plate(2, _block((26, 36), (20, 30)))
    world = _world(a, b)
    count_before = a.node_count()

    a.deform(world, [b], 1_000_000, 0.0)

    assert a.node_count() > count_before
    a._validate_leaf_topology()
    new = a.collect("node_created_years") == world.elapsed_years
    assert np.any(new)
    # Every new cell is fresh crust: no neighbour's territory, and clear of its nodes.
    new_points = a.all_points_and_elevation()[0][new]
    assert not np.any(b.contains_batch(new_points))
    distances = np.arccos(np.clip(new_points @ b.all_points_and_elevation()[0].T, -1.0, 1.0)).min(axis=1)
    assert np.all(distances > EXTEND_THRESHOLD_MULTIPLIER * SPACING)
    # Growth is isotropic -- away from the neighbour too, not only along one lattice axis.
    i, j = _columns(a)
    assert i.min() < 10 and j.min() < 20 and j.max() > 29
    assert i.max() - 19 <= quad_tectonics.MAX_ADVANCE_LAYERS_PER_STEP


def test_two_plates_close_a_gap_between_them_from_both_sides():
    a = _plate(1, _block((10, 20), (20, 30)))
    b = _plate(2, _block((30, 40), (20, 30)))
    world = _world(a, b)

    a.deform(world, [b], 1_000_000, 0.0)
    b.deform(world, [a], 1_000_000, 0.0)

    ia, ja = _columns(a)
    ib, jb = _columns(b)
    row = (ja == 25)
    gap_cells = ib[jb == 25].min() - ia[row].max() - 1
    assert gap_cells <= 2


def test_new_cells_thin_the_cells_behind_them():
    # Continental, so the thinned margin stays above the rift-melting threshold instead of
    # erupting straight back to a fresh reference column.
    a = _plate(1, _block((10, 20), (20, 30)), "continental")
    world = _world(a)
    hc_before = dict(zip(map(int, a.cell_keys), a.collect("crustal_thickness_m")))

    a.deform(world, [], 1_000_000, 0.0)

    hc_after = dict(zip(map(int, a.cell_keys), a.collect("crustal_thickness_m")))
    edge = int(pack_cell_keys(0, 19, 25))
    interior = int(pack_cell_keys(0, 15, 25))
    assert hc_after[edge] < hc_before[edge]
    assert hc_after[interior] == pytest.approx(hc_before[interior])


def _volume(plate) -> float:
    return float(np.sum(plate.node_areas_m2() * plate.collect("crustal_thickness_m")))


def _grow(plate, world, neighbours, layers):
    n = plate.node_count()
    return quad_tectonics.grow_frontier(
        plate, world, np.ones(n, dtype=bool), np.zeros(n, dtype=bool), neighbours, SPACING, layers, 10_000
    )


def test_rift_opening_with_nothing_in_view_stretches_crust_without_losing_volume():
    a = _plate(1, _block((10, 20), (20, 30)), "continental")
    world = _world(a)
    volume_before = _volume(a)
    old = set(map(int, a.cell_keys))

    assert _grow(a, world, [], 2) > 0

    # Pure stretching: every new cell's crust came from the margin behind it.
    assert _volume(a) == pytest.approx(volume_before, rel=1e-9)
    hc = dict(zip(map(int, a.cell_keys), a.collect("crustal_thickness_m")))
    new = np.array([k not in old for k in map(int, a.cell_keys)])
    codes = a.collect("crust_type_code")
    assert np.all(codes[new] == CRUST_TYPE_CONTINENTAL)
    # A stretched margin, thinning outward: edge < one row in < untouched interior.
    reference = lithosphere.REFERENCE_HC_CONTINENTAL_M
    edge, inner, interior = (hc[int(pack_cell_keys(0, i, 25))] for i in (19, 18, 15))
    assert edge < inner < interior == pytest.approx(reference)
    assert hc[int(pack_cell_keys(0, 20, 25))] < edge


def test_rift_opening_never_thins_donor_hm_below_its_floor():
    """Issue #256: donors behind a rifted cell thin by `ratio` in both Hc and Hm. A thick
    crust over a thin mantle lid never melts through, so Hm used to drop past
    `MIN_MANTLE_LITHOSPHERE_THICKNESS_M` with nothing to reset it."""
    keys = _block((10, 20), (20, 30))
    thin_hm = 1.05 * lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M
    a = _plate(1, keys, "continental", mantle_lithosphere_thickness_m=np.full(len(keys), thin_hm))
    world = _world(a)

    assert _grow(a, world, [], 2) > 0

    hc, hm = a.collect("crustal_thickness_m"), a.collect("mantle_lithosphere_thickness_m")
    old = np.isin(a.cell_keys, keys)
    assert np.any(old & (hc < lithosphere.REFERENCE_HC_CONTINENTAL_M)), "some donor must have thinned"
    assert np.all(hm >= lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M)
    density = lithosphere.node_crust_density(a.collect("crust_type_code"), a.crust_type)
    np.testing.assert_allclose(a.collect("elevation")[old], lithosphere.isostatic_elevation(hc, hm, density)[old])


def test_ocean_ridge_between_separating_plates_accretes_fresh_crust():
    # Oceanic crust already near the rift threshold breaks up rather than stretching on.
    keys = _block((10, 20), (20, 30))
    a = _plate(1, keys, crustal_thickness_m=np.full(len(keys), 5_200.0))
    b = _plate(2, _block((30, 40), (20, 30)))
    world = _world(a, b)
    old = set(map(int, a.cell_keys))

    _grow(a, world, [b], 1)

    new = np.array([k not in old for k in map(int, a.cell_keys)])
    assert np.any(new)
    codes = a.collect("crust_type_code")
    hc = a.collect("crustal_thickness_m")
    assert np.all(codes[new] == CRUST_TYPE_OCEANIC)
    # The margin facing b stretched through the threshold and erupted a fresh column.
    erupted = a.collect("is_volcano") & ~new
    assert np.any(erupted)
    np.testing.assert_allclose(hc[erupted], lithosphere.REFERENCE_HC_OCEANIC_M)


def test_growth_across_the_separation_direction_is_mostly_magmatic():
    # b sits beyond a's +i edge, so cells a grows off its +/-j sides step across the
    # separation direction: little stretch share, mostly fresh ocean floor.
    a = _plate(1, _block((10, 20), (20, 30)), "continental")
    b = _plate(2, _block((24, 34), (20, 30)), "continental")
    world = _world(a, b)
    old = set(map(int, a.cell_keys))

    _grow(a, world, [b], 1)

    hc = dict(zip(map(int, a.cell_keys), a.collect("crustal_thickness_m")))
    codes = dict(zip(map(int, a.cell_keys), a.collect("crust_type_code")))
    side = int(pack_cell_keys(0, 15, 30))
    toward_b = int(pack_cell_keys(0, 20, 25))
    assert side not in old and toward_b not in old
    assert codes[side] == CRUST_TYPE_OCEANIC
    assert codes[toward_b] == CRUST_TYPE_CONTINENTAL
    assert hc[toward_b] > hc[side]


def test_over_budget_continental_plate_does_not_grow():
    thin = lithosphere.REFERENCE_HC_OCEANIC_M
    keys = _block((10, 20), (20, 30))
    a = _plate(1, keys, "continental", crustal_thickness_m=np.full(len(keys), thin))
    world = _world(a)

    a.deform(world, [], 1_000_000, 0.0)

    np.testing.assert_array_equal(a.cell_keys, keys)


# --- Retreat ------------------------------------------------------------------------------


def test_oceanic_plate_subducts_its_contested_overlap_layer_by_layer():
    a = _plate(1, _block((10, 24), (20, 30)))
    b = _plate(2, _block((20, 34), (20, 30)))
    world = _world(a, b)

    a.deform(world, [b], 1_000_000, 10 * SPACING)

    i, j = _columns(a)
    overlap_rows = (j >= 20) & (j < 30)
    assert i[overlap_rows].max() < 20
    assert not np.any(b.contains_batch(a.all_points_and_elevation()[0]))
    assert len(world.removed_points_log) == 40


def test_retreat_depth_is_capped_by_this_steps_displacement():
    a = _plate(1, _block((10, 24), (20, 30)))
    b = _plate(2, _block((20, 34), (20, 30)))
    world = _world(a, b)

    a.deform(world, [b], 1_000_000, 1.5 * SPACING)

    i, j = _columns(a)
    assert i[(j >= 20) & (j < 30)].max() == 22


def test_continental_suture_retreat_conserves_crustal_volume():
    # Thin enough that no survivor -- even one nearest several donors at the overlap's
    # corners -- reaches SUTURE_ACCRETION_MAX_HC_M, past which the excess delaminates.
    keys = _block((10, 24), (20, 30))
    a = _plate(1, keys, "continental", crustal_thickness_m=np.full(len(keys), 20_000.0))
    b = _plate(2, _block((20, 34), (20, 30)), "continental")
    world = _world(a, b)
    ctx = boundary_context(
        world, a, [b], 1_000_000,
        lambda contested: quad_tectonics.components_of_at_least(a, contested, CONTINENTAL_CONTESTED_RETREAT_MIN_RUN),
        node_weight=a.node_areas_m2() / lithosphere.node_area_m2(SPACING),
    )
    volume_before = float(np.sum(a.node_areas_m2() * a.collect("crustal_thickness_m")))
    elevation_before = a.collect("elevation").max()

    survivors = quad_tectonics._retreat(a, world, ctx, 1.5 * SPACING, 1000)

    assert not np.all(survivors)
    volume_after = float(np.sum(a.node_areas_m2() * a.collect("crustal_thickness_m")))
    assert volume_after == pytest.approx(volume_before, rel=1e-9)
    assert a.collect("elevation").max() > elevation_before


def test_suture_accretion_carries_cap_overflow_into_later_bands():
    keys = _block((10, 20), (20, 21))
    near_cap = SUTURE_ACCRETION_MAX_HC_M - 4_000.0
    hc = np.full(len(keys), near_cap)
    hc[0] = 20_000.0
    a = _plate(1, keys, "continental", crustal_thickness_m=hc)
    donors = np.zeros(len(keys), dtype=bool)
    donors[0] = True
    survivors = ~donors
    areas = a.node_areas_m2()
    expected = float(np.sum(hc[donors] * areas[donors]) + np.sum(hc[survivors] * areas[survivors]))

    quad_tectonics._accrete_onto_survivors(a, donors, survivors)

    after = a.collect("crustal_thickness_m")
    actual = float(np.sum(after[survivors] * areas[survivors]))
    assert actual == pytest.approx(expected, rel=1e-11)
    # The ordinary three-hop band only has 12 km of room, so most of the donor must reach
    # farther inland instead of disappearing at the cap.
    assert np.any(after[SUTURE_ACCRETION_SPREAD_NODES + 1 :] > near_cap)
    assert np.all(after[survivors] <= SUTURE_ACCRETION_MAX_HC_M)


def _strip_with_full_belts(length: int):
    """A one-cell-wide continental strip: a 30 km donor at one end, near-cap belts out to
    `SUTURE_ACCRETION_MAX_HOPS`, and reference crust with room beyond them."""
    keys = _block((10, 10 + length), (20, 21))
    hc = np.full(len(keys), lithosphere.REFERENCE_HC_CONTINENTAL_M)
    hc[1 : quad_tectonics.SUTURE_ACCRETION_MAX_HOPS + 1] = SUTURE_ACCRETION_MAX_HC_M - 1_500.0
    hc[0] = 30_000.0
    a = _plate(1, keys, "continental", crustal_thickness_m=hc)
    donors = np.zeros(len(keys), dtype=bool)
    donors[0] = True
    return a, donors


def test_a_cap_hit_alone_does_not_delaminate_suture_crust():
    a, donors = _strip_with_full_belts(30)
    world = _world(a)
    survivors = ~donors
    areas = a.node_areas_m2()
    hc = a.collect("crustal_thickness_m")
    expected = float(hc @ areas)

    # No time has passed, so no root can founder: everything the full belts turn away
    # reaches the reference crust beyond them.
    quad_tectonics._accrete_onto_survivors(a, donors, survivors, world)

    after = a.collect("crustal_thickness_m")
    assert float(after[survivors] @ areas[survivors]) == pytest.approx(expected, rel=1e-11)
    assert world.continental_material_ledger["delaminated_lower_crust_m3"] == 0.0
    budget = world.orogenic_relief_budget
    assert budget["far_field_placed_m3"] > 0.0
    assert budget["delamination_completed_m3"] == 0.0
    assert budget["no_outlet_delaminated_m3"] == 0.0


def test_suture_delamination_needs_an_eligible_root_and_is_bounded():
    a, donors = _strip_with_full_belts(30)
    world = _world(a)
    continental_ledger.ensure_initialized(world)
    survivors = ~donors
    areas = a.node_areas_m2()
    hc = a.collect("crustal_thickness_m")
    hm = a.collect("mantle_lithosphere_thickness_m")
    # Belts half continental-derived, donor wholly: the delaminated material's fraction then
    # tells whose crust sank.
    belts = np.zeros(len(hc), dtype=bool)
    belts[1 : quad_tectonics.SUTURE_ACCRETION_MAX_HOPS + 1] = True
    material = np.where(belts, 0.5 * hc, hc)
    a.set_fields_on_plate(continental_material_m=material)
    donated = float(hc[0] * areas[0])
    before = float(hc[survivors] @ areas[survivors])
    capacity = float(orogeny.delamination_capacity_m(hc[belts], hm[belts], 1.0) @ areas[belts])
    assert capacity > 0.0  # near-cap belts are hot and carry an eclogitic root

    quad_tectonics._accrete_onto_survivors(a, donors, survivors, world, years=1_000_000.0)

    lost = donated - (float(a.collect("crustal_thickness_m")[survivors] @ areas[survivors]) - before)
    budget = world.orogenic_relief_budget
    assert lost == pytest.approx(budget["delamination_completed_m3"], rel=1e-9)
    assert 0.0 < lost <= min(capacity, quad_tectonics.SUTURE_ACCRETION_MAX_DELAMINATION_FRACTION * donated) * (1 + 1e-9)
    # The belts' own roots sank, carrying their half-continental provenance; the donor's
    # crust and all its material were placed.
    delaminated = world.continental_material_ledger["delaminated_lower_crust_m3"]
    assert delaminated == pytest.approx(0.5 * lost, rel=1e-9)
    material_after = float(a.collect("continental_material_m")[survivors] @ areas[survivors])
    assert material_after + delaminated == pytest.approx(float(material @ areas), rel=1e-9)


def test_full_belts_without_an_eclogitic_root_never_delaminate():
    # A belt full at a 55 km cap has only 5 km of crust below the eclogite transition, too
    # thin to founder, so however much time passes the overflow must go elsewhere.
    a, donors = _strip_with_full_belts(30)
    hc = a.collect("crustal_thickness_m")
    hm = a.collect("mantle_lithosphere_thickness_m")
    cap = 55_000.0
    hc[1 : quad_tectonics.SUTURE_ACCRETION_MAX_HOPS + 1] = cap
    areas = a.node_areas_m2()
    root_capacity = orogeny.delamination_capacity_m(hc, hm, 100.0) * areas
    assert not np.any(root_capacity)
    volume = float(hc[0] * areas[0])

    _, stages = quad_tectonics._place_suture_crust(
        hc, areas, quad_tectonics._adjacency_matrix(a), a.surface_nodes().local_xyz,
        np.array([0]), ~donors, volume, cap, None, root_capacity,
    )

    assert stages["delamination_completed_m3"] == 0.0
    assert stages["far_field_placed_m3"] == pytest.approx(volume, rel=1e-9)


def test_suture_roots_shed_restite_first():
    a, donors = _strip_with_full_belts(30)
    belts = np.zeros(len(donors), dtype=bool)
    belts[1 : quad_tectonics.SUTURE_ACCRETION_MAX_HOPS + 1] = True
    restite = np.where(belts, 15_000.0, 0.0)
    a.set_fields_on_plate(restite_m=restite)
    world = _world(a)
    continental_ledger.ensure_initialized(world)
    areas = a.node_areas_m2()

    quad_tectonics._accrete_onto_survivors(a, donors, ~donors, world, years=1_000_000.0)

    budget = world.orogenic_relief_budget
    assert budget["delamination_completed_m3"] > 0.0
    assert budget["restite_delaminated_m3"] == pytest.approx(budget["delamination_completed_m3"], rel=1e-9)
    left = float(a.collect("restite_m")[belts] @ areas[belts])
    assert left == pytest.approx(float(restite @ areas) - budget["restite_delaminated_m3"], rel=1e-9)


def test_belts_whose_moho_still_lags_have_no_root_to_shed():
    a, donors = _strip_with_full_belts(30)
    a.set_fields_on_plate(moho_thermal_lag_c=np.full(len(donors), 400.0))
    world = _world(a)

    quad_tectonics._accrete_onto_survivors(a, donors, ~donors, world, years=1_000_000.0)

    assert world.orogenic_relief_budget["delamination_completed_m3"] == 0.0
    assert world.orogenic_relief_budget["far_field_placed_m3"] > 0.0


def test_accreted_crust_buries_its_receivers_moho():
    a, donors = _strip_with_full_belts(30)
    world = _world(a)
    hc_before = a.collect("crustal_thickness_m")
    moho_before = orogeny.moho_temperature_c(hc_before, a.collect("mantle_lithosphere_thickness_m"))

    quad_tectonics._accrete_onto_survivors(a, donors, ~donors, world)

    hc = a.collect("crustal_thickness_m")
    moho = orogeny.moho_temperature_c(hc, a.collect("mantle_lithosphere_thickness_m"), a.collect("moho_thermal_lag_c"))
    survivors = ~donors
    np.testing.assert_allclose(moho[survivors], moho_before[survivors], atol=1e-9)
    # The far field took stacked crust, so it now sits colder than its new steady state.
    far = survivors & (hc > hc_before + 1_000.0) & (hc_before < 40_000.0)
    assert np.any(far) and np.all(a.collect("moho_thermal_lag_c")[far] > 0.0)


def _capped_strip(length: int, overflowing: int = 2):
    keys = _block((10, 10 + length), (20, 21))
    hc = np.full(len(keys), lithosphere.REFERENCE_HC_CONTINENTAL_M)
    hc[:overflowing] = lithosphere.MAX_CRUSTAL_THICKNESS_M
    a = _plate(1, keys, "continental", crustal_thickness_m=hc, continental_material_m=hc.copy())
    overflow = np.zeros(len(keys))
    overflow[:overflowing] = 4_000.0
    return a, overflow


def test_ceiling_overflow_places_its_melt_and_books_its_residue():
    a, overflow = _capped_strip(30)
    world = _world(a)
    continental_ledger.ensure_initialized(world)
    areas = a.node_areas_m2()
    hc0 = a.collect("crustal_thickness_m")
    volume = float(overflow @ areas)
    fraction = orogeny.melt_fraction(orogeny.moho_temperature_c(hc0, a.collect("mantle_lithosphere_thickness_m")))
    melt = float((overflow * fraction) @ areas)
    assert 0.0 < melt < volume
    material = float(a.collect("continental_material_m") @ areas)

    quad_tectonics._place_ceiling_overflow(a, world, overflow, None)

    hc = a.collect("crustal_thickness_m")
    assert float(hc @ areas) == pytest.approx(float(hc0 @ areas) + melt, rel=1e-12)
    assert np.all(hc <= lithosphere.MAX_CRUSTAL_THICKNESS_M + 1e-6)
    assert hc[2] > hc0[2]
    # Shortening was never in the material tracer, and melt intrudes hot.
    assert float(a.collect("continental_material_m") @ areas) == pytest.approx(material, rel=1e-12)
    assert not np.any(a.collect("moho_thermal_lag_c"))
    budget = world.orogenic_relief_budget
    assert budget["ceiling_overflow_m3"] == pytest.approx(volume)
    assert budget["ceiling_overflow_melt_placed_m3"] == pytest.approx(melt)
    assert budget["ceiling_overflow_residue_m3"] == pytest.approx(volume - melt)
    assert budget["ceiling_overflow_no_outlet_m3"] == 0.0
    continental_ledger.assert_closed(world)


def test_a_cold_overflowing_column_has_no_melt_to_place():
    a, overflow = _capped_strip(30)
    a.set_fields_on_plate(moho_thermal_lag_c=np.full(a.node_count(), 600.0))
    world = _world(a)
    hc0 = a.collect("crustal_thickness_m")

    quad_tectonics._place_ceiling_overflow(a, world, overflow, None)

    np.testing.assert_array_equal(a.collect("crustal_thickness_m"), hc0)
    budget = world.orogenic_relief_budget
    assert budget["ceiling_overflow_residue_m3"] == pytest.approx(budget["ceiling_overflow_m3"])
    assert budget["ceiling_overflow_melt_placed_m3"] == 0.0


def _converging_ctx(n: int, convergent: np.ndarray, divergent: np.ndarray | None = None, arc: bool = False):
    from types import SimpleNamespace

    zeros = np.zeros(n, dtype=bool)
    divergent = zeros if divergent is None else divergent
    rate = 0.1 / (365.25 * 86400.0)
    return SimpleNamespace(
        convergent=convergent, divergent=divergent, transform=zeros,
        closing_rate=np.where(convergent, rate, np.where(divergent, -rate, 0.0)),
        inputs=SimpleNamespace(neighbor_is_oceanic=zeros),
        arc_band=convergent if arc else zeros, arc_intensity=np.where(convergent & arc, 1.0, 0.0),
        fault_influence=np.ones(n), orogen_dilation_nodes=3, orogen_contested_strength=1.0, orogen_amount=1.0,
        fault_noise=None, own_points=np.zeros((n, 3)),
    )


def test_quad_column_pass_reports_ceiling_overflow_instead_of_intruding_melt():
    from app.lithosphere_plate import COLUMN_FIELDS, deform_columns

    a, _ = _capped_strip(12)
    world = _world(a)
    continental_ledger.ensure_initialized(world)
    n = a.node_count()
    convergent = np.zeros(n, dtype=bool)
    convergent[:2] = True
    near_field_dist = quad_tectonics.hop_distance(a, convergent, 3)

    def run(overflow):
        fields = {name: a.collect(name) for name in COLUMN_FIELDS}
        return deform_columns(
            world, a, _converging_ctx(n, convergent), slice(None), fields, near_field_dist, lambda: None,
            a.node_areas_m2(), 0, 1_000_000.0, ceiling_overflow=overflow,
        )

    line_engine = run(None)
    overflow = np.zeros(n)
    quad = run(overflow)

    assert np.all(overflow[:2] > 0.0) and not np.any(overflow[2:])
    ring = (near_field_dist > 0) & (near_field_dist <= 3)
    # The line engine still intrudes part of the overflow as melt into the ring; the quad
    # engine leaves the ring to its own shortening and places the overflow itself.
    assert np.all(line_engine["crustal_thickness_m"][ring] > quad["crustal_thickness_m"][ring])
    np.testing.assert_array_equal(line_engine["crustal_thickness_m"][convergent], quad["crustal_thickness_m"][convergent])


def _thermal_state_after_columns(a, ctx):
    from app.lithosphere_plate import COLUMN_FIELDS, deform_columns

    world = _world(a)
    continental_ledger.ensure_initialized(world)
    fields = {name: a.collect(name) for name in COLUMN_FIELDS}
    strained: dict = {}
    columns = deform_columns(
        world, a, ctx, slice(None), fields, None, lambda: None, a.node_areas_m2(), 0, 1_000_000.0,
        ceiling_overflow=np.zeros(a.node_count()), strained=strained,
    )
    return fields, strained, columns, quad_tectonics._column_thermal_state(a, fields, strained, columns)


def test_shortening_buries_the_moho_but_arc_magma_on_top_does_not():
    keys = _block((10, 16), (20, 21))
    a = _plate(1, keys, "continental", restite_m=np.full(len(keys), 2_000.0))
    convergent = np.zeros(len(keys), dtype=bool)
    convergent[:3] = True

    before, strained, columns, thermal = _thermal_state_after_columns(a, _converging_ctx(len(keys), convergent, arc=True))

    hc0, hm0 = before["crustal_thickness_m"], before["mantle_lithosphere_thickness_m"]
    hc = columns["crustal_thickness_m"]
    shortened = strained["crustal_thickness_m"]
    assert np.all(shortened[convergent] > hc0[convergent])
    assert np.all(hc[convergent] > shortened[convergent])  # arc magma on top
    lag = thermal["moho_thermal_lag_c"]
    np.testing.assert_allclose(
        lag, orogeny.bury_moho(np.zeros(len(keys)), hc0, hm0, shortened, strained["mantle_lithosphere_thickness_m"])
    )
    assert np.all(lag[convergent] > 0.0) and not np.any(lag[~convergent])
    # The hot arc crust leaves the Moho warmer than shortening alone would.
    moho = orogeny.moho_temperature_c(hc, columns["mantle_lithosphere_thickness_m"], lag)
    assert np.all(moho[convergent] > orogeny.moho_temperature_c(hc0, hm0)[convergent])
    # Restite thickens with the strain, not with the arc crust on top.
    np.testing.assert_allclose(thermal["restite_m"], 2_000.0 * shortened / hc0)


def test_rift_thinning_exhumes_the_moho_into_a_negative_lag():
    keys = _block((10, 16), (20, 21))
    a = _plate(1, keys, "continental", restite_m=np.full(len(keys), 2_000.0))
    divergent = np.zeros(len(keys), dtype=bool)
    divergent[:3] = True

    before, strained, columns, thermal = _thermal_state_after_columns(
        a, _converging_ctx(len(keys), np.zeros(len(keys), dtype=bool), divergent)
    )

    hc0, hm0 = before["crustal_thickness_m"], before["mantle_lithosphere_thickness_m"]
    thinned_hc, thinned_hm = strained["crustal_thickness_m"], strained["mantle_lithosphere_thickness_m"]
    assert np.all(thinned_hc[divergent] < hc0[divergent])
    lag = thermal["moho_thermal_lag_c"]
    assert np.all(lag[divergent] < 0.0) and not np.any(lag[~divergent])
    np.testing.assert_allclose(
        orogeny.moho_temperature_c(thinned_hc, thinned_hm, lag)[divergent], orogeny.moho_temperature_c(hc0, hm0)[divergent]
    )
    assert np.all(thermal["restite_m"][divergent] < 2_000.0)


def test_tectonic_escape_moves_crust_along_strike_not_inland():
    # A short suture front along j at i = 10; the overriding plate lies toward -i.
    keys = _block((10, 40), (0, 50))
    a = _plate(1, keys, "continental")
    i, j = _columns(a)
    donors = (i == 10) & (j >= 23) & (j < 28)
    hops = quad_tectonics.hop_distance(a, donors, 40)
    hc = a.collect("crustal_thickness_m")
    hc[(hops > 0) & (hops <= quad_tectonics.SUTURE_ACCRETION_MAX_HOPS)] = SUTURE_ACCRETION_MAX_HC_M
    hc[donors] = 60_000.0
    a.set_fields_on_plate(crustal_thickness_m=hc)
    points = a.surface_nodes().local_xyz
    toward_overrider = geometry.normalize(points[(i == 10) & (j == 25)][0] - points[(i == 11) & (j == 25)][0])
    convergence = np.tile(toward_overrider, (a.node_count(), 1))
    world = _world(a)
    areas = a.node_areas_m2()
    expected = float(hc @ areas)

    quad_tectonics._accrete_onto_survivors(a, donors, ~donors, world, convergence_xyz=convergence)

    after = a.collect("crustal_thickness_m")
    assert float(after[~donors] @ areas[~donors]) == pytest.approx(expected, rel=1e-11)
    gained = (after > hc + 1e-6) & ~donors
    assert world.orogenic_relief_budget["escape_placed_m3"] > 0.0
    assert np.any(gained)
    # Every receiver lies within 45 degrees of strike as seen from the front's centre, and
    # the cells straight inland of the front receive nothing.
    assert np.all(np.abs(j[gained] - 25) >= i[gained] - 10)
    assert not np.any(gained & (j >= 23) & (j < 28))


def test_a_saturated_plate_books_its_suture_remainder_as_terminal_delamination():
    keys = _block((10, 20), (20, 21))
    hc = np.full(len(keys), SUTURE_ACCRETION_MAX_HC_M)
    hc[0] = 30_000.0
    a = _plate(1, keys, "continental", crustal_thickness_m=hc)
    world = _world(a)
    continental_ledger.ensure_initialized(world)
    donors = np.zeros(len(keys), dtype=bool)
    donors[0] = True
    areas = a.node_areas_m2()

    quad_tectonics._accrete_onto_survivors(a, donors, ~donors, world)

    donated = float(hc[0] * areas[0])
    assert np.all(a.collect("crustal_thickness_m")[~donors] == SUTURE_ACCRETION_MAX_HC_M)
    assert world.orogenic_relief_budget["no_outlet_delaminated_m3"] == pytest.approx(donated)
    assert world.continental_material_ledger["delaminated_lower_crust_m3"] == pytest.approx(donated)


def test_a_saturated_continent_spills_suture_crust_onto_its_oceanic_margin():
    keys = _block((10, 24), (20, 21))
    a = _plate(1, keys, "continental")
    hc = a.collect("crustal_thickness_m")
    codes = a.collect("crust_type_code")
    hc[:10] = SUTURE_ACCRETION_MAX_HC_M
    hc[0] = 40_000.0
    codes[10:] = CRUST_TYPE_OCEANIC  # a drowned margin beyond the full continent
    hc[10:] = lithosphere.REFERENCE_HC_OCEANIC_M
    a.set_fields_on_plate(crustal_thickness_m=hc, crust_type_code=codes, continental_material_m=np.where(codes == CRUST_TYPE_OCEANIC, 0.0, hc))
    world = _world(a)
    continental_ledger.ensure_initialized(world)
    donors = np.zeros(len(keys), dtype=bool)
    donors[0] = True
    areas = a.node_areas_m2()
    material_before = float(a.collect("continental_material_m") @ areas)

    quad_tectonics._accrete_onto_survivors(a, donors, ~donors, world)

    after = a.collect("crustal_thickness_m")
    assert float(after[~donors] @ areas[~donors]) == pytest.approx(float(hc @ areas), rel=1e-11)
    budget = world.orogenic_relief_budget
    assert budget["foreland_spill_placed_m3"] == pytest.approx(float(hc[0] * areas[0]), rel=1e-9)
    assert budget["no_outlet_delaminated_m3"] == 0.0
    assert world.continental_material_ledger["delaminated_lower_crust_m3"] == 0.0
    material = a.collect("continental_material_m")
    assert float(material[~donors] @ areas[~donors]) == pytest.approx(material_before, rel=1e-9)
    # The nearest margin cell took most of it and now carries a continental column.
    assert a.collect("crust_type_code")[10] == CRUST_TYPE_CONTINENTAL


def test_a_consumed_saturated_plate_accretes_its_crust_onto_the_overriding_plate():
    # A small continent with every cell at the cap, overridden along its i = 13 edge.
    a = _plate(1, _block((10, 14), (20, 24)), "continental",
               crustal_thickness_m=np.full(16, SUTURE_ACCRETION_MAX_HC_M))
    b = _plate(2, _block((13, 40), (10, 34)), "continental")
    world = _world(a, b)
    continental_ledger.ensure_initialized(world)
    i, _ = _columns(a)
    donors = i == 13
    a_areas, b_areas = a.node_areas_m2(), b.node_areas_m2()
    donated = float(a.collect("crustal_thickness_m")[donors] @ a_areas[donors])
    b_before = float(b.collect("crustal_thickness_m") @ b_areas)
    material_before = continental_ledger.surface_volume_m3(world)

    quad_tectonics._accrete_onto_survivors(a, donors, ~donors, world, overriders=[b])

    gained = float(b.collect("crustal_thickness_m") @ b_areas) - b_before
    assert gained == pytest.approx(donated, rel=1e-9)
    budget = world.orogenic_relief_budget
    assert budget["overrider_placed_m3"] == pytest.approx(donated, rel=1e-9)
    assert budget["no_outlet_delaminated_m3"] == 0.0
    assert world.continental_material_ledger["delaminated_lower_crust_m3"] == 0.0
    # The donor cells are still on `a` here (retreat removes them next), so the material
    # they hand over is counted twice until then: once on `a`, once on `b`.
    handed_material = float(a.collect("continental_material_m")[donors] @ a_areas[donors])
    assert continental_ledger.surface_volume_m3(world) == pytest.approx(material_before + handed_material, rel=1e-9)
    # It lands next to the suture, not across the plate.
    bi, _ = _columns(b)
    grew = b.collect("crustal_thickness_m") > lithosphere.REFERENCE_HC_CONTINENTAL_M + 1.0
    assert bi[grew].min() <= 13 + quad_tectonics.SUTURE_ACCRETION_MAX_HOPS


def test_the_overriding_plate_sheds_its_own_roots_not_the_handed_crust():
    a = _plate(1, _block((10, 14), (20, 24)), "continental",
               crustal_thickness_m=np.full(16, SUTURE_ACCRETION_MAX_HC_M))
    # A hot overrider whose whole continent is near the cap, half continental-derived.
    b_hc = np.full(27 * 24, SUTURE_ACCRETION_MAX_HC_M - 500.0)
    b = _plate(2, _block((13, 40), (10, 34)), "continental", crustal_thickness_m=b_hc)
    b.set_fields_on_plate(continental_material_m=0.5 * b_hc)
    world = _world(a, b)
    continental_ledger.ensure_initialized(world)
    i, _ = _columns(a)
    donors = i == 13
    a_areas, b_areas = a.node_areas_m2(), b.node_areas_m2()
    donated = float(a.collect("crustal_thickness_m")[donors] @ a_areas[donors])
    donor_material = float(a.collect("continental_material_m")[donors] @ a_areas[donors])
    b_material_before = float(b.collect("continental_material_m") @ b_areas)
    b_hc_before = float(b.collect("crustal_thickness_m") @ b_areas)

    quad_tectonics._accrete_onto_survivors(a, donors, ~donors, world, years=1_000_000.0, overriders=[b])

    budget = world.orogenic_relief_budget
    handed = budget["overrider_placed_m3"]
    assert handed > 0.0
    # Whatever `b` didn't net-gain of what it was handed is root it shed to make room.
    b_shed = handed - (float(b.collect("crustal_thickness_m") @ b_areas) - b_hc_before)
    assert b_shed > 0.0
    # `b` lost half-continental root and kept all the handed crust's material (wholly
    # continental, in proportion to what was handed).
    handed_material = donor_material * handed / donated
    b_material_after = float(b.collect("continental_material_m") @ b_areas)
    assert b_material_after == pytest.approx(b_material_before + handed_material - 0.5 * b_shed, rel=1e-9)


def test_suture_accretion_conserves_volume_across_unequal_area_cells():
    keys = _block((10, 30), (20, 22))
    a = _plate(1, keys, "continental")
    a.refine_cells(a.cell_keys[10:14])
    i, j = _columns(a)
    hc = a.collect("crustal_thickness_m")
    donors = i == 10
    hc[~donors] = SUTURE_ACCRETION_MAX_HC_M - 5_000.0
    a.set_fields_on_plate(crustal_thickness_m=hc, continental_material_m=hc.copy())
    world = _world(a)
    areas = a.node_areas_m2()
    assert np.ptp(areas) > 0.5 * areas.max()

    quad_tectonics._accrete_onto_survivors(a, donors, ~donors, world)

    after = a.collect("crustal_thickness_m")
    assert float(after[~donors] @ areas[~donors]) == pytest.approx(float(hc @ areas), rel=1e-11)
    material = a.collect("continental_material_m")
    assert float(material[~donors] @ areas[~donors]) == pytest.approx(float(hc @ areas), rel=1e-9)


def test_suture_accretion_reaches_a_distant_same_type_survivor_when_local_band_is_empty():
    keys = _block((10, 30), (20, 21))
    a = _plate(1, keys, "continental")
    donors = np.zeros(len(keys), dtype=bool)
    donors[0] = True
    survivors = np.zeros(len(keys), dtype=bool)
    survivors[-1] = True  # 19 hops away, beyond SUTURE_ACCRETION_MAX_HOPS.
    areas = a.node_areas_m2()
    hc_before = a.collect("crustal_thickness_m")
    hm_before = a.collect("mantle_lithosphere_thickness_m")
    expected = float(hc_before[donors] @ areas[donors] + hc_before[survivors] @ areas[survivors])
    expected_hm = float(hm_before[donors] @ areas[donors] + hm_before[survivors] @ areas[survivors])

    quad_tectonics._accrete_onto_survivors(a, donors, survivors)

    actual = float(a.collect("crustal_thickness_m")[survivors] @ areas[survivors])
    assert actual == pytest.approx(expected, rel=1e-11)
    actual_hm = float(a.collect("mantle_lithosphere_thickness_m")[survivors] @ areas[survivors])
    assert actual_hm == pytest.approx(expected_hm, rel=1e-11)


def test_accretion_falls_back_to_any_type_when_no_same_type_survivor_exists():
    keys = _block((10, 20), (20, 21))
    a = _plate(1, keys, "continental")
    codes = a.collect("crust_type_code")
    codes[0] = CRUST_TYPE_OCEANIC
    a.set_fields_on_plate(crust_type_code=codes)
    donors = np.zeros(len(keys), dtype=bool)
    donors[0] = True
    survivors = ~donors
    areas = a.node_areas_m2()
    hc_before = a.collect("crustal_thickness_m")
    hm_before = a.collect("mantle_lithosphere_thickness_m")
    expected = float(np.dot(hc_before, areas))
    expected_hm = float(np.dot(hm_before, areas))

    quad_tectonics._accrete_onto_survivors(a, donors, survivors)

    actual = float(a.collect("crustal_thickness_m")[survivors] @ areas[survivors])
    assert actual == pytest.approx(expected, rel=1e-11)
    actual_hm = float(a.collect("mantle_lithosphere_thickness_m")[survivors] @ areas[survivors])
    assert actual_hm == pytest.approx(expected_hm, rel=1e-11)


def test_oceanic_plate_retreat_accretes_continental_terrane_onto_terrane_survivors():
    keys = _block((10, 24), (20, 30))
    a = _plate(1, keys, "oceanic")
    i, _ = _columns(a)
    codes = a.collect("crust_type_code")
    codes[i >= 18] = CRUST_TYPE_CONTINENTAL
    a.set_fields_on_plate(crust_type_code=codes)
    b = _plate(2, _block((20, 34), (20, 30)), "continental")
    world = _world(a, b)
    ctx = boundary_context(
        world,
        a,
        [b],
        1_000_000,
        lambda contested: quad_tectonics.components_of_at_least(a, contested, CONTINENTAL_CONTESTED_RETREAT_MIN_RUN),
        node_weight=a.node_areas_m2() / lithosphere.node_area_m2(SPACING),
    )
    continental = a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL
    volume_before = float(np.sum((a.collect("crustal_thickness_m") * a.node_areas_m2())[continental]))

    survivors = quad_tectonics._retreat(a, world, ctx, 1.5 * SPACING, 1_000)

    assert not np.all(survivors)
    continental = a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL
    volume_after = float(np.sum((a.collect("crustal_thickness_m") * a.node_areas_m2())[continental]))
    assert volume_after == pytest.approx(volume_before, rel=1e-9)


def test_fully_consumed_terrane_relocates_without_relabeling_oceanic_volume():
    keys = _block((10, 20), (20, 21))
    a = _plate(1, keys, "oceanic")
    codes = a.collect("crust_type_code")
    codes[0] = CRUST_TYPE_CONTINENTAL
    a.set_fields_on_plate(crust_type_code=codes)
    donors = np.zeros(len(keys), dtype=bool)
    donors[0] = True
    survivors = ~donors
    areas = a.node_areas_m2()
    hc_before = a.collect("crustal_thickness_m")
    continental_before = float(hc_before[0] * areas[0])
    total_before = float(np.dot(hc_before, areas))

    quad_tectonics._accrete_onto_survivors(a, donors, survivors)

    hc_after = a.collect("crustal_thickness_m")
    continental = a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL
    assert float(np.dot(hc_after[continental & survivors], areas[continental & survivors])) == pytest.approx(
        continental_before, rel=1e-11
    )
    assert float(np.dot(hc_after[survivors], areas[survivors])) == pytest.approx(total_before, rel=1e-11)


def test_disconnected_terrane_relocation_takes_only_the_footprint_it_needs():
    keys = np.sort(np.concatenate([_block((10, 11), (10, 11)), _block((30, 35), (30, 31))]))
    a = _plate(1, keys, "oceanic")
    codes = a.collect("crust_type_code")
    donor = int(a._index_of_keys(_block((10, 11), (10, 11)))[0])
    codes[donor] = CRUST_TYPE_CONTINENTAL
    a.set_fields_on_plate(crust_type_code=codes)
    donors = np.zeros(len(keys), dtype=bool)
    donors[donor] = True
    survivors = ~donors
    areas = a.node_areas_m2()
    total_before = float(np.dot(a.collect("crustal_thickness_m"), areas))
    hm_before = float(np.dot(a.collect("mantle_lithosphere_thickness_m"), areas))

    quad_tectonics._accrete_onto_survivors(a, donors, survivors)

    continental = a.collect("crust_type_code") == CRUST_TYPE_CONTINENTAL
    assert np.count_nonzero(continental & survivors) == 1
    assert np.any(~continental & survivors)
    total_after = float(a.collect("crustal_thickness_m")[survivors] @ areas[survivors])
    assert total_after == pytest.approx(total_before, rel=1e-11)
    hm_after = float(a.collect("mantle_lithosphere_thickness_m")[survivors] @ areas[survivors])
    assert hm_after == pytest.approx(hm_before, rel=1e-11)


def _plate_with_overlay(crust_type):
    a = _plate(1, _block((10, 30), (10, 30)), crust_type)
    # b has slid onto a just inside a's edge: a's own edge rows in front of it stay uncovered,
    # so no contested cell of a touches open ground.
    b = _plate(2, _block((12, 16), (18, 22)), "continental")
    return a, b


def test_oceanic_plate_carves_out_an_interior_overlap_the_edge_peel_cannot_reach():
    a, b = _plate_with_overlay("oceanic")
    world = _world(a, b)
    covered = b.contains_batch(a.all_points_and_elevation()[0])
    assert covered.sum() == 16

    a.deform(world, [b], 1_000_000, SPACING)

    assert not np.any(b.contains_batch(a.all_points_and_elevation()[0]))
    assert len(world.removed_points_log) == 16
    a._validate_leaf_topology()
    # The carve leaves a hole: a now has an inner boundary loop.
    assert len(a.boundary_loops_world()) == 2


def test_continental_plate_is_never_carved_mid_plate():
    a, b = _plate_with_overlay("continental")
    world = _world(a, b)
    covered = b.contains_batch(a.all_points_and_elevation()[0])

    a.deform(world, [b], 1_000_000, SPACING)

    assert b.contains_batch(a.all_points_and_elevation()[0]).sum() == covered.sum() == 16


def _interior_patch():
    """A 10x10 oceanic plate with a 4x4 contested patch two cells in from every edge, and
    an `open_half` with nothing peeled."""
    a = _plate(1, _block((10, 20), (10, 20)))
    i, j = _columns(a)
    patch = (i >= 13) & (i < 17) & (j >= 13) & (j < 17)
    open_half = a._probe_neighbour_indices() < 0
    return a, patch, open_half


def test_interior_carve_opens_a_connected_hole_in_a_patch_larger_than_its_budget():
    a, patch, open_half = _interior_patch()

    carved = quad_tectonics._carve_interior(a, patch, open_half, 5)

    assert carved.sum() == 5
    assert not np.any(carved & ~patch)
    assert len(quad_tectonics.components_of_at_least(a, carved, 5).nonzero()[0]) == 5


def test_interior_carve_treats_a_half_open_side_as_unreachable():
    a, patch, open_half = _interior_patch()
    # One probe on one side of a patch cell opens onto nothing -- a coarse/fine partial
    # boundary. The peel needs a wholly open side, so it can't take this cell either.
    cell = np.flatnonzero(patch)[0]
    open_half[cell, 0, 0] = True

    carved = quad_tectonics._carve_interior(a, patch, open_half, 100)

    np.testing.assert_array_equal(carved, patch)


def test_a_lone_contested_continental_cell_does_not_retreat():
    a = _plate(1, _block((10, 20), (20, 30)), "continental")
    # One cell of b pokes onto a's edge.
    b_keys = np.concatenate([_block((20, 30), (20, 30)), pack_cell_keys(np.array([0]), np.array([19]), np.array([25]))])
    b = _plate(2, np.sort(b_keys), "continental")
    world = _world(a, b)

    a.deform(world, [b], 1_000_000, 10 * SPACING)

    assert int(pack_cell_keys(0, 19, 25)) in set(map(int, a.cell_keys))


# --- Other step passes ----------------------------------------------------------------------


def test_quad_volcano_erupts_through_the_surface_field_api():
    keys = _block((10, 14), (20, 24))
    count = len(keys)
    is_volcano = np.zeros(count, dtype=bool)
    is_volcano[5] = True
    a = _plate(1, keys, is_volcano=is_volcano, volcano_active_years_remaining=np.where(is_volcano, 5e6, 0.0))
    world = _world(a)
    world.volcanism_multiplier = 1000.0
    hc_before = a.collect("crustal_thickness_m")

    volcanism.apply_volcanic_activity(world, 1_000_000)

    hc_after = a.collect("crustal_thickness_m")
    assert hc_after[5] > hc_before[5]
    assert a.collect("mineral_deposit_m")[5] > 0.0
    assert a.collect("volcano_active_years_remaining")[5] == pytest.approx(4e6)


def test_gap_fill_grows_an_adjacent_quad_plate_into_the_gap():
    whole = _plate(1, np.concatenate([_block((0, N), (0, N), face) for face in range(6)]))
    points = whole.all_points_and_elevation()[0]
    hole = points @ np.array([1.0, 0.0, 0.0]) > np.cos(12 * SPACING)
    whole.remove_cells(hole)
    world = _world(whole)
    count_before = whole.node_count()

    gaps.fill_gaps_by_growing_neighbours(world)

    assert len(world.plates) == 1
    assert whole.node_count() > count_before
    whole._validate_leaf_topology()


def test_gap_spawn_in_a_quad_world_makes_a_quad_plate():
    plate = new_plate(
        5, np.eye(3), "oceanic", SPACING, seed=0,
        is_owned=lambda pts: pts @ np.array([1.0, 0.0, 0.0]) > np.cos(5 * SPACING),
        node_is_continental=lambda pts: pts[:, 1] > 0.0,
        surface="quad",
    )

    assert isinstance(plate, PlateWithSparseQuadPatch)
    assert plate.node_count() > 0
    codes = plate.collect("crust_type_code")
    assert np.any(codes == CRUST_TYPE_CONTINENTAL) and np.any(codes != CRUST_TYPE_CONTINENTAL)
    assert np.all(plate.collect("crustal_thickness_m") > 0.0)


def test_a_quad_pair_past_its_forced_merge_time_comes_due():
    # Quad merge (quad_merge.py) lifted the guard that kept quad collisions from ending.
    a = _plate(1, _block((10, 20), (20, 30)), "continental")
    b = _plate(2, _block((20, 30), (20, 30)), "continental")
    world = _world(a, b)
    world.overlap_progress[(1, 2)] = merge_split.FORCED_MERGE_SUSTAINED_YEARS * 2

    assert merge_split._supports_merge(world, 1, 2)
    forced = merge_split.pop_ready_forced_merge(world, can_merge=lambda x, y: merge_split._supports_merge(world, x, y))

    assert set(forced) == {1, 2}


def test_mostly_magmatic_new_cells_are_volcanic_vents_and_stretched_ones_are_not():
    # Same layout as the test above: the +/-j side cells are fresh magma, the cells toward b
    # are stretched continental margin.
    a = _plate(1, _block((10, 20), (20, 30)), "continental")
    b = _plate(2, _block((24, 34), (20, 30)), "continental")
    world = _world(a, b)

    _grow(a, world, [b], 1)

    index = {k: n for n, k in enumerate(map(int, a.cell_keys))}
    side, toward_b = index[int(pack_cell_keys(0, 15, 30))], index[int(pack_cell_keys(0, 20, 25))]
    is_volcano = a.collect("is_volcano")
    remaining = a.collect("volcano_active_years_remaining")
    assert is_volcano[side] and remaining[side] > 0.0
    assert not is_volcano[toward_b]
    # Igniting leaves the column alone: the vent keeps the fresh magmatic crust it was grown with.
    assert a.collect("crust_type_code")[side] == CRUST_TYPE_OCEANIC
    assert a.collect("crustal_thickness_m")[side] > 5_000.0


def test_quad_eruption_draw_is_keyed_by_cell_not_node_order():
    keys = _block((10, 14), (20, 24))
    count = len(keys)
    vent = 9
    is_volcano = np.zeros(count, dtype=bool)
    is_volcano[vent] = True
    fields = dict(is_volcano=is_volcano, volcano_active_years_remaining=np.where(is_volcano, 5e8, 0.0))
    a = _plate(1, keys, **fields)
    # The same plate with extra cells ahead of the vent in node order.
    b = _plate(1, keys, **fields)
    b.insert_cells(_block((10, 14), (18, 20)))
    world_a, world_b = _world(a), _world(b)
    vent_key = int(keys[vent])
    outcomes = []
    for step in range(40):
        for world in (world_a, world_b):
            world.elapsed_years = step * 1_000_000.0
            # ~50% per step: p = 1 - exp(-rate * multiplier * 1 Myr).
            world.volcanism_multiplier = np.log(2.0) / volcanism.ERUPTION_RATE_PER_MYR
        hc_a = a.collect("crustal_thickness_m")[a._index_of_keys(np.array([vent_key]))[0]]
        hc_b = b.collect("crustal_thickness_m")[b._index_of_keys(np.array([vent_key]))[0]]
        volcanism.apply_volcanic_activity(world_a, 1_000_000)
        volcanism.apply_volcanic_activity(world_b, 1_000_000)
        erupted_a = a.collect("crustal_thickness_m")[a._index_of_keys(np.array([vent_key]))[0]] > hc_a
        erupted_b = b.collect("crustal_thickness_m")[b._index_of_keys(np.array([vent_key]))[0]] > hc_b
        assert erupted_a == erupted_b
        outcomes.append(erupted_a)
    assert any(outcomes) and not all(outcomes)


def _misaligned_frame() -> np.ndarray:
    """A frame whose face-0 cells sit half a cell off the identity frame's along both axes, so
    each node is ~0.7 spacings from the nearest node of an identity-frame plate: farther than
    the node-proximity overlap tolerance, though well inside that plate's cells."""
    half_cell = 0.5 * (np.pi / 2.0) / N
    z = geometry.rotation_matrix(np.array([0.0, 0.0, 1.0]), half_cell)
    y = geometry.rotation_matrix(np.array([0.0, 1.0, 0.0]), half_cell)
    return z @ y


def test_quad_overlap_is_read_by_containment():
    a = _plate(1, _block((10, 30), (20, 40)))
    b = _plate(2, _block((26, 40), (20, 40)), frame=_misaligned_frame())
    world = _world(a, b)
    tol = plates.OVERLAP_TOLERANCE_MULT * SPACING

    overlap = plates.compute_node_overlap(world.plates, tol)

    # The node-proximity reading line plates use misses this overlap entirely.
    points = [p.all_points_and_elevation()[0] for p in (a, b)]
    assert np.arccos(np.clip(points[0] @ points[1].T, -1.0, 1.0)).min() > tol
    for plate, other in ((a, b), (b, a)):
        inside = other.contains_batch(plate.all_points_and_elevation()[0])
        assert np.any(inside)
        np.testing.assert_array_equal(overlap[plate.plate_id]["overlap_mask"], inside)
        np.testing.assert_array_equal(overlap[plate.plate_id]["cover_count"], inside.astype(int))
        # Both plates are oceanic, so neither contributes to a continental cover count.
        assert not overlap[plate.plate_id]["continental_cover_count"].any()
        assert overlap[plate.plate_id]["by_partner"] == {other.plate_id: int(inside.sum())}


def test_quad_overlap_sees_a_plate_buried_inside_another():
    big = _plate(1, _block((0, N), (0, N)), "continental")
    buried = _plate(2, _block((N // 2 - 3, N // 2 + 3), (N // 2 - 3, N // 2 + 3)), "continental", frame=_misaligned_frame())
    world = _world(big, buried)

    overlap = plates.compute_node_overlap(world.plates, plates.OVERLAP_TOLERANCE_MULT * SPACING)

    assert overlap[2]["overlap_mask"].all()
    assert overlap[1]["by_partner"][2] > 0
    assert (overlap[2]["continental_cover_count"] == 1).all()


def test_overlap_tracking_stamps_onset_on_quad_plates():
    a = _plate(1, _block((10, 30), (20, 40)))
    b = _plate(2, _block((26, 40), (20, 40)), frame=_misaligned_frame())
    world = _world(a, b)
    world.elapsed_years = 7e6

    merge_split.update_overlap_tracking(world, 1_000_000)

    inside = a.contains_batch(b.all_points_and_elevation()[0])
    onset = b.collect("overlap_onset_years")
    np.testing.assert_array_equal(onset > 0.0, inside)
    assert np.all(onset[inside] == 7e6)


def _seam_world():
    """Two plates one empty column apart -- the seam boundary advance leaves, since a cell
    there is within EXTEND_THRESHOLD_MULTIPLIER spacings of the neighbour's nodes -- and a
    third plate covering the rest of the sphere, so the seam is the only uncovered ground."""
    a = _plate(1, _block((10, 20), (20, 40)))
    b = _plate(2, _block((21, 31), (20, 40)))
    seam = _block((20, 21), (20, 40))
    everything = np.concatenate([_block((0, N), (0, N), face) for face in range(6)])
    rest = _plate(3, everything[~np.isin(everything, np.concatenate([a.cell_keys, b.cell_keys, seam]))])
    return a, b, _world(a, b, rest), seam


def test_boundary_advance_leaves_a_one_cell_seam():
    a, b, world, seam = _seam_world()
    a, b = (_plate(p.plate_id, p.cell_keys) for p in (a, b))
    world = _world(a, b)

    a.deform(world, [b], 1_000_000, 0.0)
    b.deform(world, [a], 1_000_000, 0.0)

    assert not np.any(a._index_of_keys(seam) >= 0) and not np.any(b._index_of_keys(seam) >= 0)


def test_gap_fill_closes_a_one_cell_seam_between_quad_plates():
    a, b, world, seam = _seam_world()
    context = gaps._existing_node_tree(world)
    # Node distance can't see the seam; containment finds it and nothing else.
    assert len(gaps._find_gap_points(context, SPACING)) == 0
    found = gaps._find_gap_points(context, SPACING, world.plates)
    assert len(found) > 0 and np.all(_plate(9, seam).contains_batch(found))

    gaps.fill_gaps_by_growing_neighbours(world)

    seam_points = geometry.to_world(np.eye(3), a.cell_centres_local(seam))
    covered = np.zeros(len(seam_points), dtype=bool)
    for plate in world.plates:
        covered |= plate.contains_batch(seam_points)
        plate._validate_leaf_topology()
    assert covered.all()
    # Filled, not overlapped: no plate took a cell another holds.
    overlap = plates.compute_node_overlap(world.plates, plates.OVERLAP_TOLERANCE_MULT * SPACING)
    assert not any(info["overlap_mask"].any() for info in overlap.values())


def test_gap_fill_requires_most_of_a_candidate_cell_to_be_uncovered():
    plate = _plate(1, _block((10, 11), (20, 21)))
    candidates = _block((11, 13), (20, 21))

    class SampleMask:
        def contains_batch(self, points):
            covered = np.zeros(len(points), dtype=bool)
            # Candidate 0 is 9/16 covered and must be rejected; candidate 1 is 7/16 covered
            # and must remain claimable. The helper batches 16 footprint probes per cell.
            covered[:9] = True
            covered[16:23] = True
            return covered

    np.testing.assert_array_equal(
        quad_tectonics._gap_cells_mostly_uncovered(plate, candidates, [SampleMask()]),
        [False, True],
    )


def test_a_small_isolated_quad_gap_is_not_spawned_into_a_plate():
    whole = _plate(1, np.concatenate([_block((0, N), (0, N), face) for face in range(6)]))
    points = whole.all_points_and_elevation()[0]
    # A hole well under MIN_GAP_NODES, deep enough that nothing lies within
    # ADJACENT_PLATE_REACH_MULT spacings of its middle... except its own rim, so the plate
    # grows into it rather than a new plate spawning there.
    whole.remove_cells(points @ np.array([1.0, 0.0, 0.0]) > np.cos(3 * SPACING))
    world = _world(whole)

    gaps.fill_gaps_by_growing_neighbours(world)

    assert len(world.plates) == 1
    assert len(gaps._find_gap_points(gaps._existing_node_tree(world), SPACING, world.plates)) == 0


def test_gap_tracks_on_a_quad_world_see_the_seam():
    _, _, world, _ = _seam_world()

    gaps.reconcile_gap_tracks(world)

    assert len(world.gap_tracks) >= 1


def test_relattice_leaves_quad_plates_alone():
    keys = _block((10, 20), (20, 30))
    a = _plate(1, keys, "continental", crustal_thickness_m=np.linspace(30_000.0, 40_000.0, len(keys)))
    world = _world(a)
    world.steps_taken = merge_split.RELATTICE_INTERVAL_STEPS
    hc_before = a.collect("crustal_thickness_m")

    merge_split.relattice_continental_plates(world)

    np.testing.assert_array_equal(a.cell_keys, keys)
    np.testing.assert_array_equal(a.collect("crustal_thickness_m"), hc_before)
