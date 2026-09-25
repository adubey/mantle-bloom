"""Per-step deformation on sparse quad plates (issue #228 Phase 4): the cell-graph topology
primitives, 2D boundary retreat/advance, and the step passes that had to learn about quad
plates (volcanism, gap filling, merge) -- see quad_tectonics.py."""

import numpy as np
import pytest

from app import gaps, lithosphere, merge_split, quad_tectonics, volcanism
from app.elevation_lines import CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC, line_spacing_rad
from app.lithosphere_plate import (
    CONTINENTAL_CONTESTED_RETREAT_MIN_RUN,
    EXTEND_THRESHOLD_MULTIPLIER,
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


def _plate(plate_id, keys, crust_type="oceanic", **fields) -> PlateWithSparseQuadPatch:
    count = len(keys)
    hc, hm = lithosphere.reference_thickness(crust_type)
    defaults = {
        "crustal_thickness_m": np.full(count, hc),
        "mantle_lithosphere_thickness_m": np.full(count, hm),
    }
    defaults.update(fields)
    plate = PlateWithSparseQuadPatch(plate_id, np.eye(3), crust_type, N, keys, fields=defaults)
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


def test_quad_plates_are_not_offered_to_merge():
    a = _plate(1, _block((10, 20), (20, 30)), "continental")
    b = _plate(2, _block((20, 30), (20, 30)), "continental")

    assert not merge_split._supports_merge(_world(a, b), 1, 2)


def test_forced_merge_keeps_an_unmergeable_pairs_overlap_time():
    a = _plate(1, _block((10, 20), (20, 30)), "continental")
    b = _plate(2, _block((20, 30), (20, 30)), "continental")
    world = _world(a, b)
    world.overlap_progress[(1, 2)] = merge_split.FORCED_MERGE_SUSTAINED_YEARS * 2

    forced = merge_split.pop_ready_forced_merge(world, can_merge=lambda x, y: merge_split._supports_merge(world, x, y))

    assert forced is None
    assert world.overlap_progress[(1, 2)] == merge_split.FORCED_MERGE_SUSTAINED_YEARS * 2
