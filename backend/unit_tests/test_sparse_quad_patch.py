"""PlateWithSparseQuadPatch-specific invariants: the cube-sphere lattice, derived views,
exact containment, persistence, and generation from the shared tiling (issue #228 Phase 2).
Behaviour shared with PlateWithLines lives in test_plate_surface_contract.py."""

import pickle

import numpy as np
import pytest
from scipy.spatial import cKDTree

from app import geometry, persistence
from app.lithosphere_plate import generate_plates
from app.plates import gather_node_positions
from app.sparse_quad_patch import (
    PLANET_RADIUS_M,
    QUAD_SURFACE_FORMAT_VERSION,
    PlateWithSparseQuadPatch,
    cell_areas_sr,
    cells_per_face_edge,
    lattice_points,
    locate_cells,
    child_cell_keys,
    pack_cell_keys,
    parent_cell_keys,
    unpack_cell_keys,
)
from app.world import generate_world, step_world

N = 12


def _all_keys(n: int = N) -> np.ndarray:
    jj, ii = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    return np.concatenate([pack_cell_keys(np.full(n * n, f), ii.ravel(), jj.ravel()) for f in range(6)])


def _plate(keys, frame=None, n: int = N, **fields) -> PlateWithSparseQuadPatch:
    return PlateWithSparseQuadPatch(1, np.eye(3) if frame is None else frame, "oceanic", n, keys, fields=fields)


def _random_unit_vectors(count: int, seed: int = 0) -> np.ndarray:
    points = np.random.default_rng(seed).normal(size=(count, 3))
    return points / np.linalg.norm(points, axis=1, keepdims=True)


def _cap_plate(min_dot: float, n: int = N) -> PlateWithSparseQuadPatch:
    """Every cell whose centre is within the cap `local . seed > min_dot`, straddling the
    face-0 seams whenever the cap is wider than ~45 degrees."""
    whole = _plate(_all_keys(n), n=n)
    return _plate(whole.cell_keys[whole.surface_nodes().local_xyz[:, 0] > min_dot], n=n)


def test_key_packing_round_trips():
    face, i, j = np.array([0, 3, 5]), np.array([0, 7, 11]), np.array([11, 2, 0])
    f, level, ii, jj = unpack_cell_keys(pack_cell_keys(face, i, j))
    np.testing.assert_array_equal(f, face)
    np.testing.assert_array_equal(level, 0)
    np.testing.assert_array_equal(ii, i)
    np.testing.assert_array_equal(jj, j)


def test_locate_inverts_lattice_points_on_every_face():
    keys = _all_keys()
    face, _, i, j = unpack_cell_keys(keys)
    centres = lattice_points(face, i + 0.5, j + 0.5, N)
    f, ii, jj = locate_cells(centres, N)
    np.testing.assert_array_equal(pack_cell_keys(f, ii, jj), keys)


def test_full_lattice_tiles_the_sphere_with_exact_areas():
    plate = _plate(_all_keys())
    areas = plate.node_areas_m2()

    np.testing.assert_allclose(areas.sum(), 4.0 * np.pi * PLANET_RADIUS_M**2, rtol=1e-12)
    # Equiangular cells stay within a modest size range -- no polar pinching.
    assert areas.max() / areas.min() < 1.5
    assert plate.boundary_loops_world() == ()
    assert np.all(plate.contains_batch(_random_unit_vectors(2000)))


def test_adjacency_is_four_connected_and_symmetric_across_cube_seams():
    plate = _plate(_all_keys())
    graph = plate.adjacency()

    np.testing.assert_array_equal(np.diff(graph.offsets), 4)
    edges = {(a, int(b)) for a in range(plate.node_count()) for b in graph.neighbours[graph.offsets[a] : graph.offsets[a + 1]]}
    assert all((b, a) in edges for a, b in edges)
    # Neighbours are genuinely adjacent: centre separation is about one cell everywhere,
    # seams included.
    centres = plate.surface_nodes().world_xyz
    a, b = np.array(sorted(edges)).T
    separation = np.arccos(np.clip(np.einsum("ij,ij->i", centres[a], centres[b]), -1.0, 1.0))
    cell = (np.pi / 2) / N
    assert separation.min() > 0.6 * cell
    assert separation.max() < 1.2 * cell


def test_neighbour_cells_across_a_seam_share_their_edge():
    plate = _plate(_all_keys())
    face, _, i, _ = unpack_cell_keys(plate.cell_keys)
    graph = plate.adjacency()
    right_edge = np.flatnonzero((face == 0) & (i == N - 1))
    for index in right_edge:
        neighbours = graph.neighbours[graph.offsets[index] : graph.offsets[index + 1]]
        across = [n for n in neighbours if unpack_cell_keys(plate.cell_keys[n])[0] == 1]
        assert len(across) == 1


@pytest.mark.parametrize("min_dot", [0.9, 0.3, -0.4])
def test_containment_is_exactly_the_active_cells(min_dot):
    plate = _cap_plate(min_dot)
    points = _random_unit_vectors(5000, seed=1)
    face, i, j = locate_cells(points, N)
    expected = np.isin(pack_cell_keys(face, i, j), plate.cell_keys)

    np.testing.assert_array_equal(plate.contains_batch(points), expected)


@pytest.mark.parametrize("min_dot", [0.9, 0.3])
def test_boundary_loops_agree_with_cell_containment(min_dot):
    # Up to a hemisphere; the winding-number test itself only supports that range.
    plate = _cap_plate(min_dot)
    points = _random_unit_vectors(5000, seed=2)

    np.testing.assert_array_equal(
        geometry.points_in_spherical_polygon(points, plate.get_bounding_polygon()), plate.contains_batch(points)
    )


def _signed_area(loop: np.ndarray) -> float:
    """Positive for a counter-clockwise loop seen from outside (small loops only)."""
    centre = loop.mean(axis=0)
    centre /= np.linalg.norm(centre)
    return float(sum(np.dot(np.cross(a - centre, b - centre), centre) for a, b in zip(loop, np.roll(loop, -1, axis=0))))


def test_outer_loops_run_counter_clockwise_and_holes_clockwise():
    keys = [pack_cell_keys(0, i, j) for j in range(3, 8) for i in range(3, 8) if (i, j) != (5, 5)]
    plate = _plate(np.array(keys))
    loops = sorted(plate.boundary_loops_world(), key=len)

    assert [len(loop) for loop in loops] == [4, 20]
    assert _signed_area(loops[0]) < 0 < _signed_area(loops[1])


def test_diagonal_contact_yields_two_simple_loops():
    plate = _plate(np.array([pack_cell_keys(0, 4, 4), pack_cell_keys(0, 5, 5)]))
    loops = plate.boundary_loops_world()

    assert sorted(len(loop) for loop in loops) == [4, 4]


def test_row_and_column_intervals_reproduce_the_active_cells():
    plate = _cap_plate(0.3)
    face, _, i, j = unpack_cell_keys(plate.cell_keys)
    for intervals, line, position in ((plate.row_intervals(), j, i), (plate.column_intervals(), i, j)):
        assert np.all(intervals.end >= intervals.start)
        run = intervals.node_run
        np.testing.assert_array_equal(intervals.face[run], face)
        np.testing.assert_array_equal(intervals.line[run], line)
        assert np.all((position >= intervals.start[run]) & (position <= intervals.end[run]))
        assert (intervals.end - intervals.start + 1).sum() == plate.node_count()


def test_rows_split_by_a_notch_carry_several_intervals():
    keys = [pack_cell_keys(0, i, 4) for i in (1, 2, 3, 6, 7)]
    rows = _plate(np.array(keys)).row_intervals()

    np.testing.assert_array_equal(rows.start, [1, 6])
    np.testing.assert_array_equal(rows.end, [3, 7])


def test_boundary_fraction_is_measured_along_each_row_run():
    keys = [pack_cell_keys(0, i, 4) for i in range(2, 7)] + [pack_cell_keys(0, 9, 4)]
    fractions = [fraction for _, _, fraction in _plate(np.array(keys)).map_world_points_on_plate()]

    np.testing.assert_allclose(fractions, [0.0, 0.25, 0.5, 0.75, 1.0, 0.5])


def test_point_views_write_through_to_bulk_fields():
    plate = _plate(np.array([pack_cell_keys(0, 4, 4), pack_cell_keys(0, 5, 4)]))
    points = list(plate)
    points[1].set_elevation(123.0)
    points[0].set_is_volcano(True)

    np.testing.assert_array_equal(plate.collect("elevation"), [0.0, 123.0])
    np.testing.assert_array_equal(plate.collect("is_volcano"), [True, False])
    assert points[1].get_theta() > points[0].get_theta()


def test_rotation_keeps_topology_caches_and_invalidates_world_ones():
    plate = _cap_plate(0.9)
    adjacency = plate.adjacency()
    loops_before = plate.boundary_loops_world()
    rotation = geometry.plate_frame_from_seed(np.array([0.0, 0.6, 0.8]))

    plate.rotate(rotation)

    assert plate.adjacency().neighbours is not None
    np.testing.assert_array_equal(plate.adjacency().neighbours, adjacency.neighbours)
    for before, after in zip(loops_before, plate.boundary_loops_world()):
        np.testing.assert_allclose(after, before @ rotation.T, atol=1e-12)
    centre = rotation @ np.array([1.0, 0.0, 0.0])
    assert plate.contains_batch(centre[None, :])[0]
    assert not plate.contains_batch(np.array([[1.0, 0.0, 0.0]]))[0]


def test_invalid_cells_are_rejected():
    with pytest.raises(ValueError):
        _plate(np.array([pack_cell_keys(0, N, 0)]))
    with pytest.raises(ValueError):
        _plate(np.array([pack_cell_keys(0, 1, 1), pack_cell_keys(0, 1, 1)]))


def test_pickle_round_trip_keeps_authoritative_state_and_drops_caches():
    plate = _cap_plate(0.3)
    plate.set_fields_on_plate(elevation=np.arange(plate.node_count(), dtype=float))
    plate.rotate(geometry.plate_frame_from_seed(np.array([0.0, 0.0, 1.0])))
    plate.adjacency()
    plate.get_bounding_polygon()

    state = plate.__getstate__()
    assert state["_surface_format_version"] == QUAD_SURFACE_FORMAT_VERSION
    assert all(state.get(name) is None for name in PlateWithSparseQuadPatch._TOPOLOGY_CACHES)
    assert all(state.get(name) is None for name in PlateWithSparseQuadPatch._GEOMETRY_CACHES)

    loaded = pickle.loads(pickle.dumps(plate))
    np.testing.assert_array_equal(loaded.cell_keys, plate.cell_keys)
    np.testing.assert_array_equal(loaded.collect("elevation"), plate.collect("elevation"))
    np.testing.assert_array_equal(loaded.surface_nodes().node_ids, plate.surface_nodes().node_ids)
    np.testing.assert_allclose(loaded.get_bounding_polygon(), plate.get_bounding_polygon())
    assert (loaded.topology_revision, loaded.geometry_revision) == (plate.topology_revision, plate.geometry_revision)


def test_pickle_from_a_newer_surface_format_is_rejected():
    state = _cap_plate(0.9).__getstate__()
    state["_surface_format_version"] = QUAD_SURFACE_FORMAT_VERSION + 1
    fresh = PlateWithSparseQuadPatch.__new__(PlateWithSparseQuadPatch)

    with pytest.raises(ValueError, match="format version"):
        fresh.__setstate__(state)


def test_lattice_resolution_tracks_line_spacing():
    assert cells_per_face_edge(np.sqrt(4 * np.pi / 6) / 40) == 40


def test_generated_quad_plates_tile_the_sphere_from_the_shared_tiling():
    quad = generate_plates(3, num_plates=6, node_density=1.0, surface="quad")
    lines = generate_plates(3, num_plates=6, node_density=1.0)

    assert [p.crust_type for p in quad] == [p.crust_type for p in lines]
    np.testing.assert_allclose([p.frame for p in quad], [p.frame for p in lines])
    total_area = sum(p.node_areas_m2().sum() for p in quad)
    np.testing.assert_allclose(total_area, 4 * np.pi * PLANET_RADIUS_M**2, rtol=0.01)

    points = _random_unit_vectors(20000, seed=4)
    coverage = sum(p.contains_batch(points).astype(int) for p in quad)
    assert np.mean(coverage == 1) > 0.97
    # Same ownership test, so each point belongs to the same plate in both worlds, up to
    # one-cell boundary jitter. Line ownership is read as "plate of the nearest line node"
    # rather than PlateWithLines.contains_batch, whose polygon fallback is unreliable for a
    # plate covering about a hemisphere (this seed has one).
    line_points, line_plates = gather_node_positions(lines)
    counts = [p.node_count() for p in line_plates]
    node_owner = np.repeat([p.plate_id for p in line_plates], counts)
    line_owner = node_owner[cKDTree(line_points).query(points)[1]]
    quad_owner = np.array([p.plate_id for p in quad])[np.argmax([p.contains_batch(points) for p in quad], axis=0)]
    covered = coverage == 1
    assert np.mean(line_owner[covered] == quad_owner[covered]) > 0.98
    # Every field the isostasy sync reads was populated.
    for plate in quad:
        assert np.all(plate.collect("crustal_thickness_m") > 0.0)
        assert np.all(plate.collect("mantle_lithosphere_thickness_m") > 0.0)


def test_quad_world_round_trips_through_the_versioned_save_format():
    world = generate_world(seed=5, num_plates=5, surface="quad")
    data = persistence.save_world_bytes(world)
    envelope = pickle.loads(data)
    assert envelope["version"] == persistence.SAVE_FORMAT_VERSION

    loaded = persistence.load_world_bytes(data)
    assert all(isinstance(p, PlateWithSparseQuadPatch) for p in loaded.plates)
    for original, restored in zip(world.plates, loaded.plates):
        np.testing.assert_array_equal(restored.cell_keys, original.cell_keys)
        for name in ("elevation", "crustal_thickness_m", "soil_depth"):
            np.testing.assert_array_equal(restored.collect(name), original.collect(name))


def test_quad_world_refuses_plate_movement_without_mutating():
    world = generate_world(seed=5, num_plates=5, surface="quad")

    with pytest.raises(NotImplementedError):
        step_world(world, 1_000_000)
    assert world.steps_taken == 0
    assert world.elapsed_years == 0.0


def test_unknown_surface_representation_is_rejected():
    with pytest.raises(ValueError):
        generate_plates(3, num_plates=4, surface="hexes")


def test_refine_and_coarsen_round_trip_geometry_identity_and_revisions():
    root = int(pack_cell_keys(0, 4, 4))
    neighbour = int(pack_cell_keys(0, 5, 4))
    plate = _plate(np.array([root, neighbour]), elevation=np.array([12.0, 20.0]))
    area_before = plate.node_areas_m2().sum()

    lineage = plate.refine_cells(np.array([root]))

    assert lineage[root] == tuple(map(int, child_cell_keys(np.array([root]))[0]))
    np.testing.assert_array_equal(parent_cell_keys(np.array(lineage[root])), root)
    assert plate.node_count() == 5
    assert (plate.topology_revision, plate.geometry_revision) == (1, 1)
    np.testing.assert_allclose(plate.node_areas_m2().sum(), area_before, rtol=1e-12)
    for child in lineage[root]:
        assert plate.node_index_for_id((child, 0)) is not None

    reverse = plate.coarsen_cells(np.array([root]))

    assert all(reverse[child] == root for child in lineage[root])
    np.testing.assert_array_equal(plate.cell_keys, [root, neighbour])
    np.testing.assert_allclose(plate.collect("elevation"), [12.0, 20.0])
    assert (plate.topology_revision, plate.geometry_revision) == (2, 2)


def test_mixed_level_adjacency_is_symmetric_and_boundary_has_no_crack():
    roots = np.array([pack_cell_keys(0, 4, 4), pack_cell_keys(0, 5, 4)])
    plate = _plate(roots)
    plate.refine_cells(roots[:1])
    graph = plate.adjacency()
    edges = {(a, int(b)) for a in range(plate.node_count()) for b in graph.neighbours[graph.offsets[a] : graph.offsets[a + 1]]}

    assert all((b, a) in edges for a, b in edges)
    # The union is still a rectangle: refinement adds hanging boundary vertices but no
    # internal loop along the coarse/fine interface.
    assert len(plate.boundary_loops_world()) == 1


def test_boundary_loop_includes_exposed_half_of_coarse_side():
    # The fine cell covers only the lower half of the coarse cell's right side. The upper
    # half must remain in the outline even though one of that side's probes found a neighbour.
    coarse = pack_cell_keys(0, 0, 0)
    fine = pack_cell_keys(0, 2, 0, level=1)
    plate = _plate(np.array([coarse, fine]), n=8)

    loops = plate.boundary_loops_world()

    assert len(loops) == 1
    assert len(loops[0]) == 7
    centres = plate.surface_nodes().world_xyz
    assert np.all(geometry.points_in_spherical_polygon(centres, loops[0]))


def test_constructor_and_loader_reject_unbalanced_leaf_layouts():
    coarse = pack_cell_keys(0, 0, 0)
    two_levels_finer = pack_cell_keys(0, 4, 0, level=2)
    invalid_keys = np.array([coarse, two_levels_finer])

    with pytest.raises(ValueError, match="2:1 balanced"):
        _plate(invalid_keys, n=8)

    state = _plate(np.array([coarse]), n=8).__getstate__()
    state["_keys"] = invalid_keys
    fresh = PlateWithSparseQuadPatch.__new__(PlateWithSparseQuadPatch)
    with pytest.raises(ValueError, match="2:1 balanced"):
        fresh.__setstate__(state)


def test_constructor_rejects_active_ancestor_and_descendant():
    root = pack_cell_keys(0, 0, 0)
    child = child_cell_keys(np.array([root]))[0, 0]

    with pytest.raises(ValueError, match="active ancestor"):
        _plate(np.array([root, child]), n=8)


def test_coarsen_applies_every_field_policy_and_conserves_extensive_integrals():
    root = int(pack_cell_keys(0, 4, 4))
    children = child_cell_keys(np.array([root]))[0]
    plate = _plate(
        children,
        crustal_thickness_m=np.array([10.0, 20.0, 30.0, 40.0]),
        crust_type_code=np.array([2, 1, 2, 1], dtype=np.int8),
        is_volcano=np.array([False, False, True, False]),
        volcano_active_years_remaining=np.array([1.0, 7.0, 3.0, 2.0]),
        node_created_years=np.array([-1.0, 9.0, 5.0, -1.0]),
    )
    before = np.sum(plate.collect("crustal_thickness_m") * plate.node_areas_m2())

    plate.coarsen_cells(np.array([root]))

    after = np.sum(plate.collect("crustal_thickness_m") * plate.node_areas_m2())
    np.testing.assert_allclose(after, before, rtol=1e-12)
    assert plate.collect("crust_type_code")[0] == 0  # winning plate type stays inherited
    assert plate.collect("is_volcano")[0]
    assert plate.collect("volcano_active_years_remaining")[0] == 7.0
    assert plate.collect("node_created_years")[0] == 5.0


def test_refinement_balances_a_coarser_neighbour():
    left_root = int(pack_cell_keys(0, 4, 4))
    right_root = int(pack_cell_keys(0, 5, 4))
    plate = _plate(np.array([left_root, right_root]))
    first_child = int(plate.refine_cells(np.array([left_root]))[left_root][1])

    lineage = plate.refine_cells(np.array([first_child]))

    assert right_root in lineage
    levels = unpack_cell_keys(plate.cell_keys)[1]
    graph = plate.adjacency()
    for cell in range(plate.node_count()):
        neighbours = graph.neighbours[graph.offsets[cell] : graph.offsets[cell + 1]]
        assert np.all(np.abs(levels[neighbours] - levels[cell]) <= 1)


def test_coarsening_preserves_inherited_crust_type_and_uses_it_for_relief():
    root = int(pack_cell_keys(0, 4, 4))
    children = child_cell_keys(np.array([root]))[0]
    plate = PlateWithSparseQuadPatch(
        1,
        np.eye(3),
        "continental",
        N,
        children,
        fields={
            "elevation": np.array([100.0, 200.0, 300.0, 400.0]),
            "crustal_thickness_m": np.full(4, 35_000.0),
            "mantle_lithosphere_thickness_m": np.full(4, 100_000.0),
            "crust_type_code": np.array([0, 2, 0, 2], dtype=np.int8),
        },
    )

    plate.coarsen_cells(np.array([root]))

    assert plate.collect("crust_type_code")[0] == 0
    # Equal columns mean preserving the area-weighted residual is exactly an area-weighted
    # elevation mean; this exercises the isostatic path while code 0 remains implicit.
    expected = np.average([100.0, 200.0, 300.0, 400.0], weights=cell_areas_sr(*unpack_cell_keys(children)[2:], N * 2))
    np.testing.assert_allclose(plate.collect("elevation")[0], expected)


def test_default_crust_code_does_not_skip_isostatic_elevation_transfer():
    from app import lithosphere

    root = int(pack_cell_keys(0, 4, 4))
    children = child_cell_keys(np.array([root]))[0]
    hc = np.array([30_000.0, 35_000.0, 40_000.0, 45_000.0])
    hm = np.full(4, 100_000.0)
    equilibrium = lithosphere.isostatic_elevation(hc, hm, np.full(4, lithosphere.RHO_CONTINENTAL_CRUST))
    plate = PlateWithSparseQuadPatch(
        1,
        np.eye(3),
        "continental",
        N,
        children,
        fields={"elevation": equilibrium + 125.0, "crustal_thickness_m": hc, "mantle_lithosphere_thickness_m": hm},
    )

    plate.coarsen_cells(np.array([root]))

    expected_equilibrium = lithosphere.isostatic_elevation(
        plate.collect("crustal_thickness_m"),
        plate.collect("mantle_lithosphere_thickness_m"),
        np.array([lithosphere.RHO_CONTINENTAL_CRUST]),
    )
    np.testing.assert_allclose(plate.collect("elevation"), expected_equilibrium + 125.0)


def test_deep_refinement_boundary_corner_identity_does_not_overflow():
    root = int(pack_cell_keys(0, 4, 4))
    plate = _plate(np.array([root]))
    leaf = root
    for _ in range(18):
        leaf = plate.refine_cells(np.array([leaf]))[leaf][0]

    loops = plate.boundary_loops_world()

    assert len(loops) == 1
    assert np.all(np.isfinite(loops[0]))
