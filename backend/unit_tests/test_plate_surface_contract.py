"""Reusable contract tests for every authoritative PlateSurface implementation.

Each shared test runs unchanged against `PlateWithLines` and `PlateWithSparseQuadPatch`, built
as equivalent small patches around the local seed; representation-specific behaviour (legacy
line mutation, coincident line nodes) is tested separately below.
"""

import numpy as np
import pytest

from app.elevation_lines import ElevationLine
from app.plates import PlateSurface, PlateWithLines
from app.sparse_quad_patch import PlateWithSparseQuadPatch, pack_cell_keys
from app.surface_fields import SURFACE_FIELDS, RemapClass

# A 157-cell cube face has ~0.01 rad cells with cell 78 centred on the seed, matching the
# 0.01 rad line spacing below.
QUAD_N = 157
QUAD_CENTRE = 78


def _line(phi, theta, elevation=None):
    theta = np.asarray(theta, dtype=float)
    return ElevationLine(phi=phi, theta=theta, elevation=np.zeros(len(theta)) if elevation is None else elevation)


def _quad(plate_id, cells, **fields) -> PlateWithSparseQuadPatch:
    """`cells` as (i, j) offsets from the seed cell, listed in canonical (row-major) order."""
    di, dj = np.asarray(cells).T
    keys = pack_cell_keys(np.zeros(len(di), dtype=int), QUAD_CENTRE + di, QUAD_CENTRE + dj)
    return PlateWithSparseQuadPatch(plate_id, np.eye(3), "continental", QUAD_N, keys, fields=fields)


def _lines_square() -> PlateSurface:
    lines = [
        _line(-0.01, [-0.01, 0.0, 0.01], np.array([1.0, 2.0, 3.0])),
        _line(0.0, [-0.01, 0.0, 0.01], np.array([4.0, 5.0, 6.0])),
        _line(0.01, [-0.01, 0.0, 0.01], np.array([7.0, 8.0, 9.0])),
    ]
    return PlateWithLines(plate_id=7, frame=np.eye(3), crust_type="continental", lines=lines)


def _quad_square() -> PlateSurface:
    cells = [(di, dj) for dj in (-1, 0, 1) for di in (-1, 0, 1)]
    return _quad(7, cells, elevation=np.arange(1.0, 10.0))


def _lines_ring() -> PlateSurface:
    return PlateWithLines(
        plate_id=8,
        frame=np.eye(3),
        crust_type="continental",
        lines=[
            _line(-0.02, np.linspace(-0.03, 0.03, 7)),
            _line(0.0, [-0.03, -0.02]),
            _line(0.0, [0.02, 0.03]),
            _line(0.02, np.linspace(-0.03, 0.03, 7)),
        ],
    )


def _quad_ring() -> PlateSurface:
    cells = [(di, dj) for dj in (-1, 0, 1) for di in range(-3, 4) if dj != 0 or abs(di) >= 2]
    return _quad(8, cells)


SQUARES = {"lines": _lines_square, "quad": _quad_square}
RINGS = {"lines": _lines_ring, "quad": _quad_ring}


@pytest.fixture(params=sorted(SQUARES))
def surface(request) -> PlateSurface:
    return SQUARES[request.param]()


@pytest.fixture(params=sorted(RINGS))
def ring(request) -> PlateSurface:
    return RINGS[request.param]()


def test_surface_bulk_view_has_one_stable_order_for_positions_and_fields(surface):
    nodes = surface.surface_nodes("elevation", "is_volcano")

    assert nodes.local_xyz.shape == nodes.world_xyz.shape == (9, 3)
    assert nodes.node_ids.shape == (9, 2)
    assert len({tuple(node_id) for node_id in nodes.node_ids}) == 9
    assert nodes.area_m2.shape == (9,)
    assert np.all(nodes.area_m2 > 0.0)
    np.testing.assert_allclose(nodes.local_xyz, nodes.world_xyz)
    np.testing.assert_array_equal(nodes.fields["elevation"], np.arange(1.0, 10.0))
    assert nodes.fields["is_volcano"].dtype == bool


def test_surface_bulk_write_back_uses_the_same_order(surface):
    values = np.linspace(-100.0, 100.0, 9)

    surface.set_fields_on_plate(elevation=values)

    np.testing.assert_array_equal(surface.surface_nodes("elevation").fields["elevation"], values)


def test_surface_adjacency_is_valid_symmetric_csr(surface):
    graph = surface.adjacency()

    assert graph.offsets.shape == (10,)
    assert graph.offsets[0] == 0
    assert graph.offsets[-1] == len(graph.neighbours)
    edges = {
        (i, int(j))
        for i in range(9)
        for j in graph.neighbours[graph.offsets[i] : graph.offsets[i + 1]]
    }
    assert edges
    assert all(i != j for i, j in edges)
    assert all((j, i) in edges for i, j in edges)


def test_surface_revisions_distinguish_field_writes_from_rigid_motion(surface):
    ids0 = surface.surface_nodes().node_ids.copy()
    topology0, geometry0 = surface.topology_revision, surface.geometry_revision

    surface.rotate(np.eye(3))  # type: ignore[attr-defined]
    assert surface.topology_revision == topology0
    assert surface.geometry_revision == geometry0 + 1
    np.testing.assert_array_equal(surface.surface_nodes().node_ids, ids0)

    surface.set_fields_on_plate(elevation=np.zeros(9))
    assert surface.topology_revision == topology0
    assert surface.geometry_revision == geometry0 + 1


def test_surface_rigid_rotation_moves_world_views_but_not_local_ones(surface):
    before = surface.surface_nodes()
    quarter_turn = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])

    surface.rotate(quarter_turn)  # type: ignore[attr-defined]
    after = surface.surface_nodes()

    np.testing.assert_allclose(after.local_xyz, before.local_xyz)
    np.testing.assert_allclose(after.world_xyz, before.world_xyz @ quarter_turn.T)
    np.testing.assert_allclose(after.area_m2, before.area_m2)
    np.testing.assert_array_equal(surface.contains_batch(np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 0.0]])), [True, False])


def test_surface_rejects_misaligned_or_unknown_field_write_before_mutating(surface):
    before = surface.surface_nodes("elevation").fields["elevation"].copy()

    for fields in ({"elevation": np.zeros(8)}, {"made_up": np.zeros(9)}):
        with pytest.raises(ValueError):
            surface.set_fields_on_plate(**fields)
        np.testing.assert_array_equal(surface.surface_nodes("elevation").fields["elevation"], before)


def test_surface_boundary_loops_and_containment_are_available_without_storage_access(surface):
    loops = surface.boundary_loops_world()

    assert len(loops) == 1
    assert loops[0].shape[1] == 3
    query = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    np.testing.assert_array_equal(surface.contains_batch(query), [True, False])


def test_surface_boundary_api_preserves_inner_loops(ring):
    loops = ring.boundary_loops_world()

    assert len(loops) == 2
    assert all(loop.shape[1] == 3 for loop in loops)
    # The hole is outside the territory, the ring itself inside.
    np.testing.assert_array_equal(ring.contains_batch(np.array([[1.0, 0.0, 0.0], [1.0, 0.025, 0.0]])), [False, True])


def test_surface_node_id_lookup_is_storage_neutral_and_valid_only_in_its_revision(surface):
    node_id = surface.surface_nodes().node_ids[3]

    assert surface.node_index_for_id(node_id) == 3
    assert surface.node_index_for_id((int(node_id[0]), int(node_id[1]))) == 3
    assert surface.node_index_for_id((0, 0)) is None


def test_every_persistent_field_has_remap_metadata(surface):
    expected = {"elevation", *ElevationLine.OPTIONAL_FIELDS}

    assert surface.field_metadata() is SURFACE_FIELDS
    assert set(surface.field_metadata()) == expected
    assert SURFACE_FIELDS["crustal_thickness_m"].remap_class is RemapClass.EXTENSIVE
    assert SURFACE_FIELDS["is_volcano"].remap_class is RemapClass.BOOLEAN_PROVENANCE
    assert SURFACE_FIELDS["node_created_years"].remap_class is RemapClass.WRITE_ONCE_HISTORY
    assert SURFACE_FIELDS["node_created_years"].sentinel == -1.0


def test_every_persistent_field_reads_back_with_its_registry_dtype_and_default(surface):
    for name, spec in SURFACE_FIELDS.items():
        if name == "elevation":
            continue
        values = surface.surface_nodes(name).fields[name]
        assert values.shape == (9,)
        assert values.dtype.kind == spec.dtype.kind, name
        np.testing.assert_array_equal(values, np.full(9, spec.default, dtype=spec.dtype), err_msg=name)


# --- Representation-specific behaviour ----------------------------------------------------


def test_line_surface_area_is_a_compatibility_estimate_and_quad_area_is_exact():
    assert not _lines_square().surface_nodes().area_is_exact
    assert _quad_square().surface_nodes().area_is_exact


def test_line_surface_topology_revision_tracks_line_geometry_changes():
    surface = _lines_square()
    assert isinstance(surface, PlateWithLines)
    topology0, geometry0 = surface.topology_revision, surface.geometry_revision

    surface.set_lines(list(surface.lines))
    assert surface.topology_revision == topology0
    assert surface.geometry_revision == geometry0

    lines = list(surface.lines)
    lines[0] = _line(-0.01, [-0.01, 0.0])
    surface.set_lines(lines)
    assert surface.topology_revision == topology0 + 1
    assert surface.geometry_revision == geometry0 + 1


def test_surface_ids_distinguish_coincident_legacy_nodes_and_area_is_not_double_counted():
    duplicate_a = _line(0.0, [0.0])
    duplicate_b = _line(0.0, [0.0])
    neighbour = _line(0.0, [0.01])
    surface: PlateSurface = PlateWithLines(9, np.eye(3), "oceanic", [duplicate_a, duplicate_b, neighbour])

    nodes = surface.surface_nodes()

    assert len({tuple(node_id) for node_id in nodes.node_ids}) == 3
    np.testing.assert_allclose(nodes.area_m2[0], nodes.area_m2[1])
    np.testing.assert_allclose(nodes.area_m2[0] + nodes.area_m2[1], nodes.area_m2[2])
