"""Reusable contract tests for every authoritative PlateSurface implementation."""

import numpy as np

from app.elevation_lines import ElevationLine
from app.plates import PlateSurface, PlateWithLines
from app.surface_fields import SURFACE_FIELDS, RemapClass


def _surface() -> PlateSurface:
    lines = [
        ElevationLine(phi=-0.01, theta=np.array([-0.01, 0.0, 0.01]), elevation=np.array([1.0, 2.0, 3.0])),
        ElevationLine(phi=0.0, theta=np.array([-0.01, 0.0, 0.01]), elevation=np.array([4.0, 5.0, 6.0])),
        ElevationLine(phi=0.01, theta=np.array([-0.01, 0.0, 0.01]), elevation=np.array([7.0, 8.0, 9.0])),
    ]
    return PlateWithLines(plate_id=7, frame=np.eye(3), crust_type="continental", lines=lines)


def test_surface_bulk_view_has_one_stable_order_for_positions_and_fields():
    surface = _surface()
    nodes = surface.surface_nodes("elevation", "is_volcano")

    assert nodes.local_xyz.shape == nodes.world_xyz.shape == (9, 3)
    assert nodes.node_ids.shape == (9, 2)
    assert len({tuple(node_id) for node_id in nodes.node_ids}) == 9
    assert nodes.area_m2.shape == (9,)
    assert np.all(nodes.area_m2 > 0.0)
    assert not nodes.area_is_exact
    np.testing.assert_allclose(nodes.local_xyz, nodes.world_xyz)
    np.testing.assert_array_equal(nodes.fields["elevation"], np.arange(1.0, 10.0))
    assert nodes.fields["is_volcano"].dtype == bool


def test_surface_bulk_write_back_uses_the_same_order():
    surface = _surface()
    values = np.linspace(-100.0, 100.0, 9)

    surface.set_fields_on_plate(elevation=values)

    np.testing.assert_array_equal(surface.surface_nodes("elevation").fields["elevation"], values)


def test_surface_adjacency_is_valid_symmetric_csr():
    surface = _surface()
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


def test_surface_revisions_distinguish_fields_rigid_motion_and_topology_change():
    surface = _surface()
    ids0 = surface.surface_nodes().node_ids.copy()
    topology0, geometry0 = surface.topology_revision, surface.geometry_revision
    surface.rotate(np.eye(3))  # type: ignore[attr-defined]

    assert surface.topology_revision == topology0
    assert surface.geometry_revision == geometry0 + 1
    np.testing.assert_array_equal(surface.surface_nodes().node_ids, ids0)

    line_surface = surface
    assert isinstance(line_surface, PlateWithLines)
    line_surface.set_fields_on_plate(elevation=np.zeros(9))
    assert surface.topology_revision == topology0
    assert surface.geometry_revision == geometry0 + 1
    line_surface.set_lines(list(line_surface.lines))
    assert surface.topology_revision == topology0
    assert surface.geometry_revision == geometry0 + 1

    lines = list(line_surface.lines)
    lines[0] = ElevationLine(phi=-0.01, theta=np.array([-0.01, 0.0]), elevation=np.zeros(2))
    line_surface.set_lines(lines)
    assert surface.topology_revision == topology0 + 1
    assert surface.geometry_revision == geometry0 + 2


def test_surface_rejects_misaligned_or_unknown_field_write_before_mutating():
    surface = _surface()
    before = surface.surface_nodes("elevation").fields["elevation"].copy()

    for fields in ({"elevation": np.zeros(8)}, {"made_up": np.zeros(9)}):
        try:
            surface.set_fields_on_plate(**fields)
        except ValueError:
            pass
        else:
            raise AssertionError("invalid field write should fail")
        np.testing.assert_array_equal(surface.surface_nodes("elevation").fields["elevation"], before)


def test_surface_boundary_loops_and_containment_are_available_without_storage_access():
    surface = _surface()
    loops = surface.boundary_loops_world()

    assert len(loops) == 1
    assert loops[0].shape[1] == 3
    query = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    np.testing.assert_array_equal(surface.contains_batch(query), [True, False])


def test_surface_boundary_api_preserves_inner_loops():
    def line(phi, theta):
        theta = np.asarray(theta, dtype=float)
        return ElevationLine(phi=phi, theta=theta, elevation=np.zeros(len(theta)))

    surface: PlateSurface = PlateWithLines(
        plate_id=8,
        frame=np.eye(3),
        crust_type="continental",
        lines=[
            line(-0.02, np.linspace(-0.03, 0.03, 7)),
            line(0.0, [-0.03, -0.02]),
            line(0.0, [0.02, 0.03]),
            line(0.02, np.linspace(-0.03, 0.03, 7)),
        ],
    )

    loops = surface.boundary_loops_world()

    assert len(loops) == 2
    assert all(loop.shape[1] == 3 for loop in loops)


def test_surface_ids_distinguish_coincident_legacy_nodes_and_area_is_not_double_counted():
    duplicate_a = ElevationLine(phi=0.0, theta=np.array([0.0]), elevation=np.zeros(1))
    duplicate_b = ElevationLine(phi=0.0, theta=np.array([0.0]), elevation=np.zeros(1))
    neighbour = ElevationLine(phi=0.0, theta=np.array([0.01]), elevation=np.zeros(1))
    surface: PlateSurface = PlateWithLines(9, np.eye(3), "oceanic", [duplicate_a, duplicate_b, neighbour])

    nodes = surface.surface_nodes()

    assert len({tuple(node_id) for node_id in nodes.node_ids}) == 3
    np.testing.assert_allclose(nodes.area_m2[0], nodes.area_m2[1])
    np.testing.assert_allclose(nodes.area_m2[0] + nodes.area_m2[1], nodes.area_m2[2])


def test_every_persistent_field_has_remap_metadata():
    expected = {"elevation", *ElevationLine.OPTIONAL_FIELDS}
    surface = _surface()

    assert surface.field_metadata() is SURFACE_FIELDS
    assert set(surface.field_metadata()) == expected
    assert SURFACE_FIELDS["crustal_thickness_m"].remap_class is RemapClass.EXTENSIVE
    assert SURFACE_FIELDS["is_volcano"].remap_class is RemapClass.BOOLEAN_PROVENANCE
    assert SURFACE_FIELDS["node_created_years"].remap_class is RemapClass.WRITE_ONCE_HISTORY
    assert SURFACE_FIELDS["node_created_years"].sentinel == -1.0


def test_surface_node_id_lookup_is_storage_neutral_and_valid_only_in_its_revision():
    surface = _surface()
    node_id = surface.surface_nodes().node_ids[3]

    assert surface.node_index_for_id(node_id) == 3
    assert surface.node_index_for_id((int(node_id[0]), int(node_id[1]))) == 3
    assert surface.node_index_for_id((0, 0)) is None
