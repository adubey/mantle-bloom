import numpy as np
import pytest
from app import geodesic
from app.world import generate_world


@pytest.mark.parametrize("frequency", [1, 2, 3, 8])
def test_subdivision_counts_match_the_geodesic_formula(frequency):
    vertices = geodesic._icosahedron_vertices()
    faces = geodesic._icosahedron_faces(vertices)
    assert vertices.shape == (12, 3)
    assert faces.shape == (20, 3)

    subdivided, triangles = geodesic._subdivide(vertices, faces, frequency)
    assert subdivided.shape == (geodesic.tile_count(frequency), 3)
    assert triangles.shape == (20 * frequency * frequency, 3)
    # Every subdivided vertex should still land exactly on the unit sphere.
    assert np.allclose(np.linalg.norm(subdivided, axis=1), 1.0, atol=1e-9)


@pytest.mark.parametrize("frequency", [1, 2, 3, 8])
def test_dual_mesh_has_exactly_twelve_pentagons_and_the_rest_hexagons(frequency):
    vertices = geodesic._icosahedron_vertices()
    faces = geodesic._icosahedron_faces(vertices)
    subdivided, triangles = geodesic._subdivide(vertices, faces, frequency)
    _, tiles = geodesic._build_dual(subdivided, triangles)

    assert len(tiles) == geodesic.tile_count(frequency)
    degrees = [len(tile.corner_vertex_ids) for tile in tiles]
    assert sum(1 for d in degrees if d == 5) == 12
    assert sum(1 for d in degrees if d == 6) == len(tiles) - 12
    assert all(d in (5, 6) for d in degrees)
    for tile in tiles:
        assert len(tile.neighbor_ids) == len(tile.corner_vertex_ids)


def test_neighbor_relation_is_symmetric():
    vertices = geodesic._icosahedron_vertices()
    faces = geodesic._icosahedron_faces(vertices)
    subdivided, triangles = geodesic._subdivide(vertices, faces, 4)
    _, tiles = geodesic._build_dual(subdivided, triangles)
    by_id = {tile.tile_id: tile for tile in tiles}

    for tile in tiles:
        assert len(set(tile.neighbor_ids)) == len(tile.neighbor_ids)  # no duplicate neighbors
        for neighbor_id in tile.neighbor_ids:
            assert tile.tile_id in by_id[neighbor_id].neighbor_ids


def test_every_tile_is_wound_outward_and_consistently():
    # Sum of signed tetrahedron volumes (origin, corner_k, corner_k+1) over every tile's fan
    # equals the sphere's own volume iff every tile is outward/CCW-wound (a winding error on
    # any single tile would flip that tile's own contribution negative, pulling the total
    # measurably away from 4/3 * pi for the unit sphere this mesh approximates).
    vertices = geodesic._icosahedron_vertices()
    faces = geodesic._icosahedron_faces(vertices)
    subdivided, triangles = geodesic._subdivide(vertices, faces, 8)
    dual_vertices, tiles = geodesic._build_dual(subdivided, triangles)

    signed_volume = 0.0
    for tile in tiles:
        center = tile.center_xyz
        n = len(tile.corner_vertex_ids)
        for k in range(n):
            p0 = dual_vertices[tile.corner_vertex_ids[k]]
            p1 = dual_vertices[tile.corner_vertex_ids[(k + 1) % n]]
            signed_volume += np.dot(center, np.cross(p0, p1)) / 6.0

    assert signed_volume == pytest.approx(4.0 / 3.0 * np.pi, rel=0.05)


def test_neighbor_ids_are_index_aligned_with_shared_corner_edges():
    vertices = geodesic._icosahedron_vertices()
    faces = geodesic._icosahedron_faces(vertices)
    subdivided, triangles = geodesic._subdivide(vertices, faces, 4)
    _, tiles = geodesic._build_dual(subdivided, triangles)
    by_id = {tile.tile_id: tile for tile in tiles}

    for tile in tiles:
        n = len(tile.corner_vertex_ids)
        for k in range(n):
            edge = {tile.corner_vertex_ids[k], tile.corner_vertex_ids[(k + 1) % n]}
            neighbor = by_id[tile.neighbor_ids[k]]
            neighbor_edges = {
                frozenset({neighbor.corner_vertex_ids[m], neighbor.corner_vertex_ids[(m + 1) % len(neighbor.corner_vertex_ids)]})
                for m in range(len(neighbor.corner_vertex_ids))
            }
            assert frozenset(edge) in neighbor_edges


def test_unknown_frequency_choice_is_not_silently_accepted():
    assert 7 not in geodesic.FREQUENCY_CHOICES  # main.py's /world/export_hexgrid rejects this


def test_export_hexgrid_against_a_real_world():
    world = generate_world(seed=3, num_plates=8)
    result = geodesic.export_hexgrid(world, frequency=8)

    assert result["num_tiles"] == geodesic.tile_count(8) == len(result["tiles"])
    assert sum(1 for t in result["tiles"] if t["is_pentagon"]) == 12

    by_id = {t["id"]: t for t in result["tiles"]}
    for tile in result["tiles"]:
        expected_corners = 5 if tile["is_pentagon"] else 6
        assert len(tile["corner_vertex_ids"]) == expected_corners
        assert len(tile["neighbor_ids"]) == expected_corners
        for corner_id in tile["corner_vertex_ids"]:
            assert 0 <= corner_id < len(result["vertices"])
        for neighbor_id in tile["neighbor_ids"]:
            assert tile["id"] in by_id[neighbor_id]["neighbor_ids"]
        assert isinstance(tile["biome"], str)
        assert isinstance(tile["is_ocean"], bool)


def test_export_hexgrid_on_an_empty_world_still_returns_every_tile():
    from app.world import World

    empty = World(seed=0)
    result = geodesic.export_hexgrid(empty, frequency=8)
    assert result["num_tiles"] == geodesic.tile_count(8)
    assert all(t["is_ocean"] for t in result["tiles"])
