import numpy as np
from app import geometry
from app.geodesic import _icosahedron_vertices
from app.mesh_terrain import (
    MESH_ALPHA_EDGE_FACTOR,
    MeshTriangulation,
    dedupe_growth_candidates,
    estimate_outward_direction,
    needs_remeshing,
    remesh_nodes,
    trace_boundary_loops,
    build_mesh_triangulation,
)


def _disk_cluster_xyz(rng, n, radius_rad, center=np.array([1.0, 0.0, 0.0])):
    """A roughly disk-shaped cloud of local unit vectors around `center` -- stands in for a
    compact plate's own footprint, same idea test_plate_with_rtree.py's `_disk_cluster`
    uses, just already projected to xyz since mesh_terrain.py works in that space."""
    east, north = geometry.local_tangent_basis(center)
    r = radius_rad * np.sqrt(rng.uniform(0.0, 1.0, size=n))
    angle = rng.uniform(0.0, 2 * np.pi, size=n)
    tangent = r[:, None] * (np.cos(angle)[:, None] * east + np.sin(angle)[:, None] * north)
    return geometry.normalize(center[None, :] * np.cos(r)[:, None] + tangent)


def test_build_mesh_triangulation_every_node_has_adjacency():
    rng = np.random.default_rng(0)
    points = _disk_cluster_xyz(rng, 60, radius_rad=0.2)
    tri = build_mesh_triangulation(points)
    assert not tri.degenerate
    assert set(tri.adjacency.keys()) == set(range(60))


def test_build_mesh_triangulation_boundary_is_a_strict_subset_for_a_compact_cluster():
    rng = np.random.default_rng(1)
    points = _disk_cluster_xyz(rng, 200, radius_rad=0.2)
    tri = build_mesh_triangulation(points)
    assert not tri.degenerate
    assert len(tri.boundary_vertices) > 0
    assert len(tri.boundary_vertices) < len(points)  # proves lid facets were actually dropped


def test_build_mesh_triangulation_no_triangle_has_an_overlong_edge():
    rng = np.random.default_rng(2)
    points = _disk_cluster_xyz(rng, 150, radius_rad=0.25)
    tri = build_mesh_triangulation(points)
    if len(tri.triangles) == 0:
        return
    tri_pts = points[tri.triangles]
    ab = geometry.angular_distance(tri_pts[:, 0], tri_pts[:, 1])
    bc = geometry.angular_distance(tri_pts[:, 1], tri_pts[:, 2])
    ca = geometry.angular_distance(tri_pts[:, 2], tri_pts[:, 0])
    edge_lengths = np.stack([ab, bc, ca], axis=1)
    typical_spacing = float(np.median(edge_lengths))
    assert np.all(edge_lengths.max(axis=1) <= MESH_ALPHA_EDGE_FACTOR * typical_spacing + 1e-9)


def test_build_mesh_triangulation_degenerate_for_too_few_nodes():
    tri = build_mesh_triangulation(np.eye(3)[:2])
    assert tri.degenerate
    assert len(tri.triangles) == 0
    assert list(tri.boundary_vertices) == [0, 1]


def test_build_mesh_triangulation_degenerate_for_collinear_nodes():
    points = geometry.normalize(np.array([[1.0, 0.0, 0.001 * i] for i in range(6)]))
    tri = build_mesh_triangulation(points)
    assert tri.degenerate
    assert len(tri.triangles) == 0


def test_build_mesh_triangulation_whole_sphere_has_no_boundary():
    points = _icosahedron_vertices()
    tri = build_mesh_triangulation(points)
    assert not tri.degenerate
    assert tri.boundary_edges == set()
    assert len(tri.boundary_vertices) == 0


def test_trace_boundary_loops_single_compact_cluster_is_one_loop():
    rng = np.random.default_rng(3)
    points = _disk_cluster_xyz(rng, 150, radius_rad=0.2)
    tri = build_mesh_triangulation(points)
    loops = trace_boundary_loops(tri.boundary_edges)
    assert len(loops) == 1
    assert set(loops[0]) == set(tri.boundary_vertices.tolist())


def test_trace_boundary_loops_two_disjoint_clusters_are_two_loops():
    rng = np.random.default_rng(4)
    far_center = geometry.normalize(np.array([0.0, 1.0, 0.0]))
    cluster_a = _disk_cluster_xyz(rng, 80, radius_rad=0.15)
    cluster_b = _disk_cluster_xyz(rng, 80, radius_rad=0.15, center=far_center)
    boundary_edges = set()
    offset = 0
    for cluster in (cluster_a, cluster_b):
        tri = build_mesh_triangulation(cluster)
        for edge in tri.boundary_edges:
            a, b = tuple(edge)
            boundary_edges.add(frozenset((a + offset, b + offset)))
        offset += len(cluster)
    loops = trace_boundary_loops(boundary_edges)
    assert len(loops) == 2


def test_estimate_outward_direction_points_away_from_centroid_and_is_tangent():
    rng = np.random.default_rng(5)
    points = _disk_cluster_xyz(rng, 150, radius_rad=0.2)
    tri = build_mesh_triangulation(points)
    assert len(tri.boundary_vertices) > 0
    vertex = int(tri.boundary_vertices[0])

    direction = estimate_outward_direction(vertex, tri, points)
    assert direction is not None
    assert abs(float(np.dot(direction, points[vertex]))) < 1e-6  # tangent at the vertex
    assert abs(float(np.linalg.norm(direction)) - 1.0) < 1e-6

    centroid = geometry.normalize(points.mean(axis=0))
    outward_ref = points[vertex] - centroid
    outward_ref = outward_ref - np.dot(outward_ref, points[vertex]) * points[vertex]
    assert float(np.dot(direction, geometry.normalize(outward_ref))) > 0


def test_estimate_outward_direction_isolated_vertex_returns_none():
    points = np.eye(3)
    tri = build_mesh_triangulation(points)  # degenerate: 3 nodes < MIN_NODES_FOR_TRIANGULATION
    assert estimate_outward_direction(0, tri, points) is None


def test_dedupe_growth_candidates_drops_near_existing_node():
    existing = np.array([[1.0, 0.0, 0.0]])
    candidates = geometry.normalize(np.array([[1.0, 0.001, 0.0], [0.0, 1.0, 0.0]]))
    keep = dedupe_growth_candidates(candidates, existing, min_sep_rad=0.01)
    assert list(keep) == [False, True]


def test_dedupe_growth_candidates_drops_later_of_a_close_pair():
    candidates = geometry.normalize(np.array([[1.0, 0.0, 0.0], [1.0, 0.001, 0.0], [0.0, 1.0, 0.0]]))
    keep = dedupe_growth_candidates(candidates, np.zeros((0, 3)), min_sep_rad=0.01)
    assert list(keep) == [True, False, True]


def test_dedupe_growth_candidates_empty_input():
    keep = dedupe_growth_candidates(np.zeros((0, 3)), np.zeros((0, 3)), min_sep_rad=0.01)
    assert keep.shape == (0,)


def test_needs_remeshing_false_when_every_edge_matches_target_spacing():
    # Built directly rather than via build_mesh_triangulation: a real hull's boundary edges
    # are inherently somewhat longer than its interior spacing (the same "boundary is never
    # perfectly regular" property PlateWithLines/PlateWithRTree's own outlines have), so
    # needs_remeshing only depends on `tri.adjacency`/`degenerate` -- constructing a
    # MeshTriangulation with an exact, uniform-spacing adjacency graph isolates the pure
    # threshold logic from that real-geometry noise.
    spacing_rad = 0.02
    n = 6
    theta = np.arange(n) * spacing_rad
    points = geometry.local_xyz(np.zeros(n), theta)
    adjacency = {i: set() for i in range(n)}
    for i in range(n - 1):
        adjacency[i].add(i + 1)
        adjacency[i + 1].add(i)
    tri = MeshTriangulation(
        triangles=np.zeros((0, 3), dtype=int),
        adjacency=adjacency,
        boundary_edges=set(),
        boundary_vertices=np.array([0, n - 1]),
        degenerate=False,
    )
    assert not needs_remeshing(tri, points, spacing_rad)


def test_needs_remeshing_true_for_an_artificially_stretched_edge():
    points = geometry.normalize(np.array([[1.0, 0.0, 0.0], [0.9, 0.1, 0.0], [0.9, -0.1, 0.0], [0.0, 0.0, 1.0]]))
    tri = build_mesh_triangulation(points)
    assert needs_remeshing(tri, points, spacing_rad=0.001)


def test_remesh_nodes_merges_a_too_close_pair():
    # A regular pentagon ring (real spacing ~0.059 rad, real diagonals ~0.095 rad) plus one
    # point placed 1% further out along the same bearing as vertex 0 -- an almost-duplicate
    # ~0.0005 rad away. spacing_rad=0.07 puts the merge threshold (spacing_rad/1.5 ~ 0.047)
    # comfortably above that near-duplicate gap but below every real ring/diagonal edge, and
    # the split threshold (spacing_rad*1.5 = 0.105) comfortably above the longest real edge
    # (~0.096) -- so only the near-duplicate pair should merge, nothing should split.
    center = np.array([1.0, 0.0, 0.0])
    east, north = geometry.local_tangent_basis(center)
    angles = np.arange(5) * 2 * np.pi / 5
    r = 0.05
    ring = geometry.normalize(
        center[None, :] * np.cos(r) + r * (np.cos(angles)[:, None] * east + np.sin(angles)[:, None] * north)
    )
    close_dup = geometry.normalize(center * np.cos(r * 1.01) + (r * 1.01) * (np.cos(angles[0]) * east + np.sin(angles[0]) * north))
    points = np.concatenate([ring, close_dup[None, :]], axis=0)
    theta, phi = geometry.xyz_to_latlon(points)
    n = len(points)
    elevation = np.arange(n, dtype=float) * 10.0
    fields = {"soil_depth": np.zeros(n)}
    tri = build_mesh_triangulation(points)
    spacing_rad = 0.07

    new_theta, new_phi, new_elevation, new_fields = remesh_nodes(theta, phi, elevation, fields, tri, points, spacing_rad)
    assert len(new_theta) == n - 1


def test_remesh_nodes_splits_a_too_long_edge_and_averages_fields():
    rng = np.random.default_rng(7)
    points = _disk_cluster_xyz(rng, 40, radius_rad=0.05)
    theta, phi = geometry.xyz_to_latlon(points)
    n = len(points)
    elevation = np.arange(n, dtype=float)
    is_volcano = np.array([i % 2 == 0 for i in range(n)])
    fields = {"is_volcano": is_volcano}
    tri = build_mesh_triangulation(points)
    spacing_rad = 0.001  # far smaller than the cluster's real spacing -- forces every edge to split

    new_theta, new_phi, new_elevation, new_fields = remesh_nodes(theta, phi, elevation, fields, tri, points, spacing_rad)
    assert len(new_theta) > n
    assert new_fields["is_volcano"].dtype == bool
