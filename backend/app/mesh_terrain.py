"""Pure array/geometry functions behind `PlateWithMesh` (plates.py): building a Delaunay
triangle mesh over a plate's own local unit-sphere node positions, tracing its boundary,
estimating an outward growth direction at a boundary vertex, deduplicating simultaneously-
grown candidate nodes, and a remeshing quality pass. Nothing here touches `Plate`/`self` --
same reasoning `elevation_lines.py` already documents for why plate-independent node/line
machinery lives in its own sibling module rather than in plates.py itself.

Triangulation is built via `scipy.spatial.ConvexHull` on a plate's own local unit-sphere
positions (mirrors geodesic.py's `_icosahedron_faces`): every node lies on the unit sphere,
so every node is automatically a hull vertex, and the resulting hull triangulation is exactly
the sphere's own Delaunay triangulation of those points -- a standard duality that sidesteps
the wraparound problems a flat 2D Delaunay over local (theta, phi) would hit for a plate
straddling its own local theta=+-pi seam.

A plate only covers a spherical cap, though (unlike geodesic.py's whole-sphere icosahedron),
and `ConvexHull` on a cap-shaped point set still produces a *closed* polyhedron: it bridges
the empty far side of the sphere with large spurious "lid" triangles connecting only the
outer rim points. This means a naive "an edge belongs to exactly one triangle => boundary
edge" test finds no boundary at all (every true rim edge is shared by one real triangle and
one lid triangle, so every edge has two incident triangles). `build_mesh_triangulation` fixes
this with an alpha-shape-style filter: drop any triangle with an edge longer than
`MESH_ALPHA_EDGE_FACTOR` times the median hull-edge length (real local edges vastly
outnumber the handful of long lid edges, so the median is a robust local-spacing estimate),
before computing adjacency/boundary from the surviving triangles only.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import ConvexHull, QhullError, cKDTree

from . import geometry

# How many times the median hull-edge length a triangle's longest edge may be before it's
# treated as a spurious "lid" facet bridging the empty far side of a plate's spherical cap
# (see module docstring) and dropped. Same "generous multiple of local spacing" spirit as
# plates.OUTLINE_NEIGHBORHOOD_SPACING_FACTOR (3.0) -- wide enough to keep real, only-mildly-
# irregular local edges, narrow enough to reject the lid triangles, which span a large
# fraction of the whole plate rather than one local gap.
#
# Known limitation, confirmed directly against small synthetic plates (~60 nodes): the
# global median-edge-length threshold this drives can end up coarser than the true boundary
# for a *small/compact* point cloud, where a real (non-lid) rim edge can itself run close to
# this factor's own multiple of the interior spacing -- `build_mesh_triangulation` then
# correctly keeps that triangle (it's real, not a lid facet), but the effect is a smaller,
# sparser `boundary_vertices` set than an intuitive "every rim node" count would suggest
# (observed as few as ~6 boundary vertices out of 61 total nodes for one such synthetic
# case). `outline_world`/`contains` stay geometrically correct either way -- the traced
# boundary loop is still a real, non-overclaiming edge loop -- but `deform()`'s own boundary
# growth/shrink has fewer eligible vertices to act on per call for a small/young plate,
# meaning slower per-step boundary evolution until the plate grows large enough for its own
# rim to be a small fraction of its total edges (confirmed clean at 200+ nodes in
# unit_tests/test_mesh_terrain.py). Accepted for now, same "bounded, documented imprecision"
# tradeoff PlateWithRTree.outline_world's own docstring already accepts for concave plates,
# rather than a per-triangle aspect-ratio filter or other local refinement.
MESH_ALPHA_EDGE_FACTOR = 4.0

# A kept mesh edge outside [spacing_rad / MESH_REMESH_TOLERANCE, spacing_rad *
# MESH_REMESH_TOLERANCE] triggers a remesh pass -- same value and role as
# elevation_lines.IRREGULARITY_TOLERANCE (the row-based `needs_regularizing`'s own tolerance).
MESH_REMESH_TOLERANCE = 1.5

# Below this many nodes there aren't enough points for a meaningful 3D convex hull -- treat
# the triangulation as degenerate (every node is its own boundary) rather than let
# ConvexHull raise. Matches plates.OUTLINE_MIN_NODES_FOR_HULL's own value/reasoning; kept as
# a separate constant (not imported) to avoid a plates.py <-> mesh_terrain.py import cycle --
# callers in plates.py that want to stay in sync with that constant can pass it explicitly.
MIN_NODES_FOR_TRIANGULATION = 4


@dataclass
class MeshTriangulation:
    """A plate's current mesh, built from its own local unit-sphere node positions.

    `triangles`: (M, 3) int vertex-index triples, outward-wound, spurious "lid" facets
    already excluded (see module docstring). `adjacency`: vertex index -> the set of vertex
    indices it shares a kept triangle edge with. `boundary_edges`: edges belonging to exactly
    one kept triangle, plus every edge touching an isolated vertex (one that ended up in zero
    kept triangles -- a genuine outlier, far from everything else; trivially its own
    boundary). `boundary_vertices`: sorted unique vertex indices touched by `boundary_edges`.
    `degenerate`: True if there were too few nodes for a real hull, or the hull was
    numerically degenerate (collinear/near-coplanar points) -- `triangles`/`adjacency`/
    `boundary_edges` are then empty and `boundary_vertices` is every node index."""

    triangles: np.ndarray
    adjacency: dict[int, set[int]]
    boundary_edges: set[frozenset[int]]
    boundary_vertices: np.ndarray
    degenerate: bool


def _degenerate_triangulation(n: int) -> MeshTriangulation:
    return MeshTriangulation(
        triangles=np.zeros((0, 3), dtype=int),
        adjacency={i: set() for i in range(n)},
        boundary_edges=set(),
        boundary_vertices=np.arange(n),
        degenerate=True,
    )


def build_mesh_triangulation(local_xyz: np.ndarray, min_nodes: int = MIN_NODES_FOR_TRIANGULATION) -> MeshTriangulation:
    """Build the Delaunay triangle mesh over `local_xyz` (a plate's own local unit-sphere
    node positions, e.g. `geometry.local_xyz(plate.phi, plate.theta)`) -- see module
    docstring for the 3D-convex-hull-on-a-sphere technique and the alpha-shape lid-facet
    filter. `degenerate=True` (every node its own boundary, no triangles) for `< min_nodes`
    nodes or a numerically degenerate (collinear/near-coplanar) point set."""
    n = len(local_xyz)
    if n < min_nodes:
        return _degenerate_triangulation(n)
    try:
        hull = ConvexHull(local_xyz)
    except QhullError:
        return _degenerate_triangulation(n)

    triangles = hull.simplices.copy()
    tri_pts = local_xyz[triangles]  # (M, 3, 3)
    normals = np.cross(tri_pts[:, 1] - tri_pts[:, 0], tri_pts[:, 2] - tri_pts[:, 0])
    centroids = tri_pts.mean(axis=1)
    flip = np.sum(normals * centroids, axis=-1) < 0
    triangles[flip, 1], triangles[flip, 2] = triangles[flip, 2].copy(), triangles[flip, 1].copy()

    tri_pts = local_xyz[triangles]  # re-read post winding-fix (vertex identities unchanged, order may have)
    ab = geometry.angular_distance(tri_pts[:, 0], tri_pts[:, 1])
    bc = geometry.angular_distance(tri_pts[:, 1], tri_pts[:, 2])
    ca = geometry.angular_distance(tri_pts[:, 2], tri_pts[:, 0])
    edge_lengths = np.stack([ab, bc, ca], axis=1)
    typical_spacing = float(np.median(edge_lengths)) if edge_lengths.size else 0.0
    max_edge = MESH_ALPHA_EDGE_FACTOR * typical_spacing
    kept_mask = edge_lengths.max(axis=1) <= max_edge if max_edge > 0 else np.zeros(len(triangles), dtype=bool)
    kept_triangles = triangles[kept_mask]

    adjacency: dict[int, set[int]] = {i: set() for i in range(n)}
    edge_count: dict[frozenset[int], int] = {}
    for a, b, c in kept_triangles:
        for u, w in ((int(a), int(b)), (int(b), int(c)), (int(c), int(a))):
            adjacency[u].add(w)
            adjacency[w].add(u)
            edge = frozenset((u, w))
            edge_count[edge] = edge_count.get(edge, 0) + 1

    boundary_edges = {edge for edge, count in edge_count.items() if count == 1}
    boundary_vertex_set: set[int] = set()
    for edge in boundary_edges:
        boundary_vertex_set.update(edge)
    for i in range(n):
        if not adjacency[i]:
            boundary_vertex_set.add(i)

    return MeshTriangulation(
        triangles=kept_triangles,
        adjacency=adjacency,
        boundary_edges=boundary_edges,
        boundary_vertices=np.array(sorted(boundary_vertex_set), dtype=int),
        degenerate=False,
    )


def trace_boundary_loops(boundary_edges: set[frozenset[int]]) -> list[list[int]]:
    """Walk `boundary_edges` (see `MeshTriangulation.boundary_edges`) into ordered cyclic
    vertex loops -- generically one loop for a simply-connected plate, but a shrinking/
    growing plate can genuinely pinch into more than one disjoint piece, so this always
    returns a list (possibly empty, for no boundary at all -- e.g. a plate covering the
    entire sphere)."""
    adjacency: dict[int, list[int]] = {}
    for edge in boundary_edges:
        a, b = tuple(edge)
        adjacency.setdefault(a, []).append(b)
        adjacency.setdefault(b, []).append(a)

    visited_edges: set[frozenset[int]] = set()
    loops: list[list[int]] = []
    for start in adjacency:
        for first_step in adjacency[start]:
            edge = frozenset((start, first_step))
            if edge in visited_edges:
                continue
            visited_edges.add(edge)
            loop = [start]
            prev, cur = start, first_step
            while cur != start:
                loop.append(cur)
                candidates = [w for w in adjacency[cur] if w != prev]
                if not candidates:
                    break  # non-manifold guard -- shouldn't happen post-filtering
                next_vertex = candidates[0]
                visited_edges.add(frozenset((cur, next_vertex)))
                prev, cur = cur, next_vertex
            loops.append(loop)
    return loops


def estimate_outward_direction(vertex: int, tri: MeshTriangulation, local_xyz: np.ndarray) -> np.ndarray | None:
    """A local outward tangent-plane unit direction at boundary vertex `vertex`, for
    `deform()`'s boundary growth -- purely local (no global boundary-loop normal or
    convexity assumption needed), so it stays correct at a concave notch.

    Primary estimate: the tangent between `vertex`'s own two boundary-loop neighbors,
    rotated 90 degrees in the tangent plane (`cross(vertex_position, tangent)`), oriented
    away from `vertex`'s own 1-ring centroid. Falls back to "away from the 1-ring centroid"
    directly if there aren't two boundary neighbors (or the tangent is degenerate), and to
    `None` (caller's own further fallback -- away from the whole plate's centroid) if there's
    no adjacency information at all (a fully isolated vertex)."""
    p_vertex = local_xyz[vertex]
    boundary_neighbours = [w for edge in tri.boundary_edges if vertex in edge for w in edge if w != vertex]
    ring = tri.adjacency.get(vertex) or set(boundary_neighbours)
    if not ring:
        return None

    direction = None
    if len(boundary_neighbours) >= 2:
        p_prev, p_next = local_xyz[boundary_neighbours[0]], local_xyz[boundary_neighbours[1]]
        tangent = p_next - p_prev
        tangent = tangent - np.dot(tangent, p_vertex) * p_vertex
        tangent_norm = float(np.linalg.norm(tangent))
        if tangent_norm > 1e-9:
            direction = np.cross(p_vertex, tangent / tangent_norm)

    ring_mean = geometry.normalize(np.mean(local_xyz[list(ring)], axis=0))
    if direction is None:
        direction = p_vertex - ring_mean

    direction = direction - np.dot(direction, p_vertex) * p_vertex
    direction_norm = float(np.linalg.norm(direction))
    if direction_norm < 1e-9:
        return None
    direction = direction / direction_norm

    if np.dot(direction, p_vertex - ring_mean) < 0:
        direction = -direction
    return direction


def dedupe_growth_candidates(candidates_xyz: np.ndarray, existing_xyz: np.ndarray, min_sep_rad: float) -> np.ndarray:
    """Boolean keep-mask over `candidates_xyz` (unit vectors): drops any candidate within
    `min_sep_rad` of an existing node, then greedily drops any later candidate within
    `min_sep_rad` of an earlier-surviving one (deterministic: lower index always wins) --
    needed because, unlike `PlateWithLines`' own single-row growth, several independent
    boundary vertices grow "in parallel" in one `deform()` call and their outward rays can
    converge on the same real-world spot."""
    n = len(candidates_xyz)
    keep = np.ones(n, dtype=bool)
    if n == 0:
        return keep
    if len(existing_xyz) > 0:
        dist_to_existing, _ = cKDTree(existing_xyz).query(candidates_xyz)
        keep &= dist_to_existing >= min_sep_rad
    alive_idx = np.nonzero(keep)[0]
    if len(alive_idx) < 2:
        return keep
    pairs = sorted(cKDTree(candidates_xyz[alive_idx]).query_pairs(min_sep_rad))
    for i, j in pairs:
        global_i, global_j = int(alive_idx[i]), int(alive_idx[j])
        if keep[global_i] and keep[global_j]:
            keep[global_j] = False
    return keep


def _unique_edges(tri: MeshTriangulation) -> set[frozenset[int]]:
    return {frozenset((u, w)) for u, neighbours in tri.adjacency.items() for w in neighbours}


def needs_remeshing(tri: MeshTriangulation, local_xyz: np.ndarray, spacing_rad: float, tolerance: float = MESH_REMESH_TOLERANCE) -> bool:
    """True if any kept mesh edge has drifted outside `[spacing_rad / tolerance, spacing_rad
    * tolerance]` -- the mesh analogue of `elevation_lines.needs_regularizing`'s per-line gap
    check, triggering `remesh_nodes` below. Mirrors the *real* trigger `PlateWithLines.deform`
    uses: inline, every `deform()` call, gated by an irregularity check -- not a periodic
    cadence (elevation_lines.regularize_world_lines/regularize_plate_lines, which use that
    framing, are themselves dead code with no live caller)."""
    if tri.degenerate:
        return False
    edges = _unique_edges(tri)
    if not edges:
        return False
    for edge in edges:
        a, b = tuple(edge)
        length = float(geometry.angular_distance(local_xyz[a], local_xyz[b]))
        if length > tolerance * spacing_rad or length < spacing_rad / tolerance:
            return True
    return False


def remesh_nodes(
    theta: np.ndarray,
    phi: np.ndarray,
    elevation: np.ndarray,
    fields: dict[str, np.ndarray],
    tri: MeshTriangulation,
    local_xyz: np.ndarray,
    spacing_rad: float,
    tolerance: float = MESH_REMESH_TOLERANCE,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    """The mesh analogue of `elevation_lines.regularize_line`: repair edge-length drift the
    per-step boundary growth/shrink doesn't fully prevent. Two passes over `tri`'s own kept
    edges (see `MeshTriangulation.adjacency`) -- merge any pair closer than `spacing_rad /
    tolerance` (deterministically dropping the higher-index endpoint), then split any
    surviving edge (both endpoints still alive) longer than `spacing_rad * tolerance` by
    inserting a midpoint node (position = normalized average; every field value = the
    endpoints' arithmetic mean, except a `bool`-dtype field like `is_volcano`, which is their
    logical AND -- a new node counts as a volcano only if both parents did). Returns the
    updated (theta, phi, elevation, fields) arrays; the caller is responsible for committing
    them (e.g. via `PlateWithMesh.set_nodes`, which also rebuilds the triangulation)."""
    n = len(theta)
    edges = _unique_edges(tri)
    if not edges:
        return theta, phi, elevation, fields

    merge_threshold = spacing_rad / tolerance
    edges_by_length = sorted(edges, key=lambda e: float(geometry.angular_distance(local_xyz[min(e)], local_xyz[max(e)])))
    alive = np.ones(n, dtype=bool)
    for edge in edges_by_length:
        a, b = min(edge), max(edge)
        length = float(geometry.angular_distance(local_xyz[a], local_xyz[b]))
        if length >= merge_threshold:
            break  # sorted ascending -- nothing further can be too-short either
        if alive[a] and alive[b]:
            alive[b] = False

    split_threshold = tolerance * spacing_rad
    new_theta: list[float] = []
    new_phi: list[float] = []
    new_elevation: list[float] = []
    new_fields: dict[str, list] = {name: [] for name in fields}
    for edge in edges:
        a, b = min(edge), max(edge)
        if not (alive[a] and alive[b]):
            continue
        length = float(geometry.angular_distance(local_xyz[a], local_xyz[b]))
        if length <= split_threshold:
            continue
        midpoint = geometry.normalize(local_xyz[a] + local_xyz[b])
        mid_phi, mid_theta = geometry.xyz_to_latlon(midpoint[None, :])
        new_theta.append(float(mid_theta[0]))
        new_phi.append(float(mid_phi[0]))
        new_elevation.append(float((elevation[a] + elevation[b]) / 2.0))
        for name, values in fields.items():
            if values.dtype == bool:
                new_fields[name].append(bool(values[a]) and bool(values[b]))
            else:
                new_fields[name].append(float((values[a] + values[b]) / 2.0))

    final_theta = np.concatenate([theta[alive], np.array(new_theta)])
    final_phi = np.concatenate([phi[alive], np.array(new_phi)])
    final_elevation = np.concatenate([elevation[alive], np.array(new_elevation)])
    final_fields = {
        name: np.concatenate([values[alive], np.array(new_fields[name], dtype=values.dtype)])
        for name, values in fields.items()
    }
    return final_theta, final_phi, final_elevation, final_fields
