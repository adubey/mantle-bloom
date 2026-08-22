"""Geodesic icosahedral hex/pentagon tiling of the sphere ("File > Export Hex Grid") --
independent of the plate/elevation simulation's own node cloud; samples the *current*
world's elevation/biome onto a separate, uniform tiling built for an external application
to consume (see docs/hex-export-format.md for the file format and the neighbor-finding
pseudocode this module's own `neighbor_ids` field backs).

Construction (the standard "hexasphere"/Goldberg-polyhedron-dual technique):
  1. A regular icosahedron (12 vertices, 20 triangular faces) -- see `_icosahedron_vertices`/
     `_icosahedron_faces`.
  2. Each face subdivided into frequency**2 smaller triangles (Class I geodesic
     subdivision) -- see `_subdivide`.
  3. The *dual* of that subdivided mesh: one tile per subdivided-mesh vertex, its polygon
     boundary the centroids of its incident triangles (5 triangles -- a pentagon -- at the
     original 12 icosahedron vertices, 6 -- a hexagon -- everywhere else) -- see
     `_build_dual`.

Tile adjacency falls directly out of the subdivided mesh's own edges (two tiles are
neighbors iff the subdivided mesh has an edge between their two originating vertices), so
it's derived exactly, not approximated by a nearest-center distance threshold (which would
misclassify the 12 pentagon sites -- see docs/hex-export-format.md).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.spatial import ConvexHull, cKDTree

from . import biomes, climate, geometry, plates
from .world import World

PLANET_RADIUS_KM = plates.PLANET_RADIUS_KM

# UI-facing choices for `frequency` (the geodesic subdivision level) -- a discrete set, not
# a free-form input, same "pick from a few sane presets" reasoning plates.NODE_DENSITY_CHOICES
# and climate.CLIMATE_DENSITY_CHOICES already use. Tile count is exactly 10*frequency**2 + 2
# (see tile_count below) -- 642 / 2562 / 10242 tiles respectively.
FREQUENCY_CHOICES = (8, 16, 32)
DEFAULT_FREQUENCY = 16


def tile_count(frequency: int) -> int:
    """Standard geodesic-icosahedron vertex count -- see module docstring / _subdivide."""
    return 10 * frequency * frequency + 2


def _icosahedron_vertices() -> np.ndarray:
    """12 unit vectors -- the standard construction: cyclic permutations of (0, +-1, +-phi)."""
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    raw = []
    for a in (-1.0, 1.0):
        for b in (-1.0, 1.0):
            raw.append((0.0, a, b * phi))
            raw.append((a, b * phi, 0.0))
            raw.append((b * phi, 0.0, a))
    return geometry.normalize(np.array(raw, dtype=float))


def _icosahedron_faces(vertices: np.ndarray) -> np.ndarray:
    """20 triangular faces (vertex-index triples), each wound outward/CCW (as seen from
    outside the sphere). `scipy.spatial.ConvexHull` on the 12 vertices gives the facets
    directly -- the icosahedron is simplicial (no coplanar quads to merge), so every facet
    is already a triangle -- but doesn't guarantee outward winding, so that's fixed up per
    face below (sign of the face normal dotted with its own centroid direction)."""
    hull = ConvexHull(vertices)
    faces = hull.simplices.copy()
    for face in faces:
        a, b, c = vertices[face]
        normal = np.cross(b - a, c - a)
        centroid = (a + b + c) / 3.0
        if np.dot(normal, centroid) < 0:
            face[1], face[2] = face[2], face[1]
    return faces


def _subdivide(vertices: np.ndarray, faces: np.ndarray, frequency: int) -> tuple[np.ndarray, np.ndarray]:
    """Class I geodesic subdivision at `frequency` -- see module docstring. Returns
    (subdivided_vertices (10*frequency**2+2, 3), triangles (20*frequency**2, 3) int indices
    into subdivided_vertices), every triangle wound the same outward/CCW sense its parent
    icosahedron face was (barycentric interpolation preserves winding)."""
    key_scale = 1.0e7  # dedup precision on unit-sphere coordinates -- far finer than the
    # floating-point gap between two independently-computed copies of the same shared-edge
    # point (two neighboring faces both interpolate toward their shared A/B/C corners).
    vertex_by_key: dict[tuple[int, int, int], int] = {}
    all_vertices: list[np.ndarray] = []

    def vertex_id(point: np.ndarray) -> int:
        key = tuple(np.round(point * key_scale).astype(np.int64).tolist())
        existing = vertex_by_key.get(key)
        if existing is not None:
            return existing
        new_id = len(all_vertices)
        vertex_by_key[key] = new_id
        all_vertices.append(point)
        return new_id

    triangles: list[tuple[int, int, int]] = []
    for ia, ib, ic in faces:
        a, b, c = vertices[ia], vertices[ib], vertices[ic]
        # idx[i][j] -- global vertex id of the barycentric grid point i steps toward B, j
        # steps toward C (i, j >= 0, i + j <= frequency; i = j = 0 is A itself).
        idx = [[-1] * (frequency - i + 1) for i in range(frequency + 1)]
        for i in range(frequency + 1):
            for j in range(frequency + 1 - i):
                k = frequency - i - j
                point = geometry.normalize((k * a + i * b + j * c) / frequency)
                idx[i][j] = vertex_id(point)
        for i in range(frequency):
            for j in range(frequency - i):
                triangles.append((idx[i][j], idx[i + 1][j], idx[i][j + 1]))
                if i + j + 1 < frequency:
                    triangles.append((idx[i + 1][j], idx[i + 1][j + 1], idx[i][j + 1]))

    return np.array(all_vertices), np.array(triangles, dtype=int)


@dataclass
class HexTile:
    tile_id: int
    is_pentagon: bool
    center_xyz: np.ndarray  # unit vector
    corner_vertex_ids: list[int]  # indices into _build_dual's returned dual-vertex array
    # Index-aligned with corner_vertex_ids: neighbor_ids[k] is the tile sharing the edge
    # between corner_vertex_ids[k] and corner_vertex_ids[(k+1) % len] -- see module
    # docstring and docs/hex-export-format.md.
    neighbor_ids: list[int]


def _build_dual(subdivided_vertices: np.ndarray, triangles: np.ndarray) -> tuple[np.ndarray, list[HexTile]]:
    """The dual mesh: one tile per subdivided-mesh vertex, its corners the (cyclically
    ordered, outward-wound) centroids of its incident triangles -- see module docstring.
    Returns (dual_vertices (num_triangles, 3) -- one per triangle, its own centroid -- and
    the per-tile list)."""
    dual_vertices = geometry.normalize(subdivided_vertices[triangles].mean(axis=1))

    incident: dict[int, list[int]] = {i: [] for i in range(len(subdivided_vertices))}
    for tri_id, (i, j, k) in enumerate(triangles):
        incident[i].append(tri_id)
        incident[j].append(tri_id)
        incident[k].append(tri_id)

    tiles: list[HexTile] = []
    for v in range(len(subdivided_vertices)):
        v_pos = subdivided_vertices[v]
        tri_ids = incident[v]
        east, north = geometry.local_tangent_basis(v_pos)

        def angle_of(tri_id: int) -> float:
            c = dual_vertices[tri_id]
            tangent = c - np.dot(c, v_pos) * v_pos
            return math.atan2(np.dot(tangent, north), np.dot(tangent, east))

        ordered = sorted(tri_ids, key=angle_of)

        # Force outward/CCW as seen from outside the sphere, checked directly (rather than
        # relied on from local_tangent_basis's own handedness convention, which this module
        # doesn't otherwise depend on): for an outward-wound polygon, the cross product of
        # the first two consecutive corner-edges should point the same way as v_pos itself.
        c0, c1 = dual_vertices[ordered[0]], dual_vertices[ordered[1]]
        if np.dot(np.cross(c0 - v_pos, c1 - v_pos), v_pos) < 0:
            ordered.reverse()

        n = len(ordered)
        neighbor_ids = []
        for pos in range(n):
            tri_a, tri_b = triangles[ordered[pos]], triangles[ordered[(pos + 1) % n]]
            shared = (set(tri_a.tolist()) & set(tri_b.tolist())) - {v}
            neighbor_ids.append(shared.pop())

        tiles.append(
            HexTile(tile_id=v, is_pentagon=(n == 5), center_xyz=v_pos, corner_vertex_ids=ordered, neighbor_ids=neighbor_ids)
        )

    return dual_vertices, tiles


def _hex_slope(centers_xyz: np.ndarray, elevation_m: np.ndarray, tiles: list[HexTile]) -> np.ndarray:
    """Dimensionless rise/run slope at each tile center -- real elevation difference to its
    steepest neighbor, divided by the real great-circle distance to that neighbor. Uses the
    dome's own tile adjacency (already exact) rather than render_image.grid_slope's
    regular-lat/lon-grid neighbor assumption, which doesn't apply to this irregular tiling.
    Feeds biomes.classify_wetland the same way render_image._biome_fields's own slope does."""
    slope = np.zeros(len(tiles))
    for tile in tiles:
        i = tile.tile_id
        neighbor_idx = np.array(tile.neighbor_ids, dtype=int)
        dist_km = geometry.angular_distance(centers_xyz[i], centers_xyz[neighbor_idx]) * PLANET_RADIUS_KM
        rise_m = np.abs(elevation_m[i] - elevation_m[neighbor_idx])
        slope[i] = np.max(rise_m / np.maximum(dist_km * 1000.0, 1.0))
    return slope


def export_hexgrid(world: World, frequency: int) -> dict:
    """Builds a `frequency`-level geodesic dome (see module docstring), samples the world's
    current elevation/ocean/biome onto each tile center using the same nearest-node/
    climate-resample techniques render_image.py's Biome/Combined views already use (see
    render_image._biome_fields), and returns the JSON-ready export dict (see
    docs/hex-export-format.md)."""
    ico_vertices = _icosahedron_vertices()
    ico_faces = _icosahedron_faces(ico_vertices)
    subdivided_vertices, triangles = _subdivide(ico_vertices, ico_faces, frequency)
    dual_vertices, tiles = _build_dual(subdivided_vertices, triangles)

    centers_xyz = subdivided_vertices
    center_lat, center_lon = geometry.xyz_to_latlon(centers_xyz)

    collected = plates.collect_all_points(world.plates)
    if collected is None:
        elevation_m = np.zeros(len(tiles))
    else:
        all_points, all_elevation, _ = collected
        _, idx = cKDTree(all_points).query(centers_xyz)
        elevation_m = all_elevation[idx]
    is_ocean = elevation_m <= world.sea_level_m

    climate_fields = climate.compute_climate_cached(world)
    climate_tree = cKDTree(climate_fields.world_xyz.reshape(-1, 3))
    _, climate_idx = climate_tree.query(centers_xyz)
    air_temp_c = climate_fields.air_temperature_c.reshape(-1)[climate_idx]
    precip_mm = climate_fields.precipitation_mm.reshape(-1)[climate_idx]

    slope = _hex_slope(centers_xyz, elevation_m, tiles)
    biome_ids = biomes.classify_biomes(air_temp_c, precip_mm, elevation_m, slope, is_ocean, world.sea_level_m)

    return {
        "planet_radius_km": PLANET_RADIUS_KM,
        "frequency": frequency,
        "num_tiles": len(tiles),
        "vertices": np.round(dual_vertices, 6).tolist(),
        "tiles": [
            {
                "id": tile.tile_id,
                "is_pentagon": tile.is_pentagon,
                "center_lat_deg": round(math.degrees(center_lat[tile.tile_id]), 6),
                "center_lon_deg": round(math.degrees(center_lon[tile.tile_id]), 6),
                "center_xyz": np.round(tile.center_xyz, 6).tolist(),
                "corner_vertex_ids": tile.corner_vertex_ids,
                "elevation_m": round(float(elevation_m[tile.tile_id]), 3),
                "is_ocean": bool(is_ocean[tile.tile_id]),
                "biome": biomes.BIOME_NAMES[int(biome_ids[tile.tile_id])],
                "neighbor_ids": tile.neighbor_ids,
            }
            for tile in tiles
        ],
    }
