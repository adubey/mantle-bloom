# Hex grid export format

`POST /world/export_hexgrid` (see [api-reference.md](api-reference.md)) tiles the sphere
into a geodesic-icosahedron hex/pentagon dome -- the standard "hexasphere"/Goldberg-
polyhedron-dual construction used by hex-globe strategy games -- independent of the plate
simulation's own node cloud, and maps the current world's elevation/ocean/biome onto each
tile. See `backend/app/geodesic.py` for the construction itself; this document is about the
*output file* an external application would actually consume.

## Shape

```json
{
  "planet_radius_km": 6371.0,
  "frequency": 16,
  "num_tiles": 2562,
  "vertices": [[x, y, z], ...],
  "tiles": [
    {
      "id": 0,
      "is_pentagon": false,
      "center_lat_deg": 12.3,
      "center_lon_deg": -45.6,
      "center_xyz": [x, y, z],
      "corner_vertex_ids": [14, 87, 203, 55, 9, 301],
      "elevation_m": 812.4,
      "is_ocean": false,
      "biome": "Humid Subtropical",
      "neighbor_ids": [1, 4, 7, 12, 18, 23]
    }
  ]
}
```

- `vertices` is a shared, globally-indexed list of tile *corner* points (unit vectors,
  `x^2+y^2+z^2 = 1`) -- multiply by `planet_radius_km` for real km. This is a **separate** id
  space from tile ids: `vertices` has `20 * frequency**2` entries (one per triangle in the
  underlying subdivided icosahedron -- see geodesic.py), while `tiles` has
  `10 * frequency**2 + 2` entries (one per tile).
- Every tile is a pentagon (`corner_vertex_ids`/`neighbor_ids` both length 5) or a hexagon
  (length 6) -- exactly the 12 original icosahedron vertices are pentagons, everything else
  is a hexagon.
- `corner_vertex_ids` is wound consistently CCW as seen from outside the sphere (outward-
  facing normal), indices into the top-level `vertices` array.
- `neighbor_ids` is **index-aligned** with `corner_vertex_ids`: `neighbor_ids[k]` is the tile
  sharing the edge between `corner_vertex_ids[k]` and `corner_vertex_ids[(k+1) % n]`.
- `elevation_m`/`is_ocean`/`biome` are sampled from the world exactly as it stood at export
  time -- `is_ocean` is `elevation_m <= ` the world's current sea level, `biome` is one of
  `backend/app/biomes.py`'s `BIOME_NAMES` (a descriptive Köppen-Geiger climate class such as
  `"Hot Desert"` or `"Subarctic (Boreal)"` for land, or a pelagic class such as
  `"Subtropical Gyre"` for ocean -- the same categories the Biome/Combined map views use).

## Finding a tile's neighbors

### Direct (the common case)

`neighbor_ids` is already precomputed at export time -- just look it up:

```
function neighborsOf(tile):
    return tile.neighbor_ids.map(id => tileById[id])
```

### Geometric re-derivation, from the dome's own shared-vertex structure

Useful if a client only kept `vertices`/`corner_vertex_ids` -- e.g. after re-processing the
export into its own data structure and dropping `neighbor_ids` -- or wants to verify/rebuild
adjacency independently. Two tiles are neighbors *iff* their corner lists share an edge (two
consecutive corner vertex ids). This is exact and pentagon-safe: unlike testing "are these
two tile centers close together" with a fixed distance threshold, it never misclassifies at
the 12 pentagon sites, where a tile's neighbors sit at a measurably different distance than a
hexagon's do.

```
edgeToTiles = {}                      // key: sorted (vertex_id, vertex_id) pair
for tile in tiles:
    n = len(tile.corner_vertex_ids)
    for k in 0..n-1:
        a, b = tile.corner_vertex_ids[k], tile.corner_vertex_ids[(k+1) % n]
        key = (min(a, b), max(a, b))
        edgeToTiles.setdefault(key, []).append(tile.id)

function neighborsOf(tile):
    result = []
    n = len(tile.corner_vertex_ids)
    for k in 0..n-1:
        a, b = tile.corner_vertex_ids[k], tile.corner_vertex_ids[(k+1) % n]
        pair = edgeToTiles[(min(a, b), max(a, b))]
        result.append(pair[0] if pair[0] != tile.id else pair[1])
    return result   // index-aligned with corner_vertex_ids, i.e. result[k] is the
                     // neighbor across the edge between corner k and corner k+1
```

Building `edgeToTiles` once up front is `O(num_tiles)` (every tile contributes exactly as
many edges as it has corners, 5 or 6), and every `neighborsOf` lookup after that is `O(1)`
per neighbor.
