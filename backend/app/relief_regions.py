"""Named real (Earth/Pangaea) and lore-described (Z&D) mountain-belt/plateau regions, for the
Premade-worlds "local relief grounded in real geography" feature -- see
`terrain_noise.ContinentalRelief`'s `belt_mask`/`plateau_mask` override hooks. Only *where* a
belt/plateau is allowed to appear becomes geography-driven this way; the ridge/terrace texture
within one is still the same algorithmic, seeded noise as always (see terrain_noise.py) --
these regions place mountains in the right rough area, they don't hand-draw individual peaks.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
from scipy.spatial import cKDTree

from . import geometry, real_plates

Ring = list[tuple[float, float]]  # (lon_deg, lat_deg) control points

# Buffer widths (radians) the belt/plateau masks fall off over, past a belt's own centerline
# or a plateau's own edge -- wide enough for a smooth transition, narrow enough that a belt
# still reads as a belt rather than a smear across half the continent.
_BELT_BUFFER_RAD = 0.12
_PLATEAU_EDGE_BUFFER_RAD = 0.18
# How many points a belt polyline segment is densified to (see _densify_ring) -- densifying
# in lon/lat before projecting to xyz, so a KDTree nearest-neighbour distance is accurate
# along a segment, not just at its endpoints.
_DENSIFY_PER_SEGMENT = 6


def _smoothstep(lo: float, hi: float, x: np.ndarray) -> np.ndarray:
    t = np.clip((x - lo) / (hi - lo), 0.0, 1.0)
    return t * t * (3.0 - 2.0 * t)


def _ring_to_xyz(ring: Ring) -> np.ndarray:
    arr = np.array(ring, dtype=float)
    return geometry.latlon_to_xyz(np.radians(arr[:, 1]), np.radians(arr[:, 0]))


def _densify_ring(ring: Ring, closed: bool = False) -> Ring:
    pts = list(ring) + ([ring[0]] if closed else [])
    out: list[tuple[float, float]] = []
    for (lon0, lat0), (lon1, lat1) in zip(pts, pts[1:]):
        for t in np.linspace(0.0, 1.0, _DENSIFY_PER_SEGMENT, endpoint=False):
            out.append((lon0 + (lon1 - lon0) * t, lat0 + (lat1 - lat0) * t))
    out.append(pts[-1])
    return out


def _chord_for_angle(rad: float) -> float:
    """Chordal (3D euclidean) distance on the unit sphere corresponding to angular distance
    `rad` -- for converting an angular buffer into the units a plain cKDTree over xyz points
    actually measures in."""
    return 2.0 * np.sin(rad / 2.0)


def build_belt_mask(belts: dict[str, Ring], buffer_rad: float = _BELT_BUFFER_RAD) -> Callable[[np.ndarray], np.ndarray]:
    """A `ContinentalRelief(belt_mask=...)`-compatible callable: 1.0 on a belt's own
    centerline, smoothly falling to 0.0 by `buffer_rad` away from it."""
    all_pts = np.concatenate([_ring_to_xyz(_densify_ring(ring)) for ring in belts.values()])
    tree = cKDTree(all_pts)
    chord_buffer = _chord_for_angle(buffer_rad)

    def mask(xyz: np.ndarray) -> np.ndarray:
        dist, _ = tree.query(np.asarray(xyz, dtype=float), workers=-1)
        return _smoothstep(chord_buffer, 0.0, dist)

    return mask


def build_plateau_mask(plateaus: dict[str, Ring], edge_buffer_rad: float = _PLATEAU_EDGE_BUFFER_RAD) -> Callable[[np.ndarray], np.ndarray]:
    """A `ContinentalRelief(plateau_mask=...)`-compatible callable: 1.0 anywhere inside one of
    `plateaus`' own polygons, falling smoothly to 0.0 by `edge_buffer_rad` past its boundary.
    An empty `plateaus` gives an always-zero mask (no plateaus at all), rather than an error."""
    if not plateaus:
        return lambda xyz: np.zeros(len(np.asarray(xyz)), dtype=float)
    polygons = [_ring_to_xyz(ring) for ring in plateaus.values()]
    boundary_pts = np.concatenate(polygons)
    tree = cKDTree(boundary_pts)
    chord_buffer = _chord_for_angle(edge_buffer_rad)

    def mask(xyz: np.ndarray) -> np.ndarray:
        xyz = np.asarray(xyz, dtype=float)
        inside = geometry.points_in_any_spherical_polygon(xyz, polygons)
        dist, _ = tree.query(xyz, workers=-1)
        outside_falloff = _smoothstep(chord_buffer, 0.0, dist)
        return np.where(inside, 1.0, outside_falloff)

    return mask


# --- Earth: real major mountain belts (same control points authored for the Earth sketch's
# own mountain ink -- see the frontend's premadeWorlds.ts generation script) and real major
# plateaus (new here; simple, coarse polygons -- these are "roughly where a plateau is", not
# a precise cadastral boundary).
EARTH_BELTS: dict[str, Ring] = {
    "Rockies+Sierra Madre": [
        (-105, 34), (-106, 40), (-113, 44), (-114, 48), (-120, 52), (-124, 56), (-138, 60), (-151, 63),
    ],
    "Appalachians": [(-86, 33), (-83, 36), (-79, 39), (-75, 42), (-70, 46)],
    "Andes": [
        (-72, 8), (-77, 2), (-78, -5), (-76, -12), (-71, -17), (-69, -22), (-70, -28), (-70, -34), (-71, -40),
        (-72, -46), (-72, -51), (-70, -55),
    ],
    "Atlas": [(-9, 30), (-6, 32), (-1, 34), (4, 36)],
    "Alps+Carpathians": [(7, 44), (9, 46), (13, 47), (16, 48), (20, 48), (24, 47.5), (26, 46)],
    "East African Rift/Ethiopian Highlands": [(38, 12), (36, 6), (35, 0), (34, -6), (33, -10), (34, -14), (35, -18)],
    "Urals": [(59, 51), (60, 56), (62, 61), (65, 66)],
    "Himalaya+Karakoram+Tian Shan": [
        (73, 36), (77, 35), (81, 32), (85, 28), (90, 27), (95, 28), (99, 30), (103, 33), (86, 40), (80, 42), (75, 42),
    ],
    "Scandinavian Mountains": [(7, 58), (9, 61), (14, 65), (18, 68)],
    "Great Dividing Range": [(145, -12), (147, -18), (148, -25), (150, -30), (149, -34), (147, -37)],
    "Zagros": [(46, 34), (49, 31), (52, 28), (56, 27)],
}

EARTH_PLATEAUS: dict[str, Ring] = {
    "Tibetan Plateau": [(78, 36), (95, 36), (98, 32), (90, 28), (82, 29), (78, 32)],
    "Colorado Plateau": [(-112, 39), (-107, 39), (-106, 35), (-111, 34), (-113, 36)],
    "Deccan Plateau": [(74, 21), (79, 21), (80, 15), (76, 11), (73, 15)],
    "Brazilian Highlands": [(-50, -14), (-42, -16), (-40, -22), (-47, -25), (-52, -20)],
    "Patagonian Plateau": [(-72, -45), (-67, -46), (-68, -51), (-72, -51)],
}

# Every Earth belt/plateau's "which Pangaea piece does it move with" -- matches
# real_plates.PANGAEA_PLATE_GROUPS' own grouping (a belt/plateau on a continent that gets
# rigidly repositioned for Pangaea moves with it; one that stays put -- Africa/Eurasia/India
# in this reconstruction -- uses "identity").
_EARTH_REGION_GROUPS = {
    "Rockies+Sierra Madre": "americas",
    "Appalachians": "americas",
    "Andes": "americas",
    "Atlas": "identity",
    "Alps+Carpathians": "identity",
    "East African Rift/Ethiopian Highlands": "identity",
    "Urals": "identity",
    "Himalaya+Karakoram+Tian Shan": "identity",
    "Scandinavian Mountains": "identity",
    "Great Dividing Range": "australia",
    "Zagros": "identity",
    "Tibetan Plateau": "identity",
    "Colorado Plateau": "americas",
    "Deccan Plateau": "identity",
    "Brazilian Highlands": "americas",
    "Patagonian Plateau": "americas",
}


def _transform_ring(ring: Ring, group: str) -> Ring:
    lon = np.array([p[0] for p in ring], dtype=float)
    lat = np.array([p[1] for p in ring], dtype=float)
    lon2, lat2 = real_plates._rotate_lonlat(lon, lat, real_plates.PANGAEA_TRANSFORMS[group])
    return list(zip(lon2.tolist(), lat2.tolist()))


def pangaea_belts() -> dict[str, Ring]:
    """Earth's own belts, each carried along by the same rigid transform its continent gets
    for the Pangaea coastline reconstruction (see real_plates.PANGAEA_TRANSFORMS), plus one
    detail real geology supports well: the Appalachian-Caledonian-Variscan-Mauritanide belt is
    one continuous orogen in the Pangaea configuration -- North America, Europe and NW Africa
    were sutured together along exactly this line -- represented here as a single belt bridging
    the (already-transformed) Appalachians into NW Africa's Atlas position."""
    out = {name: _transform_ring(ring, _EARTH_REGION_GROUPS[name]) for name, ring in EARTH_BELTS.items()}
    out["Central Pangean Suture"] = out["Appalachians"] + out["Atlas"]
    return out


def pangaea_plateaus() -> dict[str, Ring]:
    return {name: _transform_ring(ring, _EARTH_REGION_GROUPS[name]) for name, ring in EARTH_PLATEAUS.items()}


# --- Z&D: hand-placed belts/regions gathered from public wiki descriptions of the setting's
# geography (combining the general layout several community maps agree on, not tracing any
# one of them) -- positioned against the same coordinates the frontend's premadeWorlds.ts
# Z&D sketch-building script uses for Westeros/Essos/Sothoryos.
GOT_BELTS: dict[str, Ring] = {
    "Mountains of the Moon": [(-33, 10), (-31, 4), (-33, 0)],
    "Red Mountains": [(-40, -1), (-36, -3), (-32.5, -6)],
    "Bone Mountains": [(56, -10), (60, -18), (58, -26)],
    "Grey Cliffs": [(-45.5, 20), (-46.5, 14), (-44, 6)],
    "Fourteen Flames": [(50, -22), (54, -24), (52, -27)],
    "Sothoryan Interior Mountains": [(10, -40), (28, -42), (44, -39)],
}
