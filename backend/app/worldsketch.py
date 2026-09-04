"""Turn a hand-drawn or loaded coastline image into generation inputs for
`lithosphere_plate.generate_plates` -- the "Human-made" Generate World tab's backing pipeline
(see `main.py`'s `/world/generate` `sketch` field).

The convention, shared by both the in-app drawing tool and an image loaded from another
program: a mostly-white canvas, coastline drawn as a near-black outline (strokes, not filled
regions), with two optional extra ink colors -- a blue-ish stroke for rivers, a brown/orange
stroke for mountains. Land vs. sea is therefore never painted directly; it's resolved by
flood-filling the white area on either side of the coastline, seeded from the four map corners
as presumed ocean. This is a deliberately simple heuristic (real hand-drawn sketches are rough,
not perfectly closed loops), so:

- small pen gaps in the coastline are closed by a light dilation before flood-filling (a gap
  wider than `_COAST_GAP_CLOSE_ITERS` iterations at the working grid resolution will still leak
  land into ocean or vice versa -- this is the tool's tolerance for a rough sketch, not a bug);
- a landmass or ocean that happens to wrap the antimeridian is kept as one connected region
  (`_label_wrapped` pads the grid across the longitude seam before labeling -- latitude/the
  poles are not periodic and are left unpadded);
- if a drawn "continent" turns out to fully enclose the globe (no ocean at all reachable from
  any corner), that's taken at face value -- an all-land world.

`parse_sketch_image` is the module's one real entry point; `sketch_plate_sites` is the other
(consumed by `lithosphere_plate.generate_plates`, see that function's own `sketch` param) --
everything else here is a helper for those two.
"""

from __future__ import annotations

import io
from dataclasses import dataclass

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.cluster.vq import kmeans2

from . import geometry

# The working equirectangular grid every sketch is resampled onto, independent of the source
# canvas/image resolution -- 0.5 degree cells, plenty of fidelity for a mouse/finger sketch and
# cheap to flood-fill/label. Row 0 is the north pole, column 0 is the antimeridian (-180),
# matching main.py's sample_at's own row/col convention.
SKETCH_GRID_W = 720
SKETCH_GRID_H = 360

# A safety cap on the source image before per-pixel color classification -- bounds classify
# cost/memory for an unexpectedly huge upload (e.g. a phone photo). Classification happens
# before downsampling to the working grid (see parse_sketch_image), so thin strokes in a
# larger-than-this source are downscaled here too, same as any other oversized upload.
_MAX_INPUT_DIM = 2200

# Dilation iterations (3x3 structuring element) applied to the coastline mask, at working-grid
# resolution, before flood-filling -- closes small pen gaps so a rough, not-quite-closed
# coastline doesn't leak ocean into land or vice versa. Each iteration closes the coastline in
# by 1 cell (~55 km) on each side of a gap, so `n` iterations close gaps up to `2*n` cells wide.
_COAST_GAP_CLOSE_ITERS = 2

# Color classification thresholds (0-255 channels). A near-black/gray stroke (low chroma, dark)
# is coastline; a stroke whose blue or red channel clearly leads the others (and isn't near-
# white) is a river/mountain stroke respectively. See module docstring for the ink convention.
_COAST_BRIGHTNESS_MAX = 140.0
_GRAY_CHROMA_MAX = 30
_HUE_MARGIN = 40
_INK_BRIGHTNESS_MAX = 235.0

# Frontend pen colors (SketchEditor.tsx) -- kept here too so a test/tooling script can draw a
# synthetic sketch with colors guaranteed to classify correctly, without duplicating the exact
# hex values by hand.
COAST_COLOR = "#000000"
RIVER_COLOR = "#1e6fd9"
MOUNTAIN_COLOR = "#8b4513"


def _classify_colors(rgb: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """`rgb` (H, W, 3) uint8 -> (coast, river, mountain) boolean masks -- see the threshold
    constants above and the module docstring's ink convention."""
    r = rgb[..., 0].astype(np.int16)
    g = rgb[..., 1].astype(np.int16)
    b = rgb[..., 2].astype(np.int16)
    brightness = (r + g + b) / 3.0
    chroma = rgb.max(axis=-1).astype(np.int16) - rgb.min(axis=-1).astype(np.int16)
    coast = (brightness < _COAST_BRIGHTNESS_MAX) & (chroma < _GRAY_CHROMA_MAX)
    river = (b - np.maximum(r, g) > _HUE_MARGIN) & (brightness < _INK_BRIGHTNESS_MAX) & ~coast
    mountain = (r - b > _HUE_MARGIN) & (r >= g) & (g >= b - 10) & (brightness < _INK_BRIGHTNESS_MAX) & ~coast
    return coast, river, mountain


def _max_downsample(mask: np.ndarray, out_h: int, out_w: int) -> np.ndarray:
    """Boolean `mask` (H, W) resampled down to (out_h, out_w) by max-pooling (any True pixel
    in a cell's source neighborhood keeps that cell True) -- preserves thin strokes far better
    than an averaging resize would, which is what a plain PIL resize of a 1-2px-wide line onto
    a much smaller grid would otherwise mostly erase."""
    h, w = mask.shape
    if (h, w) == (out_h, out_w):
        return mask.copy()
    fh = max(1, int(np.ceil(h / out_h)))
    fw = max(1, int(np.ceil(w / out_w)))
    filtered = ndimage.maximum_filter(mask, size=(fh, fw), mode="nearest")
    row_idx = np.clip((np.arange(out_h) * h / out_h).astype(int), 0, h - 1)
    col_idx = np.clip((np.arange(out_w) * w / out_w).astype(int), 0, w - 1)
    return filtered[np.ix_(row_idx, col_idx)]


def _label_wrapped(mask: np.ndarray) -> tuple[np.ndarray, int]:
    """`scipy.ndimage.label` (4-connectivity) over `mask` (H, W), treating the left/right edges
    as adjacent (longitude wraps; latitude/the poles do not, so rows are left unpadded) --
    padding the array across the seam before labeling means a component that crosses the
    antimeridian is connected the same way any other adjacent pair of True cells would be, no
    separate post-hoc label-merging needed."""
    padded = np.pad(mask, ((0, 0), (1, 1)), mode="wrap")
    structure = ndimage.generate_binary_structure(2, 1)
    labels_padded, n = ndimage.label(padded, structure=structure)
    return labels_padded[:, 1:-1], n


def _flood_fill_land(coast: np.ndarray) -> np.ndarray:
    """`coast` (H, W boolean, already gap-closed) -> a land mask: every non-coast pixel
    connected (via `_label_wrapped`) to any of the four grid corners is ocean; everything else
    (including the coastline pixels themselves) is land. See module docstring for why the
    corners are the seed -- a hand-drawn or loaded world map conventionally leaves its corners
    as open ocean."""
    non_coast = ~coast
    labels, _n = _label_wrapped(non_coast)
    h, w = coast.shape
    corner_labels = {int(labels[0, 0]), int(labels[0, w - 1]), int(labels[h - 1, 0]), int(labels[h - 1, w - 1])}
    corner_labels.discard(0)  # 0 == background (coast) label, never a real ocean component
    ocean = np.isin(labels, list(corner_labels)) if corner_labels else np.zeros_like(coast)
    return ~ocean


@dataclass
class SketchMasks:
    """The parsed result of `parse_sketch_image`: land/mountain/river as boolean grids on the
    working `SKETCH_GRID_H` x `SKETCH_GRID_W` equirectangular lattice, plus world-xyz sampling
    for arbitrary points (what `lithosphere_plate.generate_plates`'s `hc_at` closures and
    `sketch_plate_sites` actually consume). `mountain`/`river` are already restricted to land
    cells (see `parse_sketch_image`) -- a stray stroke that lands in the ocean has no effect."""

    land: np.ndarray
    mountain: np.ndarray
    river: np.ndarray

    @property
    def height(self) -> int:
        return self.land.shape[0]

    @property
    def width(self) -> int:
        return self.land.shape[1]

    def _rowcol(self, xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        lat_rad, lon_rad = geometry.xyz_to_latlon(xyz)
        lat_deg = np.degrees(lat_rad)
        lon_deg = np.degrees(lon_rad)
        h, w = self.height, self.width
        row = np.clip(np.round((90.0 - lat_deg) * h / 180.0 - 0.5), 0, h - 1).astype(int)
        col = np.mod(np.round((lon_deg + 180.0) * w / 360.0 - 0.5), w).astype(int)
        return row, col

    def sample_land(self, xyz: np.ndarray) -> np.ndarray:
        row, col = self._rowcol(xyz)
        return self.land[row, col]

    def sample_mountain(self, xyz: np.ndarray) -> np.ndarray:
        row, col = self._rowcol(xyz)
        return self.mountain[row, col]

    def sample_river(self, xyz: np.ndarray) -> np.ndarray:
        row, col = self._rowcol(xyz)
        return self.river[row, col]

    def cell_centers_xyz(self) -> np.ndarray:
        """(H, W, 3) world unit vectors for every grid-cell center. Cached on first use."""
        cached = getattr(self, "_cell_xyz", None)
        if cached is None:
            h, w = self.height, self.width
            lat_deg = 90.0 - (np.arange(h) + 0.5) * 180.0 / h
            lon_deg = -180.0 + (np.arange(w) + 0.5) * 360.0 / w
            lat_grid, lon_grid = np.meshgrid(lat_deg, lon_deg, indexing="ij")
            cached = geometry.latlon_to_xyz(np.radians(lat_grid), np.radians(lon_grid))
            self._cell_xyz = cached
        return cached

    def cell_weight(self) -> np.ndarray:
        """(H, W) per-cell area weight (~cos(lat)), for turning a cell count into a share of
        actual sphere area -- an equirectangular grid's cells shrink toward the poles."""
        cached = getattr(self, "_cell_weight", None)
        if cached is None:
            lat_deg = 90.0 - (np.arange(self.height) + 0.5) * 180.0 / self.height
            cached = np.repeat(np.cos(np.radians(lat_deg))[:, None], self.width, axis=1)
            self._cell_weight = cached
        return cached


def parse_sketch_image(image_bytes: bytes) -> SketchMasks:
    """Decode `image_bytes` (a PNG, from either the in-app drawing tool or a file the user
    picked) into `SketchMasks`. Raises `ValueError` if the bytes aren't a decodable image --
    the caller (main.py) turns that into a 400. Any source size/aspect ratio is accepted and
    stretched onto the working `SKETCH_GRID_W` x `SKETCH_GRID_H` grid, matching the assumption
    the in-app drawing tool's own fixed-aspect canvas already guarantees."""
    try:
        img = Image.open(io.BytesIO(image_bytes)).convert("RGB")
    except Exception as exc:
        raise ValueError("could not decode sketch image") from exc
    if max(img.size) > _MAX_INPUT_DIM:
        img = img.copy()
        img.thumbnail((_MAX_INPUT_DIM, _MAX_INPUT_DIM), Image.LANCZOS)
    rgb = np.asarray(img)

    coast, river, mountain = _classify_colors(rgb)
    coast = _max_downsample(coast, SKETCH_GRID_H, SKETCH_GRID_W)
    river = _max_downsample(river, SKETCH_GRID_H, SKETCH_GRID_W)
    mountain = _max_downsample(mountain, SKETCH_GRID_H, SKETCH_GRID_W)

    coast_closed = ndimage.binary_dilation(coast, structure=np.ones((3, 3), dtype=bool), iterations=_COAST_GAP_CLOSE_ITERS)
    land = _flood_fill_land(coast_closed)
    # Rivers/mountains only ever matter on land (see hc_at's own gate in lithosphere_plate.py)
    # -- restricting them here keeps SketchMasks' own fields consistent with what applies.
    return SketchMasks(land=land, mountain=mountain & land, river=river & land)


@dataclass
class Landmass:
    """One connected landmass from a parsed sketch -- its member grid cells' world positions
    and a total area weight (see `SketchMasks.cell_weight`), used by `sketch_plate_sites` to
    decide how many continental plates to split it across."""

    xyz: np.ndarray  # (n, 3)
    area: float


def landmasses(masks: SketchMasks) -> list[Landmass]:
    """Every connected component of `masks.land` (wraparound-aware, see `_label_wrapped`),
    largest first."""
    labels, n = _label_wrapped(masks.land)
    if n == 0:
        return []
    cell_xyz = masks.cell_centers_xyz()
    cell_weight = masks.cell_weight()
    result = []
    for label_id in range(1, n + 1):
        member = labels == label_id
        if not member.any():
            continue
        result.append(Landmass(xyz=cell_xyz[member], area=float(cell_weight[member].sum())))
    result.sort(key=lambda lm: -lm.area)
    return result


def _distribute_counts(total: int, weights: list[float]) -> list[int]:
    """`total` plate slots divided across `len(weights)` landmasses proportional to `weights`
    (largest-remainder rounding), guaranteeing every landmass at least 1 plate when there are
    enough slots to go around -- so a small island never silently vanishes into whichever
    oceanic plate's Voronoi cell happens to cover it. When `total < len(weights)`, only the
    `total` largest landmasses get a plate (the rest fall inside a neighboring plate's
    territory, an acceptable degrade for "more islands than requested continental plates")."""
    n = len(weights)
    if n == 0:
        return []
    total = max(0, total)
    if total < n:
        # Drawn land should never end up with zero continental plates just because a low
        # continental_fraction rounded the target down to 0 -- always give at least the
        # largest landmass one, even if `total` itself was 0.
        order = np.argsort(weights)[::-1]
        counts = [0] * n
        for i in order[: max(total, 1)]:
            counts[int(i)] = 1
        return counts
    counts = [1] * n
    remaining = total - n
    if remaining > 0:
        total_weight = sum(weights) or 1.0
        raw_extra = [remaining * w / total_weight for w in weights]
        extra = [int(np.floor(x)) for x in raw_extra]
        leftover = remaining - sum(extra)
        order = sorted(range(n), key=lambda i: raw_extra[i] - extra[i], reverse=True)
        for i in order[:leftover]:
            extra[i] += 1
        counts = [c + e for c, e in zip(counts, extra)]
    return counts


def _kmeans_sites(xyz: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    """`k` seed sites (unit vectors) covering `xyz` (n, 3) -- a weighted centroid for `k <= 1`,
    else `scipy.cluster.vq.kmeans2` cluster centers (re-normalized back onto the sphere; the
    landmass's own angular extent is always small enough for Euclidean k-means on its xyz
    points to be a fine stand-in for a true spherical clustering)."""
    if len(xyz) == 0:
        return np.zeros((0, 3))
    k = max(1, min(k, len(xyz)))
    if k == 1:
        return geometry.normalize(xyz.mean(axis=0, keepdims=True))
    seed_val = int(rng.integers(0, 2**31 - 1))
    centers, _ = kmeans2(xyz, k, minit="++", seed=seed_val)
    return geometry.normalize(centers)


def _farthest_point_sites(candidates_xyz: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    """Greedy farthest-point sampling of `k` sites from `candidates_xyz` (n, 3, unit vectors)
    -- repeatedly picks whichever remaining candidate is farthest (in angular terms) from every
    site already chosen, so oceanic plate sites spread out over open ocean instead of
    clustering. Angular distance is monotonic in dot product, so "farthest" is just "smallest
    max similarity to anything chosen so far" -- no arccos needed."""
    if k <= 0 or len(candidates_xyz) == 0:
        return np.zeros((0, 3))
    k = min(k, len(candidates_xyz))
    start = int(rng.integers(0, len(candidates_xyz)))
    chosen = [start]
    best_similarity = candidates_xyz @ candidates_xyz[start]
    for _ in range(1, k):
        next_idx = int(np.argmin(best_similarity))
        chosen.append(next_idx)
        best_similarity = np.maximum(best_similarity, candidates_xyz @ candidates_xyz[next_idx])
    return candidates_xyz[chosen]


def sketch_plate_sites(
    masks: SketchMasks, num_plates: int, num_continents: int, rng: np.random.Generator
) -> tuple[np.ndarray, list[str]]:
    """The "plate boundaries that fit the map" step: `num_plates` primary seed sites (unit
    vectors) plus a parallel crust-type list, biased so plate boundaries fall in open ocean and
    across drawn landmasses rather than scattered uniformly at random (see
    `lithosphere_plate.generate_plates`'s `sketch` param, which calls this).

    `num_continents` plate slots are divided across `landmasses(masks)` by area share
    (`_distribute_counts`), each landmass's share placed via `_kmeans_sites` -- so a large
    hand-drawn continent assigned several plates gets pre-partitioned into plausible
    sub-regions before Voronoi tiling ever runs, the same way Earth's own continents straddle
    more than one plate. The remaining slots go to oceanic sites, farthest-point-sampled over
    every non-land grid cell so they spread out across open ocean.

    The returned site count may exceed `num_plates` (every landmass is guaranteed at least one
    plate even if that overshoots the target, see `_distribute_counts`) -- the caller
    (`generate_plates`) uses `len(site_xyz)` as the real final plate count, the same tolerance
    it already has for bumping plate count up to satisfy `MIN_OCEANIC_PLATES`."""
    lands = landmasses(masks)
    num_continents = max(0, min(num_continents, num_plates)) if lands else 0
    counts = _distribute_counts(num_continents, [lm.area for lm in lands]) if lands else []

    continental_chunks = [
        _kmeans_sites(lm.xyz, k, rng) for lm, k in zip(lands, counts) if k > 0
    ]
    continental_sites = np.concatenate(continental_chunks, axis=0) if continental_chunks else np.zeros((0, 3))

    num_oceanic = max(num_plates - len(continental_sites), 0)
    ocean_xyz = masks.cell_centers_xyz()[~masks.land]
    oceanic_sites = _farthest_point_sites(ocean_xyz, num_oceanic, rng)

    site_xyz = np.concatenate([continental_sites, oceanic_sites], axis=0)
    crust_types = ["continental"] * len(continental_sites) + ["oceanic"] * len(oceanic_sites)
    return site_xyz, crust_types
