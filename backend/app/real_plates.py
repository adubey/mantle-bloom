"""Real-world plate geometry and motion, for the Premade-worlds "predetermined plates" and
"mantle centers worked backwards from known plate motion" features (see App.tsx's Generate
World dialog) -- Earth uses this data directly; Pangaea reuses the same plates rigidly
rotated by the identical per-continent transform the coastline reconstruction already applies
(see the frontend's premadeWorlds.ts generation script and PANGAEA_PLATE_TRANSFORMS below).

Data source: NNR-MORVEL56 (Argus, D. F., R. G. Gordon, and C. DeMets, 2011, "Geologically
current motion of 56 plates relative to the no-net-rotation reference frame", Geochemistry,
Geophysics, Geosystems), as bundled by the `platemotion` PyPI package -- both the absolute
Euler pole table and the matching plate boundary outlines, trimmed to the 16 plates that
matter at this app's usual plate-count scale (see data/major_plates.json). Ridge/trench
candidate positions for the mantle-center fit come from PB2002 (Bird, 2002), via the
`fraxen/tectonicplates` GeoJSON export, clustered down from individual boundary steps to a
handful of representative points per major ridge/trench system (see data/ridge_trench_centers.json).
`data/plate_owner_grid.npz` (which of those same 16 polygons contains each cell of a 0.25-
degree grid) is baked from that same boundary data. `data/real_land_mask.npz` (real coastline,
Natural Earth `ne_10m_land` polygons rasterized onto a 0.05-degree grid) is the one exception
to "boundary/motion data" this module otherwise deals in -- see the "Exact Earth coastline/
plates" section below for why coastline needed its own real-data source, not just plates. All
four files are baked in at generation time by one-off scripts, not fetched at runtime.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree

from . import geometry, mantle

# Plain __file__-relative path (matching desktop.py's own convention) rather than
# importlib.resources -- the PyInstaller onedir build's `datas` list (see
# bin/mantle-bloom.spec) copies this directory to the same relative location, but doesn't
# make it importable as a package resource.
_DATA_DIR = Path(__file__).resolve().parent / "data"

# Degrees/Ma -> radians/year, matching mantle.py's own angular-rate convention.
_DEG_PER_MYR_TO_RAD_PER_YR = (np.pi / 180.0) / 1e6

# Fixed falloff for every ridge/trench candidate center (radians) -- picked to comfortably
# span the ~22 degree cluster cell (see build_ridge_trench.py) each candidate represents,
# in the same spirit as mantle.generate_convection_centers' own 0.6-1.4 rad random range.
_CANDIDATE_FALLOFF_RAD = 0.55

# A candidate ridge/trench center's initial strength guess, before the least-squares fit
# below refines it -- sign only matters here (positive = upwelling/ridge, negative =
# downwelling/trench), the fit finds the actual magnitude.
_INITIAL_STRENGTH = mantle.MANTLE_FLOW_REFERENCE_RATE

# Sample points per plate when building the mantle-center fit's target (point, velocity)
# pairs -- a handful spread across each plate's own polygon, not the full sketch grid (the
# fit only needs to be representative, not exhaustive; see fit_mantle_centers).
_FIT_SAMPLES_PER_PLATE = 24


@dataclass
class RealPlate:
    name: str
    code: str
    boundary_xyz: np.ndarray  # (n, 3) unit vectors, ordered ring
    omega: np.ndarray  # (3,) rad/yr angular velocity (axis = Euler pole, magnitude = rate)


def _pole_to_omega(lat_deg: float, lon_deg: float, rate_deg_per_myr: float) -> np.ndarray:
    axis = geometry.latlon_to_xyz(np.radians(np.array([lat_deg])), np.radians(np.array([lon_deg])))[0]
    return axis * (rate_deg_per_myr * _DEG_PER_MYR_TO_RAD_PER_YR)


# `geometry.points_in_spherical_polygon` builds an (n_points, n_vertices) intermediate for
# every call -- fine at the "a plate's own few hundred near-boundary nodes vs. a neighbour's
# handful-to-~100-vertex outline" scale it was built for, but the real NNR-MORVEL56 outlines
# run to 1000+ vertices, and real_plate_sites/_target_samples both test the *entire* sketch
# grid (259200 cells) against every plate -- (259200, 1200)-shaped intermediates measured at
# well over a minute per large plate. Decimating to this many vertices keeps the outline's
# overall shape (plenty for classifying which 0.5-degree grid cell belongs to which plate)
# while cutting that cost by an order of magnitude.
_MAX_BOUNDARY_VERTICES = 120


def _decimate_ring(xyz: np.ndarray, max_vertices: int) -> np.ndarray:
    if len(xyz) <= max_vertices:
        return xyz
    stride = math.ceil(len(xyz) / max_vertices)
    return xyz[::stride]


def load_major_plates() -> list[RealPlate]:
    """The 16 major real plates (NNR-MORVEL56), boundary + absolute Euler pole. Cheap enough
    (a few hundred KB of JSON, parsed once) to not bother caching beyond Python's own import
    caching of whatever calls this once per generation."""
    raw = json.loads((_DATA_DIR / "major_plates.json").read_text())
    plates = []
    for p in raw["plates"]:
        boundary_lonlat = np.array(p["boundary"], dtype=float)
        boundary_xyz = geometry.latlon_to_xyz(np.radians(boundary_lonlat[:, 1]), np.radians(boundary_lonlat[:, 0]))
        boundary_xyz = _decimate_ring(boundary_xyz, _MAX_BOUNDARY_VERTICES)
        pole = p["euler_pole"]
        omega = _pole_to_omega(pole["lat"], pole["lon"], pole["rate_deg_per_myr"])
        plates.append(RealPlate(name=p["name"], code=p["code"], boundary_xyz=boundary_xyz, omega=omega))
    return plates


# --- Exact Earth coastline/plates ---------------------------------------------------------
#
# "Present-day Earth" no longer decides land/sea from the drawn/traced sketch PNG at all (see
# worldsketch.py's module docstring for why that path is lossy: a hand-tolerant coastline
# stroke needs gap-closing dilation to flood-fill cleanly, and that same dilation swallows any
# real strait/sea narrower than it -- the Mediterranean via Gibraltar, the Korea Strait,
# Torres Strait all measured as casualties of the *same* dilation amount that keeps ordinary
# coastlines closed, at every resolution the baked sketch PNG offers). `real_land_mask.npz`
# instead is a direct rasterization (see the one-off bake, not checked in -- Natural Earth's
# `ne_10m_land` polygons filled straight onto a 0.05-degree grid, no outline-then-flood-fill
# step to lose a gap in) of the actual coastline, at 10x the old sketch grid's resolution.
# Gibraltar itself is still narrower than even this grid resolves (confirmed by inspection:
# the source polygons touch with zero gap at Gibraltar, a generalization of the underlying
# cartographic data, not a bug in this bake) -- see hydrology.py's own Gibraltar note for how
# the Mediterranean still ends up counted as ocean despite that.
#
# Plate *membership* similarly no longer comes from Voronoi tiling over a handful of sites
# (`real_plate_sites` below, still used for Pangaea/GoT) -- `plate_owner_grid.npz` (baked the
# same way, from the unrotated `major_plates.json` boundaries) directly assigns every 0.25-
# degree cell to whichever of the 16 real plates' polygon contains it (covering ~93% of the
# globe; the rest is real microplates this app's simplified 16-plate set omits entirely).
# `build_exact_earth_plates` below turns that into actual `LithospherePlate` territories: one
# oceanic sub-plate per real plate (its crust type is "monolithic" in this engine, see
# `generate_plates`'s own docstring, so a real plate's onshore/offshore crust still has to be
# two separate LithospherePlate objects here, same as the old Voronoi path already did) plus
# one continental sub-plate per land component *above `_MIN_LAND_COMPONENT_KM2`* within it --
# smaller than that is folded into the surrounding ocean rather than getting its own plate.
# That floor exists only because a `LithospherePlate` is a real, persistent per-step cost
# (its own KD-tree, its own deform() pass, ...) and this app has never been run with more than
# a few dozen of them (32 in its own stress tests) -- every real island get its own exact
# plate would mean several hundred, untested territory for the whole engine, not just
# generation. `_MIN_LAND_COMPONENT_KM2` = 50,000 keeps the total near that tested range while
# still covering every landmass big enough to plausibly matter (all of Indonesia's/Japan's/
# the Mediterranean's/the Caribbean's major islands clear it) -- a smaller real island
# (Cyprus, Crete, Puerto Rico, ~8-10k km2 each) is the actual, deliberate cost of that
# tradeoff: it renders as open ocean, not a bug, and the floor is one constant to lower if
# that trade is ever worth revisiting against the performance it costs.
_EXACT_GRID_H = 720
_EXACT_GRID_W = 1440
_MIN_LAND_COMPONENT_KM2 = 50_000.0
_EARTH_RADIUS_KM = 6371.0


def _exact_grid_rowcol(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Same row/col convention as `worldsketch.SketchMasks._rowcol`, standalone here so
    `build_exact_earth_plates`'s label grid can be looked up without constructing a
    `SketchMasks` for it (its H/W differ from the sketch's own working grid)."""
    lat_rad, lon_rad = geometry.xyz_to_latlon(xyz)
    lat_deg = np.degrees(lat_rad)
    lon_deg = np.degrees(lon_rad)
    row = np.clip(np.round((90.0 - lat_deg) * _EXACT_GRID_H / 180.0 - 0.5), 0, _EXACT_GRID_H - 1).astype(int)
    col = np.mod(np.round((lon_deg + 180.0) * _EXACT_GRID_W / 360.0 - 0.5), _EXACT_GRID_W).astype(int)
    return row, col


_real_land_mask_cache: np.ndarray | None = None


def load_real_land_mask() -> np.ndarray:
    """The fine (0.05-degree) real-coastline raster -- see the module comment above. Cached
    at module scope since it's ~26M booleans decompressed and every call this generation (or
    load, if a saved exact-Earth world is later reopened) wants the identical array."""
    global _real_land_mask_cache
    if _real_land_mask_cache is None:
        with np.load(_DATA_DIR / "real_land_mask.npz") as data:
            _real_land_mask_cache = data["land"]
    return _real_land_mask_cache


_plate_owner_grid_cache: tuple[np.ndarray, list[str]] | None = None


def load_plate_owner_grid() -> tuple[np.ndarray, list[str]]:
    """`(owner, codes)` -- `owner` ((_EXACT_GRID_H, _EXACT_GRID_W) int8, -1 or an index into
    `codes`) is which real plate's polygon (unrotated `major_plates.json` boundaries) contains
    each cell, baked the same way `real_land_mask.npz` was. Cached like that mask."""
    global _plate_owner_grid_cache
    if _plate_owner_grid_cache is None:
        with np.load(_DATA_DIR / "plate_owner_grid.npz") as data:
            _plate_owner_grid_cache = (data["owner"], [str(c) for c in data["codes"]])
    return _plate_owner_grid_cache


@dataclass
class ExactEarthPlates:
    """`build_exact_earth_plates`'s result: `label_grid` ((_EXACT_GRID_H, _EXACT_GRID_W)
    int32, -1 or an index into `crust_types`/`frame_xyz`) assigns every cell to one final
    LithospherePlate-to-be; `crust_types[i]`/`frame_xyz[i]` (a representative unit vector,
    that sub-plate's own member-cell centroid) describe sub-plate `i`."""

    label_grid: np.ndarray
    crust_types: list[str]
    frame_xyz: list[np.ndarray]

    def is_owned(self, sub_plate_id: int, world_pts: np.ndarray) -> np.ndarray:
        row, col = _exact_grid_rowcol(world_pts)
        return self.label_grid[row, col] == sub_plate_id


_exact_earth_plates_cache: ExactEarthPlates | None = None


def build_exact_earth_plates() -> ExactEarthPlates:
    """The real-plate-grounded partition "Present-day Earth" builds its `LithospherePlate`s
    from -- see the module comment above for the oceanic/continental split and why
    `_MIN_LAND_COMPONENT_KM2` exists. Deterministic (no seed/rng: both input rasters are
    static data), so cached at module scope like the rasters it's built from."""
    global _exact_earth_plates_cache
    if _exact_earth_plates_cache is not None:
        return _exact_earth_plates_cache

    from . import worldsketch  # local import: real_plates.py is the newer, more special-purpose module

    owner, codes = load_plate_owner_grid()
    fine_land = load_real_land_mask()
    fh, fw = fine_land.shape[0] // _EXACT_GRID_H, fine_land.shape[1] // _EXACT_GRID_W
    coarse_land = fine_land[: fh * _EXACT_GRID_H, : fw * _EXACT_GRID_W].reshape(_EXACT_GRID_H, fh, _EXACT_GRID_W, fw).max(axis=(1, 3))

    lat_deg = 90.0 - (np.arange(_EXACT_GRID_H) + 0.5) * 180.0 / _EXACT_GRID_H
    lon_deg = -180.0 + (np.arange(_EXACT_GRID_W) + 0.5) * 360.0 / _EXACT_GRID_W
    lat_grid, lon_grid = np.meshgrid(np.radians(lat_deg), np.radians(lon_deg), indexing="ij")
    cell_xyz = geometry.latlon_to_xyz(lat_grid, lon_grid)
    cell_km2 = np.cos(lat_grid) * (np.radians(180.0 / _EXACT_GRID_H) * np.radians(360.0 / _EXACT_GRID_W)) * _EARTH_RADIUS_KM**2

    label_grid = np.full((_EXACT_GRID_H, _EXACT_GRID_W), -1, dtype=np.int32)
    crust_types: list[str] = []
    frame_xyz: list[np.ndarray] = []

    def add_sub_plate(mask: np.ndarray, crust_type: str) -> None:
        label_grid[mask] = len(crust_types)
        crust_types.append(crust_type)
        frame_xyz.append(geometry.normalize(cell_xyz[mask].mean(axis=0)))

    # The 16 real plates: exactly one oceanic sub-plate per plate (its own onshore/offshore
    # crust already differs in reality; splitting the *sea* side further by its own internal
    # connectivity would multiply plate count for no coastline-fidelity benefit, since an
    # oceanic sub-plate's hc_at never reads the land mask either way -- see the module
    # comment), plus one continental sub-plate per land component clearing the area floor.
    for i in range(len(codes)):
        region_mask = owner == i
        if not np.any(region_mask):
            continue
        land_mask = region_mask & coarse_land
        labeled, n_components = worldsketch._label_wrapped(land_mask)
        kept_land = np.zeros_like(land_mask)
        for comp in range(1, n_components + 1):
            comp_mask = labeled == comp
            if cell_km2[comp_mask].sum() >= _MIN_LAND_COMPONENT_KM2:
                add_sub_plate(comp_mask, "continental")
                kept_land |= comp_mask
        sea_mask = region_mask & ~kept_land
        if np.any(sea_mask):
            add_sub_plate(sea_mask, "oceanic")

    # Leftover (~7% of the globe, real microplates this app's 16-plate set omits): each
    # connected patch -- land and sea together, since these are already small -- becomes its
    # own plate, continental if it holds any land clearing the same floor, else oceanic. This
    # is the "give leftover area its own small plates" choice, not folding it into whichever
    # of the 16 plates happens to be nearest.
    leftover_mask = owner == -1
    labeled, n_components = worldsketch._label_wrapped(leftover_mask)
    for comp in range(1, n_components + 1):
        comp_mask = labeled == comp
        if cell_km2[comp_mask].sum() < _MIN_LAND_COMPONENT_KM2:
            continue  # too small a scrap to bother with its own plate -- absorbed into its nearest neighbor below
        has_land = bool(coarse_land[comp_mask].any())
        add_sub_plate(comp_mask, "continental" if has_land else "oceanic")

    # Every cell above belongs to some sub-plate now except the leftover scraps just skipped
    # for being too small -- unlike a too-small *land* component within a real plate's own
    # region (which already lands in that region's own oceanic sea_mask above, never
    # unclaimed), a skipped leftover scrap has no such natural home, so it would otherwise
    # sit outside every plate's territory and never get a node at all. Absorb each into
    # whichever already-decided sub-plate has the nearest member cell -- same "too small to
    # matter, folds into its surroundings" degrade as everywhere else in this bake, just
    # applied after the fact instead of by skipping it up front.
    unclaimed = label_grid == -1
    if np.any(unclaimed) and not np.all(unclaimed):
        claimed_tree = cKDTree(cell_xyz[~unclaimed])
        claimed_labels = label_grid[~unclaimed]
        _, nearest = claimed_tree.query(cell_xyz[unclaimed])
        label_grid[unclaimed] = claimed_labels[nearest]

    _exact_earth_plates_cache = ExactEarthPlates(label_grid=label_grid, crust_types=crust_types, frame_xyz=frame_xyz)
    return _exact_earth_plates_cache


def exact_earth_masks(approx_sketch: "worldsketch.SketchMasks") -> "worldsketch.SketchMasks":
    """`approx_sketch` (the ordinary sketch-PNG parse the frontend still sends, so it still
    supplies mountain/river painting the same as before -- the reported mismatch was
    coastlines/plates, not those) with its `land` grid swapped for the exact raster (see the
    module comment above). `mountain`/`river` are nearest-upsampled from the sketch's own
    coarser grid onto the exact grid's resolution -- `SketchMasks` indexes all three fields
    with the same row/col math, so they have to share a shape."""
    from . import worldsketch  # local import: real_plates.py is the newer, more special-purpose module

    land = load_real_land_mask()
    fh = land.shape[0] // approx_sketch.height
    fw = land.shape[1] // approx_sketch.width
    mountain = np.repeat(np.repeat(approx_sketch.mountain, fh, axis=0), fw, axis=1)[: land.shape[0], : land.shape[1]]
    river = np.repeat(np.repeat(approx_sketch.river, fh, axis=0), fw, axis=1)[: land.shape[0], : land.shape[1]]
    return worldsketch.SketchMasks(land=land, mountain=mountain & land, river=river & land)


@dataclass
class RidgeTrenchCandidate:
    position: np.ndarray  # unit vector
    sign: float  # +1 upwelling (ridge/rift), -1 downwelling (trench)
    weight: float  # relative prominence (source boundary-step count in its cluster cell)


def load_ridge_trench_candidates() -> list[RidgeTrenchCandidate]:
    raw = json.loads((_DATA_DIR / "ridge_trench_centers.json").read_text())
    out = []
    for c in raw["centers"]:
        position = geometry.latlon_to_xyz(np.radians(np.array([c["lat"]])), np.radians(np.array([c["lon"]])))[0]
        out.append(RidgeTrenchCandidate(position=position, sign=float(c["sign"]), weight=float(c["weight"])))
    return out


def real_plate_sites(
    sketch,  # worldsketch.SketchMasks -- typed loosely to avoid a circular import
    real_plates: list[RealPlate],
    num_plates: int,
    num_continents: int,
    rng: np.random.Generator,
    pooled_oceanic: bool = False,
):
    """`worldsketch.sketch_plate_sites`'s own algorithm, scoped per real plate instead of per
    drawn landmass: each real plate's boundary polygon is intersected with the sketch's own
    land mask, continental plate slots are distributed across every plate's *land*
    intersection by area (`worldsketch.distribute_counts`) and placed by k-means within it
    (`worldsketch.kmeans_sites`) -- so the resulting Voronoi tiling approximates each real
    plate's actual shape (continental and oceanic portions can end up as separate
    LithospherePlate objects, same as a real plate's onshore/offshore crust already differs in
    this engine's monolithic-crust-type-per-plate model) rather than being seeded from the
    drawn coastline's connected components alone.

    `pooled_oceanic` (Pangaea's own use -- see `pangaea_real_plates`) treats oceanic slots the
    plain `worldsketch.sketch_plate_sites` way instead: farthest-point-sampled across every
    sea cell regardless of which real plate (if any) it falls in. `real_plates` there is only
    today's *continental* plates, rigidly rotated into the Pangaea fit -- there's no
    trustworthy 200 Mya analogue for the oceanic plates (Panthalassa bore no resemblance to
    today's Pacific/Nazca/Cocos/... footprints), so their sites fall back to today's ordinary
    random placement rather than pretending to a precision the data doesn't support."""
    from . import worldsketch  # local import: real_plates.py is the newer, more special-purpose module

    cell_xyz = sketch.cell_centers_xyz().reshape(-1, 3)
    cell_weight = sketch.cell_weight().reshape(-1)
    land_flat = sketch.land.reshape(-1)

    land_groups: list[tuple[np.ndarray, float]] = []
    sea_groups: list[tuple[np.ndarray, float]] = []
    for rp in real_plates:
        in_plate = geometry.points_in_spherical_polygon(cell_xyz, rp.boundary_xyz)
        land_mask = in_plate & land_flat
        if np.any(land_mask):
            land_groups.append((cell_xyz[land_mask], float(cell_weight[land_mask].sum())))
        if not pooled_oceanic:
            sea_mask = in_plate & ~land_flat
            if np.any(sea_mask):
                sea_groups.append((cell_xyz[sea_mask], float(cell_weight[sea_mask].sum())))

    num_continents = max(0, min(num_continents, num_plates))
    continental_counts = worldsketch.distribute_counts(num_continents, [a for _, a in land_groups]) if land_groups else []
    continental_chunks = [
        worldsketch.kmeans_sites(pts, k, rng) for (pts, _), k in zip(land_groups, continental_counts) if k > 0
    ]
    continental_sites = np.concatenate(continental_chunks, axis=0) if continental_chunks else np.zeros((0, 3))

    num_oceanic = max(num_plates - len(continental_sites), 0)
    if pooled_oceanic:
        ocean_xyz = cell_xyz[~land_flat]
        oceanic_sites = worldsketch.farthest_point_sites(ocean_xyz, num_oceanic, rng)
    else:
        oceanic_counts = worldsketch.distribute_counts(num_oceanic, [a for _, a in sea_groups]) if sea_groups else []
        oceanic_chunks = [
            worldsketch.farthest_point_sites(pts, k, rng) for (pts, _), k in zip(sea_groups, oceanic_counts) if k > 0
        ]
        oceanic_sites = np.concatenate(oceanic_chunks, axis=0) if oceanic_chunks else np.zeros((0, 3))

    site_xyz = np.concatenate([continental_sites, oceanic_sites], axis=0)
    crust_types = ["continental"] * len(continental_sites) + ["oceanic"] * len(oceanic_sites)

    # Rescue pass: `distribute_counts`'s "only the `num_continents` largest landmasses get a
    # plate" degrade (its own docstring calls this "an acceptable loss for more islands than
    # requested continental plates") is fine when the dropped land ends up nearest an already-
    # placed *continental* site regardless of which real plate that site nominally belongs to
    # (`hc_at` reads the sketch per-node the same way for every continental plate, see
    # `lithosphere_plate.generate_plates`'s own `sketch` docstring) -- a sliver of land near a
    # real plate boundary inherits its continental neighbor's site just like it would in
    # reality. It's not fine when the nearest site turns out to be *oceanic* instead: that
    # land renders as open ocean regardless of what was drawn, whether it's a whole isolated
    # continent with no continental neighbor to inherit from (Australia, Antarctica) or part
    # of one that a real plate's own internal boundary (e.g. Nubia/Somalia's East African Rift
    # split) happens to cut off from its nearest surviving continental site. Rather than
    # guessing which real plate's land deserves a guaranteed site up front -- by real plate
    # (over-guarantees every oceanic plate's boundary-noise sliver of coastal land too, at the
    # oceanic site budget's expense) or even by whole land group (too coarse when only part of
    # a large group ends up orphaned) -- check every drawn land cell directly against the
    # sites actually placed, and add one rescue site per *connected patch* of cells that ended
    # up oceanic-nearest (`ndimage.label`, same 8-connectivity `worldsketch._flood_fill_land`
    # uses), sized independently of real-plate identity. Added on top rather than carved out
    # of `num_oceanic`, so this can only ever add a handful of plates for genuinely orphaned
    # land, never shrink oceanic coverage -- and can only ever turn oceanic-nearest land
    # continental, never the reverse, so it can't introduce a new mismatch of its own.
    if len(site_xyz) > 0 and np.any(land_flat):
        site_tree = cKDTree(site_xyz)
        is_oceanic_site = np.array([c == "oceanic" for c in crust_types])
        land_idx = np.flatnonzero(land_flat)
        _, nearest_site = site_tree.query(cell_xyz[land_idx])
        orphaned = is_oceanic_site[nearest_site]
        if np.any(orphaned):
            orphaned_mask = np.zeros(land_flat.shape, dtype=bool)
            orphaned_mask[land_idx[orphaned]] = True
            orphaned_grid = orphaned_mask.reshape(sketch.height, sketch.width)
            labeled, n_components = ndimage.label(orphaned_grid, structure=np.ones((3, 3)))
            labeled_flat = labeled.reshape(-1)
            rescued_sites = np.concatenate(
                [worldsketch.kmeans_sites(cell_xyz[labeled_flat == comp], 1, rng) for comp in range(1, n_components + 1)],
                axis=0,
            )
            site_xyz = np.concatenate([site_xyz, rescued_sites], axis=0)
            crust_types = crust_types + ["continental"] * len(rescued_sites)

    return site_xyz, crust_types


# --- Pangaea: today's continental plates, rigidly repositioned -----------------------------
# The identical planar (pivot, rotate, translate) transforms the frontend's Pangaea coastline
# script applies to each continent's real coastline (premadeWorlds.ts's generation script,
# derived there from a 2-point rigid-registration solve -- see that script's own comments),
# reused here so a plate's boundary agrees with where its coastline actually is. A crude
# planar approximation of a true spherical rotation, not scientific paleogeographic
# reconstruction -- fine for "roughly the right shape in roughly the right place", not for
# precision. Oceanic plates (Panthalassa-era) have no real analogue and are left out; India is
# kept attached at today's position (`"identity"`, same as Eurasia) for consistency with the
# coastline, which -- inheriting Natural Earth's single fused Africa+Eurasia landmass -- never
# separately relocated it either.
PANGAEA_TRANSFORMS = {
    "identity": {"pivot": (0.0, 0.0), "angle_deg": 0.0, "translate": (0.0, 0.0)},
    "americas": {"pivot": (-53.0, 47.0), "angle_deg": 8.12261955682186, "translate": (44.0, -7.0)},
    "australia": {"pivot": (115.0, -34.0), "angle_deg": -21.43152092332241, "translate": (-40.0, 4.0)},
}
PANGAEA_PLATE_GROUPS = {
    "North America": "americas",
    "South America": "americas",
    "Africa": "identity",
    "Eurasia": "identity",
    "India": "identity",
    "Australia": "australia",
}

# India's drift before its Cenozoic collision with Asia was famously much faster than any
# plate moves today (a well-documented average around 15-18 cm/yr, several times the present
# rate) -- reflected as a rate multiplier on its own real absolute pole, direction unchanged.
# A directional stand-in, not a claimed paleo-precise value (see the module docstring).
_INDIA_PANGAEA_RATE_MULTIPLIER = 3.0


def _rotate_lonlat(lon: np.ndarray, lat: np.ndarray, transform: dict) -> tuple[np.ndarray, np.ndarray]:
    px, py = transform["pivot"]
    theta = np.radians(transform["angle_deg"])
    x, y = lon - px, lat - py
    xr = x * np.cos(theta) - y * np.sin(theta)
    yr = x * np.sin(theta) + y * np.cos(theta)
    return xr + px + transform["translate"][0], yr + py + transform["translate"][1]


def pangaea_real_plates() -> list[RealPlate]:
    """Today's continental real plates (`PANGAEA_PLATE_GROUPS`), boundary rigidly repositioned
    by `PANGAEA_TRANSFORMS`; motion direction kept as each plate's own real absolute pole
    (the sense of e.g. Atlantic opening hasn't reversed since Pangaea began breaking up), with
    India's rate boosted (`_INDIA_PANGAEA_RATE_MULTIPLIER`) for its well-documented historically
    fast drift. See the module docstring for why this is directional/approximate, not to
    NNR-MORVEL56 precision."""
    by_name = {p.name: p for p in load_major_plates()}
    out = []
    for name, group in PANGAEA_PLATE_GROUPS.items():
        rp = by_name[name]
        lat_rad, lon_rad = geometry.xyz_to_latlon(rp.boundary_xyz)
        lon2, lat2 = _rotate_lonlat(np.degrees(lon_rad), np.degrees(lat_rad), PANGAEA_TRANSFORMS[group])
        boundary_xyz = geometry.latlon_to_xyz(np.radians(lat2), np.radians(lon2))
        omega = rp.omega * (_INDIA_PANGAEA_RATE_MULTIPLIER if name == "India" else 1.0)
        out.append(RealPlate(name=rp.name, code=rp.code, boundary_xyz=boundary_xyz, omega=omega))
    return out


# How close two Pangaea plates' boundary points need to be (radians) to count as sitting on
# the same rift/suture line, for pangaea_ridge_trench_candidates' upwelling half.
_PANGAEA_SUTURE_REACH_RAD = 0.12
# Minimum angular spacing (radians) kept between clustered suture candidates -- greedy
# farthest-point thinning, so a long suture becomes a handful of representative points
# instead of every close boundary-point pair along its whole length.
_PANGAEA_SUTURE_CLUSTER_RAD = 0.35
_PANGAEA_DOWNWELLING_CANDIDATES = 20
_PANGAEA_INDIA_RIFT_CANDIDATES = 6


def _cluster_positions(positions: np.ndarray, min_spacing_rad: float) -> np.ndarray:
    """Greedy farthest-point thinning: keep the first position, then only keep each next one
    if it's more than `min_spacing_rad` from every position kept so far."""
    if len(positions) == 0:
        return positions
    kept = [positions[0]]
    for p in positions[1:]:
        if all(np.arccos(np.clip(p @ k, -1.0, 1.0)) > min_spacing_rad for k in kept):
            kept.append(p)
    return np.array(kept)


def pangaea_ridge_trench_candidates(sketch, rng: np.random.Generator) -> list[RidgeTrenchCandidate]:
    """Candidate mantle-center positions for Pangaea: there's no real 200 Mya ridge/trench
    map to draw on (unlike Earth's own PB2002-derived `load_ridge_trench_candidates`), so this
    is a directional stand-in, not real data. Upwelling candidates sit where two of the
    repositioned continental plates' own boundaries now run close together, clustered down to
    a handful of representative points per suture (`_cluster_positions`) -- the Americas/
    Africa-Eurasia suture is exactly the Central Atlantic rift that historically opened first,
    so marking it as upwelling is the one part of this that *is* well-grounded. A few more sit
    directly along India's own boundary, since India rifting away from the rest of Gondwana
    (well before its later collision with Asia) is itself well-documented, just not by a
    precise ridge geometry. Downwelling candidates are spread across the open Panthalassa
    ocean (farthest-point sampled from the sketch's own sea cells) -- subduction has to be
    consuming old ocean crust *somewhere* to balance new crust forming at the rifts, but
    there's no basis for saying precisely where, so this doesn't pretend to more precision
    than that."""
    from . import worldsketch

    plates = pangaea_real_plates()
    by_name = {p.name: p for p in plates}
    americas = np.concatenate([by_name["North America"].boundary_xyz, by_name["South America"].boundary_xyz])
    rest = np.concatenate([by_name[n].boundary_xyz for n in ("Africa", "Eurasia", "Australia")])

    dots = np.clip(americas @ rest.T, -1.0, 1.0)
    angular = np.arccos(dots)
    close_pairs = np.argwhere(angular < _PANGAEA_SUTURE_REACH_RAD)
    midpoints = np.array([geometry.normalize(americas[i] + rest[j]) for i, j in close_pairs])
    suture_candidates = _cluster_positions(midpoints, _PANGAEA_SUTURE_CLUSTER_RAD)
    upwelling = [RidgeTrenchCandidate(position=p, sign=1.0, weight=1.0) for p in suture_candidates]

    india_xyz = by_name["India"].boundary_xyz
    india_idx = rng.choice(len(india_xyz), size=min(_PANGAEA_INDIA_RIFT_CANDIDATES, len(india_xyz)), replace=False)
    upwelling += [RidgeTrenchCandidate(position=p, sign=1.0, weight=1.0) for p in india_xyz[india_idx]]

    cell_xyz = sketch.cell_centers_xyz().reshape(-1, 3)
    ocean_xyz = cell_xyz[~sketch.land.reshape(-1)]
    downwelling_xyz = worldsketch.farthest_point_sites(ocean_xyz, _PANGAEA_DOWNWELLING_CANDIDATES, rng)
    downwelling = [RidgeTrenchCandidate(position=p, sign=-1.0, weight=1.0) for p in downwelling_xyz]

    return upwelling + downwelling


def _target_samples(real_plates: list[RealPlate], rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """`_FIT_SAMPLES_PER_PLATE` (point, target_velocity) pairs per plate: points drawn from a
    Fibonacci-sphere lattice filtered to the plate's own polygon (deterministic, independent
    of `rng`'s draw order elsewhere), velocities from `omega x point`."""
    n = 20000
    i = np.arange(n) + 0.5
    z = 1.0 - 2.0 * i / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, 1.0))
    phi = np.pi * (1.0 + np.sqrt(5.0)) * i
    lattice = np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=-1)

    points_chunks, velocity_chunks = [], []
    for rp in real_plates:
        inside = geometry.points_in_spherical_polygon(lattice, rp.boundary_xyz)
        idx = np.flatnonzero(inside)
        if len(idx) == 0:
            continue
        take = idx if len(idx) <= _FIT_SAMPLES_PER_PLATE else rng.choice(idx, size=_FIT_SAMPLES_PER_PLATE, replace=False)
        pts = lattice[take]
        vel = np.cross(np.broadcast_to(rp.omega, pts.shape), pts)
        points_chunks.append(pts)
        velocity_chunks.append(vel)
    return np.concatenate(points_chunks, axis=0), np.concatenate(velocity_chunks, axis=0)


def fit_mantle_centers(
    real_plates: list[RealPlate],
    candidates: list[RidgeTrenchCandidate],
    rng: np.random.Generator,
) -> list[mantle.ConvectionCenter]:
    """Mantle convection centers "worked backwards" from known (real or, for Pangaea,
    directionally-approximated) plate motion: candidate center *positions* are fixed (real
    ridge/trench locations, or Pangaea's directional stand-ins), and since
    `mantle.flow_at` is linear in each center's `strength` for a fixed position/falloff, the
    strength that best reproduces every plate's own target velocity samples (`omega x point`,
    see `_target_samples`) at every sample point is an ordinary linear least-squares solve --
    `np.linalg.lstsq`'s own minimum-norm behaviour keeps weakly-determined centers (too far
    from any sample point to matter) near zero rather than blowing up."""
    points, targets = _target_samples(real_plates, rng)
    n = len(points)
    m = len(candidates)
    a = np.zeros((n * 3, m))
    for j, c in enumerate(candidates):
        centers = [mantle.ConvectionCenter(position=c.position, strength=1.0, falloff=_CANDIDATE_FALLOFF_RAD)]
        a[:, j] = mantle.flow_at(points, centers).reshape(-1)
    b = targets.reshape(-1)
    strengths, *_ = np.linalg.lstsq(a, b, rcond=None)

    # Sanity clamp: a handful of centers with almost no sample-point leverage can otherwise
    # fit to an implausibly large strength (competing against another to cancel out at the
    # few points that do see both) -- cap at a generous multiple of the reference rate,
    # matching the spread `generate_convection_centers`' own random draw already allows.
    cap = 4.0 * mantle.MANTLE_FLOW_REFERENCE_RATE
    strengths = np.clip(strengths, -cap, cap)

    return [
        mantle.ConvectionCenter(position=c.position, strength=float(s), falloff=_CANDIDATE_FALLOFF_RAD)
        for c, s in zip(candidates, strengths)
    ]
