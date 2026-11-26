"""Thin Python bridge and immutable flat-array snapshot for Rust map rendering (#1039)."""

from __future__ import annotations

import logging
import os
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from scipy.spatial import cKDTree

from . import geometry, plates, projections
from .elevation_lines import line_spacing_rad
from .hydroclimate import hydrology, ice_sheet

if TYPE_CHECKING:
    from .world import World

_logger = logging.getLogger(__name__)
_DISABLE_RUST = os.environ.get("MANTLE_BLOOM_DISABLE_RUST") == "1"

_REQUIRED_RUST_EXPORTS = (
    "render_pool_threads",
    "grid_edt",
    "smooth_biome_field_blend_rs",
    "biome_relative_shade_factor_rs",
    "fill_rects_rs",
    "spill_path_rs",
    "render_biome_or_combined",
    "render_elevation",
    "render_constants_rs",
)


def _resolve_rust_module(mod):
    if mod is None:
        return None
    missing = [name for name in _REQUIRED_RUST_EXPORTS if not hasattr(mod, name)]
    if missing:
        _logger.warning(
            "mantle_bloom_rust is missing render exports %s (stale build?); falling back to Python. "
            "Rebuild with `backend/.venv/bin/python bin/build_rust.py` or `./bin/setup.sh`.",
            missing,
        )
        return None
    return mod


if _DISABLE_RUST:
    _imported_rust = None
else:
    try:
        import mantle_bloom_rust as _imported_rust
    except ImportError:
        _imported_rust = None

_rust = _resolve_rust_module(_imported_rust)
RUST_AVAILABLE: bool = _rust is not None

RUST_RENDER_VIEWS = frozenset({"combined", "biome", "elevation"})
RUST_RENDER_PROJECTIONS = frozenset({"behrmann", "eckert4", "globe"})

# Capture original helper identities at import time so tests that monkeypatch any of
# these functions transparently fall back to the pure-Python path.
_ORIG_SAMPLE_IS_OCEAN = hydrology.sample_is_ocean
_ORIG_SAMPLE_IS_SEA = hydrology.sample_is_sea
_ORIG_SAMPLE_PERENNIAL_LAKE_ICE = hydrology.sample_perennial_lake_ice
_ORIG_SPILL_PASS_LEVEL = hydrology._spill_pass_level
_ORIG_FEWEST_HOPS_PATH = hydrology._fewest_hops_path


def _f64(arr) -> np.ndarray:
    return np.ascontiguousarray(arr, dtype=np.float64)


def _i64(arr) -> np.ndarray:
    return np.ascontiguousarray(arr, dtype=np.int64)


def _u8(arr) -> np.ndarray:
    return np.ascontiguousarray(arr, dtype=np.uint8)


def _bool(arr) -> np.ndarray:
    return np.ascontiguousarray(arr, dtype=np.bool_)


def render_pool_threads() -> int:
    """Number of worker threads in the dedicated Rust rendering pool."""
    if _rust is None:
        return 0
    return int(_rust.render_pool_threads())


def grid_edt(mask: np.ndarray) -> np.ndarray:
    """2D Euclidean distance transform of a boolean grid via Rust."""
    return _rust.grid_edt(_bool(mask))


def smooth_biome_field_blend(
    temperature_c: np.ndarray,
    precipitation_mm: np.ndarray,
    is_ocean: np.ndarray,
    lat_deg: np.ndarray,
    axial_tilt_deg: float,
    has_sea_ice: np.ndarray | None = None,
    sea_ice_weight: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run smoothed biome classification + inline top-2 soft blend in Rust."""
    lat = np.asarray(lat_deg, dtype=np.float64)
    if lat.ndim == 2:
        lat = lat[:, 0]
    return _rust.smooth_biome_field_blend_rs(
        _f64(temperature_c),
        _f64(precipitation_mm),
        _bool(is_ocean),
        _f64(lat),
        float(axial_tilt_deg),
        None if has_sea_ice is None else _bool(has_sea_ice),
        None if sea_ice_weight is None else _f64(sea_ice_weight),
    )


def biome_relative_shade_factor(biome_ids: np.ndarray, elevation_m: np.ndarray) -> np.ndarray:
    """Compute per-class rank-based relief shade factor in Rust."""
    ids = np.asarray(biome_ids)
    elev = np.asarray(elevation_m)
    out = _rust.biome_relative_shade_factor_rs(
        _i64(ids.reshape(-1)),
        _f64(elev.reshape(-1)),
    )
    return out.reshape(ids.shape)


def fill_rects(
    pixels: np.ndarray,
    x0: np.ndarray,
    x1: np.ndarray,
    y0: np.ndarray,
    y1: np.ndarray,
    colors: np.ndarray,
) -> None:
    """Fill axis-aligned rectangles into an existing contiguous uint8 (H, W, C) array."""
    _rust.fill_rects_rs(
        pixels,
        _i64(x0),
        _i64(x1),
        _i64(y0),
        _i64(y1),
        _u8(colors),
    )


def spill_path_can_use_rust() -> bool:
    """True when Rust is enabled and spill-path Numba helpers are not monkeypatched."""
    return (
        RUST_AVAILABLE
        and hydrology._spill_pass_level is _ORIG_SPILL_PASS_LEVEL
        and hydrology._fewest_hops_path is _ORIG_FEWEST_HOPS_PATH
    )


def spill_path(
    indptr: np.ndarray,
    indices: np.ndarray,
    elevation: np.ndarray,
    is_ocean: np.ndarray,
    source: int,
    target: int,
    ceiling: float,
    max_nodes: int,
) -> list[int]:
    """Compute minimax + fewest-hop spill path in Rust without entering Numba."""
    return list(
        _rust.spill_path_rs(
            _i64(indptr),
            _i64(indices),
            _f64(elevation),
            _bool(is_ocean),
            int(source),
            int(target),
            float(ceiling),
            int(max_nodes),
        )
    )


def can_render_rust(
    world: "World",
    projection: str,
    view: str,
    *,
    orig_render_map_png=None,
    current_render_map_png=None,
) -> bool:
    """Return True when the Rust render pipeline can serve this render safely."""
    if not RUST_AVAILABLE or _rust is None:
        return False
    if view not in RUST_RENDER_VIEWS or projection not in RUST_RENDER_PROJECTIONS:
        return False
    if not getattr(world, "plates", None):
        return False
    if getattr(world, "node_cloud_resample_mode", "kdtree") != "kdtree":
        return False
    if (
        hydrology.sample_is_ocean is not _ORIG_SAMPLE_IS_OCEAN
        or hydrology.sample_is_sea is not _ORIG_SAMPLE_IS_SEA
        or hydrology.sample_perennial_lake_ice is not _ORIG_SAMPLE_PERENNIAL_LAKE_ICE
    ):
        return False
    if (
        orig_render_map_png is not None
        and current_render_map_png is not None
        and current_render_map_png is not orig_render_map_png
    ):
        return False
    return True


@dataclass(frozen=True)
class RenderSnapshot:
    """Immutable flat-array snapshot owned by a background Rust render worker."""

    view: str
    projection: str
    width: int
    height: int
    view_rotation: np.ndarray
    stroke_scale_override: float
    grid_h: int
    grid_w: int
    grid_spacing_rad: float
    sea_level_m: float
    axial_tilt_deg: float
    solar_multiplier: float
    smooth_sea_ice: bool
    soft_sea_ice_edge: bool
    show_mountains: bool
    show_plains_plateaus: bool
    tectonics_bbox: tuple[float, float, float, float] | None
    points: np.ndarray
    elevation: np.ndarray
    lake_depth: np.ndarray
    lake_ice: np.ndarray
    glacier_depth: np.ndarray
    hillshade: np.ndarray | None
    channel_depth: np.ndarray | None
    channel_width: np.ndarray | None
    terrain_relief: np.ndarray | None
    hydro_points_same: bool
    hydro_points: np.ndarray | None
    hydro_is_ocean: np.ndarray | None
    hydro_is_sea: np.ndarray | None
    hydro_perennial_lake_ice: np.ndarray | None
    clim_air_temp: np.ndarray
    clim_ocean_temp: np.ndarray
    clim_precip: np.ndarray
    clim_lat_deg: np.ndarray
    clim_elev_m: np.ndarray | None
    clim_is_ocean: np.ndarray | None
    clim_h: int
    clim_w: int
    river_from_xyz: np.ndarray
    river_to_xyz: np.ndarray
    river_width_frac: np.ndarray
    river_color_alpha: np.ndarray
    river_canyon_fill: np.ndarray


def _collect_tectonics_bbox(
    world: "World",
    projection: str,
    view_rotation: np.ndarray,
) -> tuple[float, float, float, float] | None:
    """Compute the projected (min_x, min_y, max_x, max_y) of plate outlines and poles."""
    target_proj = "behrmann" if projection == projections.GLOBE else projection
    chunks: list[np.ndarray] = []
    for plate in world.plates:
        outline = plate.get_bounding_polygon()
        if len(outline) > 0:
            rot_pts = outline @ view_rotation.T
            lat, lon = geometry.xyz_to_latlon(rot_pts)
            x, y = projections.project(target_proj, lat, lon)
            chunks.append(np.stack([x, y], axis=-1))
        speed = float(np.linalg.norm(plate.omega))
        if speed > 1e-15:
            pole_xyz = (plate.omega / speed)[None, :] @ view_rotation.T
            lat, lon = geometry.xyz_to_latlon(pole_xyz)
            x, y = projections.project(target_proj, lat, lon)
            chunks.append(np.stack([x, y], axis=-1))
    if not chunks:
        return None
    all_xy = np.concatenate(chunks, axis=0)
    return (
        float(all_xy[:, 0].min()),
        float(all_xy[:, 1].min()),
        float(all_xy[:, 0].max()),
        float(all_xy[:, 1].max()),
    )


def _extract_river_arrays(world: "World") -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Extract compact flat river segment arrays without running Pillow or Numba."""
    from . import render_image

    empty_pts = np.empty((0, 3), dtype=np.float64)
    empty_f = np.empty(0, dtype=np.float64)
    empty_i = np.empty(0, dtype=np.int32)
    hydro = world.hydrology_cache
    if hydro is None:
        return empty_pts, empty_pts, empty_f, empty_f, empty_i

    style = render_image.RIVER_STYLE
    river_idx, color_alpha, width_frac = render_image._rivers_to_draw(
        world, hide_short_networks=True, style=style
    )
    if len(river_idx) == 0:
        return empty_pts, empty_pts, empty_f, empty_f, empty_i

    water_frac = render_image._river_water_fraction(hydro, river_idx, style)
    from_idx, target_idx, owner, _ = hydrology.flow_segments(hydro, river_idx)
    if len(from_idx) == 0:
        return empty_pts, empty_pts, empty_f, empty_f, empty_i

    color_alpha_seg = np.ascontiguousarray(color_alpha[owner] * water_frac[owner], dtype=np.float64)
    width_frac_seg = np.ascontiguousarray(width_frac[owner], dtype=np.float64)
    canyon_fill_seg = np.ascontiguousarray(np.round(255.0 * (1.0 - water_frac[owner])), dtype=np.int32)
    course = render_image._river_course_points(hydro, river_idx, style, line_spacing_rad(world.node_density))
    from_xyz = np.array(course[from_idx], dtype=np.float64, copy=True)
    to_xyz = np.array(course[target_idx], dtype=np.float64, copy=True)
    return from_xyz, to_xyz, width_frac_seg, color_alpha_seg, canyon_fill_seg


def build_render_snapshot(
    world: "World",
    projection: str,
    view: str,
    width: int,
    height: int,
    view_rotation: np.ndarray | None,
    *,
    show_mountains: bool = False,
    show_plains_plateaus: bool = False,
    smooth_sea_ice: bool = True,
    soft_sea_ice_edge: bool = True,
    stroke_scale_override: float = 1.0,
) -> RenderSnapshot | None:
    """Copy only the flat arrays needed by `view` into an immutable `RenderSnapshot`."""
    from . import render_image
    from .hydroclimate import climate

    rot = np.eye(3, dtype=np.float64) if view_rotation is None else np.array(view_rotation, dtype=np.float64, copy=True)
    collected = plates.collect_all_points(world.plates)
    if collected is None:
        return None
    all_points, all_elevation, _owner = collected
    all_points = np.ascontiguousarray(all_points, dtype=np.float64)
    all_elevation = np.ascontiguousarray(all_elevation, dtype=np.float64)

    lake_depth = np.ascontiguousarray(plates.collect_all_lake_depth(world.plates), dtype=np.float64)
    lake_ice = np.ascontiguousarray(plates.collect_all_lake_ice(world.plates), dtype=np.float64)

    raw_glacier = plates.collect_all_glacier_depth(world.plates)
    if np.any(raw_glacier > 0.0):
        ocean = hydrology.sample_is_ocean(world, all_points, all_elevation <= world.sea_level_m)
        unloaded_bed = all_elevation - plates.collect_all_ice_load_deflection(world.plates)
        grounded = ice_sheet.is_grounded(raw_glacier, unloaded_bed, world.sea_level_m, ocean)
        glacier_depth = np.ascontiguousarray(np.where(grounded, raw_glacier, 0.0), dtype=np.float64)
    else:
        glacier_depth = np.ascontiguousarray(raw_glacier, dtype=np.float64)

    hillshade = None
    channel_depth = None
    channel_width = None
    if view in ("combined", "elevation"):
        cached_hs = getattr(world, "node_hillshade_cache", None)
        if cached_hs is not None and len(cached_hs) == len(all_points):
            hillshade = np.array(cached_hs, dtype=np.float64, copy=True)
        else:
            channel_depth = np.ascontiguousarray(plates.collect_all_channel_depth(world.plates), dtype=np.float64)
            channel_width = np.ascontiguousarray(plates.collect_all_channel_width(world.plates), dtype=np.float64)

    terrain_relief = None
    if view == "elevation" and (show_mountains or show_plains_plateaus):
        cached_kdtree = getattr(world, "node_kdtree_cache", None)
        if cached_kdtree is not None and len(cached_kdtree[0]) == len(all_points):
            tree = cached_kdtree[3]
        else:
            shared = getattr(world, "node_position_tree_cache", None)
            if shared is not None and shared[0].shape[0] == all_points.shape[0]:
                tree = shared[1]
            else:
                tree = cKDTree(all_points)
        terrain_relief = np.ascontiguousarray(
            render_image._classify_terrain_relief(all_points, all_elevation, tree, world.sea_level_m),
            dtype=np.uint8,
        )

    tectonics_bbox = _collect_tectonics_bbox(world, projection, rot) if view == "elevation" else None

    hydro = getattr(world, "hydrology_cache", None)
    hydro_points_same = False
    hydro_points = None
    hydro_is_ocean = None
    hydro_is_sea = None
    hydro_perennial_lake_ice = None
    if hydro is not None and len(hydro.points) > 0:
        if len(hydro.points) == len(all_points) and np.array_equal(hydro.points, all_points):
            hydro_points_same = True
        else:
            hydro_points = np.array(hydro.points, dtype=np.float64, copy=True)
        if len(hydro.is_ocean) == len(hydro.points):
            hydro_is_ocean = np.array(hydro.is_ocean, dtype=np.bool_, copy=True)
        h_sea = getattr(hydro, "is_sea", None)
        if h_sea is not None and len(h_sea) == len(hydro.points):
            hydro_is_sea = np.array(h_sea, dtype=np.bool_, copy=True)
        h_lake_ice = getattr(hydro, "perennial_lake_ice", None)
        if h_lake_ice is not None and h_lake_ice.shape == (len(hydro.points),):
            hydro_perennial_lake_ice = np.array(h_lake_ice, dtype=np.bool_, copy=True)

    fields = climate.compute_climate_cached(world)
    clim_h, clim_w = fields.ocean_temperature_c.shape
    clim_air_temp = np.ascontiguousarray(fields.air_temperature_c, dtype=np.float64).ravel().copy()
    clim_ocean_temp = np.ascontiguousarray(fields.ocean_temperature_c, dtype=np.float64).ravel().copy()
    clim_precip = np.ascontiguousarray(fields.precipitation_mm, dtype=np.float64).ravel().copy()
    clim_lat_deg = np.array(fields.lat_deg, dtype=np.float64, copy=True)
    need_lake_fallback = hydro_perennial_lake_ice is None and bool(np.any(lake_ice > 0.0))
    clim_elev_m = (
        np.ascontiguousarray(fields.elevation_m, dtype=np.float64).ravel().copy()
        if need_lake_fallback
        else None
    )
    clim_is_ocean = (
        np.ascontiguousarray(fields.is_ocean, dtype=np.bool_).ravel().copy()
        if need_lake_fallback
        else None
    )

    river_from, river_to, river_wf, river_ca, river_cf = _extract_river_arrays(world)
    grid_h, grid_w = render_image.biome_grid_dimensions(world.climate_density)
    grid_spacing = min(render_image.GRID_SPACING_RAD, line_spacing_rad(world.node_density))

    return RenderSnapshot(
        view=view,
        projection=projection,
        width=int(width),
        height=int(height),
        view_rotation=rot,
        stroke_scale_override=float(stroke_scale_override),
        grid_h=int(grid_h),
        grid_w=int(grid_w),
        grid_spacing_rad=float(grid_spacing),
        sea_level_m=float(world.sea_level_m),
        axial_tilt_deg=float(world.axial_tilt_deg),
        solar_multiplier=float(world.solar_multiplier),
        smooth_sea_ice=bool(smooth_sea_ice),
        soft_sea_ice_edge=bool(soft_sea_ice_edge),
        show_mountains=bool(show_mountains),
        show_plains_plateaus=bool(show_plains_plateaus),
        tectonics_bbox=tectonics_bbox,
        points=all_points,
        elevation=all_elevation,
        lake_depth=lake_depth,
        lake_ice=lake_ice,
        glacier_depth=glacier_depth,
        hillshade=hillshade,
        channel_depth=channel_depth,
        channel_width=channel_width,
        terrain_relief=terrain_relief,
        hydro_points_same=hydro_points_same,
        hydro_points=hydro_points,
        hydro_is_ocean=hydro_is_ocean,
        hydro_is_sea=hydro_is_sea,
        hydro_perennial_lake_ice=hydro_perennial_lake_ice,
        clim_air_temp=clim_air_temp,
        clim_ocean_temp=clim_ocean_temp,
        clim_precip=clim_precip,
        clim_lat_deg=clim_lat_deg,
        clim_elev_m=clim_elev_m,
        clim_is_ocean=clim_is_ocean,
        clim_h=int(clim_h),
        clim_w=int(clim_w),
        river_from_xyz=river_from,
        river_to_xyz=river_to,
        river_width_frac=river_wf,
        river_color_alpha=river_ca,
        river_canyon_fill=river_cf,
    )


def render_snapshot(
    snap: RenderSnapshot,
    *,
    want_png: bool = True,
    want_rgb: bool = False,
) -> tuple[bytes | None, np.ndarray | None]:
    """Execute a `RenderSnapshot` in Rust on the dedicated render pool (GIL released)."""
    if snap.view in ("combined", "biome"):
        return _rust.render_biome_or_combined(
            snap.view,
            snap.projection,
            snap.width,
            snap.height,
            snap.view_rotation,
            snap.grid_h,
            snap.grid_w,
            snap.sea_level_m,
            snap.axial_tilt_deg,
            snap.solar_multiplier,
            snap.smooth_sea_ice,
            snap.soft_sea_ice_edge,
            snap.points,
            snap.elevation,
            snap.lake_depth,
            snap.lake_ice,
            snap.glacier_depth,
            snap.hillshade,
            snap.channel_depth,
            snap.channel_width,
            snap.hydro_points_same,
            snap.hydro_points,
            snap.hydro_is_ocean,
            snap.hydro_is_sea,
            snap.hydro_perennial_lake_ice,
            snap.clim_air_temp,
            snap.clim_ocean_temp,
            snap.clim_precip,
            snap.clim_lat_deg,
            snap.clim_elev_m,
            snap.clim_is_ocean,
            snap.clim_h,
            snap.clim_w,
            snap.river_from_xyz,
            snap.river_to_xyz,
            snap.river_width_frac,
            snap.river_color_alpha,
            snap.river_canyon_fill,
            want_png,
            want_rgb,
            snap.stroke_scale_override,
        )
    if snap.view == "elevation":
        return _rust.render_elevation(
            snap.projection,
            snap.width,
            snap.height,
            snap.view_rotation,
            snap.grid_spacing_rad,
            snap.sea_level_m,
            snap.axial_tilt_deg,
            snap.solar_multiplier,
            snap.smooth_sea_ice,
            snap.soft_sea_ice_edge,
            snap.show_mountains,
            snap.show_plains_plateaus,
            snap.tectonics_bbox,
            snap.points,
            snap.elevation,
            snap.lake_depth,
            snap.lake_ice,
            snap.glacier_depth,
            snap.hillshade,
            snap.channel_depth,
            snap.channel_width,
            snap.terrain_relief,
            snap.hydro_points_same,
            snap.hydro_points,
            snap.hydro_is_ocean,
            snap.hydro_is_sea,
            snap.hydro_perennial_lake_ice,
            snap.clim_air_temp,
            snap.clim_ocean_temp,
            snap.clim_lat_deg,
            snap.clim_elev_m,
            snap.clim_is_ocean,
            snap.clim_h,
            snap.clim_w,
            snap.river_from_xyz,
            snap.river_to_xyz,
            snap.river_width_frac,
            snap.river_color_alpha,
            snap.river_canyon_fill,
            want_png,
            want_rgb,
            snap.stroke_scale_override,
        )
    raise ValueError(f"Unsupported Rust render view {snap.view!r}")
