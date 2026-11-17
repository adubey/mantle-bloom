#!/usr/bin/env python3
"""RFC 0019 G3 climate-only comparison for registered tuning saves.

Usage:
  backend/.venv/bin/python bin/debug/rfc0019_g3_compare.py <save>... --out <private.json>

Runs the current climate and an offline conservative moisture solve on each requested save.
Use only the registered tuning saves until parameters are frozen; never pass held-out saves
during calibration. The source conversion and acceptance limits are in
analysis/rfc0017/g3-preregistration.md. This script does not step or mutate production worlds.
"""
from __future__ import annotations

import argparse
from copy import copy
import hashlib
from itertools import product
import json
import resource
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from PIL import Image, ImageDraw
from scipy.spatial import cKDTree

from rfc0019_transport_prototypes import transition_rates, solve_stationary

EARTH_RADIUS_M = 6_371_000.0
AIR_DENSITY_KG_M3 = 1.2
SECONDS_PER_YEAR = 365.25 * 86_400.0
RAIN_DAYS = 9.1
BAND_EDGES_DEG = np.arange(0.0, 90.1, 15.0)
ERA5_CLASS_RAINOUT_DAYS = {
    "land": np.array([8.12, 11.56, 8.79, 5.61, 4.84, 3.79]),
    "ocean": np.array([9.08, 12.65, 6.22, 3.84, 3.46, 5.88]),
}
WARMUP_COUNT = 1
TIMING_REPEATS = 3


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("saves", nargs="+", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument(
        "--map-dir", type=Path,
        help="optional private directory for common-scale rain-map quicklooks",
    )
    parser.add_argument("--density", choices=("1", "4", "both"), default="both")
    parser.add_argument("--wind", choices=("diagnostic", "cfd", "both"), default="both")
    parser.add_argument("--solver", choices=("finite_volume", "wind_arrows", "both"), default="both")
    parser.add_argument("--rain-days", type=float, default=RAIN_DAYS)
    parser.add_argument("--effective-depth-km", type=float, default=2.0)
    parser.add_argument("--wind-multiplier", type=float, default=1.0)
    parser.add_argument(
        "--land-proxy-multiplier", type=float, default=1.0,
        help="scale lake/river/vegetation source proxies together; ocean target stays unchanged",
    )
    parser.add_argument(
        "--ocean-source-mapping", choices=("humidity_proxy", "bulk_flux"),
        default="humidity_proxy",
        help="use the calibrated humidity proxy or a bulk sea-air humidity-deficit flux",
    )
    parser.add_argument("--ocean-transfer-coefficient", type=float, default=0.003)
    parser.add_argument(
        "--spatial-era5-rainout", action="store_true",
        help="use G1 class-mean column-water/precip ratios by 15-degree band and land/ocean",
    )
    parser.add_argument(
        "--saturation-convection", action="store_true",
        help="add a conservative rain sink for column water above a saturation threshold",
    )
    parser.add_argument(
        "--orographic-lift", action="store_true",
        help="add the preregistered bounded windward-lift rain sink for issue #913",
    )
    parser.add_argument(
        "--lift-comparison", action="store_true",
        help="run uniform-sink control and lift-aware treatment on identical inputs",
    )
    parser.add_argument("--relative-humidity-threshold", type=float, default=0.8)
    parser.add_argument("--convection-days", type=float, default=1.0,
                        help="e-folding time for rainout of water above the saturation threshold")
    parser.add_argument(
        "--sweep-ocean-coefficient", action="store_true",
        help="test the 0.001/0.003/0.005 bulk-transfer bracket on tuning saves",
    )
    parser.add_argument(
        "--sweep-land-proxy", action="store_true",
        help="test land-proxy multipliers 1, 5, and 10 on tuning saves",
    )
    parser.add_argument(
        "--sweep-tuning", action="store_true",
        help="run preregistered 1.5/2.25/3 km, 5/9.1/18 day, and 0.5/1/2 wind sensitivities",
    )
    return parser.parse_args()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _grid_area(height: int, width: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return climate-grid latitudes/edges and exact spherical cell areas."""
    lat_edges = np.linspace(np.pi / 2, -np.pi / 2, height + 1)
    lat = 0.5 * (lat_edges[:-1] + lat_edges[1:])
    dlon = 2.0 * np.pi / width
    row_area = EARTH_RADIUS_M**2 * dlon * (
        np.sin(lat_edges[:-1]) - np.sin(lat_edges[1:])
    )
    area = np.broadcast_to(row_area[:, None], (height, width)).copy()
    return lat, lat_edges, area


def _solver_order_rainout_sink(rainout_days: np.ndarray) -> np.ndarray:
    """Return a north-to-south map's sink in the solver's south-to-north order. #913"""
    return 1.0 / (np.flip(rainout_days, axis=0).ravel() * 86_400.0)


def _check_asymmetric_rainout_fixture() -> None:
    """Keep land/ocean rainout classes aligned across the latitude flip. #913"""
    north_south_ocean = np.array([[False, True, False], [True, False, True]])
    north_south_days = np.where(north_south_ocean, 9.0, 2.0)
    actual = _solver_order_rainout_sink(north_south_days)
    expected = 1.0 / (
        np.where(np.flip(north_south_ocean, axis=0), 9.0, 2.0).ravel() * 86_400.0
    )
    if not np.array_equal(actual, expected):
        raise AssertionError("ERA5 land/ocean sinks do not follow the solver latitude order")


def _orographic_lift_rate(fields, area: np.ndarray,
                          wind_u: np.ndarray | None = None,
                          wind_v: np.ndarray | None = None,
                          smoothing_sigma_km: float = 250.0,
                          moisture_scale_m: float = 2250.0,
                          cap_days: float = 5.0) -> tuple[np.ndarray, dict]:
    """Return bounded sink from positive ascent over physically smoothed terrain. #913"""
    xyz = fields.world_xyz.reshape(-1, 3)
    points = xyz * EARTH_RADIUS_M
    sigma_m = smoothing_sigma_km * 1000.0
    neighbors = cKDTree(points).query_ball_point(points, r=3.0 * sigma_m)
    elevation = fields.elevation_m.ravel()
    smooth = np.empty_like(elevation)
    for index, member_ids in enumerate(neighbors):
        members = np.asarray(member_ids, dtype=int)
        cosine = np.clip(xyz[members] @ xyz[index], -1.0, 1.0)
        distance_m = EARTH_RADIUS_M * np.arccos(cosine)
        weights = np.exp(-0.5 * (distance_m / sigma_m) ** 2)
        smooth[index] = np.dot(weights, elevation[members]) / weights.sum()
    smooth = smooth.reshape(fields.elevation_m.shape)

    latitude = np.deg2rad(fields.lat_deg)
    dz_dnorth = np.gradient(smooth, latitude, axis=0) / EARTH_RADIUS_M
    longitude_width = 2.0 * np.pi / smooth.shape[1]
    east_spacing = (
        EARTH_RADIUS_M * np.maximum(np.cos(latitude), 1e-6)[:, None] * longitude_width
    )
    dz_deast = (np.roll(smooth, -1, axis=1) - np.roll(smooth, 1, axis=1)) / (
        2.0 * east_spacing
    )
    wind_u = fields.wind_u if wind_u is None else wind_u
    wind_v = fields.wind_v if wind_v is None else wind_v
    ascent = np.maximum(wind_u * dz_deast + wind_v * dz_dnorth, 0.0)
    cap_rate = 1.0 / (cap_days * 86_400.0)
    lift_rate = np.minimum(ascent / moisture_scale_m, cap_rate)
    total_area = max(float(area.sum()), 1e-300)
    return lift_rate, {
        "smoothing_sigma_km": smoothing_sigma_km,
        "moisture_scale_m": moisture_scale_m,
        "cap_days": cap_days,
        "area_fraction_with_lift": float(area[lift_rate > 0.0].sum() / total_area),
        "area_fraction_at_cap": float(area[lift_rate >= cap_rate].sum() / total_area),
        "area_weighted_mean_rate_s-1": _weighted_mean(lift_rate, area),
        "area_weighted_p90_rate_s-1": _weighted_quantiles(lift_rate, area)["p90"],
    }


def _cache_on_grid(cache, shape: tuple[int, int]):
    """Area-average the lagged climate inputs when diagnosing a different grid density."""
    if cache is None:
        return None, False
    source_h, source_w = cache.precipitation_mm.shape
    target_h, target_w = shape
    if (source_h, source_w) == shape:
        return cache, False
    if source_h % target_h or source_w % target_w:
        raise ValueError(f"cannot area-average cache {(source_h, source_w)} onto {shape}")
    factor_h, factor_w = source_h // target_h, source_w // target_w
    edges = np.linspace(np.pi / 2, -np.pi / 2, source_h + 1)
    row_area = np.sin(edges[:-1]) - np.sin(edges[1:])

    def average(values: np.ndarray) -> np.ndarray:
        blocks = np.asarray(values).reshape(target_h, factor_h, target_w, factor_w)
        weights = row_area.reshape(target_h, factor_h, 1, 1)
        return np.sum(blocks * weights, axis=(1, 3)) / (
            np.sum(row_area.reshape(target_h, factor_h), axis=1)[:, None] * factor_w
        )

    ocean_fraction = average(cache.is_ocean.astype(float))
    lagged = SimpleNamespace(
        is_ocean=ocean_fraction >= 0.5,
        ocean_temperature_c=average(cache.ocean_temperature_c),
        air_temperature_c=average(cache.air_temperature_c),
        precipitation_mm=average(cache.precipitation_mm),
        lat_deg=90.0 - (np.arange(target_h) + 0.5) * (180.0 / target_h),
    )
    return lagged, True


def _surface_source_proxies(fields, world, climate, hydrology,
                            sampling_times: list[float] | None = None) -> dict[str, np.ndarray]:
    sample_start = time.perf_counter()
    elevation, is_ocean, lake_depth, channel_depth, is_sea = climate._sample_elevation_and_crust(
        world, fields.world_xyz
    )
    if sampling_times is not None:
        sampling_times.append(time.perf_counter() - sample_start)
    vegetation = climate._vegetation_transpiration_source(world, elevation, is_ocean)
    temperature = fields.air_temperature_c
    water_humidity = climate.open_water_humidity_g_kg(temperature)
    sea_fraction = np.where(
        is_sea, climate.SEA_EVAPORATION_FRACTION, climate.LAKE_EVAPORATION_FRACTION
    )
    lake_fraction = np.clip(lake_depth / climate.LAKE_EVAPORATION_REFERENCE_DEPTH_M, 0.0, 1.0)
    river_fraction = np.clip(channel_depth / climate.RIVER_EVAPORATION_REFERENCE_DEPTH_M, 0.0, 1.0)
    lake_sea = lake_fraction * sea_fraction * water_humidity
    river = river_fraction * climate.RIVER_EVAPORATION_HUMIDITY_G_KG
    unfrozen = temperature >= hydrology.FREEZE_POINT_C
    lake_sea = np.where(unfrozen & ~is_ocean, lake_sea, 0.0)
    river = np.where(unfrozen & ~is_ocean, river, 0.0)
    vegetation = np.where(~is_ocean, vegetation, 0.0)
    ocean_target = climate.open_water_humidity_g_kg(fields.ocean_temperature_c)
    ocean_target = np.where(is_ocean, ocean_target, 0.0)
    return {
        "ocean_target": ocean_target,
        "lake_sea_proxy": lake_sea,
        "river_proxy": river,
        "vegetation_proxy": vegetation,
        "is_ocean": is_ocean,
    }


def _compute_candidate_fields(world, height, width, climate, biomes):
    """Time the zero-rain biome pass that production would replace with candidate rain."""
    classification_times = []
    original = biomes.smooth_biome_field

    def timed_classification(*args, **kwargs):
        start = time.perf_counter()
        try:
            return original(*args, **kwargs)
        finally:
            classification_times.append(time.perf_counter() - start)

    biomes.smooth_biome_field = timed_classification
    try:
        fields = climate.compute_climate(world, height, width, skip_moisture=True)
    finally:
        biomes.smooth_biome_field = original
    if len(classification_times) != 1:
        raise AssertionError("expected one zero-rain biome classification")
    return fields, classification_times[0]


def _flux_sources(
    proxies: dict[str, np.ndarray], area: np.ndarray, depth_m: float, tau_s: float,
    land_proxy_multiplier: float = 1.0, ocean_bulk_flux: np.ndarray | None = None,
):
    scale = AIR_DENSITY_KG_M3 * depth_m / (1000.0 * tau_s)
    tagged = {
        name: values * scale * area * (1.0 if name == "ocean_target" else land_proxy_multiplier)
        for name, values in proxies.items()
        if name.endswith("_proxy") or name == "ocean_target"
        if name != "ocean_target" or ocean_bulk_flux is None
    }
    if ocean_bulk_flux is not None:
        tagged["ocean_bulk_flux"] = ocean_bulk_flux * area
    return tagged, sum(tagged.values())


def _solve_saturation_convection(
    source: np.ndarray,
    area: np.ndarray,
    saturation_water: np.ndarray,
    background_tau_s: float,
    convection_tau_s: float,
    rates,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Solve a stationary budget with background rain plus rain from excess column water.

    The sink is `W/tau_background + max(W - Wsat, 0)/tau_convection`. A fixed active set
    makes it a linear solve; the additional `Wsat/tau_convection` term is moved to an
    algebraic effective source, then the physical rain and ledger are reconstructed from
    the solved column water.
    """
    background_rate = 1.0 / background_tau_s
    convective_rate = 1.0 / convection_tau_s
    active = np.ones(source.size, dtype=bool)
    effective_source = source
    mass = None
    solve_metrics = None
    max_active_iterations = 50
    for iteration in range(1, max_active_iterations + 1):
        sink = background_rate + convective_rate * active
        effective_source = source + convective_rate * active * saturation_water * area
        mass, _, solve_metrics = solve_stationary(effective_source, sink, rates)
        column_water = mass / area
        updated = column_water > saturation_water
        if np.array_equal(updated, active):
            break
        active = updated
    else:
        raise ValueError("saturation-triggered rainout active set did not converge")

    assert mass is not None and solve_metrics is not None
    rain = background_rate * mass + convective_rate * active * np.maximum(
        mass - saturation_water * area, 0.0
    )
    inflow = np.zeros(source.size)
    outflow = np.zeros(source.size)
    for origin, destination, rate in rates:
        flux = rate * mass[origin]
        outflow[origin] += flux
        inflow[destination] += flux
    residual = source + inflow - outflow - rain
    total_source = max(float(source.sum()), 1e-300)
    if float(np.max(np.abs(residual))) / total_source > 1e-8:
        raise AssertionError("saturation-triggered physical budget failed cell closure")
    rain += residual
    if np.any(rain < 0.0):
        raise AssertionError("saturation-triggered ledger correction made rain negative")
    closure = abs(float(rain.sum()) - float(source.sum())) / total_source
    if closure > 1e-8:
        raise AssertionError("saturation-triggered physical budget failed global closure")
    metrics = dict(solve_metrics)
    metrics.update({
        "global_closure_relative": closure,
        "max_cell_residual_relative": float(np.max(np.abs(residual))) / total_source,
        "saturation_active_fraction": float(np.mean(active)),
        "saturation_active_set_iterations": iteration,
    })
    return mass, rain, metrics


def _weighted_mean(values: np.ndarray, weights: np.ndarray) -> float | None:
    total = float(weights.sum())
    return None if total <= 0.0 else float(np.sum(values * weights) / total)


def _rain_metrics(fields, rain_mm_yr: np.ndarray, area: np.ndarray, biomes,
                  biome_ids: np.ndarray | None = None) -> dict:
    land = ~fields.is_ocean
    rain_mass = rain_mm_yr * area
    total = float(rain_mass.sum())
    land_total = float(rain_mass[land].sum())
    ocean_total = float(rain_mass[~land].sum())
    band_means = []
    absolute_lat = np.abs(fields.lat_deg)
    for lo, hi in zip(BAND_EDGES_DEG[:-1], BAND_EDGES_DEG[1:]):
        band = (absolute_lat >= lo) & (absolute_lat < hi)
        row_area = area[band]
        band_land = land[band]
        band_ocean = ~band_land
        band_means.append({
            "lat_abs_deg": [float(lo), float(hi)],
            "land_mm_yr": _weighted_mean(rain_mm_yr[band][band_land], row_area[band_land]),
            "ocean_mm_yr": _weighted_mean(rain_mm_yr[band][band_ocean], row_area[band_ocean]),
        })
    codes = biomes.KOPPEN_CODES
    land_area = float(area[land].sum())
    koppen = {}
    for group in "ABCDE":
        group_mask = np.zeros(fields.biome_ids.shape, dtype=bool)
        for biome_id, code in enumerate(codes):
            if code.startswith(group):
                group_mask |= (fields.biome_ids if biome_ids is None else biome_ids) == biome_id
        koppen[group] = float(area[group_mask & land].sum() / max(land_area, 1e-300))
    return {
        "land_rain_mm_yr": _weighted_mean(rain_mm_yr[land], area[land]),
        "ocean_rain_mm_yr": _weighted_mean(rain_mm_yr[~land], area[~land]),
        "land_rain_share": land_total / max(total, 1e-300),
        "ocean_rain_share": ocean_total / max(total, 1e-300),
        "global_mean_rain_mm_yr": total / float(area.sum()),
        "latitude_bands": band_means,
        "koppen_group_shares": koppen,
    }


def _weighted_quantiles(values: np.ndarray, weights: np.ndarray) -> dict[str, float | None]:
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(valid):
        return {"p10": None, "p50": None, "p90": None}
    values = values[valid].ravel()
    weights = weights[valid].ravel()
    order = np.argsort(values)
    values = values[order]
    cumulative = np.cumsum(weights[order])
    total = cumulative[-1]
    return {
        name: float(values[min(np.searchsorted(cumulative, fraction * total), values.size - 1)])
        for name, fraction in (("p10", 0.1), ("p50", 0.5), ("p90", 0.9))
    }


def _area_average(values: np.ndarray, area: np.ndarray,
                  target_shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray]:
    target_h, target_w = target_shape
    height, width = values.shape
    if height % target_h or width % target_w:
        raise ValueError(f"cannot area-average {(height, width)} onto {target_shape}")
    factor_h, factor_w = height // target_h, width // target_w
    blocks = values.reshape(target_h, factor_h, target_w, factor_w)
    area_blocks = area.reshape(target_h, factor_h, target_w, factor_w)
    block_area = area_blocks.sum(axis=(1, 3))
    mean = np.sum(blocks * area_blocks, axis=(1, 3)) / np.maximum(block_area, 1e-300)
    return mean, block_area


def _rms_pair(reference: np.ndarray, candidate: np.ndarray, area: np.ndarray) -> dict[str, float]:
    norm = max(float(np.sum(area * reference**2)), 1e-300)
    raw = float(np.sqrt(np.sum(area * (candidate - reference) ** 2) / norm))
    scale = float(np.sum(area * candidate * reference)
                  / max(float(np.sum(area * candidate**2)), 1e-300))
    fit = float(np.sqrt(np.sum(area * (scale * candidate - reference) ** 2) / norm))
    return {"raw_rel_rms": raw, "bestfit_scale": scale, "bestfit_rel_rms": fit}


def _rain_map_color(values: np.ndarray, maximum: float) -> np.ndarray:
    palette = np.array([
        [68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37],
    ], dtype=float)
    maximum = max(maximum, 1e-12)
    scaled = np.clip(np.log1p(np.maximum(values, 0.0)) / np.log1p(maximum), 0.0, 1.0)
    position = scaled * (len(palette) - 1)
    lower = np.floor(position).astype(int)
    upper = np.minimum(lower + 1, len(palette) - 1)
    weight = (position - lower)[..., None]
    return np.rint(palette[lower] * (1.0 - weight) + palette[upper] * weight).astype(np.uint8)


def _difference_map_color(values: np.ndarray, maximum: float) -> np.ndarray:
    normalized = np.clip(values / max(maximum, 1e-300), -1.0, 1.0)
    neutral = np.full((*values.shape, 3), 245.0)
    negative = np.array([35.0, 85.0, 180.0])
    positive = np.array([190.0, 45.0, 45.0])
    color = neutral.copy()
    below = normalized < 0.0
    above = normalized > 0.0
    color[below] = neutral[below] * (1.0 + normalized[below, None]) - negative * normalized[below, None]
    color[above] = neutral[above] * (1.0 - normalized[above, None]) + positive * normalized[above, None]
    return np.clip(np.rint(color), 0, 255).astype(np.uint8)


def _save_common_scale_rain_map(
    map_dir: Path, save: Path, density: int, wind: str, solver: str, arm: str,
    current_rain: np.ndarray, candidate_rain: np.ndarray, area: np.ndarray,
) -> str:
    target_shape = (current_rain.shape if density == 1 else
                    (current_rain.shape[0] // 4, current_rain.shape[1] // 4))
    current, _ = _area_average(current_rain, area, target_shape)
    candidate, _ = _area_average(candidate_rain, area, target_shape)
    difference = candidate - current
    shared_max = float(np.percentile(np.concatenate((current.ravel(), candidate.ravel())), 99.5))
    difference_max = float(np.percentile(np.abs(difference), 99.5))
    zoom = 4
    panel_width = target_shape[1] * zoom
    panel_height = target_shape[0] * zoom
    gutter, header, footer = 18, 34, 28
    image = Image.new("RGB", (gutter + 3 * (panel_width + gutter), header + panel_height + footer),
                      (248, 248, 248))
    draw = ImageDraw.Draw(image)
    labels = ("Current source", "Conservative candidate", "Candidate − current")
    fields = (current, candidate, difference)
    for index, (label, field) in enumerate(zip(labels, fields)):
        left = gutter + index * (panel_width + gutter)
        draw.text((left, 8), label, fill=(25, 25, 25))
        colors = (_rain_map_color(field, shared_max) if index < 2
                  else _difference_map_color(field, difference_max))
        raster = Image.fromarray(colors, mode="RGB").resize(
            (panel_width, panel_height), resample=Image.Resampling.NEAREST
        )
        image.paste(raster, (left, header))
        draw.rectangle((left, header, left + panel_width - 1, header + panel_height - 1),
                       outline=(70, 70, 70))
        draw.text((left, header + panel_height + 6), "90°N", fill=(25, 25, 25))
        draw.text((left + panel_width // 2 - 12, header + panel_height + 6), "0°", fill=(25, 25, 25))
        draw.text((left + panel_width - 26, header + panel_height + 6), "90°S", fill=(25, 25, 25))
    draw.text((gutter, image.height - 13),
              f"Rain (mm/yr), shared log scale to p99.5={shared_max:.0f}; "
              f"difference scale ±{difference_max:.0f}; longitude 0° → 360°",
              fill=(25, 25, 25))
    map_dir.mkdir(parents=True, exist_ok=True)
    name = f"{save.stem}-d{density}-{wind}-{solver}-{arm}-common-scale.png"
    image.save(map_dir / name, optimize=True)
    return name


def _surface_quantiles(rain: np.ndarray, area: np.ndarray,
                       is_ocean: np.ndarray) -> dict[str, dict[str, float | None]]:
    return {
        "all": _weighted_quantiles(rain, area),
        "land": _weighted_quantiles(rain, area * ~is_ocean),
        "ocean": _weighted_quantiles(rain, area * is_ocean),
    }


def _hemisphere_bands(fields, rain: np.ndarray, area: np.ndarray) -> dict[str, list[dict]]:
    absolute_lat = np.abs(fields.lat_deg)[:, None]
    hemisphere_masks = {
        "north": np.broadcast_to(fields.lat_deg[:, None] >= 0.0, rain.shape),
        "south": np.broadcast_to(fields.lat_deg[:, None] < 0.0, rain.shape),
    }
    output = {}
    for name, hemisphere in hemisphere_masks.items():
        bands = []
        for lower, upper in zip(BAND_EDGES_DEG[:-1], BAND_EDGES_DEG[1:]):
            band = hemisphere & (absolute_lat >= lower) & (absolute_lat < upper)
            land = band & ~fields.is_ocean
            ocean = band & fields.is_ocean
            bands.append({
                "lat_abs_deg": [float(lower), float(upper)],
                "land_mm_yr": _weighted_mean(rain[land], area[land]),
                "ocean_mm_yr": _weighted_mean(rain[ocean], area[ocean]),
            })
        output[name] = bands
    return output


def _windward_leeward_masks(fields, area: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    """Select extreme 10% terrain-ascent/descent land masks from the same wind and terrain."""
    lat = np.deg2rad(fields.lat_deg)
    dz_dnorth = np.gradient(fields.elevation_m, lat, axis=0) / EARTH_RADIUS_M
    dlon = 2.0 * np.pi / fields.elevation_m.shape[1]
    east_spacing = EARTH_RADIUS_M * np.maximum(np.cos(lat), 1e-6)[:, None] * dlon
    dz_deast = (np.roll(fields.elevation_m, -1, axis=1)
                - np.roll(fields.elevation_m, 1, axis=1)) / (2.0 * east_spacing)
    ascent = fields.wind_u * dz_deast + fields.wind_v * dz_dnorth
    land = ~fields.is_ocean
    positive = land & (ascent > 0.0)
    negative = land & (ascent < 0.0)
    up_cut = _weighted_quantiles(ascent[positive], area[positive])["p90"]
    down_cut = _weighted_quantiles(ascent[negative], area[negative])["p10"]
    windward = positive & (ascent >= up_cut) if up_cut is not None else positive
    leeward = negative & (ascent <= down_cut) if down_cut is not None else negative
    return windward, leeward, {
        "windward_ascent_cut_m_s": up_cut,
        "leeward_descent_cut_m_s": down_cut,
        "windward_area_fraction_of_land": float(area[windward].sum()
                                                  / max(area[land].sum(), 1e-300)),
        "leeward_area_fraction_of_land": float(area[leeward].sum()
                                                 / max(area[land].sum(), 1e-300)),
    }


def _extended_rain_diagnostics(fields, rain: np.ndarray, area: np.ndarray,
                               windward: np.ndarray, leeward: np.ndarray,
                               koppen_dry_fraction: float) -> dict:
    high_lat_land = ~fields.is_ocean & (np.abs(fields.lat_deg)[:, None] >= 60.0)
    windward_mean = _weighted_mean(rain[windward], area[windward])
    leeward_mean = _weighted_mean(rain[leeward], area[leeward])
    return {
        "intensity_quantiles_mm_yr": _surface_quantiles(rain, area, fields.is_ocean),
        "hemisphere_bands": _hemisphere_bands(fields, rain, area),
        "high_latitude_land_mean_mm_yr": _weighted_mean(rain[high_lat_land], area[high_lat_land]),
        "koppen_dry_fraction": koppen_dry_fraction,
        "windward_mean_mm_yr": windward_mean,
        "leeward_mean_mm_yr": leeward_mean,
        "windward_leeward_area_weighted_pair": (
            None if windward_mean is None or leeward_mean is None
            else windward_mean / max(leeward_mean, 1e-300)
        ),
    }


def _compare_one(world, save: Path, density: int, wind_mode: str, solver_names: tuple[str, ...],
                 args, climate, hydrology, biomes, original_cache) -> dict:
    height, width = climate.grid_dimensions(float(density))
    world.climate_density = float(density)
    world.climate_cache, cache_resampled = _cache_on_grid(
        original_cache, (height, width)
    )
    requested_wind = wind_mode
    world.wind_model = wind_mode
    wind_available = wind_mode != "cfd" or world.atmosphere_cfd_state is not None
    if not wind_available:
        raise ValueError(f"{save.name}: requested CFD wind but save has no atmosphere_cfd_state")

    current_times = []
    current = None
    for repeat in range(WARMUP_COUNT + TIMING_REPEATS):
        start = time.perf_counter()
        current = climate.compute_climate(world, height, width)
        elapsed = time.perf_counter() - start
        if repeat >= WARMUP_COUNT:
            current_times.append(elapsed)
    assert current is not None
    _, _, area = _grid_area(height, width)
    current_metrics = _rain_metrics(current, current.precipitation_mm, area, biomes)
    windward, leeward, exposure = _windward_leeward_masks(current, area)
    current_extended = _extended_rain_diagnostics(
        current, current.precipitation_mm, area, windward, leeward,
        current_metrics["koppen_group_shares"]["B"],
    )
    scheme_results = {}
    for solver_name in solver_names:
        finite_volume = solver_name == "finite_volume"
        solver_times = []
        candidate_update_times = []
        component_samples = []
        saturation_metrics = None
        lift_metrics = None
        lift_sink_solver = None
        mass = rain = closure = None
        candidate_fields = tagged_sources = source = u = v = area = None
        for repeat in range(WARMUP_COUNT + TIMING_REPEATS):
            update_start = time.perf_counter()
            climate_start = time.perf_counter()
            candidate_fields, zero_rain_biome_elapsed = _compute_candidate_fields(
                world, height, width, climate, biomes
            )
            climate_elapsed = time.perf_counter() - climate_start
            lat, lat_edges, area = _grid_area(height, width)
            proxy_start = time.perf_counter()
            duplicate_surface_samples = []
            proxies = _surface_source_proxies(
                candidate_fields, world, climate, hydrology, duplicate_surface_samples
            )
            ocean_bulk_flux = None
            if args.ocean_source_mapping == "bulk_flux":
                saturation = climate.saturation_specific_humidity_g_kg(
                    candidate_fields.ocean_temperature_c
                ) / 1000.0
                ambient = current.humidity_g_kg / 1000.0
                humidity_deficit = np.maximum(saturation - ambient, 0.0)
                wind_speed = np.hypot(candidate_fields.wind_u, candidate_fields.wind_v)
                wind_speed *= args.wind_multiplier
                ocean_bulk_flux = (
                    AIR_DENSITY_KG_M3 * args.ocean_transfer_coefficient
                    * wind_speed * humidity_deficit
                )
                ocean_bulk_flux = np.where(candidate_fields.is_ocean, ocean_bulk_flux, 0.0)
            tagged_sources, source = _flux_sources(
                proxies, area, args.effective_depth_km * 1000.0, 9.1 * 86_400.0,
                args.land_proxy_multiplier, ocean_bulk_flux,
            )
            proxy_elapsed = time.perf_counter() - proxy_start
            if args.spatial_era5_rainout:
                rainout_days = np.full(source.shape, args.rain_days, dtype=float)
                absolute_lat = np.abs(candidate_fields.lat_deg)
                for index, (lower, upper) in enumerate(
                    zip(BAND_EDGES_DEG[:-1], BAND_EDGES_DEG[1:])
                ):
                    in_band = (absolute_lat >= lower) & (absolute_lat < upper)
                    rainout_days[in_band] = np.where(
                        candidate_fields.is_ocean[in_band],
                        ERA5_CLASS_RAINOUT_DAYS["ocean"][index],
                        ERA5_CLASS_RAINOUT_DAYS["land"][index],
                    )
                sink = _solver_order_rainout_sink(rainout_days)
            else:
                sink = np.full(source.size, 1.0 / (args.rain_days * 86_400.0))
            lift_rate_start = time.perf_counter()
            if args.orographic_lift:
                lift_rate_map, lift_metrics = _orographic_lift_rate(
                    candidate_fields, area,
                    wind_u=candidate_fields.wind_u * args.wind_multiplier,
                    wind_v=candidate_fields.wind_v * args.wind_multiplier,
                )
                lift_sink_solver = np.flip(lift_rate_map, axis=0).ravel()
                sink = sink + lift_sink_solver
            lift_rate_elapsed = time.perf_counter() - lift_rate_start
            u = candidate_fields.wind_u * args.wind_multiplier
            v = candidate_fields.wind_v * args.wind_multiplier
            solver_lat = lat[::-1]
            solver_edges = lat_edges[::-1]
            solver_area = np.flip(area, axis=0)
            solver_u = np.flip(u, axis=0)
            solver_v = np.flip(v, axis=0)
            solver_source = np.flip(source, axis=0).ravel()
            rate_start = time.perf_counter()
            rates = transition_rates(
                solver_u, solver_v, solver_lat, solver_edges, solver_area, finite_volume
            )
            rate_elapsed = time.perf_counter() - rate_start
            solve_start = time.perf_counter()
            if args.saturation_convection:
                surface_temperature = np.where(
                    candidate_fields.is_ocean, candidate_fields.ocean_temperature_c,
                    candidate_fields.air_temperature_c,
                )
                saturation_q = climate.saturation_specific_humidity_g_kg(
                    surface_temperature
                ) / 1000.0
                saturation_water = (
                    args.relative_humidity_threshold * AIR_DENSITY_KG_M3
                    * args.effective_depth_km * 1000.0 * saturation_q
                )
                mass, rain, solve_metrics = _solve_saturation_convection(
                    solver_source, solver_area.ravel(),
                    np.flip(saturation_water, axis=0).ravel(),
                    args.rain_days * 86_400.0,
                    args.convection_days * 86_400.0, rates,
                )
                closure = solve_metrics["global_closure_relative"]
                saturation_metrics = solve_metrics
            else:
                mass, rain, closure = solve_stationary(solver_source, sink, rates)
            solve_elapsed = time.perf_counter() - solve_start
            mass_grid = np.flip(mass.reshape((height, width)), axis=0)
            candidate_rain = (
                np.flip(rain.reshape((height, width)), axis=0) / area * SECONDS_PER_YEAR
            )
            current_temperature = np.where(
                candidate_fields.is_ocean, candidate_fields.ocean_temperature_c,
                candidate_fields.air_temperature_c,
            )
            slope = biomes.grid_slope(candidate_fields.elevation_m, candidate_fields.lat_deg)
            biome_start = time.perf_counter()
            candidate_biomes = biomes.smooth_biome_field(
                current_temperature, candidate_rain, candidate_fields.elevation_m, slope,
                candidate_fields.is_ocean, world.sea_level_m, lat_deg=candidate_fields.lat_deg,
                axial_tilt_deg=world.axial_tilt_deg,
            )
            candidate_biome_elapsed = time.perf_counter() - biome_start
            update_elapsed = time.perf_counter() - update_start
            if repeat >= WARMUP_COUNT:
                solver_times.append(solve_elapsed)
                candidate_update_times.append(update_elapsed)
                component_samples.append({
                    "climate_without_moisture_seconds": climate_elapsed,
                    "zero_rain_biome_classification_seconds": zero_rain_biome_elapsed,
                    "proxy_and_source_seconds": proxy_elapsed,
                    "duplicate_surface_sampling_seconds": duplicate_surface_samples[-1],
                    "source_proxy_work_excluding_resample_seconds": max(
                        0.0, proxy_elapsed - duplicate_surface_samples[-1]
                    ),
                    "lift_rate_build_seconds": lift_rate_elapsed,
                    "transition_rate_build_seconds": rate_elapsed,
                    "linear_solve_and_ledger_seconds": solve_elapsed,
                    "candidate_biome_classification_seconds": candidate_biome_elapsed,
                    "postsolve_work_excluding_biome_seconds": max(
                        0.0, update_elapsed - climate_elapsed - proxy_elapsed
                        - lift_rate_elapsed - rate_elapsed - solve_elapsed
                        - candidate_biome_elapsed
                    ),
                    # Remove the discarded zero-rain biome pass and repeat surface sample. #913
                    "production_shaped_update_seconds": max(
                        0.0, update_elapsed - zero_rain_biome_elapsed
                        - duplicate_surface_samples[-1]
                    ),
                    "full_candidate_update_seconds": update_elapsed,
                })
        assert mass is not None and rain is not None and closure is not None
        assert candidate_fields is not None and source is not None and area is not None
        assert tagged_sources is not None and u is not None and v is not None
        solver_seconds = float(np.median(solver_times))
        candidate_metrics = _rain_metrics(
            candidate_fields, candidate_rain, area, biomes, biome_ids=candidate_biomes
        )
        candidate_extended = _extended_rain_diagnostics(
            candidate_fields, candidate_rain, area, windward, leeward,
            candidate_metrics["koppen_group_shares"]["B"],
        )
        scale_comparison = {
            "native": _rms_pair(current.precipitation_mm, candidate_rain, area),
        }
        if density == 4:
            common_shape = (height // 4, width // 4)
            current_common, common_area = _area_average(
                current.precipitation_mm, area, common_shape
            )
            candidate_common, candidate_common_area = _area_average(
                candidate_rain, area, common_shape
            )
            if not np.array_equal(common_area, candidate_common_area):
                raise AssertionError("density-4 maps were not reduced with the same area weights")
            scale_comparison["area_average_to_density_1"] = _rms_pair(
                current_common, candidate_common, common_area
            )
        map_name = None
        if args.map_dir is not None:
            map_name = _save_common_scale_rain_map(
                args.map_dir, save, density, wind_mode, solver_name,
                "lift" if args.orographic_lift else "uniform",
                current.precipitation_mm, candidate_rain, area,
            )
        baseline_rain = current.precipitation_mm
        baseline_norm = max(float(np.sum(area * baseline_rain**2)), 1e-300)
        rel_rms = float(np.sqrt(np.sum(area * (candidate_rain - baseline_rain) ** 2)
                                / baseline_norm))
        bestfit_scale = float(
            np.sum(area * candidate_rain * baseline_rain)
            / max(float(np.sum(area * candidate_rain**2)), 1e-300)
        )
        bestfit_rms = float(np.sqrt(np.sum(
            area * (bestfit_scale * candidate_rain - baseline_rain) ** 2
        ) / baseline_norm))
        source_s = {name: float(values.sum()) for name, values in tagged_sources.items()}
        wind_flux_u = mass_grid * u
        wind_flux_v = mass_grid * v
        scheme_results[solver_name] = {
            "metrics": candidate_metrics,
            "extended_metrics": candidate_extended,
            "scale_comparison": scale_comparison,
            "common_scale_map": map_name,
            "precipitation_rel_rms_to_current": rel_rms,
            "bestfit_uniform_precipitation_scale": bestfit_scale,
            "bestfit_effective_depth_km": args.effective_depth_km * bestfit_scale,
            "bestfit_precipitation_rel_rms": bestfit_rms,
            "closure": closure,
            "saturation_convection_metrics": saturation_metrics,
            "orographic_lift_metrics": lift_metrics,
            "lift_rain_share": (
                None if lift_sink_solver is None
                else float(np.sum(lift_sink_solver * mass) / max(float(rain.sum()), 1e-300))
            ),
            "minimum_column_water_kg_m2": float(np.min(mass_grid / area)),
            "source_kg_s_by_tag": source_s,
            "source_kg_s_total": float(source.sum()),
            "rain_kg_s_total": float(rain.sum()),
            "wind_transport_area_mean_kg_m_s": {
                "eastward": _weighted_mean(wind_flux_u, area),
                "northward": _weighted_mean(wind_flux_v, area),
            },
            "solver_seconds": solver_seconds,
            "climate_update_median_seconds": float(np.median(candidate_update_times)),
            "current_source_update_median_seconds": float(np.median(current_times)),
            "candidate_update_components_median_seconds": {
                name: float(np.median([sample[name] for sample in component_samples]))
                for name in component_samples[-1]
            },
        }
    return {
        "save": save.name,
        "save_sha256": _sha256(save),
        "density": density,
        "shape": [height, width],
        "lagged_climate_cache_resampled": cache_resampled,
        "requested_wind_mode": requested_wind,
        "effective_wind_mode": "cfd" if wind_available and requested_wind == "cfd" else "diagnostic",
        "solver_parameters": {
            "wind_multiplier": args.wind_multiplier,
            "effective_depth_km": args.effective_depth_km,
            "source_timescale_days": 9.1,
            "rainout_days": args.rain_days,
            "orographic_lift": args.orographic_lift,
            "land_proxy_multiplier": args.land_proxy_multiplier,
            "ocean_source_mapping": args.ocean_source_mapping,
            "ocean_transfer_coefficient": args.ocean_transfer_coefficient,
            "spatial_era5_rainout": args.spatial_era5_rainout,
            "saturation_convection": args.saturation_convection,
            "relative_humidity_threshold": args.relative_humidity_threshold,
            "convection_days": args.convection_days,
        },
        "current_source_metrics": current_metrics,
        "current_source_extended_metrics": current_extended,
        "solvers": scheme_results,
        "windward_leeward_mask_definition": exposure,
        "peak_rss_mb_process": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6,
    }


def main() -> None:
    _check_asymmetric_rainout_fixture()
    args = _arguments()
    if args.lift_comparison and args.orographic_lift:
        raise ValueError("choose --lift-comparison or --orographic-lift, not both")
    if (args.orographic_lift or args.lift_comparison) and (
        args.saturation_convection or args.spatial_era5_rainout
    ):
        raise ValueError("orographic lift arm requires the uniform background sink")
    backend = Path(__file__).resolve().parents[2] / "backend"
    sys.path.insert(0, str(backend))
    from app import persistence
    from app.hydroclimate import climate, hydrology, biomes

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=backend, check=True,
                            capture_output=True, text=True).stdout.strip()
    densities = (1, 4) if args.density == "both" else (int(args.density),)
    winds = ("diagnostic", "cfd") if args.wind == "both" else (args.wind,)
    solvers = (("finite_volume", "wind_arrows") if args.solver == "both"
               else (args.solver,))
    output = {
        "rfc": "0019",
        "gate": "G3 tuning comparison",
        "commit": commit,
        "parameters": vars(args) | {"saves": [str(path) for path in args.saves],
                                    "out": str(args.out),
                                    "map_dir": None if args.map_dir is None else str(args.map_dir)},
        "results": [],
    }
    for save in args.saves:
        data = save.read_bytes()
        world = persistence.load_world_bytes(data)
        saved_density = world.climate_density
        cache_shape = None if world.climate_cache is None else list(world.climate_cache.precipitation_mm.shape)
        original_cache = world.climate_cache
        for density in densities:
            for wind in winds:
                if args.sweep_tuning:
                    settings = product((1.5, 2.25, 3.0), (5.0, 9.1, 18.0), (0.5, 1.0, 2.0))
                else:
                    settings = [(args.effective_depth_km, args.rain_days, args.wind_multiplier)]
                coefficients = (0.001, 0.003, 0.005) if args.sweep_ocean_coefficient else (
                    args.ocean_transfer_coefficient,
                )
                land_multipliers = (1.0, 5.0, 10.0) if args.sweep_land_proxy else (
                    args.land_proxy_multiplier,
                )
                for depth, rain_days, wind_multiplier in settings:
                    for ocean_coefficient in coefficients:
                        for land_multiplier in land_multipliers:
                            arms = (("uniform", "lift_aware") if args.lift_comparison else
                                    (("lift_aware",) if args.orographic_lift else ("uniform",)))
                            for arm in arms:
                                run_args = copy(args)
                                run_args.effective_depth_km = depth
                                run_args.rain_days = rain_days
                                run_args.wind_multiplier = wind_multiplier
                                run_args.ocean_transfer_coefficient = ocean_coefficient
                                run_args.land_proxy_multiplier = land_multiplier
                                run_args.orographic_lift = arm == "lift_aware"
                                result = _compare_one(world, save, density, wind, solvers, run_args,
                                                      climate, hydrology, biomes, original_cache)
                                result["rainout_arm"] = arm
                                result["saved_climate_density"] = saved_density
                                result["loaded_climate_cache_shape"] = cache_shape
                                output["results"].append(result)
                                args.out.parent.mkdir(parents=True, exist_ok=True)
                                args.out.write_text(json.dumps(output, indent=2))
                                print(json.dumps({key: result[key] for key in (
                                    "save", "density", "requested_wind_mode", "peak_rss_mb_process")
                                } | {"rainout_arm": arm, "effective_depth_km": depth,
                                     "rain_days": rain_days, "wind_multiplier": wind_multiplier,
                                     "ocean_transfer_coefficient": ocean_coefficient,
                                     "land_proxy_multiplier": land_multiplier}), flush=True)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
