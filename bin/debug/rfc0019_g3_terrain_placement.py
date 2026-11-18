"""Test wind-aligned fine rain placement using the stored terrain on the pinned tuning saves.

Usage: backend/.venv/bin/python bin/debug/rfc0019_g3_terrain_placement.py SAVE BOUNDED_CAPTURE
Only the cached elevation is read from the save. The matched coarse climate supplies
the rain budget, wind, and surface classes. This is not a density-4 climate update.
"""

from __future__ import annotations

import argparse
import json
import sys
import zipfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.ndimage import gaussian_filter1d, zoom

import rfc0019_g3_compare as g3
from rfc0019_coarse_moisture import aggregate_inputs, distribute_rain
from rfc0019_g3_archive_cap import _compare_projection
from rfc0019_g3_bounded_column import _downslope_pairs
from rfc0019_g3_mixing import GRID_SHAPE, TUNING_SAVES

FACTOR = 4
FINE_SHAPE = tuple(axis * FACTOR for axis in GRID_SHAPE)
SMOOTH_M = 250_000.0
TERRAIN_WEIGHT = 0.35


def _elevation(save: Path) -> np.ndarray:
    with zipfile.ZipFile(save) as archive:
        manifest = json.loads(archive.read("manifest.json"))
        elevation_key = manifest["world"]["climate_cache"]["fields"]["elevation_m"]["__npy__"]
        with archive.open(f"a/{elevation_key}.npy") as stored:
            native = np.load(stored, allow_pickle=False)
    _, _, area = g3._grid_area(*native.shape)
    return g3._area_average(native, area, FINE_SHAPE)[0]


def _upscale_wind(coarse: np.ndarray) -> np.ndarray:
    padded = np.pad(coarse, ((0, 0), (1, 1)), mode="wrap")
    return zoom(padded, FACTOR, order=1, mode="nearest", grid_mode=True)[:, FACTOR:-FACTOR]


def _signed_slope(fields) -> np.ndarray:
    height, width = fields.elevation_m.shape
    latitude = np.deg2rad(fields.lat_deg)
    meridional_spacing = g3.EARTH_RADIUS_M * np.pi / height
    zonal_spacing = 2 * np.pi * g3.EARTH_RADIUS_M * np.cos(latitude) / width
    elevation = gaussian_filter1d(fields.elevation_m, SMOOTH_M / meridional_spacing, axis=0, mode="nearest")
    for index, spacing in enumerate(zonal_spacing):
        elevation[index] = gaussian_filter1d(elevation[index], SMOOTH_M / spacing, mode="wrap")
    north = np.gradient(elevation, latitude, axis=0) / g3.EARTH_RADIUS_M
    east = (np.roll(elevation, -1, axis=1) - np.roll(elevation, 1, axis=1)) / (2 * zonal_spacing[:, None])
    speed = np.hypot(fields.wind_u, fields.wind_v)
    slope = fields.wind_u * east + fields.wind_v * north
    return np.divide(slope, speed, out=np.zeros_like(slope), where=speed > 0)


def _terrain_weights(fields, area: np.ndarray) -> np.ndarray:
    slopes = _signed_slope(fields)
    land = ~fields.is_ocean
    rms = np.sqrt(np.average(slopes[land]**2, weights=area[land]))
    if rms == 0:
        return np.ones_like(area)
    return 1.0 + TERRAIN_WEIGHT * np.where(land, np.tanh(slopes / rms), 0.0)


def _pair_metrics(uniform, lift, up, down, include):
    uniform_delta = (uniform[up] - uniform[down])[include]
    lift_delta = (lift[up] - lift[down])[include]
    if not lift_delta.size:
        return {"pairs": 0}
    mean_rain = (0.5 * (uniform[up] + uniform[down]))[include]
    significant = 1e-12 * max(1.0, float(np.max(mean_rain)))
    return {
        "pairs": int(lift_delta.size),
        "uniform_pair_median_rain_mm_yr": float(np.median(mean_rain)),
        "uniform_median_mm_yr": float(np.median(uniform_delta)),
        "terrain_median_mm_yr": float(np.median(lift_delta)),
        "uniform_windward_wetter_share": float(np.mean(uniform_delta > significant)),
        "terrain_windward_wetter_share": float(np.mean(lift_delta > significant)),
    }


def _place_rain(fields, elevation, coarse_rain):
    fine_sea = fields.is_ocean.repeat(FACTOR, axis=0).repeat(FACTOR, axis=1)
    grid = aggregate_inputs(fine_sea[::-1], {}, factor=FACTOR)
    fine = SimpleNamespace(
        is_ocean=fine_sea, elevation_m=elevation, lat_deg=np.linspace(89.75, -89.75, FINE_SHAPE[0]),
        wind_u=_upscale_wind(fields.wind_u), wind_v=_upscale_wind(fields.wind_v),
    )
    area = grid.native_area[::-1]
    budgets = coarse_rain[::-1].ravel() * grid.area / g3.SECONDS_PER_YEAR
    uniform = distribute_rain(grid, budgets, np.ones(FINE_SHAPE))
    weighted = distribute_rain(grid, budgets, _terrain_weights(fine, area)[::-1])
    max_residual = 0.0
    for deposition in (uniform, weighted):
        deposited = np.bincount(grid.native_ids.ravel(), weights=deposition.ravel())
        np.testing.assert_allclose(deposited, budgets, rtol=1e-12, atol=0.0)
        max_residual = max(max_residual, float(np.max(np.abs(deposited - budgets))))
    uniform_rate = uniform[::-1] / area * g3.SECONDS_PER_YEAR
    weighted_rate = weighted[::-1] / area * g3.SECONDS_PER_YEAR
    up, down = _downslope_pairs(fine, g3._windward_leeward_masks(fine, area))
    ids = grid.native_ids[::-1]
    within = ids[up] == ids[down]
    return {
        "within_one_budget": _pair_metrics(uniform_rate, weighted_rate, up, down, within),
        "across_budgets": _pair_metrics(uniform_rate, weighted_rate, up, down, ~within),
        "maximum_reservoir_residual_of_total_rain": max_residual / float(budgets.sum()),
    }


def _sample(capture, elevation, climate, freeze_point):
    latitude, _, xyz = climate._build_grid(*GRID_SHAPE)
    result = {}
    for mode in ("diagnostic", "cfd"):
        with np.load(capture / f"{mode}-columns.npz", allow_pickle=False) as archive:
            fields = SimpleNamespace(
                lat_deg=latitude, world_xyz=xyz, is_ocean=archive["is_ocean"],
                air_temperature_c=archive["air_temperature_c"], elevation_m=archive["elevation_m"],
                wind_u=archive["wind_u_m_s"], wind_v=archive["wind_v_m_s"],
            )
            coarse = aggregate_inputs(fields.is_ocean[::-1], {})
            _, maps = _compare_projection(
                fields, coarse, archive["source_kg_s"], archive["control_mass_kg"],
                archive["reference_precipitation_mm_yr"], climate, freeze_point, screen="placement",
            )
        result[mode] = _place_rain(fields, elevation, maps["broad_transport_half_spread"])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("save", type=Path)
    parser.add_argument("capture", type=Path)
    args = parser.parse_args()
    provenance = json.loads((args.capture / "results.json").read_text())
    digest = TUNING_SAVES.get(args.save.name)
    if not digest or not digest == provenance["save_sha256"] == g3._sha256(args.save):
        parser.error("expected the pinned save matching the captured climate")

    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root / "backend"))
    from app.hydroclimate import climate, hydrology

    report = {"save_sha256": digest, "fine_shape": FINE_SHAPE, "smooth_m": SMOOTH_M, "terrain_weight": TERRAIN_WEIGHT}
    report["results"] = _sample(args.capture, _elevation(args.save), climate, hydrology.FREEZE_POINT_C)
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
