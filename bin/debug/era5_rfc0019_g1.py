#!/usr/bin/env python3
"""RFC 0019 G1: sample simultaneous ERA5 variables on the WeatherBench2 1.5° grid.

Usage: /path/to/python bin/debug/era5_rfc0019_g1.py --out analysis/rfc0017/g1-era5.json
Requires numpy and numcodecs; reads public ERA5 Zarr chunks via HTTPS. Samples two days
near the start of each month of 2000, not a climatology. Flux is a 13-pressure-level
approximation, with surface extrapolation; this is not ERA5's native integrated flux.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import urllib.request
from pathlib import Path

import numpy as np
from numcodecs import Blosc

STORE = 'https://storage.googleapis.com/weatherbench2/datasets/era5/1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr'
SURFACE = ('2m_dewpoint_temperature', '2m_temperature', 'surface_pressure', 'total_column_water_vapour', 'total_precipitation_6hr', '10m_u_component_of_wind', '10m_v_component_of_wind')
LEVEL_FIELDS = ('specific_humidity', 'u_component_of_wind', 'v_component_of_wind')
G = 9.80665


def fetch(path: str) -> bytes:
    with urllib.request.urlopen(f'{STORE}/{path}', timeout=120) as response:
        return response.read()


def meta(name: str) -> dict:
    return json.loads(fetch(f'{name}/.zarray'))


def read(name: str, key: str, metadata: dict, hashes: dict) -> np.ndarray:
    raw = fetch(f'{name}/{key}')
    hashes[f'{name}/{key}'] = hashlib.sha256(raw).hexdigest()
    return np.frombuffer(Blosc().decode(raw), dtype=np.dtype(metadata['dtype'])).reshape(metadata['chunks'])


def humidity_and_density(dewpoint, temperature, pressure):
    # Magnus water-vapour pressure over liquid water; approximation below freezing.
    dew_c = dewpoint - 273.15
    vapour_pa = 611.2 * np.exp(17.67 * dew_c / (dew_c + 243.5))
    humidity = 0.622 * vapour_pa / (pressure - 0.378 * vapour_pa)
    density = pressure / (287.05 * temperature * (1.0 + 0.61 * humidity))
    return humidity, density


def integrated_flux(q, u, v, levels_pa, pressure, surface_q, surface_u, surface_v):
    """Trapezoidal pressure integration, excluding below-surface pressure levels."""
    last_p = np.zeros_like(pressure, dtype=np.float64)
    last_x = (q[:, 0] * u[:, 0]).astype(np.float64)
    last_y = (q[:, 0] * v[:, 0]).astype(np.float64)
    last_q = q[:, 0].astype(np.float64)
    integral_q = np.zeros_like(last_p)
    integral_x = np.zeros_like(last_p)
    integral_y = np.zeros_like(last_p)
    for i, level in enumerate(levels_pa):
        valid = (pressure >= level) & np.isfinite(q[:, i]) & np.isfinite(u[:, i]) & np.isfinite(v[:, i])
        x = q[:, i] * u[:, i]
        y = q[:, i] * v[:, i]
        dp = np.where(valid, level - last_p, 0.0)
        integral_x += np.where(valid, dp * (last_x + x) * 0.5, 0.0)
        integral_y += np.where(valid, dp * (last_y + y) * 0.5, 0.0)
        integral_q += np.where(valid, dp * (last_q + q[:, i]) * 0.5, 0.0)
        last_p = np.where(valid, level, last_p)
        last_x = np.where(valid, x, last_x)
        last_y = np.where(valid, y, last_y)
        last_q = np.where(valid, q[:, i], last_q)
    dp = np.maximum(pressure - last_p, 0.0)
    integral_x += dp * (last_x + surface_q * surface_u) * 0.5
    integral_y += dp * (last_y + surface_q * surface_v) * 0.5
    integral_q += dp * (last_q + surface_q) * 0.5
    return integral_x / G, integral_y / G, integral_q / G


def weighted_quantiles(values, weights, quantiles=(0.1, 0.5, 0.9)):
    if values.size == 0:
        return None
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return [float(x) for x in np.interp(np.asarray(quantiles) * cumulative[-1], cumulative, values[order])]


def summarize(fields, area, mask):
    selected = mask & np.isfinite(fields['height']) & np.isfinite(fields['flux_x']) & np.isfinite(fields['flux_y']) & (fields['proxy_speed'] > 5.0)
    weights = area[selected]
    def mean(name, scale=1.0):
        valid = mask & np.isfinite(fields[name])
        return float(np.average(fields[name][valid], weights=area[valid]) * scale) if valid.any() else None
    dot = fields['proxy_x'] * fields['flux_x'] + fields['proxy_y'] * fields['flux_y']
    squared = fields['proxy_x'] ** 2 + fields['proxy_y'] ** 2
    ratio = dot[selected] / squared[selected]
    return {
        'area_fraction': float(area[mask].sum() / area.sum()),
        'column_water_kg_m2': mean('water'),
        'pressure_integrated_water_kg_m2': mean('integrated_water'),
        'surface_specific_humidity_g_kg': mean('humidity', 1000),
        'surface_wind_vector_m_s': [mean('wind_x'), mean('wind_y')],
        'precipitation_mm_day': mean('precip'),
        'ivt_vector_kg_m_s': [mean('flux_x'), mean('flux_y')],
        'effective_depth_m_quantiles_10_50_90': weighted_quantiles(fields['height'][selected], weights),
        'parallel_transport_ratio_quantiles_10_50_90': weighted_quantiles(ratio, weights),
        'signed_vector_regression_multiplier': (
            float(np.sum(area[selected] * dot[selected]) / np.sum(area[selected] * squared[selected]))
            if selected.any() else None
        ),
        'selected_fraction_of_class_area': float(weights.sum() / area[mask].sum()) if area[mask].sum() else 0.0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--months', type=int, nargs='+', default=list(range(1, 13)))
    args = parser.parse_args()
    names = (*SURFACE, *LEVEL_FIELDS, 'latitude', 'level', 'land_sea_mask')
    metadata = {name: meta(name) for name in names}
    hashes = {}
    lat = read('latitude', '0', metadata['latitude'], hashes)
    levels = read('level', '0', metadata['level'], hashes).astype(float) * 100.0
    land_fraction = read('land_sea_mask', '0.0', metadata['land_sea_mask'], hashes)
    weights = np.broadcast_to(np.cos(np.radians(lat))[None, :], land_fraction.shape).copy()
    fields = {name: np.zeros(land_fraction.shape, dtype=np.float64) for name in ('water','humidity','wind_x','wind_y','precip','proxy_x','proxy_y','flux_x','flux_y','rho_q','integrated_water')}
    starts = []
    for month in args.months:
        first = dt.datetime(2000, month, 1)
        steps = int((first - dt.datetime(1959, 1, 1)).total_seconds() / 21600)
        chunk = steps // 8
        starts.append((dt.datetime(1959, 1, 1) + dt.timedelta(hours=chunk * 48)).isoformat() + 'Z')
        surface = {name: read(name, f'{chunk}.0.0', metadata[name], hashes).astype(np.float64) for name in SURFACE}
        q, u, v = (read(name, f'{chunk}.0.0.0', metadata[name], hashes) for name in LEVEL_FIELDS)
        humidity, density = humidity_and_density(surface[SURFACE[0]], surface[SURFACE[1]], surface[SURFACE[2]])
        fx, fy, integrated_water = integrated_flux(q, u, v, levels, surface[SURFACE[2]], humidity, surface[SURFACE[5]], surface[SURFACE[6]])
        values = {'water': surface[SURFACE[3]], 'humidity': humidity, 'wind_x': surface[SURFACE[5]], 'wind_y': surface[SURFACE[6]], 'precip': surface[SURFACE[4]] * 4000.0, 'proxy_x': surface[SURFACE[3]] * surface[SURFACE[5]], 'proxy_y': surface[SURFACE[3]] * surface[SURFACE[6]], 'flux_x': fx, 'flux_y': fy, 'rho_q': density * humidity, 'integrated_water': integrated_water}
        for name, value in values.items():
            fields[name] += value.mean(axis=0) / len(args.months)
        print(f'processed {starts[-1]}', flush=True)
    fields['height'] = fields['water'] / fields['rho_q']
    fields['proxy_speed'] = np.hypot(fields['proxy_x'], fields['proxy_y'])
    land = land_fraction >= 0.5
    lat_grid = np.broadcast_to(lat[None, :], land.shape)
    classes = {'global': np.ones_like(land, dtype=bool), 'land': land, 'ocean': ~land}
    for lo in range(0, 90, 15):
        band = (np.abs(lat_grid) >= lo) & (np.abs(lat_grid) < lo + 15 if lo < 75 else np.abs(lat_grid) <= 90)
        classes[f'abs_lat_{lo}_{lo+15}_land'] = band & land
        classes[f'abs_lat_{lo}_{lo+15}_ocean'] = band & ~land
    result = {'source': STORE, 'sampling': {'year': 2000, 'months': args.months, 'block_starts_utc': starts, 'samples_per_block': 8, 'spacing_hours': 6, 'grid_degrees': 1.5}, 'method': {'specific_humidity': 'Magnus water saturation from 2m dew point and surface pressure; below-freezing approximation', 'precipitation': '6-hour accumulation multiplied by 4, m to mm/day', 'ivt': '13 pressure levels 50-1000 hPa; trapezoidal q*u/g and q*v/g; below-ground levels excluded; top held constant to 0 hPa; surface interpolated from 2m humidity and 10m wind', 'depth': 'ratio of sample-mean column vapour to sample-mean near-surface vapour density, not measured mixing depth', 'transport': 'signed vector projection of mean pressure-integrated IVT onto mean instantaneous W times 10m vector; mean ratios use proxy norm >5 kg m-1 s-1; sample is not a climatology'}, 'compressed_chunk_sha256': hashes, 'classes': {name: summarize(fields, weights, mask) for name, mask in classes.items()}}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(f'wrote {args.out}', flush=True)


if __name__ == '__main__':
    main()
