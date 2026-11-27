"""Build a G3b calibration packet from an already-exposed density-4 capture.

Usage: backend/.venv/bin/python bin/debug/rfc0019_g3b_calibration_packet.py CAPTURE OUT

Re-solves the frozen ridge_025 and paired control for both transport solvers and
wind fields. Produces native-grid candidate arrays, fixed linear-scale maps,
co-registered geography/wind maps, latitude profiles and geography-stratified
coastal comparisons. No world is loaded or stepped. This is a calibration tool,
not a frozen evaluation protocol; do not use on fresh evaluation captures.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import numpy as np
from PIL import Image, ImageDraw

import rfc0019_g3_ridge_rainout as ridge

SEED = 'seed910211954-221.6myr.bloom'
SAVE_SHA256 = '1f167a8d9e0edf06c3c52b23598dc2f776adc2927f843beb992b5193ad30fd67'
SOLVERS = ('finite_volume', 'wind_arrows')
WINDS = ('diagnostic', 'cfd')


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_capture(capture):
    meta = json.loads((capture / 'results.json').read_text())
    if (Path(meta['save']).name != SEED or meta.get('heldout') is not False
            or meta.get('save_sha256') != SAVE_SHA256 or meta.get('density') != 4):
        raise ValueError('only the already-exposed pinned development capture is permitted')
    return meta


def land_exposure(ocean, u, v, latitude):
    """Classify land by an upwind ray to 500 km, using winds and geography only.

    Cells whose ray crosses a pole or whose wind is calm are unclassified. Sampling
    is 50 km along a locally straight bearing on a sphere, wrapping longitude.
    This is an exposure proxy, not an orographic rain-shadow test.
    """
    h, w = ocean.shape
    rows, cols = np.indices((h, w))
    speed = np.hypot(u, v)
    valid = speed > 1e-6
    ocean_seen = np.zeros_like(ocean)
    active = valid & ~ocean
    lat = latitude[:, None]
    for km in range(50, 501, 50):
        dr = np.rint((v / np.maximum(speed, 1e-30)) * km / 111.195 * h / 180).astype(int)
        coslat = np.maximum(np.cos(np.radians(lat)), 1e-6)
        dc = np.rint((-u / np.maximum(speed, 1e-30)) * km / (111.195 * coslat) * w / 360).astype(int)
        rr = rows + dr
        in_bounds = (rr >= 0) & (rr < h)
        active &= in_bounds
        ocean_seen |= active & ocean[np.clip(rr, 0, h - 1), (cols + dc) % w]
    # -1 unknown, 0 no ocean along ray, 1 ocean along ray
    return np.where(active, ocean_seen.astype(np.int8), -1)


def weighted_mean(values, area, mask):
    total = float(area[mask].sum())
    return {'area_km2': total / 1e6, 'cells': int(mask.sum()),
            'rain_mm_yr': float(np.sum(values[mask] * area[mask]) / total) if total else None}


def summaries(rain, area, ocean, distance, exposure, lat):
    lat2 = np.broadcast_to(lat[:, None], rain.shape)
    profiles = []
    for lo in range(-90, 90, 5):
        band = (lat2 >= lo) & (lat2 < lo + 5)
        profiles.append({'latitude_deg': [lo, lo + 5], **{
            name: weighted_mean(rain, area, band & mask)
            for name, mask in [('all', np.ones_like(ocean)), ('land', ~ocean), ('ocean', ocean)]}})
    comparisons = []
    for lo in range(-90, 90, 15):
        band = (lat2 >= lo) & (lat2 < lo + 15) & ~ocean
        for exp, label in [(1, 'ocean_upwind_500km'), (0, 'no_ocean_upwind_500km'), (-1, 'unknown')]:
            for low, high, coast in [(0, 500, 'near'), (500, 1500, 'middle'), (1500, np.inf, 'interior')]:
                mask = band & (exposure == exp) & (distance >= low) & (distance < high)
                comparisons.append({'latitude_deg': [lo, lo + 15], 'exposure': label,
                                    'distance': coast, **weighted_mean(rain, area, mask)})
    return profiles, comparisons


def rain_colors(rain, maximum):
    palette = np.array([[68, 1, 84], [59, 82, 139], [33, 145, 140], [94, 201, 98], [253, 231, 37]])
    pos = np.clip(rain / maximum, 0, 1) * (len(palette) - 1)
    low = np.floor(pos).astype(int)
    color = palette[low] * (1 - (pos - low)[..., None]) + palette[np.minimum(low + 1, 4)] * (pos - low)[..., None]
    color[rain > maximum] = (255, 0, 255)
    return color.astype(np.uint8)


def terrain_colors(ocean, elevation):
    land_height = np.clip(elevation / 4000, 0, 1)
    terrain = np.empty((*ocean.shape, 3), dtype=np.uint8)
    terrain[:] = (180, 210, 232)
    terrain[~ocean] = np.stack((110 + 140 * land_height, 135 + 100 * land_height,
                                 85 + 130 * land_height), axis=-1)[~ocean].astype(np.uint8)
    return terrain


def draw_map(path, rain, ocean, elevation, u, v, maximum, label):
    h, w = ocean.shape
    panel_h, panel_w = h, w
    image = Image.new('RGB', (2 * panel_w + 30, panel_h + 90), 'white')
    draw = ImageDraw.Draw(image)
    image.paste(Image.fromarray(rain_colors(rain, maximum)), (10, 35))
    x0 = panel_w + 20
    image.paste(Image.fromarray(terrain_colors(ocean, elevation)), (x0, 35))
    for row in range(15, h, 30):
        for col in range(15, w, 30):
            speed = np.hypot(u[row, col], v[row, col])
            if speed < 1e-6:
                continue
            dx, dy = 8 * u[row, col] / speed, -8 * v[row, col] / speed
            x, y = x0 + col, 35 + row
            draw.line((x - dx, y - dy, x + dx, y + dy), fill='black', width=1)
            draw.ellipse((x + dx - 2, y + dy - 2, x + dx + 2, y + dy + 2), fill='black')
    draw.text((10, 5), label + ' | rain', fill='black')
    draw.text((x0, 5), 'Terrain and winds (dot is downwind)', fill='black')
    gradient = np.tile(np.linspace(0, maximum, 256), (10, 1))
    image.paste(Image.fromarray(rain_colors(gradient, maximum)), (10, panel_h + 40))
    draw.text((275, panel_h + 39), f'0 to {maximum:.1f} mm/yr (linear); magenta > maximum', fill='black')
    draw.text((10, panel_h + 58), 'North up; longitude -180 to +180. Ocean blue; land 0 to 4000 m green to pale.', fill='black')
    image.save(path)


def draw_profiles(path, profiles, label):
    image = Image.new('RGB', (900, 500), 'white')
    draw = ImageDraw.Draw(image)
    left, top, width, height = 65, 45, 800, 380
    vals = [row[k]['rain_mm_yr'] for row in profiles for k in ('all', 'land', 'ocean')]
    ymax = max(x for x in vals if x is not None) * 1.05
    draw.rectangle((left, top, left + width, top + height), outline='black')
    for degrees in (-90, -60, -30, 0, 30, 60, 90):
        x = left + (degrees + 90) / 180 * width
        draw.line((x, top + height, x, top + height + 5), fill='black')
        draw.text((x - 8, top + height + 8), str(degrees), fill='black')
    for fraction in (0, 0.5, 1):
        y = top + height * (1 - fraction)
        draw.line((left - 5, y, left, y), fill='black')
        draw.text((4, y - 5), f'{ymax * fraction:.0f}', fill='black')
    for name, color in [('all', 'black'), ('land', 'saddlebrown'), ('ocean', 'blue')]:
        previous = None
        for i, row in enumerate(profiles):
            value = row[name]['rain_mm_yr']
            if value is None:
                previous = None
                continue
            point = (left + (i + .5) * width / len(profiles), top + height * (1 - value / ymax))
            if previous is not None:
                draw.line((*previous, *point), fill=color, width=2)
            previous = point
    draw.text((left, 10), label + ' | area-weighted rain, 5-degree latitude bands', fill='black')
    draw.text((left, 450), 'Latitude (degrees, south to north); vertical axis: mm/yr', fill='black')
    draw.text((left, 470), 'Black: all  Brown: land  Blue: ocean', fill='black')
    image.save(path)


def write_strata(path, comparisons, area):
    planet_area = float(area.sum()) / 1e6
    with path.open('w', newline='') as stream:
        names = ['lat_south', 'lat_north', 'exposure', 'distance', 'area_km2', 'planet_area_fraction',
                 'cells', 'rain_mm_yr', 'area_at_least_0_1_percent_planet']
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        for row in comparisons:
            fraction = row['area_km2'] / planet_area
            writer.writerow({'lat_south': row['latitude_deg'][0], 'lat_north': row['latitude_deg'][1],
                             'exposure': row['exposure'], 'distance': row['distance'],
                             'area_km2': row['area_km2'], 'planet_area_fraction': fraction,
                             'cells': row['cells'], 'rain_mm_yr': row['rain_mm_yr'],
                             'area_at_least_0_1_percent_planet': fraction >= 0.001})


def write_ridge_profiles(path, map_path, profiles, elevation, ocean):
    h, w = ocean.shape
    image = Image.fromarray(terrain_colors(ocean, elevation)).resize((w * 2, h * 2), Image.Resampling.NEAREST)
    draw = ImageDraw.Draw(image)
    with path.open('w', newline='') as stream:
        names = ['latitude_deg', 'longitude_deg', 'axis', 'downwind_step_solver_order',
                 'crest_coarse_cell', 'windward_mean_mm_yr', 'lee_mean_mm_yr', 'contrast_mm_yr']
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        for p in profiles:
            samples = p['samples_by_offset_native_cells']
            sides = [[samples[str(offset)] for offset in offsets]
                     for offsets in ((-3, -2, -1), (1, 2, 3))]
            values = [[s['rain_mm_yr'] for s in side if s is not None] for side in sides]
            if min(map(len, values)) < 2:
                continue
            windward, lee = (float(np.mean(side)) for side in values)
            row = h - 1 - p['native_peak_row_solver_order']
            col = p['native_peak_col']
            writer.writerow({'latitude_deg': 90 - (row + 0.5) * 180 / h,
                             'longitude_deg': -180 + (col + 0.5) * 360 / w,
                             'axis': p['axis'], 'downwind_step_solver_order': p['downwind_step'],
                             'crest_coarse_cell': p['crest_coarse_cell'],
                             'windward_mean_mm_yr': windward, 'lee_mean_mm_yr': lee,
                             'contrast_mm_yr': windward - lee})
            color = (200, 0, 0) if windward > lee else (0, 40, 200)
            x, y = col * 2, row * 2
            draw.ellipse((x - 2, y - 2, x + 2, y + 2), fill=color)
    image.save(map_path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('capture', type=Path)
    parser.add_argument('out', type=Path)
    args = parser.parse_args()
    meta = validate_capture(args.capture)
    args.out.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'backend'))
    from app.hydroclimate.hydrology import FREEZE_POINT_C

    archives = {w: np.load(args.capture / f'{w}-columns.npz', allow_pickle=False) for w in WINDS}
    positive = np.concatenate([a['reference_precipitation_mm_yr'].ravel() for a in archives.values()])
    maximum = float(np.percentile(positive[positive > 0], 99))
    result = {'capture': SEED, 'save_sha256': SAVE_SHA256, 'heldout': False, 'density': 4,
              'code_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
              'generator_sha256': sha256(Path(__file__)),
              'source_sha256': {name: sha256(Path(name)) for name in ridge.SOURCE_FILES},
              'capture_sha256': {f'{w}-columns.npz': sha256(args.capture / f'{w}-columns.npz') for w in WINDS},
              'scale': {'type': 'linear', 'maximum_mm_yr': maximum,
                        'basis': 'unweighted pooled p99 of positive current-climate native cells, both winds'},
              'projection': 'equirectangular native grid, north up, longitude -180 to +180',
              'exposure': 'upwind ray sampled every 50 km to 500 km; ocean encountered or not; calm/pole crossing unknown',
              'capture_grid': meta['grid_shape'], 'results': {}}
    for wind, a in archives.items():
        lat = a['lat_deg_north_to_south']
        fields = SimpleNamespace(is_ocean=a['is_ocean'], lat_deg=lat, air_temperature_c=a['air_temperature_c'],
                                 elevation_m=a['elevation_m'], wind_u=a['wind_u_m_s'], wind_v=a['wind_v_m_s'])
        exposure = land_exposure(fields.is_ocean, fields.wind_u, fields.wind_v, lat)
        draw_map(args.out / f'{wind}-reference-map.png', a['reference_precipitation_mm_yr'],
                 fields.is_ocean, fields.elevation_m, fields.wind_u, fields.wind_v, maximum,
                 f'{wind} current climate (reference)')
        result['results'][wind] = {}
        for solver in SOLVERS:
            print(f'Solving {wind} {solver}', flush=True)
            metrics, maps, area = ridge._solve_arms(fields, a, FREEZE_POINT_C,
                                                    {'control': 0.0, 'ridge_025': 0.25}, solver=solver)
            rain = maps['ridge_025']
            stem = f'{wind}-{solver}'
            np.savez_compressed(args.out / f'{stem}-arrays.npz', rain_mm_yr=rain,
                                reference_mm_yr=a['reference_precipitation_mm_yr'],
                                area_m2=area, exposure=exposure)
            draw_map(args.out / f'{stem}-map.png', rain, fields.is_ocean, fields.elevation_m,
                     fields.wind_u, fields.wind_v, maximum, stem)
            draw_map(args.out / f'{stem}-unsaturated.png', rain, fields.is_ocean, fields.elevation_m,
                     fields.wind_u, fields.wind_v, float(np.max(rain)), stem + ' diagnostic scale')
            profiles, comparisons = summaries(rain, area, fields.is_ocean, a['coast_distance_km'], exposure, lat)
            draw_profiles(args.out / f'{stem}-latitude.png', profiles, stem)
            write_strata(args.out / f'{stem}-coast-strata.csv', comparisons, area)
            row = metrics['ridge_025']
            write_ridge_profiles(args.out / f'{stem}-ridge-profiles.csv', args.out / f'{stem}-ridges.png',
                                 row['ridge_native_cell_profiles']['profiles'], fields.elevation_m,
                                 fields.is_ocean)
            result['results'][wind][solver] = {'latitude_profiles': profiles, 'coast_strata': comparisons,
                'geographic_support': row['geographic_support'], 'ridge_transects': row['ridge_transects'],
                'ridge_native_cell_profiles': row['ridge_native_cell_profiles'], 'hard_gates': row['hard_gates'],
                'global_closure': row['global_closure'], 'summed_cell_closure': row['summed_cell_closure'],
                'over_scale_cells': int(np.count_nonzero(rain > maximum))}
            print(f'Finished {stem}', flush=True)
    (args.out / 'packet.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    for a in archives.values():
        a.close()


if __name__ == '__main__':
    main()
