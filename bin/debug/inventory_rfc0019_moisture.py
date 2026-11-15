#!/usr/bin/env python3
"""Inventory RFC 0019 G0 climate-grid proxies and precipitation on one saved world.

Usage: backend/.venv/bin/python bin/debug/inventory_rfc0019_moisture.py SAVE --out OUTPUT
The humidity proxies are area-integrated g/kg, not atmospheric water fluxes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'backend'))
from app import persistence, plates  # noqa: E402
from app.hydroclimate import climate  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('save', type=Path)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    data = args.save.read_bytes()
    world = persistence.load_world_bytes(data)
    fields = climate.compute_climate(world, *climate.grid_dimensions(world.climate_density))
    height, width = fields.is_ocean.shape
    edges = np.linspace(np.pi / 2, -np.pi / 2, height + 1)
    row_area = (plates.PLANET_RADIUS_KM * 1000) ** 2 * (2 * np.pi / width) * (np.sin(edges[:-1]) - np.sin(edges[1:]))
    area = np.broadcast_to(row_area[:, None], (height, width))
    ocean = fields.is_ocean
    land = ~ocean
    _, _, lake, river, sea = climate._sample_elevation_and_crust(world, fields.world_xyz)
    zero = np.zeros_like(lake)
    vegetation = climate._vegetation_transpiration_source(world, fields.elevation_m, ocean)
    lake_source = climate._land_moisture_source(fields.air_temperature_c, lake, zero, zero, is_sea=sea)
    river_source = climate._land_moisture_source(fields.air_temperature_c, zero, river, zero, is_sea=sea)
    ocean_target = climate.open_water_humidity_g_kg(fields.ocean_temperature_c)

    def area_integral(value, mask):
        return float(np.sum(value[mask] * area[mask]))

    result = {
        'save': args.save.name,
        'save_sha256': hashlib.sha256(data).hexdigest(),
        'grid_shape': [height, width],
        'land_area_km2': float(area[land].sum() / 1e6),
        'ocean_area_km2': float(area[ocean].sum() / 1e6),
        'precipitation_km3_per_year': {
            'land': area_integral(fields.precipitation_mm, land) / 1e12,
            'ocean': area_integral(fields.precipitation_mm, ocean) / 1e12,
        },
        'humidity_proxy_area_integral_g_per_kg_m2': {
            'ocean_target': area_integral(ocean_target, ocean),
            'lake_and_inland_sea_input': area_integral(lake_source, land),
            'river_input': area_integral(river_source, land),
            'vegetation_input': area_integral(vegetation, land),
        },
        'note': 'Ocean target is a reset, not an additive flux; land values are pre-cap per-sweep inputs. None is water mass or an evaporation rate.',
    }
    args.out.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
