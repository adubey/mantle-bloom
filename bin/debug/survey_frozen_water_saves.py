#!/usr/bin/env python3
"""RFC 0015 seed sets: summarise candidate saves for the diagnosis and held-out sets (#548).

Reads each `.mbworld` save without stepping it and prints one JSON line per save: seed, age,
node density, axial tilt, solar multiplier, ice-age period and wind model, the SHA-256, and
what the gates need a save to have: grounded ice volume and visible ice area (plate nodes,
exact areas), and the share of the ocean poleward of 60 degrees in each hemisphere (climate
grid). No gate is measured here, so surveying a save doesn't spend it for the held-out set.
Saves are loaded one at a time (peak RSS about 2 GB).

Usage:

    cd backend
    .venv/bin/python ../bin/debug/survey_frozen_water_saves.py <save>.mbworld [...] --out <scratch>/survey.jsonl
"""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import numpy as np
from replay_lib import load_replay

POLAR_LAT_DEG = 60.0


def survey(world) -> dict:
    from app import plates as plates_mod
    from app.elevation_lines import line_spacing_rad
    from app.hydroclimate import climate, hydrology

    area = plates_mod.collect_all_accounting_areas_m2(world.plates, line_spacing_rad(world.node_density))
    depth = np.clip(plates_mod.collect_all_glacier_depth(world.plates), 0.0, None)
    lat_deg, _lon, world_xyz = climate._build_grid(climate.GRID_HEIGHT, climate.GRID_WIDTH)
    _elev, is_ocean, *_rest = climate._sample_elevation_and_crust(world, world_xyz)
    w = np.broadcast_to(np.cos(np.radians(lat_deg)).reshape(-1, 1), is_ocean.shape)
    lat = np.broadcast_to(lat_deg.reshape(-1, 1), is_ocean.shape)
    ocean_area = float(w[is_ocean].sum())
    return {
        "seed": int(world.seed),
        "elapsed_myr": round(world.elapsed_years / 1e6, 3),
        "node_density": world.node_density,
        "axial_tilt_deg": float(world.axial_tilt_deg),
        "solar_multiplier": float(world.solar_multiplier),
        "ice_age_period_years": float(getattr(world, "ice_age_period_years", 0.0) or 0.0),
        "wind_model": getattr(world, "wind_model", None),
        "ice_volume_km3": float(np.dot(depth, area)) / 1e9,
        "visible_ice_area_km2": float(area[depth >= hydrology.GLACIER_VISIBLE_DEPTH_M].sum()) / 1e6,
        "ocean_share_of_globe": ocean_area / float(w.sum()),
        "ocean_poleward_60N_share_of_ocean": float(w[is_ocean & (lat >= POLAR_LAT_DEG)].sum()) / ocean_area,
        "ocean_poleward_60S_share_of_ocean": float(w[is_ocean & (lat <= -POLAR_LAT_DEG)].sum()) / ocean_area,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("saves", type=Path, nargs="+")
    ap.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w") as out:
        for save in args.saves:
            args.save, args.step_years = save, None
            world, header = load_replay(args, issue=548)
            row = {"save": save.name, "save_sha256": header["save_sha256"], "commit": header["commit"], **survey(world)}
            out.write(json.dumps(row) + "\n")
            out.flush()
            print(json.dumps(row), flush=True)
            del world
            gc.collect()


if __name__ == "__main__":
    main()
