#!/usr/bin/env python3
"""Issue #758, post hoc: why nearby land got drier when D's basin became ocean.

Replays a save a few steps and, after each, reads the climate cache over two regions fixed at
the start: the save's largest enclosed below-sea-level body (climate cells within
`--exclude-km` of its cells) and the land ring around it (land cells within `--near-km` but
beyond `--exclude-km`, as `replay_enclosed_ocean.py` samples). Over each it reports the share
of ocean and visible-lake cells, and cos-latitude-weighted mean elevation, air temperature,
ocean temperature (ocean cells), humidity and precipitation.

Run from `backend/` so the worktree's environment is active:

    .venv/bin/python ../bin/debug/probe_enclosed_ocean_climate.py <save> --myr 0.3 \\
        --step-years 100000 [--backend <baseline>/backend] --out <scratch>/probe.json
"""

from __future__ import annotations

import json

import numpy as np
from replay_lib import load_replay, parse_replay_args
from scipy.spatial import cKDTree


def _arguments(ap):
    ap.add_argument("--near-km", type=float, default=1_000.0)
    ap.add_argument("--exclude-km", type=float, default=100.0)


args = parse_replay_args(__doc__, _arguments)
world, result = load_replay(args, issue=758)

from app import plates  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.elevation_lines import PLANET_RADIUS_KM  # noqa: E402
from app.hydroclimate import hydrology  # noqa: E402

world.ice_age_period_years = 0.0


def chord(km):
    return 2.0 * np.sin(km / PLANET_RADIUS_KM / 2.0)


hydro = world.hydrology_cache
elevation = plates.collect_all_elevation(hydro.plates_in_order)
below = elevation <= world.sea_level_m
labels = hydrology._below_sea_components(below, hydro.neighbor_idx)
counts = np.bincount(labels[below])
largest = int(counts.argmax())
counts[largest] = 0
body_points = hydro.points[below & (labels == int(counts.argmax()))]
fields = world.climate_cache
dist, _ = cKDTree(body_points).query(fields.world_xyz.reshape(-1, 3), distance_upper_bound=chord(args.near_km), workers=-1)
dist = dist.reshape(fields.precipitation_mm.shape)
regions = {"basin": dist <= chord(args.exclude_km), "ring": (dist > chord(args.exclude_km)) & np.isfinite(dist)}


def mean(values, mask, weight):
    w = weight * mask
    return float((values * w).sum() / w.sum()) if w.sum() > 0 else None


def read(step):
    f = world.climate_cache
    weight = np.cos(np.radians(f.lat_deg))[:, None] * np.ones_like(f.precipitation_mm)
    lake = ~f.is_ocean & (f.lake_depth_m >= hydrology.LAKE_MIN_VISIBLE_DEPTH_M)
    row = {"step": step, "elapsed_myr": world.elapsed_years / 1e6, "sea_level_m": world.sea_level_m}
    for name, region in regions.items():
        land = region & ~f.is_ocean & ~lake
        row[name] = {
            "cells": int(region.sum()),
            "ocean_share": mean(f.is_ocean.astype(float), region, weight),
            "lake_share": mean(lake.astype(float), region, weight),
            "elevation_m": mean(f.elevation_m, region, weight),
            "air_c": mean(f.air_temperature_c, region, weight),
            "ocean_c": mean(f.ocean_temperature_c, region & f.is_ocean, weight),
            "humidity_g_kg": mean(f.humidity_g_kg, region, weight),
            "precipitation_mm": mean(f.precipitation_mm, region, weight),
            "land_humidity_g_kg": mean(f.humidity_g_kg, land, weight),
            "land_precipitation_mm": mean(f.precipitation_mm, land, weight),
            "land_air_c": mean(f.air_temperature_c, land, weight),
        }
    return row


steps = int(round(args.myr * 1e6 / args.step_years))
result["steps"] = [read(0)]
for i in range(1, steps + 1):
    world_mod.step_world(world, args.step_years)
    result["steps"].append(read(i))
    print(json.dumps(result["steps"][-1]), flush=True)
args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(result, indent=1))
print(f"wrote {args.out}")
