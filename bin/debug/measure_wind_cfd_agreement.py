#!/usr/bin/env python3
"""How closely the diagnostic wind model tracks the shallow-water CFD wind (#526).

Steps two copies of one generated world the same way, one under `World.wind_model == "cfd"` and
one under `"diagnostic"`, using the setup of stress_tests/test_climate.py's
test_diagnostic_wind_model_tracks_the_cfd_biome_map, and reports wind-speed rel-RMS, zonal-mean
ocean speeds for both models (the diagnostic one is the 10 m wind; the CFD layer wind runs
weaker, #530), land-biome agreement and precipitation rel-RMS. Both worlds are masked with the
CFD world's land/ocean, as in the stress test, though their terrain drifts apart over the run. It imports the
backend of the worktree it is run from (the current directory), not the one it lives in, so the
same copy measures a baseline worktree too.

Usage, from a worktree root:
    backend/.venv/bin/python <path>/bin/debug/measure_wind_cfd_agreement.py [--seeds 7] [--density 2.0]
    [--steps 10] [--out analysis/issue526/cfd_agreement.json]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.cwd() / "backend"))
from app.hydroclimate import climate  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402

BAND_DEG = 15


def _area_mean(field: np.ndarray, weight: np.ndarray, mask: np.ndarray) -> float:
    return float((field * weight)[mask].sum() / weight[mask].sum())


def _rel_rms(reference: np.ndarray, other: np.ndarray) -> float:
    return float(np.sqrt(np.mean((reference - other) ** 2) / np.mean(reference ** 2)))


def _measure(seed: int, density: float, steps: int) -> dict:
    kwargs = dict(num_plates=10, continental_fraction=0.7, land_fraction=0.29, climate_density=density, fluid_density=density)
    cfd_world = generate_world(seed, **kwargs)
    cfd_world.wind_model = "cfd"  # not the default: World.wind_model is "diagnostic"
    abl_world = generate_world(seed, **kwargs)
    abl_world.wind_model = "diagnostic"
    for _ in range(steps):
        step_world(cfd_world, years=2_000_000)
        step_world(abl_world, years=2_000_000)
    cfd = climate.compute_climate_cached(cfd_world)
    abl = climate.compute_climate_cached(abl_world)

    weight = np.cos(np.radians(cfd.lat_deg))[:, None] * np.ones_like(cfd.precipitation_mm)
    lat = np.broadcast_to(cfd.lat_deg[:, None], weight.shape)
    land, ocean = ~cfd.is_ocean, cfd.is_ocean
    cfd_speed, abl_speed = np.hypot(cfd.wind_u, cfd.wind_v), np.hypot(abl.wind_u, abl.wind_v)
    bands = []
    for lo in range(0, 90, BAND_DEG):
        in_band = (np.abs(lat) >= lo) & (np.abs(lat) < lo + BAND_DEG) & ocean
        if in_band.any():
            bands.append({
                "abs_lat": f"{lo}-{lo + BAND_DEG}",
                "cfd_ocean_m_s": round(_area_mean(cfd_speed, weight, in_band), 2),
                "diagnostic_ocean_m_s": round(_area_mean(abl_speed, weight, in_band), 2),
            })
    return {
        "seed": seed,
        "density": density,
        "steps": steps,
        "wind_speed_rel_rms": round(_rel_rms(cfd_speed, abl_speed), 3),
        "wind_vector_rel_rms": round(float(np.sqrt(
            np.mean((cfd.wind_u - abl.wind_u) ** 2 + (cfd.wind_v - abl.wind_v) ** 2) / np.mean(cfd_speed ** 2)
        )), 3),
        "cfd_ocean_m_s": round(_area_mean(cfd_speed, weight, ocean), 2),
        "diagnostic_ocean_m_s": round(_area_mean(abl_speed, weight, ocean), 2),
        "cfd_land_m_s": round(_area_mean(cfd_speed, weight, land), 2),
        "diagnostic_land_m_s": round(_area_mean(abl_speed, weight, land), 2),
        "land_biome_agreement": round(float(np.mean(cfd.biome_ids[land] == abl.biome_ids[land])), 3),
        "precipitation_rel_rms": round(_rel_rms(cfd.precipitation_mm, abl.precipitation_mm), 3),
        "ocean_bands": bands,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--seeds", default="7", help="comma-separated world seeds")
    parser.add_argument("--density", type=float, default=2.0, help="climate_density and fluid_density")
    parser.add_argument("--steps", type=int, default=10, help="2 Myr steps per world")
    parser.add_argument("--out", help="write the JSON results here")
    args = parser.parse_args()

    results = [_measure(int(seed), args.density, args.steps) for seed in args.seeds.split(",")]
    print(json.dumps(results, indent=1))
    if args.out:
        with open(args.out, "w") as handle:
            json.dump(results, handle, indent=1)
            handle.write("\n")


if __name__ == "__main__":
    main()
