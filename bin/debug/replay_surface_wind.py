#!/usr/bin/env python3
"""Replay a saved world and track the wind's downstream outcomes (#526): land weathering volume,
10 m wind speeds, mean ocean-current speed, precipitation, land area and the land-biome map.

Weathering volume is the demand `erosion.apply_erosion` forms from WEATHERING_COEFFICIENT
(`ErosionResult.weathering`, m per node, before the per-node erosion cap) integrated over node
area on land nodes, captured by wrapping apply_erosion so both arms of an A/B measure it the same
way. It imports the backend of the worktree it is run from (the current directory), so one copy
measures a baseline worktree too. Climate is read from `World.climate_cache`, the snapshot
erosion computed that step. `--baseline` takes a previous run's `--out` and adds land-biome
agreement against its maps at each checkpoint (both runs must share save and checkpoints).

Usage, from a worktree root:
    backend/.venv/bin/python <path>/bin/debug/replay_surface_wind.py <local-save> --myr 10 \
        --checkpoint-myr 2.5 --out <scratch-output>/branch [--baseline <scratch-output>/main]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path.cwd() / "backend"))
from app import persistence, world as world_mod  # noqa: E402
from app.hydroclimate import erosion  # noqa: E402


def _install_weathering_probe(totals: dict) -> None:
    original = erosion.apply_erosion

    def probed(world, years, node_cloud=None):
        result = original(world, years, node_cloud=node_cloud)
        if result is not None and node_cloud is not None:
            area = erosion._gather_areas(world, node_cloud[1])
            land = result.elevation > world.sea_level_m
            totals["weathering_m3"] += float(np.sum(result.weathering[land] * area[land]))
        return result

    # world.step_world calls it as erosion.apply_erosion, so patching the module is enough.
    erosion.apply_erosion = probed


def _peak_rss_mb() -> float:
    # ru_maxrss is bytes on macOS and kilobytes on Linux.
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / 2**20 if sys.platform == "darwin" else peak / 2**10


def _snapshot(world, weathering_m3: float, myr: float) -> tuple[dict, np.ndarray]:
    fields = world.climate_cache
    weight = np.cos(np.radians(fields.lat_deg))[:, None] * np.ones_like(fields.precipitation_mm)
    land, ocean = ~fields.is_ocean, fields.is_ocean

    def mean(field, mask):
        return round(float((field * weight)[mask].sum() / weight[mask].sum()), 4)

    wind = np.hypot(fields.wind_u, fields.wind_v)
    current = np.hypot(fields.current_u, fields.current_v)
    return {
        "myr": round(myr, 3),
        "weathering_m3_since_last": weathering_m3,
        "wind_land_m_s": mean(wind, land),
        "wind_ocean_m_s": mean(wind, ocean),
        "current_ocean_mean": mean(current, ocean),
        "precipitation_mm": mean(fields.precipitation_mm, land | ocean),
        "precipitation_land_mm": mean(fields.precipitation_mm, land),
        "land_area_fraction": mean(land.astype(float), land | ocean),
        "land_mean_elevation_m": mean(fields.elevation_m, land),
    }, np.where(land, fields.biome_ids, -1)  # -1 marks ocean, so agreement is land-only


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("save", type=Path)
    parser.add_argument("--myr", type=float, required=True, help="how far to replay")
    parser.add_argument("--step-years", type=float, default=100_000.0, help="UI default step")
    parser.add_argument("--checkpoint-myr", type=float, default=2.5)
    parser.add_argument("--out", type=Path, required=True, help="writes <out>.json and <out>_biomes.npz")
    parser.add_argument("--baseline", type=Path, help="an earlier run's --out, for biome agreement")
    args = parser.parse_args()

    data = args.save.read_bytes()
    world = persistence.load_world_bytes(data)
    totals = {"weathering_m3": 0.0}
    _install_weathering_probe(totals)

    steps = int(round(args.myr * 1e6 / args.step_years))
    per_checkpoint = max(1, int(round(args.checkpoint_myr * 1e6 / args.step_years)))
    baseline_biomes = np.load(f"{args.baseline}_biomes.npz") if args.baseline else None
    checkpoints, biomes = [], {}
    started = time.time()
    for step in range(1, steps + 1):
        world_mod.step_world(world, args.step_years)
        if step % per_checkpoint == 0 or step == steps:
            row, biome_ids = _snapshot(world, totals["weathering_m3"], step * args.step_years / 1e6)
            totals["weathering_m3"] = 0.0
            key = f"step{step}"
            biomes[key] = biome_ids
            checkpoints.append(row)
            print(json.dumps(row), flush=True)

    if baseline_biomes is not None:
        for row, key in zip(checkpoints, biomes):
            if key in baseline_biomes.files:
                ours, theirs = biomes[key], baseline_biomes[key]
                land = theirs != -1
                row["land_biome_agreement_vs_baseline"] = round(float(np.mean(ours[land] == theirs[land])), 4)

    commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    summary = {
        "save_sha256": hashlib.sha256(data).hexdigest(),
        "commit": commit,
        "wind_model": world.wind_model,
        "myr": args.myr,
        "step_years": args.step_years,
        "runtime_s": round(time.time() - started, 1),
        "peak_rss_mb": round(_peak_rss_mb(), 1),
        "checkpoints": checkpoints,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    Path(f"{args.out}.json").write_text(json.dumps(summary, indent=1) + "\n")
    np.savez_compressed(f"{args.out}_biomes.npz", **biomes)


if __name__ == "__main__":
    main()
