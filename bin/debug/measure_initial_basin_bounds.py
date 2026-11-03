#!/usr/bin/env python3
"""Measure isolated below-sea continental components in fresh synthetic-relief worlds.

Run from the repo root:
    backend/.venv/bin/python bin/debug/measure_initial_basin_bounds.py
    backend/.venv/bin/python bin/debug/measure_initial_basin_bounds.py --seeds 7 11 2024

The defaults and generation settings match the Generate dialog's cratons-first defaults.
This is an initial-world diagnostic only; it does not advance the simulation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad
from app.hydroclimate import hydrology
from app.plates import collect_all_accounting_areas_m2
from app.synthetic_relief import MAX_ISOLATED_HOLE_AREA_KM2, MAX_ISOLATED_HOLE_DEPTH_M
from app.world import generate_world

DEFAULT_SEEDS = (528019447, 7, 11, 2024, 123456789, 987654321, 42, 314159265, 876543210, 135792468, 528019447)
DEPTH_TOLERANCE_M = 1e-6


def measure(seed: int) -> dict[str, object]:
    world = generate_world(
        seed=seed,
        continental_fraction=0.70,
        land_fraction=0.29,
        node_density=4.0,
        climate_density=4.0,
        fluid_density=2.0,
        cratons_first=True,
        synthetic_relief=True,
    )
    points_chunks, elevation_chunks, continental_chunks = [], [], []
    for plate in world.plates:
        points, elevation = plate.all_points_and_elevation()
        points_chunks.append(points)
        elevation_chunks.append(elevation)
        continental_chunks.append(
            effective_is_continental_from_codes(
                plate.collect("crust_type_code"), plate.crust_type == "continental"
            )
        )
    points = np.concatenate(points_chunks)
    elevation = np.concatenate(elevation_chunks)
    continental = np.concatenate(continental_chunks)
    areas_m2 = collect_all_accounting_areas_m2(world.plates, line_spacing_rad(world.node_density))
    neighbours = hydrology._build_neighbor_graph(points, world=world)
    below = elevation <= world.sea_level_m
    ocean = hydrology.connected_ocean_mask(points, elevation, world.sea_level_m, neighbours)
    holes = below & ~ocean & continental
    labels = hydrology._below_sea_components(holes, neighbours)
    components = []
    for label in np.unique(labels[holes]):
        mask = holes & (labels == label)
        components.append(
            {
                "area_km2": float(areas_m2[mask].sum() / 1e6),
                "node_count": int(mask.sum()),
                "floor_m": float(elevation[mask].min()),
            }
        )
    components.sort(key=lambda item: item["area_km2"], reverse=True)
    land = elevation > world.sea_level_m
    return {
        "seed": seed,
        "sea_level_m": float(world.sea_level_m),
        "land_area_fraction": float(areas_m2[land].sum() / areas_m2.sum()),
        "land_node_fraction": float(land.mean()),
        "elevation_sha256": hashlib.sha256(elevation.astype("<f8").tobytes()).hexdigest(),
        "hole_count": len(components),
        "total_hole_area_km2": float(sum(item["area_km2"] for item in components)),
        "largest_hole_area_km2": components[0]["area_km2"] if components else 0.0,
        "deepest_hole_floor_m": min((item["floor_m"] for item in components), default=0.0),
        "all_holes_within_bounds": all(
            item["area_km2"] <= MAX_ISOLATED_HOLE_AREA_KM2
            and item["floor_m"] >= -MAX_ISOLATED_HOLE_DEPTH_M - DEPTH_TOLERANCE_M
            for item in components
        ),
        "components": components,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seeds", type=int, nargs="+", default=DEFAULT_SEEDS)
    args = parser.parse_args()
    results = [measure(seed) for seed in args.seeds]
    first_hash: dict[int, str] = {}
    for result in results:
        seed = int(result["seed"])
        digest = str(result["elevation_sha256"])
        result["deterministic_repeat"] = digest == first_hash[seed] if seed in first_hash else None
        first_hash.setdefault(seed, digest)
    print(json.dumps(results, separators=(",", ":")))


if __name__ == "__main__":
    main()
