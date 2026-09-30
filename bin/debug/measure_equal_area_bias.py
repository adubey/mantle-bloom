#!/usr/bin/env python3
"""Issue #249: how far the "every node has the nominal area" assumption is off on line vs quad
worlds.

Several whole-world calculations -- eustasy's water budget, stats' land area and continental
volume, torque's asthenosphere drag, magma transport's deposit thickness -- weight every node
by `lithosphere.node_area_m2(spacing)`. This reports, per surface and age, each node's actual
area (`SurfaceNodes.area_m2`) against that nominal value, and the resulting error in ocean
volume and land area.

    cd backend
    .venv/bin/python ../bin/debug/measure_equal_area_bias.py --seed 7 --steps 0,10
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import eustasy, lithosphere  # noqa: E402
from app.elevation_lines import PLANET_RADIUS_KM, line_spacing_rad  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402


def report(world, surface: str) -> str:
    live = [p for p in world.plates if p.node_count()]
    area = np.concatenate([p.surface_nodes().area_m2 for p in live])
    elevation = eustasy.all_elevations(world)
    nominal = lithosphere.node_area_m2(line_spacing_rad(world.node_density))
    depth = np.clip(world.sea_level_m - elevation, 0.0, None)
    land = elevation > world.sea_level_m
    sphere = 4.0 * np.pi * (PLANET_RADIUS_KM * 1000.0) ** 2
    ratio = np.percentile(area / nominal, [5, 50, 95])
    return (
        f"| {surface} | {world.elapsed_years / 1e6:g} | {len(area)} | {ratio[0]:.3f} / {ratio[1]:.3f} / {ratio[2]:.3f} | "
        f"{area.sum() / sphere:.4f} | {np.sum(depth * area) / (np.sum(depth) * nominal):.4f} | "
        f"{area[land].sum() / (land.sum() * nominal):.4f} |"
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--steps", default="0,10", help="comma-separated step counts (1 Myr each) to report at")
    args = parser.parse_args()
    report_at = sorted(int(s) for s in args.steps.split(","))
    print("| surface | Myr | nodes | node area / nominal, p05 / p50 / p95 | Σ area / sphere | ocean volume, real / equal-area | land area, real / equal-area |")
    print("|---|---:|---:|---|---:|---:|---:|")
    for surface in ("lines", "quad"):
        world = generate_world(seed=args.seed, node_density=args.density, surface=surface)
        for step in range(report_at[-1] + 1):
            if step in report_at:
                print(report(world, surface), flush=True)
            if step < report_at[-1]:
                step_world(world, 1e6)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
