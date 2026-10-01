"""Issue #259: per-step churn in the area a quad surface represents, and the sea-level jitter it
drives.

Steps one world and, after every step, sums `SurfaceNodes.area_m2` over all plates (as a
fraction of the sphere) and reads `world.sea_level_m`. Reports the std of the per-step change
in each over the second half of the run (the S2 `stability_sea_level_m` window) and their
correlation.

    cd backend && .venv/bin/python ../bin/debug/measure_quad_area_churn.py --seed 1 --steps 40
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app.sparse_quad_patch import PLANET_RADIUS_M  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402

SPHERE_M2 = 4.0 * np.pi * PLANET_RADIUS_M**2


def represented_fraction(world) -> float:
    return sum(float(np.sum(p.surface_nodes().area_m2)) for p in world.plates if p.node_count()) / SPHERE_M2


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--steps", type=int, default=40)
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--surface", choices=("quad", "lines"), default="quad")
    parser.add_argument("--years", type=float, default=1_000_000.0)
    args = parser.parse_args()

    world = generate_world(seed=args.seed, node_density=args.density, surface=args.surface)
    area, sea = [], []
    for step in range(args.steps):
        step_world(world, years=args.years)
        area.append(represented_fraction(world))
        sea.append(float(world.sea_level_m))
        print(f"step {step + 1:3d}  area/sphere {area[-1]:.5f}  sea level {sea[-1]:8.2f} m", file=sys.stderr)

    half = args.steps // 2
    d_area, d_sea = np.diff(area[half:]), np.diff(sea[half:])
    print(
        json.dumps(
            {
                "seed": args.seed,
                "surface": args.surface,
                "area_range": [round(min(area), 5), round(max(area), 5)],
                "std_d_area": round(float(np.std(d_area)), 6),
                "std_d_sea_level_m": round(float(np.std(d_sea)), 2),
                "corr_d_area_d_sea": round(float(np.corrcoef(d_area, d_sea)[0, 1]), 4),
            }
        )
    )


if __name__ == "__main__":
    main()
