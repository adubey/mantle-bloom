"""Issue #228 Phase 4: measure rift opening and interior subduction on quad vs line worlds.

Steps one world per surface from the same seed and reports, per step window:

- total crustal volume (sum Hc x exact node area) and land fraction;
- the per-phase Hc budget for boundary growth (`boundary_advance` on quad plates,
  `line_growth_shrink`/`row_claim`/... on line plates);
- on quad plates, how many contested (`shrinkable`) cells lie *off* the plate edge -- the
  interior-subduction backlog that edge-only retreat can never reach.

    cd backend && ../.venv/bin/python ../bin/debug/measure_quad_rift_carve.py --seed 7 --steps 60
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import lithosphere, quad_tectonics  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app.lithosphere_plate import boundary_context  # noqa: E402
from app.sparse_quad_patch import PlateWithSparseQuadPatch  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402


def node_areas(world, plate) -> np.ndarray:
    if hasattr(plate, "node_areas_m2"):
        return plate.node_areas_m2()
    return np.full(plate.node_count(), lithosphere.node_area_m2(line_spacing_rad(world.node_density)))


def totals(world) -> dict:
    volume = land_area = area = 0.0
    nodes = 0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        a = node_areas(world, plate)
        volume += float(np.sum(plate.collect("crustal_thickness_m") * a))
        land_area += float(np.sum(a[plate.collect("elevation") > world.sea_level_m]))
        area += float(a.sum())
        nodes += plate.node_count()
    return {"volume_km3": round(volume / 1e9), "land_fraction": round(land_area / area, 4), "nodes": nodes, "plates": len(world.plates)}


def interior_contested(world) -> dict:
    """Contested cells with no exposed side, over every quad plate."""
    edge_total = 0
    interior_total = 0
    for plate in world.plates:
        if not isinstance(plate, PlateWithSparseQuadPatch) or plate.node_count() == 0:
            continue
        others = [p for p in world.plates if p is not plate]
        ctx = boundary_context(world, plate, others, 1.0, lambda c: quad_tectonics.components_of_at_least(plate, c, 3))
        probes = plate._probe_neighbour_indices()
        exposed = np.any(np.all(probes < 0, axis=2), axis=1)
        edge_total += int(np.count_nonzero(ctx.shrinkable & exposed))
        interior_total += int(np.count_nonzero(ctx.shrinkable & ~exposed))
    return {"contested_edge": edge_total, "contested_interior": interior_total}


def budget_summary(world) -> dict:
    out = {}
    for phase, entry in world.phase_budget.items():
        s = entry["scopes"]["all"]
        out[phase] = {
            "d_nodes": s["count_after"] - s["count_before"],
            "d_sum_hc_km": (s["sum_hc_after"] - s["sum_hc_before"]) / 1e3,
        }
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=60)
    parser.add_argument("--years", type=float, default=1_000_000.0)
    parser.add_argument("--every", type=int, default=10)
    parser.add_argument("--surfaces", default="quad,lines")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    report = {}
    for surface in args.surfaces.split(","):
        world = generate_world(seed=args.seed, node_density=args.density, surface=surface)
        world.debug_diagnostics = True
        rows = [{"step": 0, **totals(world)}]
        start = time.perf_counter()
        for step in range(1, args.steps + 1):
            step_world(world, args.years)
            if step % args.every == 0 or step == args.steps:
                row = {"step": step, **totals(world), "elapsed_s": round(time.perf_counter() - start, 1)}
                if surface == "quad":
                    row.update(interior_contested(world))
                rows.append(row)
                print(surface, json.dumps(row), flush=True)
        report[surface] = {"rows": rows, "phase_budget": budget_summary(world)}
        for phase, v in sorted(report[surface]["phase_budget"].items()):
            print(f"  {surface:5s} {phase:32s} d_nodes={v['d_nodes']:>9.0f}  d_sum_hc_km={v['d_sum_hc_km']:>12.1f}")
    if args.out:
        args.out.write_text(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
