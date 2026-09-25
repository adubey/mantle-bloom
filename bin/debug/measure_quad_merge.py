"""Issue #228 Phase 4: measure cross-plate merge on quad vs line worlds.

Steps one world per variant from the same seed:

- `quad`: quad plates with quad merge (quad_merge.py);
- `quad-nomerge`: quad plates with merge disabled, as before quad merge landed;
- `lines`: the line engine.

Every `merge_plates` call is wrapped to book the fusing pair's crustal volume before and the
surviving plate's after, both as sum Hc x exact node area (quad) or x nominal node area
(lines, whose merge resamples onto a constant-area lattice), plus the call's runtime. Every
`--every` steps it also reports total crustal volume, land fraction, node and plate counts.

    cd backend && .venv/bin/python ../bin/debug/measure_quad_merge.py --seed 7 --steps 150
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import lithosphere, merge_split  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app.world import generate_world, step_world  # noqa: E402


def node_areas(world, plate) -> np.ndarray:
    if hasattr(plate, "node_areas_m2"):
        return plate.node_areas_m2()
    return np.full(plate.node_count(), lithosphere.node_area_m2(line_spacing_rad(world.node_density)))


def plate_volume(world, plate) -> float:
    return float(np.sum(plate.collect("crustal_thickness_m") * node_areas(world, plate)))


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
    return {"volume_1e9_km3": round(volume / 1e18, 4), "land_fraction": round(land_area / area, 4), "nodes": nodes, "plates": len(world.plates)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--steps", type=int, default=150)
    parser.add_argument("--years", type=float, default=1_000_000.0)
    parser.add_argument("--every", type=int, default=25)
    parser.add_argument("--variants", default="quad,quad-nomerge,lines")
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    original_merge = merge_split.merge_plates
    original_supports = merge_split._supports_merge
    report = {}
    for variant in args.variants.split(","):
        surface = "lines" if variant == "lines" else "quad"
        world = generate_world(seed=args.seed, node_density=args.density, surface=surface)
        merges: list[dict] = []

        def measured_merge(world, id_keep, id_absorb, merges=merges):
            keep = next(p for p in world.plates if p.plate_id == id_keep)
            absorb = next(p for p in world.plates if p.plate_id == id_absorb)
            before = plate_volume(world, keep) + plate_volume(world, absorb)
            cells_before = keep.node_count() + absorb.node_count()
            start = time.perf_counter()
            original_merge(world, id_keep, id_absorb)
            seconds = time.perf_counter() - start
            after = plate_volume(world, keep)
            merges.append({
                "step": world.steps_taken, "keep": id_keep, "absorb": id_absorb,
                "nodes_before": cells_before, "nodes_after": keep.node_count(),
                "volume_rel_change": (after - before) / before, "seconds": round(seconds, 3),
            })
            print(variant, "merge", json.dumps(merges[-1]), flush=True)

        merge_split.merge_plates = measured_merge
        merge_split._supports_merge = (lambda *_: False) if variant == "quad-nomerge" else original_supports
        try:
            rows = [{"step": 0, **totals(world)}]
            start = time.perf_counter()
            for step in range(1, args.steps + 1):
                step_world(world, args.years)
                if step % args.every == 0 or step == args.steps:
                    row = {"step": step, **totals(world), "elapsed_s": round(time.perf_counter() - start, 1)}
                    rows.append(row)
                    print(variant, json.dumps(row), flush=True)
        finally:
            merge_split.merge_plates = original_merge
            merge_split._supports_merge = original_supports
        report[variant] = {"rows": rows, "merges": merges}
    if args.out:
        args.out.write_text(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
