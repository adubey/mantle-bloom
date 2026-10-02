#!/usr/bin/env python3
"""Issue #249: where the continental crust quad boundary retreat removes goes.

`attribute_continental_budget.py` shows quad `_retreat` as the dominant continental-crust sink
on quad worlds. This splits every retreat call's removed continental volume (effective
crust type continental, sum(Hc * area)) by what happens to it:

- `accreted_kept`: suture donors (`ctx.accrete`) plus continental terranes riding oceanic
  plates whose volume `_accrete_onto_survivors` put back onto surviving cells;
- `accretion_cap_loss`: donor volume clipped by `SUTURE_ACCRETION_MAX_HC_M`;
- `subducted_continental_plate`: removed non-donor cells of continental plates (the
  arc-magmatism-budgeted override retreat);
- `subducted_oceanic_plate`: removed continental cells of oceanic plates (terranes riding an
  overridden oceanic plate, and interior carve).

    cd backend
    .venv/bin/python ../bin/debug/probe_quad_retreat_continental.py --seed 1 --steps 240
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import elevation_lines, quad_tectonics, world as world_mod  # noqa: E402


def continental_volume(plate, mask=None) -> float:
    hc = plate.collect("crustal_thickness_m")
    areas = plate.node_areas_m2()
    continental = elevation_lines.effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
    if mask is not None:
        continental = continental & mask
    return float(np.sum((hc * areas)[continental]) / 1e9)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--steps", type=int, default=240)
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--every", type=int, default=20)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    totals = defaultdict(float)
    original_retreat = quad_tectonics._retreat
    original_accrete = quad_tectonics._accrete_onto_survivors
    current = {}

    def accrete(plate, donors, survivors):
        before = continental_volume(plate, survivors)
        donated = continental_volume(plate, donors)
        original_accrete(plate, donors, survivors)
        gained = continental_volume(plate, survivors) - before
        current["accreted_kept"] = gained
        current["accretion_cap_loss"] = donated - gained

    def retreat(plate, world, ctx, max_distance, max_cells):
        current.clear()
        n = plate.node_count()
        snapshot = (plate.collect("crustal_thickness_m"), plate.node_areas_m2(), plate.collect("crust_type_code"))
        survivors = original_retreat(plate, world, ctx, max_distance, max_cells)
        if len(survivors) != n or np.all(survivors):
            return survivors
        hc, areas, codes = snapshot
        continental = elevation_lines.effective_is_continental_from_codes(codes, plate.crust_type == "continental")
        removed = ~survivors
        volume = (hc * areas / 1e9)[removed & continental]
        terrane = continental if plate.crust_type == "oceanic" else np.zeros_like(continental)
        donors = (removed & (ctx.accrete | terrane))[removed & continental]
        totals["removed_continental"] += float(volume.sum())
        if plate.crust_type == "continental":
            totals["subducted_continental_plate"] += float(volume[~donors].sum())
        else:
            totals["subducted_oceanic_plate"] += float(volume[~donors].sum())
            # Oceanic-plate donors (rare) count as accreted below like any other.
        totals["accreted_kept"] += current.get("accreted_kept", 0.0)
        totals["accretion_cap_loss"] += current.get("accretion_cap_loss", 0.0)
        totals["calls"] += 1
        return survivors

    quad_tectonics._retreat = retreat
    quad_tectonics._accrete_onto_survivors = accrete
    rows = []
    try:
        world = world_mod.generate_world(seed=args.seed, node_density=args.density, surface="quad")
        start = sum(continental_volume(p) for p in world.plates if p.node_count())
        for step in range(1, args.steps + 1):
            world_mod.step_world(world, 1e6)
            if step % args.every == 0:
                row = {"step": step, "continental_now": sum(continental_volume(p) for p in world.plates if p.node_count()) / start}
                row.update({k: (v / start if k != "calls" else v) for k, v in totals.items()})
                rows.append(row)
                print(json.dumps({k: round(v, 4) for k, v in row.items()}), flush=True)
    finally:
        quad_tectonics._retreat = original_retreat
        quad_tectonics._accrete_onto_survivors = original_accrete
    if args.out:
        args.out.write_text(json.dumps({"seed": args.seed, "density": args.density, "start_km3": start, "rows": rows}, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
