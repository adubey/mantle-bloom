#!/usr/bin/env python3
"""Issue #290: attribute every `delaminated_lower_crust_m3` increment and every suture
no-outlet remainder to the call site and front that produced it.

    cd backend
    .venv/bin/python ../bin/debug/attribute_suture_delamination.py --seed 3 --steps 100
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import numpy as np


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    parser.add_argument("--seed", type=int, default=3)
    parser.add_argument("--steps", type=int, default=100)
    parser.add_argument("--years", type=float, default=1_000_000.0)
    args = parser.parse_args()

    sys.path.insert(0, str(args.backend.resolve()))
    import traceback

    from app import continental_ledger, quad_tectonics  # noqa: E402
    from app.world import generate_world, step_world  # noqa: E402

    by_site: collections.Counter = collections.Counter()
    fronts: list[dict] = []
    original_record = continental_ledger.record

    def record(world, account, volume_m3):
        if account == "delaminated_lower_crust_m3" and volume_m3 > 0.0:
            frame = traceback.extract_stack(limit=3)[0]
            by_site[f"{Path(frame.filename).name}:{frame.name}"] += volume_m3
        return original_record(world, account, volume_m3)

    continental_ledger.record = record

    original_place = quad_tectonics._place_suture_crust
    original_accrete = quad_tectonics._accrete_onto_survivors
    current: dict = {}

    def accrete(plate, donors, survivors, world=None, **kwargs):
        current.update(plate_id=plate.plate_id, step=world.steps_taken if world is not None else -1)
        return original_accrete(plate, donors, survivors, world, **kwargs)

    def place(thickness, areas, adjacency, points, front, eligible, volume, cap, strike, root_capacity, shed=None):
        room = float(np.dot(np.maximum(cap - thickness[eligible], 0.0), areas[eligible]))
        changed, stages = original_place(
            thickness, areas, adjacency, points, front, eligible, volume, cap, strike, root_capacity, shed
        )
        if stages["no_outlet_subducted_m3"] > 0.0:
            fronts.append({
                **current,
                "spill_call": root_capacity is None and strike is None,
                "continental_front": root_capacity is not None,
                "front_cells": int(len(front)),
                "eligible_cells": int(eligible.sum()),
                "plate_cells": int(len(thickness)),
                "volume_m3": volume,
                "room_m3": room,
                "no_outlet_m3": stages["no_outlet_subducted_m3"],
            })
        return changed, stages

    quad_tectonics._place_suture_crust = place
    quad_tectonics._accrete_onto_survivors = accrete

    world = generate_world(seed=args.seed)
    for step in range(1, args.steps + 1):
        step_world(world, args.years)
    print(json.dumps({"delaminated_by_site_m3": {k: float(f"{v:.4g}") for k, v in by_site.most_common()}}, indent=2))
    cont = [f for f in fronts if f["continental_front"]]
    print(f"no-outlet fronts: {len(fronts)} ({len(cont)} continental)")
    for key in (True, False):
        sel = [f for f in fronts if f["continental_front"] == key]
        if sel:
            print(f"  continental={key}: no_outlet={sum(f['no_outlet_m3'] for f in sel):.3g} m3")
    per_plate: dict = collections.defaultdict(lambda: {"no_outlet_m3": 0.0, "steps": set(), "cells": []})
    for f in fronts:
        if f["spill_call"] or not f["continental_front"]:
            entry = per_plate[f["plate_id"]]
            entry["no_outlet_m3"] += f["no_outlet_m3"] if f["spill_call"] else 0.0
            entry["steps"].add(f["step"])
            entry["cells"].append(f["plate_cells"])
    print("terminal (post-spill) no-outlet by plate:")
    for pid, e in sorted(per_plate.items(), key=lambda kv: -kv[1]["no_outlet_m3"])[:10]:
        steps = sorted(e["steps"])
        print(f"  plate {pid}: {e['no_outlet_m3']:.3g} m3 over {len(steps)} steps ({steps[0]}-{steps[-1]}), cells {e['cells'][0]} -> {e['cells'][-1]}")
    for f in sorted(fronts, key=lambda f: -f["no_outlet_m3"])[:15]:
        print(json.dumps({k: (float(f"{v:.3g}") if isinstance(v, float) else v) for k, v in f.items()}))


if __name__ == "__main__":
    main()
