#!/usr/bin/env python3
"""Issue #275: measure how closely one `erosion.apply_erosion` call conserves material.

Generates a world per seed, runs erosion alone a few times and, around every call, measures

- total crustal (Hc) volume, sum(Hc * area) over every node with a column;
- the continental-material tracer's live volume (`continental_ledger.surface_volume_m3`);
- the ledger's declared sinks, so a tracer drop can be split into "declared" and "silent";
- how much tracer volume sits on oceanic host columns;
- wall time.

Areas are each plate's own `accounting_areas_m2` (exact cell areas on quad worlds).

    cd backend
    .venv/bin/python ../bin/debug/measure_erosion_budget.py --seeds 3 4 --steps 5 --surface quad
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import continental_ledger, erosion  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app.world import generate_world  # noqa: E402

SINKS = ("delaminated_lower_crust_m3", "deeply_subducted_m3", "numerical_unplaced_m3")


def hc_volume_m3(world) -> float:
    spacing = line_spacing_rad(world.node_density)
    return sum(float(np.sum(p.collect("crustal_thickness_m") * p.accounting_areas_m2(spacing))) for p in world.plates)


def snapshot(world) -> dict[str, float]:
    inv = continental_ledger.inventories(world)
    return {
        "hc_m3": hc_volume_m3(world),
        "tracer_m3": inv["surface_continental_derived_m3"],
        "on_ocean_m3": inv["continental_sediment_on_oceanic_hosts_m3"],
        "sinks_m3": sum(inv[k] for k in SINKS),
    }


def run(seed: int, steps: int, years: float, surface: str, num_plates: int) -> list[dict[str, float]]:
    world = generate_world(seed=seed, num_plates=num_plates, surface=surface)
    rows = []
    for step in range(steps):
        before = snapshot(world)
        start = time.perf_counter()
        result = erosion.apply_erosion(world, years)
        elapsed = time.perf_counter() - start
        after = snapshot(world)
        tracer_change = after["tracer_m3"] - before["tracer_m3"]
        declared = after["sinks_m3"] - before["sinks_m3"]
        budget = getattr(result, "budget", None) or {}
        rows.append(
            {
                "seed": seed,
                "step": step,
                "seconds": round(elapsed, 3),
                "hc_change_m3": after["hc_m3"] - before["hc_m3"],
                "tracer_change_m3": tracer_change,
                "declared_sinks_m3": declared,
                "silent_tracer_change_m3": tracer_change + declared,
                "silent_relative": (tracer_change + declared) / max(before["tracer_m3"], 1.0),
                "on_ocean_m3": after["on_ocean_m3"],
                **{f"budget.{k}": v for k, v in budget.items()},
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", type=int, nargs="+", default=[3])
    parser.add_argument("--steps", type=int, default=3)
    parser.add_argument("--years", type=float, default=5_000_000.0)
    parser.add_argument("--surface", default="quad")
    parser.add_argument("--num-plates", type=int, default=8)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()

    rows = [row for seed in args.seeds for row in run(seed, args.steps, args.years, args.surface, args.num_plates)]
    for row in rows:
        print(json.dumps({k: (f"{v:.4g}" if isinstance(v, float) else v) for k, v in row.items()}))
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
