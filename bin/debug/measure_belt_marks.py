#!/usr/bin/env python3
"""Replay A/B for RFC 0009's belt marks (issue #461): belt area, craton formation and ledgers
on default worlds, and step time at density 4.

`run` generates a default world per seed (today's plates-first layout) and steps it to
--years, writing a JSONL row every --every years: belt area (0 on a tree without belt marks),
craton area and volume, cumulative formed volume, genuine continental area, the craton and
continental-material ledger balances, and the formation audit -- craton volume and cells that
formed this run on columns marked as belt, which must stay 0. `bench` times --steps steps at
--density after --warmup steps, --repeats times. Each row records which backend ran
(`model_file`): the script puts its own tree's backend first on sys.path, so a baseline run
must use a copy of this file inside the baseline worktree.

Usage:
    backend/.venv/bin/python bin/debug/measure_belt_marks.py run --seeds 7,11,2024 \
        --years 1e9 --every 100e6 --out /tmp/belts-branch.jsonl
    backend/.venv/bin/python bin/debug/measure_belt_marks.py bench --seed 7 --density 4 \
        --out /tmp/bench-branch.jsonl
    backend/.venv/bin/python bin/debug/measure_belt_marks.py summarize /tmp/belts-main.jsonl /tmp/belts-branch.jsonl
"""

from __future__ import annotations

import argparse
import json
import resource
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

import numpy as np  # noqa: E402

from app import continental_ledger, cratons  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad  # noqa: E402

HAS_BELTS = hasattr(cratons, "belt_mask")


class FormationAudit:
    """Wraps `cratons.update` to count craton formed on belt-marked columns."""

    def __init__(self) -> None:
        self.cells = 0
        self.volume_m3 = 0.0
        self._update = cratons.update
        cratons.update = self._audited

    def _audited(self, world, years: float) -> None:
        before = {}
        if HAS_BELTS:
            for plate in world.plates:
                if plate.node_count():
                    before[plate.plate_id] = (plate.collect("craton_crust_m") <= 0.0) & cratons.belt_mask(plate, world.elapsed_years)
        self._update(world, years)
        for plate in world.plates:
            mask = before.get(plate.plate_id)
            if mask is None or len(mask) != plate.node_count():
                continue
            formed = mask & (plate.collect("craton_crust_m") > 0.0)
            self.cells += int(formed.sum())
            areas = plate.accounting_areas_m2(line_spacing_rad(world.node_density))
            self.volume_m3 += float(plate.collect("craton_crust_m")[formed] @ areas[formed])


def _snapshot(world, audit: FormationAudit) -> dict:
    spacing = line_spacing_rad(world.node_density)
    genuine = 0.0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        continental = effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        hc = plate.collect("crustal_thickness_m")
        genuine += float(plate.accounting_areas_m2(spacing)[continental & (hc >= cratons.CRATON_HOST_MIN_HC_M)].sum())
    craton = cratons.diagnostics(world)
    return {
        "myr": world.elapsed_years / 1e6,
        "belt_area_km2": craton.get("craton_belt_area_m2", 0.0) / 1e6,
        "core_count": craton.get("craton_core_count"),
        "craton_area_km2": craton["craton_area_m2"] / 1e6,
        "craton_volume_km3": craton["craton_volume_m3"] / 1e9,
        "formed_km3": world.craton_ledger.get("formed_m3", 0.0) / 1e9,
        "genuine_continental_area_km2": genuine / 1e6,
        "craton_balance_error_km3": craton["craton_balance_error_m3"] / 1e9,
        "continental_balance_error_km3": continental_ledger.balance_error_m3(world) / 1e9,
        "formed_on_belt_cells": audit.cells,
        "formed_on_belt_km3": audit.volume_m3 / 1e9,
        "plate_count": len(world.plates),
    }


def run(args) -> None:
    audit = FormationAudit()
    with open(args.out, "a") as out:
        for seed in (int(s) for s in args.seeds.split(",")):
            audit.cells, audit.volume_m3 = 0, 0.0
            world = world_mod.generate_world(seed)
            started = time.time()
            rows = [_snapshot(world, audit)]
            next_checkpoint = args.every
            while world.elapsed_years < args.years - 1.0:
                world_mod.step_world(world, min(args.step, args.years - world.elapsed_years))
                if world.elapsed_years >= next_checkpoint - 1.0:
                    rows.append(_snapshot(world, audit))
                    next_checkpoint += args.every
                    print(json.dumps({"seed": seed, **rows[-1]}), flush=True)
            for row in rows:
                row.update(seed=seed, model_file=cratons.__file__, wall_s=time.time() - started)
                out.write(json.dumps(row) + "\n")
            out.flush()


def bench(args) -> None:
    with open(args.out, "a") as out:
        for repeat in range(args.repeats):
            world = world_mod.generate_world(args.seed, node_density=args.density)
            world.simulate_climate_biomes = False
            for _ in range(args.warmup):
                world_mod.step_world(world, args.step)
            times = []
            for _ in range(args.steps):
                started = time.perf_counter()
                world_mod.step_world(world, args.step)
                times.append(time.perf_counter() - started)
            row = {
                "seed": args.seed, "density": args.density, "repeat": repeat, "step_s": times,
                "median_step_s": float(np.median(times)), "model_file": cratons.__file__,
                "max_rss_mb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20,
            }
            out.write(json.dumps(row) + "\n")
            out.flush()
            print(json.dumps({k: row[k] for k in ("repeat", "median_step_s")}), flush=True)


def summarize(args) -> None:
    runs = {}
    for path in args.paths:
        rows = [json.loads(line) for line in Path(path).read_text().splitlines() if line.strip()]
        runs[Path(path).stem] = rows
    for name, rows in runs.items():
        print(f"== {name} ({rows[0].get('model_file')})")
        if "median_step_s" in rows[0]:
            medians = [r["median_step_s"] for r in rows]
            print(f"  median step {np.median(medians):.3f} s over repeats {', '.join(f'{m:.3f}' for m in medians)}")
            continue
        for row in rows:
            if round(row["myr"]) % 100 == 0:
                print(
                    f"  seed {row['seed']:5d} {row['myr']:6.0f} Myr belt {row['belt_area_km2'] / 1e6:6.2f} M km2"
                    f"  craton {row['craton_area_km2'] / 1e6:6.2f} M km2  formed {row['formed_km3'] / 1e6:7.2f} M km3"
                    f"  continental {row['genuine_continental_area_km2'] / 1e6:6.1f} M km2"
                    f"  on-belt {row['formed_on_belt_cells']}  balance {row['craton_balance_error_km3']:.2g}/"
                    f"{row['continental_balance_error_km3']:.2g} km3"
                )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("run")
    p.add_argument("--seeds", default="7,11,2024")
    p.add_argument("--years", type=float, default=1e9)
    p.add_argument("--every", type=float, default=100e6)
    p.add_argument("--step", type=float, default=5e6)
    p.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("bench")
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--density", type=float, default=4.0)
    p.add_argument("--step", type=float, default=2e6)
    p.add_argument("--warmup", type=int, default=2)
    p.add_argument("--steps", type=int, default=8)
    p.add_argument("--repeats", type=int, default=3)
    p.add_argument("--out", type=Path, required=True)
    p = sub.add_parser("summarize")
    p.add_argument("paths", nargs="+")
    args = parser.parse_args()
    {"run": run, "bench": bench, "summarize": summarize}[args.command](args)


if __name__ == "__main__":
    main()
