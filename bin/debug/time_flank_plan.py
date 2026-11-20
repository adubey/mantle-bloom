#!/usr/bin/env python3
"""Time RFC 0018's flank-recipient search (`rift_strength.flank_plan`) as a rift's recipients recede (#987).

Generates a plates-only world with `split_method = "rift_lines"` and continental braking off, then:

- `--steps N` steps it N times (`--step-myr` each), timing every `rift_damage.damage_pass` and `flank_plan` call and
  recording each step's flank `hops_max`, held donors and broken cells.
- For the `--plates` largest plates with a rift line, it builds the line's band as `rift_damage._line_problem` does,
  marks every unbroken plate cell within K hops of the band broken for each K in `--extra-hops` (the issue's
  "widening broken corridor"), and times one `flank_plan` call with the damage pass's crust-type groups.

Every plan and the stepped world's rift and column state get a SHA-256 digest, so running the same command in an
origin/main worktree gives the before/after timings and shows whether the plans and states are bit-identical.

Usage (density 4 peaks near 3 GB; one process at a time, see the compute-budget skill):
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/time_flank_plan.py --seed 106533960 --node-density 4 \\
      --steps 10 --params '{"k_cap": 0.14}' --extra-hops 0 5 10 20 40 80 --out <report.json>
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

from app import merge_split, rift_damage, rift_lines, rift_strength
from app import world as world_module

STATE_FIELDS = ("rift_strain", "rift_warm_strain", "rift_broken", "crustal_thickness_m", "mantle_lithosphere_thickness_m")


def plan_digest(plan) -> str:
    reach = plan.reach.tocsr()
    reach.sort_indices()
    h = hashlib.sha256()
    for array in (plan.band, plan.held, plan.hops, reach.indptr, reach.indices):
        h.update(np.ascontiguousarray(array).astype(np.int64).tobytes())
    return h.hexdigest()[:16]


def state_digest(world) -> str:
    h = hashlib.sha256()
    for plate in sorted(world.plates, key=lambda p: p.plate_id):
        h.update(np.int64(plate.plate_id).tobytes())
        for name in STATE_FIELDS:
            h.update(np.ascontiguousarray(plate.collect(name), dtype=np.float64).tobytes())
    return h.hexdigest()[:16]


def make_world(args):
    world = world_module.generate_world(args.seed, node_density=args.node_density)
    world.simulate_climate_biomes = False
    world.contact_drag_multiplier = 0.0
    world.craton_drag_multiplier = 0.0
    world.split_method = "rift_lines"
    world.rift_line_params = dict(args.params)
    merge_split.check_split_settings(world)
    return world


def step_with_timers(world, args) -> list[dict]:
    """Step the world, timing the damage pass and every flank plan inside it."""
    original_pass, original_plan = rift_damage.damage_pass, rift_strength.flank_plan
    current: dict = {}

    def timed_pass(world, years):
        started = time.perf_counter()
        try:
            return original_pass(world, years)
        finally:
            current["damage_s"] = current.get("damage_s", 0.0) + time.perf_counter() - started

    def timed_plan(*a, **k):
        started = time.perf_counter()
        plan = original_plan(*a, **k)
        current["flank_s"] = current.get("flank_s", 0.0) + time.perf_counter() - started
        current["hops_max"] = max(current.get("hops_max", 0), int(plan.hops.max(initial=0)))
        current["held"] = current.get("held", 0) + int(plan.held.sum())
        return plan

    rows = []
    rift_damage.damage_pass, rift_strength.flank_plan = timed_pass, timed_plan
    try:
        for k in range(1, args.steps + 1):
            current.clear()
            started = time.perf_counter()
            world_module.step_world(world, args.step_myr * 1e6)
            broken = sum(int(p.collect("rift_broken").astype(bool).sum()) for p in world.plates)
            rows.append({"step": k, "step_s": round(time.perf_counter() - started, 3), "damage_s": round(current.get("damage_s", 0.0), 4),
                         "flank_s": round(current.get("flank_s", 0.0), 4), "hops_max": current.get("hops_max", 0),
                         "held": current.get("held", 0), "broken": broken, "plates": len(world.plates)})
            print(json.dumps(rows[-1]), flush=True)
    finally:
        rift_damage.damage_pass, rift_strength.flank_plan = original_pass, original_plan
    return rows


def dilate(graph: csr_matrix, mask: np.ndarray, hops: int) -> np.ndarray:
    pattern = (graph > 0).astype(np.int8)
    out = mask.copy()
    for _ in range(hops):
        out = out | (pattern @ out.astype(np.int8) > 0)
    return out


def bench_plates(world, args) -> list[dict]:
    """One `flank_plan` per plate and corridor width, on the line the damage pass would choose."""
    p = rift_strength.RiftLineParams.with_overrides(world.rift_line_params)
    rows = []
    for plate in sorted(world.plates, key=lambda q: -q.node_count())[: args.plates]:
        prep = rift_damage.prepare_plate(world, plate, p, args.step_myr * 1e6)
        if prep.line is None:
            continue
        band, _, _ = rift_lines.band_reach(prep.points, prep.line.chain, prep.eps, p)
        in_band = np.zeros(len(prep.points), dtype=bool)
        in_band[band] = True
        for extra in args.extra_hops:
            broken = prep.broken | (dilate(prep.graph, in_band, extra) & ~in_band)
            started = time.perf_counter()
            plan = rift_strength.flank_plan(prep.graph, in_band, broken, p.flank_rows, group=prep.continental)
            seconds = time.perf_counter() - started
            settled = plan.hops[~plan.held]
            rows.append({"plate": int(plate.plate_id), "cells": len(prep.points), "donors": len(plan.band), "extra_hops": extra,
                         "broken": int(broken.sum()), "flank_s": round(seconds, 4), "hops_max": int(plan.hops.max(initial=0)),
                         "hops_mean": round(float(settled.mean()), 1) if len(settled) else 0.0, "held": int(plan.held.sum()),
                         "digest": plan_digest(plan)})
            print(json.dumps(rows[-1]), flush=True)
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--node-density", type=float, default=4.0)
    parser.add_argument("--params", type=json.loads, default={}, help="World.rift_line_params, as JSON")
    parser.add_argument("--steps", type=int, default=0)
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--plates", type=int, default=3)
    parser.add_argument("--extra-hops", type=int, nargs="*", default=[0, 5, 10, 20, 40, 80])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    world = make_world(args)
    report = {"seed": args.seed, "node_density": args.node_density, "params": args.params, "generate_s": round(time.perf_counter() - started, 1)}
    report["steps"] = step_with_timers(world, args)
    report["state_digest"] = state_digest(world)
    report["bench"] = bench_plates(world, args)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=1) + "\n")
    print(json.dumps({"state_digest": report["state_digest"], "out": str(args.out)}))


if __name__ == "__main__":
    main()
