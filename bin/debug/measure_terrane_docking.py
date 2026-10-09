#!/usr/bin/env python3
"""Replay a saved world and measure what happens to continental terranes on oceanic plates
(issue #321).

Runs on any branch, so a baseline and a docking branch can be compared from the same save:

- **Terrane survival**: continental-coded area and Hc volume on nominally oceanic plates.
- **Continental inventory**: continental area and Hc volume on every plate, land area.
- **Substrate removal**: Hm booked to the oceanic subduction sink, by node type, and the
  oceanic-coded area left on oceanic plates.
- **Paths**: relocation calls (`quad_tectonics._relocate_terrane_column`, tapped here, so it
  counts on both branches) and, where the branch has them, the docking counters in
  `World.suture_transfer_stats["terranes"]` and the polarity fallback tiers.
- **Ledgers**: continental-material, craton and mobile-cover balance errors.
- **Runtime** per step.

Usage::

    PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_terrane_docking.py \\
        ~/Downloads/mantle-bloom-seed936513024-90000000y.mbworld --myr 20 --out docking.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import continental_ledger, cratons, mobile_cover, persistence, quad_tectonics, world as world_mod  # noqa: E402
from app.elevation_lines import effective_is_continental_from_codes  # noqa: E402

KM2 = 1e6
KM3 = 1e9


def metrics(world) -> dict:
    out = {
        "elapsed_myr": world.elapsed_years / 1e6,
        "terrane_area_mkm2": 0.0,
        "terrane_hc_mkm3": 0.0,
        "terrane_cells": 0,
        "terrane_craton_mkm3": 0.0,
        "craton_mkm3": 0.0,
        "oceanic_plate_oceanic_area_mkm2": 0.0,
        "continental_area_mkm2": 0.0,
        "continental_hc_mkm3": 0.0,
        "land_area_mkm2": 0.0,
    }
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        areas = plate.node_areas_m2()
        hc = plate.collect("crustal_thickness_m")
        continental = effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        land = plate.collect("elevation") > world.sea_level_m
        out["continental_area_mkm2"] += float(areas[continental].sum()) / KM2 / 1e6
        out["continental_hc_mkm3"] += float(hc[continental] @ areas[continental]) / KM3 / 1e6
        out["land_area_mkm2"] += float(areas[land].sum()) / KM2 / 1e6
        craton = plate.collect("craton_crust_m") * areas
        out["craton_mkm3"] += float(craton.sum()) / KM3 / 1e6
        if plate.crust_type == "oceanic":
            out["terrane_craton_mkm3"] += float(craton[continental].sum()) / KM3 / 1e6
            out["terrane_area_mkm2"] += float(areas[continental].sum()) / KM2 / 1e6
            out["terrane_hc_mkm3"] += float(hc[continental] @ areas[continental]) / KM3 / 1e6
            out["terrane_cells"] += int(np.count_nonzero(continental))
            out["oceanic_plate_oceanic_area_mkm2"] += float(areas[~continental].sum()) / KM2 / 1e6
    sink = world.hm_source_sink_ledger.get("oceanic_and_deep_subduction", {}).get("scopes", {})
    out["oceanic_subduction_hm_mkm3"] = {
        scope: sink.get(scope, {}).get("sink_m3", 0.0) / KM3 / 1e6 for scope in ("continental_node", "oceanic_node")
    }
    out["craton_ledger_km3"] = {k: v / KM3 for k, v in (getattr(world, "craton_ledger", {}) or {}).items()}
    out["ledger_errors_km3"] = {
        "continental_material": continental_ledger.balance_error_m3(world) / KM3,
        "craton": cratons.balance_error_m3(world) / KM3,
        "mobile_cover": mobile_cover.balance_error_m3(world) / KM3,
    }
    stats = getattr(world, "suture_transfer_stats", {}) or {}
    out["terranes"] = dict(stats.get("terranes", {}))
    polarity = getattr(world, "collision_polarity_stats", {}) or {}
    out["polarity_fallbacks"] = {k: v for k, v in polarity.items() if k.startswith("fallback_")}
    budget = getattr(world, "orogenic_relief_budget", {}) or {}
    out["orogeny_km3"] = {
        k: budget.get(k, 0.0) / KM3
        for k in ("suture_donated_m3", "suture_lower_crust_subducted_m3", "no_outlet_subducted_m3", "terrane_returned_m3")
    }
    return out


class RelocationTap:
    """Counts `_relocate_terrane_column` calls and the Hc they actually moved."""

    def __init__(self):
        self.calls = 0
        self.relocated = 0
        self.hc_m3 = 0.0

    def __enter__(self):
        self.original = quad_tectonics._relocate_terrane_column

        def wrapped(*args, **kwargs):
            out = self.original(*args, **kwargs)
            self.calls += 1
            if np.any(out):
                self.relocated += 1
                self.hc_m3 += float(args[9])
            return out

        quad_tectonics._relocate_terrane_column = wrapped
        return self

    def __exit__(self, *exc):
        quad_tectonics._relocate_terrane_column = self.original


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("save", type=Path)
    parser.add_argument("--myr", type=float, default=20.0)
    parser.add_argument("--step-years", type=float, default=100_000.0)
    parser.add_argument("--checkpoint-myr", type=float, default=2.0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    commit = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=False).stdout.strip()
    world = persistence.load_world_bytes(args.save.read_bytes())
    world.debug_diagnostics = True
    world.reset_phase_budget()
    continental_ledger.ensure_initialized(world)
    steps = int(round(args.myr * 1e6 / args.step_years))
    every = max(1, int(round(args.checkpoint_myr * 1e6 / args.step_years)))
    result = {"save": str(args.save), "commit": commit, "myr": args.myr, "checkpoints": [metrics(world)], "step_seconds": []}
    args.out.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with RelocationTap() as tap:
        for i in range(1, steps + 1):
            t = time.perf_counter()
            world_mod.step_world(world, args.step_years)
            result["step_seconds"].append(time.perf_counter() - t)
            if i % every == 0 or i == steps:
                row = metrics(world)
                row["relocation_calls"] = tap.calls
                row["relocations"] = tap.relocated
                row["relocated_hc_km3"] = tap.hc_m3 / KM3
                result["checkpoints"].append(row)
                print(
                    f"{i}/{steps} {row['elapsed_myr']:.1f} Myr  terrane {row['terrane_area_mkm2']:.3f} M km2  "
                    f"cont {row['continental_area_mkm2']:.2f}  reloc {tap.relocated}  {row['terranes']}",
                    flush=True,
                )
                args.out.write_text(json.dumps(result, indent=1))
    seconds = np.asarray(result["step_seconds"])
    result["runtime"] = {"total_s": time.time() - started, "mean_step_s": float(seconds.mean()), "p95_step_s": float(np.percentile(seconds, 95))}
    args.out.write_text(json.dumps(result, indent=1))
    print(json.dumps(result["runtime"]))


if __name__ == "__main__":
    main()
