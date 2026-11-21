#!/usr/bin/env python3
"""Craton budget and plate speed through time on a freshly generated world (#1001).

Generates `--seed` at node density 1, steps it `--steps` times by `--step-myr`, and every
`--snapshot-every` steps appends one JSON line to `--out` with:

- the live craton inventory (volume, area, share of genuine continental area, core count);
- the cumulative craton ledger: formed, each sink account, the formed and decratonized area
  counters, and the balance error;
- area-weighted plate speed |omega x r| over continental and oceanic nodes and overall, in
  cm/yr (the measure #942 uses), plus the plate count.

`--formation-multiplier` and `--decratonization-multiplier` set
`World.craton_formation_multiplier` and `World.craton_decratonization_multiplier` after
generation (or after `--resume`); left out, the model defaults (1.0) stand. The glacier march
is off unless `--glacier-march`, as on a world the UI creates: marched steps cost minutes each.

`--save-at-myr` saves the world partway, and `--resume` continues a save, refusing one whose
seed differs. With `--require-no-formation`, it also refuses a save that has formed craton
(`formed_m3 > 0`): arms that differ only in formation step identically until a host first
forms, so they can share a prefix up to such a save. Run each arm in its own process, one
after another, for example:

    PYTHONPATH=backend:bin/debug backend/.venv/bin/python bin/debug/craton_budget_trajectory.py \\
      --seed 326494039 --steps 500 --out <scratch>/seed326494039-default.jsonl

A `.summary.json` beside `--out` holds the final ledger and the runtime.
"""

from __future__ import annotations

import argparse
import json
import resource
import time
from pathlib import Path

import numpy as np

from app import cratons, mantle, persistence
from app.world import generate_world, step_world
from g2_replay_common import add_replay_arguments, prepare_replay_output, write_replay_summary


def plate_speeds_cm_per_yr(world) -> dict[str, float | None]:
    """Area-weighted mean |omega x r| over continental, oceanic and all nodes."""
    sums = {"continental": [0.0, 0.0], "oceanic": [0.0, 0.0], "all": [0.0, 0.0]}
    for plate in world.plates:
        if not plate.node_count():
            continue
        points = plate.all_points_and_elevation()[0]
        points = points / np.linalg.norm(points, axis=1, keepdims=True)
        speed = mantle.rad_per_yr_to_cm_per_yr(1.0) * np.linalg.norm(np.cross(np.asarray(plate.omega, dtype=float), points), axis=1)
        areas = cratons._areas(world, plate)
        continental = cratons._continental(plate)
        for key, mask in (("continental", continental), ("oceanic", ~continental), ("all", np.ones_like(continental))):
            sums[key][0] += float(np.dot(speed[mask], areas[mask]))
            sums[key][1] += float(areas[mask].sum())
    return {f"speed_{key}_cm_per_yr": (num / den if den > 0.0 else None) for key, (num, den) in sums.items()}


def snapshot(world, step: int) -> dict:
    diag = cratons.diagnostics(world)
    ledger = {key: diag[f"craton_{key}"] / 1e15 for key in cratons.CRATON_SOURCE_ACCOUNTS + cratons.CRATON_SINK_ACCOUNTS}
    area_counters = {key: diag.get(f"craton_{key}", 0.0) / 1e12 for key in getattr(cratons, "CRATON_AREA_COUNTERS", ())}
    return {
        "step": step,
        "elapsed_myr": world.elapsed_years / 1e6,
        "craton_formation_multiplier": getattr(world, "craton_formation_multiplier", None),
        "craton_decratonization_multiplier": getattr(world, "craton_decratonization_multiplier", None),
        "plate_count": len(world.plates),
        "craton_volume_mkm3": diag["craton_volume_m3"] / 1e15,
        "craton_area_mkm2": diag["craton_area_m2"] / 1e12,
        "craton_share_of_genuine_continental": diag["craton_area_fraction_of_continental"],
        "craton_core_count": diag["craton_core_count"],
        "ledger_mkm3": ledger,
        "area_counters_mkm2": area_counters,
        "destroyed_mkm3": sum(ledger[key] for key in cratons.CRATON_SINK_ACCOUNTS),
        "balance_error_mkm3": diag["craton_balance_error_m3"] / 1e15,
        "peak_rss_gb": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9,
        **plate_speeds_cm_per_yr(world),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--formation-multiplier", type=float, help="World.craton_formation_multiplier (default: model's)")
    parser.add_argument("--decratonization-multiplier", type=float, help="World.craton_decratonization_multiplier (default: model's)")
    parser.add_argument("--require-no-formation", action="store_true", help="refuse a --resume save that has formed craton")
    parser.add_argument("--glacier-march", action="store_true", help="march glaciers each step (model default; slow)")
    parser.add_argument("--save-at-myr", type=float, help="save the world once elapsed time reaches this")
    parser.add_argument("--save-path", type=Path, help="where --save-at-myr writes the .mbworld")
    parser.add_argument("--resume", type=Path, help="continue a save written by --save-at-myr instead of generating")
    add_replay_arguments(parser, default_steps=500)
    args = parser.parse_args()
    prepare_replay_output(parser, args)
    args.out.write_text("")

    if (args.save_at_myr is None) != (args.save_path is None):
        parser.error("--save-at-myr and --save-path go together")
    started = time.perf_counter()
    if args.resume:
        world = persistence.load_world_bytes(args.resume.read_bytes())
        if world.seed != args.seed:
            parser.error(f"--resume holds seed {world.seed}, not {args.seed}")
        if args.require_no_formation and world.craton_ledger.get("formed_m3", 0.0) != 0.0:
            parser.error("--resume save already formed craton: the arms diverged before it, so it can't be shared")
    else:
        world = generate_world(seed=args.seed, node_density=1.0)
    world.simulate_glacier_march = args.glacier_march
    if args.formation_multiplier is not None:
        world.craton_formation_multiplier = args.formation_multiplier
    if args.decratonization_multiplier is not None:
        world.craton_decratonization_multiplier = args.decratonization_multiplier
    cratons.ensure_initialized(world)
    start_step = round(world.elapsed_years / (args.step_myr * 1e6))
    rows = [snapshot(world, start_step)]
    with args.out.open("a") as handle:
        handle.write(json.dumps(rows[-1], separators=(",", ":")) + "\n")
    for index in range(args.steps):
        step_world(world, args.step_myr * 1e6)
        step = start_step + index + 1
        if args.save_at_myr is not None and args.save_path and world.elapsed_years >= args.save_at_myr * 1e6:
            args.save_path.write_bytes(persistence.save_world_bytes(world))
            print(f"saved {world.elapsed_years / 1e6:g} Myr, formed_m3={world.craton_ledger.get('formed_m3', 0.0):g}", flush=True)
            args.save_path = None
        if step % args.snapshot_every == 0 or index + 1 == args.steps:
            rows.append(snapshot(world, step))
            with args.out.open("a") as handle:
                handle.write(json.dumps(rows[-1], separators=(",", ":")) + "\n")
            row = rows[-1]
            print(
                f"seed={args.seed} step={step} myr={row['elapsed_myr']:g} "
                f"craton={row['craton_volume_mkm3']:.0f} Mkm3 share={row['craton_share_of_genuine_continental']:.3f} "
                f"speed={row['speed_all_cm_per_yr']:.2f} cm/yr runtime_min={(time.perf_counter() - started) / 60:.1f}",
                flush=True,
            )
    write_replay_summary(
        args.out, seed=args.seed, steps=args.steps, step_myr=args.step_myr, rows=rows, started=started,
        ledger=world.craton_ledger,
        craton_formation_multiplier=getattr(world, "craton_formation_multiplier", None),
        craton_decratonization_multiplier=getattr(world, "craton_decratonization_multiplier", None),
        peak_rss_gb=rows[-1]["peak_rss_gb"],
    )


if __name__ == "__main__":
    main()
