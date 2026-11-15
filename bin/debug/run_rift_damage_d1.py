"""One RFC 0018 D1 run (#888): a generated development world stepped plates-only with the damage-band prototype.

Usage (repo root, this worktree's venv; follow the compute-budget skill first, one heavy run at a time):
    backend/.venv/bin/python bin/debug/run_rift_damage_d1.py --seed 106533960 --out <dir> \
        [--setting '{"eps_dot0": 8e-8, "k_cap": 1.0}'] [--years 200] [--step-myr 1] [--node-density 1.0] \
        [--shadow-every 10] [--snapshot-every 5] [--max-steps N]

Writes ``<dir>/run.json`` (parameters, seed, timings, counters, ledger closure, the recorder's score and every split
event) and, with ``--snapshot-every``, ``<dir>/snapshots/t<Myr>.npz`` of per-cell arrays for the map sequences.
``--shadow-every`` solves the thin shell every that many Myr and records proxy C's D0 statistics (stage 2 reporting).
A held-out or reserve seed from ``analysis/issue797/g2-seed-preregistration.json`` is refused.
"""

import argparse
import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

sys.path[:0] = [str(Path(__file__).resolve().parent)]

import rift_damage_mechanism as mech  # noqa: E402
import rift_damage_recorder as rec  # noqa: E402
import rift_damage_step as step_mod  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


def forbidden_seeds() -> set[int]:
    """The G2 preregistration's held-out and reserve seeds: D1 must never touch them."""
    data = json.loads((REPO / "analysis/issue797/g2-seed-preregistration.json").read_text())
    return {int(seed) for key in ("heldout_seeds", "reserve_seeds") for seed in data[key]}


def snapshot(world, path: Path) -> None:
    arrays = {k: [] for k in ("lon", "lat", "elevation", "craton", "belt", "eps", "broken", "plate", "area")}
    for plate in world.plates:
        pts, _ = plate.all_points_and_elevation()
        arrays["lon"].append(np.degrees(np.arctan2(pts[:, 1], pts[:, 0])).astype(np.float32))
        arrays["lat"].append(np.degrees(np.arcsin(np.clip(pts[:, 2], -1, 1))).astype(np.float32))
        arrays["elevation"].append(plate.collect("elevation").astype(np.float32))
        arrays["craton"].append((plate.collect("craton_crust_m") > 0).astype(np.int8))
        arrays["belt"].append(step_mod.cratons.belt_mask(plate, world.elapsed_years).astype(np.int8))
        arrays["eps"].append(plate.collect("rift_strain").astype(np.float32))
        arrays["broken"].append(plate.collect("rift_broken").astype(np.int8))
        arrays["plate"].append(np.full(len(pts), plate.plate_id, np.int32))
        arrays["area"].append(plate.node_areas_m2().astype(np.float32))
    np.savez_compressed(path, **{k: np.concatenate(v) for k, v in arrays.items()})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--setting", type=json.loads, default={})
    parser.add_argument("--years", type=float, default=200.0, help="Myr")
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--node-density", type=float, default=1.0)
    parser.add_argument("--shadow-every", type=float, default=0.0, help="Myr; 0 = never")
    parser.add_argument("--snapshot-every", type=float, default=0.0, help="Myr; 0 = never")
    parser.add_argument("--max-steps", type=int)
    parser.add_argument("--production-rate", type=float, help="the production splitter's splits per world per 100 Myr on this world")
    args = parser.parse_args()
    if args.seed in forbidden_seeds():
        sys.exit(f"seed {args.seed} is in the G2 preregistration (held-out or reserve); D1 must not use it")
    args.out.mkdir(parents=True, exist_ok=True)
    params = mech.Params(**args.setting)
    state = step_mod.install(params)
    from app import world as world_module

    started = time.perf_counter()
    world = world_module.generate_world(args.seed, node_density=args.node_density)
    world.simulate_climate_biomes = False
    generated_s = time.perf_counter() - started
    recorder = rec.Recorder(params)
    steps = round(args.years / args.step_myr)
    if args.max_steps:
        steps = min(steps, args.max_steps)
    shadow, plate_counts = [], []
    stepping = time.perf_counter()
    if args.snapshot_every:
        (args.out / "snapshots").mkdir(exist_ok=True)
        snapshot(world, args.out / "snapshots" / "t000.npz")
    for k in range(1, steps + 1):
        world_module.step_world(world, args.step_myr * 1e6)
        recorder.observe(world, state, args.step_myr * 1e6)
        myr = k * args.step_myr
        plate_counts.append((myr, len(world.plates)))
        if args.snapshot_every and round(myr / args.snapshot_every) * args.snapshot_every == myr:
            snapshot(world, args.out / "snapshots" / f"t{int(myr):03d}.npz")
        if args.shadow_every and round(myr / args.shadow_every) * args.shadow_every == myr:
            import shadow_shell_check

            shadow.append({"myr": myr, **shadow_shell_check.shadow_stats(world)})
        if k % 10 == 0:
            print(json.dumps({"myr": myr, "plates": len(world.plates), "splits": len(state.events), "counters": state.counters,
                              "wall_s": round(time.perf_counter() - stepping, 1)}), flush=True)
    recorder.finish(world)
    score = rec.score(list(recorder.candidates.values()), state.events, recorder.last_years, recorder.super_sea_max_width_m,
                      recorder.super_sea_violations, params.a_min, args.production_rate, steps * args.step_myr)
    result = {
        "seed": args.seed, "params": asdict(params), "years_myr": steps * args.step_myr, "step_myr": args.step_myr,
        "node_density": args.node_density, "generate_s": generated_s, "step_s": time.perf_counter() - stepping,
        "line": {"strain_on_line_share": state.on_line_strain / max(state.total_strain, 1e-30), "growing_samples": len(state.lam_samples),
                 "open_samples": sum(o for _, o, _ in state.lam_samples),
                 "lam_quartiles": [float(x) for x in np.percentile([lam for lam, _, _ in state.lam_samples], [25, 50, 75])] if state.lam_samples else None},
        "seconds_damage": state.seconds_damage, "seconds_split_decision": state.seconds_split, "counters": state.counters,
        "ledger_relative_gap_max": state.ledger_gap, "plates_final": len(world.plates), "plate_counts": plate_counts[::10],
        "score": score, "events": [{k: v for k, v in e.items() if k != "path_keys"} for e in state.events], "shadow_shell": shadow,
        "candidates": [{"id": c.ident, "plate": c.plate_id, "onset_myr": c.onset_years / 1e6, "onset_width_km": c.onset_width_m / 1e3,
                        "status": c.status, "from_start": c.from_start, "max_cells": c.max_cells, "end_myr": None if c.end_years is None else c.end_years / 1e6,
                        "ever_failed": c.ever_failed, "samples": [[round(s[0] / 1e6, 2), s[1], round(s[4] / 1e3, 1), s[5], round(s[6], 2)] for s in c.samples[::5]]}
                       for c in recorder.candidates.values() if c.max_cells >= rec.MIN_CELLS_TO_COUNT],
    }
    (args.out / "run.json").write_text(json.dumps(result, indent=1, default=float) + "\n")
    print(json.dumps({"done": True, "targets_passed": score["targets_passed"], "splits": score["splits"]["value"], "step_s": result["step_s"]}))


if __name__ == "__main__":
    main()
