"""One RFC 0018 D1 run on the production split methods (#888, #846): a generated development world stepped
plates-only, with `World.split_method = "rift_lines"` or the great-circle splitter.

Usage (repo root, this worktree's venv; follow the compute-budget skill first; about 5-7 minutes per world):
    backend/.venv/bin/python bin/debug/run_rift_lines_d1.py --seed 106533960 --out <dir> \
        [--split-method rift_lines|great_circle] [--params '{"k_cap": 0.14, "mark_power": 2.8}'] \
        [--world-knobs '{"contact_drag_multiplier": 0.0}'] [--years 200] [--step-myr 1] [--snapshot-every 5] [--max-steps N]

``--params`` becomes `World.rift_line_params` (checked by `RiftLineParams.with_overrides`); ``--world-knobs`` sets
other tuning attributes on the generated world, such as the contact brake.

Writes ``<dir>/run.json``: the settings, timings, every split (time, plate, smaller daughter's share) and the plate
count. For rift lines it also has the method's counters, its ledger closure, the run's own target score, the split
events, the candidates, and every 10 Myr the area-weighted speed of continental and oceanic crust and how much of the
step's rift lines is oceanic. ``<dir>/recorder.pkl`` holds the recorder's candidates and events for
`score_rift_lines.py`. With ``--snapshot-every``, ``<dir>/snapshots/t<Myr>.npz`` holds per-cell arrays for maps.

A held-out or reserve seed from ``analysis/issue797/g2-seed-preregistration.json`` is refused: D1 must never touch them.
"""

import argparse
import json
import pickle
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

sys.path[:0] = [str(Path(__file__).resolve().parent), str(Path(__file__).resolve().parents[2] / "backend")]

import rift_line_recorder as rec  # noqa: E402
from app import cratons, merge_split, rift_damage  # noqa: E402
from app import world as world_module  # noqa: E402
from app.rift_strength import RADIUS_M  # noqa: E402

REPO = Path(__file__).resolve().parents[2]


def forbidden_seeds() -> set[int]:
    """The G2 preregistration's held-out and reserve seeds."""
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
        arrays["belt"].append(cratons.belt_mask(plate, world.elapsed_years).astype(np.int8))
        arrays["eps"].append(plate.collect("rift_strain").astype(np.float32))
        arrays["broken"].append(plate.collect("rift_broken").astype(np.int8))
        arrays["plate"].append(np.full(len(pts), plate.plate_id, np.int32))
        arrays["area"].append(plate.node_areas_m2().astype(np.float32))
    np.savez_compressed(path, **{k: np.concatenate(v) for k, v in arrays.items()})


def kinematics(world) -> dict:
    """Area-weighted mean speed (cm/yr) of continental and oceanic crust over the plates the damage pass saw, the
    continental share of their area, and the length-weighted oceanic share and count of the step's rift lines (#888:
    does the brake keep continents slow and push the lines into the ocean)."""
    state = rift_damage.runtime(world)
    moved = {"continental": [0.0, 0.0], "oceanic": [0.0, 0.0]}
    lines, line_oceanic, line_weight, growing = 0, 0.0, 0.0, 0
    by_id = {plate.plate_id: plate for plate in world.plates}
    for plate_id, step in state.steps.items():
        plate = by_id.get(plate_id)
        if plate is None or len(step.points) != len(step.oceanic):
            continue
        speed = np.linalg.norm(np.cross(plate.omega, step.points), axis=1) * RADIUS_M * 100.0
        for name, mask in (("continental", ~step.oceanic), ("oceanic", step.oceanic)):
            moved[name][0] += float(np.dot(speed[mask], step.areas[mask]))
            moved[name][1] += float(step.areas[mask].sum())
        if step.line is not None:
            ds = step.line.geometry.ds
            lines += 1
            line_oceanic += float(ds[step.oceanic[step.line.chain]].sum())
            line_weight += float(ds.sum())
            growing += step.line.is_open or step.line.lam > state.params.lam_grow
    total = moved["continental"][1] + moved["oceanic"][1]
    return {"speed_cm_per_yr": {k: v[0] / v[1] if v[1] else None for k, v in moved.items()},
            "continental_area_share": moved["continental"][1] / total if total else None, "lines": lines, "lines_growing": int(growing),
            "line_oceanic_share": line_oceanic / line_weight if line_weight else None}


def candidate_rows(recorder: rec.Recorder) -> list[dict]:
    return [{"id": c.ident, "plate": c.plate_id, "onset_myr": c.onset_years / 1e6, "onset_width_km": c.onset_width_m / 1e3,
             "status": c.status, "from_start": c.from_start, "max_cells": c.max_cells, "end_myr": None if c.end_years is None else c.end_years / 1e6,
             "ever_failed": c.ever_failed, "samples": [[round(s[0] / 1e6, 2), s[1], round(s[4] / 1e3, 1), s[5], round(s[6], 2)] for s in c.samples[::5]]}
            for c in recorder.candidates.values() if c.max_cells >= rec.MIN_CELLS_TO_COUNT]


def make_world(args):
    world = world_module.generate_world(args.seed, node_density=args.node_density)
    world.simulate_climate_biomes = False
    for name, value in args.world_knobs.items():
        if not hasattr(world, name):
            sys.exit(f"the world has no knob {name!r}")
        setattr(world, name, value)
    world.split_method = args.split_method
    world.rift_line_params = dict(args.params)
    merge_split.check_split_settings(world)  # a bad --params fails here, before any stepping
    return world


MERGES: list[dict] = []  # #995: every merge applied during `step_all` (the splits list has the splits); cleared per call


def step_all(world, args, out: Path, recorder, on_step=None):
    """Step the world, recording every split the splitter applies. Returns (splits, plate counts, kinematics).
    `on_step(world, myr)`, if given, runs after each step's own recording (RFC 0021's runner, #958)."""
    splits, plate_counts, motion = [], [], []
    original = merge_split.maybe_split_plate
    original_merge = merge_split.merge_plates
    MERGES.clear()

    def counting_merge(world, id_keep, id_absorb):
        total = sum(float(p.node_areas_m2().sum()) for p in world.plates)
        areas = {p.plate_id: float(p.node_areas_m2().sum()) for p in world.plates if p.plate_id in (id_keep, id_absorb)}
        original_merge(world, id_keep, id_absorb)
        MERGES.append({"myr": world.elapsed_years / 1e6, "keep": int(id_keep), "absorb": int(id_absorb),
                       "keep_share": areas.get(id_keep, 0.0) / total, "absorb_share": areas.get(id_absorb, 0.0) / total})

    def counting(world, plate, step_years=0.0):
        result = original(world, plate, step_years)
        if result is not None:
            areas = [float(d.node_areas_m2().sum()) for d in result]
            daughters = {id(d) for d in result}
            others = [float(o.node_areas_m2().sum()) for o in world.plates if id(o) not in daughters and o is not plate]
            parent = sum(areas)
            splits.append({"myr": world.elapsed_years / 1e6, "plate": int(plate.plate_id), "small_share": min(areas) / sum(areas),
                           # #995: the splitting plate's size against the rest, to ask whether the big plates rift
                           "parent_area_m2": parent, "parent_rank": 1 + sum(a > parent for a in others), "n_plates_before": len(others) + 1,
                           "largest_other_m2": max(others, default=0.0), "total_area_m2": parent + sum(others),
                           "daughters": [int(d.plate_id) for d in result]})
        return result

    steps = round(args.years / args.step_myr)
    if args.max_steps:
        steps = min(steps, args.max_steps)
    merge_split.maybe_split_plate = counting
    merge_split.merge_plates = counting_merge
    started = time.perf_counter()
    try:
        for k in range(1, steps + 1):
            world_module.step_world(world, args.step_myr * 1e6)
            myr = k * args.step_myr
            plate_counts.append((myr, len(world.plates)))
            if recorder is not None:
                recorder.observe(world, args.step_myr * 1e6)
            if on_step is not None:
                on_step(world, myr)
            if args.snapshot_every and round(myr / args.snapshot_every) * args.snapshot_every == myr:
                snapshot(world, out / "snapshots" / f"t{int(myr):03d}.npz")
            if k % 10 == 0:
                if recorder is not None:
                    motion.append({"myr": myr, **kinematics(world)})
                print(json.dumps({"myr": myr, "plates": len(world.plates), "splits": len(splits), "wall_s": round(time.perf_counter() - started, 1)}), flush=True)
    finally:
        merge_split.maybe_split_plate = original
        merge_split.merge_plates = original_merge
    return splits, plate_counts, motion, steps * args.step_myr, time.perf_counter() - started


def rift_line_result(world, recorder: rec.Recorder, years_myr: float, out: Path) -> dict:
    state = rift_damage.runtime(world)
    recorder.finish(world)
    score = rec.score(list(recorder.candidates.values()), state.events, recorder.last_years, recorder.super_sea_max_width_m,
                      recorder.super_sea_violations, state.params.a_min, None, years_myr)
    with open(out / "recorder.pkl", "wb") as handle:  # full candidates and events, for scoring several runs together
        pickle.dump({"candidates": list(recorder.candidates.values()), "events": state.events, "last_years": recorder.last_years,
                     "super_sea_max_m": recorder.super_sea_max_width_m, "super_sea_violations": recorder.super_sea_violations,
                     "years_myr": years_myr}, handle)
    return {"params": asdict(state.params), "counters": state.counters, "ledger_relative_gap_max": state.ledger_gap, "score": score,
            "events": [{k: v for k, v in e.items() if k != "path_keys"} for e in state.events], "candidates": candidate_rows(recorder)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--split-method", choices=rift_damage.SPLIT_METHODS, default="rift_lines")
    parser.add_argument("--params", type=json.loads, default={}, help="World.rift_line_params, as JSON")
    parser.add_argument("--world-knobs", type=json.loads, default={}, help="other tuning attributes set on the generated world, as JSON")
    parser.add_argument("--years", type=float, default=200.0, help="Myr")
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--node-density", type=float, default=1.0)
    parser.add_argument("--snapshot-every", type=float, default=0.0, help="Myr; 0 = never")
    parser.add_argument("--max-steps", type=int)
    args = parser.parse_args()
    if args.seed in forbidden_seeds():
        sys.exit(f"seed {args.seed} is in the G2 preregistration (held-out or reserve); D1 must not use it")
    args.out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    world = make_world(args)
    generated_s = time.perf_counter() - started
    lines_on = args.split_method == "rift_lines"
    recorder = rec.Recorder(rift_damage.runtime(world).params.eps_cand) if lines_on else None
    if args.snapshot_every:
        (args.out / "snapshots").mkdir(exist_ok=True)
        snapshot(world, args.out / "snapshots" / "t000.npz")
    splits, plate_counts, motion, years_myr, step_s = step_all(world, args, args.out, recorder)
    result = {"seed": args.seed, "split_method": args.split_method, "rift_line_params": args.params, "world_knobs": args.world_knobs,
              "years_myr": years_myr, "step_myr": args.step_myr, "node_density": args.node_density, "generate_s": generated_s, "step_s": step_s,
              "n_splits": len(splits), "splits": splits, "plates_final": len(world.plates), "plate_counts": plate_counts[::10], "kinematics": motion}
    if lines_on:
        result.update(rift_line_result(world, recorder, years_myr, args.out))
    (args.out / "run.json").write_text(json.dumps(result, indent=1, default=float) + "\n")
    print(json.dumps({"done": True, "splits": len(splits), "step_s": round(step_s),
                      **({"targets_passed": result["score"]["targets_passed"]} if lines_on else {})}))


if __name__ == "__main__":
    main()
