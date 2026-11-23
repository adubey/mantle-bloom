"""Why do the biggest plates rarely rift? Per-plate rift-line state over time (#995).

Steps one generated world plates-only under one RFC 0021 setting (the `plate_contact_graph.py record` set-up) and, every
`--every` Myr, records for every plate with a rift line (`rift_damage.runtime(world).steps`):

- area share of the planet, craton share of its area, age (the split age gate is 2 Myr);
- the line's pull `P` and strength `S` at the step end, the ratio P/S (the breakthrough test is `P >= S`), the broken
  length share (needs `f_break` = 0.5), whether the plate has a line at all;
- which gate fails: no line / not enough of the line broken / pull below strength / passes both.

At the end it writes the counters (`break_tests`, `split_rejected` and why) and the split events.

Usage (repo root; compute-budget first; about 1 GB, plates-only):
    backend/.venv/bin/python bin/debug/big_plate_rift_probe.py --settings analysis/issue995/settings.json \\
        --setting brake-o0.25-c1 --seed 591453549 --years 200 --every 5 --out analysis/issue995/rift-state-<seed>.json
"""

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace


HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "backend")]


def gate(step, f_break: float) -> str:
    if step is None or step.line is None:
        return "no_line"
    if step.broken_share < f_break:
        return "too_little_broken"
    if not (step.pull_end > 0.0 and step.pull_end >= step.strength_end):
        return "pull_below_strength"
    return "passes"


def sample(world, myr: float) -> dict:
    from app import rift_damage

    state = rift_damage.runtime(world)
    plates = [p for p in world.plates if p.node_count() > 0]
    areas = {p.plate_id: float(p.node_areas_m2().sum()) for p in plates}
    total = sum(areas.values())
    rows = []
    for plate in plates:
        step = state.steps.get(plate.plate_id)
        node_area = plate.node_areas_m2()
        craton = float(node_area[plate.collect("craton_crust_m") > 0].sum() / node_area.sum())
        row = {"plate": plate.plate_id, "area_share": areas[plate.plate_id] / total, "craton_share": craton,
               "age_myr": plate.age_years / 1e6, "nodes": plate.node_count(), "gate": gate(step, state.params.f_break)}
        if step is not None and step.line is not None:
            row.update(pull_n=float(step.pull_end), strength_n=float(step.strength_end), broken_share=float(step.broken_share),
                       line_open=bool(step.line.is_open), peak_share=float(step.line.peak_share))
        rows.append(row)
    return {"myr": myr, "plates": rows}


def install_ramp(early: float, late: float, tau_myr: float, big: float = 0.10) -> dict:
    """Time-dependent strength for big plates (a #995 probe, not a model): the breakthrough test sees the line strength times
    `late + (early - late) * exp(-t / tau)`, with t the Myr the plate has held at least `big` of the planet's area.
    Plates under `big` are untouched, and damage accumulation is untouched (only the P >= S test is scaled). Returns the
    `big_since` dict, which the caller updates each step with `track_big`."""
    import math

    from app import rift_damage

    big_since: dict[int, float] = {}
    original = rift_damage.breakthrough

    def ramped(world, plate):
        step = rift_damage.runtime(world).steps.get(plate.plate_id)
        t0 = big_since.get(plate.plate_id)
        if step is None or t0 is None:
            return original(world, plate)
        saved = step.strength_end
        step.strength_end = saved * (late + (early - late) * math.exp(-(world.elapsed_years / 1e6 - t0) / tau_myr))
        try:
            return original(world, plate)
        finally:
            step.strength_end = saved

    rift_damage.breakthrough = ramped
    return big_since


def track_big(world, big_since: dict, myr: float, big: float = 0.10) -> None:
    areas = {p.plate_id: float(p.node_areas_m2().sum()) for p in world.plates if p.node_count() > 0}
    total = sum(areas.values())
    for pid in list(big_since):
        if areas.get(pid, 0.0) / total < big:
            del big_since[pid]
    for pid, a in areas.items():
        if a / total >= big:
            big_since.setdefault(pid, myr)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--settings", type=Path, required=True)
    ap.add_argument("--setting", required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--years", type=float, default=200.0)
    ap.add_argument("--step-myr", type=float, default=1.0)
    ap.add_argument("--every", type=float, default=5.0)
    ap.add_argument("--ramp", type=float, nargs=3, metavar=("EARLY", "LATE", "TAU_MYR"), help="strength factor ramp for big plates")
    ap.add_argument("--scratch", required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    import calibrate_plate_forces as cal
    import rift_line_recorder as rec
    import run_rift_lines_d1 as d1
    from app import rift_damage

    setting = cal.load_setting(args.settings, args.setting)
    if args.seed in cal.refused_seeds():
        sys.exit("held-out seed")
    world = d1.make_world(SimpleNamespace(seed=args.seed, node_density=1.0, world_knobs=setting["world_knobs"],
                                          split_method=setting["split_method"], params=setting["params"]))
    cal.apply_module_knobs(setting["module_knobs"])
    recorder = rec.Recorder(rift_damage.runtime(world).params.eps_cand)
    samples = []

    big_since = install_ramp(*args.ramp) if args.ramp else None

    def on_step(world, myr: float) -> None:
        if big_since is not None:
            track_big(world, big_since, myr)
        if round(myr / args.every) * args.every == myr:
            samples.append(sample(world, myr))

    scratch = Path(args.scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    splits, *_ = d1.step_all(world, SimpleNamespace(years=args.years, step_myr=args.step_myr, max_steps=None, snapshot_every=0.0), scratch,
                             recorder, on_step)
    state = rift_damage.runtime(world)
    args.out.write_text(json.dumps({"seed": args.seed, "setting": args.setting, "counters": state.counters, "splits": splits, "merges": list(d1.MERGES), "ramp": args.ramp,
                                    "samples": samples}, default=float) + "\n")


if __name__ == "__main__":
    main()
