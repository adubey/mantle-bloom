"""Run and score RFC 0021 settings: plate forces, RFC 0008's collision measures and RFC 0018's rifting score (#958).

One run steps a generated development world plates-only under one setting and records, in ``<out>/<setting>/<seed>/``:

- ``run.json`` and ``recorder.pkl``: `run_rift_lines_d1.py`'s split record and, for rift lines, RFC 0018's recorder;
- ``breakdown.json``: every ``--every`` Myr, `plate_torque_breakdown.sample` (force shares along each plate's axis,
  speeds, the solver's closure and RFC 0021's G0 measures), and every step each plate's ``|omega| R`` and area, for
  T6's pre- and post-contact speeds;
- ``contacts.jsonl``: `probe_craton_contacts.py`'s per-step contact rows (RFC 0008: jams, collision-time drops,
  fusions and their craton overlap, plate count, the largest plate). With ``--no-contacts`` it is not written.

Scoring pools one setting's worlds into one JSON: T1-T5 and T7 from the breakdowns, T8's measures from the contact
rows, and T9 from RFC 0018's `rift_line_recorder.pool`, against the production splitter's rate on the same worlds.

Usage (repo root, this worktree's venv; follow the compute-budget skill; 5-7 minutes and about 0.5 GB per world):
    backend/.venv/bin/python bin/debug/calibrate_plate_forces.py queue --settings analysis/issue958/settings.json \\
        --setting rifting-baseline-brake-off --seeds-file analysis/issue958/development-seeds.json \\
        --out ~/mantle-bloom-captures/issue958 [--step-myr 1] [--years 200] [--every 10] [--no-contacts] [--dry-run]
    backend/.venv/bin/python bin/debug/calibrate_plate_forces.py run --settings <file> --setting <name> --seed <seed> --out <dir>
    backend/.venv/bin/python bin/debug/calibrate_plate_forces.py score --runs ~/mantle-bloom-captures/issue958/<setting> \\
        [--production-runs ~/mantle-bloom-captures/issue958/great-circle | --production-rate 7.7] [--t2-share 0.05] \\
        [--out analysis/issue958/<setting>.json]
    backend/.venv/bin/python bin/debug/calibrate_plate_forces.py t2-census --seeds-file analysis/issue958/development-seeds.json \\
        --out analysis/issue958/t2-census.json

A settings file maps a setting name to ``split_method`` (default ``rift_lines``), ``params`` (`World.rift_line_params`)
and ``world_knobs``, with ``common_params`` merged under every rift-line setting's own. ``queue`` runs each seed in its
own process, one at a time, and skips a seed whose ``run.json`` exists, so it resumes. The G2 preregistration's
held-out and reserve seeds and RFC 0021's held-out seeds are refused unless ``--heldout`` (G3 only).
"""

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
sys.path[:0] = [str(HERE), str(REPO / "backend")]

HELDOUT_FILE = REPO / "analysis/issue958/heldout-seeds.json"
T2_SHARE_FILE = REPO / "analysis/issue958/t2-threshold.json"
T5_MIN_FLOOR_AGE_MYR = 50.0  # RFC 0021 T5: B.3's column reaches 81 of its 84 km at 50 Myr


def load_setting(path: Path, name: str) -> dict:
    spec = json.loads(path.read_text())
    setting = spec["settings"][name]
    method = setting.get("split_method", "rift_lines")
    params = {**spec.get("common_params", {}), **setting.get("params", {})} if method == "rift_lines" else {}
    return {"name": name, "split_method": method, "params": params, "world_knobs": setting.get("world_knobs", {})}


def refused_seeds() -> set[int]:
    import run_rift_lines_d1 as d1

    seeds = d1.forbidden_seeds()
    if HELDOUT_FILE.exists():
        seeds |= {int(s) for s in json.loads(HELDOUT_FILE.read_text())["heldout_seeds"]}
    return seeds


# --- one run ---------------------------------------------------------------------------------------------------------


class _Taps:
    """The contact probe's tapped functions also see the breakdown's read-only solve; this mutes them during it."""

    def __init__(self, probe) -> None:
        self.probe = probe

    def __enter__(self):
        if self.probe is not None:
            self.probe.on_slab_pull = lambda *a, **k: None
            self.probe.on_collision_drag = lambda *a, **k: None
        return self

    def __exit__(self, *exc) -> None:
        if self.probe is not None:
            del self.probe.on_slab_pull, self.probe.on_collision_drag


def _plate_speeds(world, radius_cm: float) -> list[list[float]]:
    return [[int(p.plate_id), float(np.linalg.norm(p.omega) * radius_cm), float(p.node_areas_m2().sum() / 1e6)]
            for p in world.plates if p.node_count() > 0]


def run_one(args: argparse.Namespace) -> None:
    import plate_torque_breakdown as breakdown
    import probe_craton_contacts as contacts
    import rift_line_recorder as rec
    import run_rift_lines_d1 as d1
    from app import cratons, lithosphere, rift_damage

    if args.seed in refused_seeds() and not args.heldout:
        sys.exit(f"seed {args.seed} is held out (G2 preregistration or RFC 0021); only G3 may run it, with --heldout")
    setting = load_setting(args.settings, args.setting)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    world = d1.make_world(SimpleNamespace(seed=args.seed, node_density=1.0, world_knobs=setting["world_knobs"],
                                          split_method=setting["split_method"], params=setting["params"]))
    recorder = rec.Recorder(rift_damage.runtime(world).params.eps_cand) if setting["split_method"] == "rift_lines" else None
    radius_cm = lithosphere.PLANET_RADIUS_M * 100.0
    samples, speeds = [], []
    probe, handle = None, None
    if not args.no_contacts:
        probe = contacts.Probe(world, 500.0)
        contacts.install(probe)
        handle = open(out / "contacts.jsonl", "w")
        handle.write(json.dumps({"seed": args.seed, "step_myr": args.step_myr, "reach_km": 500.0, "start_craton_m3": cratons.live_volume_m3(world)}) + "\n")

    def on_step(world, myr: float) -> None:
        if probe is not None:
            for row in probe.rows:
                handle.write(json.dumps(row) + "\n")
            probe.rows = []
        speeds.append([myr, _plate_speeds(world, radius_cm)])
        if round(myr / args.every) * args.every == myr:
            with _Taps(probe):
                samples.append(breakdown.sample(world))

    stepping = SimpleNamespace(years=args.years, step_myr=args.step_myr, max_steps=None, snapshot_every=0.0)
    splits, plate_counts, motion, years_myr, step_s = d1.step_all(world, stepping, out, recorder, on_step)
    if probe is not None:
        probe.finish()
        for row in probe.rows:
            handle.write(json.dumps(row) + "\n")
        handle.write(json.dumps({"end_craton_m3": cratons.live_volume_m3(world), "ledger": dict(world.craton_ledger), "wall_s": time.perf_counter() - started}) + "\n")
        handle.close()
    meta = {"seed": args.seed, "setting": setting, "years": years_myr, "step_myr": args.step_myr, "every": args.every}
    (out / "breakdown.json").write_text(json.dumps({**meta, "world_knobs": setting["world_knobs"], "samples": samples, "plate_speeds": speeds}, default=float) + "\n")
    result = {**meta, "split_method": setting["split_method"], "years_myr": years_myr, "step_s": step_s, "wall_s": time.perf_counter() - started,
              "n_splits": len(splits), "splits": splits, "plates_final": len(world.plates), "plate_counts": plate_counts[::10], "kinematics": motion,
              "ledger_relative_gap_max": {}}
    if recorder is not None:
        result.update(d1.rift_line_result(world, recorder, years_myr, out))
    (out / "run.json").write_text(json.dumps(result, indent=1, default=float) + "\n")
    print(json.dumps({"seed": args.seed, "splits": len(splits), "step_s": round(step_s)}), flush=True)


def queue(args: argparse.Namespace) -> None:
    seeds = args.seeds or json.loads(args.seeds_file.read_text())["seeds"]
    failed = []
    for seed in seeds:
        target = args.out / args.setting / str(seed)
        if (target / "run.json").exists():
            continue
        command = [sys.executable, __file__, "run", "--settings", str(args.settings), "--setting", args.setting, "--seed", str(seed),
                   "--out", str(target), "--years", str(args.years), "--step-myr", str(args.step_myr), "--every", str(args.every)]
        command += ["--no-contacts"] * args.no_contacts + ["--heldout"] * args.heldout
        if args.dry_run:
            print(" ".join(command[1:]))
            continue
        target.mkdir(parents=True, exist_ok=True)
        with open(target / "run.log", "w") as log:
            code = subprocess.call(command, stdout=log, stderr=subprocess.STDOUT)
        with open(args.out / "queue.log", "a") as log:
            log.write(f"{args.setting}/{seed} exit {code} {time.strftime('%H:%M')}\n")
        if code:
            failed.append(seed)
    print(json.dumps({"setting": args.setting, "seeds": len(seeds), "failed": failed}))
    sys.exit(1 if failed else 0)


# --- scoring ---------------------------------------------------------------------------------------------------------


def _weighted(pairs: list[tuple[float, float]]) -> float | None:
    weight = sum(w for _, w in pairs)
    return sum(v * w for v, w in pairs) / weight if weight else None


def slab_attached(row: dict, share: float) -> bool:
    """T2's definition (the owner's, 2026-11-18, after B.4b's generation census failed): the plate has a slab of its
    own, oceanic floor it lost under its neighbours in the last 10 Myr of at least `share` of its area."""
    return row["g0"]["consumed_10myr_km2"] >= share * row["area_km2"]


def force_targets(rows: list[dict], share: float) -> dict:
    """T1-T5 pooled over plate-samples, area-weighted (T5 length-weighted)."""
    import summarize_torque_breakdown as summary

    pooled = summary.summarize(rows)["all"]
    attached = [(r["speed_cm_per_yr"]["all"], r["area_km2"]) for r in rows if slab_attached(r, share)]
    other = [(r["speed_cm_per_yr"]["all"], r["area_km2"]) for r in rows if not slab_attached(r, share)]
    fast, slow = _weighted(attached), _weighted(other)
    slabless = [(r["g0"]["basal_traction_pa"], r["area_km2"]) for r in rows if not r["g0"]["has_subducting_node"]]
    old = [(r["g0"]["ridge_push_n_per_m"], r["g0"]["ridge_push_km"]) for r in rows
           if r["g0"]["ridge_push_n_per_m"] is not None and (r["g0"]["oceanic_age_max_myr"] or 0.0) >= T5_MIN_FLOOR_AGE_MYR]
    return {
        "T1_speed_cm_per_yr": _weighted([(r["speed_cm_per_yr"]["all"], r["area_km2"]) for r in rows]),
        # In the no-net-rotation frame, as T1's Earth figure is (`earth_plate_speed.py`); older runs lack it.
        "T1_speed_nnr_cm_per_yr": _weighted([(r["speed_nnr_cm_per_yr"], r["area_km2"]) for r in rows if "speed_nnr_cm_per_yr" in r]),
        "T2_slab_attached_ratio": fast / slow if fast is not None and slow else None,
        "T2_slab_attached_area_share": sum(w for _, w in attached) / sum(r["area_km2"] for r in rows) if rows else None,
        "T3_slab_pull_driving_share": pooled["driving_share"].get("slab_pull"),
        "T4_slabless_basal_traction_mpa": None if _weighted(slabless) is None else _weighted(slabless) / 1e6,
        "T5_old_ridge_push_n_per_m": _weighted(old),
        "breakdown": pooled,
    }


def _world_speeds(runs: list[dict]) -> dict:
    speeds = {str(run["seed"]): _weighted([(r["speed_cm_per_yr"]["all"], r["area_km2"]) for s in run["samples"] for r in s["plates"]]) for run in runs}
    values = [v for v in speeds.values() if v]
    return {"per_world_cm_per_yr": speeds, "T7_fastest_over_slowest": max(values) / min(values) if values else None}


def _closure(runs: list[dict]) -> dict:
    ratios = [s["solver"]["residual_torque_norm_nm"] / s["solver"]["driving_torque_norm_nm"]
              for run in runs for s in run["samples"] if s["solver"]["driving_torque_norm_nm"]]
    return {"residual_over_driving_max": max(ratios, default=None),
            "solver_converged_all": all(s["solver"]["solver_converged"] for run in runs for s in run["samples"])}


def contact_measures(directories: list[Path]) -> dict:
    """T8: RFC 0008's measures pooled over worlds, from the probe's per-seed summary."""
    import probe_craton_contacts as contacts

    seeds = [contacts.summarize_seed(d / "contacts.jsonl") for d in directories if (d / "contacts.jsonl").exists()]
    if not seeds:
        return {}
    return {
        "worlds": len(seeds),
        "collision_time_drops_jammed": sum(s["collision_time_drops_jammed"] for s in seeds),
        "craton_overlap_at_fusion_km3": sum(s["fusions"]["craton_overlap_km3"] for s in seeds),
        "largest_plate_share_mean": float(np.mean([s["plates"]["largest_share_mean"] for s in seeds])),
        "plate_count_mean": float(np.mean([s["plates"]["count_mean"] for s in seeds])),
        "fusions": sum(s["fusions"]["roll"] + s["fusions"]["forced"] for s in seeds),
        "per_world": seeds,
    }


def production_rate(directory: Path) -> float:
    runs = [json.loads((d / "run.json").read_text()) for d in sorted(directory.iterdir()) if (d / "run.json").exists()]
    return sum(r["n_splits"] for r in runs) / sum(r["years_myr"] for r in runs) * 100.0


def rifting_score(directory: Path, rate: float | None) -> dict | None:
    import rift_line_recorder as rec
    import score_rift_lines
    from app.rift_strength import RiftLineParams

    if not all((d / "recorder.pkl").exists() for d in directory.iterdir() if d.is_dir()):
        return None
    return score_rift_lines.row(directory.name, rec.pool(directory, rate, RiftLineParams().a_min))


def t2_census(args: argparse.Namespace) -> None:
    """Per plate of each generated world (t = 0, no stepping), the share of its boundary-band length that B.4b's test
    calls slab-attached, with its area and closing length at any rate: the distribution T2's threshold is fixed on."""
    import plate_torque_breakdown as breakdown
    from app import world as world_module

    seeds = json.loads(args.seeds_file.read_text())["seeds"]
    worlds = {}
    for seed in seeds:
        if seed in refused_seeds():
            sys.exit(f"seed {seed} is held out")
        world = world_module.generate_world(seed)
        world.simulate_climate_biomes = False
        rows = breakdown.sample(world)["plates"]
        worlds[str(seed)] = [{"plate_id": r["plate_id"], "area_km2": round(r["area_km2"]), "continental_share": round(r["continental_area_share"], 3),
                              "band_km": round(r["boundary_km"]["band_total"]), "geometric_slab_km": round(r["g0"]["geometric_slab_km"]),
                              "closing_down_going_km": round(r["boundary_km"]["down_going_closing_any_rate"]),
                              "past_gate_km": round(r["boundary_km"]["down_going_closing_past_gate"]),
                              "rate_cm_per_yr": round(r["rate_cm_per_yr"], 3)} for r in rows]
        print(seed, len(rows), flush=True)
    text = json.dumps({"seeds": seeds, "worlds": worlds}, indent=1) + "\n"
    args.out.write_text(text)


def score(args: argparse.Namespace) -> None:
    directories = sorted(d for d in args.runs.iterdir() if (d / "breakdown.json").exists())
    runs = [json.loads((d / "breakdown.json").read_text()) for d in directories]
    share = args.t2_share if args.t2_share is not None else json.loads(T2_SHARE_FILE.read_text())["share"] if T2_SHARE_FILE.exists() else None
    if share is None:
        sys.exit(f"no T2 share: pass --t2-share or fix it in {T2_SHARE_FILE}")
    rows = [row for run in runs for sample in run["samples"] for row in sample["plates"]]
    rate = args.production_rate if args.production_rate is not None else production_rate(args.production_runs) if args.production_runs else None
    result = {
        "setting": runs[0]["setting"] if runs else None, "seeds": [run["seed"] for run in runs], "t2_share": share, "production_rate_per_world_100myr": rate,
        **force_targets(rows, share), **_world_speeds(runs), "closure": _closure(runs),
        "T8_contacts": contact_measures(directories), "T9_rifting": rifting_score(args.runs, rate),
    }
    text = json.dumps(result, indent=1, default=float) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text)
    print(text, end="")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "queue"):
        p = sub.add_parser(name)
        p.add_argument("--settings", type=Path, required=True)
        p.add_argument("--setting", required=True)
        p.add_argument("--out", type=Path, required=True)
        p.add_argument("--years", type=float, default=200.0, help="Myr")
        p.add_argument("--step-myr", type=float, default=1.0)
        p.add_argument("--every", type=float, default=10.0, help="Myr between breakdown samples")
        p.add_argument("--no-contacts", action="store_true")
        p.add_argument("--heldout", action="store_true", help="allow held-out seeds (G3 only)")
    sub.choices["run"].add_argument("--seed", type=int, required=True)
    sub.choices["queue"].add_argument("--seeds", type=int, nargs="+")
    sub.choices["queue"].add_argument("--seeds-file", type=Path)
    sub.choices["queue"].add_argument("--dry-run", action="store_true")
    p = sub.add_parser("score")
    p.add_argument("--runs", type=Path, required=True)
    p.add_argument("--production-runs", type=Path)
    p.add_argument("--production-rate", type=float)
    p.add_argument("--t2-share", type=float)
    p.add_argument("--out", type=Path)
    p = sub.add_parser("t2-census")
    p.add_argument("--seeds-file", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    {"run": run_one, "queue": queue, "score": score, "t2-census": t2_census}[args.command](args)


if __name__ == "__main__":
    main()
