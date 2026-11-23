"""Scaffolding the #998 sampling scripts share: `plate_speed_samples.py`, `downwelling_check.py` and
`contact_flow_coherence.py`.

Each script steps one seeded world plates-only and records something at chosen times. This module holds what they have
in common: building the world from a setting, the boundary-band context of a world, stepping to the sample times, running
seeds in series as separate processes, and the `run` / `queue` command line. Not a script itself.
"""

import argparse
import copy
import subprocess
import sys
import time
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "backend")]


def build_world(args: argparse.Namespace) -> tuple[dict, object]:
    """The setting's world for `args.seed` with its module knobs applied; refuses a held-out seed."""
    import calibrate_plate_forces as cal
    import run_rift_lines_d1 as d1

    if args.seed in cal.refused_seeds():
        sys.exit(f"seed {args.seed} is held out")
    setting = cal.load_setting(args.settings, args.setting)
    world = d1.make_world(SimpleNamespace(seed=args.seed, node_density=1.0, world_knobs=setting["world_knobs"],
                                          split_method=setting["split_method"], params=setting["params"]))
    cal.apply_module_knobs(setting["module_knobs"])
    return setting, world


def reference_centers(world, scale: float) -> list:
    """Copies of the world's generated convection centres with strength x `scale`, taken before any flow intervention."""
    centers = copy.deepcopy(world.mantle_centers)
    for center in centers:
        center.strength *= scale
    return centers


def boundary_context(world) -> SimpleNamespace:
    """Spacing, force reach, the plates with nodes (by id), the RFC 0021 force forms and the boundary search cache
    (built here when no solve has made it yet, as at t=0)."""
    from app import torque
    from app.elevation_lines import line_spacing_rad

    spacing = line_spacing_rad(world.node_density)
    cache = world.boundary_search_cache
    if cache is None:
        cache = torque.BoundarySearchCache()
        world.boundary_search_cache = cache
    plates = sorted((p for p in world.plates if p.node_count() > 0), key=lambda p: p.plate_id)
    return SimpleNamespace(spacing=spacing, reach=torque.boundary_force_reach_rad(spacing), plates=plates, forms=torque.force_forms(world), cache=cache)


def plate_boundary(world, ctx: SimpleNamespace, plate):
    """(neighbours, inputs, band) for one plate: its boundary inputs with the force forms attached, and the mask of nodes
    inside the contact band that have a neighbour."""
    from app import torque

    neighbours = ctx.cache.neighbours(plate, ctx.plates, ctx.reach)
    inputs = torque.gather_boundary_force_inputs(plate, neighbours, ctx.spacing, ctx.reach, ctx.cache)
    torque.apply_force_forms(world, plate, neighbours, inputs, ctx.forms)
    band = (inputs.dist_to_neighbor <= ctx.reach) & (inputs.neighbor_plate_id >= 0)
    return neighbours, inputs, band


def step_and_sample(world, setting: dict, args: argparse.Namespace, take) -> None:
    """Call `take(myr)` at t=0 (before any step, when 0 is in `--at`) and after each step that lands on a time in `--at`."""
    import run_rift_lines_d1 as d1

    wanted = sorted(set(args.at))
    if wanted and wanted[0] == 0:
        take(0.0)
    final = max(wanted)
    scratch = args.out.parent / f"{args.out.stem}-scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    stepping = SimpleNamespace(years=final, step_myr=args.step_myr, max_steps=None, snapshot_every=0.0)
    recorder = None
    if setting["split_method"] == "rift_lines":
        import rift_line_recorder as rec
        from app import rift_damage

        recorder = rec.Recorder(rift_damage.runtime(world).params.eps_cand)
    if final > 0:
        d1.step_all(world, stepping, scratch, recorder, lambda w, myr: take(myr) if round(myr, 6) in wanted else None)


def queue_seeds(args: argparse.Namespace, script: str, forward: tuple[str, ...]) -> None:
    """Run the `run` command of `script` for each seed in series, one process each, skipping finished ones. `forward` names
    the parsed arguments (underscored) passed on as `--flag value`."""
    import json

    seeds = args.seeds or json.loads(args.seeds_file.read_text())["seeds"]
    args.out.mkdir(parents=True, exist_ok=True)
    failed = []
    for seed in seeds:
        target = args.out / f"{seed}.json"
        if target.exists():
            continue
        command = [sys.executable, script, "run", "--settings", str(args.settings), "--setting", args.setting, "--seed", str(seed),
                   "--step-myr", str(args.step_myr), "--out", str(target)]
        for name in forward:
            command += [f"--{name.replace('_', '-')}", str(getattr(args, name))]
        command += ["--at", *map(str, args.at)]
        with open(args.out / f"{seed}.log", "w") as log:
            code = subprocess.call(command, stdout=log, stderr=subprocess.STDOUT)
        print(seed, "exit", code, time.strftime("%H:%M"), flush=True)
        if code:
            failed.append(seed)
    sys.exit(1 if failed else 0)


def build_parser(doc: str, script: str, run_one, *, forward: tuple[str, ...], extra=lambda p: None, default_at=(0, 5, 20, 50), aggregate=None, reference_scale=False) -> argparse.ArgumentParser:
    """The `run` and `queue` commands with the arguments every sampler takes; `extra(parser)` adds a script's own flags
    (name them in `forward` so `queue` passes them on). `reference_scale` adds `--reference-scale` (the scale of the original
    centres that downwelling potential is measured against). `aggregate(args)` adds an `aggregate --runs DIR --out FILE` command."""
    ap = argparse.ArgumentParser(description=doc, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "queue"):
        p = sub.add_parser(name)
        p.add_argument("--settings", type=Path, required=True)
        p.add_argument("--setting", required=True)
        p.add_argument("--at", type=float, nargs="+", default=list(default_at))
        p.add_argument("--step-myr", type=float, default=1.0)
        p.add_argument("--flow-scale", type=float, default=1.0, help="multiplier on every mantle convection centre's strength (default 1; 0 = off, negative flips the sign)")
        p.add_argument("--out", type=Path, required=True)
        if reference_scale:
            p.add_argument("--reference-scale", type=float, default=4.0, help="scale of the original centres that D is measured against")
        extra(p)
        if name == "run":
            p.add_argument("--seed", type=int, required=True)
            p.set_defaults(fn=run_one)
        else:
            p.add_argument("--seeds", type=int, nargs="*")
            p.add_argument("--seeds-file", type=Path, default=Path("analysis/issue958/development-seeds.json"))
            p.set_defaults(fn=lambda a: queue_seeds(a, script, forward))
    if aggregate is not None:
        agg = sub.add_parser("aggregate")
        agg.add_argument("--runs", type=Path, required=True)
        agg.add_argument("--out", type=Path, required=True)
        agg.set_defaults(fn=aggregate)
    return ap


def scale_flow(world, factor: float) -> None:
    """Scale every convection centre's strength after generation, so every seed starts from the same world and only the
    stepped flow changes (exploration of #998 only)."""
    for center in world.mantle_centers:
        center.strength *= factor


__all__ = ["boundary_context", "build_parser", "build_world", "plate_boundary", "queue_seeds", "reference_centers", "scale_flow", "step_and_sample"]
