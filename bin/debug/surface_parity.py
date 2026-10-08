#!/usr/bin/env python3
"""Issues #247/#249: plate-surface audit runs and their quality gates.

    cd backend
    .venv/bin/python ../bin/debug/surface_parity.py audit --preset smoke --seeds 3 --out ../analysis/parity-smoke

`audit` runs every seed -- in parallel with `--jobs` -- then judges them: `comparison.json`
holds every gate result and `report.md` a readable summary. The exit status is 1 when the
verdict is `fail`, so a script can gate on it. `run` does one seed's world, optionally
continuing a saved one (`--from-world`); `judge` re-judges an existing output directory, e.g.
after combining runs made on several machines.

Presets (`surface_parity.PRESETS`): `smoke` (CI-sized, density 0.5, 4 Myr), `standard`
(density 1, 120 Myr), `long` (density 1, 400 Myr), `issue147` (the issue #147 profile world:
density 4, climate 4, fluid 2, 100 kyr steps to 60 Myr). Any preset field can be overridden.
See docs/surface-parity.md.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import platform
import shlex
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

import numpy as np  # noqa: E402
import scipy  # noqa: E402

from app import surface_parity, surface_parity_gates  # noqa: E402


def config_from_args(args: argparse.Namespace) -> surface_parity.RunConfig:
    config = surface_parity.PRESETS[args.preset]
    overrides = {}
    if args.checkpoints_myr:
        overrides["checkpoints_myr"] = tuple(float(v) for v in args.checkpoints_myr.split(","))
    for name in ("node_density", "step_years", "audit_every", "samples", "climate_density", "fluid_density"):
        value = getattr(args, name)
        if value is not None:
            overrides[name] = value
    if args.render:
        overrides["render"] = True
    if args.no_load_checks:
        overrides["load_checks"] = False
    if overrides:
        overrides.setdefault("name", config.name)
        config = dataclasses.replace(config, **overrides)
    config.checkpoint_steps  # validates checkpoint ages against the step length
    return config


def provenance(argv: list[str]) -> dict:
    def git(*command: str) -> str:
        return subprocess.run(["git", *command], cwd=REPO, capture_output=True, text=True).stdout.strip()

    return {
        "command": " ".join(shlex.quote(a) for a in ["bin/debug/surface_parity.py", *argv]),
        "git_commit": git("rev-parse", "HEAD"),
        "git_dirty": bool(git("status", "--porcelain", "--", "backend/app", "bin/debug")),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "machine": platform.machine(),
        "platform": platform.platform(),
    }


def _run_one(config: surface_parity.RunConfig, seed: int, out: Path, initial_world: Path | None = None) -> str:
    surface_parity.run_world(config, seed, out, log=lambda message: print(message, flush=True), initial_world=initial_world)
    return surface_parity.run_name(seed)


def judge(out: Path, commands: list[str] | None = None) -> int:
    runs = surface_parity.load_results(out)
    if not runs:
        print(f"no runs under {out}", file=sys.stderr)
        return 2
    evaluation = surface_parity_gates.evaluate(runs)
    if commands is None:
        meta = out / "provenance.json"
        commands = [json.loads(meta.read_text())["command"]] if meta.exists() else []
    (out / "comparison.json").write_text(surface_parity.dumps(evaluation))
    (out / "report.md").write_text(surface_parity_gates.render_report(evaluation, runs, commands))
    counts = {}
    for gate in evaluation["gates"].values():
        counts[gate["status"]] = counts.get(gate["status"], 0) + 1
    print(f"verdict {evaluation['verdict']}: {counts} -- {out / 'report.md'}")
    for key, gate in evaluation["gates"].items():
        if gate["status"] in (surface_parity_gates.FAIL, surface_parity_gates.WARN):
            print(f"  {gate['status']:5} {key}: {gate['worst']}")
    return 1 if evaluation["verdict"] == surface_parity_gates.FAIL else 0


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    for name in ("run", "audit"):
        p = sub.add_parser(name)
        p.add_argument("--preset", choices=sorted(surface_parity.PRESETS), default="smoke")
        p.add_argument("--out", type=Path, required=True)
        p.add_argument("--checkpoints-myr", help="comma-separated ages; overrides the preset's")
        p.add_argument("--node-density", type=float)
        p.add_argument("--step-years", type=float)
        p.add_argument("--audit-every", type=int)
        p.add_argument("--samples", type=int)
        p.add_argument("--climate-density", type=float)
        p.add_argument("--fluid-density", type=float)
        p.add_argument("--render", action="store_true", help="write elevation/platesDetail PNGs at each checkpoint")
        p.add_argument("--no-load-checks", action="store_true")
        if name == "run":
            p.add_argument("--seed", type=int, required=True)
            p.add_argument("--from-world", type=Path, help="continue this .mbworld instead of generating (its seed must match)")
        else:
            p.add_argument("--seeds", default="3", help="comma-separated seeds")
            p.add_argument("--jobs", type=int, default=1, help="parallel worlds; timings are only comparable at --jobs 1")
    p = sub.add_parser("judge")
    p.add_argument("out", type=Path)
    args = parser.parse_args(argv)

    if args.command == "judge":
        return judge(args.out)

    config = config_from_args(args)
    args.out.mkdir(parents=True, exist_ok=True)
    meta = provenance(argv) | {"config": config.to_json()}
    if args.command == "run":
        (args.out / f"provenance-{surface_parity.run_name(args.seed)}.json").write_text(json.dumps(meta, indent=2) + "\n")
        _run_one(config, args.seed, args.out, args.from_world)
        return 0

    meta["jobs"] = args.jobs
    (args.out / "provenance.json").write_text(json.dumps(meta, indent=2) + "\n")
    seeds = [int(seed) for seed in args.seeds.split(",")]
    if args.jobs > 1:
        with ProcessPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(_run_one, config, seed, args.out) for seed in seeds]
            for future in as_completed(futures):
                print("finished", future.result(), flush=True)
    else:
        for seed in seeds:
            _run_one(config, seed, args.out)
    return judge(args.out, [meta["command"]])


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
