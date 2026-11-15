"""The stage 1 queue of RFC 0018's D1 search (#888): every predeclared setting on every stage 1 world, one run at a time.

Usage (repo root, this worktree's venv; detach it, for example with ``subprocess.Popen(..., start_new_session=True)``):
    backend/.venv/bin/python bin/debug/run_stage1.py --out ~/mantle-bloom-captures/issue888/stage1 [--limit N] [--settings 0 3 ...]

Runs ``run_rift_damage_d1.py`` for each (setting, seed) pair of the predeclaration (``stage1_settings.py``, the five seeds below)
in setting order, skipping pairs whose ``run.json`` exists, so it resumes. Before each run it waits until the compute
budget says RUN (``bin/compute_budget.py check --run sim``) and until none of this task's other heavy runs is alive; it
never overrides a WAIT. ``--limit`` stops after that many runs (the first one is timed before the rest are released).
Pooled scoring is in ``score_stage1.py``.
"""

import argparse
import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE)]

import stage1_settings  # noqa: E402

REPO = HERE.parents[1]
STAGE1_SEEDS = [106533960, 1097624328, 1397817354, 1629738962, 591453549]
OWN_RUNS = ("run_production_baseline.py", "check_rift_lines_world_d1.py", "run_rift_damage_d1.py", "diagnose_line_placement.py")


def own_run_alive() -> bool:
    out = subprocess.run(["pgrep", "-f", "|".join(OWN_RUNS)], capture_output=True, text=True).stdout.split()
    return bool([pid for pid in out if int(pid) != os.getpid()])


def budget_says_run() -> bool:
    return subprocess.run([sys.executable, str(REPO / "bin/compute_budget.py"), "check", "--run", "sim"], capture_output=True).returncode == 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--settings", type=int, nargs="*")
    parser.add_argument("--years", type=float, default=200.0)
    args = parser.parse_args()
    table = stage1_settings.settings()

    done = 0
    for index, setting in enumerate(table):
        if args.settings is not None and index not in args.settings:
            continue
        for seed in STAGE1_SEEDS:
            target = args.out / f"setting{index:02d}" / str(seed)
            if (target / "run.json").exists():
                continue
            while own_run_alive() or not budget_says_run():
                time.sleep(60)
            target.mkdir(parents=True, exist_ok=True)
            started = time.perf_counter()
            command = [sys.executable, str(HERE / "run_rift_damage_d1.py"), "--seed", str(seed), "--out", str(target), "--years", str(args.years),
                       "--setting", json.dumps(setting)]
            subprocess.run(command, cwd=REPO / "backend", stdout=open(target / "run.log", "w"), stderr=subprocess.STDOUT, check=False)
            print(index, seed, round(time.perf_counter() - started), (target / "run.json").exists(), flush=True)
            done += 1
            if args.limit and done >= args.limit:
                return



if __name__ == "__main__":
    main()
