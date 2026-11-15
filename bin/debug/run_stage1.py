"""The stage 1 queue of RFC 0018's D1 search (#888): every predeclared setting on every stage 1 world, one run at a time.

Usage (repo root, this worktree's venv; detach it, for example with ``subprocess.Popen(..., start_new_session=True)``):
    backend/.venv/bin/python bin/debug/run_stage1.py --out ~/mantle-bloom-captures/issue888/stage1 [--limit N] [--settings 0 3 ...]

Runs ``run_rift_damage_d1.py`` for each (setting, seed) pair of the predeclaration (``stage1_settings.py``, the five seeds below)
in setting order, skipping pairs whose ``run.json`` exists, so it resumes. With ``--workers N`` it keeps up to N heavy
runs of this task alive at once (the baseline runs count), starting a new one only while free-plus-inactive memory is
above ``--min-free-gb``. A run of this task uses about 0.4 GB, far below the 6 GB the compute-budget script reserves, so the
owner approved parallel runs on 2026-11-15 when no other agent is active; with other agents active, use the
compute-budget skill instead. ``--limit`` stops after that many runs.
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


def alive_runs() -> int:
    """How many heavy runs of this task are alive (this queue's workers count, its own process does not)."""
    out = subprocess.run(["pgrep", "-f", "|".join(OWN_RUNS)], capture_output=True, text=True).stdout.split()
    return len([pid for pid in out if int(pid) != os.getpid()])


def memory_ok(minimum_gb: float) -> bool:
    """Free plus inactive memory is at least ``minimum_gb`` (macOS ``vm_stat``)."""
    text = subprocess.run(["vm_stat"], capture_output=True, text=True).stdout
    pages = {line.split(":")[0].strip(): int(line.split(":")[1].strip().rstrip(".")) for line in text.splitlines() if ":" in line and "Pages" in line}
    return (pages.get("Pages free", 0) + pages.get("Pages inactive", 0)) * 16384 / 2**30 >= minimum_gb


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--settings", type=int, nargs="*")
    parser.add_argument("--years", type=float, default=200.0)
    parser.add_argument("--workers", type=int, default=1, help="heavy runs of this task alive at once (the owner approved 6 on 2026-11-15)")
    parser.add_argument("--min-free-gb", type=float, default=2.5)
    args = parser.parse_args()
    table = stage1_settings.settings()
    jobs = [(i, s, seed) for i, s in enumerate(table) for seed in STAGE1_SEEDS if args.settings is None or i in args.settings]
    jobs = [(i, s, seed) for i, s, seed in jobs if not (args.out / f"setting{i:02d}" / str(seed) / "run.json").exists()]
    if args.limit:
        jobs = jobs[:args.limit]
    running: list = []
    while jobs or running:
        for proc, started, label in list(running):
            if proc.poll() is not None:
                print(label, round(time.perf_counter() - started), proc.returncode, flush=True)
                running.remove((proc, started, label))
        if jobs and alive_runs() < args.workers and memory_ok(args.min_free_gb):
            index, setting, seed = jobs.pop(0)
            target = args.out / f"setting{index:02d}" / str(seed)
            target.mkdir(parents=True, exist_ok=True)
            command = [sys.executable, str(HERE / "run_rift_damage_d1.py"), "--seed", str(seed), "--out", str(target), "--years", str(args.years),
                       "--setting", json.dumps(setting)]
            proc = subprocess.Popen(command, cwd=REPO / "backend", stdout=open(target / "run.log", "w"), stderr=subprocess.STDOUT)
            running.append((proc, time.perf_counter(), f"setting{index:02d}/{seed}"))
            time.sleep(5)
            continue
        time.sleep(20)


if __name__ == "__main__":
    main()
