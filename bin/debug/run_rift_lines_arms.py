"""Run every (arm, seed) pair of an arms file with `run_rift_lines_d1.py`, a few at a time (RFC 0018 D1, #888).

Usage (repo root, this worktree's venv; follow the compute-budget skill first, and count the whole queue as one heavy
run with ``--jobs`` workers):
    backend/.venv/bin/python bin/debug/run_rift_lines_arms.py --arms analysis/issue888/rift-lines-baseline/arms.json \
        --out ~/mantle-bloom-captures/issue888/rift-lines-baseline [--jobs 3] [--only base deep-marks] [--dry-run]

The arms file holds ``seeds``, ``years`` and ``arms``: a map from arm name to ``split_method`` (default
``rift_lines``), ``params`` (`World.rift_line_params`) and ``world_knobs``, with ``common_params`` and
``common_world_knobs`` applied under every rift-line arm's own and every arm's own, respectively. Each run writes
``<out>/<arm>/<seed>/``, with its log as ``run.log``; a pair whose ``run.json`` exists is skipped, so the queue resumes.
``<out>/queue.log`` gets one line per finished run. Score with `score_rift_lines.py`.
"""

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent


def jobs(spec: dict, only: list[str] | None) -> list[tuple[str, int, list[str]]]:
    """(arm, seed, command) for every pair of the arms file, in file order."""
    out = []
    for name, arm in spec["arms"].items():
        if only and name not in only:
            continue
        method = arm.get("split_method", "rift_lines")
        params = {**spec.get("common_params", {}), **arm.get("params", {})} if method == "rift_lines" else {}
        knobs = {**spec.get("common_world_knobs", {}), **arm.get("world_knobs", {})}
        for seed in spec["seeds"]:
            out.append((name, seed, [sys.executable, str(HERE / "run_rift_lines_d1.py"), "--seed", str(seed), "--split-method", method,
                                     "--params", json.dumps(params), "--world-knobs", json.dumps(knobs), "--years", str(spec.get("years", 200.0))]))
    return out


def run_one(out: Path, name: str, seed: int, command: list[str]) -> int:
    target = out / name / str(seed)
    if (target / "run.json").exists():
        return 0
    target.mkdir(parents=True, exist_ok=True)
    with open(target / "run.log", "w") as log:
        code = subprocess.call([*command, "--out", str(target)], stdout=log, stderr=subprocess.STDOUT)
    with open(out / "queue.log", "a") as queue:
        queue.write(f"{name}/{seed} exit {code} {time.strftime('%H:%M')}\n")
    return code


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--arms", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--only", nargs="+", help="arm names to run")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    queue = jobs(json.loads(args.arms.read_text()), args.only)
    if args.dry_run:
        for name, seed, command in queue:
            print(name, seed, " ".join(command[1:]))
        return
    args.out.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(args.jobs) as pool:
        codes = list(pool.map(lambda job: run_one(args.out, *job), queue))
    failed = [f"{n}/{s}" for (n, s, _), c in zip(queue, codes) if c]
    print(json.dumps({"runs": len(queue), "failed": failed}))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
