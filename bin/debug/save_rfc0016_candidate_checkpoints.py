#!/usr/bin/env python3
"""Save reproducible source-world checkpoints for RFC 0016 window screens (#871).

Advance a catalog save at the registered 100 kyr schedule and write exact production
checkpoints after the requested number of steps. Candidate event at step k uses the
four outputs k-2 through k+1, so its checkpoint is after step k-3.

Run from backend/:

    .venv/bin/python ../bin/debug/save_rfc0016_candidate_checkpoints.py <save> \
        --after-steps 34 43 --step-years 100000 --out <captures>/checkpoints
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import persistence
from app import world as world_mod


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("save", type=Path)
    parser.add_argument("--after-steps", type=int, nargs="+", required=True)
    parser.add_argument("--step-years", type=float, default=100_000.0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    wanted = sorted(set(args.after_steps))
    if wanted[0] < 1:
        parser.error("--after-steps values must be positive")
    args.out.mkdir(parents=True, exist_ok=True)
    source = args.save.read_bytes()
    world = persistence.load_world_bytes(source)
    world.ice_age_period_years = 0.0
    records = []
    started = time.perf_counter()
    for step in range(1, wanted[-1] + 1):
        world_mod.step_world(world, args.step_years)
        if step in wanted:
            payload = persistence.save_world_bytes(world)
            checkpoint = args.out / f"after-step-{step:03d}.bloom"
            checkpoint.write_bytes(payload)
            record = {
                "after_step": step,
                "elapsed_years": world.elapsed_years,
                "sha256": hashlib.sha256(payload).hexdigest(),
                "bytes": len(payload),
            }
            records.append(record)
            print(json.dumps(record), flush=True)
    summary = {
        "source_save_sha256": hashlib.sha256(source).hexdigest(),
        "source_commit": __import__("subprocess").check_output(
            ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[2], text=True
        ).strip(),
        "step_years": args.step_years,
        "elapsed_s": time.perf_counter() - started,
        "checkpoints": records,
    }
    (args.out / "checkpoints.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
