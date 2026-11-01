#!/usr/bin/env python3
"""Save exact G8 lake-event checkpoints from a registered diagnosis world.

Replay the source save with the intended step schedule and save its exact state at each
requested model age. The ice-age control is disabled to match the RFC 0014 G6/G7 diagnosis
replays. This script does not perturb the world or evaluate a mechanism.

Example (run from the repository root):

    NUMBA_CACHE_DIR=/tmp/issue680-numba-cache backend/.venv/bin/python \\
      bin/debug/save_lake_g8_checkpoints.py <seed910-save.mbworld> \\
      --myr 1.7 --step-years 100000 --target-myr 221.3 --target-myr 222.1 \\
      --out /tmp/issue680-seed910-checkpoints.json

The output JSON records the source save checksum, source commit, target ages, checkpoint
checksums, and relative checkpoint filenames. Save files remain in the output file's directory.
"""

from __future__ import annotations

import hashlib
import json
import time

from replay_lib import load_replay, parse_replay_args


def add_arguments(parser) -> None:
    parser.add_argument(
        "--target-myr",
        action="append",
        type=float,
        required=True,
        help="model age at which to save a checkpoint; repeat for multiple checkpoints",
    )


args = parse_replay_args(__doc__, add_arguments)
world, result = load_replay(args, issue=680)
world.ice_age_period_years = 0.0

from app import persistence, world as world_mod  # noqa: E402

step_years = int(round(args.step_years))
target_years = sorted({int(round(myr * 1e6)) for myr in args.target_myr})
if len(target_years) != len(args.target_myr):
    raise SystemExit("target ages must be distinct")
if step_years <= 0 or any((target - world.elapsed_years) % step_years for target in target_years):
    raise SystemExit("every target age must align with the replay step schedule")
if any(target <= world.elapsed_years for target in target_years):
    raise SystemExit("every target age must be after the source save")
if any(target > world.elapsed_years + args.myr * 1e6 for target in target_years):
    raise SystemExit("every target age must fall within --myr")

args.out.parent.mkdir(parents=True, exist_ok=True)
checkpoints = []
next_target = 0
total_steps = max(1, int(round(args.myr * 1e6 / args.step_years)))
for step in range(1, total_steps + 1):
    started = time.perf_counter()
    world_mod.step_world(world, args.step_years)
    print(json.dumps({
        "step": step,
        "elapsed_years": world.elapsed_years,
        "runtime_s": time.perf_counter() - started,
    }), flush=True)
    if next_target < len(target_years) and round(world.elapsed_years) == target_years[next_target]:
        target_myr = target_years[next_target] / 1e6
        save_name = f"seed{world.seed}-{target_myr:g}myr.mbworld"
        checkpoint_path = args.out.parent / save_name
        if checkpoint_path.exists():
            raise SystemExit(f"refusing to overwrite existing checkpoint: {checkpoint_path}")
        payload = persistence.save_exact_bytes(world)
        checkpoint_path.write_bytes(payload)
        checkpoints.append({
            "elapsed_years": world.elapsed_years,
            "file": save_name,
            "sha256": hashlib.sha256(payload).hexdigest(),
            "bytes": len(payload),
        })
        manifest = {
            **{key: value for key, value in result.items() if key != "checkpoints"},
            "forcing": {"ice_age_period_years": 0.0},
            "step_years": args.step_years,
            "checkpoints": checkpoints,
        }
        args.out.write_text(json.dumps(manifest, indent=2) + "\n")
        next_target += 1

if next_target != len(target_years):
    raise SystemExit("replay ended before all target checkpoints were saved")
