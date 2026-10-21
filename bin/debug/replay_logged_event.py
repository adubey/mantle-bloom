#!/usr/bin/env python3
"""Replay one step a server's replay log recorded, and check it reproduces (#442).

`bin/restart.sh --replay-log` saves the world before any step whose health metrics go askew
(backend/app/replay_log.py). This loads a record's `before.mbworld`, runs the recorded step
(its controls are in the save), and compares the result with the post-step fingerprint,
plate count and craton ledger the server logged. It exits 1 when they don't match: a build
or dependency change since the record was written can do that, and then the record only
reproduces under its own commit (`git` in record.json).

A record of a step that raised replays as a match when the replay raises too.

Usage:
    backend/.venv/bin/python bin/debug/replay_logged_event.py \\
        ~/.mantle-bloom/replay-logs/8442/<record> [--save-after /tmp/after.mbworld]

To dig further, hand the pre-step save to a probe, for example for a craton clip:

    backend/.venv/bin/python bin/debug/probe_craton_clips.py <record>/before.mbworld \\
        --from-myr <before> --to-myr <after> --step-years <step_years>

The script prints that command with the record's numbers filled in.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import persistence, replay_log  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("record", type=Path, help="a record directory (holds before.mbworld and record.json)")
    parser.add_argument("--save-after", type=Path, help="also save the replayed post-step world here")
    args = parser.parse_args()

    result = replay_log.replay(args.record)
    sidecar = result.sidecar
    before_myr = sidecar["elapsed_years_before"] / 1e6
    after_myr = (sidecar["elapsed_years_before"] + sidecar["step_years"]) / 1e6
    print(f"record:   {args.record}")
    print(f"seed {sidecar['seed']}, {before_myr:.3f} -> {after_myr:.3f} Myr, git {sidecar.get('git')}")
    print("triggers:")
    for name, values in sidecar["triggers"].items():
        print(f"  {name}: {json.dumps(values, default=str)[:400]}")
    if result.error is not None:
        print(f"replay raised {type(result.error).__name__}: {result.error}")
    print(f"plate count matches: {result.plate_count_matches}")
    print(f"craton ledger matches: {result.craton_ledger_matches}")
    print(f"fingerprint mismatches: {', '.join(result.mismatched) or 'none'}")
    print("REPRODUCED" if result.matches else "NOT REPRODUCED")
    if args.save_after is not None:
        args.save_after.write_bytes(persistence.save_world_bytes(result.world))
        print(f"saved the replayed world to {args.save_after}")
    print("\nprobe the step:")
    print(
        f"  backend/.venv/bin/python bin/debug/probe_craton_clips.py {args.record / replay_log.BEFORE_FILE} "
        f"--from-myr {before_myr:g} --to-myr {after_myr:g} --step-years {sidecar['step_years']:g}"
    )
    return 0 if result.matches else 1


if __name__ == "__main__":
    sys.exit(main())
