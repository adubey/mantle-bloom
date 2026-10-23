"""Shared CLI and output helpers for G2 replay diagnostics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import time


def add_replay_arguments(parser: argparse.ArgumentParser, *, default_steps: int) -> None:
    parser.add_argument("--steps", type=int, default=default_steps)
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--snapshot-every", type=int, default=10)
    parser.add_argument("--out", type=Path, required=True)


def prepare_replay_output(parser: argparse.ArgumentParser, args: argparse.Namespace) -> None:
    if args.steps < 1 or args.step_myr <= 0 or args.snapshot_every < 1:
        parser.error("steps, step-myr, and snapshot-every must be positive")
    args.out.parent.mkdir(parents=True, exist_ok=True)


def write_replay_summary(
    out: Path,
    *,
    seed: int,
    steps: int,
    step_myr: float,
    rows: list[dict],
    started: float,
    ledger: dict,
    **experiment_fields,
) -> None:
    summary = {
        "done": True,
        "seed": seed,
        "steps": steps,
        "step_myr": step_myr,
        "elapsed_myr": steps * step_myr,
        "runtime_min": (time.perf_counter() - started) / 60,
        "checkpoints": len(rows),
        "ledger": ledger,
        **experiment_fields,
    }
    out.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True), flush=True)
