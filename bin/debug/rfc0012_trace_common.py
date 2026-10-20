"""Shared command-line setup for RFC 0012 replay trace tools."""

from __future__ import annotations

import argparse
from collections.abc import Callable
import json
from pathlib import Path


def run_trace_cli(description: str, run: Callable[..., None]) -> None:
    """Parse the common replay arguments and invoke a trace runner."""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--model-root", type=Path, default=Path.cwd())
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    if args.steps < 1 or args.step_myr <= 0:
        parser.error("steps and step-myr must be positive")
    run(args.seed, args.steps, args.step_myr, args.out, args.model_root.resolve())


def write_trace_summary(output: Path, summary: dict) -> None:
    """Append and print the final row shared by RFC 0012 trace replays."""
    with output.open("a") as handle:
        handle.write(json.dumps(summary, separators=(",", ":"), allow_nan=False) + "\n")
    print(json.dumps(summary, sort_keys=True), flush=True)
