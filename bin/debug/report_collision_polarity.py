#!/usr/bin/env python3
"""Issue #318: replay diagnostics for continental collision polarity (collision_polarity.py).

Loads a save (or generates a world from `--seed`), steps it, and at every checkpoint records:

- which evidence decided each new front (`consumption`, `slab`, `arc`) and how often the
  labelled fallback had to (its rate, why: no evidence or ambiguous evidence, and which tier
  broke the tie: motion, size or plate id);
- ambiguity: contradictory physical evidence, opposing arcs, and arc cues the physical
  evidence overruled;
- record churn: fronts created, split, merged (and merges whose records disagreed), expired,
  dropped or fused by topology changes, plus the live record and evidence-bin counts;
- cost: prepass seconds per step against the whole step, the boundary searches the prepass
  ran, how many of deform()'s searches the shared cache answered, and the net extra
  searches per step.

A loaded save carries no history, so every front already in contact when it was written
decides by fallback; `--warmup-myr` steps past that before the counters are reset.

    cd backend
    .venv/bin/python ../bin/debug/report_collision_polarity.py \\
        ~/Downloads/mantle-bloom-seed997271774-87400000y.mbworld --myr 43.1 \\
        --out ../analysis/issue318/polarity-seed997271774.json
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import collision_polarity, persistence, world as world_mod  # noqa: E402


def _load(args) -> "world_mod.World":
    if args.save is not None:
        return persistence.load_world_bytes(args.save.read_bytes())
    return world_mod.generate_world(seed=args.seed)


def _checkpoint(world, step_seconds: float, steps: int) -> dict:
    report = collision_polarity.summary(world)
    report["elapsed_myr"] = world.elapsed_years / 1e6
    report["step_seconds_mean"] = step_seconds / steps if steps else 0.0
    report["prepass_share_of_step"] = report["prepass_seconds_per_step"] / report["step_seconds_mean"] if steps else 0.0
    return report


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("save", type=Path, nargs="?", help="a .mbworld save; omit to generate from --seed")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--myr", type=float, default=10.0)
    ap.add_argument("--warmup-myr", type=float, default=0.0, help="step this long first, then reset the counters")
    ap.add_argument("--step-years", type=float, default=100_000.0, help="UI default step")
    ap.add_argument("--checkpoint-myr", type=float, default=5.0)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    world = _load(args)
    for _ in range(int(round(args.warmup_myr * 1e6 / args.step_years))):
        world_mod.step_world(world, args.step_years)
    world.collision_polarity_stats = {}

    steps = int(round(args.myr * 1e6 / args.step_years))
    per_checkpoint = max(1, int(round(args.checkpoint_myr * 1e6 / args.step_years)))
    result = {"save": str(args.save), "seed": args.seed, "step_years": args.step_years, "checkpoints": []}
    step_seconds = 0.0
    for k in range(1, steps + 1):
        start = time.perf_counter()
        world_mod.step_world(world, args.step_years)
        step_seconds += time.perf_counter() - start
        if k % per_checkpoint == 0 or k == steps:
            report = _checkpoint(world, step_seconds, k)
            result["checkpoints"].append(report)
            print(
                f"{report['elapsed_myr']:7.1f} Myr  fronts {report['fronts_active']:3d}  "
                f"fallback {report['fallback_rate']:.0%}  ambiguous {report['ambiguity_rate']:.0%}  "
                f"decided {report['decision_share']}  churn {report['churn']}  "
                f"prepass {report['prepass_share_of_step']:.1%} of step, "
                f"deform reuse {report['deform_search_reuse']:.0%}",
                flush=True,
            )
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=1, default=float))


if __name__ == "__main__":
    main()
