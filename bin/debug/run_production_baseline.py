"""One production-splitter run per development world: the rate the rift-line arm is compared with (RFC 0018 D1, #888).

Usage (repo root, this worktree's venv; follow the compute-budget skill, about 6 minutes per world):
    backend/.venv/bin/python bin/debug/run_production_baseline.py --out ~/mantle-bloom-captures/issue888/production-baseline \
        [--seeds 106533960 ...] [--years 200]

Generates each development world with production defaults, steps it plates-only (``simulate_climate_biomes =
False``) at 1 Myr for ``--years`` Myr, and records every split ``merge_split.maybe_split_plate`` applies: the time,
the plate and the smaller daughter's area share (a share under 7% is a tiny breakup). Writes ``<out>/<seed>.json`` per
world (existing files are kept, so a run resumes) and ``summary.json`` with the pooled rate in splits per world per
100 Myr, the figure the "splits happen" target compares the rift-line arm with. A held-out or reserve seed from
``analysis/issue797/g2-seed-preregistration.json`` is refused.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path[:0] = [str(Path(__file__).resolve().parent)]

from check_rift_lines_world_d1 import DEV_SEEDS  # noqa: E402
from run_rift_damage_d1 import forbidden_seeds  # noqa: E402

TINY_SHARE = 0.07


def run_world(seed: int, years: float, out: Path) -> dict:
    target = out / f"{seed}.json"
    if target.exists():
        return json.loads(target.read_text())
    from app import merge_split
    from app import world as world_module

    events = []
    original = merge_split.maybe_split_plate

    def counting(world, plate, step_years=0.0):
        result = original(world, plate, step_years)
        if result is not None:
            a, b = result
            areas = [float(d.node_areas_m2().sum()) for d in (a, b)]
            events.append({"myr": world.elapsed_years / 1e6, "plate": int(plate.plate_id), "small_share": min(areas) / sum(areas)})
        return result

    merge_split.maybe_split_plate = counting
    try:
        started = time.perf_counter()
        world = world_module.generate_world(seed)
        world.simulate_climate_biomes = False
        counts = []
        for step in range(round(years)):
            world_module.step_world(world, 1e6)
            counts.append(len(world.plates))
        result = {"seed": seed, "years": years, "splits": events, "n_splits": len(events), "tiny": sum(e["small_share"] < TINY_SHARE for e in events),
                  "plates_final": len(world.plates), "plates_min_max": [int(min(counts)), int(max(counts))], "wall_s": time.perf_counter() - started}
    finally:
        merge_split.maybe_split_plate = original
    target.write_text(json.dumps(result, indent=1) + "\n")
    print(seed, result["n_splits"], round(result["wall_s"]), flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=DEV_SEEDS)
    parser.add_argument("--years", type=float, default=200.0)
    args = parser.parse_args()
    if set(args.seeds) & forbidden_seeds():
        sys.exit("a held-out or reserve seed was given")
    args.out.mkdir(parents=True, exist_ok=True)
    results = [run_world(seed, args.years, args.out) for seed in args.seeds]
    splits = sum(r["n_splits"] for r in results)
    summary = {"worlds": len(results), "years": args.years, "splits": splits, "splits_per_world_100myr": splits / len(results) / args.years * 100.0,
               "per_world": {str(r["seed"]): r["n_splits"] for r in results}, "tiny_breakups": sum(r["tiny"] for r in results),
               "small_share_median": float(np.median([e["small_share"] for r in results for e in r["splits"]])) if splits else None}
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
