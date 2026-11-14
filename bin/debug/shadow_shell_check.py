"""Shadow thin-shell check of proxy C on a world (D1 reporting, #888; also the cherry-picked D0 proxy's smoke test).

Usage (repo root, this worktree's venv; one shell solve of about 1.3 s per plate):
    backend/.venv/bin/python bin/debug/shadow_shell_check.py --seeds 106533960 1097624328 [--steps 0] [--out stats.json]

For each world (generated from the seed and stepped plates-only ``--steps`` Myr), solves the thin shell on every
plate, scores the line-cut proxy against it with the D0 statistics at ``k_cut = 1.33`` (sign agreement, direction,
Spearman, top-decile ratios), and prints them. Plates with no connected thin-shell mesh are skipped and counted.
``shadow_stats(world)`` is what the D1 driver calls every 10 Myr on its stage-2 runs.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path[:0] = [str(Path(__file__).resolve().parent)]

import calibrate_rift_load_proxy as d0  # noqa: E402

K_CUT = 1.33


def shadow_stats(world) -> dict:
    records, skipped = [], 0
    for plate in world.plates:
        record = d0.compute_record(world, plate)
        skipped += record is None
        if record is not None:
            d0.prepare(record, "cut")
            records.append(record)
    if not records:
        return {"plates": 0, "skipped_no_mesh": skipped}
    result = d0.score(d0.pooled_metrics(records, "cut", (K_CUT,)))
    return {"plates": len(records), "skipped_no_mesh": skipped, **result, "passes": d0.passes(result)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seeds", type=int, nargs="+", required=True)
    parser.add_argument("--steps", type=int, default=0)
    parser.add_argument("--out", type=Path)
    args = parser.parse_args()
    from app import world as world_module

    out = {}
    for seed in args.seeds:
        world = world_module.generate_world(seed)
        world.simulate_climate_biomes = False
        for _ in range(args.steps):
            world_module.step_world(world, 1e6)
        out[str(seed)] = shadow_stats(world)
        print(seed, json.dumps(out[str(seed)]), flush=True)
    if args.out:
        args.out.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
