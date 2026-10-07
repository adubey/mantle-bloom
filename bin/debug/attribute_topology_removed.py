#!/usr/bin/env python3
"""Issue #305: what is still booked as `topology_removed_m3`, and why.

Steps generated worlds and classifies every continental material booking the topology
cleanup makes (merge_split.defragment_plates / remove_defunct_plates) by the piece dropped:
whether it carries continental crust, its node count, and its nearest distance (in node
spacings) to any live plate. Prints one JSON summary per seed.

    cd backend
    .venv/bin/python ../bin/debug/attribute_topology_removed.py --seeds 1 2 3 --steps 200
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    parser.add_argument("--seeds", type=int, nargs="+", required=True)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--years", type=float, default=1_000_000.0)
    args = parser.parse_args()

    sys.path.insert(0, str(args.backend.resolve()))
    import numpy as np
    from scipy.spatial import cKDTree

    from app import continental_ledger, merge_split
    from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad
    from app.world import generate_world, step_world

    for seed in args.seeds:
        world = generate_world(seed=seed)
        continental_ledger.ensure_initialized(world)
        start = continental_ledger.surface_volume_m3(world)
        spacing = line_spacing_rad(world.node_density)
        buckets: dict[str, float] = defaultdict(float)
        counts: dict[str, int] = defaultdict(int)
        current: dict[str, object] = {}

        original_accrete = merge_split._accrete_stranded_terrane

        def accrete(world_, fragment, radius, contacts):
            receiver = original_accrete(world_, fragment, radius, contacts)
            if receiver is None:
                current["fragment"] = fragment
            return receiver

        original_record = continental_ledger.record

        def record(world_, account, volume):
            if account == "topology_removed_m3" and volume > 0.0:
                fragment = current.pop("fragment", None)
                if fragment is None:
                    key = "unknown"
                else:
                    others = [p.all_points_and_elevation()[0] for p in world_.plates if p.plate_id != fragment.plate_id]
                    pts = fragment.all_points_and_elevation()[0]
                    gap = float(cKDTree(np.concatenate(others)).query(pts)[0].min()) / spacing if others else float("inf")
                    codes = fragment.collect("crust_type_code")
                    continental = effective_is_continental_from_codes(codes, fragment.crust_type == "continental")
                    kind = "continental" if np.any(continental) else "oceanic"
                    size = "1-3" if fragment.node_count() < 4 else ("4-49" if fragment.node_count() < 50 else "50+")
                    reach = "gap<=2.5" if gap <= 2.5 else ("gap<=5" if gap <= 5 else "gap>5")
                    key = f"{kind}/{size}/{reach}"
                buckets[key] += volume
                counts[key] += 1
            original_record(world_, account, volume)

        merge_split._accrete_stranded_terrane = accrete
        continental_ledger.record = record
        try:
            for _ in range(args.steps):
                step_world(world, args.years)
        finally:
            merge_split._accrete_stranded_terrane = original_accrete
            continental_ledger.record = original_record
        summary = {
            "seed": seed,
            "topology_removed_vs_start": world.continental_material_ledger["topology_removed_m3"] / start,
            "by_kind": {k: {"vs_start": v / start, "events": counts[k]} for k, v in sorted(buckets.items(), key=lambda kv: -kv[1])},
        }
        print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
