#!/usr/bin/env python3
"""Issue #318: what are the front merges whose records disagree on polarity? Replays a save
and, whenever front matching absorbs a record whose lower plate differs from the record that
survives, records both: how each was decided (source, scope, ambiguous, fallback basis), age,
contact steps, split lineage, stored size, and the share of the merged front each covered.

    cd backend
    .venv/bin/python ../bin/debug/diagnose_polarity_merge_conflicts.py \\
        ~/Downloads/mantle-bloom-seed997271774-87400000y.mbworld --myr 43.1 \\
        --out ../analysis/issue318/merge-conflicts-seed997271774.json
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import collision_polarity as cp, geometry, persistence, world as world_mod  # noqa: E402


def _describe(record: cp.CollisionFront, now: float) -> dict:
    return {
        "front_id": record.front_id,
        "lower": record.lower_plate_id,
        "source": record.source,
        "scope": getattr(record, "scope", None),
        "ambiguous": record.ambiguous,
        "fallback_basis": record.fallback_basis,
        "age_myr": (now - record.established_years) / 1e6,
        "since_contact_myr": (now - record.last_contact_years) / 1e6,
        "contact_steps": record.contact_steps,
        "parent_id": record.parent_id,
        "stored_points": {str(k): int(len(v)) for k, v in record.side_points.items()},
        "votes": record.votes,
    }


def _coverage(record: cp.CollisionFront, front: dict, contacts: dict, tol: float) -> float:
    """Share of the front's nodes within `tol` of the record's stored nodes."""
    hit = total = 0
    for side, idx in front.items():
        stored = record.side_points.get(side)
        total += len(idx)
        if stored is None or not len(stored) or not len(idx):
            continue
        stored_world = geometry.to_world(contacts[side].plate.frame, stored)
        hit += int(np.isfinite(cKDTree(stored_world).query(contacts[side].points[idx], distance_upper_bound=tol)[0]).sum())
    return hit / total if total else 0.0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("save", type=Path)
    ap.add_argument("--myr", type=float, default=43.1)
    ap.add_argument("--step-years", type=float, default=100_000.0)
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    world = persistence.load_world_bytes(args.save.read_bytes())
    conflicts: list[dict] = []
    merges = {"total": 0}
    match = cp._match_fronts

    def match_and_record(w, pair, fronts, contacts, spacing_rad, now):
        before = {r.front_id: copy.deepcopy(r) for r in w.collision_fronts if r.plate_ids == pair}
        out = match(w, pair, fronts, contacts, spacing_rad, now)
        after = {r.front_id for r in w.collision_fronts if r.plate_ids == pair}
        absorbed = [r for fid, r in before.items() if fid not in after]
        merges["total"] += len(absorbed)
        tol = cp.FRONT_MATCH_SPACINGS * spacing_rad
        for gone in absorbed:
            # The front it folded into: the one its stored nodes cover most.
            covers = [_coverage(gone, f, contacts, tol) for f in fronts]
            k = int(np.argmax(covers))
            survivor = out[k]
            if survivor.lower_plate_id == gone.lower_plate_id:
                continue
            prior = before.get(survivor.front_id)
            conflicts.append(
                {
                    "elapsed_myr": now / 1e6,
                    "pair": list(pair),
                    "merged_front_nodes": {str(s): int(len(i)) for s, i in fronts[k].items()},
                    "survivor": _describe(prior if prior is not None else survivor, now),
                    "survivor_coverage": _coverage(prior if prior is not None else survivor, fronts[k], contacts, tol),
                    "absorbed": _describe(gone, now),
                    "absorbed_coverage": covers[k],
                    "fronts_this_step": len(fronts),
                }
            )
        return out

    cp._match_fronts = match_and_record
    steps = int(round(args.myr * 1e6 / args.step_years))
    for k in range(steps):
        world_mod.step_world(world, args.step_years)
        if (k + 1) % 50 == 0:
            print(f"{world.elapsed_years / 1e6:.1f} Myr: {merges['total']} merges, {len(conflicts)} conflicts", flush=True)
    if args.out is not None:
        args.out.write_text(json.dumps({"save": str(args.save), "merges": merges["total"], "conflicts": conflicts}, indent=1, default=float))


if __name__ == "__main__":
    main()
