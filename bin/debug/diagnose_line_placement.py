"""Why rift lines fail the D1 line check's validity test (#888, after its first result of 66% valid).

Usage (repo root, this worktree's venv; follow the compute-budget skill, about 90 s per world):
    backend/.venv/bin/python bin/debug/diagnose_line_placement.py --out ~/mantle-bloom-captures/issue888/line-diagnosis \
        [--seeds 106533960 ...] [--years 50]

For each development world stepped plates-only to ``--years`` Myr by the production splitter, takes every eligible
plate's rift line (no straining) and records: the sizes of the two largest components left when the line is removed
and of every stray piece; whether the RFC's vote rule for the line cells (each joins the side it touches most) leaves
both daughters connected; and the outcome of an alternative rule, a breadth-first assignment in which each line cell
joins the side of a neighbour that already has a side (connected by construction), at stray limits of 1%, 2%, 5%
and none. Nothing is changed in the model: this measures how much of the failure is the assignment rule, the stray
limit or the geometry of the line itself. A held-out or reserve seed is refused.
"""

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import connected_components

sys.path[:0] = [str(Path(__file__).resolve().parent)]

import check_rift_lines_world_d1 as line_check  # noqa: E402
import rift_damage_mechanism as mech  # noqa: E402
import rift_damage_step as st  # noqa: E402
import rift_line as rl  # noqa: E402
from run_rift_damage_d1 import forbidden_seeds  # noqa: E402


def bfs_sides(graph, chain, side):
    """Give every line cell the side of an already-labelled neighbour, in rounds, so each side stays connected."""
    pattern = (graph > 0).astype(np.int8).tocsr()
    open_cells = list(chain)
    while open_cells:
        assigned = {}
        for cell in open_cells:
            votes = rl.side_votes(pattern, side, cell)
            if votes.sum():
                assigned[cell] = int(np.argmax(votes))
        if not assigned:
            return False
        for cell, label in assigned.items():
            side[cell] = label
        open_cells = [c for c in open_cells if c not in assigned]
    return True


def alternative(graph, areas, chain, limit, p):
    """The breadth-first assignment with strays attached up to ``limit`` of the area; returns (valid, why, smaller share)."""
    n = graph.shape[0]
    total = float(areas.sum())
    keep = np.ones(n, bool)
    keep[chain] = False
    count, labels = connected_components(graph[keep][:, keep], directed=False)
    full = np.full(n, -1)
    full[keep] = labels
    if count < 2:
        return False, "no_separation", None
    sizes = np.bincount(labels, weights=areas[keep], minlength=count)
    order = np.argsort(-sizes)
    side = np.full(n, -1)
    for k, comp in enumerate(order[:2]):
        side[full == comp] = k
    if any(sizes[c] > limit * total for c in order[2:]):
        return False, "stray_too_large", None
    if not bfs_sides(graph, chain, side):
        return False, "unassigned", None
    pattern = (graph > 0).astype(np.int8).tocsr()
    for comp in order[2:]:
        members = np.flatnonzero(full == comp)
        side[members] = int(np.argmax(rl.side_votes(pattern, side, members)))
    mask = side == 0
    for part in (mask, ~mask):
        if connected_components(graph[part][:, part], directed=False)[0] != 1:
            return False, "daughter_disconnected", None
    small = float(min(areas[mask].sum(), areas[~mask].sum()) / total)
    return small >= p.a_min, None if small >= p.a_min else "daughter_small", small


def diagnose(world, p) -> list[dict]:
    from app import merge_split

    min_nodes = max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density))
    rows = []
    for plate in world.plates:
        if plate.node_count() < min_nodes:
            continue
        prep = st.prepare_plate(world, plate, p, 1e6)
        if prep.line is None:
            rows.append({"cells": plate.node_count(), "line": False, "corridors": prep.candidates})
            continue
        chain, graph, areas = prep.line.chain, prep.graph, prep.areas
        keep = np.ones(len(areas), bool)
        keep[chain] = False
        count, labels = connected_components(graph[keep][:, keep], directed=False)
        sizes = np.sort(np.bincount(labels, weights=areas[keep], minlength=count))[::-1] / areas.sum() if count else np.array([])
        mask, why = rl.split_masks(graph, areas, chain, p)
        row = {"cells": plate.node_count(), "line": True, "line_cells": len(chain), "components": int(count),
               "main_shares": [float(x) for x in sizes[:2]], "stray_shares": [float(x) for x in sizes[2:]],
               "rfc_rule": "valid" if mask is not None else why}
        for limit in (0.01, 0.02, 0.05, 1.0):
            valid, reason, small = alternative(graph, areas, chain, limit, p)
            row[f"bfs_stray{limit:g}"] = "valid" if valid else reason
        rows.append(row)
    return rows


def summarize(rows) -> dict:
    lines = [r for r in rows if r["line"]]
    out = {"plates": len(rows), "no_line": len(rows) - len(lines), "lines": len(lines)}
    for key in ("rfc_rule", "bfs_stray0.01", "bfs_stray0.02", "bfs_stray0.05", "bfs_stray1"):
        tally: dict = {}
        for r in lines:
            tally[r[key]] = tally.get(r[key], 0) + 1
        out[key] = {k: v for k, v in sorted(tally.items(), key=lambda kv: -kv[1])}
        out[key]["valid_share_of_lines"] = tally.get("valid", 0) / max(len(lines), 1)
    strays = [r["stray_shares"][0] if r["stray_shares"] else 0.0 for r in lines]
    out["largest_stray_share_quantiles"] = [float(x) for x in np.percentile(strays, [50, 75, 90, 99])]
    out["lines_with_any_stray"] = int(sum(bool(r["stray_shares"]) for r in lines))
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--seeds", type=int, nargs="+", default=line_check.DEV_SEEDS[:8])
    parser.add_argument("--years", type=float, default=50.0)
    args = parser.parse_args()
    if set(args.seeds) & forbidden_seeds():
        sys.exit("a held-out or reserve seed was given")
    args.out.mkdir(parents=True, exist_ok=True)
    mech.register_fields()
    p = replace(mech.Params())
    from app import world as world_module

    rows = []
    for seed in args.seeds:
        target = args.out / f"{seed}.json"
        if target.exists():
            rows += json.loads(target.read_text())
            continue
        world = world_module.generate_world(seed)
        world.simulate_climate_biomes = False
        for _ in range(round(args.years)):
            world_module.step_world(world, 1e6)
        got = diagnose(world, p)
        target.write_text(json.dumps(got) + "\n")
        rows += got
        print(seed, len(got), flush=True)
    summary = summarize(rows)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
