"""Which plate-level features predict "no structurally valid candidate" before any eigensolve? (#848)

Usage (repo root, this worktree's venv; follow the compute-budget skill first):
    backend/.venv/bin/python bin/debug/study_rfc0018_skip_predictors.py collect \
        --save <world.mbworld> --g0-script <measure_rfc0018_g0.py> --out rows.jsonl \
        [--warmup-steps 3] [--steps 10]
    backend/.venv/bin/python bin/debug/study_rfc0018_skip_predictors.py analyze rows.jsonl [rows2.jsonl ...] \
        [--out analysis/issue848/skip-predictors.json]

``collect`` loads a checkpoint, steps it plates-only with the production splitter, and at each measured
step runs the rebuilt mass-balanced generator (``rebuild_rfc0018_mass_balanced.propose``) on every split
opportunity, recording the candidate outcome next to features computed BEFORE any eigensolve: plate size and
age, internal stress, flow-fit residual against its threshold, the balance, area share, spatial contiguity and
cut length of the two flow clusters, and steps since the plate last ran ``quad_tectonics._open_rift``.

Under the revised event model (RFC 0018 revision, #857) a plate with no structurally valid candidate has
``q0 = 1``: no path and no geological mark. Skipping its search is therefore exact. ``analyze`` reports, per
feature, the best threshold that skips opportunities with zero false negatives (no skipped opportunity had a
valid candidate), fit on all-but-one world and scored on the held-out one. A zero-false-negative threshold on
a sample is evidence, not a proof; only a bound derived from the model could make a skip exact.
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

import rebuild_rfc0018_mass_balanced as rebuild  # noqa: E402
from scipy.cluster.vq import kmeans2  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import connected_components  # noqa: E402

from app import mantle, merge_split, persistence, quad_tectonics  # noqa: E402
from app import world as world_mod  # noqa: E402

FEATURES = ["nodes", "age_myr", "internal_stress", "residual_ratio", "small_cluster_node_share",
            "small_cluster_area_share", "cluster_components", "cluster_edge_share", "cluster_cut_per_node",
            "steps_since_open_rift"]


def cheap_features(g0, world, plate, since_open_rift) -> dict:
    """Plate-level features computed without any Fiedler solve."""
    points, _ = plate.all_points_and_elevation()
    flow = mantle.flow_at(points, world.mantle_centers)
    _, labels = kmeans2(flow, k=2, minit="++", seed=plate.plate_id)
    i, j, lengths = g0._shared_edges(plate)
    n = len(points)
    areas = plate.node_areas_m2()
    sizes = np.array([(labels == k).sum() for k in (0, 1)])
    small = int(np.argmin(sizes))
    unlike = labels[i] != labels[j]
    graph = coo_matrix((np.ones(int((~unlike).sum())), (i[~unlike], j[~unlike])), shape=(n, n))
    components = int(connected_components(graph, directed=False)[0])  # parts of the like-labelled graph
    _, radius = merge_split.geometry.bounding_sphere(points)
    relax = max(min(1.0, radius / merge_split.SPLIT_SIZE_CERTAIN_RIFT_RAD),
                min(1.0, plate.internal_stress / merge_split.STRESS_CERTAIN_RIFT))
    threshold = merge_split.SPLIT_RMS_RESIDUAL_THRESHOLD * (1.0 - relax)
    residual = merge_split._fit_residual_rms(points, flow, plate.omega)
    return {
        "nodes": int(n), "age_myr": float(plate.age_years / 1e6), "internal_stress": float(plate.internal_stress),
        "residual_ratio": float(residual / threshold) if threshold > 0 else float("inf"),
        "small_cluster_node_share": float(sizes[small] / n),
        "small_cluster_area_share": float(areas[labels == small].sum() / areas.sum()),
        "cluster_components": components,
        "cluster_edge_share": float(unlike.mean()),
        "cluster_cut_per_node": float(lengths[unlike].sum() / max(median_length(lengths), 1.0) / n),
        "steps_since_open_rift": since_open_rift,
    }


def median_length(lengths) -> float:
    return float(np.median(lengths))


def collect(args) -> None:
    g0 = rebuild.load_g0(args.g0_script.expanduser())
    with open(args.save.expanduser(), "rb") as handle:
        world = persistence.load_world_file(handle)
    world.simulate_climate_biomes = False
    state = {"step": 0, "diagnosed": None}
    last_open_rift = {}
    original_open_rift, original_split = quad_tectonics._open_rift, merge_split.maybe_split_plate
    out = open(args.out.expanduser(), "w")

    def open_rift_spy(plate, *a, **kw):
        last_open_rift[plate.plate_id] = state["step"]
        return original_open_rift(plate, *a, **kw)

    def split_spy(current_world, plate, step_years=0.0):
        if state["step"] > args.warmup_steps and state["diagnosed"] != state["step"]:
            state["diagnosed"] = state["step"]
            for probe in list(current_world.plates):
                if not g0._is_opportunity(current_world, probe):
                    continue
                since = last_open_rift.get(probe.plate_id)
                row = {"save": args.save.name, "seed": int(current_world.seed), "step": state["step"] - args.warmup_steps,
                       "elapsed_myr": current_world.elapsed_years / 1e6, "plate_id": int(probe.plate_id)}
                row.update(cheap_features(g0, current_world, probe, None if since is None else state["step"] - since))
                started = time.perf_counter()
                try:
                    candidates = rebuild.propose(g0, current_world, probe)
                    prior = rebuild.breakup_prior(candidates)
                    row.update({"status": "ok", "candidates": len(candidates), "p_break": prior["P_break"],
                                "valid": sum(c["structurally_valid"] for c in candidates)})
                except ValueError as exc:  # the generator raises when no candidate can be proposed
                    row.update({"status": str(exc), "candidates": 0, "valid": 0, "p_break": 0.0})
                row["search_s"] = time.perf_counter() - started
                out.write(json.dumps(row) + "\n")
                out.flush()
        return original_split(current_world, plate, step_years)

    quad_tectonics._open_rift, merge_split.maybe_split_plate = open_rift_spy, split_spy
    try:
        for _ in range(args.warmup_steps + args.steps):
            state["step"] += 1
            world_mod.step_world(world, 1e6)
            print(json.dumps({"step": state["step"], "plates": len(world.plates)}), flush=True)
    finally:
        quad_tectonics._open_rift, merge_split.maybe_split_plate = original_open_rift, original_split
        out.close()


def best_zero_fn_threshold(values: np.ndarray, valid: np.ndarray):
    """The rule 'skip if feature <= t' or '>= t' that skips the most opportunities with no valid one skipped."""
    best = {"skipped": 0}
    for direction in ("le", "ge"):
        positives = values[valid]
        if not len(positives):
            continue
        # skip only values strictly beyond every opportunity that has a valid candidate
        t = positives.min() if direction == "le" else positives.max()
        mask = values < t if direction == "le" else values > t
        if mask.sum() > best["skipped"]:
            best = {"direction": direction, "threshold": float(t), "skipped": int(mask.sum())}
    return best


def analyze(args) -> None:
    rows = []
    for index, path in enumerate(args.rows):
        for line in Path(path).read_text().splitlines():
            if line.strip():
                row = json.loads(line)
                # Rows written before the seed field existed: the files come in (50, 100 Myr) pairs per world.
                row["world"] = str(row.get("seed", f"files{index // 2}"))
                rows.append(row)
    worlds = sorted({r["world"] for r in rows})
    valid = np.array([r["valid"] > 0 for r in rows])
    cost = np.array([r["search_s"] for r in rows])
    p_break = np.array([r["p_break"] for r in rows])
    summary = {"opportunities": len(rows), "worlds": worlds, "with_valid": int(valid.sum()),
               "share_without_valid": float(1 - valid.mean()), "search_s_total": float(cost.sum()),
               "search_s_without_valid": float(cost[~valid].sum()), "p_break_mean": float(np.mean([r["p_break"] for r in rows])),
               "features": {}}
    for name in FEATURES:
        raw = [r[name] for r in rows]
        values = np.array([np.nan if v is None else v for v in raw], float)
        if np.isnan(values).any():  # a missing value (never stretched) counts as very large
            values = np.where(np.isnan(values), 1e9, values)
        fit = best_zero_fn_threshold(values, valid)
        held_out = []
        for world in worlds:
            test = np.array([r["world"] == world for r in rows])
            rule = best_zero_fn_threshold(values[~test], valid[~test])
            base = {"world": world, "opportunities": int(test.sum()), "p_break_mass": float(p_break[test].sum()),
                    "search_s": float(cost[test].sum())}
            if "direction" not in rule:
                held_out.append({**base, "skipped": 0, "false_negatives": 0, "p_break_mass_lost": 0.0, "search_s_saved": 0.0})
                continue
            mask = values[test] < rule["threshold"] if rule["direction"] == "le" else values[test] > rule["threshold"]
            held_out.append({**base, "skipped": int(mask.sum()), "false_negatives": int((mask & valid[test]).sum()),
                             "p_break_mass_lost": float(p_break[test][mask].sum()),
                             "search_s_saved": float(cost[test][mask].sum())})
        summary["features"][name] = {"all_data_rule": fit, "all_data_share_skipped": fit["skipped"] / len(rows),
                                     "leave_one_world_out": held_out}
    text = json.dumps(summary, indent=1)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text(text + "\n")
    print(text)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = parser.add_subparsers(dest="mode", required=True)
    c = sub.add_parser("collect")
    c.add_argument("--save", required=True, type=Path)
    c.add_argument("--g0-script", required=True, type=Path)
    c.add_argument("--out", required=True, type=Path)
    c.add_argument("--warmup-steps", type=int, default=3)
    c.add_argument("--steps", type=int, default=10)
    a = sub.add_parser("analyze")
    a.add_argument("rows", nargs="+")
    a.add_argument("--out")
    args = parser.parse_args()
    collect(args) if args.mode == "collect" else analyze(args)


if __name__ == "__main__":
    main()
