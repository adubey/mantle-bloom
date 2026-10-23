#!/usr/bin/env python3
"""Trace why continental plate pairs do or do not reach a merge in an RFC 0012 replay.

This opt-in diagnostic replaces ``find_continental_collision_pairs`` with an equivalent
instrumented implementation and wraps the existing collision/overlap timers. It records
filter-stage counts, nearby pair geometry and closing rates, sustained-collision progress,
timer/probability rolls, forced-overlap progress, and the merge actually selected. It does
not change model state. Run matched baseline/pilot seeds serially:

    PYTHONPATH=<model-root>/backend NUMBA_CACHE_DIR=/private/tmp/rfc0012-merge-numba \
      <model-root>/backend/.venv/bin/python \
      bin/debug/trace_rfc0012_merge_opportunities.py --model-root <model-root> \
      --seed 326494039 --steps 500 --step-myr 1 \
      --out /private/tmp/rfc0012-merge-opportunities/pilot-seed326494039.jsonl

The output is JSONL, one row per completed step plus a compact run summary. Full replay
outputs belong in ``/private/tmp``; do not check them into the repository.
"""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.spatial import cKDTree

from app import geometry, merge_split, plates, world as world_module
from app.boundary import MERGE_THRESHOLD_RAD, TRANSFORM_RATE_THRESHOLD, closing_rate
from app.elevation_lines import TARGET_LINE_SPACING_RAD, line_spacing_rad
from rfc0012_trace_common import run_trace_cli, write_trace_summary

EARTH_RADIUS_M = 6_371_000.0
SECONDS_PER_YEAR = 31_557_600.0


def _json(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(key): _json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(item) for item in value]
    return value


def _cm_per_year(rate_rad_per_year: float) -> float:
    return rate_rad_per_year * EARTH_RADIUS_M * 100.0


def _instrument_collision_search(world, runtime: dict) -> list[tuple[int, int]]:
    """Mirror the production filters while retaining each pair's rejection stage."""
    contact_rad = MERGE_THRESHOLD_RAD * (line_spacing_rad(world.node_density) / TARGET_LINE_SPACING_RAD)
    points = {}
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        mask = plates.node_is_continental(plate)
        if np.any(mask):
            points[plate.plate_id] = plate.all_points_and_elevation()[0][mask]

    continental = [plate for plate in world.plates if plate.plate_id in points]
    spheres = {pid: geometry.bounding_sphere(cloud) for pid, cloud in points.items()}
    trees = {}
    pairs = []
    details = []
    stages = Counter()

    for index, plate_a in enumerate(continental):
        for plate_b in continental[index + 1 :]:
            pid_a, pid_b = plate_a.plate_id, plate_b.plate_id
            cloud_a, cloud_b = points[pid_a], points[pid_b]
            stages["continental_pairs"] += 1
            center_a, radius_a = spheres[pid_a]
            center_b, radius_b = spheres[pid_b]
            center_distance = float(geometry.angular_distance(center_a, center_b))
            sphere_gap = center_distance - radius_a - radius_b
            if sphere_gap > contact_rad:
                stages["bounding_sphere_separated"] += 1
                continue
            stages["bounding_sphere_survivors"] += 1
            if pid_b not in trees:
                trees[pid_b] = cKDTree(cloud_b)
            candidate_mask = geometry.angular_distance(cloud_a, center_b) <= radius_b + contact_rad
            candidate_count = int(candidate_mask.sum())
            if candidate_count == 0:
                stages["no_point_candidates"] += 1
                details.append({
                    "pair": [pid_a, pid_b], "stage": "no_point_candidates",
                    "sphere_gap_rad": sphere_gap, "continental_nodes": [len(cloud_a), len(cloud_b)],
                })
                continue
            candidate_idx = np.flatnonzero(candidate_mask)
            distances, neighbor_idx = trees[pid_b].query(
                cloud_a[candidate_idx], workers=plates.query_workers(len(candidate_idx))
            )
            close = distances < contact_rad
            close_count = int(close.sum())
            detail = {
                "pair": [pid_a, pid_b],
                "continental_nodes": [len(cloud_a), len(cloud_b)],
                "sphere_gap_rad": sphere_gap,
                "candidate_points": candidate_count,
                "close_contact_points": close_count,
                "relative_rotation_cm_per_yr": _cm_per_year(float(np.linalg.norm(plate_a.omega - plate_b.omega))),
            }
            if close_count < merge_split.MERGE_MIN_CONTACT_NODES:
                stages["too_few_close_contacts"] += 1
                detail["stage"] = "too_few_close_contacts"
                details.append(detail)
                continue

            close_a = cloud_a[candidate_idx[close]]
            close_b = cloud_b[neighbor_idx[close]]
            rates = closing_rate(close_a, plate_a.omega, plate_b.omega, close_b)
            convergent = rates > TRANSFORM_RATE_THRESHOLD
            convergent_count = int(convergent.sum())
            detail.update({
                "convergent_contacts": convergent_count,
                "convergence_threshold_cm_per_yr": _cm_per_year(TRANSFORM_RATE_THRESHOLD),
                "closing_cm_per_yr": {
                    "mean": _cm_per_year(float(np.mean(rates))),
                    "max": _cm_per_year(float(np.max(rates))),
                    "convergent_mean": _cm_per_year(float(np.mean(rates[convergent]))) if convergent_count else None,
                },
            })
            if convergent_count < merge_split.MERGE_MIN_CONTACT_NODES:
                stages["insufficient_convergent_contacts"] += 1
                detail["stage"] = "insufficient_convergent_contacts"
                details.append(detail)
                continue

            stages["convergent_collision_candidates"] += 1
            detail["stage"] = "convergent_collision_candidate"
            details.append(detail)
            pairs.append(tuple(sorted((pid_a, pid_b))))

    runtime["collision_search"] = {
        "stages": dict(stages),
        "pairs": pairs.copy(),
        "pair_details": details,
        "contact_distance_rad": contact_rad,
        "min_contact_nodes": merge_split.MERGE_MIN_CONTACT_NODES,
        "plate_count": len(world.plates),
        "continental_plate_count": len(continental),
    }
    return pairs


def _install(world, output: Path, model_root: Path, seed: int, step_myr: float):
    runtime = {"step": 0, "counters": Counter()}
    original_finder = merge_split.find_continental_collision_pairs
    original_collision_update = merge_split.update_collision_progress
    original_overlap_update = merge_split.update_overlap_progress
    original_forced_pop = merge_split.pop_ready_forced_merge
    original_merge = merge_split._merge_plates_with_craton_audit
    original_step = world_module.step_world
    runtime["commit"] = subprocess.check_output(["git", "-C", str(model_root), "rev-parse", "HEAD"], text=True).strip()

    def traced_finder(current_world):
        return _instrument_collision_search(current_world, runtime)

    def traced_collision_update(current_world, years):
        progress_before = dict(current_world.collision_progress)
        ready = original_collision_update(current_world, years)
        search = runtime.get("collision_search", {})
        current_pairs = set(search.get("pairs", []))
        jammed = merge_split.find_jammed_pairs(current_world)
        live = {
            plate.plate_id for plate in current_world.plates
            if plate.node_count() > 0 and np.any(plates.node_is_continental(plate))
        }
        retained_jam_pairs = {
            pair for pair in progress_before
            if pair in jammed and pair not in current_pairs and pair[0] in live and pair[1] in live
        }
        tracked_pairs = sorted(current_pairs | retained_jam_pairs)
        decisions = []
        for pair in tracked_pairs:
            accumulated = progress_before.get(pair, 0.0) + years
            threshold = merge_split._collision_threshold_years(current_world.seed, pair)
            timer_ready = accumulated >= threshold
            probability = merge_split._merge_probability(current_world, pair) if timer_ready else None
            draw = None
            passed = None
            if timer_ready:
                rng = np.random.default_rng((current_world.seed, pair[0], pair[1], round(current_world.elapsed_years)))
                draw = float(rng.random())
                passed = draw < probability
            decisions.append({
                "pair": pair,
                "source": "convergent_contact" if pair in current_pairs else "retained_craton_jam",
                "jammed_contact_fraction": jammed.get(pair),
                "progress_before_years": progress_before.get(pair, 0.0),
                "progress_after_years": current_world.collision_progress.get(pair),
                "timer_threshold_years": threshold,
                "timer_ready": timer_ready,
                "merge_probability": probability,
                "random_draw": draw,
                "probability_pass": passed,
                "returned_ready": pair in ready,
            })
            runtime["counters"]["collision_timer_ready"] += int(timer_ready)
            runtime["counters"]["collision_probability_pass"] += int(passed is True)
            runtime["counters"]["collision_probability_fail"] += int(passed is False)
        expired = [
            {"pair": pair, "progress_lost_years": progress}
            for pair, progress in progress_before.items()
            if pair not in current_world.collision_progress
        ]
        runtime["collision_progress"] = {
            "pairs": decisions,
            "expired_pairs": expired,
            "ready_pairs": ready.copy(),
            "selected_by_one_merge_limit": ready[0] if ready else None,
            "other_ready_deferred": ready[1:],
            "ready_forced_overlap_pairs": [
                pair for pair, progress in current_world.overlap_progress.items()
                if progress >= merge_split.FORCED_MERGE_SUSTAINED_YEARS
                and pair[0] in live and pair[1] in live
            ],
        }
        forced_ready = runtime["collision_progress"]["ready_forced_overlap_pairs"]
        runtime["forced_merge"] = {
            "ready_pairs": forced_ready,
            "selected_pair": None,
            "deferred_by_collision_merge": bool(ready and forced_ready),
        }
        runtime["counters"]["collision_pairs_current"] += len(current_pairs)
        runtime["counters"]["collision_pairs_jam_retained"] += len(retained_jam_pairs)
        runtime["counters"]["collision_pairs_expired"] += len(expired)
        runtime["counters"]["ready_collision_pairs"] += len(ready)
        predicted_ready = {item["pair"] for item in decisions if item["probability_pass"] is True}
        runtime["counters"]["collision_ready_decision_mismatches"] += int(predicted_ready != set(ready))
        return ready

    def traced_overlap_update(current_world, years, overlap):
        before = dict(current_world.overlap_progress)
        fractions = merge_split._deep_continental_overlap_fractions(current_world, overlap)
        original_overlap_update(current_world, years, overlap)
        rows = []
        for pair, fraction in fractions.items():
            if fraction > 0.0 or pair in before:
                after = current_world.overlap_progress.get(pair)
                rows.append({
                    "pair": pair,
                    "deep_overlap_fraction": fraction,
                    "required_fraction": merge_split.FORCED_MERGE_OVERLAP_FRACTION,
                    "progress_before_years": before.get(pair, 0.0),
                    "progress_after_years": after,
                    "qualified_this_step": fraction >= merge_split.FORCED_MERGE_OVERLAP_FRACTION,
                    "timer_ready_after_update": bool(after is not None and after >= merge_split.FORCED_MERGE_SUSTAINED_YEARS),
                })
        reset = [pair for pair in before if pair not in current_world.overlap_progress]
        runtime["overlap_progress_update"] = {"pairs": rows, "reset_pairs": reset}
        runtime["counters"]["forced_overlap_qualified_steps"] += sum(row["qualified_this_step"] for row in rows)
        runtime["counters"]["forced_overlap_resets"] += len(reset)

    def traced_forced_pop(current_world):
        progress_before = dict(current_world.overlap_progress)
        live = {
            plate.plate_id for plate in current_world.plates
            if plate.node_count() > 0 and np.any(plates.node_is_continental(plate))
        }
        ready = [
            pair for pair, years in progress_before.items()
            if years >= merge_split.FORCED_MERGE_SUSTAINED_YEARS and pair[0] in live and pair[1] in live
        ]
        selected = original_forced_pop(current_world)
        runtime["forced_merge"] = {
            "ready_pairs": ready,
            "selected_pair": selected,
            "deferred_by_collision_merge": False,
        }
        runtime["counters"]["forced_merge_ready_pairs"] += len(ready)
        runtime["counters"]["forced_merge_events"] += int(selected is not None)
        return selected

    def traced_merge(current_world, id_keep, id_absorb):
        original_merge(current_world, id_keep, id_absorb)
        pair = tuple(sorted((id_keep, id_absorb)))
        collision_selected = runtime.get("collision_progress", {}).get("selected_by_one_merge_limit")
        runtime["actual_merge"] = {
            "keep_plate_id": id_keep,
            "absorbed_plate_id": id_absorb,
            "path": "sustained_collision" if pair == collision_selected else "forced_overlap",
        }
        runtime["counters"]["actual_merges"] += 1

    def traced_step(current_world, years, **kwargs):
        runtime["step"] = current_world.steps_taken + 1
        runtime["collision_search"] = {}
        runtime["collision_progress"] = {}
        runtime["overlap_progress_update"] = {}
        runtime["forced_merge"] = {}
        runtime["actual_merge"] = None
        started = time.perf_counter()
        original_step(current_world, years, **kwargs)
        row = {
            "kind": "step",
            "model": str(model_root),
            "commit": runtime["commit"],
            "seed": seed,
            "step": current_world.steps_taken,
            "elapsed_myr": current_world.elapsed_years / 1e6,
            "runtime_seconds": time.perf_counter() - started,
            "collision_search": runtime["collision_search"],
            "collision_progress": runtime["collision_progress"],
            "overlap_progress_update": runtime["overlap_progress_update"],
            "forced_merge": runtime["forced_merge"],
            "actual_merge": runtime["actual_merge"],
            "live_plate_count": len(current_world.plates),
        }
        with output.open("a") as handle:
            handle.write(json.dumps(_json(row), separators=(",", ":"), allow_nan=False) + "\n")
        return None

    merge_split.find_continental_collision_pairs = traced_finder
    merge_split.update_collision_progress = traced_collision_update
    merge_split.update_overlap_progress = traced_overlap_update
    merge_split.pop_ready_forced_merge = traced_forced_pop
    merge_split._merge_plates_with_craton_audit = traced_merge
    world_module.step_world = traced_step

    def restore():
        merge_split.find_continental_collision_pairs = original_finder
        merge_split.update_collision_progress = original_collision_update
        merge_split.update_overlap_progress = original_overlap_update
        merge_split.pop_ready_forced_merge = original_forced_pop
        merge_split._merge_plates_with_craton_audit = original_merge
        world_module.step_world = original_step

    return runtime, restore


def run(seed: int, steps: int, step_myr: float, output: Path, model_root: Path) -> None:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite trace: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    world = world_module.generate_world(seed=seed, node_density=1.0)
    runtime, restore = _install(world, output, model_root, seed, step_myr)
    started = time.perf_counter()
    try:
        for index in range(steps):
            world_module.step_world(world, step_myr * 1_000_000.0)
            if (index + 1) % 25 == 0:
                print(
                    f"seed={seed} step={index + 1}/{steps} "
                    f"elapsed_myr={(index + 1) * step_myr:g} "
                    f"runtime_min={(time.perf_counter() - started) / 60:.1f}",
                    flush=True,
                )
    finally:
        restore()
    summary = {
        "kind": "summary",
        "seed": seed,
        "steps": steps,
        "step_myr": step_myr,
        "elapsed_myr": steps * step_myr,
        "model": str(model_root),
        "commit": runtime["commit"],
        "runtime_seconds": time.perf_counter() - started,
        "counters": dict(runtime["counters"]),
        "final_collision_progress_pairs": len(world.collision_progress),
        "final_overlap_progress_pairs": len(world.overlap_progress),
    }
    write_trace_summary(output, summary)


def main() -> None:
    run_trace_cli(__doc__, run)


if __name__ == "__main__":
    main()
