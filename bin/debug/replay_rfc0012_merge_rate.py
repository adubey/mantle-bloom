#!/usr/bin/env python3
"""Run one RFC 0012 #547 merge-only replay from a generated world.

This is one arm of the preregistered baseline/candidate matrix in
`analysis/issue547/README.md`. Invoke it once per seed and step size, serially across
baseline and candidate worktrees. It writes timestamped progress and error logs, checkpoint
summaries, cumulative motion/deformation/ledger measures, and a restartable world save at
30 Myr by default. Resume from a summary checkpoint with `--resume-from`.

    PYTHONHASHSEED=0 backend/.venv/bin/python bin/debug/replay_rfc0012_merge_rate.py \
        --seed 326494039 --step-myr 0.1 --myr 50 --out <scratch-output>.json
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import resource
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

backend_root = Path(
    os.environ.get("MANTLE_BLOOM_REPLAY_BACKEND", str(Path(__file__).resolve().parents[2] / "backend"))
).resolve()
sys.path.insert(0, str(backend_root))

from app import continental_ledger, cratons, elevation_lines, hm_ledger, merge_split, persistence, world as world_mod
from app.hydroclimate import mobile_cover


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    if not len(values):
        return float("nan")
    order = np.argsort(values)
    cumulative = np.cumsum(weights[order])
    return float(values[order[np.searchsorted(cumulative, q * cumulative[-1])]])


def _snapshot(world) -> dict:
    spacing = elevation_lines.line_spacing_rad(world.node_density)
    hc, hm, area, land = [], [], [], []
    for plate in world.plates:
        if not plate.node_count():
            continue
        hc.append(np.asarray(plate.collect("crustal_thickness_m"), dtype=float))
        hm.append(np.asarray(plate.collect("mantle_lithosphere_thickness_m"), dtype=float))
        area.append(np.asarray(plate.accounting_areas_m2(spacing), dtype=float))
        land.append(np.asarray(plate.collect("elevation"), dtype=float) > world.sea_level_m)
    hc, hm, area, land = (np.concatenate(parts) for parts in (hc, hm, area, land))
    total_area = float(area.sum())
    return {
        "myr": world.elapsed_years / 1e6,
        "plate_count": len(world.plates),
        "land_fraction": float(area[land].sum() / total_area) if total_area else 0.0,
        "hc_area_weighted_km": {
            f"p{int(q * 100)}": _weighted_quantile(hc, area, q) / 1e3 for q in (0.5, 0.9, 0.99)
        },
        "hm_area_weighted_km": {
            f"p{int(q * 100)}": _weighted_quantile(hm, area, q) / 1e3 for q in (0.5, 0.9, 0.99)
        },
        "ledger_balance_error_m3": {
            "continental": continental_ledger.balance_error_m3(world),
            "craton": cratons.balance_error_m3(world),
            "mobile_cover": mobile_cover.balance_error_m3(world),
            "hm_by_scope": {},
        },
        "hm_inventory_m3": hm_ledger.inventory_scopes_m3(world),
        "hm_accounts_m3": hm_ledger.cumulative_scopes(world),
    }


def _phase_totals(budget: dict) -> dict:
    totals = {}
    for phase, entry in budget.items():
        scopes = entry.get("scopes")
        if scopes and "all" in scopes:
            row = scopes["all"]
            totals[phase] = {
                "hc_m3": row["hc_volume_after_m3"] - row["hc_volume_before_m3"],
                "hm_m3": row["hm_volume_after_m3"] - row["hm_volume_before_m3"],
            }
    return totals


def _hm_closure(snapshot: dict, initial_inventory: dict) -> dict:
    return {
        scope: snapshot["hm_inventory_m3"][scope]
        - initial_inventory[scope]
        - snapshot["hm_accounts_m3"][scope]["source_m3"]
        + snapshot["hm_accounts_m3"][scope]["sink_m3"]
        - snapshot["hm_accounts_m3"][scope]["reclassification_m3"]
        for scope in initial_inventory
    }


def _configure_logging(path: Path) -> logging.Logger:
    path.parent.mkdir(parents=True, exist_ok=True)
    root = logging.getLogger()
    root.setLevel(logging.INFO)
    root.handlers.clear()
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    for handler in (logging.FileHandler(path, mode="a"), logging.StreamHandler()):
        handler.setFormatter(formatter)
        root.addHandler(handler)
    logging.captureWarnings(True)
    return logging.getLogger("rfc0012_merge_replay")


def _write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=1, default=dict) + "\n")
    temporary.replace(path)


def _peak_rss_mib() -> float:
    # getrusage reports bytes on macOS/BSD and KiB on Linux.
    divisor = 2**10 if sys.platform == "linux" else 2**20
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / divisor


def _topology_events_for_previous_step(topology_events: list[list[Any]], elapsed_years: float, step_years: float) -> list[str]:
    previous_step_years = elapsed_years - step_years
    return [
        message
        for event_years, message in topology_events
        if np.isclose(event_years, previous_step_years, rtol=0.0, atol=1.0)
    ]


def _load_resume(path: Path, args) -> tuple[object, dict]:
    checkpoint = json.loads(path.read_text())
    for key, expected in (("seed", args.seed), ("step_myr", args.step_myr)):
        if checkpoint.get(key) != expected:
            raise ValueError(f"Resume checkpoint {key}={checkpoint.get(key)!r}, expected {expected!r}")
    save_path = checkpoint.get("world_save_path")
    if not save_path or not Path(save_path).is_file():
        raise ValueError(f"Resume checkpoint does not name an existing world save: {save_path!r}")
    world = persistence.load_world_bytes(Path(save_path).read_bytes())
    world.debug_diagnostics = True
    save_myr = float(checkpoint["world_save_myr"])
    save_index = next(
        (i for i, row in enumerate(checkpoint["checkpoints"]) if abs(float(row["myr"]) - save_myr) < 1e-9),
        None,
    )
    if save_index is None:
        raise ValueError(f"Summary has no checkpoint matching world save at {save_myr:g} Myr")
    # A crash can leave the JSON several checkpoints ahead of its one restart save. Roll
    # cumulative diagnostics back to the saved checkpoint before continuing from that world.
    checkpoint["checkpoints"] = checkpoint["checkpoints"][: save_index + 1]
    checkpoint["topology_events"] = [event for event in checkpoint["topology_events"] if event[0] <= world.elapsed_years]
    saved_row = checkpoint["checkpoints"][-1]
    if "motion_cumulative" not in saved_row:
        raise ValueError("Resume summary predates cumulative checkpoint metrics; restart the replay to preserve ledger accuracy")
    checkpoint["motion"] = saved_row.get("motion_cumulative", {})
    deformation = defaultdict(lambda: defaultdict(float))
    for row in checkpoint["checkpoints"]:
        for phase, measures in row.get("deformation_since_checkpoint", {}).items():
            for measure, value in measures.items():
                deformation[phase][measure] += value
    checkpoint["deformation"] = {phase: dict(values) for phase, values in deformation.items()}
    checkpoint["runtime_seconds"] = saved_row.get("runtime_seconds", 0.0)
    checkpoint["peak_rss_mib"] = saved_row.get("peak_rss_mib", 0.0)
    expected_years = save_myr * 1e6
    if abs(world.elapsed_years - expected_years) > 1.0:
        raise ValueError(
            f"Saved world is at {world.elapsed_years / 1e6:g} Myr, summary is at {expected_years / 1e6:g} Myr"
        )
    result = dict(checkpoint)
    result["motion"] = defaultdict(float, checkpoint["motion"])
    result["deformation"] = defaultdict(
        lambda: defaultdict(float),
        {phase: defaultdict(float, totals) for phase, totals in checkpoint["deformation"].items()},
    )
    result["deformation_window"] = defaultdict(lambda: defaultdict(float))
    result["resumed_from"] = str(path)
    return world, result


def _validate_intervals(parser, args) -> int:
    save_step = round(args.world_save_myr / args.step_myr)
    if abs(save_step * args.step_myr - args.world_save_myr) > 1e-9:
        parser.error("--world-save-myr must fall on a timestep boundary")
    checkpoint_steps = round(args.checkpoint_myr / args.step_myr)
    if args.checkpoint_myr <= 0 or abs(checkpoint_steps * args.step_myr - args.checkpoint_myr) > 1e-9:
        parser.error("--checkpoint-myr must be a positive multiple of --step-myr")
    if abs(args.world_save_myr / args.checkpoint_myr - round(args.world_save_myr / args.checkpoint_myr)) > 1e-9:
        parser.error("--world-save-myr must fall on a checkpoint boundary")
    if args.world_save_myr > args.myr:
        parser.error("--world-save-myr must not exceed --myr")
    return checkpoint_steps


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--step-myr", type=float, choices=(0.0125, 0.025, 0.05, 0.1, 1.0, 2.0), required=True)
    parser.add_argument("--myr", type=float, default=50.0)
    parser.add_argument("--checkpoint-myr", type=float, default=10.0)
    parser.add_argument("--world-save-myr", type=float, default=30.0)
    parser.add_argument("--world-save-path", type=Path)
    parser.add_argument("--resume-from", type=Path, help="resume from a checkpoint JSON with its paired .mbworld save")
    parser.add_argument("--log", type=Path, help="append timestamped progress and errors here; defaults beside --out")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    logger = _configure_logging(args.log or args.out.with_suffix(".log"))
    checkpoint_steps = _validate_intervals(parser, args)
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=backend_root.parent, check=True, capture_output=True, text=True
    ).stdout.strip()
    seconds = args.step_myr * 1e6
    steps = round(args.myr / args.step_myr)
    if args.resume_from:
        world, result = _load_resume(args.resume_from, args)
        if args.myr < world.elapsed_years / 1e6:
            raise ValueError("--myr must be at or beyond the saved world's current time")
        initial_hm = result["checkpoints"][0]["hm_inventory_m3"]
        start_index = round(world.elapsed_years / seconds)
        world_save_path = Path(result["world_save_path"])
        if not world_save_path.is_file():
            raise ValueError(f"Saved world disappeared before resume: {world_save_path}")
        logger.info("Resuming from %.6f Myr checkpoint %s", world.elapsed_years / 1e6, args.resume_from)
    else:
        world = world_mod.generate_world(seed=args.seed, node_density=4.0)
        world.debug_diagnostics = True
        cratons.ensure_initialized(world)
        continental_ledger.ensure_initialized(world)
        mobile_cover.ensure_ledger(world)
        initial = _snapshot(world)
        initial_hm = initial["hm_inventory_m3"]
        result = {
            "commit": commit,
            "seed": args.seed,
            "node_density": 4.0,
            "step_myr": args.step_myr,
            "horizon_myr": args.myr,
            "checkpoints": [initial],
            "topology_events": [],
            "motion": defaultdict(float),
            "deformation": defaultdict(lambda: defaultdict(float)),
            "deformation_window": defaultdict(lambda: defaultdict(float)),
            "world_save_path": None,
            "world_save_myr": None,
        }
        start_index = 0
        logger.info("Generated seed %d at node density 4; model commit %s", args.seed, commit)
    previous_phases = _phase_totals(world.phase_budget)
    original_topology = merge_split.apply_topology_changes

    def capture_topology_events(step_world, years):
        for message in original_topology(step_world, years):
            result["topology_events"].append([step_world.elapsed_years, message])
            yield message

    merge_split.apply_topology_changes = capture_topology_events
    started = time.time()
    runtime_before = float(result.get("runtime_seconds", 0.0))
    current_step = start_index
    try:
        for index in range(start_index + 1, steps + 1):
            current_step = index
            motion = []
            world_mod.step_world(world, seconds, motion_diagnostics=motion)
            if motion:
                diag = motion[-1]
                result["motion"]["accumulated_driving_torque_nm_year"] += diag["driving_torque_norm_nm"] * seconds
                result["motion"]["accumulated_residual_torque_nm_year"] += diag["residual_torque_norm_nm"] * seconds
                rates = [row["rate_cm_per_yr"] for row in diag["plate_rates_cm_per_yr"]]
                result["motion"]["rate_sum_cm_per_yr"] += sum(rates)
                result["motion"]["rate_samples"] += len(rates)
                result["motion"]["peak_solved_rate_cm_per_yr"] = max(
                    result["motion"]["peak_solved_rate_cm_per_yr"], max(rates, default=0.0)
                )
                contacts = diag["blocked_contact_responses"]
                result["motion"]["contact_samples"] += len(contacts)
                result["motion"]["absolute_closing_sum_m_per_s"] += sum(
                    abs(c["closing_rate_m_per_s"]) for c in contacts
                )
                for contact in contacts:
                    closing = abs(contact["closing_rate_m_per_s"])
                    if closing > result["motion"]["peak_absolute_closing_m_per_s"]:
                        result["motion"]["peak_absolute_closing_m_per_s"] = closing
                        result["motion"]["peak_closing_contact"] = {
                            "step_index": index,
                            "step_myr": world.elapsed_years / 1e6,
                            "signed_closing_rate_m_per_s": contact["closing_rate_m_per_s"],
                            "absolute_closing_rate_m_per_s": closing,
                            "lower_plate_id": contact["lower_plate_id"],
                            "upper_plate_id": contact["upper_plate_id"],
                            "force_n_per_m": contact["force_n_per_m"],
                            "capacity_n_per_m": contact["capacity_n_per_m"],
                            "yielded": contact["yielded"],
                            "plate_count": len(world.plates),
                            "topology_events_previous_step": _topology_events_for_previous_step(
                                result["topology_events"], world.elapsed_years, seconds
                            ),
                        }
                result["motion"]["contact_force_sum_n_per_m"] += sum(
                    c["force_n_per_m"] for c in contacts
                )
            phases = _phase_totals(world.phase_budget)
            for phase in phases.keys() | previous_phases.keys():
                for measure in ("hc_m3", "hm_m3"):
                    delta = phases.get(phase, {}).get(measure, 0.0) - previous_phases.get(phase, {}).get(measure, 0.0)
                    result["deformation"][phase][measure] += delta
                    result["deformation_window"][phase][measure] += delta
            previous_phases = phases
            logger.info(
                "step %d/%d complete: %.6f Myr, %d plates, peak closing %.6e m/s, RSS %.0f MiB",
                index,
                steps,
                world.elapsed_years / 1e6,
                len(world.plates),
                result["motion"].get("peak_absolute_closing_m_per_s", 0.0),
                _peak_rss_mib(),
            )
            if index % checkpoint_steps == 0 or index == steps:
                row = _snapshot(world)
                row["hm_ledger_closure_error_m3"] = _hm_closure(row, initial_hm)
                row["topology_events_so_far"] = len(result["topology_events"])
                row["deformation_since_checkpoint"] = {
                    phase: dict(values) for phase, values in result["deformation_window"].items()
                }
                row["motion_cumulative"] = dict(result["motion"])
                result["deformation_window"] = defaultdict(lambda: defaultdict(float))
                result["checkpoints"].append(row)
                args.out.parent.mkdir(parents=True, exist_ok=True)
                result["runtime_seconds"] = runtime_before + time.time() - started
                result["peak_rss_mib"] = _peak_rss_mib()
                row["runtime_seconds"] = result["runtime_seconds"]
                row["peak_rss_mib"] = result["peak_rss_mib"]
                superseded_save_path = None
                if world.elapsed_years / 1e6 >= args.world_save_myr:
                    previous_save = result.get("world_save_path")
                    base_save_path = args.world_save_path or args.out.with_suffix("")
                    save_path = base_save_path.with_name(
                        f"{base_save_path.stem}-{world.elapsed_years / 1e6:g}myr.mbworld"
                    )
                    save_path.parent.mkdir(parents=True, exist_ok=True)
                    temporary = save_path.with_suffix(save_path.suffix + ".tmp")
                    temporary.write_bytes(persistence.save_world_bytes(world))
                    temporary.replace(save_path)
                    result["world_save_path"] = str(save_path)
                    result["world_save_myr"] = world.elapsed_years / 1e6
                    if previous_save and Path(previous_save) != save_path:
                        superseded_save_path = Path(previous_save)
                    logger.info("Saved restartable world at %.6f Myr to %s", world.elapsed_years / 1e6, save_path)
                args.out.parent.mkdir(parents=True, exist_ok=True)
                _write_json(args.out, result)
                if superseded_save_path is not None:
                    superseded_save_path.unlink(missing_ok=True)
                logger.info("Wrote summary checkpoint at %.6f Myr to %s", world.elapsed_years / 1e6, args.out)
    except BaseException:
        logger.exception(
            "Replay interrupted at step %d/%d (%.6f Myr); last world save: %s",
            current_step,
            steps,
            world.elapsed_years / 1e6,
            result.get("world_save_path"),
        )
        raise
    motion = result["motion"]
    rate_samples = max(motion.pop("rate_samples"), 1)
    contact_samples = max(motion.pop("contact_samples"), 1)
    motion["mean_solved_rate_cm_per_yr"] = motion.pop("rate_sum_cm_per_yr") / rate_samples
    motion["mean_absolute_closing_m_per_s"] = motion.pop("absolute_closing_sum_m_per_s") / contact_samples
    motion["mean_contact_force_n_per_m"] = motion.pop("contact_force_sum_n_per_m") / contact_samples
    result["runtime_seconds"] = runtime_before + time.time() - started
    result["peak_rss_mib"] = _peak_rss_mib()
    result["motion"] = dict(motion)
    result["deformation"] = {phase: dict(values) for phase, values in result["deformation"].items()}
    _write_json(args.out, result)
    logger.info("Replay complete at %.6f Myr; final summary written to %s", world.elapsed_years / 1e6, args.out)


if __name__ == "__main__":
    main()
