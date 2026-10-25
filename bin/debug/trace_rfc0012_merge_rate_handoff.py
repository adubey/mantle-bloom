#!/usr/bin/env python3
"""Trace the first matched merge's stored-rate handoff and following motion solve.

Instrumentation only: wraps the existing merge and quasi-static solve, recording state
without changing model calculations. Run the same command serially with ``--model-root`` set
to the baseline and candidate checkouts. Raw traces belong in scratch storage.

    PYTHONHASHSEED=0 PYTHONPATH=<model-root>/backend \
      <model-root>/backend/.venv/bin/python \
      <trace-script> --model-root <model-root> --seed 580014012 --step-myr 2 \
      --stop-myr 18 --keep-plate 8 --absorb-plate 13 --out <scratch-output>.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any

import numpy as np

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("--model-root", type=Path, required=True)
parser.add_argument("--seed", type=int, default=580014012)
parser.add_argument("--step-myr", type=float, default=2.0)
parser.add_argument("--stop-myr", type=float, default=18.0)
parser.add_argument("--keep-plate", type=int, default=8)
parser.add_argument("--absorb-plate", type=int, default=13)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()

model_root = args.model_root.resolve()
sys.path.insert(0, str(model_root / "backend"))

from app import cratons, hm_ledger, merge_split, quad_merge, torque, world as world_mod  # noqa: E402
from app import continental_ledger  # noqa: E402
from app.hydroclimate import mobile_cover  # noqa: E402
from app.surface_fields import SURFACE_FIELDS  # noqa: E402


def _array_hash_update(digest: Any, values: np.ndarray) -> None:
    array = np.ascontiguousarray(values)
    digest.update(str(array.dtype).encode())
    digest.update(str(array.shape).encode())
    digest.update(array.tobytes())


def _initial_fingerprint(world) -> str:
    digest = hashlib.sha256()
    digest.update(f"{world.seed}:{world.node_density}".encode())
    for plate in sorted(world.plates, key=lambda item: item.plate_id):
        digest.update(str(plate.plate_id).encode())
        _array_hash_update(digest, plate._keys)
        _array_hash_update(digest, plate._frame)
        _array_hash_update(digest, np.asarray(plate.omega, dtype=float))
        for name in sorted(SURFACE_FIELDS):
            digest.update(name.encode())
            _array_hash_update(digest, plate.collect(name))
    return digest.hexdigest()


def _omega(plate) -> list[float]:
    return np.asarray(plate.omega, dtype=float).tolist()


def _ledger_snapshot(world, initial_hm: dict[str, float]) -> dict[str, Any]:
    inventory = hm_ledger.inventory_scopes_m3(world)
    accounts = hm_ledger.cumulative_scopes(world)
    hm_closure = {
        scope: inventory[scope]
        - initial_hm[scope]
        - accounts[scope]["source_m3"]
        + accounts[scope]["sink_m3"]
        - accounts[scope]["reclassification_m3"]
        for scope in initial_hm
    }
    return {
        "continental_balance_error_m3": continental_ledger.balance_error_m3(world),
        "craton_balance_error_m3": cratons.balance_error_m3(world),
        "mobile_cover_balance_error_m3": mobile_cover.balance_error_m3(world),
        "hm_closure_error_m3": hm_closure,
    }


def main() -> None:
    if args.step_myr <= 0 or args.stop_myr <= 0 or not np.isclose(args.stop_myr / args.step_myr, round(args.stop_myr / args.step_myr)):
        parser.error("stop-myr must be a positive multiple of step-myr")
    if args.out.exists():
        raise FileExistsError(f"refusing to overwrite trace: {args.out}")

    commit = subprocess.check_output(["git", "-C", str(model_root), "rev-parse", "HEAD"], text=True).strip()
    world = world_mod.generate_world(seed=args.seed, node_density=4.0)
    world.debug_diagnostics = True
    cratons.ensure_initialized(world)
    continental_ledger.ensure_initialized(world)
    mobile_cover.ensure_ledger(world)
    initial_hm = hm_ledger.inventory_scopes_m3(world)
    initial_fingerprint = _initial_fingerprint(world)
    original_merge = quad_merge.merge
    original_solve = torque.solve_plates_quasistatic
    original_topology = merge_split.apply_topology_changes
    trace: dict[str, Any] = {
        "commit": commit,
        "seed": args.seed,
        "node_density": 4.0,
        "step_myr": args.step_myr,
        "stop_myr": args.stop_myr,
        "initial_state_sha256": initial_fingerprint,
        "merge": None,
        "next_solve": None,
        "topology_events": [],
    }

    def traced_merge(keep, absorb, other_points_xyz, other_plates=(), world=None):
        matching = keep.plate_id == args.keep_plate and absorb.plate_id == args.absorb_plate
        before = {"keep_omega": _omega(keep), "absorb_omega": _omega(absorb)} if matching else None
        original_merge(keep, absorb, other_points_xyz, other_plates, world)
        if matching:
            trace["merge"] = {
                "step": world.steps_taken if world is not None else None,
                "elapsed_myr": world.elapsed_years / 1e6 if world is not None else None,
                "keep_plate_id": keep.plate_id,
                "absorb_plate_id": absorb.plate_id,
                "before": before,
                "survivor_omega_after_merge": _omega(keep),
                "survivor_node_count_after_merge": keep.node_count(),
            }

    def traced_solve(current_world, *, diagnostics=None):
        pre_solve = {plate.plate_id: _omega(plate) for plate in current_world.plates if plate.node_count()}
        solve_diagnostics: dict[str, Any] = {}
        solved = original_solve(current_world, diagnostics=solve_diagnostics)
        if current_world.steps_taken == round(16.0 / args.step_myr) + 1:
            trace["next_solve"] = {
                "step": current_world.steps_taken,
                "elapsed_myr_before_step": current_world.elapsed_years / 1e6,
                "pre_solve_omega": pre_solve,
                "solved_omega": {str(pid): np.asarray(rate, dtype=float).tolist() for pid, rate in solved.items()},
                "global": {
                    key: solve_diagnostics[key]
                    for key in (
                        "plates", "blocked_contacts", "yielded_contacts", "solver_converged",
                        "driving_torque_norm_nm", "residual_torque_norm_nm", "contact_force_max_n_per_m",
                    ) if key in solve_diagnostics
                },
                "plates": solve_diagnostics.get("plate_rates_cm_per_yr", []),
                "blocked_contact_responses": solve_diagnostics.get("blocked_contact_responses", []),
            }
        if diagnostics is not None:
            diagnostics.update(solve_diagnostics)
        return solved

    def traced_topology(current_world, years):
        for message in original_topology(current_world, years):
            trace["topology_events"].append({
                "step": current_world.steps_taken,
                "elapsed_myr": current_world.elapsed_years / 1e6,
                "message": message,
            })
            yield message

    quad_merge.merge = traced_merge
    torque.solve_plates_quasistatic = traced_solve
    merge_split.apply_topology_changes = traced_topology
    try:
        for _ in range(round(args.stop_myr / args.step_myr)):
            world_mod.step_world(world, args.step_myr * 1e6)
            if world.elapsed_years >= 16e6:
                # Keep the exact post-topology survivor rate at the event for comparison.
                if trace["merge"] is not None and "survivor_omega_after_step" not in trace["merge"]:
                    survivor = next((p for p in world.plates if p.plate_id == args.keep_plate), None)
                    trace["merge"]["survivor_omega_after_step"] = _omega(survivor) if survivor else None
    finally:
        quad_merge.merge = original_merge
        torque.solve_plates_quasistatic = original_solve
        merge_split.apply_topology_changes = original_topology

    trace["final"] = {
        "elapsed_myr": world.elapsed_years / 1e6,
        "plate_count": len(world.plates),
        "ledger": _ledger_snapshot(world, initial_hm),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(trace, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"out": str(args.out), "commit": commit, "merge_found": trace["merge"] is not None, "next_solve_found": trace["next_solve"] is not None}), flush=True)


if __name__ == "__main__":
    main()
