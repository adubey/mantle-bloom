#!/usr/bin/env python3
"""Find the step function that leaves craton outside its column (issue #340).

A phase audit (`cratons.PhaseAudit.settle`) clips craton that sits above its column's Hc or
on a non-continental cell, and books it to a guessed account, so a one-off clip in
`split_craton_rework.py`'s per-site split doesn't say which writer caused it. This replays a
save to `--from-myr`, then steps through to `--to-myr` with taps on every audit's settle and
on the topology functions, printing each one that changes the craton outside columns
(`excess`) or the craton ledger by more than 1 km^3.

The first stretch is the slow part, so `--snapshot` saves the world at `--from-myr`, and a
later run with the same path starts there. `--capture` pickles the inputs of the first plate
fusion (`merge_split._fuse_plates`) that leaves craton outside its column, as
`(keep, absorb, other_points)`, so that fusion can be replayed alone with
`keep.merge_with(absorb, other_points)`.

    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_clips.py \\
        ~/Downloads/mantle-bloom-seed910211954-220200000y.mbworld --from-myr 238.4 --to-myr 240.45 \\
        --snapshot /tmp/world-238.mbworld --capture /tmp/fuse.pkl
"""

from __future__ import annotations

import argparse
import copy
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import continental_ledger, cratons, merge_split, persistence, world as world_mod  # noqa: E402
from app.elevation_lines import effective_is_continental_from_codes  # noqa: E402

KM3 = 1e9
TOPOLOGY_FUNCTIONS = (
    "apply_topology_changes",
    "_merge_plates_with_craton_audit",
    "remove_defunct_plates",
    "_accrete_stranded_terrane",
    "maybe_split_plate",
    "defragment_plates",
)


def excess_m3(world) -> float:
    """Craton a phase audit would clip: above its column's Hc or on a non-continental cell."""
    total = 0.0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        craton = plate.collect("craton_crust_m")
        if not np.any(craton > 0.0):
            continue
        hc = plate.collect("crustal_thickness_m")
        continental = effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        kept = np.where(continental, np.clip(craton, 0.0, np.maximum(hc, 0.0)), 0.0)
        total += float((craton - kept) @ plate.node_areas_m2())
    return total


def tap(module, name: str) -> None:
    """Print `module.name`'s effect on the craton outside columns and on the ledger."""
    original = getattr(module, name)

    def tapped(world, *args, **kwargs):
        before, ledger = excess_m3(world), dict(world.craton_ledger)
        result = original(world, *args, **kwargs)
        after = excess_m3(world)
        moved = {k: (v - ledger.get(k, 0.0)) / KM3 for k, v in world.craton_ledger.items() if abs(v - ledger.get(k, 0.0)) > KM3}
        if abs(after - before) > KM3 or moved:
            print(
                f"  {name} at {world.elapsed_years / 1e6:.1f} Myr: outside columns {before / KM3:,.0f} -> "
                f"{after / KM3:,.0f} km3, ledger {moved}",
                flush=True,
            )
        return result

    setattr(module, name, tapped)


def tap_settle() -> None:
    original = cratons.PhaseAudit.settle

    def settle(self, clip_account, residual_account="unattributed_m3"):
        excess = excess_m3(self.world)
        if excess > KM3:
            print(f"  settle({clip_account}) clips {excess / KM3:,.0f} km3 at {self.world.elapsed_years / 1e6:.1f} Myr", flush=True)
        return original(self, clip_account, residual_account)

    cratons.PhaseAudit.settle = settle


def tap_capture(path: Path) -> None:
    original = merge_split._fuse_plates

    def fuse(world, keep, absorb, phase):
        if path.exists():
            return original(world, keep, absorb, phase)
        others = [p.all_points_and_elevation()[0] for p in world.plates if p is not keep and p is not absorb]
        snapshot = (copy.deepcopy(keep), copy.deepcopy(absorb), np.concatenate(others) if others else np.zeros((0, 3)))
        before = excess_m3(world)
        result = original(world, keep, absorb, phase)
        if excess_m3(world) > before + KM3:
            path.write_bytes(pickle.dumps(snapshot))
            print(f"  captured {phase} of plate {absorb.plate_id} into {keep.plate_id} -> {path}", flush=True)
        return result

    merge_split._fuse_plates = fuse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("save", type=Path)
    parser.add_argument("--from-myr", type=float, required=True)
    parser.add_argument("--to-myr", type=float, required=True)
    parser.add_argument("--step-years", type=float, default=100_000.0)
    parser.add_argument("--snapshot", type=Path, help="world saved at --from-myr; reused when it exists")
    parser.add_argument("--capture", type=Path, help="pickle the first fusion that leaves craton outside its column")
    args = parser.parse_args()

    if args.snapshot is not None and args.snapshot.exists():
        world = persistence.load_world_bytes(args.snapshot.read_bytes())
    else:
        world = persistence.load_world_bytes(args.save.read_bytes())
        world.reset_phase_budget()
        continental_ledger.ensure_initialized(world)
        cratons.ensure_initialized(world)
        while world.elapsed_years < args.from_myr * 1e6 - 0.5 * args.step_years:
            world_mod.step_world(world, args.step_years)
        if args.snapshot is not None:
            args.snapshot.write_bytes(persistence.save_world_bytes(world))
    print(f"probing from {world.elapsed_years / 1e6:.2f} Myr", flush=True)

    tap_settle()
    for name in TOPOLOGY_FUNCTIONS:
        tap(merge_split, name)
    tap(world_mod, "_clamp_column_caps")
    if args.capture is not None:
        tap_capture(args.capture)
    while world.elapsed_years < args.to_myr * 1e6:
        world_mod.step_world(world, args.step_years)
    print(f"done at {world.elapsed_years / 1e6:.2f} Myr", flush=True)


if __name__ == "__main__":
    main()
