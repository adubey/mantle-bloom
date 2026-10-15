#!/usr/bin/env python3
"""Find the step function that leaves craton outside its column (issue #340).

A phase audit (`cratons.PhaseAudit.settle`) clips craton that sits above its column's Hc or
on a non-continental cell, and books it to a guessed account, so a one-off clip in
`split_craton_rework.py`'s per-site split doesn't say which writer caused it. This replays a
save to `--from-myr`, then steps through to `--to-myr` with taps on every audit's settle and
on the topology functions, printing each one that changes the craton outside columns
(`excess`) or the craton ledger by more than 1 km^3.

The first stretch is the slow part, so `--snapshot` saves the world at `--from-myr`, and a
later run with the same path starts there. `--seed` generates a world instead of loading a
save (`world.generate_world(seed)`, every other argument default: the #403 controls) and
starts probing at 0.

`--capture-dir` (RFC 0001's diagnosis gate, #440) saves every collision merge
(`merge_split._merge_plates_with_craton_audit`) with RFC 0001's signature: it books more
than `--min-share` of the live craton and also clips continental material above its column
(`collision_subducted_m3`). Each capture is a pickle holding the whole world just before the
merge (as `.mbworld` bytes) and the pair's ids, so the merge sees exactly the third plates
(`others`) and their footprints it saw in the run; `measure_merge_placement.py` replays it.
`--contacts` also writes RFC 0008's per-step contact rows (`probe_craton_contacts.Probe`),
with the craton-on-craton overlap volume of every overlapping pair, for the history before
each captured merge.

    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_clips.py \\
        ~/Downloads/mantle-bloom-seed910211954-220200000y.mbworld --from-myr 238.4 --to-myr 240.45 \\
        --snapshot /tmp/world-238.mbworld
    backend/.venv/bin/python bin/debug/probe_craton_clips.py --seed 7 --to-myr 400 \\
        --step-years 5e6 --capture-dir /tmp/rfc0001 --contacts /tmp/rfc0001/contacts-7.jsonl
"""

from __future__ import annotations

import argparse
import json
import pickle
import subprocess
import sys
import time
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


def _commit() -> str:
    root = Path(__file__).resolve().parents[2]
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    dirty = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", "backend/app"], cwd=root).returncode != 0
    return head + ("-dirty" if dirty else "")


def tap_capture(directory: Path, min_share: float, limit: int) -> None:
    """Save each collision merge with RFC 0001's signature (see the module docstring)."""
    original = merge_split._merge_plates_with_craton_audit
    directory.mkdir(parents=True, exist_ok=True)
    commit = _commit()
    saved = []

    def merge(world, id_keep: int, id_absorb: int) -> None:
        if len(saved) >= limit:
            return original(world, id_keep, id_absorb)
        data = persistence.save_world_bytes(world)
        live = cratons.live_volume_m3(world)
        sinks = cratons.booked_sinks_m3(world)
        material = world.continental_material_ledger.get("collision_subducted_m3", 0.0)
        original(world, id_keep, id_absorb)
        booked = cratons.booked_sinks_m3(world) - sinks
        clipped = world.continental_material_ledger.get("collision_subducted_m3", 0.0) - material
        if live > 0.0 and booked > min_share * live and clipped > 0.0:
            path = directory / f"seed{world.seed}-{world.elapsed_years / 1e6:.0f}myr-{id_keep}-{id_absorb}.pkl"
            meta = {"seed": world.seed, "myr": world.elapsed_years / 1e6, "keep": id_keep, "absorb": id_absorb,
                    "live_m3": live, "booked_m3": booked, "material_clipped_m3": clipped, "commit": commit}
            path.write_bytes(pickle.dumps({**meta, "world": data}))
            saved.append(path)
            print(f"  captured {json.dumps(meta)} -> {path}", flush=True)

    merge_split._merge_plates_with_craton_audit = merge


def craton_overlap_m3(world, a, b) -> float:
    """Craton on craton where `b`'s cratonic nodes sit inside `a`: the lesser of the two
    columns' cratons, by `b`'s cell area (RFC 0008's measure at fusion, each step here)."""
    craton_b = b.collect("craton_crust_m")
    cratonic = np.flatnonzero(craton_b > 0.0)
    if not len(cratonic) or not np.any(a.collect("craton_crust_m") > 0.0):
        return 0.0
    points = b.all_points_and_elevation()[0][cratonic]
    inside = a.contains_batch(points)
    if not np.any(inside):
        return 0.0
    _, idx = a.get_node_kdtree().query(points[inside])
    stacked = np.minimum(craton_b[cratonic[inside]], a.collect("craton_crust_m")[idx])
    return float(stacked @ b.node_areas_m2()[cratonic[inside]])


class ContactLog:
    """RFC 0008's contact rows each step, with each overlapping pair's craton-on-craton volume."""

    def __init__(self, world, path: Path) -> None:
        import probe_craton_contacts

        self.probe = probe_craton_contacts.Probe(world, 500.0)
        probe_craton_contacts.install(self.probe)
        path.parent.mkdir(parents=True, exist_ok=True)
        self.handle = path.open("w")

    def write(self, world) -> None:
        plates = {p.plate_id: p for p in world.plates}
        for row in self.probe.rows:
            for key, pair in row.get("pairs", {}).items():
                a, b = (plates.get(int(pid)) for pid in key.split("-"))
                if pair["overlap"] > 0.0 and a is not None and b is not None:
                    pair["craton_overlap_km3"] = round((craton_overlap_m3(world, a, b) + craton_overlap_m3(world, b, a)) / KM3, 1)
            self.handle.write(json.dumps(row) + "\n")
        self.handle.flush()
        self.probe.rows = []


def _start(args):
    """The world to probe from: a generated seed, the snapshot, or the save stepped to --from-myr."""
    if args.seed is not None:
        return world_mod.generate_world(args.seed)
    if args.snapshot is not None and args.snapshot.exists():
        return persistence.load_world_bytes(args.snapshot.read_bytes())
    world = persistence.load_world_bytes(args.save.read_bytes())
    world.reset_phase_budget()
    continental_ledger.ensure_initialized(world)
    cratons.ensure_initialized(world)
    while world.elapsed_years < args.from_myr * 1e6 - 0.5 * args.step_years:
        world_mod.step_world(world, args.step_years)
    if args.snapshot is not None:
        args.snapshot.write_bytes(persistence.save_world_bytes(world))
    return world


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("save", type=Path, nargs="?")
    parser.add_argument("--seed", type=int, help="generate this seed instead of loading a save")
    parser.add_argument("--from-myr", type=float, default=0.0)
    parser.add_argument("--to-myr", type=float, required=True)
    parser.add_argument("--step-years", type=float, default=100_000.0)
    parser.add_argument("--snapshot", type=Path, help="world saved at --from-myr; reused when it exists")
    parser.add_argument("--capture-dir", type=Path, help="save each merge with RFC 0001's signature here")
    parser.add_argument("--min-share", type=float, default=0.01, help="signature: craton booked / live craton")
    parser.add_argument("--max-captures", type=int, default=20)
    parser.add_argument("--contacts", type=Path, help="write RFC 0008's per-step contact rows here (JSON lines)")
    args = parser.parse_args()
    if (args.save is None) == (args.seed is None):
        parser.error("give a save or --seed, not both")

    world = _start(args)
    print(f"probing from {world.elapsed_years / 1e6:.2f} Myr", flush=True)
    tap_settle()
    for name in TOPOLOGY_FUNCTIONS:
        tap(merge_split, name)
    tap(world_mod, "_clamp_column_caps")
    contacts = ContactLog(world, args.contacts) if args.contacts is not None else None
    if args.capture_dir is not None:
        tap_capture(args.capture_dir, args.min_share, args.max_captures)
    started = time.time()
    while world.elapsed_years < args.to_myr * 1e6 - 1.0:
        world_mod.step_world(world, min(args.step_years, args.to_myr * 1e6 - world.elapsed_years))
        if contacts is not None:
            contacts.write(world)
    print(f"done at {world.elapsed_years / 1e6:.2f} Myr in {time.time() - started:.0f} s", flush=True)


if __name__ == "__main__":
    main()
