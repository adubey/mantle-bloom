#!/usr/bin/env python3
"""Attribute divergent and stretch craton loss to the split cut that made the boundary
(RFC 0003 phase 1b, issue #388).

Phase 1 (#378) split `rifted_m3` by path, but couldn't say whether divergent thinning and
stretch donation run on boundaries a craton-crossing cut created. This replays a seed with
probe-only surface fields, registered in this process alone, so no saved state changes:

- at each successful split, the nodes of both daughters within the cut's band are marked
  "through" if the band held craton before the split, else "clear". A later cut over the same
  node overwrites its mark. The marks are intensive fractions in [0, 1], so a remap that blends
  a marked cell with an unmarked one keeps the marked share of its area;
- every later divergent (`lithosphere_plate.deform_columns`) and stretch
  (`quad_tectonics._open_rift`) craton loss is booked per node by those fractions: cut through
  craton, cut clear of craton, or never in a cut band (the rest).

The band is the failed-rift band phase 1 counted with (`FAILED_RIFT_BAND_MULT` spacings, the
"band" scheme). Two sensitivity schemes ride along: "ever", where a "through" mark is never
overwritten by a later clear cut (an upper bound), and "wide", the same as "band" but twice
as wide, for boundary loss that eats past the original band.

F, the share of `rifted_m3` steering the cut (solution 1) could remove, is the failed-rift
loss plus the divergent and stretch loss marked "through". The run starts at year 0 and goes
to --years (1 Gyr in the RFC), so every mark is younger than 1 Gyr.

    backend/.venv/bin/python bin/debug/attribute_rift_boundary_loss.py --seed 7 --years 1e9 \\
        --out /tmp/issue388-seed7.jsonl
    backend/.venv/bin/python bin/debug/attribute_rift_boundary_loss.py --summarize /tmp/issue388-seed*.jsonl
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

# The probe fields must be registered before any module that reads the registry at import.
from app import surface_fields  # noqa: E402

SCHEMES = {"band": 1.0, "ever": 1.0, "wide": 2.0}  # band width, in FAILED_RIFT_BAND_MULT spacings
MARKS = ("through", "clear")
PATHS = ("divergent", "stretch")
PROBE_FIELDS = {scheme: {mark: f"probe388_{scheme}_{mark}" for mark in MARKS} for scheme in SCHEMES}
for _names in PROBE_FIELDS.values():
    for _name in _names.values():
        surface_fields.SURFACE_FIELDS[_name] = surface_fields._field(float, 0.0, surface_fields.RemapClass.INTENSIVE)

from app import cratons, merge_split, quad_tectonics  # noqa: E402
from app import world as world_mod  # noqa: E402

KM3 = 1e9


class Probe:
    def __init__(self) -> None:
        # booked[scheme][path][mark] in m^3, mark in ("through", "clear", "never").
        self.booked = {s: {p: {m: 0.0 for m in (*MARKS, "never")} for p in PATHS} for s in SCHEMES}
        self.splits = {s: {m: 0 for m in MARKS} for s in SCHEMES}
        self._pending_cut: tuple | None = None

    def book(self, plate, path: str, lost_m3: np.ndarray) -> None:
        """Book per-node craton loss `lost_m3` on `plate` (current node order) by mark."""
        total = float(lost_m3.sum())
        if total <= 0.0:
            return
        for scheme, names in PROBE_FIELDS.items():
            shares = {mark: float(np.dot(np.clip(plate.collect(name), 0.0, 1.0), lost_m3)) for mark, name in names.items()}
            if scheme == "ever":
                # A sticky through mark can coexist with a later clear mark; through wins.
                shares["clear"] = min(shares["clear"], total - shares["through"])
            for mark, value in shares.items():
                self.booked[scheme][path][mark] += value
            self.booked[scheme][path]["never"] += total - sum(shares.values())

    def note_cut(self, world, plate, points, cut_normal, band_sin, kind) -> None:
        if kind != "split":
            self._pending_cut = None
            return
        distance = np.abs(points @ np.asarray(cut_normal, dtype=float))
        craton = plate.collect("craton_crust_m")
        held = {}
        for scheme, mult in SCHEMES.items():
            sin = float(np.sin(mult * np.arcsin(band_sin)))
            held[scheme] = (sin, bool(np.any(craton[distance < sin] > 0.0)))
        self._pending_cut = (np.asarray(cut_normal, dtype=float).copy(), held)

    def mark_daughters(self, daughters) -> None:
        normal, held = self._pending_cut
        self._pending_cut = None
        for scheme, (sin, through) in held.items():
            mark = "through" if through else "clear"
            self.splits[scheme][mark] += 1
            names = PROBE_FIELDS[scheme]
            for plate in daughters:
                points, _ = plate.all_points_and_elevation()
                in_band = np.abs(points @ normal) < sin
                if not np.any(in_band):
                    continue
                fields = {name: plate.collect(name) for name in names.values()}
                if scheme == "ever":
                    # Sticky: a clear cut only marks the share no through cut has.
                    kept = 1.0 if through else fields[names["through"]][in_band]
                    fields[names["through"]][in_band] = kept
                    fields[names["clear"]][in_band] = 1.0 - kept
                else:
                    fields[names["through"]][in_band] = float(through)
                    fields[names["clear"]][in_band] = float(not through)
                plate.set_fields_on_plate(**fields)


def install(probe: Probe) -> None:
    record_rift_cut = cratons.record_rift_cut
    maybe_split_plate = merge_split.maybe_split_plate
    deform_columns = quad_tectonics.deform_columns
    open_rift = quad_tectonics._open_rift

    def tapped_record_rift_cut(world, plate, points, cut_normal, band_sin, kind):
        probe.note_cut(world, plate, points, cut_normal, band_sin, kind)
        return record_rift_cut(world, plate, points, cut_normal, band_sin, kind)

    def tapped_maybe_split_plate(world, plate):
        probe._pending_cut = None
        result = maybe_split_plate(world, plate)
        if result is not None:
            probe.mark_daughters(result)
        return result

    def tapped_deform_columns(world, plate, ctx, fields, local_xyz, node_area_m2, *args, **kwargs):
        before = fields["craton_crust_m"].copy()
        columns = deform_columns(world, plate, ctx, fields, local_xyz, node_area_m2, *args, **kwargs)
        # The same per-node loss deform_columns books to rifted_divergent_m3.
        probe.book(plate, "divergent", np.where(columns["craton_crust_m"] < before, before - columns["craton_crust_m"], 0.0) * node_area_m2)
        return columns

    def tapped_open_rift(plate, world, *args, **kwargs):
        before = plate.collect("craton_crust_m")
        open_rift(plate, world, *args, **kwargs)
        after = plate.collect("craton_crust_m")
        probe.book(plate, "stretch", np.maximum(before - after, 0.0) * plate.node_areas_m2())

    cratons.record_rift_cut = tapped_record_rift_cut
    merge_split.maybe_split_plate = tapped_maybe_split_plate
    quad_tectonics.deform_columns = tapped_deform_columns
    quad_tectonics._open_rift = tapped_open_rift


def _row(world, probe: Probe, start_craton_m3: float, wall_s: float) -> dict:
    ledger = world.craton_ledger
    return {
        "seed": world.seed,
        "myr": world.elapsed_years / 1e6,
        "wall_s": round(wall_s, 1),
        "start_craton_km3": start_craton_m3 / KM3,
        "craton_km3": cratons.diagnostics(world)["craton_volume_m3"] / KM3,
        "ledger_km3": {key: ledger.get(key, 0.0) / KM3 for key in ("rifted_m3", *cratons.CRATON_RIFT_PATH_ACCOUNTS)},
        "split_cuts": ledger.get("split_cuts", 0.0),
        "splits_marked": probe.splits,
        "booked_km3": {s: {p: {m: v / KM3 for m, v in marks.items()} for p, marks in paths.items()} for s, paths in probe.booked.items()},
    }


def run(seed: int, years: float, every: float, step: float, out: Path) -> None:
    probe = Probe()
    install(probe)
    world = world_mod.generate_world(seed)
    cratons.ensure_initialized(world)
    start = cratons.diagnostics(world)["craton_volume_m3"]
    started = time.time()
    next_checkpoint = every
    with out.open("w") as handle:
        while world.elapsed_years < years - 1.0:
            world_mod.step_world(world, min(step, years - world.elapsed_years))
            if world.elapsed_years >= next_checkpoint - 1.0:
                row = _row(world, probe, start, time.time() - started)
                handle.write(json.dumps(row) + "\n")
                handle.flush()
                print(json.dumps({"seed": seed, "myr": row["myr"], "F_band": _f(row, "band")}), flush=True)
                next_checkpoint += every


def _f(row: dict, scheme: str) -> float:
    rifted = row["ledger_km3"]["rifted_m3"]
    if rifted <= 0.0:
        return 0.0
    booked = row["booked_km3"][scheme]
    return (row["ledger_km3"]["rifted_failed_m3"] + sum(booked[p]["through"] for p in PATHS)) / rifted


def summarize(paths: list[Path]) -> dict:
    """The last checkpoint of each file: F per scheme, the marks' shares of divergent and
    stretch loss, and how well the per-node booking matches the ledger's path accounts."""
    result = {}
    for path in paths:
        rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
        row = rows[-1]
        ledger = row["ledger_km3"]
        rifted = ledger["rifted_m3"]
        band = row["booked_km3"]["band"]
        boundary = sum(sum(band[p].values()) for p in PATHS)
        entry = {
            "myr": row["myr"],
            "rifted_share_of_start_craton": round(rifted / row["start_craton_km3"], 4),
            "rifted_km3": round(rifted),
            "path_share_of_rifted": {
                "divergent": round(ledger["rifted_divergent_m3"] / rifted, 3),
                "stretch": round(ledger["rifted_stretch_m3"] / rifted, 3),
                "failed": round(ledger["rifted_failed_m3"] / rifted, 3),
            },
            "booking_vs_ledger": {
                p: round(sum(band[p].values()) / ledger[f"rifted_{p}_m3"], 6) if ledger[f"rifted_{p}_m3"] else None for p in PATHS
            },
            "boundary_loss_share_by_mark": {
                scheme: {m: round(sum(row["booked_km3"][scheme][p][m] for p in PATHS) / boundary, 3) for m in (*MARKS, "never")}
                for scheme in SCHEMES
            } if boundary else {},
            "boundary_loss_by_path_and_mark": {
                p: {m: round(band[p][m] / ledger[f"rifted_{p}_m3"], 3) if ledger[f"rifted_{p}_m3"] else None for m in (*MARKS, "never")}
                for p in PATHS
            },
            "splits_marked": row["splits_marked"]["band"],
            "splits_marked_wide": row["splits_marked"]["wide"],
            "F": {scheme: round(_f(row, scheme), 3) for scheme in SCHEMES},
        }
        result[str(row["seed"])] = entry
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--years", type=float, default=1e9)
    parser.add_argument("--every", type=float, default=100e6)
    parser.add_argument("--step", type=float, default=5e6, help="step_world years per call (compare_cratons.py's default)")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--summarize", type=Path, nargs="+", help="JSONL files to summarize instead of running")
    args = parser.parse_args()
    if args.summarize:
        print(json.dumps(summarize(args.summarize), indent=1))
        return
    if args.seed is None or args.out is None:
        parser.error("--seed and --out are required to run")
    run(args.seed, args.years, args.every, args.step, args.out)


if __name__ == "__main__":
    main()
