#!/usr/bin/env python3
"""Issue #747: give the starting rock in an existing save the dates new worlds get.

Saves from before #744 hold the -1.0 "predates tracking" sentinel in `node_created_years` on
every starting cell. This writes a copy of a save in which each such cell takes a negative
starting date from `initial_rock_ages.assign_initial_rock_ages`: ocean ages away from the
save's ridges and continental crust away from its cratons, with area-weighted rank mapped
onto the Earth-like age distributions. Cells with a real creation date keep it and take no
part in the ranking, and nothing else in the world changes.

This approximates what generation would have assigned, because ridges are read from the
plate motions stored at the save's time, not those of t=0. A cell that blended starting rock
with dated rock already reads as the dated rock, and that can't be recovered.
`World.initial_rock_ages` keys on the save's plate ids, and each record carries
`dated_at_years` to tell it from a record written at generation.

The copy is written in the exact profile, which keeps everything any other profile keeps,
so it continues just as the source does. The source is never overwritten. A save with no
undated cells writes nothing, so a second run changes nothing. A line-backed save converts
on load like any other, and the copy is the converted quad world.

Usage (from repo root, this repo's venv):
    backend/.venv/bin/python bin/debug/date_saved_world.py SAVE.mbworld \\
        [--out DATED.mbworld] [--force]

--out defaults to <save>-dated.mbworld beside the source. --force replaces an existing --out.
Prints a JSON summary of what changed.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "backend"))

from app import elevation_lines, initial_rock_ages, npz_save, persistence  # noqa: E402
from app.surface_fields import SURFACE_FIELDS  # noqa: E402
from app.world import World  # noqa: E402

UNDATED = SURFACE_FIELDS["node_created_years"].sentinel


def _undated_cells(world: World) -> int:
    return sum(int(np.count_nonzero(plate.collect("node_created_years") == UNDATED)) for plate in world.plates)


def date_world(world: World) -> dict:
    """Date every undated cell of `world` in place and return a summary of the change."""
    undated = _undated_cells(world)
    omegas = {plate.plate_id: plate.omega for plate in world.plates}
    records = initial_rock_ages.assign_initial_rock_ages(world, omegas, undated_only=True)
    world.initial_rock_ages = {
        **world.initial_rock_ages,
        **{plate_id: {**record, "dated_at_years": world.elapsed_years} for plate_id, record in records.items()},
    }
    return {"elapsed_myr": world.elapsed_years / 1e6, "plates": len(world.plates), "plates_dated": len(records),
            "cells": sum(plate.node_count() for plate in world.plates), "cells_dated": undated,
            **_starting_rock_ages(world)}


def _starting_rock_ages(world: World) -> dict:
    """Mean age now, in Myr, of the rock dated before t=0, by crust type."""
    ages = {True: [], False: []}
    for plate in world.plates:
        created = plate.collect("node_created_years")
        continental = elevation_lines.effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental")
        for is_continental, rows in ages.items():
            rows.append(created[(created < 0.0) & (continental == is_continental)])
    summary = {}
    for is_continental, rows in ages.items():
        created = np.concatenate(rows) if rows else np.zeros(0)
        kind = "continental" if is_continental else "oceanic"
        summary[f"{kind}_mean_age_myr"] = (
            round(float((world.elapsed_years - created.mean()) / 1e6), 1) if len(created) else None)
    return summary


def _write(data: bytes, out: Path) -> None:
    # Through a sibling temporary file, so a failed write can't leave a truncated save.
    partial = out.with_name(f"{out.name}.partial")
    partial.write_bytes(data)
    os.replace(partial, out)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("save", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    out = args.out or args.save.with_name(f"{args.save.stem}-dated{args.save.suffix}")
    if out.resolve() == args.save.resolve():
        parser.error("--out must differ from the source save")
    if out.exists() and not args.force:
        parser.error(f"{out} exists; pass --force to replace it")

    data = args.save.read_bytes()
    world = persistence.load_world_bytes(data)
    if _undated_cells(world) == 0:
        print(json.dumps({"written": None, "reason": "no undated cells"}))
        return 0
    summary = date_world(world)
    _write(npz_save.save_world_v5(world, profile="exact", name=persistence.load_world_name(data)), out)
    print(json.dumps({"written": str(out), **summary}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
