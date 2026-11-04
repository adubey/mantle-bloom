#!/usr/bin/env python3
"""Capture identical pre-split G0 worlds for isolated RFC 0018 explorations.

Capture each preregistered development world at the first age-qualified split decision,
before the production splitter mutates it. The production split function still runs after
the snapshot is saved; these files are read-only inputs to separate diagnostic branches.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python -u \
        bin/debug/capture_rfc0018_g0_snapshots.py \
        --out-dir ~/mantle-bloom-captures/rfc0018-g0-snapshots
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import time

from app import merge_split, persistence, world as world_mod

PRIMARY_SEEDS = [1629738962, 1097710202, 591453549, 969845080, 763144949, 870996356,
                 494568917, 1854897345, 954941394, 600370317, 1643146570, 146046339,
                 1397817354, 1630166320, 604997962, 704718994, 1097624328, 487506346,
                 106533960, 1439649077]
REPLACEMENTS = [1227913777, 303885793, 2045331463, 157959777, 1208093067]
CAPTURE_YEAR = 2_000_000.0
STEP_YEARS = 1_000_000.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir = args.out_dir.expanduser()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.out_dir / "manifest.json"
    manifest = {
        "protocol": "RFC0018-G0-snapshot-v1",
        "source_commit": "63290b053de15b72280fe28fe0086ff37affdba2",
        "capture_year": CAPTURE_YEAR,
        "step_years": STEP_YEARS,
        "world_defaults": True,
        "primary_seeds": PRIMARY_SEEDS,
        "replacement_seeds": REPLACEMENTS,
        "worlds": [],
    }
    original = merge_split.maybe_split_plate
    for seed in PRIMARY_SEEDS + REPLACEMENTS:
        if sum(bool(row.get("captured")) for row in manifest["worlds"]) >= 20:
            break
        row = {"seed": seed, "captured": False}
        started = time.perf_counter()
        try:
            world = world_mod.generate_world(seed)
            if world is None or not world.plates:
                raise RuntimeError("generation returned no world or no plates")
        except Exception as exc:
            row["generation_error"] = f"{type(exc).__name__}: {exc}"
            row["runtime_s"] = time.perf_counter() - started
            manifest["worlds"].append(row)
            manifest_path.write_text(json.dumps(manifest, separators=(",", ":")) + "\n")
            continue

        captured = False

        def capture_before_production_split(current_world, plate, step_years=0.0):
            nonlocal captured
            if not captured and current_world.elapsed_years >= CAPTURE_YEAR:
                filename = f"g0-snapshot-seed-{seed}.mbworld"
                data = persistence.save_world_bytes(
                    current_world,
                    name=f"RFC0018 G0 frozen snapshot seed {seed}",
                )
                path = args.out_dir / filename
                path.write_bytes(data)
                loaded = persistence.load_world_bytes(data)
                if loaded.seed != seed or loaded.elapsed_years != CAPTURE_YEAR:
                    raise RuntimeError("snapshot failed seed/time validation")
                row.update({
                    "captured": True,
                    "file": filename,
                    "sha256": hashlib.sha256(data).hexdigest(),
                    "plate_count": len(loaded.plates),
                    "opportunity_step_elapsed_myr": current_world.elapsed_years / 1e6,
                })
                captured = True
                manifest["worlds"].append(row)
                manifest_path.write_text(json.dumps(manifest, separators=(",", ":")) + "\n")
                print(json.dumps({"seed": seed, "captured": True,
                                  "plates": row["plate_count"], "sha256": row["sha256"]}), flush=True)
            return original(current_world, plate, step_years)

        merge_split.maybe_split_plate = capture_before_production_split
        try:
            while world.elapsed_years < CAPTURE_YEAR:
                world_mod.step_world(world, STEP_YEARS)
        finally:
            merge_split.maybe_split_plate = original
        if not captured:
            row["runtime_s"] = time.perf_counter() - started
            row["capture_error"] = "no split-decision callback occurred by the capture year"
            manifest["worlds"].append(row)
            manifest_path.write_text(json.dumps(manifest, separators=(",", ":")) + "\n")
        else:
            row["runtime_s"] = time.perf_counter() - started
            manifest_path.write_text(json.dumps(manifest, separators=(",", ":")) + "\n")

    complete = sum(bool(row.get("captured")) for row in manifest["worlds"])
    sums = "".join(
        f"{row['sha256']}  {row['file']}\n"
        for row in manifest["worlds"] if row.get("captured")
    )
    (args.out_dir / "SHA256SUMS").write_text(sums)
    print(json.dumps({"captured": complete, "world_rows": len(manifest["worlds"]),
                      "manifest": str(manifest_path)}, indent=2))


if __name__ == "__main__":
    main()
