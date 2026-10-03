#!/usr/bin/env python3
"""Long-run, multi-seed comparison of worlds with and without cratons (GitHub issue #274).

Each (seed, variant) job generates a default quad world and steps it to --years, recording a
checkpoint every --every years: genuine continental area (Hc >= cratons.CRATON_HOST_MIN_HC_M),
node land fraction, the largest plate's genuine continental area, craton extent and ledger,
and how much tectonics is still happening (plate count, merges/splits/rift events so far, and
mean plate rotation rate). The "off" variant disables craton seeding and formation, which
leaves every craton resistance inert (strength 0 everywhere) -- the pre-#274 physics.

Usage:
    backend/.venv/bin/python bin/debug/compare_cratons.py --seeds 7,11,2024 --years 500e6 \
        --out /tmp/cratons.jsonl --workers 6
    backend/.venv/bin/python bin/debug/compare_cratons.py --summarize /tmp/cratons.jsonl
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

import numpy as np  # noqa: E402

from app import cratons  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad  # noqa: E402

VARIANTS = ("off", "on")
EVENT_KEYS = {"merged": ("merged", "fused"), "split": ("split into",), "failed_rift": ("rift failed",), "craton": ("Craton destruction",)}


def _snapshot(world, events_seen: dict[str, int]) -> dict:
    spacing = line_spacing_rad(world.node_density)
    genuine_area = 0.0
    largest = 0.0
    land_nodes = 0
    nodes = 0
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        hc = plate.collect("crustal_thickness_m")
        areas = plate.accounting_areas_m2(spacing)
        continental = effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        genuine = float(areas[continental & (hc >= cratons.CRATON_HOST_MIN_HC_M)].sum())
        genuine_area += genuine
        largest = max(largest, genuine)
        land_nodes += int(np.count_nonzero(plate.collect("elevation") > world.sea_level_m))
        nodes += plate.node_count()
    craton = cratons.diagnostics(world)
    return {
        "myr": world.elapsed_years / 1e6,
        "genuine_continental_area_km2": genuine_area / 1e6,
        "largest_plate_continental_area_km2": largest / 1e6,
        "land_fraction_node": land_nodes / max(nodes, 1),
        "plate_count": len(world.plates),
        "mean_rotation_rate": float(np.mean([np.linalg.norm(p.omega) for p in world.plates])),
        "craton_area_km2": craton["craton_area_m2"] / 1e6,
        "craton_continental_fraction": craton["craton_area_fraction_of_continental"],
        "craton_volume_km3": craton["craton_volume_m3"] / 1e9,
        "craton_balance_error_km3": craton["craton_balance_error_m3"] / 1e9,
        "craton_ledger_km3": {key: value / 1e9 for key, value in world.craton_ledger.items()},
        "events": dict(events_seen),
    }


_DEFAULTS = {name: getattr(cratons, name) for name in dir(cratons) if name.startswith("CRATON_") and isinstance(getattr(cratons, name), float)}


def run_job(seed: int, variant: str, years: float, every: float, step: float, overrides: dict[str, float]) -> list[dict]:
    # Pool workers are reused across jobs: always reset, then apply overrides / disable for "off".
    for name, value in _DEFAULTS.items():
        setattr(cratons, name, value)
    for name, value in overrides.items():
        setattr(cratons, name, value)
    if variant == "off":
        cratons.CRATON_SEED_MARGIN_KM = 1e6  # wider than any plate: no interior qualifies
        cratons.CRATON_FORMATION_MYR = float("inf")
    world = world_mod.generate_world(seed)
    events_seen: dict[str, int] = defaultdict(int)
    rows = [{"seed": seed, "variant": variant, **_snapshot(world, events_seen)}]
    started = time.time()
    next_checkpoint = every
    while world.elapsed_years < years - 1.0:
        previous = world.elapsed_years
        world_mod.step_world(world, min(step, years - world.elapsed_years))
        # Topology and craton events are logged after the step's clock advances; the log is
        # capped, so select by timestamp rather than position.
        for logged_at, message in world.events:
            if logged_at > previous:
                for key, needles in EVENT_KEYS.items():
                    if any(needle in message for needle in needles):
                        events_seen[key] += 1
        if world.elapsed_years >= next_checkpoint - 1.0:
            row = {"seed": seed, "variant": variant, **_snapshot(world, events_seen), "wall_s": time.time() - started}
            rows.append(row)
            print(json.dumps({k: row[k] for k in ("seed", "variant", "myr", "genuine_continental_area_km2", "craton_continental_fraction")}), flush=True)
            next_checkpoint += every
    return rows


def summarize(path: Path) -> None:
    rows = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    by = defaultdict(dict)
    for row in rows:
        by[(row["seed"], row["variant"])][round(row["myr"])] = row
    seeds = sorted({seed for seed, _ in by})
    checkpoints = sorted({myr for runs in by.values() for myr in runs})
    print("Genuine continental area (10^6 km^2) / largest plate's / land fraction (node) / craton share of continental")
    for seed in seeds:
        print(f"seed {seed}")
        for myr in checkpoints:
            cells = []
            for variant in VARIANTS:
                row = by.get((seed, variant), {}).get(myr)
                if row is None:
                    cells.append(f"{variant}: --")
                    continue
                cells.append(
                    f"{variant}: {row['genuine_continental_area_km2'] / 1e6:6.1f} / {row['largest_plate_continental_area_km2'] / 1e6:5.1f}"
                    f" / {row['land_fraction_node']:.3f} / {row['craton_continental_fraction']:.2f}"
                )
            print(f"  {myr:4d} Myr  " + "   ".join(cells))
    print("\nFinal craton ledger (10^6 km^3) and topology events")
    for (seed, variant), runs in sorted(by.items()):
        last = runs[max(runs)]
        ledger = {k.removesuffix("_m3"): round(v / 1e6, 2) for k, v in last["craton_ledger_km3"].items() if v}
        print(f"  {seed} {variant}: {ledger} balance={last['craton_balance_error_km3']:.3g} km^3 events={last['events']}"
              f" plates={last['plate_count']} mean_omega={last['mean_rotation_rate']:.3g}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default="7,11,2024")
    parser.add_argument("--variants", default=",".join(VARIANTS))
    parser.add_argument("--years", type=float, default=500e6)
    parser.add_argument("--every", type=float, default=50e6)
    parser.add_argument("--step", type=float, default=5e6)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--set", action="append", default=[], metavar="CRATON_NAME=VALUE", help="override a cratons.py constant")
    parser.add_argument("--out", type=Path)
    parser.add_argument("--summarize", type=Path)
    args = parser.parse_args()
    if args.summarize:
        summarize(args.summarize)
        return
    if args.out is None:
        parser.error("--out is required unless --summarize is given")
    overrides = {name: float(value) for name, value in (item.split("=", 1) for item in args.set)}
    unknown = sorted(set(overrides) - set(_DEFAULTS))
    if unknown:
        parser.error(f"unknown cratons constant(s): {unknown}")
    jobs = [(int(seed), variant) for seed in args.seeds.split(",") for variant in args.variants.split(",")]
    with ProcessPoolExecutor(max_workers=args.workers) as pool, args.out.open("a") as out:
        futures = {pool.submit(run_job, seed, variant, args.years, args.every, args.step, overrides): (seed, variant) for seed, variant in jobs}
        for future in as_completed(futures):
            for row in future.result():
                out.write(json.dumps(row) + "\n")
            out.flush()
    summarize(args.out)


if __name__ == "__main__":
    main()
