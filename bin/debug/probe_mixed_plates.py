#!/usr/bin/env python3
"""RFC 0009 gate G2: do mixed plates work? (#438)

The cratons-first generator makes most plates mixed (continental and oceanic crust) from
t = 0, where today's generator makes every plate wholly one type. This probe asks whether the
runtime behaves the same on such plates, before that generator is built.

Per seed it runs three worlds over the same starting crust, to --years with 5 Myr steps:
- **control**: `generate_world(seed)` unchanged (today's single-type plates);
- **mixed A** and **mixed B**: the same finished world, re-tiled. A new `build_plate_tiling`
  (same plate count and extra sites) is drawn from its own RNG stream, keyed on
  (seed, TILING_TAG, run), and each new plate's cells copy every surface field from the
  nearest control node: Hc, Hm, elevation, crust type, the seeded cratons and the ledger
  tracer. A plate's nominal `crust_type` is the type of the majority of its area and its
  minority-type cells carry an explicit `crust_type_code`, as RFC 0009's generator would
  build them. The re-tiled world then runs `finish_generation` like any new world, with the
  copied craton volume booked as `initial_m3` so no second seeding runs.

At t = 0 and each --every checkpoint it records genuine continental area (`compare_cratons`'
measure) and plate overlap (`stats`), both as a share of the sphere, active plates, the share
of plates mixed by continental-crust area (5-95%), cumulative merges and splits (logged
events, as `compare_cratons` counts them) and the continental-material and craton ledger
balance errors. To diagnose a failure it also records each continental-material ledger
account and `measure_plate_layout`'s share of the continent's edge within 300 km of a plate
boundary (the margins boundary sinks can reach). `--summarize` applies the RFC's decision rule to a finished result:
noise floor = max(|mixed A - mixed B|, floor), floor 1 plate, 1 event or 1 percentage point;
a metric fails when |mean(mixed) - control| > 2 x noise floor on at least two seeds.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_mixed_plates.py \
        --out /tmp/g2.json
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_mixed_plates.py \
        --seeds 7 --years 1e7 --every 5e6 --out /tmp/g2-smoke.json
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_mixed_plates.py \
        --summarize /tmp/g2.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from collections import defaultdict
from pathlib import Path

import compare_cratons
import measure_plate_layout
import numpy as np
from scipy.spatial import cKDTree

from app import continental_ledger, cratons, geometry, lithosphere_plate, quad_merge, stats
from app import world as world_mod
from app.elevation_lines import CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC, effective_is_continental_from_codes, line_spacing_rad
from app.sparse_quad_patch import PlateWithSparseQuadPatch
from app.surface_fields import SURFACE_FIELDS

TILING_TAG = 438
RUNS = ("control", "mixed_a", "mixed_b")
MIXED_SHARE = (0.05, 0.95)
# The RFC's floors, in each metric's units (percent of the sphere for areas).
FLOORS = {"genuine_pct": 1.0, "overlap_pct": 1.0, "active_plates": 1.0, "merges": 1.0, "splits": 1.0}


def _continental(plate) -> np.ndarray:
    return effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")


def retile(world, run: int):
    """A new World over `world`'s crust on an independent tiling -- see the module docstring."""
    spacing = line_spacing_rad(world.node_density)
    live = [p for p in world.plates if p.node_count() > 0]
    tree = cKDTree(np.concatenate([p.all_points_and_elevation()[0] for p in live]))
    fields = {name: np.concatenate([p.collect(name) for p in live]) for name in SURFACE_FIELDS}
    fields["crust_type_code"] = np.where(np.concatenate([_continental(p) for p in live]), CRUST_TYPE_CONTINENTAL, CRUST_TYPE_OCEANIC).astype(
        fields["crust_type_code"].dtype
    )
    tiling = lithosphere_plate.build_plate_tiling(np.random.default_rng((world.seed, TILING_TAG, run)), len(live))
    owner = cKDTree(tiling.site_xyz)
    plates = []
    for i in range(tiling.num_plates):
        frame = geometry.plate_frame_from_seed(tiling.primary_site(i))

        def is_owned(pts: np.ndarray, _i: int = i) -> np.ndarray:
            return tiling.site_plate[owner.query(pts)[1]] == _i

        cells = PlateWithSparseQuadPatch.from_lattice(i, frame, "oceanic", spacing, is_owned)
        source = tree.query(cells.all_points_and_elevation()[0])[1]
        explicit = fields["crust_type_code"][source]
        areas = cells.node_areas_m2()
        majority = "continental" if areas[explicit == CRUST_TYPE_CONTINENTAL].sum() > areas[explicit != CRUST_TYPE_CONTINENTAL].sum() else "oceanic"
        own = {name: values[source] for name, values in fields.items()}
        own["crust_type_code"] = quad_merge._relative_codes(explicit, majority)
        plates.append(PlateWithSparseQuadPatch(i, frame, majority, cells.cells_per_edge, cells.cell_keys, fields=own))
    mixed = world_mod.World(
        seed=world.seed,
        plates=plates,
        mantle_centers=world.mantle_centers,
        next_plate_id=len(plates),
        axial_tilt_deg=world.axial_tilt_deg,
        node_density=world.node_density,
        climate_density=world.climate_density,
        fluid_density=world.fluid_density,
    )
    # The copied cratons are this world's seeded cratons: book them so ensure_initialized
    # doesn't seed again from the new plates' interiors.
    cratons.ensure_ledger(mixed)
    mixed.craton_ledger["initial_m3"] = cratons.live_volume_m3(mixed)
    n_continental = sum(1 for p in plates if p.crust_type == "continental")
    world_mod.finish_generation(mixed, f"World re-tiled with {len(plates)} mixed plates ({n_continental} majority continental).")
    return mixed


def _snapshot(world, events_seen: dict[str, int]) -> dict:
    base = compare_cratons._snapshot(world, events_seen)
    live = [p for p in world.plates if p.node_count() > 0]
    mixed = 0
    for plate in live:
        areas = plate.node_areas_m2()
        share = float(areas[_continental(plate)].sum() / areas.sum())
        mixed += MIXED_SHARE[0] <= share <= MIXED_SHARE[1]
    sphere_km2 = stats.SPHERE_AREA_M2 / 1e6
    return {
        "myr": base["myr"],
        "genuine_pct": 100.0 * base["genuine_continental_area_km2"] / sphere_km2,
        "overlap_pct": 100.0 * stats._overlap_once_totals(world)["plate_overlap_area_fraction"],
        "active_plates": len(live),
        "mixed_plate_share": mixed / max(len(live), 1),
        "merges": events_seen.get("merged", 0),
        "splits": events_seen.get("split", 0),
        "failed_rifts": events_seen.get("failed_rift", 0),
        "craton_mkm2": base["craton_area_km2"] / 1e6,
        "continental_material_km3": continental_ledger.surface_volume_m3(world) / 1e9,
        "continental_balance_error_rel": continental_ledger.balance_error_m3(world)
        / max(world.continental_material_ledger["initial_continental_m3"], 1.0),
        "craton_balance_error_rel": cratons.balance_error_m3(world) / max(world.craton_ledger.get("initial_m3", 0.0), 1.0),
        "continent_edge_on_boundary_share": measure_plate_layout._measure(measure_plate_layout._grid(world), "")["continent_edge_on_boundary_share"],
        "continental_ledger_km3": {key: value / 1e9 for key, value in world.continental_material_ledger.items()},
    }


def _run(seed: int, run: str, years: float, every: float, step: float) -> dict:
    started = time.perf_counter()
    world = world_mod.generate_world(seed)
    if run != "control":
        world = retile(world, RUNS.index(run))
    generation_s = time.perf_counter() - started
    events_seen: dict[str, int] = defaultdict(int)
    checkpoints = [_snapshot(world, events_seen)]
    next_checkpoint = every
    while world.elapsed_years < years - 1.0:
        previous = world.elapsed_years
        world_mod.step_world(world, min(step, years - world.elapsed_years))
        for logged_at, message in world.events:
            if logged_at > previous:
                for key, needles in compare_cratons.EVENT_KEYS.items():
                    if any(needle in message for needle in needles):
                        events_seen[key] += 1
        if world.elapsed_years >= next_checkpoint - 1.0:
            checkpoints.append({**_snapshot(world, events_seen), "wall_s": round(time.perf_counter() - started, 1)})
            print(json.dumps({"seed": seed, "run": run, **checkpoints[-1]}), flush=True)
            next_checkpoint += every
    return {"seed": seed, "run": run, "generation_s": generation_s, "checkpoints": checkpoints}


def _metrics(run: dict) -> dict[str, float]:
    """The gate's metrics: genuine area and cumulative events at the end, plate count and
    overlap at every checkpoint after t = 0."""
    last = run["checkpoints"][-1]
    out = {"genuine_pct@end": last["genuine_pct"], "merges@end": last["merges"], "splits@end": last["splits"]}
    for c in run["checkpoints"][1:]:
        for key in ("active_plates", "overlap_pct"):
            out[f"{key}@{round(c['myr'])}"] = c[key]
    return out


def summarize(result: dict) -> dict:
    """Apply the RFC's decision rule per metric (see the module docstring)."""
    by_seed: dict[int, dict[str, dict]] = defaultdict(dict)
    for run in result["runs"]:
        by_seed[run["seed"]][run["run"]] = _metrics(run)
    rows: dict[str, dict] = {}
    for seed, runs in sorted(by_seed.items()):
        if set(runs) != set(RUNS):
            continue
        for metric, control in runs["control"].items():
            a, b = runs["mixed_a"][metric], runs["mixed_b"][metric]
            noise = max(abs(a - b), FLOORS[metric.split("@")[0]])
            delta = abs((a + b) / 2 - control)
            row = rows.setdefault(metric, {"seeds": {}, "failing_seeds": 0})
            row["seeds"][seed] = {"control": control, "mixed_a": a, "mixed_b": b, "delta": delta, "noise": noise, "fails": delta > 2 * noise}
            row["failing_seeds"] += delta > 2 * noise
    for row in rows.values():
        row["fails"] = row["failing_seeds"] >= 2
    ledgers = max(
        max(abs(c["continental_balance_error_rel"]), abs(c["craton_balance_error_rel"])) for run in result["runs"] for c in run["checkpoints"]
    )
    return {"metrics": rows, "failing": sorted(m for m, r in rows.items() if r["fails"]), "max_ledger_error_rel": ledgers}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default="7,11,2024")
    parser.add_argument("--runs", default=",".join(RUNS))
    parser.add_argument("--years", type=float, default=500e6)
    parser.add_argument("--every", type=float, default=100e6)
    parser.add_argument("--step", type=float, default=5e6)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--summarize", type=Path, help="print the decision rule for a finished --out file and exit")
    args = parser.parse_args()
    if args.summarize:
        print(json.dumps(summarize(json.loads(args.summarize.read_text())), indent=1))
        return
    if args.out is None:
        parser.error("--out is required unless --summarize is given")
    app_root = Path(world_mod.__file__).resolve().parents[2]
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=app_root, text=True).strip()
    result = {"model_commit": commit, "years": args.years, "step": args.step, "every": args.every, "runs": []}
    for seed in map(int, args.seeds.split(",")):
        for run in args.runs.split(","):
            result["runs"].append(_run(seed, run, args.years, args.every, args.step))
            result["summary"] = summarize(result)
            args.out.parent.mkdir(parents=True, exist_ok=True)
            args.out.write_text(json.dumps(result, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
