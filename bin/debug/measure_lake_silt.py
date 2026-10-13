#!/usr/bin/env python3
"""Issue #363: replay a saved world and measure its lakes and lake silt, on `main` or on the branch.

Replays `--myr` from the save in `--step-years` steps and records, at every checkpoint:

- land fraction (by area), lake count (connected bodies of visible lake water), lake area
  (km^2) and lake water volume (km^3);
- silt laid on lake beds, summed over the steps since the last checkpoint (km^3): `applied`
  is `HydrologyFields.silt_deposited` over cell area, after each lake's room clamp;
- the same window's silt *budget* before that clamp, in both units, from a tap on
  `lakes._water_balance`: `legacy` is the model-unit rule main used before #363 (0.0012, 1/100
  of the retired lake fill rate) times the lake's `water_deposited` (stacked mm/yr) times step
  Myr times the reference cell area, reported only where `_water_balance` still takes
  `water_deposited`; `real_per_unit` is the lake's real inflow
  (`LakeClimate.inflow_m3_per_yr`) times step years, i.e. the silt at a coefficient of 1.
  Their ratio calibrates the real-unit coefficient. Each is split by whether the lake ended
  the step spilling;
- step runtime and peak RSS.

`--backend` points at the backend to import, so the same script measures both trees:

    cd backend
    .venv/bin/python ../bin/debug/measure_lake_silt.py ~/Downloads/<save>.mbworld --myr 10 --out <scratch>/branch.json
    .venv/bin/python ../bin/debug/measure_lake_silt.py ~/Downloads/<save>.mbworld --myr 10 \\
        --backend ../../mantle-bloom-origin-main/backend --out <scratch>/main.json
"""

from __future__ import annotations

import argparse
import hashlib
import inspect
import json
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

# The pre-#363 SILT_ACCUMULATION_COEFFICIENT (1/100 of the retired lake fill rate 0.12).
LEGACY_SILT_COEFFICIENT = 0.12 / 100.0


def main() -> None:  # noqa: C901
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("save", type=Path)
    ap.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
    ap.add_argument("--myr", type=float, required=True)
    ap.add_argument("--step-years", type=float, default=100_000.0)
    ap.add_argument("--checkpoint-myr", type=float, default=1.0)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    backend = args.backend.resolve()
    sys.path.insert(0, str(backend))
    from app import persistence, plates as plates_mod, world as world_mod  # noqa: E402
    from app.elevation_lines import line_spacing_rad  # noqa: E402
    from app.hydroclimate import hydrology, lakes  # noqa: E402

    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=backend, check=True, capture_output=True, text=True).stdout.strip()
    data = args.save.read_bytes()
    world = persistence.load_world_bytes(data)
    result = {
        "issue": 363,
        "backend": str(backend),
        "commit": commit,
        "save": str(args.save),
        "save_sha256": hashlib.sha256(data).hexdigest(),
        "step_years": args.step_years,
        "node_density": world.node_density,
        "silt_coefficient": lakes.SILT_ACCUMULATION_COEFFICIENT,
        "checkpoints": [],
    }
    del data

    def empty_window() -> dict:
        budget = {f"{unit}_{kind}": 0.0 for unit in ("legacy_km3", "real_per_unit_km3") for kind in ("closed", "spilling")}
        return {"seconds": 0.0, "applied_km3": 0.0, **budget}

    window = empty_window()
    original_balance = lakes._water_balance
    balance_signature = inspect.signature(original_balance)

    def tapped_balance(*a, **k):
        level = original_balance(*a, **k)
        bound = balance_signature.bind(*a, **k)
        bound.apply_defaults()
        p = bound.arguments
        lake, members = p["lake"], p["lake"].members
        hyps = p["hypsometry"] or lakes.CellHypsometry.whole_cells(p["elevation"])
        # A frozen lake returns its dry level but still silts at its balanced level, so redo
        # the balance (same cap as `_water_balance`) to see whether it holds water.
        dry = hyps.dry_level(members)
        rim = (lake.max_depth - lake.floor_elevation) if lake.max_depth is not None else None
        cap_depth = p["tier_max_depth"] if rim is None else min(rim, p["tier_max_depth"])
        cap = dry if lake.drains_at_floor else max(lake.floor_elevation + cap_depth, dry)
        wet_level = lakes.balanced_level(lake, hyps, p["climate"], dry, cap)
        if not (hyps.wet_fraction(wet_level, members) > 0.0).any():
            return level
        kind = "spilling" if lake.max_depth is not None and wet_level >= lake.max_depth - 1e-6 else "closed"
        if "water_deposited" in p:
            legacy = LEGACY_SILT_COEFFICIENT * float(p["water_deposited"][members].sum()) * p["years_myr"] * hyps.reference_area_m2
            window[f"legacy_km3_{kind}"] += legacy / 1e9
        if p["climate"] is not None:
            years = p["years_myr"] * 1e6
            window[f"real_per_unit_km3_{kind}"] += float(p["climate"].inflow_m3_per_yr[members].sum()) * years / 1e9
        return level

    lakes._water_balance = tapped_balance

    def lake_bodies(hydro, visible: np.ndarray) -> int:
        """Connected components of visibly wet land nodes."""
        seen = np.zeros(len(visible), dtype=bool)
        count = 0
        for start in np.flatnonzero(visible):
            if seen[start]:
                continue
            count += 1
            seen[start] = True
            stack = [start]
            while stack:
                node = stack.pop()
                for nb in hydro.neighbor_idx[node]:
                    if nb >= 0 and visible[nb] and not seen[nb]:
                        seen[nb] = True
                        stack.append(nb)
        return count

    def snapshot() -> dict:
        hydro = world.hydrology_cache
        if hydro is None:
            return {"elapsed_myr": world.elapsed_years / 1e6}
        areas = plates_mod.collect_all_accounting_areas_m2(hydro.plates_in_order, line_spacing_rad(world.node_density))
        land = ~hydro.is_ocean
        visible = land & (hydro.lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M)
        return {
            "elapsed_myr": world.elapsed_years / 1e6,
            "land_pct": 100.0 * float(areas[land].sum() / areas.sum()),
            "lake_count": lake_bodies(hydro, visible),
            "lake_area_km2": float(areas[visible].sum()) / 1e6,
            "sea_tier_area_km2": float(areas[visible & hydro.is_sea].sum()) / 1e6,
            "lake_volume_km3": float(np.dot(np.where(land, hydro.lake_depth, 0.0), areas)) / 1e9,
        }

    steps = int(round(args.myr * 1e6 / args.step_years))
    per_checkpoint = max(1, int(round(args.checkpoint_myr * 1e6 / args.step_years)))
    result["checkpoints"].append(snapshot())
    print(json.dumps(result["checkpoints"][-1]), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    for i in range(1, steps + 1):
        started = time.perf_counter()
        world_mod.step_world(world, args.step_years)
        window["seconds"] += time.perf_counter() - started
        hydro = world.hydrology_cache
        if hydro is not None:
            areas = plates_mod.collect_all_accounting_areas_m2(hydro.plates_in_order, line_spacing_rad(world.node_density))
            window["applied_km3"] += float(np.dot(hydro.silt_deposited, areas)) / 1e9
        if i % per_checkpoint == 0 or i == steps:
            snap = snapshot()
            seconds = window.pop("seconds")
            snap["silt"] = {k: round(v, 6) for k, v in window.items()}
            snap["step_seconds_mean"] = seconds / per_checkpoint
            snap["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
            result["checkpoints"].append(snap)
            print(json.dumps(snap), flush=True)
            args.out.write_text(json.dumps(result, indent=1))
            window = empty_window()
    args.out.write_text(json.dumps(result, indent=1))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
