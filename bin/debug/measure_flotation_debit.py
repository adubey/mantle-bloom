#!/usr/bin/env python3
"""RFC 0015 gate G4 (#912): ice above flotation and the eustatic debit.

Replays a save with `step_world` and, after each step, on the step's final state:

- classifies every ice cell (`glacier_depth > 0`) at the replay's own sea level as **land**
  (outside the connected-ocean mask, so water depth `D` = 0), **grounded marine** (in the
  mask, `0.917 H > D`) or **floating** (in the mask, `0.917 H <= D`), with each class's
  volume and area;
- solves three sea levels against the same hypsometry and water budget:
  - **T, today**: debit `0.917 H` on every cell (re-solved here; the replay's own level is
    recorded beside it);
  - **S, D6's simple form**: debit `max(0, 0.917 H - D)` with `D` taken at the level known
    when `update_sea_level` runs, the previous step's closing level. The arm carries its own
    level from step to step, starting from the save's;
  - **G, generalised** (Goelzer et al., 2020, Eqs. 5, 13, 14 with the reference level set to
    the level being solved): the liquid-water balance
    `sum_mask(h) A * max(0, (h - z) - 0.917 H) = V_total - L - sum A * 0.917 H`;
- reports `h_T - h_S` (today against D6), `h_S - h_G` (the generalised correction, the gate's
  statistic) and Goelzer's density term (Eq. 10) at an ocean density of 1028 kg/m^3.

Arms S and G are diagnostics: the replay runs today's rule, so ice, isostasy and coastlines
follow `h_T`. The method was fixed in #912's first comment before any run.

    cd backend
    PYTHONPATH=. .venv/bin/python ../bin/debug/measure_flotation_debit.py \
        <save> --steps 10 --myr 0.25 [--ice-age-period-years 8e6] --out <captures>/<name>.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

BACKEND = Path(os.environ.get("MANTLE_BLOOM_BACKEND", str(Path(__file__).resolve().parents[2] / "backend")))
sys.path.insert(0, str(BACKEND))

from app import persistence  # noqa: E402
from app import plates as plates_mod  # noqa: E402
from app import world as world_mod  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app.hydroclimate import eustasy, ice_sheet  # noqa: E402

RHO_ICE = 917.0
RHO_WATER = 1000.0
RHO_OCEAN = 1028.0  # Goelzer et al. (2020); the model itself has one water density
_TOLERANCE_M3 = 1.0e3  # far below a millimetre of sea level over any ocean
_MAX_ITERS = 200


def _bisect(volume_at, target_m3: float, lo: float, hi: float) -> float:
    """The `h` in `[lo, hi]` with `volume_at(h) == target_m3`, for a non-decreasing
    `volume_at`. Tighter than `eustasy._bisect_sea_level` (1 m over a mean cell), so
    differences between arms aren't solver tolerance."""
    if target_m3 <= 0.0:
        return lo
    for _ in range(_MAX_ITERS):
        mid = 0.5 * (lo + hi)
        v = volume_at(mid)
        if abs(v - target_m3) < _TOLERANCE_M3 or hi - lo < 1e-9:
            return mid
        if v < target_m3:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


class State:
    """One step's final fields, in eustasy's node order, and the solves on them."""

    def __init__(self, world) -> None:
        self.z, self.areas = eustasy.hypsometry(world)
        spacing = line_spacing_rad(world.node_density)
        acct = plates_mod.collect_all_accounting_areas_m2(world.plates, spacing)
        assert np.array_equal(acct, self.areas)
        self.ice_m = np.clip(plates_mod.collect_all_glacier_depth(world.plates), 0.0, None)
        self.w = eustasy._ICE_WATER_EQUIVALENT * self.ice_m  # water-equivalent column, m
        self.budget = float(world.ocean_water_volume_m3)
        self.ice_debit_today = float(np.dot(self.areas, self.w))
        # Everything else `trapped_water_volume_m3` counts: lakes, lake ice, river water.
        self.other_trapped = eustasy.trapped_water_volume_m3(world) - self.ice_debit_today
        anchor = eustasy._seed_anchor(world, self.z, None, None)
        self.seeded = anchor is not None
        # The level at which each node joins the connected ocean: the mask at `h` is
        # `levels <= h`, as in `_solve_sea_level_connected`. With no seed, eustasy falls
        # back to other masks; the gate saves all have one, and the output records it.
        self.levels = eustasy._seeded_flood_levels(self.z, *anchor) if self.seeded else self.z
        self.lo = float(self.z.min())
        self.hi = float(self.z.max()) + self.budget / float(self.areas.sum()) + 1.0

    def mask(self, h: float) -> np.ndarray:
        return self.levels <= h

    def depth(self, h: float) -> np.ndarray:
        """`D`: water depth over cells of the connected ocean at `h`, 0 elsewhere."""
        return np.where(self.mask(h), np.clip(h - self.z, 0.0, None), 0.0)

    def capacity(self, h: float) -> float:
        """Eustasy's `C(h)`: geometric capacity of the connected ocean, ice or not."""
        return float(np.dot(self.areas, self.depth(h)))

    def liquid(self, h: float) -> float:
        """The liquid ocean at `h`: capacity the ice doesn't occupy (Goelzer Eq. 5)."""
        return float(np.dot(self.areas, np.clip(self.depth(h) - self.w, 0.0, None)))

    def flotation_debit(self, h_ref: float) -> float:
        """D6's debit with `D` taken at `h_ref`: `sum A * max(0, 0.917 H - D)`."""
        return float(np.dot(self.areas, np.clip(self.w - self.depth(h_ref), 0.0, None)))

    def solve_capacity(self, ice_debit: float) -> float:
        return _bisect(self.capacity, self.budget - self.other_trapped - ice_debit, self.lo, self.hi)

    def solve_generalised(self) -> float:
        return _bisect(self.liquid, self.budget - self.other_trapped - self.ice_debit_today, self.lo, self.hi)

    def classes(self, h: float) -> dict:
        ice = self.ice_m > 0.0
        d = self.depth(h)
        marine = ice & self.mask(h) & (d > 0.0)
        grounded = marine & (self.w > d)
        floating = marine & ~grounded
        land = ice & ~marine
        out = {}
        for name, sel in (("land", land), ("grounded_marine", grounded), ("floating", floating)):
            vol = float(np.dot(self.areas[sel], self.ice_m[sel]))
            out[name] = {
                "cells": int(sel.sum()),
                "area_km2": float(self.areas[sel].sum()) / 1e6,
                "ice_km3": vol / 1e9,
                "water_equivalent_km3": eustasy._ICE_WATER_EQUIVALENT * vol / 1e9,
            }
        marine_area = float(self.areas[marine].sum())
        out["marine_mean_depth_m"] = float(np.dot(self.areas[marine], d[marine]) / marine_area) if marine_area else 0.0
        # Water equivalent of grounded marine ice held below flotation: what today over-debits.
        out["grounded_below_flotation_km3"] = float(np.dot(self.areas[grounded], d[grounded])) / 1e9
        return out


def measure(world, h_s_prev: float) -> dict:
    s = State(world)
    h_model = float(world.sea_level_m)
    h_t = s.solve_capacity(s.ice_debit_today)
    debit_s = s.flotation_debit(h_s_prev)
    h_s = s.solve_capacity(debit_s)
    h_g = s.solve_generalised()
    ocean_area = float(s.areas[s.mask(h_t) & (s.z < h_t)].sum())
    v_den = float(np.dot(s.areas, s.ice_m)) * (RHO_ICE / RHO_WATER - RHO_ICE / RHO_OCEAN)
    hydro = getattr(world, "hydrology_cache", None)
    land_ledger = getattr(hydro, "land_water_ledger", None) if hydro is not None else None
    glacier_ledger = getattr(hydro, "glacier_water_ledger", None) if hydro is not None else None
    return {
        "age_myr": float(world.elapsed_years) / 1e6,
        "seeded_mask": s.seeded,
        "classes": s.classes(h_t),
        "total_ice_km3": float(np.dot(s.areas, s.ice_m)) / 1e9,
        "land_water_residual_m3": float(land_ledger.residual_m3) if land_ledger is not None else None,
        "glacier_water_residual_m3": float(glacier_ledger.residual_m3) if glacier_ledger is not None else None,
        "ocean_area_km2": ocean_area / 1e6,
        "debit_today_km3": s.ice_debit_today / 1e9,
        "debit_d6_km3": debit_s / 1e9,
        "debit_d6_at_solved_level_km3": s.flotation_debit(h_g) / 1e9,
        "debit_difference_km3": (s.ice_debit_today - debit_s) / 1e9,
        "debit_difference_over_ocean_m": (s.ice_debit_today - debit_s) / ocean_area if ocean_area else 0.0,
        "sea_level_model_m": h_model,
        "h_today_m": h_t,
        "h_d6_simple_m": h_s,
        "h_generalised_m": h_g,
        "today_minus_d6_m": h_t - h_s,
        "simple_minus_generalised_m": h_s - h_g,
        "resolve_minus_model_m": h_t - h_model,
        "density_term_km3": v_den / 1e9,
        "density_term_over_ocean_m": v_den / ocean_area if ocean_area else 0.0,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("save", type=Path)
    ap.add_argument("--steps", type=int, default=10)
    ap.add_argument("--myr", type=float, default=0.25)
    ap.add_argument("--ice-age-period-years", type=float, help="replay with ice ages at this period")
    ap.add_argument("--out", type=Path, required=True)
    ns = ap.parse_args()

    raw = ns.save.read_bytes()
    world = persistence.load_world_bytes(raw)
    if ns.ice_age_period_years is not None:
        world.ice_age_period_years = ns.ice_age_period_years
    preset = next((k for k, v in ice_sheet.MARCH_COST_PRESETS.items() if v == ice_sheet.MARCH_COST), "custom")
    commit = subprocess.run(["git", "rev-parse", "HEAD"], cwd=BACKEND, check=True, capture_output=True, text=True).stdout.strip()

    h_s = float(world.sea_level_m)  # arm S starts from the save's level
    steps = []
    # The opening state: classes and debits before any step, with no lag to measure.
    opening = measure(world, h_s) if getattr(world, "hydrology_cache", None) is not None else None
    for i in range(ns.steps):
        start = time.perf_counter()
        world_mod.step_world(world, ns.myr * 1e6)
        wall = time.perf_counter() - start
        rec = measure(world, h_s)
        h_s = rec["h_d6_simple_m"]
        rec["step"] = i
        rec["wall_s"] = round(wall, 1)
        steps.append(rec)
        c = rec["classes"]
        print(
            f"step {i}: {wall:.0f}s age={rec['age_myr']:.2f} ice land/grounded/floating km3="
            f"{c['land']['ice_km3']:.0f}/{c['grounded_marine']['ice_km3']:.1f}/{c['floating']['ice_km3']:.1f} "
            f"h_T={rec['h_today_m']:.3f} T-S={rec['today_minus_d6_m']:.4f} S-G={rec['simple_minus_generalised_m']:.4f}",
            flush=True,
        )

    gaps = [abs(r["simple_minus_generalised_m"]) for r in steps]
    ns.out.write_text(json.dumps({
        "issue": 912,
        "commit": commit,
        "save": ns.save.name,
        "sha256": hashlib.sha256(raw).hexdigest(),
        "march_preset": preset,
        "ice_age_period_years": float(world.ice_age_period_years),
        "steps": ns.steps,
        "myr_per_step": ns.myr,
        "peak_rss_gb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9, 2),
        "max_abs_simple_minus_generalised_m": max(gaps) if gaps else None,
        "max_abs_today_minus_d6_m": max((abs(r["today_minus_d6_m"]) for r in steps), default=None),
        "opening": opening,
        "per_step": steps,
    }, indent=1) + "\n")
    print(f"wrote {ns.out} max|S-G|={max(gaps):.4f} m")


if __name__ == "__main__":
    main()
