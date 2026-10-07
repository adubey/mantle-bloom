#!/usr/bin/env python3
"""Issue #310: attribute the mantle-lithosphere (Hm) ratchet toward the Hm cap to step phases.

Replays a saved world with `World.debug_diagnostics` on, so `phase_budget` (issue #216) books
every instrumented phase's area-weighted Hm volume change, and adds what that budget can't see:

- `convergent_clip`: Hm that `rheology.apply_convergent_deformation` grew and then clipped at
  `MAX_MANTLE_LITHOSPHERE_THICKNESS_M` -- the implicit "delamination" sink the issue questions;
- `unattributed`: the step's whole-world Hm volume change minus the sum of every recorded
  phase's change -- writers no phase brackets (gap fill, magma transport, splits, ...);
- cap-saturation metrics at each checkpoint (global, cratonic/non-cratonic, per plate).

Also prints the same metrics for `--compare` saves (e.g. the original 130.5 Myr checkpoint),
so the replay can be checked against the run that produced the issue's numbers.

    cd backend
    .venv/bin/python ../bin/debug/attribute_hm_ratchet.py \\
        ~/Downloads/mantle-bloom-seed997271774-87400000y.mbworld --myr 43.1 \\
        --compare ~/Downloads/mantle-bloom-seed997271774-130500000y.mbworld \\
        --out ../analysis/issue310/hm-budget-seed997271774.json
"""

from __future__ import annotations

import argparse
import copy
import inspect
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import elevation_lines, lithosphere, persistence, rheology, world as world_mod  # noqa: E402

HM_CAP = lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M
HC_CAP = lithosphere.MAX_CRUSTAL_THICKNESS_M
CAP_TOLERANCE_M = 1.0  # a node within 1 m of the cap counts as capped
CAP_FRACTIONS = (0.8, 0.9, 0.99)


def load(path: Path):
    return persistence.load_world_bytes(path.read_bytes())


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    if values.size == 0:
        return float("nan")
    order = np.argsort(values)
    cum = np.cumsum(weights[order])
    return float(values[order][np.searchsorted(cum, q * cum[-1])])


def metrics(world) -> dict:
    """Area-weighted Hm/Hc saturation metrics, globally and per plate (exact quad cell areas)."""
    spacing = elevation_lines.line_spacing_rad(world.node_density)
    parts = defaultdict(list)
    per_plate = {}
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        hc = np.asarray(plate.collect("crustal_thickness_m"), dtype=float)
        hm = np.asarray(plate.collect("mantle_lithosphere_thickness_m"), dtype=float)
        area = np.asarray(plate.accounting_areas_m2(spacing), dtype=float)
        cont = elevation_lines.effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        craton = np.asarray(plate.collect("craton_crust_m"), dtype=float) > 0.0
        land = np.asarray(plate.collect("elevation"), dtype=float) > world.sea_level_m
        for name, arr in (("hc", hc), ("hm", hm), ("area", area), ("cont", cont), ("craton", craton), ("land", land)):
            parts[name].append(arr)
        cont_area = float(area[cont].sum())
        if cont_area > 0.0:
            per_plate[str(plate.plate_id)] = {
                "continental_area_km2": cont_area / 1e6,
                "hm_cap_fraction_of_continental": float(area[cont & (hm >= HM_CAP - CAP_TOLERANCE_M)].sum()) / cont_area,
            }
    hc, hm, area, cont, craton, land = (np.concatenate(parts[k]) for k in ("hc", "hm", "area", "cont", "craton", "land"))
    capped = hm >= HM_CAP - CAP_TOLERANCE_M
    cont_area = float(area[cont].sum())

    def group(mask: np.ndarray) -> dict:
        a = area[mask]
        total = float(a.sum())
        if total == 0.0:
            return {"area_km2": 0.0}
        h = hm[mask]
        out = {
            "area_km2": total / 1e6,
            "hm_cap_area_km2": float(a[capped[mask]].sum()) / 1e6,
            "hm_cap_fraction": float(a[capped[mask]].sum()) / total,
            "hm_mean_km": float(np.dot(h, a)) / total / 1e3,
            "hm_p50_km": _weighted_quantile(h, a, 0.50) / 1e3,
            "hm_p95_km": _weighted_quantile(h, a, 0.95) / 1e3,
            "hm_p99_km": _weighted_quantile(h, a, 0.99) / 1e3,
            "hm_volume_km3": float(np.dot(h, a)) / 1e9,
        }
        for f in CAP_FRACTIONS:
            out[f"hm_above_{f:g}_cap_fraction"] = float(a[h >= f * HM_CAP].sum()) / total
        return out

    return {
        "elapsed_myr": world.elapsed_years / 1e6,
        "plate_count": len(world.plates),
        "sea_level_m": world.sea_level_m,
        "land_area_km2": float(area[land].sum()) / 1e6,
        "continental_area_km2": cont_area / 1e6,
        "hc_cap_area_km2": float(area[hc >= HC_CAP - CAP_TOLERANCE_M].sum()) / 1e6,
        "hm_cap_area_km2": float(area[capped].sum()) / 1e6,
        "hm_volume_all_km3": float(np.dot(hm, area)) / 1e9,
        "continental": group(cont),
        "continental_craton": group(cont & craton),
        "continental_non_craton": group(cont & ~craton),
        "oceanic": group(~cont),
        "per_plate": dict(sorted(per_plate.items(), key=lambda kv: -kv[1]["hm_cap_fraction_of_continental"] * kv[1]["continental_area_km2"])),
    }


class ClipTap:
    """Wraps `rheology.apply_convergent_deformation` to book the Hm it clips at the cap.
    Areas and node types come from `lithosphere_plate.deform_columns`' own locals (its only
    caller): `thicken` selects the slice passed in, `node_area_m2` is the slice's per-cell
    area, `codes0`/`plate` resolve continental vs oceanic, `fields` carries craton share."""

    def __init__(self):
        self.original = rheology.apply_convergent_deformation
        self.totals = defaultdict(float)

    def __enter__(self):
        tap = self

        def wrapped(hc_m, hm_m, closing_rate, years_myr, fault_factor, strength=1.0):
            new_hc, new_hm, overflow_hc = tap.original(hc_m, hm_m, closing_rate, years_myr, fault_factor, strength=strength)
            caller = inspect.currentframe().f_back.f_locals
            thicken = caller["thicken"]
            area = np.broadcast_to(caller["node_area_m2"], thicken.shape)[thicken]
            cont = elevation_lines.effective_is_continental_from_codes(caller["codes0"], caller["plate"].crust_type == "continental")[thicken]
            craton = np.asarray(caller["fields"].get("craton_crust_m", np.zeros(thicken.shape)))[thicken] > 0.0
            fractional = np.clip(rheology.plastic_strain_rate_per_myr(closing_rate), 0.0, None) * years_myr * fault_factor * strength
            uncapped = hm_m * (1.0 + fractional)
            clipped = (uncapped - new_hm) * area / 1e9
            added = (new_hm - hm_m) * area / 1e9
            for name, mask in (("continental", cont), ("continental_craton", cont & craton), ("oceanic", ~cont)):
                tap.totals[f"{name}_hm_added_km3"] += float(added[mask].sum())
                tap.totals[f"{name}_hm_clipped_km3"] += float(clipped[mask].sum())
                tap.totals[f"{name}_hc_clipped_km3"] += float((overflow_hc * area / 1e9)[mask].sum())
            return new_hc, new_hm, overflow_hc

        rheology.apply_convergent_deformation = wrapped
        return self

    def __exit__(self, *exc):
        rheology.apply_convergent_deformation = self.original


def hm_inventory_km3(world) -> dict:
    """Whole-world area-weighted Hm volume, all nodes and continental nodes."""
    spacing = elevation_lines.line_spacing_rad(world.node_density)
    out = {"all": 0.0, "continental_node": 0.0, "oceanic_node": 0.0}
    for plate in world.plates:
        if plate.node_count() == 0:
            continue
        vol = np.asarray(plate.collect("mantle_lithosphere_thickness_m"), dtype=float) * plate.accounting_areas_m2(spacing) / 1e9
        cont = elevation_lines.effective_is_continental_from_codes(plate.collect("crust_type_code"), plate.crust_type == "continental")
        out["all"] += float(vol.sum())
        out["continental_node"] += float(vol[cont].sum())
        out["oceanic_node"] += float(vol[~cont].sum())
    return out


def phase_deltas_km3(budget: dict) -> dict:
    """{phase: {scope: Hm volume after - before, km^3}} from a `World.phase_budget`."""
    return {
        phase: {scope: (t["hm_volume_after_m3"] - t["hm_volume_before_m3"]) / 1e9 for scope, t in totals["scopes"].items()}
        for phase, totals in budget.items()
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("save", type=Path)
    ap.add_argument("--myr", type=float, required=True, help="how far to replay")
    ap.add_argument("--step-years", type=float, default=100_000.0, help="UI default step")
    ap.add_argument("--checkpoint-myr", type=float, default=5.0)
    ap.add_argument("--compare", type=Path, nargs="*", default=[])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    result = {"save": str(args.save), "step_years": args.step_years, "compare": {}, "checkpoints": []}
    for path in args.compare:
        result["compare"][str(path)] = metrics(load(path))
        print(f"compare {path.name}: {json.dumps({k: v for k, v in result['compare'][str(path)].items() if k != 'per_plate'})}", flush=True)

    world = load(args.save)
    world.debug_diagnostics = True
    world.reset_phase_budget()
    steps = int(round(args.myr * 1e6 / args.step_years))
    per_checkpoint = max(1, int(round(args.checkpoint_myr * 1e6 / args.step_years)))
    inventory0 = hm_inventory_km3(world)
    unattributed = defaultdict(float)
    result["checkpoints"].append({"metrics": metrics(world), "inventory_km3": inventory0})

    args.out.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with ClipTap() as tap:
        for i in range(1, steps + 1):
            before_inv = hm_inventory_km3(world)
            before_phase = copy.deepcopy(phase_deltas_km3(world.phase_budget))
            world_mod.step_world(world, args.step_years)
            after_inv = hm_inventory_km3(world)
            after_phase = phase_deltas_km3(world.phase_budget)
            for scope_inv, scope_phase in (("all", "all"), ("continental_node", "continental_node"), ("oceanic_node", "oceanic_node")):
                recorded = sum(s.get(scope_phase, 0.0) for s in after_phase.values()) - sum(s.get(scope_phase, 0.0) for s in before_phase.values())
                unattributed[scope_inv] += (after_inv[scope_inv] - before_inv[scope_inv]) - recorded
            if i % per_checkpoint == 0 or i == steps:
                entry = {
                    "step": i,
                    "metrics": metrics(world),
                    "inventory_km3": after_inv,
                    "phase_hm_delta_km3": after_phase,
                    "convergent_clip_km3": dict(tap.totals),
                    "unattributed_hm_delta_km3": dict(unattributed),
                }
                result["checkpoints"].append(entry)
                m = entry["metrics"]
                print(
                    f"step {i}/{steps} t={m['elapsed_myr']:.1f} Myr  hm_cap={m['hm_cap_area_km2'] / 1e6:.2f} M km2 "
                    f"({m['continental']['hm_cap_fraction'] * 100:.2f}% cont)  cont={m['continental_area_km2'] / 1e6:.1f} M km2  "
                    f"clipped_cont={tap.totals['continental_hm_clipped_km3']:.3g} km3  elapsed={time.time() - started:.0f}s",
                    flush=True,
                )
                args.out.write_text(json.dumps(result, indent=1))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
