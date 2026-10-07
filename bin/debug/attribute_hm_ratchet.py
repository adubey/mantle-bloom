#!/usr/bin/env python3
"""Issue #310: attribute the mantle-lithosphere (Hm) ratchet toward the Hm cap to step phases.

Replays a saved world with `World.debug_diagnostics` on, so `phase_budget` (issue #216) books
every instrumented phase's area-weighted Hm volume change, and adds what that budget can't see:

- `convergent_clip`: Hm that `rheology.apply_convergent_deformation` grew and then clipped at
  `MAX_MANTLE_LITHOSPHERE_THICKNESS_M` -- the implicit "delamination" sink the issue questions;
- `unattributed`: the step's whole-world Hm volume change minus the sum of every recorded
  phase's change -- writers no phase brackets (gap fill, magma transport, splits, ...);
- cap-saturation metrics at each checkpoint (global, cratonic/non-cratonic, per plate);
- land, emergent continental area and hypsometry, and issue #314's collisional-shortening
  ledger (`phase_budget.SHORTENING_PHASE`: demanded = absorbed + returned to the boundary)
  where the engine books one.

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
import functools
import hashlib
import inspect
import json
import resource
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import elevation_lines, hm_ledger, lithosphere, persistence, plates, quad_tectonics, rheology, world as world_mod  # noqa: E402

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
        elevation = np.asarray(plate.collect("elevation"), dtype=float) - world.sea_level_m
        land = elevation > 0.0
        for name, arr in (("hc", hc), ("hm", hm), ("area", area), ("cont", cont), ("craton", craton), ("land", land), ("elevation", elevation)):
            parts[name].append(arr)
        cont_area = float(area[cont].sum())
        if cont_area > 0.0:
            per_plate[str(plate.plate_id)] = {
                "continental_area_km2": cont_area / 1e6,
                "hm_cap_fraction_of_continental": float(area[cont & (hm >= HM_CAP - CAP_TOLERANCE_M)].sum()) / cont_area,
            }
    hc, hm, area, cont, craton, land, elevation = (
        np.concatenate(parts[k]) for k in ("hc", "hm", "area", "cont", "craton", "land", "elevation")
    )
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
        "emergent_continental_area_km2": float(area[cont & land].sum()) / 1e6,
        "land_elevation_m": {f"p{q}": _weighted_quantile(elevation[land], area[land], q / 100.0) for q in (10, 50, 90, 99)},
        "continental_hc_km": {f"p{q}": _weighted_quantile(hc[cont], area[cont], q / 100.0) / 1e3 for q in (50, 90, 99)},
        "continental_hc_above_50km_area_km2": float(area[cont & (hc >= 50_000.0)].sum()) / 1e6,
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


def hm_state(world) -> dict:
    """Lightweight per-step live inventory and already-capped area, split by host state."""
    spacing = elevation_lines.line_spacing_rad(world.node_density)
    volume = defaultdict(float)
    capped_area = defaultdict(float)
    for plate in world.plates:
        hm = np.asarray(plate.collect("mantle_lithosphere_thickness_m"), dtype=float)
        if len(hm) == 0:
            continue
        area = np.asarray(plate.accounting_areas_m2(spacing), dtype=float)
        cont = elevation_lines.effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental"
        )
        craton = np.asarray(plate.collect("craton_crust_m"), dtype=float) > 0.0
        capped = hm >= HM_CAP - CAP_TOLERANCE_M
        for name, mask in (
            ("all", np.ones(len(hm), dtype=bool)),
            ("continental_node", cont),
            ("oceanic_node", ~cont),
            ("continental_craton_node", cont & craton),
            ("continental_non_craton_node", cont & ~craton),
        ):
            volume[name] += float(np.dot(hm[mask], area[mask])) / 1e9
            capped_area[name] += float(area[mask & capped].sum()) / 1e6
    return {"live_hm_km3": dict(volume), "already_capped_area_km2": dict(capped_area)}


def cap_transition_totals_km2(budget: dict) -> dict:
    out = {}
    for phase, totals in budget.items():
        cap = totals.get("hm_cap_transitions")
        if not cap:
            continue
        out[phase] = {
            scope: {key.replace("_m2", "_km2"): value / 1e6 for key, value in row.items()}
            for scope, row in cap["scopes"].items()
        }
    return out


def _nested_delta(after: dict, before: dict) -> dict:
    out = {}
    for outer, after_row in after.items():
        before_row = before.get(outer, {})
        rows = {}
        for inner, values in after_row.items():
            if not isinstance(values, dict):
                continue
            delta = {key: float(value) - float(before_row.get(inner, {}).get(key, 0.0)) for key, value in values.items()}
            if any(value != 0.0 for value in delta.values()):
                rows[inner] = delta
        if rows:
            out[outer] = rows
    return out


def account_totals_km3(world) -> dict:
    return {
        account: {
            scope: {key.replace("_m3", "_km3"): value / 1e9 for key, value in row.items()}
            for scope, row in entry["scopes"].items()
        }
        for account, entry in world.hm_source_sink_ledger.items()
    }


def suture_totals_km3(world) -> dict:
    budget = world.hm_suture_budget
    return {
        "fronts": budget.get("fronts", 0),
        "donor_hm_km3": budget.get("donor_hm_m3", 0.0) / 1e9,
        "placed_hm_km3": budget.get("placed_hm_m3", 0.0) / 1e9,
        "unplaced_hm_km3": budget.get("unplaced_hm_m3", 0.0) / 1e9,
        "by_pair": {
            key: {
                "fronts": row["fronts"],
                "donor_hm_km3": row["donor_hm_m3"] / 1e9,
                "placed_hm_km3": row["placed_hm_m3"] / 1e9,
                "unplaced_hm_km3": row["unplaced_hm_m3"] / 1e9,
            }
            for key, row in budget.get("by_pair", {}).items()
        },
    }


class RuntimeProfile:
    """Low-overhead debug-only timers for issue #317's replay cost comparison."""

    def __init__(self):
        self.stats = defaultdict(lambda: {"calls": 0, "seconds": 0.0})
        self.originals = []

    def _wrap(self, obj, name: str, label: str, *, only_when=None):
        original = getattr(obj, name)
        self.originals.append((obj, name, original))

        @functools.wraps(original)
        def wrapped(*args, **kwargs):
            enabled = only_when is None or only_when(*args, **kwargs)
            started = time.perf_counter() if enabled else 0.0
            try:
                return original(*args, **kwargs)
            finally:
                if enabled:
                    row = self.stats[label]
                    row["calls"] += 1
                    row["seconds"] += time.perf_counter() - started

        setattr(obj, name, wrapped)

    def __enter__(self):
        self._wrap(quad_tectonics, "boundary_context", "boundary_classification")
        self._wrap(quad_tectonics, "_spread_accretion_volume", "hm_spreading")
        self._wrap(quad_tectonics, "_place_suture_crust", "graph_traversal")
        self._wrap(quad_tectonics, "_adjacency_matrix", "spatial_graph_builds")
        self._wrap(
            plates.Plate,
            "get_node_kdtree",
            "spatial_index_builds",
            only_when=lambda plate, *args, **kwargs: getattr(plate, "_node_kdtree_cache", None) is None,
        )
        return self

    def __exit__(self, *exc):
        for obj, name, original in reversed(self.originals):
            setattr(obj, name, original)


def phase_deltas_km3(budget: dict) -> dict:
    """{phase: {scope: Hm volume after - before, km^3}} from a `World.phase_budget`."""
    return {
        phase: {scope: (t["hm_volume_after_m3"] - t["hm_volume_before_m3"]) / 1e9 for scope, t in totals["scopes"].items()}
        for phase, totals in budget.items()
        if "scopes" in totals
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

    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=Path(__file__).resolve().parents[2], check=True, capture_output=True, text=True
    ).stdout.strip()
    save_sha256 = hashlib.sha256(args.save.read_bytes()).hexdigest()
    result = {
        "issue": 317,
        "commit": commit,
        "save": str(args.save),
        "save_sha256": save_sha256,
        "step_years": args.step_years,
        "compare": {},
        "steps": [],
        "checkpoints": [],
    }
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
    result["checkpoints"].append({"metrics": metrics(world), "inventory_km3": inventory0, "hm_state": hm_state(world)})

    args.out.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with ClipTap() as tap, RuntimeProfile() as profile:
        for i in range(1, steps + 1):
            step_started = time.perf_counter()
            before_inv = hm_inventory_km3(world)
            before_phase = copy.deepcopy(phase_deltas_km3(world.phase_budget))
            before_caps = cap_transition_totals_km2(world.phase_budget)
            before_accounts = account_totals_km3(world)
            before_suture = suture_totals_km3(world)
            world_mod.step_world(world, args.step_years)
            after_inv = hm_inventory_km3(world)
            after_phase = phase_deltas_km3(world.phase_budget)
            after_caps = cap_transition_totals_km2(world.phase_budget)
            after_accounts = account_totals_km3(world)
            after_suture = suture_totals_km3(world)
            step_seconds = time.perf_counter() - step_started
            step_residual = {}
            ledger_residual = {}
            for scope_inv, scope_phase in (("all", "all"), ("continental_node", "continental_node"), ("oceanic_node", "oceanic_node")):
                recorded = sum(s.get(scope_phase, 0.0) for s in after_phase.values()) - sum(s.get(scope_phase, 0.0) for s in before_phase.values())
                residual = (after_inv[scope_inv] - before_inv[scope_inv]) - recorded
                unattributed[scope_inv] += residual
                step_residual[scope_inv] = residual
                scope_accounts = [row.get(scope_inv, {}) for row in after_accounts.values()]
                prior_scope_accounts = [row.get(scope_inv, {}) for row in before_accounts.values()]
                source_delta = sum(row.get("source_km3", 0.0) for row in scope_accounts) - sum(row.get("source_km3", 0.0) for row in prior_scope_accounts)
                sink_delta = sum(row.get("sink_km3", 0.0) for row in scope_accounts) - sum(row.get("sink_km3", 0.0) for row in prior_scope_accounts)
                reclass_delta = sum(row.get("reclassification_km3", 0.0) for row in scope_accounts) - sum(row.get("reclassification_km3", 0.0) for row in prior_scope_accounts)
                ledger_residual[scope_inv] = (after_inv[scope_inv] - before_inv[scope_inv]) - source_delta + sink_delta - reclass_delta
            result["steps"].append(
                {
                    "step": i,
                    "elapsed_myr": world.elapsed_years / 1e6,
                    "runtime_seconds": step_seconds,
                    **hm_state(world),
                    "newly_capped_by_phase_km2": _nested_delta(after_caps, before_caps),
                    "hm_source_sink_delta_km3": _nested_delta(after_accounts, before_accounts),
                    "hm_ledger_residual_km3": ledger_residual,
                    "signed_residual_km3": step_residual,
                    "suture_delta": {
                        "fronts": after_suture["fronts"] - before_suture["fronts"],
                        "donor_hm_km3": after_suture["donor_hm_km3"] - before_suture["donor_hm_km3"],
                        "placed_hm_km3": after_suture["placed_hm_km3"] - before_suture["placed_hm_km3"],
                        "unplaced_hm_km3": after_suture["unplaced_hm_km3"] - before_suture["unplaced_hm_km3"],
                    },
                }
            )
            if i % per_checkpoint == 0 or i == steps:
                entry = {
                    "step": i,
                    "metrics": metrics(world),
                    "inventory_km3": after_inv,
                    "phase_hm_delta_km3": after_phase,
                    "convergent_clip_km3": dict(tap.totals),
                    "unattributed_hm_delta_km3": dict(unattributed),
                    "shortening": dict(world.phase_budget.get("convergent_shortening", {}).get("shortening", {})),
                    "hm_source_sink_accounts_km3": after_accounts,
                    "hm_cap_transitions_km2": after_caps,
                    "hm_suture_budget_km3": after_suture,
                    "runtime_profile": dict(profile.stats),
                }
                result["checkpoints"].append(entry)
                m = entry["metrics"]
                print(
                    f"step {i}/{steps} t={m['elapsed_myr']:.1f} Myr  hm_cap={m['hm_cap_area_km2'] / 1e6:.2f} M km2 "
                    f"({m['continental']['hm_cap_fraction'] * 100:.2f}% cont)  cont={m['continental_area_km2'] / 1e6:.1f} M km2  "
                    f"land={m['land_area_km2'] / 1e6:.1f} M km2  clipped_cont={tap.totals['continental_hm_clipped_km3']:.3g} km3  elapsed={time.time() - started:.0f}s",
                    flush=True,
                )
                args.out.write_text(json.dumps(result, indent=1))
    result["runtime_seconds"] = time.time() - started
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    result["peak_memory_mb"] = peak / (1024.0 * 1024.0) if sys.platform == "darwin" else peak / 1024.0
    result["runtime_profile"] = dict(profile.stats)
    closure_accounts = hm_ledger.cumulative_scopes(world)
    final_inventory = hm_inventory_km3(world)
    result["hm_closure_km3"] = {
        scope: {
            "live_change_km3": final_inventory[scope] - inventory0[scope],
            "source_km3": closure_accounts[scope]["source_m3"] / 1e9,
            "sink_km3": closure_accounts[scope]["sink_m3"] / 1e9,
            "reclassification_km3": closure_accounts[scope]["reclassification_m3"] / 1e9,
            "ledger_residual_km3": (final_inventory[scope] - inventory0[scope])
            - closure_accounts[scope]["source_m3"] / 1e9
            + closure_accounts[scope]["sink_m3"] / 1e9
            - closure_accounts[scope]["reclassification_m3"] / 1e9,
            "phase_budget_residual_km3": unattributed[scope],
        }
        for scope in ("all", "continental_node", "oceanic_node")
    }
    args.out.write_text(json.dumps(result, indent=1))
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
