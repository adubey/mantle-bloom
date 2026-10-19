#!/usr/bin/env python3
"""Replay one RFC 0012 seed and trace subduction, merges, quiet gates, and unloading.

The trace is opt-in and installed with runtime wrappers; it does not change model state or
save files. Run in series for the four RFC seeds, once against each worktree:

    PYTHONPATH=<worktree>/backend <worktree>/backend/.venv/bin/python \
      bin/debug/trace_rfc0012_ledger_gates.py --seed 326494039 --steps 500 \
      --step-myr 1 --out /private/tmp/rfc0012-ledger-gates/pilot-seed326494039.jsonl

Each JSONL row is one completed step. It includes ledger deltas and closure residuals,
per-front subduction candidates/fates, merge inventories/capacity, per-plate quiet-clock
eligibility/reset reasons, local boundary-rate/strain proxies, and erosion/unloading totals.
"""

from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import subprocess
import time

import numpy as np

from app import cratons, crust_transfer, lithosphere, merge_split, quad_tectonics
from app.elevation_lines import effective_is_continental_from_codes, line_spacing_rad
from app import world as world_module


SECONDS_PER_YEAR = 31_557_600.0


def _belt_mask(plate, elapsed_years: float) -> np.ndarray:
    belt_reader = getattr(cratons, "belt_mask", None)
    return belt_reader(plate, elapsed_years) if belt_reader else np.zeros(plate.node_count(), dtype=bool)


def _json(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): _json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json(v) for v in value]
    return value


def _ledgers(world) -> dict:
    return {
        name: dict(value)
        for name, value in vars(world).items()
        if isinstance(value, dict) and ("ledger" in name or "budget" in name)
    }


def _delta(after: dict, before: dict) -> dict:
    result = {}
    for ledger, values in after.items():
        old = before.get(ledger, {})
        keys = set(values) | set(old)
        changes = {key: float(values.get(key, 0.0)) - float(old.get(key, 0.0)) for key in keys}
        result[ledger] = {key: value for key, value in changes.items() if value != 0.0}
    return result


def _weighted_stats(values: np.ndarray, weights: np.ndarray, mask: np.ndarray) -> dict:
    selected = np.asarray(values)[mask]
    areas = np.asarray(weights)[mask]
    if not selected.size:
        return {"count": 0, "area_m2": 0.0}
    order = np.argsort(selected)
    sorted_values, sorted_areas = selected[order], areas[order]
    midpoint = np.cumsum(sorted_areas) - sorted_areas / 2.0
    median = float(np.interp(sorted_areas.sum() / 2.0, midpoint, sorted_values))
    return {
        "count": int(selected.size),
        "area_m2": float(areas.sum()),
        "min": float(selected.min()),
        "median_area_weighted": median,
        "mean_area_weighted": float(np.average(selected, weights=areas)),
        "max": float(selected.max()),
    }


def _plate_inventory(world, plate) -> dict:
    area = plate.accounting_areas_m2(line_spacing_rad(world.node_density))
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    craton = plate.collect("craton_crust_m")
    return {
        "plate_id": int(plate.plate_id),
        "nodes": int(plate.node_count()),
        "area_m2": float(area.sum()),
        "hc_m3": float(hc @ area),
        "hm_m3": float(hm @ area),
        "craton_m3": float(craton @ area),
        "merge_hc_headroom_m3": float(np.maximum(lithosphere.MAX_CRUSTAL_THICKNESS_M - hc, 0.0) @ area),
    }


def _gate_snapshot(world, plate, years_myr: float) -> dict:
    area = cratons._areas(world, plate)
    hc = plate.collect("crustal_thickness_m")
    hm = plate.collect("mantle_lithosphere_thickness_m")
    elevation = plate.collect("elevation")
    craton = plate.collect("craton_crust_m")
    stable = plate.collect("stable_continental_myr")
    host = cratons._genuine_continental(plate)
    spacing = line_spacing_rad(world.node_density)
    hops = cratons._hops(cratons.CRATON_FORMATION_MARGIN_KM, spacing)
    interior = cratons._interior_hops(plate, host, hops)
    if interior is None:
        interior = np.zeros(plate.node_count(), dtype=int)
    belt = _belt_mask(plate, world.elapsed_years)
    reason_masks = {
        "not_genuine_continental": ~host,
        "host_too_thick": host & (hc > cratons.CRATON_HOST_MAX_HC_M),
        "outside_interior_margin": host & (interior < hops),
        "belt_exclusion": host & (interior >= hops) & (hc <= cratons.CRATON_HOST_MAX_HC_M) & belt,
        "already_cratonized": host & (interior >= hops) & (hc <= cratons.CRATON_HOST_MAX_HC_M) & ~belt & (craton > 0.0),
    }
    quiet = host & (hc <= cratons.CRATON_HOST_MAX_HC_M) & (interior >= hops) & ~belt
    erosion_trace = _erosion_plate_trace(world, plate, area)
    return {
        "plate_id": int(plate.plate_id),
        "nodes": int(plate.node_count()),
        "quiet_area_m2": float(area[quiet].sum()),
        "eligible_uncratonized_area_m2": float(area[quiet & (craton <= 0.0)].sum()),
        "mature_area_m2": float(area[quiet & (craton <= 0.0) & (stable + years_myr >= cratons.CRATON_FORMATION_MYR)].sum()),
        "reset_reason": {
            reason: {
                "count": int(mask.sum()),
                "area_m2": float(area[mask].sum()),
                "prior_clock_area_myr_m2": float((stable[mask] * area[mask]).sum()),
            }
            for reason, mask in reason_masks.items()
        },
        "quiet_hc_km": _weighted_stats(hc / 1000.0, area, quiet),
        "quiet_hm_km": _weighted_stats(hm / 1000.0, area, quiet),
        "quiet_elevation_km": _weighted_stats(elevation / 1000.0, area, quiet),
        "quiet_clock_myr": _weighted_stats(stable, area, quiet),
        "erosion_unloading": erosion_trace,
    }


def _erosion_plate_trace(world, plate, area) -> dict | None:
    result = getattr(world, "erosion_cache", None)
    if result is None:
        return None
    if len(result.points) != sum(p.node_count() for p in world.plates):
        return None
    offset = 0
    found = False
    for current in world.plates:
        n = current.node_count()
        if current.plate_id == plate.plate_id:
            found = True
            break
        offset += n
    if not found:
        return None
    sl = slice(offset, offset + plate.node_count())
    net = result.net_elevation_change_m[sl]
    load_change = result.ice_load_change_pa[sl] if result.ice_load_change_pa is not None else np.zeros_like(net)
    deflection = result.ice_deflection_change_m[sl] if result.ice_deflection_change_m is not None else np.zeros_like(net)
    return {
        "net_elevation_change_volume_m3": float(net @ area),
        "net_lowered_area_m2": float(area[net < 0.0].sum()),
        "net_raised_area_m2": float(area[net > 0.0].sum()),
        "ice_unloading_area_m2": float(area[load_change < 0.0].sum()),
        "ice_load_change_pa_area_mean": float(np.average(load_change, weights=area)),
        "ice_deflection_change_m_area_mean": float(np.average(deflection, weights=area)),
        "ice_deflection_change_m_min": float(deflection.min(initial=0.0)),
        "ice_deflection_change_m_max": float(deflection.max(initial=0.0)),
    }


def _install_traces(world, output: Path, model_root: Path, seed: int, step_myr: float):
    pending: list[dict] = []
    runtime = {"step": 0}

    def emit(kind: str, payload: dict) -> None:
        pending.append({"kind": kind, **_json(payload)})

    original_retreat = quad_tectonics._retreat
    original_context = quad_tectonics.boundary_context
    original_deform_columns = quad_tectonics.deform_columns

    def traced_deform_columns(
        current_world, plate, ctx, fields, local_xyz, node_area_m2, rng_index, years,
        ceiling_overflow=None, strained=None, accommodate=None,
    ):
        result = original_deform_columns(
            current_world, plate, ctx, fields, local_xyz, node_area_m2, rng_index, years,
            ceiling_overflow=ceiling_overflow, strained=strained, accommodate=accommodate,
        )
        if strained is not None:
            hc0 = np.asarray(fields["crustal_thickness_m"])
            hm0 = np.asarray(fields["mantle_lithosphere_thickness_m"])
            hc1 = np.asarray(strained["crustal_thickness_m"])
            hm1 = np.asarray(strained["mantle_lithosphere_thickness_m"])
            area = np.broadcast_to(np.asarray(node_area_m2, dtype=float), hc0.shape)
            host = effective_is_continental_from_codes(
                fields["crust_type_code"], plate.crust_type == "continental"
            ) & (hc0 >= cratons.CRATON_HOST_MIN_HC_M) & (hc0 <= cratons.CRATON_HOST_MAX_HC_M)
            interior = cratons._interior_hops(
                plate, host,
                cratons._hops(cratons.CRATON_FORMATION_MARGIN_KM, line_spacing_rad(current_world.node_density)),
            )
            if interior is not None and len(interior) == len(host):
                quiet_host = host & (interior >= cratons._hops(
                    cratons.CRATON_FORMATION_MARGIN_KM, line_spacing_rad(current_world.node_density)
                )) & ~_belt_mask(plate, current_world.elapsed_years)
            else:
                quiet_host = np.zeros(len(host), dtype=bool)
            emit("local_lithosphere_deformation", {
                "step": runtime["step"],
                "elapsed_myr": current_world.elapsed_years / 1e6,
                "plate_id": plate.plate_id,
                "cells_with_hc_change": int(np.count_nonzero(hc1 != hc0)),
                "cells_with_hm_change": int(np.count_nonzero(hm1 != hm0)),
                "melted_cells": int(np.asarray(strained["melted"], dtype=bool).sum()),
                "hc_delta_m3": float((hc1 - hc0) @ area),
                "hm_delta_m3": float((hm1 - hm0) @ area),
                "quiet_host_hc_delta_m3": float(((hc1 - hc0) * quiet_host) @ area),
                "quiet_host_hm_delta_m3": float(((hm1 - hm0) * quiet_host) @ area),
                "quiet_host_cells": int(quiet_host.sum()),
                "quiet_host_area_m2": float(area[quiet_host].sum()),
                "quiet_host_hc_change_m": _weighted_stats(hc1 - hc0, area, quiet_host),
                "quiet_host_hm_change_m": _weighted_stats(hm1 - hm0, area, quiet_host),
                "process_note": "strain-only Hc/Hm deltas before arc/rift magma and thermal evolution; captures the physical quiet-gate history directly",
            })
        return result

    quad_tectonics.deform_columns = traced_deform_columns

    def traced_context(current_world, plate, other_plates, years, continental_retreat_runs, node_weight=1.0):
        ctx = original_context(current_world, plate, other_plates, years, continental_retreat_runs, node_weight)
        area = plate.node_areas_m2()
        if len(area) != len(ctx.closing_rate):
            area = np.ones(len(ctx.closing_rate), dtype=float)
        hc = np.asarray(ctx.inputs.own_hc)
        host = effective_is_continental_from_codes(
            ctx.inputs.own_crust_type_codes, plate.crust_type == "continental"
        ) & (hc >= cratons.CRATON_HOST_MIN_HC_M)
        interior = cratons._interior_hops(plate, host, cratons._hops(
            cratons.CRATON_FORMATION_MARGIN_KM, line_spacing_rad(current_world.node_density)
        ))
        if interior is None or len(interior) != len(host):
            interior_mask = np.zeros(len(host), dtype=bool)
        else:
            interior_mask = interior >= cratons._hops(
                cratons.CRATON_FORMATION_MARGIN_KM, line_spacing_rad(current_world.node_density)
            )
        quiet_host = (
            host & (hc <= cratons.CRATON_HOST_MAX_HC_M) & interior_mask
            & ~_belt_mask(plate, current_world.elapsed_years)
        )
        rate = np.asarray(ctx.closing_rate) * 100.0 * SECONDS_PER_YEAR
        local_strain = np.divide(
            np.abs(ctx.closing_rate) * years * SECONDS_PER_YEAR,
            max(float(ctx.spacing_rad) * 6_371_000.0, 1.0),
        )
        emit("boundary_deformation", {
            "step": runtime["step"],
            "elapsed_myr": current_world.elapsed_years / 1e6,
            "plate_id": plate.plate_id,
            "candidate_contact_cells": int(np.asarray(ctx.contested).sum()),
            "convergent_cells": int(np.asarray(ctx.convergent).sum()),
            "divergent_cells": int(np.asarray(ctx.divergent).sum()),
            "transform_cells": int(np.asarray(ctx.transform).sum()),
            "quiet_host_contact_cells": int((quiet_host & ctx.contested).sum()),
            "contact_closing_cm_per_yr": _weighted_stats(rate, area, np.asarray(ctx.contested, dtype=bool)),
            "quiet_host_contact_closing_cm_per_yr": _weighted_stats(rate, area, quiet_host & ctx.contested),
            "contact_strain_per_step_proxy": _weighted_stats(local_strain, area, np.asarray(ctx.contested, dtype=bool)),
            "quiet_host_contact_strain_per_step_proxy": _weighted_stats(local_strain, area, quiet_host & ctx.contested),
            "proxy_note": "relative normal displacement divided by nominal cell spacing; this is a kinematic proxy, not constitutive strain",
        })
        return ctx

    quad_tectonics.boundary_context = traced_context

    def traced_retreat(plate, current_world, ctx, max_distance, max_cells, years=0.0):
        # The context is captured before per-column deformation; use its aligned arrays here
        # rather than collecting plate fields after deform_columns has written its result.
        hc = np.asarray(ctx.inputs.own_hc)
        hm = np.asarray(ctx.inputs.own_hm)
        craton = np.asarray(ctx.inputs.own_craton_m)
        continental = effective_is_continental_from_codes(
            ctx.inputs.own_crust_type_codes, plate.crust_type == "continental"
        )
        fate = quad_tectonics._retreat_disposition(ctx, continental)
        candidates = np.asarray(ctx.shrinkable, dtype=bool)
        before = _ledgers(current_world)
        survivors = original_retreat(plate, current_world, ctx, max_distance, max_cells, years)
        removed = ~np.asarray(survivors, dtype=bool)
        areas = np.ones(len(removed), dtype=float)
        selected_neighbours = np.asarray(ctx.inputs.neighbor_plate_id)[candidates]
        neighbor_counts = Counter(int(x) for x in selected_neighbours if x >= 0)
        closing_cm_yr = np.asarray(ctx.closing_rate) * 100.0 * SECONDS_PER_YEAR
        strain = np.divide(
            np.abs(ctx.closing_rate) * years * SECONDS_PER_YEAR,
            np.maximum(float(ctx.spacing_rad) * 6_371_000.0, 1.0),
        )
        emit("subduction_front", {
            "step": runtime["step"],
            "elapsed_myr": current_world.elapsed_years / 1e6,
            "plate_id": plate.plate_id,
            "candidate_cells": int(candidates.sum()),
            "retreatable_noncraton_cells": int((candidates & (craton <= 0.0)).sum()),
            "eligible_donor_cells": int((candidates & fate.donor).sum()),
            "eligible_budgeted_cells": int((candidates & fate.budgeted).sum()),
            "removed_cells": int(removed.sum()),
            "removed_donor_cells": int((removed & fate.donor).sum()),
            "removed_deep_subduction_cells": int((removed & ~fate.donor).sum()),
            "budget_blocked_cells": int((candidates & fate.budgeted & ~removed).sum()),
            "removed_hc_km_cell_sum": float(hc[removed].sum() / 1000.0),
            "removed_craton_km_cell_sum": float(craton[removed].sum() / 1000.0),
            "removed_hm_km_cell_sum": float(hm[removed].sum() / 1000.0),
            "volume_note": "thickness sums are per-cell diagnostic proxies; ledger deltas below are exact volumes",
            "neighbor_candidate_cell_counts": dict(neighbor_counts),
            "removed_contact_closing_cm_per_yr": _weighted_stats(closing_cm_yr, areas, removed),
            "removed_contact_strain_proxy": _weighted_stats(strain, areas, removed),
            "ledger_delta": _delta(_ledgers(current_world), before),
        })
        return survivors

    quad_tectonics._retreat = traced_retreat

    original_book = crust_transfer._book

    def traced_book(current_world, column, upper, result, thrust_up):
        before = _ledgers(current_world)
        original_book(current_world, column, upper, result, thrust_up)
        emit("suture_transfer", {
            "step": runtime["step"],
            "elapsed_myr": current_world.elapsed_years / 1e6,
            "lower_plate_id": column.plate_id,
            "upper_plate_id": upper.plate_id,
            "requested_hc_m3": column.hc_m3,
            "requested_craton_m3": column.craton_m3,
            "requested_material_m3": column.material_m3,
            "underthrust_m3": result.underthrust_m3,
            "scraped_m3": result.scraped_m3,
            "upper_placed_m3": result.upper_placed_m3,
            "handed_m3": result.handed_m3,
            "lost_m3": result.lost_m3,
            "unplaced_m3": result.unplaced_m3,
            "ledger_delta": _delta(_ledgers(current_world), before),
        })

    crust_transfer._book = traced_book

    original_merge = merge_split._merge_plates_with_craton_audit

    def traced_merge(current_world, id_keep, id_absorb):
        keep = next(p for p in current_world.plates if p.plate_id == id_keep)
        absorb = next(p for p in current_world.plates if p.plate_id == id_absorb)
        before_ledgers = _ledgers(current_world)
        pair_before = [_plate_inventory(current_world, keep), _plate_inventory(current_world, absorb)]
        pair_progress = {
            "collision_progress_years": current_world.collision_progress.get(tuple(sorted((id_keep, id_absorb)))),
            "overlap_progress_years": current_world.overlap_progress.get(tuple(sorted((id_keep, id_absorb)))),
        }
        before_craton = cratons.live_volume_m3(current_world)
        original_merge(current_world, id_keep, id_absorb)
        survivor = next(p for p in current_world.plates if p.plate_id == id_keep)
        after_craton = cratons.live_volume_m3(current_world)
        emit("plate_merge", {
            "step": runtime["step"],
            "elapsed_myr": current_world.elapsed_years / 1e6,
            "keep_plate_id": id_keep,
            "absorbed_plate_id": id_absorb,
            "sustained_pair_history": pair_progress,
            "pair_before": pair_before,
            "survivor_after": _plate_inventory(current_world, survivor),
            "craton_live_before_m3": before_craton,
            "craton_live_after_m3": after_craton,
            "craton_live_delta_m3": after_craton - before_craton,
            "merge_capped_delta_m3": _delta(_ledgers(current_world), before_ledgers).get("craton_ledger", {}).get("merge_capped_m3", 0.0),
            "ledger_delta": _delta(_ledgers(current_world), before_ledgers),
        })

    merge_split._merge_plates_with_craton_audit = traced_merge

    original_update = cratons.update

    def traced_update(current_world, years):
        cratons.ensure_initialized(current_world)
        gate_rows = [_gate_snapshot(current_world, plate, years / 1e6) for plate in current_world.plates if plate.node_count()]
        before_by_plate = {
            plate.plate_id: _plate_inventory(current_world, plate)["craton_m3"]
            for plate in current_world.plates if plate.node_count()
        }
        before = float(current_world.craton_ledger.get("formed_m3", 0.0))
        original_update(current_world, years)
        formed = float(current_world.craton_ledger.get("formed_m3", 0.0)) - before
        for row in gate_rows:
            plate = next((p for p in current_world.plates if p.plate_id == row["plate_id"]), None)
            formed_on_plate = 0.0
            if plate is not None:
                formed_on_plate = max(
                    _plate_inventory(current_world, plate)["craton_m3"] - before_by_plate[row["plate_id"]],
                    0.0,
                )
            row["step"] = runtime["step"]
            row["elapsed_myr"] = current_world.elapsed_years / 1e6
            row["formed_this_step_m3"] = formed_on_plate
            emit("quiet_gate", row)
        if formed > 0.0:
            emit("formation_total", {"step": runtime["step"], "formed_this_step_m3": formed})

    cratons.update = traced_update

    original_step = world_module.step_world

    def traced_step(current_world, years, **kwargs):
        runtime["step"] = current_world.steps_taken + 1
        before_ledgers = _ledgers(current_world)
        started = time.perf_counter()
        original_step(current_world, years, **kwargs)
        geology = _world_snapshot(current_world)
        erosion_world = _erosion_world_trace(current_world)
        row = {
            "kind": "step",
            "step": current_world.steps_taken,
            "elapsed_myr": current_world.elapsed_years / 1e6,
            "runtime_seconds": time.perf_counter() - started,
            "model": str(model_root),
            "commit": _commit(model_root),
            "seed": seed,
            "geology": geology,
            "ledger_delta": _delta(_ledgers(current_world), before_ledgers),
            "craton_balance_error_m3": cratons.balance_error_m3(current_world),
            "continental_material_balance_error_m3": _continental_balance(current_world),
            "erosion_unloading": erosion_world,
            "topology_lineage": _json(getattr(current_world, "topology_lineage", None)),
            "events": pending.copy(),
        }
        with output.open("a") as handle:
            handle.write(json.dumps(row, separators=(",", ":"), allow_nan=False) + "\n")
        pending.clear()
        return None

    world_module.step_world = traced_step
    return {
        "restore": lambda: _restore(
            original_step, original_update, original_merge, original_book, original_retreat,
            original_context, original_deform_columns,
        ),
    }


def _restore(step, update, merge, book, retreat, context, deform_columns):
    world_module.step_world = step
    cratons.update = update
    merge_split._merge_plates_with_craton_audit = merge
    crust_transfer._book = book
    quad_tectonics._retreat = retreat
    quad_tectonics.boundary_context = context
    quad_tectonics.deform_columns = deform_columns


def _world_snapshot(world) -> dict:
    area = np.concatenate([p.accounting_areas_m2(line_spacing_rad(world.node_density)) for p in world.plates if p.node_count()])
    elevation = np.concatenate([p.collect("elevation") for p in world.plates if p.node_count()])
    hc = np.concatenate([p.collect("crustal_thickness_m") for p in world.plates if p.node_count()])
    hm = np.concatenate([p.collect("mantle_lithosphere_thickness_m") for p in world.plates if p.node_count()])
    return {
        "plates": len(world.plates),
        "nodes": int(len(area)),
        "land_fraction": float(area[elevation > 0.0].sum() / area.sum()),
        "mean_hc_km": float(np.average(hc / 1000.0, weights=area)),
        "mean_hm_km": float(np.average(hm / 1000.0, weights=area)),
    }


def _erosion_world_trace(world) -> dict | None:
    result = getattr(world, "erosion_cache", None)
    if result is None:
        return None
    areas = np.concatenate([p.node_areas_m2() for p in world.plates if p.node_count()])
    net = result.net_elevation_change_m
    load = result.ice_load_change_pa if result.ice_load_change_pa is not None else np.zeros_like(net)
    deflection = result.ice_deflection_change_m if result.ice_deflection_change_m is not None else np.zeros_like(net)
    if len(areas) != len(net):
        return {"alignment_error": True, "areas": len(areas), "erosion_nodes": len(net)}
    return {
        "net_geomorphic_elevation_change_m3": float(net @ areas),
        "net_lowered_area_m2": float(areas[net < 0.0].sum()),
        "net_raised_area_m2": float(areas[net > 0.0].sum()),
        "ice_unloading_area_m2": float(areas[load < 0.0].sum()),
        "ice_load_change_pa_area_mean": float(np.average(load, weights=areas)),
        "ice_deflection_change_m_area_mean": float(np.average(deflection, weights=areas)),
        "ice_deflection_change_m_min": float(deflection.min(initial=0.0)),
        "ice_deflection_change_m_max": float(deflection.max(initial=0.0)),
    }


def _continental_balance(world) -> float:
    from app import continental_ledger
    return continental_ledger.balance_error_m3(world)


def _commit(root: Path) -> str:
    return subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()


def run(seed: int, steps: int, step_myr: float, output: Path, model_root: Path) -> None:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite trace: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    world = world_module.generate_world(seed=seed, node_density=1.0)
    _install_traces(world, output, model_root, seed, step_myr)
    dt_years = step_myr * 1_000_000.0
    started = time.perf_counter()
    try:
        for step in range(steps):
            world_module.step_world(world, dt_years)
            if (step + 1) % 25 == 0:
                print(f"seed={seed} step={step + 1}/{steps} elapsed_myr={(step + 1) * step_myr:g} runtime_min={(time.perf_counter() - started) / 60:.1f}", flush=True)
    finally:
        # Each run is a fresh interpreter, but restore wrappers to keep smoke imports reusable.
        pass
    summary = {
        "kind": "summary",
        "seed": seed,
        "steps": steps,
        "step_myr": step_myr,
        "elapsed_myr": steps * step_myr,
        "commit": _commit(model_root),
        "runtime_seconds": time.perf_counter() - started,
        "final_geology": _world_snapshot(world),
        "ledgers": _ledgers(world),
        "craton_balance_error_m3": cratons.balance_error_m3(world),
        "continental_material_balance_error_m3": _continental_balance(world),
    }
    with output.open("a") as handle:
        handle.write(json.dumps(summary, separators=(",", ":"), allow_nan=False) + "\n")
    print(json.dumps(summary, sort_keys=True), flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--model-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    if args.steps < 1 or args.step_myr <= 0:
        parser.error("steps and step-myr must be positive")
    run(args.seed, args.steps, args.step_myr, args.out, args.model_root.resolve())


if __name__ == "__main__":
    main()
