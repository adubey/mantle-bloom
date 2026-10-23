#!/usr/bin/env python3
"""Measure one G2 run under a selected glacial-erosion multiplier.

Run each seed/multiplier pair as a separate process to keep memory bounded. For example:

    PYTHONPATH=backend backend/.venv/bin/python bin/debug/g2_glacial_erosion_sensitivity.py \
      --seed 326494039 --multiplier 2 --steps 250 --step-myr 1 \
      --out /tmp/g2-glacier-326494039-2x.jsonl

The output records craton-formation gates, cumulative applied glacier erosion carried by
stable node IDs where the backend exposes per-node erosion, glacier-covered mountain area,
and local relief summaries. The checked-in backend currently records aggregate glacier
erosion only; per-node lineage requires the temporary diagnostic overlay used for RFC 0012's
matched replay. It is a sensitivity replay, not a calibrated glacier mass-balance experiment.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from app import cratons
from app.elevation_lines import line_spacing_rad
from app.world import generate_world, step_world
from g2_replay_common import add_replay_arguments, prepare_replay_output, write_replay_summary


def _weighted(values: np.ndarray, area: np.ndarray, mask: np.ndarray) -> float | None:
    if not np.any(mask):
        return None
    return float(np.dot(values[mask], area[mask]) / area[mask].sum())


def _node_id_keys(ids: np.ndarray) -> np.ndarray:
    keys = np.empty(len(ids), dtype=[("a", "<u8"), ("b", "<u8")])
    keys["a"] = ids[:, 0]
    keys["b"] = ids[:, 1]
    return keys


def _capture(world, elapsed_myr: float, step: int, cumulative_glacial_m3: float,
             cumulative_by_node_m: np.ndarray, step_glacier_m: np.ndarray) -> dict:
    spacing = line_spacing_rad(world.node_density)
    areas = []
    hcs = []
    hms = []
    elevations = []
    glaciers = []
    hosts = []
    interiors = []
    cratons_live = []
    stable_clocks = []
    hops = cratons._hops(cratons.CRATON_FORMATION_MARGIN_KM, spacing)
    for plate in world.plates:
        count = plate.node_count()
        if not count:
            continue
        area = cratons._areas(world, plate)
        host = cratons._genuine_continental(plate)
        interior = cratons._interior_hops(plate, host, hops)
        if interior is None:
            interior = np.zeros(count, dtype=int)
        areas.append(area)
        hcs.append(plate.collect("crustal_thickness_m"))
        hms.append(plate.collect("mantle_lithosphere_thickness_m"))
        elevations.append(plate.collect("elevation"))
        glaciers.append(plate.collect("glacier_depth"))
        hosts.append(host)
        interiors.append(interior)
        cratons_live.append(plate.collect("craton_crust_m") > 0.0)
        stable_clocks.append(plate.collect("stable_continental_myr"))

    if not areas:
        return {}
    area = np.concatenate(areas)
    hc = np.concatenate(hcs)
    hm = np.concatenate(hms)
    elevation = np.concatenate(elevations)
    glacier = np.concatenate(glaciers)
    host = np.concatenate(hosts)
    interior = np.concatenate(interiors)
    live_craton = np.concatenate(cratons_live)
    stable = np.concatenate(stable_clocks)
    eligible = host & (hc <= cratons.CRATON_HOST_MAX_HC_M) & (interior >= hops)
    forming_pool = eligible & ~live_craton

    # Hydrology's neighbor graph is gathered in the same plate/node order as the erosion
    # cache. The range of neighboring elevations is a simple, resolution-scale relief proxy.
    hydro = world.hydrology_cache
    if hydro is not None and len(hydro.neighbor_idx) == len(elevation):
        neighbor_elevation = elevation[hydro.neighbor_idx]
        local_relief_m = np.ptp(neighbor_elevation, axis=1)
    else:
        local_relief_m = np.zeros_like(elevation)
    glaciated = glacier >= 10.0
    high_relief = local_relief_m >= 1000.0
    mountain = elevation >= 2000.0
    return {
        "step": step,
        "elapsed_myr": elapsed_myr,
        "plate_count": len(world.plates),
        "cumulative_glacial_erosion_m3": cumulative_glacial_m3,
        "formed_total_m3": float(world.craton_ledger.get("formed_m3", 0.0)),
        "initial_craton_m3": float(world.craton_ledger.get("initial_m3", 0.0)),
        "live_craton_m3": sum(
            float(np.dot(plate.collect("craton_crust_m"), cratons._areas(world, plate)))
            for plate in world.plates if plate.node_count()
        ),
        "gates": {
            "eligible_uncratonized_km2": float(area[forming_pool].sum() / 1e6),
            "eligible_uncratonized_hc_mean_km": _weighted(hc / 1000.0, area, forming_pool),
            "eligible_uncratonized_hm_mean_km": _weighted(hm / 1000.0, area, forming_pool),
            "eligible_uncratonized_glacier_covered_km2": float(area[forming_pool & glaciated].sum() / 1e6),
            "eligible_uncratonized_relief_mean_m": _weighted(local_relief_m, area, forming_pool),
            "eligible_uncratonized_glacial_erosion_mean_m": _weighted(cumulative_by_node_m, area, forming_pool),
            "clock_ge_180_myr_km2": float(area[forming_pool & (stable >= 180.0)].sum() / 1e6),
            "clock_ge_190_myr_km2": float(area[forming_pool & (stable >= 190.0)].sum() / 1e6),
        },
        "mountains": {
            "area_ge_2km_km2": float(area[mountain].sum() / 1e6),
            "glaciated_area_ge_2km_km2": float(area[mountain & glaciated].sum() / 1e6),
            "glaciated_high_relief_area_km2": float(area[glaciated & high_relief].sum() / 1e6),
            "local_relief_mean_m": _weighted(local_relief_m, area, mountain),
            "glaciated_relief_mean_m": _weighted(local_relief_m, area, mountain & glaciated),
            "glacial_erosion_mean_m": _weighted(cumulative_by_node_m, area, mountain),
            "glacial_erosion_p90_m": float(np.quantile(cumulative_by_node_m[mountain], 0.9)) if np.any(mountain) else None,
        },
        "last_step_glacial_erosion_m3": float(
            world.erosion_cache.budget.get("glacial_erosion_m3", 0.0)
            if world.erosion_cache is not None else 0.0
        ),
        "last_step_glacial_erosion_by_node_mean_m": _weighted(step_glacier_m, area, np.ones(len(area), dtype=bool)),
    }


def run(seed: int, multiplier: float, steps: int, step_myr: float,
        snapshot_every: int, out: Path) -> None:
    world = generate_world(seed=seed, node_density=1.0)
    world.glacier_erosion_multiplier = multiplier
    original_update = cratons.update
    rows = []
    cumulative_glacial_m3 = 0.0
    prior_ids = None
    prior_cumulative = np.zeros(0)
    step_glacier = np.zeros(0)
    started = time.perf_counter()

    def wrapped_update(current_world, years):
        nonlocal cumulative_glacial_m3, prior_ids, prior_cumulative, step_glacier
        erosion_result = current_world.erosion_cache
        if erosion_result is not None and erosion_result.budget:
            cumulative_glacial_m3 += float(erosion_result.budget.get("glacial_erosion_m3", 0.0))
            applied = np.asarray(erosion_result.budget.get("glacial_erosion_m_by_node", []))
            ids = np.concatenate([
                plate.surface_nodes().node_ids for plate in current_world.plates if plate.node_count()
            ]) if current_world.plates else np.zeros((0, 2), dtype=np.uint64)
            step_glacier = applied if len(applied) == len(ids) else np.zeros(len(ids))
            if len(step_glacier) == len(ids):
                if prior_ids is None:
                    cumulative = step_glacier.copy()
                else:
                    old_keys = _node_id_keys(prior_ids)
                    order = np.argsort(old_keys)
                    sorted_keys = old_keys[order]
                    current_keys = _node_id_keys(ids)
                    positions = np.searchsorted(sorted_keys, current_keys)
                    valid = positions < len(sorted_keys)
                    if len(sorted_keys):
                        valid &= sorted_keys[np.minimum(positions, len(sorted_keys) - 1)] == current_keys
                    cumulative = step_glacier.copy()
                    cumulative[valid] += prior_cumulative[order[positions[valid]]]
                prior_ids = ids.copy()
                prior_cumulative = cumulative

        before = float(current_world.craton_ledger.get("formed_m3", 0.0))
        original_update(current_world, years)
        formed_this_step = float(current_world.craton_ledger.get("formed_m3", 0.0)) - before
        step = int(round(current_world.elapsed_years / years)) if years else 0
        if step % snapshot_every == 0 or formed_this_step > 0.0 or step == 1 or step == steps:
            row = _capture(
                current_world, current_world.elapsed_years / 1e6, step,
                cumulative_glacial_m3, prior_cumulative, step_glacier,
            )
            row["formed_this_step_m3"] = formed_this_step
            rows.append(row)
            with out.open("a") as handle:
                handle.write(json.dumps(row, separators=(",", ":")) + "\n")

    cratons.update = wrapped_update
    try:
        for index in range(steps):
            step_world(world, step_myr * 1e6)
            if (index + 1) % 25 == 0:
                print(
                    f"seed={seed} multiplier={multiplier:g} step={index + 1}/{steps} "
                    f"myr={(index + 1) * step_myr:g} "
                    f"runtime_min={(time.perf_counter() - started) / 60:.1f}",
                    flush=True,
                )
    finally:
        cratons.update = original_update
    write_replay_summary(
        out,
        seed=seed,
        steps=steps,
        step_myr=step_myr,
        rows=rows,
        started=started,
        ledger=world.craton_ledger,
        glacier_erosion_multiplier=multiplier,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--multiplier", type=float, choices=(1.0, 2.0), required=True)
    add_replay_arguments(parser, default_steps=250)
    args = parser.parse_args()
    prepare_replay_output(parser, args)
    run(args.seed, args.multiplier, args.steps, args.step_myr, args.snapshot_every, args.out)


if __name__ == "__main__":
    main()
