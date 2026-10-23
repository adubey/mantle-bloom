#!/usr/bin/env python3
"""Measure divergent-thinning width and rift-ledger closure in a generated world.

The probe wraps the production column-deformation call without changing its inputs or
outputs. For each connected patch of cells that actually loses crust to divergent
deformation, it estimates area-equivalent width as 2 * area / perimeter on the plate mesh.
The result is the spatial support of the current model's active thinning in one step, not
the full sedimentary-basin or tectonic-province width. Checkpoint rows also audit the rift
path subaccounts and the continental-material and craton balances.

Run one replay at a time; do not run alongside unit or stress tests:

    PYTHONHASHSEED=0 PYTHONPATH=backend NUMBA_CACHE_DIR=<private-cache> \\
      backend/.venv/bin/python bin/debug/probe_rfc0012_world_rifts.py \\
      --seed 326494039 --steps 500 --step-myr 1 --every-myr 10 \\
      --out <scratch-output.jsonl>

Raw per-step/component traces stay in private scratch; check in only compact summaries.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from app import continental_ledger, cratons, quad_tectonics, world as world_module
from app.sparse_quad_patch import lattice_points, unpack_cell_keys

EARTH_RADIUS_M = 6_371_000.0


def _component_geometry(plate, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return area and perimeter for edge-connected components of an active cell mask."""
    indices = np.flatnonzero(mask)
    if not len(indices):
        return np.zeros(0), np.zeros(0)

    graph = plate.adjacency()
    rows = np.repeat(np.arange(plate.node_count()), np.diff(graph.offsets))
    adjacency = csr_matrix(
        (np.ones(len(graph.neighbours), dtype=np.int8), (rows, graph.neighbours)),
        shape=(plate.node_count(), plate.node_count()),
    )[indices][:, indices]
    component_count, labels = connected_components(adjacency, directed=False)
    areas = plate.node_areas_m2()[indices]
    component_areas = np.bincount(labels, weights=areas, minlength=component_count)

    keys = plate._keys[indices]
    face, level, i, j = unpack_cell_keys(keys)
    scale = np.left_shift(1, level).astype(float)
    x0, x1 = i / scale, (i + 1) / scale
    y0, y1 = j / scale, (j + 1) / scale
    xm, ym = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    # The eight half-edges match _probe_neighbour_indices' (side, half) order.
    starts_a = np.stack((x0, xm, x1, x1, xm, x0, x0, x0), axis=1).reshape(-1)
    starts_b = np.stack((y0, y0, y0, ym, y1, y1, ym, y0), axis=1).reshape(-1)
    ends_a = np.stack((xm, x1, x1, x1, x0, xm, x0, x0), axis=1).reshape(-1)
    ends_b = np.stack((y0, y0, ym, y1, y1, y1, y0, ym), axis=1).reshape(-1)
    edge_faces = np.repeat(face, 8)
    edge_component = np.repeat(labels, 8)
    probes = plate._probe_neighbour_indices()[indices].reshape(-1)
    neighbor_active = np.zeros(len(probes), dtype=bool)
    present = probes >= 0
    neighbor_active[present] = mask[probes[present]]
    boundary = ~neighbor_active

    start_xyz = lattice_points(edge_faces[boundary], starts_a[boundary], starts_b[boundary], plate._n)
    end_xyz = lattice_points(edge_faces[boundary], ends_a[boundary], ends_b[boundary], plate._n)
    arc_m = np.arccos(np.clip(np.einsum("ij,ij->i", start_xyz, end_xyz), -1.0, 1.0)) * EARTH_RADIUS_M
    component_perimeters = np.bincount(
        edge_component[boundary], weights=arc_m, minlength=component_count
    )
    return component_areas, component_perimeters


def _wrap_deformation(runtime: dict) -> None:
    original = quad_tectonics.deform_columns

    def traced(world, plate, ctx, fields, local_xyz, node_area_m2, *args, **kwargs):
        columns = original(world, plate, ctx, fields, local_xyz, node_area_m2, *args, **kwargs)
        hc_before = np.asarray(fields["crustal_thickness_m"], dtype=float)
        hc_after = np.asarray(columns["crustal_thickness_m"], dtype=float)
        areas = np.broadcast_to(np.asarray(node_area_m2, dtype=float), hc_before.shape)
        loss = np.maximum(hc_before - hc_after, 0.0)
        thinned = ctx.divergent & (loss > 0.0)
        if np.any(thinned):
            component_areas, perimeters = _component_geometry(plate, thinned)
            widths_km = np.divide(
                2.0 * component_areas,
                perimeters,
                out=np.zeros_like(component_areas),
                where=perimeters > 0.0,
            ) / 1000.0
            runtime["components"].extend(
                {"area_m2": float(area), "width_km": float(width)}
                for area, width in zip(component_areas, widths_km, strict=True)
                if area > 0.0
            )
            runtime["divergent_loss_m3"] += float(np.dot(loss[thinned], areas[thinned]))
            runtime["thinned_nodes"] += int(thinned.sum())
        return columns

    quad_tectonics.deform_columns = traced
    runtime["restore"] = lambda: setattr(quad_tectonics, "deform_columns", original)


def _summarize_step(runtime: dict) -> dict:
    components = runtime["components"]
    if not components:
        return {
            "components": 0,
            "active_area_km2": 0.0,
            "area_equivalent_width_km": None,
            "component_width_km_p50": None,
            "component_width_km_p90": None,
            "maximum_component_width_km": None,
            "active_area_share_ge_100km": 0.0,
            "active_area_share_ge_300km": 0.0,
            "component_share_ge_100km": 0.0,
            "component_share_ge_300km": 0.0,
        }
    areas = np.asarray([item["area_m2"] for item in components])
    widths = np.asarray([item["width_km"] for item in components])
    order = np.argsort(widths)
    cumulative = np.cumsum(areas[order])

    def weighted_quantile(q: float) -> float:
        return float(widths[order][np.searchsorted(cumulative, q * cumulative[-1])])

    perimeter_half_m = sum(item["area_m2"] / (item["width_km"] * 1000.0) for item in components)
    return {
        "components": len(components),
        "active_area_km2": float(areas.sum() / 1e6),
        "area_equivalent_width_km": float(areas.sum() / perimeter_half_m / 1000.0),
        "component_width_km_p50": weighted_quantile(0.5),
        "component_width_km_p90": weighted_quantile(0.9),
        "maximum_component_width_km": float(widths.max()),
        "active_area_share_ge_100km": float(areas[widths >= 100.0].sum() / areas.sum()),
        "active_area_share_ge_300km": float(areas[widths >= 300.0].sum() / areas.sum()),
        "component_share_ge_100km": float(np.mean(widths >= 100.0)),
        "component_share_ge_300km": float(np.mean(widths >= 300.0)),
    }


def run(seed: int, steps: int, step_myr: float, every_myr: float, out: Path) -> None:
    world = world_module.generate_world(seed)
    cratons.ensure_initialized(world)
    started = time.perf_counter()
    model_commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    runtime = {"components": [], "divergent_loss_m3": 0.0, "thinned_nodes": 0}
    _wrap_deformation(runtime)
    checkpoint_steps = max(1, round(every_myr / step_myr))
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        with out.open("w", encoding="utf-8") as handle:
            for step in range(1, steps + 1):
                runtime["components"] = []
                runtime["divergent_loss_m3"] = 0.0
                runtime["thinned_nodes"] = 0
                world_module.step_world(world, step_myr * 1e6)
                if step % checkpoint_steps and step != steps:
                    continue
                continental_ledger.assert_closed(world)
                craton = cratons.diagnostics(world)
                craton_ledger = world.craton_ledger
                rift_paths = sum(craton_ledger.get(key, 0.0) for key in cratons.CRATON_RIFT_PATH_ACCOUNTS)
                row = {
                    "seed": seed,
                    "step": step,
                    "myr": world.elapsed_years / 1e6,
                    "model_commit": model_commit,
                    "rift_width": _summarize_step(runtime),
                    "divergent_crust_loss_m3": runtime["divergent_loss_m3"],
                    "divergent_thinned_nodes": runtime["thinned_nodes"],
                    "rifted_craton_m3": craton_ledger.get("rifted_m3", 0.0),
                    "rift_path_sum_m3": rift_paths,
                    "rift_path_residual_m3": rift_paths - craton_ledger.get("rifted_m3", 0.0),
                    "craton_balance_error_m3": craton["craton_balance_error_m3"],
                    "continental_material_balance_error_m3": continental_ledger.balance_error_m3(world),
                    "wall_seconds": time.perf_counter() - started,
                }
                handle.write(json.dumps(row, allow_nan=False) + "\n")
                handle.flush()
                print(json.dumps({key: row[key] for key in ("step", "myr", "rift_width", "wall_seconds")}), flush=True)
    finally:
        runtime["restore"]()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=326494039)
    parser.add_argument("--steps", type=int, default=500)
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--every-myr", type=float, default=10.0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    run(args.seed, args.steps, args.step_myr, args.every_myr, args.out)


if __name__ == "__main__":
    main()
