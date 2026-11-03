#!/usr/bin/env python3
"""Issue #758: list every below-sea-level body in saved worlds, ocean or not.

For each connected below-sea-level component (the components `connected_ocean_mask` labels),
reports node count, area, share of the largest component, oceanic-crust share by area,
median rock age (`node_created_years`), median floor elevation, its sill (the lowest level at
which it would join the largest component, by a minimax flood over the hydrology graph) and
whether `connected_ocean_mask` currently calls it ocean. Used to place the absolute-area rule
for an enclosed ocean basin in the gap between real lakes and ocean-scale bodies.

Run from `backend/` so the worktree's environment is active:

    .venv/bin/python ../bin/debug/survey_enclosed_basins.py <save> [<save> ...] \\
        --min-nodes 20 --out <scratch>/survey.json
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import sys
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("saves", type=Path, nargs="+")
ap.add_argument("--backend", type=Path, default=Path(__file__).resolve().parents[2] / "backend")
ap.add_argument("--min-nodes", type=int, default=20)
ap.add_argument("--out", type=Path, required=True)
args = ap.parse_args()
sys.path.insert(0, str(args.backend.resolve()))

from app import elevation_lines, persistence, plates  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402
from app.hydroclimate import hydrology  # noqa: E402


def sill_levels(elevation: np.ndarray, neighbor_idx: np.ndarray, seed: np.ndarray) -> np.ndarray:
    """Minimax path elevation from `seed` to every node: the water level at which each node
    joins the seed's body."""
    n = len(elevation)
    level = np.full(n, np.inf)
    heap = [(float(elevation[i]), int(i)) for i in np.flatnonzero(seed)]
    for lv, i in heap:
        level[i] = lv
    heapq.heapify(heap)
    while heap:
        lv, i = heapq.heappop(heap)
        if lv > level[i]:
            continue
        for j in neighbor_idx[i]:
            cand = max(lv, float(elevation[j]))
            if cand < level[j]:
                level[j] = cand
                heapq.heappush(heap, (cand, int(j)))
    return level


def survey(path: Path) -> dict:
    data = path.read_bytes()
    world = persistence.load_world_bytes(data)
    points, order = plates.gather_node_positions(world.plates)
    elevation = plates.collect_all_elevation(order)
    areas = plates.collect_all_accounting_areas_m2(order, line_spacing_rad(world.node_density))
    oceanic = np.concatenate([
        ~elevation_lines.effective_is_continental_from_codes(p.collect("crust_type_code"), p.crust_type == "continental")
        for p in order
    ])
    created = plates.collect_all_node_created_years(order)
    age_myr = (world.elapsed_years - created) / 1e6
    lake_depth = plates.collect_all_lake_depth(order)
    sea = world.sea_level_m
    neighbor_idx = hydrology._build_neighbor_graph(points)
    below = elevation <= sea
    labels = hydrology._below_sea_components(below, neighbor_idx)
    is_ocean = hydrology.connected_ocean_mask(points, elevation, sea, neighbor_idx)
    counts = np.bincount(labels[below])
    largest = int(counts.argmax())
    sill = sill_levels(elevation, neighbor_idx, labels == largest)
    bodies = []
    for lab in np.flatnonzero(counts >= args.min_nodes):
        m = below & (labels == lab)
        a = areas[m]
        bodies.append({
            "nodes": int(m.sum()),
            "area_km2": float(a.sum() / 1e6),
            "share_of_largest": float(m.sum() / counts[largest]),
            "oceanic_area_share": float(a[oceanic[m]].sum() / a.sum()),
            "median_rock_age_myr": float(np.median(age_myr[m])),
            "p10_rock_age_myr": float(np.percentile(age_myr[m], 10)),
            "median_floor_m": float(np.median(elevation[m])),
            "floor_below_sea_m": float(sea - np.median(elevation[m])),
            "sill_above_sea_m": float(sill[m].min() - sea) if lab != largest else 0.0,
            "lake_volume_km3": float(np.dot(a, lake_depth[m]) / 1e9),
            "volume_to_sea_level_km3": float(np.dot(a, sea - elevation[m]) / 1e9),
            "largest": bool(lab == largest),
            "is_ocean_now": bool(is_ocean[m].all()),
        })
    bodies.sort(key=lambda b: -b["nodes"])
    return {
        "save": path.name,
        "sha256": hashlib.sha256(data).hexdigest(),
        "elapsed_myr": world.elapsed_years / 1e6,
        "sea_level_m": sea,
        "nodes": int(len(elevation)),
        "mean_node_area_km2": float(areas.mean() / 1e6),
        "bodies": bodies,
    }


results = []
for save in args.saves:
    r = survey(save)
    results.append(r)
    print(f"{r['save']}  {r['elapsed_myr']:.1f} Myr  sea {r['sea_level_m']:.0f} m  {r['nodes']} nodes")
    for b in r["bodies"][:12]:
        print(
            f"  {b['nodes']:7d} n {b['area_km2']:12.0f} km2 {100 * b['share_of_largest']:6.2f}%  "
            f"oc {b['oceanic_area_share']:.2f} age {b['median_rock_age_myr']:6.1f} "
            f"depth {b['floor_below_sea_m']:6.0f} sill {b['sill_above_sea_m']:6.0f} "
            f"ocean={b['is_ocean_now']}"
        )
args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(results, indent=1))
