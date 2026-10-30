#!/usr/bin/env python3
"""Select the preregistered Gate 3b continental component and check mesh eligibility.

Usage:
    PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
        backend/.venv/bin/python bin/debug/probe_gate3b_rift_potential.py \
        --out analysis/issue699/rift-potential-seed11-density4-preflight.json

This selection-only preflight stops before force scoring if the existing whole-production-plate
shell mesh cannot represent the largest connected continental component.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from app import cratons, elevation_lines, world as world_mod


def run(seed: int, density: float) -> dict:
    world = world_mod.generate_world(
        seed, node_density=density, cratons_first=True, synthetic_relief=False
    )
    world.simulate_climate_biomes = False
    cratons.ensure_initialized(world)
    components: list[tuple[float, object, np.ndarray]] = []
    for plate in world.plates:
        continental = elevation_lines.effective_is_continental_from_codes(
            plate.collect("crust_type_code"), plate.crust_type == "continental"
        )
        indices = np.flatnonzero(continental)
        if not len(indices):
            continue
        adjacency = plate.adjacency()
        rows = np.repeat(np.arange(plate.node_count(), dtype=np.int64), np.diff(adjacency.offsets))
        cols = adjacency.neighbours.astype(np.int64, copy=False)
        local_index = np.full(plate.node_count(), -1, dtype=np.int64)
        local_index[indices] = np.arange(len(indices), dtype=np.int64)
        keep = continental[rows] & continental[cols]
        graph = csr_matrix(
            (np.ones(int(keep.sum())), (local_index[rows[keep]], local_index[cols[keep]])),
            shape=(len(indices), len(indices)),
        )
        count, labels = connected_components(graph, directed=False)
        area_by_node = plate.node_areas_m2()[indices]
        for label in range(count):
            selected = indices[labels == label]
            area = float(area_by_node[labels == label].sum())
            components.append((area, plate, selected))

    if not components:
        raise ValueError("generated world has no edge-connected continental component")
    area_m2, plate, selected = max(components, key=lambda row: row[0])
    result = {
        "schema_version": 1,
        "title": "RFC 0010 Gate 3b rift-potential accounting preflight",
        "seed": seed,
        "node_density": density,
        "elapsed_years": float(world.elapsed_years),
        "initial_state": "Generated t=0; cratons-first; climate/biomes disabled; no startup step",
        "plate_selection": "Largest edge-connected continental component by represented area, before inspecting response",
        "plate_id": int(plate.plate_id),
        "production_plate_node_count": int(plate.node_count()),
        "production_plate_area_km2": float(plate.node_areas_m2().sum() / 1e6),
        "selected_continental_component_node_count": int(len(selected)),
        "selected_continental_component_area_km2": area_m2 / 1e6,
        "other_continental_component_count": int(len(components) - 1),
    }
    if len(selected) != plate.node_count():
        result.update({
            "eligibility": "ineligible_before_force_transfer",
            "reason": (
                "The selected continental component is a proper subset of its production plate. "
                "Existing coarse-shell mesh validation uses whole-plate containment and cannot "
                "honor this component boundary; no force, capacity, or event score was computed."
            ),
            "next_requirement": (
                "Build shell mesh and force aggregation that preserve the selected continental-"
                "component boundary; rerun frozen accounting checks before the persistence arm."
            ),
            "history_eligible": False,
            "interpretation": "Diagnostic eligibility failure, not evidence for or against a persistent rift-potential event.",
        })
    else:
        result.update({
            "eligibility": "component_matches_production_plate",
            "reason": "Selection matches the whole production plate; force-transfer accounting may proceed.",
        })
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--density", type=float, default=4.0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = run(args.seed, args.density)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"out": str(args.out), **report}, sort_keys=True))


if __name__ == "__main__":
    main()
