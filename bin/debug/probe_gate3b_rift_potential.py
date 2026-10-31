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

from app import cratons, world as world_mod
import gate3b_diagnostic_common as shared


def run(seed: int, density: float) -> dict:
    world = world_mod.generate_world(
        seed, node_density=density, cratons_first=True, synthetic_relief=False
    )
    world.simulate_climate_biomes = False
    cratons.ensure_initialized(world)
    components = shared.edge_connected_continental_components(world)

    if not components:
        raise ValueError("generated world has no edge-connected continental component")
    area_m2, plate, selected = components[0]
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
    shared.write_json_report(args.out, report)
    print(json.dumps({"out": str(args.out), **report}, sort_keys=True))


if __name__ == "__main__":
    main()
