#!/usr/bin/env python3
"""Compare two source-backed river velocity estimates at a saved checkpoint.

Reads the exact hydrology cache retained in a v5 checkpoint and summarizes reach length,
discharge, drainage area, slope, velocity, travel time and steady-flow inventory. The Jobson
comparison currently sets mean annual discharge equal to the checkpoint's interval-mean rate;
that is an explicit sensitivity assumption, not a validated annual mean. No simulation step is
run and no world state is changed.

Example:

    backend/.venv/bin/python bin/debug/measure_river_velocity_candidates.py \\
        <checkpoint.mbworld> --backend <source-worktree>/backend --out <summary.json>
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("save", type=Path)
parser.add_argument("--backend", type=Path, default=Path("backend"))
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()

backend = args.backend.resolve()
sys.path.insert(0, str(backend))
from app import elevation_lines, geometry, persistence, plates  # noqa: E402
from app.hydroclimate import hydrology  # noqa: E402

payload = args.save.read_bytes()
world = persistence.load_world_bytes(payload)
fields = world.hydrology_cache
if fields is None:
    raise SystemExit("checkpoint has no saved hydrology cache")
if not hasattr(fields, "discharge_m3_s") or fields.discharge_m3_s is None:
    raise SystemExit("checkpoint hydrology cache has no physical discharge")

areas = plates.collect_all_accounting_areas_m2(
    fields.plates_in_order, elevation_lines.line_spacing_rad(world.node_density)
)
target = np.asarray(fields.flow_target, dtype=np.int64)
q = np.asarray(fields.discharge_m3_s, dtype=float)
n = len(target)
if any(len(values) != n for values in (areas, q, fields.elevation, fields.is_ocean)):
    raise SystemExit("saved hydrology cache is not aligned with checkpoint surface nodes")

node_area = np.asarray(areas, dtype=float)
has_target = target >= 0
edge_sources = np.flatnonzero(has_target)
graph = coo_matrix(
    (np.ones(edge_sources.size), (edge_sources, target[has_target])), shape=(n, n)
).tocsr()
component_count, component = connected_components(graph, directed=True, connection="strong")
component_area = np.bincount(component, weights=node_area, minlength=component_count)

# Condense directed flow cycles before accumulation. Within a cycle, every node sees the
# cycle's whole contributing area; downstream catchments then receive that area once.
component_edges = np.column_stack((component[edge_sources], component[target[has_target]]))
component_edges = np.unique(component_edges[component_edges[:, 0] != component_edges[:, 1]], axis=0)
outgoing: list[list[int]] = [[] for _ in range(component_count)]
indegree = np.zeros(component_count, dtype=np.int64)
for source_component, target_component in component_edges:
    outgoing[int(source_component)].append(int(target_component))
    indegree[int(target_component)] += 1
ready = list(np.flatnonzero(indegree == 0))
cursor = 0
while cursor < len(ready):
    source_component = ready[cursor]
    cursor += 1
    for target_component in outgoing[source_component]:
        component_area[target_component] += component_area[source_component]
        indegree[target_component] -= 1
        if indegree[target_component] == 0:
            ready.append(target_component)
if cursor != component_count:
    raise SystemExit("cycle condensation did not produce an acyclic routing graph")
contributing_area = component_area[component]

safe_target = np.where(target >= 0, target, 0)
length_m = (
    geometry.angular_distance(fields.points, fields.points[safe_target])
    * hydrology.PLANET_RADIUS_KM * 1000.0
)
length_m = np.maximum(length_m, 1.0)
slope = hydrology._slope_to_flow_target(fields.points, fields.elevation, target)
lake_depth = np.asarray(getattr(fields, "lake_depth", np.zeros(n)), dtype=float)
lake_nodes = lake_depth > hydrology.LAKE_MIN_VISIBLE_DEPTH_M
eligible = (~np.asarray(fields.is_ocean, dtype=bool)) & ~lake_nodes & (target >= 0)
positive = eligible & np.isfinite(q) & (q > 0.0)

# Raymond et al. (2012) SI hydraulic-geometry fit: U = 0.19 Q^0.29.
raymond_u = 0.19 * np.power(q[positive], 0.29)
raymond_tau = length_m[positive] / raymond_u

# Jobson (1996) national peak-tracer regression. With Qa=Q, Q'a=1; this is a
# sensitivity to a missing long-term mean-flow input, not a calibrated mean velocity.
qa = q[positive]
da = np.maximum(contributing_area[positive], 1.0)
da_prime = np.power(da, 1.25) * np.sqrt(9.80665) / qa
jobson_u = 0.094 + 0.0143 * np.power(da_prime, 0.919) * np.power(1.0, -0.469) * np.power(
    np.maximum(slope[positive], 0.0), 0.159
) * (q[positive] / da)
jobson_tau = length_m[positive] / jobson_u


def summary(values: np.ndarray) -> dict[str, float | list[float]]:
    return {
        "min": float(np.min(values)),
        "p10": float(np.percentile(values, 10)),
        "median": float(np.median(values)),
        "p90": float(np.percentile(values, 90)),
        "max": float(np.max(values)),
    }


def inventory(q_values: np.ndarray, tau: np.ndarray) -> float:
    return float(np.sum(q_values * tau))


def reach_summary(mask: np.ndarray) -> dict[str, float | int | dict]:
    selected_q = q[mask]
    selected_l = length_m[mask]
    selected_areas = contributing_area[mask]
    selected_slope = slope[mask]
    ray = 0.19 * np.power(selected_q, 0.29)
    da_p = np.power(np.maximum(selected_areas, 1.0), 1.25) * np.sqrt(9.80665) / selected_q
    job = 0.094 + 0.0143 * np.power(da_p, 0.919) * np.power(
        np.maximum(selected_slope, 0.0), 0.159
    ) * (selected_q / np.maximum(selected_areas, 1.0))
    ray_tau = selected_l / ray
    job_tau = selected_l / job
    return {
        "reach_count": int(mask.sum()),
        "discharge_m3_s": summary(selected_q),
        "edge_length_m": summary(selected_l),
        "drainage_area_m2": summary(selected_areas),
        "slope": summary(selected_slope),
        "raymond_2012": {
            "velocity_m_s": summary(ray),
            "residence_time_s": summary(ray_tau),
            "steady_inventory_m3": inventory(selected_q, ray_tau),
        },
        "jobson_1996_peak_tracer_with_qa_eq_q_sensitivity": {
            "velocity_m_s": summary(job),
            "residence_time_s": summary(job_tau),
            "steady_inventory_m3": inventory(selected_q, job_tau),
            "fraction_outside_reported_0p01_to_1p51_m_s": float(
                np.mean((job < 0.01) | (job > 1.51))
            ),
        },
    }


all_reaches = reach_summary(positive)
lake_outlets = positive & lake_nodes[safe_target]
result = {
    "backend_commit": subprocess.run(
        ["git", "rev-parse", "HEAD"], cwd=backend, check=True, capture_output=True, text=True
    ).stdout.strip(),
    "checkpoint_sha256": hashlib.sha256(payload).hexdigest(),
    "elapsed_years": float(world.elapsed_years),
    "node_count": n,
    "eligible_routed_edges": int(eligible.sum()),
    "positive_discharge_edges": int(positive.sum()),
    "zero_or_nonpositive_discharge_edges": int((eligible & ~positive).sum()),
    "routing_cycle_components": int(np.sum(np.bincount(component) > 1)),
    "jobson_assumption": "checkpoint interval-mean discharge equals mean annual discharge (Qa=Q)",
    "positive_flow_reaches": all_reaches,
    "positive_flow_reaches_entering_lake": reach_summary(lake_outlets) if np.any(lake_outlets) else None,
}
args.out.parent.mkdir(parents=True, exist_ok=True)
args.out.write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
