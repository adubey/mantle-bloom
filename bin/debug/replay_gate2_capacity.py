#!/usr/bin/env python3
"""Compare Gate 2 strength alternatives on repeatable Gate 1 rift-cut candidates.

Recreate the cratons-first baseline for seeds 7, 11, and 2024 at 2 Myr steps. At each
production cut crossing craton, reconstruct the event-local spatial force field, solve the
free-edge shell diagnostic once on the clipped plate mesh, then score that same stress and
cut band against physical column-capacity scenarios and the labeled belt-weakening rule.
Capacity scores do not feed back into the world, so both alternatives see identical loads,
cut events, and states. This is a diagnostic replay, not an A/B of production behavior.

Run one seed at a time, never alongside unit/stress tests or another replay:

    PYTHONHASHSEED=0 PYTHONPATH=backend backend/.venv/bin/python \
        bin/debug/replay_gate2_capacity.py --seed 7 --out /tmp/issue487-g2-seed7.json

Use ``--years 1e8`` for a short smoke replay. The full baseline is 1e9 years.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "backend"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import validate_thin_shell_plate_meshes as mesh_validation  # noqa: E402
import prototype_mincut_split as mincut  # noqa: E402
from app import cratons, geometry, initial_terrain, lithosphere, merge_split, surface_fields, thin_shell, torque, world as world_mod  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402

STEP_YEARS = 2e6
EFFECTIVE_MODULUS_PA = 70e9
PHENOMENOLOGICAL_WEAKENING = 0.5
PHENOMENOLOGICAL_REFERENCE_FACTORS = (0.5, 1.0, 2.0)
GATE1_EXPECTED_CUT_COUNTS = {7: 144, 11: 111, 2024: 125}
GATE1_EXPECTED_CRATON_CROSSING_CUTS = {7: 44, 11: 42, 2024: 44}
GATE1_BASELINE_COMMIT = "0118752d0ffe75060464204d5b94f81a60fe9a07"


def _physical_capacities() -> dict[str, list[float]]:
    data = json.loads((ROOT / "analysis/issue487/gate2-column-envelope.json").read_text())
    summaries = data["sensitivity"]["summary_by_crust_type"]
    result = {}
    for crust_type in ("continental", "oceanic_young_ridge", "oceanic_mature"):
        key = "continental" if crust_type == "continental" else "oceanic_young_ridge" if crust_type == "oceanic_young_ridge" else "oceanic_mature"
        marginal = summaries[key]["marginal_by_parameter"]["base_temperature_c"]["1330.0"]
        result[crust_type] = [float(value) * 1e12 for value in marginal["capacity_tn_per_m_p10_median_p90"]]
    return result


def _event_force_field(world, plate):
    spacing = line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    cache = getattr(world, "boundary_search_cache", None)
    neighbours = plate.get_neighbours(world.plates)
    inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
    _convergent, divergent, _transform, contested = torque.classify_boundary_nodes(
        plate, neighbours, inputs, reach, cache
    )
    collision = torque.continental_contact_mask(plate, inputs, contested)
    subducting = torque.subducting_boundary_mask(plate, inputs, reach)
    craton_contact = torque.contact_craton_m(plate, neighbours, inputs)
    jam = torque.contact_jam_factor(plate, world, neighbours, inputs, reach)
    return torque.gather_per_node_force_field(
        plate, world, inputs, spacing, subducting, divergent, collision, craton_contact, jam
    )


def _topology_audit(plate, node_count: int, triangles: np.ndarray) -> dict:
    """Compare exact source-cell components with components present in accepted elements."""
    adjacency = plate.adjacency()
    source_graph = csr_matrix(
        (np.ones(len(adjacency.neighbours)), adjacency.neighbours, adjacency.offsets),
        shape=(node_count, node_count),
    )
    source_count, source_labels = connected_components(source_graph, directed=False)
    edges = np.concatenate((triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]))
    mesh_graph = coo_matrix(
        (np.ones(2 * len(edges)),
         (np.r_[edges[:, 0], edges[:, 1]], np.r_[edges[:, 1], edges[:, 0]])),
        shape=(node_count, node_count),
    ).tocsr()
    _, mesh_labels = connected_components(mesh_graph, directed=False)
    active = np.zeros(node_count, dtype=bool)
    active[np.unique(triangles)] = True
    source_components = []
    for component_id in range(source_count):
        nodes = source_labels == component_id
        _, mesh_sizes = np.unique(mesh_labels[nodes & active], return_counts=True)
        source_components.append({
            "source_component_id": int(component_id),
            "source_node_count": int(nodes.sum()),
            "meshed_node_count": int((nodes & active).sum()),
            "unmeshed_node_count": int((nodes & ~active).sum()),
            "mesh_component_sizes": sorted((int(size) for size in mesh_sizes), reverse=True),
        })
    mesh_components = []
    for component_id in np.unique(mesh_labels[active]):
        nodes = (mesh_labels == component_id) & active
        mesh_components.append({
            "mesh_component_id": int(component_id),
            "mesh_node_count": int(nodes.sum()),
            "source_component_ids": [int(label) for label in np.unique(source_labels[nodes])],
        })
    return {
        "source_component_count": int(source_count),
        "source_components": source_components,
        "mesh_components": mesh_components,
    }


def _capacity_scores(resultant: np.ndarray, band: np.ndarray, continental: np.ndarray,
                     young_oceanic: np.ndarray, belt_marked: np.ndarray, capacities: dict) -> dict:
    in_band = band
    scores: dict[str, object] = {}
    physical_cases = {}
    for class_name, class_mask in (
        ("continental", continental),
        ("oceanic_young_ridge", young_oceanic & ~continental),
        ("oceanic_mature", ~young_oceanic & ~continental),
    ):
        selected = in_band & class_mask
        if not np.any(selected):
            continue
        caps = capacities[class_name]
        physical_cases[class_name] = {
            "nodes": int(selected.sum()),
            "capacity_p10_median_p90_n_per_m": caps,
            "failed_fraction_p10_median_p90": [
                float(np.mean(resultant[selected] >= cap)) for cap in caps
            ],
        }
    scores["physical"] = physical_cases
    weak = belt_marked & in_band
    reference = np.empty(len(resultant), dtype=float)
    for class_name, class_mask in (
        ("continental", continental),
        ("oceanic_young_ridge", young_oceanic & ~continental),
        ("oceanic_mature", ~young_oceanic & ~continental),
    ):
        reference[class_mask] = capacities[class_name][1]
    phenom = []
    for factor in PHENOMENOLOGICAL_REFERENCE_FACTORS:
        capacity = reference * factor * np.where(belt_marked, 1.0 - PHENOMENOLOGICAL_WEAKENING, 1.0)
        phenom.append({
            "multiplier_of_physical_median_reference": factor,
            "failed_fraction_in_cut_band": float(np.mean(resultant[in_band] >= capacity[in_band])),
            "failed_fraction_in_marked_belt": (
                float(np.mean(resultant[weak] >= capacity[weak])) if np.any(weak) else None
            ),
            "marked_belt_nodes_in_cut_band": int(weak.sum()),
        })
    scores["phenomenological"] = {
        "weak_zone": "existing belt_formed_years marks; lambda=0.5; class-matched physical median used only as the reference scale",
        "reference_capacity_sweep": phenom,
    }
    return scores


def _score_event(world, plate, points, cut_normal, kind, capacities) -> dict | None:
    points = np.asarray(points, dtype=float)
    normal = np.asarray(cut_normal, dtype=float)
    spacing = line_spacing_rad(world.node_density)
    band_chord = 2.0 * np.sin(0.5 * merge_split.FAILED_RIFT_BAND_MULT * spacing)
    graph_i, graph_j, _length = mincut._graph(points, spacing)
    flow_mask = points @ normal > 0.0
    band = mincut._boundary_band(points, graph_i, graph_j, flow_mask, band_chord)
    craton = plate.collect("craton_crust_m") > 0.0
    if not np.any(band & craton):
        return None

    field = _event_force_field(world, plate)
    source_points = np.asarray(field.points_xyz, dtype=float)
    mesh_method = "clipped_gnomonic_delaunay"
    try:
        candidate_triangles, _rejected, _total = mesh_validation._candidate_triangles(plate, source_points)
    except ValueError as error:
        if "does not fit in one tangent hemisphere" not in str(error):
            raise
        # Broad plates may exceed one gnomonic patch; clip their spherical hull by the true
        # territory outline so the diagnostic still uses the original curved footprint.
        hull = thin_shell._spherical_surface_triangles(source_points)
        centroids = geometry.normalize(source_points[hull].mean(axis=1))
        candidate_triangles = hull[plate.contains_batch(centroids)]
        mesh_method = "spherical_hull_clipped_to_plate_outline"
    if len(candidate_triangles) == 0:
        return {"skipped": "no clipped mesh triangles", "plate": plate.plate_id, "myr": world.elapsed_years / 1e6}
    nondegenerate = []
    for triangle in candidate_triangles:
        element = source_points[triangle]
        radial = geometry.normalize(element.sum(axis=0))
        edge = element[1] - element[0]
        tangent = edge - np.dot(edge, radial) * radial
        tangent_norm = float(np.linalg.norm(tangent))
        if tangent_norm <= np.finfo(float).eps:
            continue
        axis_x = tangent / tangent_norm
        axis_y = np.cross(radial, axis_x)
        xy = np.column_stack((element @ axis_x, element @ axis_y)) * lithosphere.PLANET_RADIUS_M
        twice_area = (
            (xy[1, 0] - xy[0, 0]) * (xy[2, 1] - xy[0, 1])
            - (xy[1, 1] - xy[0, 1]) * (xy[2, 0] - xy[0, 0])
        )
        if abs(twice_area) * 0.5 > np.finfo(float).eps * lithosphere.PLANET_RADIUS_M**2:
            nondegenerate.append(triangle)
    candidate_triangles = np.asarray(nondegenerate, dtype=np.int64)
    if len(candidate_triangles) == 0:
        return {"skipped": "all clipped mesh triangles are degenerate", "plate": plate.plate_id, "myr": world.elapsed_years / 1e6}
    topology_audit = _topology_audit(plate, len(source_points), candidate_triangles)
    if topology_audit["source_component_count"] > 1:
        return {
            "skipped": "source plate adjacency has disconnected cell regions",
            "plate": int(plate.plate_id),
            "myr": float(world.elapsed_years / 1e6),
            "mesh_method": mesh_method,
            "mesh_node_count": int(len(np.unique(candidate_triangles))),
            "mesh_triangle_count": int(len(candidate_triangles)),
            "topology_audit": topology_audit,
        }
    used = np.unique(candidate_triangles)
    remap = np.full(len(source_points), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    triangles = remap[candidate_triangles]
    band = band[used]
    if not np.any(band):
        return {"skipped": "cut band has no nodes on clipped mesh", "plate": plate.plate_id, "myr": world.elapsed_years / 1e6}

    hc = plate.collect("crustal_thickness_m")[used]
    hm = plate.collect("mantle_lithosphere_thickness_m")[used]
    codes = plate.collect("crust_type_code")[used]
    continental = np.asarray(
        __import__("app.elevation_lines", fromlist=["effective_is_continental_from_codes"])
        .effective_is_continental_from_codes(codes, plate.crust_type == "continental"), dtype=bool
    )
    # Active divergent oceanic nodes are the only runtime proxy for young-ridge columns;
    # this remains a coarse classification, as the saved state has no oceanic crust-age field.
    neighbours = plate.get_neighbours(world.plates)
    reach = torque.boundary_force_reach_rad(spacing)
    inputs = torque.gather_boundary_force_inputs(
        plate, neighbours, spacing, reach, getattr(world, "boundary_search_cache", None)
    )
    _conv, divergent, _trans, _contested = torque.classify_boundary_nodes(
        plate, neighbours, inputs, reach, getattr(world, "boundary_search_cache", None)
    )
    young_oceanic = divergent[used] & ~continental
    belt = plate.collect("belt_formed_years")[used] != surface_fields.BELT_UNMARKED_YEARS
    try:
        result = thin_shell.solve_thin_shell(
            source_points[used], field.areas_m2[used], field.total_n[used],
            thickness_m=np.maximum(hc + hm, 1.0), young_modulus_pa=EFFECTIVE_MODULUS_PA,
            extensional_capacity_resultant_n_per_m=np.full(len(used), 1e10),
            mesh_triangles=triangles,
        )
    except ValueError as error:
        if "disconnected element components" not in str(error):
            raise
        return {
            "skipped": "thin-shell mesh is disconnected",
            "mesh_error": str(error),
            "plate": int(plate.plate_id),
            "myr": float(world.elapsed_years / 1e6),
            "mesh_node_count": int(len(used)),
            "mesh_triangle_count": int(len(triangles)),
            "mesh_method": mesh_method,
            "topology_audit": topology_audit,
        }
    resultant = result.max_incremental_extensional_resultant_n_per_m
    scores = _capacity_scores(resultant, band, continental, young_oceanic, belt, capacities)
    in_band_craton = band & (plate.collect("craton_crust_m")[used] > 0.0)
    return {
        "kind": kind,
        "myr": float(world.elapsed_years / 1e6),
        "plate": int(plate.plate_id),
        "node_count": int(len(source_points)),
        "mesh_node_count": int(len(used)),
        "mesh_method": mesh_method,
        "topology_audit": topology_audit,
        "band_node_count": int(band.sum()),
        "band_craton_node_count": int(in_band_craton.sum()),
        "band_max_resultant_n_per_m": float(np.max(resultant[band])),
        "band_p50_resultant_n_per_m": float(np.median(resultant[band])),
        "band_p95_resultant_n_per_m": float(np.quantile(resultant[band], 0.95)),
        "solver_stop_code": int(result.solver_stop_code),
        "equilibrium_residual_n": float(result.equilibrium_residual_n),
        "equilibrium_residual_over_load": float(
            result.equilibrium_residual_n / max(np.linalg.norm(result.equilibrated_force_n), 1.0)
        ),
        "capacity_scores": scores,
    }


def _run(seed: int, years: float, out_path: Path) -> dict:
    for name, value in __import__("compare_cratons")._DEFAULTS.items():
        setattr(cratons, name, value)
    cratons.CRATON_FORMATION_MYR = float("inf")
    capacities = _physical_capacities()
    events: list[dict] = []
    original_cut = cratons.record_rift_cut
    original_lowland_relief = initial_terrain.apply_lowland_spherical_relief

    def tapped(world, plate, points, cut_normal, band_sin, kind):
        started = time.perf_counter()
        row = _score_event(world, plate, points, cut_normal, kind, capacities)
        if row is not None:
            row["event_index"] = len(events)
            row["diagnostic_wall_s"] = time.perf_counter() - started
            events.append(row)
            out_path.write_text(json.dumps({"seed": seed, "years": years, "events": events}, indent=2) + "\n")
            print(json.dumps({"seed": seed, "event": len(events), "myr": row.get("myr"),
                              "plate": row.get("plate"), "diag_s": row["diagnostic_wall_s"]}), flush=True)
        original_cut(world, plate, points, cut_normal, band_sin, kind)

    cratons.record_rift_cut = tapped
    # The repeatable Gate 1 run predates #490; suppress that later initial-terrain addition
    # so this replay starts from the exact accepted baseline state.
    initial_terrain.apply_lowland_spherical_relief = lambda _world: None
    started = time.perf_counter()
    try:
        world = world_mod.generate_world(seed, cratons_first=True)
        cratons.ensure_initialized(world)
        while world.elapsed_years < years - 1.0:
            world_mod.step_world(world, min(STEP_YEARS, years - world.elapsed_years))
            if world.steps_taken % 500 == 0:
                print(json.dumps({"seed": seed, "myr": world.elapsed_years / 1e6,
                                  "events_scored": len(events), "wall_s": time.perf_counter() - started}), flush=True)
    finally:
        cratons.record_rift_cut = original_cut
        initial_terrain.apply_lowland_spherical_relief = original_lowland_relief
    return {
        "seed": seed,
        "years": years,
        "step_years": STEP_YEARS,
        "craton_formation_off": True,
        "gate1_baseline_commit": GATE1_BASELINE_COMMIT,
        "events_scored": len(events),
        "gate1_craton_crossing_cut_count_expected": GATE1_EXPECTED_CRATON_CROSSING_CUTS[seed],
        "gate1_craton_crossing_cut_count_matches": len(events) == GATE1_EXPECTED_CRATON_CROSSING_CUTS[seed],
        "production_cut_count": int(round(
            world.craton_ledger.get("split_cuts", 0.0)
            + world.craton_ledger.get("failed_rift_cuts", 0.0)
        )),
        "gate1_cut_count_expected": GATE1_EXPECTED_CUT_COUNTS[seed],
        "gate1_cut_count_matches": int(round(
            world.craton_ledger.get("split_cuts", 0.0)
            + world.craton_ledger.get("failed_rift_cuts", 0.0)
        )) == GATE1_EXPECTED_CUT_COUNTS[seed],
        "wall_seconds": time.perf_counter() - started,
        "physical_capacity_n_per_m_by_class_p10_median_p90": capacities,
        "events": events,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, choices=(7, 11, 2024), required=True)
    parser.add_argument("--years", type=float, default=1e9)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    report = _run(args.seed, args.years, args.out)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in ("seed", "years", "events_scored", "wall_seconds")} ), flush=True)


if __name__ == "__main__":
    main()
