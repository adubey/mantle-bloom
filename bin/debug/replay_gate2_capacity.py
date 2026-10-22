#!/usr/bin/env python3
"""Compare Gate 2 strength alternatives on repeatable Gate 1 rift-cut candidates.

Replay the pinned Gate 1 cratons-first model for seeds 7, 11, and 2024 at 2 Myr steps. The
backend source must be the Gate 1 baseline commit with this diagnostic and the thin-shell
per-node force-field helpers overlaid; running against later production model code can change
the candidate population. Always check the Gate 1 count-match fields in the output. At each
production cut crossing craton, reconstruct the event-local spatial force field, solve the
free-edge shell diagnostic once on the clipped plate mesh, then score that same stress and
cut band against physical column-capacity scenarios and the labeled belt-weakening rule.
Capacity scores do not feed back into the world, so both alternatives see identical loads,
cut events, and states. This is a diagnostic replay, not an A/B of production behavior.

Run at most two seeds at once, never alongside unit/stress tests:

    PYTHONHASHSEED=0 PYTHONPATH=backend backend/.venv/bin/python \
        bin/debug/replay_gate2_capacity.py --seed 7 --out /tmp/issue487-g2-seed7.json

Use ``--years 1e8`` for a short smoke replay. The full baseline is 1e9 years.
"""

from __future__ import annotations

import argparse
import functools
import json
import math
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
from diagnose_gate2_column_envelope import Scenario, evaluate_case  # noqa: E402
from app import cratons, geometry, initial_terrain, lithosphere, merge_split, surface_fields, thin_shell, torque, world as world_mod  # noqa: E402
from app import orogeny  # noqa: E402
from app.elevation_lines import line_spacing_rad  # noqa: E402

STEP_YEARS = 2e6
EFFECTIVE_MODULUS_PA = 70e9
MIN_DIAGNOSTIC_BAND_COVERAGE = 0.95
MIN_REPAIRED_TRIANGLE_FRACTION = 0.99
PHENOMENOLOGICAL_WEAKENING = 0.5
PHENOMENOLOGICAL_REFERENCE_FACTORS = (0.5, 1.0, 2.0)
GATE1_EXPECTED_CUT_COUNTS = {7: 144, 11: 111, 2024: 125}
GATE1_EXPECTED_CRATON_CROSSING_CUTS = {7: 44, 11: 42, 2024: 44}
GATE1_BASELINE_COMMIT = "0118752d0ffe75060464204d5b94f81a60fe9a07"
COLUMN_SCENARIO = Scenario(
    base_temperature_c=orogeny.LITHOSPHERE_BASE_TEMPERATURE_C,
    strain_rate_s=1e-15,
    pore_pressure_ratio=0.0,
    pore_pressure_mode="hydrostatic_water",
    k0=0.8,
    cohesion_pa=20e6,
    friction_angle_deg=math.degrees(math.atan(0.6)),
    lower_crust_law="plagioclase_an75",
    mantle_law="dry_olivine",
)
COLUMN_DZ_M = 125.0


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


@functools.lru_cache(maxsize=100_000)
def _column_capacity(crust_type: str, hc_m: float, hm_m: float, lag_c: float) -> float:
    """Evaluate one local column; round cache keys below in the caller for tractability."""
    result = evaluate_case(
        crust_type,
        COLUMN_SCENARIO,
        COLUMN_DZ_M,
        column_crust_thickness_m=hc_m,
        column_mantle_lithosphere_thickness_m=hm_m,
        moho_thermal_lag_c=lag_c,
    )
    return float(result["integrated_additional_strength_n_per_m"])


def _spatial_capacity_scores(
    resultant: np.ndarray,
    in_band: np.ndarray,
    continental: np.ndarray,
    young_oceanic: np.ndarray,
    belt_marked: np.ndarray,
    hc_m: np.ndarray,
    hm_m: np.ndarray,
    moho_lag_c: np.ndarray,
    craton_hc_m: np.ndarray,
) -> dict:
    """Compare event-local N_ext with mapped Hc/Hm column capacities, without fitting labels."""
    classes = (
        ("continental", continental),
        ("oceanic_young_ridge", young_oceanic & ~continental),
        ("oceanic_mature", ~young_oceanic & ~continental),
    )
    craton_fraction = np.divide(
        craton_hc_m, hc_m, out=np.zeros_like(hc_m, dtype=float), where=hc_m > 0.0
    )
    craton_fraction = np.clip(craton_fraction, 0.0, 1.0)
    strata = {
        "all": np.ones(len(resultant), dtype=bool),
        "craton_none": craton_fraction <= 1e-9,
        "craton_partial": (craton_fraction > 1e-9) & (craton_fraction < 0.5),
        "craton_majority": craton_fraction >= 0.5,
        "inherited_belt": belt_marked,
        "no_inherited_belt": ~belt_marked,
    }
    for name, mask in classes:
        strata[f"crust_{name}"] = mask

    mode_lags = {
        "steady_conductive": np.zeros(len(resultant), dtype=float),
        "moho_lag_interpolated": moho_lag_c,
    }
    summary = {}
    for mode, lags in mode_lags.items():
        capacities_by_node = np.empty(len(resultant), dtype=float)
        crust_labels = np.empty(len(resultant), dtype=object)
        for name, mask in classes:
            crust_labels[mask] = name
        for index in range(len(resultant)):
            # State fields are quantized to 100 m / 1 C solely to reuse column solves;
            # these increments are finer than the physical/rheology uncertainty here.
            capacities_by_node[index] = _column_capacity(
                str(crust_labels[index]),
                round(float(hc_m[index]) / 100.0) * 100.0,
                round(float(hm_m[index]) / 100.0) * 100.0,
                round(float(lags[index])),
            )
        by_stratum = {}
        for name, mask in strata.items():
            selected = in_band & mask & np.isfinite(capacities_by_node) & (capacities_by_node > 0.0)
            if not np.any(selected):
                continue
            caps = capacities_by_node[selected]
            ratios = resultant[selected] / caps
            by_stratum[name] = {
                "nodes": int(selected.sum()),
                "capacity_n_per_m_p10_median_p90": [
                    float(value) for value in np.percentile(caps, (10, 50, 90))
                ],
                "n_ext_over_capacity_median_p95": [
                    float(value) for value in np.percentile(ratios, (50, 95))
                ],
                "over_capacity_nodes": int(np.count_nonzero(ratios >= 1.0)),
                "over_capacity_fraction": float(np.mean(ratios >= 1.0)),
                "hc_median_m": float(np.median(hc_m[selected])),
                "hm_median_m": float(np.median(hm_m[selected])),
                "moho_lag_median_c": float(np.median(moho_lag_c[selected])),
            }
        summary[mode] = {
            "capacity_scenario": {
                "base_temperature_c": COLUMN_SCENARIO.base_temperature_c,
                "strain_rate_s-1": COLUMN_SCENARIO.strain_rate_s,
                "pore_pressure_ratio": "hydrostatic_water_depth_profile",
                "k0": COLUMN_SCENARIO.k0,
                "cohesion_pa": COLUMN_SCENARIO.cohesion_pa,
                "friction_angle_deg": COLUMN_SCENARIO.friction_angle_deg,
                "lower_crust_law": COLUMN_SCENARIO.lower_crust_law,
                "mantle_law": COLUMN_SCENARIO.mantle_law,
                "thermal_mode": mode,
                "lag_profile_assumption": (
                    "subtract lag linearly from zero at surface to the Moho, then linearly back to zero at lithosphere base"
                    if mode == "moho_lag_interpolated" else "steady conductive profile; lag ignored"
                ),
                "column_grid_spacing_m": COLUMN_DZ_M,
                "capacity_cache_rounding": {"Hc_m": 100, "Hm_m": 100, "Moho_lag_c": 1},
            },
            "by_stratum": by_stratum,
        }
    return summary


def _nondegenerate_triangles(points_xyz: np.ndarray, triangles: np.ndarray) -> np.ndarray:
    """Drop projected zero-area elements before connectivity checks or solving."""
    keep = []
    for triangle in triangles:
        element = points_xyz[triangle]
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
            keep.append(triangle)
    return np.asarray(keep, dtype=np.int64).reshape(-1, 3)


def _remap_triangle_nodes(triangles: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Compact triangle indices and return the original node IDs used by the mesh."""
    used = np.unique(triangles)
    remap = np.full(int(triangles.max()) + 1, -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    return used, remap[triangles]


def _triangle_mesh_component_sizes(triangles: np.ndarray) -> np.ndarray:
    """Count full-edge-connected triangle components independent of unused plate nodes."""
    if len(triangles) == 0:
        return np.empty(0, dtype=np.int64)
    used, compact = _remap_triangle_nodes(triangles)
    return thin_shell._mesh_component_sizes(len(used), compact)


def _hull_triangles_inside_outline(plate, points: np.ndarray, weights: np.ndarray) -> tuple[np.ndarray, int]:
    """Clip curved-hull triangles with a stated set of interior outline probes."""
    hull = thin_shell._spherical_surface_triangles(points)
    probes = np.sum(points[hull][:, None, :, :] * weights[None, :, :, None], axis=2)
    probes = geometry.normalize(probes.reshape(-1, 3))
    inside = plate.contains_batch(probes).reshape(len(hull), len(weights)).all(axis=1)
    selected = _nondegenerate_triangles(points, hull[inside])
    return selected, int((~inside).sum())


def _repair_connected_source_mesh(plate, points: np.ndarray, triangles: np.ndarray) -> tuple[np.ndarray, dict | None]:
    """Retry an edge-disconnected clipped mesh with less strict outline probes.

    Exact source-cell adjacency has already established that the plate is one component.
    The strict eight-probe clip can nevertheless leave tiny detached triangulation scraps at
    concave edges. Retriangulate the spherical hull, first requiring the centroid and all three
    edge midpoints inside the plate, then trying centroid-only clipping as a recorded fallback.
    """
    if len(_triangle_mesh_component_sizes(triangles)) <= 1:
        return triangles, None
    midpoint_weights = np.array([
        [1 / 3, 1 / 3, 1 / 3],
        [0.5, 0.5, 0.0],
        [0.5, 0.0, 0.5],
        [0.0, 0.5, 0.5],
    ])
    centroid_weights = midpoint_weights[:1]
    original_components = _triangle_mesh_component_sizes(triangles).tolist()
    for policy, weights in (("centroid_and_edge_midpoints", midpoint_weights), ("centroid_only", centroid_weights)):
        candidate, rejected = _hull_triangles_inside_outline(plate, points, weights)
        if len(_triangle_mesh_component_sizes(candidate)) == 1:
            return candidate, {
                "policy": policy,
                "original_triangle_component_sizes": original_components,
                "repaired_triangle_count": int(len(candidate)),
                "strict_outline_rejected_triangle_count": rejected,
            }

    labels = thin_shell._mesh_component_labels(triangles)
    sizes = np.bincount(labels)
    largest_component = int(np.argmax(sizes))
    retained_fraction = float(sizes[largest_component] / len(triangles))
    if retained_fraction >= MIN_REPAIRED_TRIANGLE_FRACTION:
        repaired = triangles[labels == largest_component]
        return repaired, {
            "policy": "retain_dominant_edge_connected_component",
            "original_triangle_component_sizes": original_components,
            "repaired_triangle_count": int(len(repaired)),
            "discarded_triangle_count": int(len(triangles) - len(repaired)),
            "retained_triangle_fraction": retained_fraction,
        }
    return triangles, None


def _skip_candidate(reason: str, world, plate, points: np.ndarray, kind: str, band: np.ndarray,
                    craton: np.ndarray, **details) -> dict:
    """Keep candidate identity and cut-band coverage metadata on every diagnostic skip."""
    return {
        "skipped": reason,
        "kind": kind,
        "plate": int(plate.plate_id),
        "myr": float(world.elapsed_years / 1e6),
        "node_count": int(len(points)),
        "candidate_band_node_count": int(band.sum()),
        "candidate_band_craton_node_count": int((band & craton).sum()),
        **details,
    }


def _score_event(
    world, plate, points, cut_normal, kind, capacities,
    capture_node_samples: bool = False, include_capacity_scores: bool = True,
) -> dict | None:
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
        candidate_triangles = thin_shell._spherical_surface_triangles(source_points)
        centroids = geometry.normalize(source_points[candidate_triangles].mean(axis=1))
        candidate_triangles = candidate_triangles[plate.contains_batch(centroids)]
        mesh_method = "spherical_hull_centroid_clipped"
    candidate_triangles = _nondegenerate_triangles(source_points, candidate_triangles)
    if len(candidate_triangles) == 0:
        return _skip_candidate("no nondegenerate clipped mesh triangles", world, plate, points, kind, band, craton)

    topology_audit = _topology_audit(plate, len(source_points), candidate_triangles)
    if topology_audit["source_component_count"] > 1:
        return _skip_candidate(
            "source plate adjacency has disconnected cell regions", world, plate, points, kind, band, craton,
            mesh_method=mesh_method,
            mesh_node_count=int(len(np.unique(candidate_triangles))),
            mesh_triangle_count=int(len(candidate_triangles)),
            topology_audit=topology_audit,
        )

    primary_triangle_count = len(candidate_triangles)
    candidate_triangles, mesh_repair = _repair_connected_source_mesh(plate, source_points, candidate_triangles)
    if mesh_repair is not None:
        mesh_method = "spherical_hull_outline_probe_repair"
    if len(_triangle_mesh_component_sizes(candidate_triangles)) > 1:
        return _skip_candidate(
            "thin-shell mesh is disconnected", world, plate, points, kind, band, craton,
            mesh_error=f"thin-shell mesh has disconnected element components with triangle counts {_triangle_mesh_component_sizes(candidate_triangles).tolist()}",
            plate_node_count=int(len(source_points)),
            mesh_node_count=int(len(np.unique(candidate_triangles))),
            mesh_triangle_count=int(len(candidate_triangles)),
            mesh_method=mesh_method,
            topology_audit=topology_audit,
        )

    used, triangles = _remap_triangle_nodes(candidate_triangles)
    in_band = band[used]
    in_band_craton = in_band & craton[used]
    band_coverage = float(in_band.sum() / max(band.sum(), 1))
    craton_band_coverage = float(in_band_craton.sum() / max((band & craton).sum(), 1))
    if band_coverage < MIN_DIAGNOSTIC_BAND_COVERAGE or craton_band_coverage < MIN_DIAGNOSTIC_BAND_COVERAGE:
        return _skip_candidate(
            "mesh does not cover at least 95% of cut-band and craton-band nodes", world, plate, points, kind, band, craton,
            mesh_method=mesh_method, mesh_repair=mesh_repair,
            mesh_band_node_count=int(in_band.sum()), mesh_craton_band_node_count=int(in_band_craton.sum()),
            mesh_band_coverage=band_coverage, mesh_craton_band_coverage=craton_band_coverage,
            topology_audit=topology_audit,
        )
    if not np.any(in_band):
        return _skip_candidate("cut band has no nodes on clipped mesh", world, plate, points, kind, band, craton,
                               mesh_method=mesh_method, topology_audit=topology_audit)
    if not np.any(in_band_craton):
        return _skip_candidate("cut band has no craton nodes on clipped mesh", world, plate, points, kind, band, craton,
                               mesh_method=mesh_method, topology_audit=topology_audit)

    hc = plate.collect("crustal_thickness_m")[used]
    hm = plate.collect("mantle_lithosphere_thickness_m")[used]
    moho_lag = plate.collect("moho_thermal_lag_c")[used]
    craton_hc = plate.collect("craton_crust_m")[used]
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
        return _skip_candidate(
            "thin-shell mesh is disconnected", world, plate, points, kind, band, craton,
            mesh_error=str(error), mesh_node_count=int(len(used)),
            mesh_triangle_count=int(len(triangles)), mesh_method=mesh_method,
            topology_audit=topology_audit,
        )
    resultant = result.max_incremental_extensional_resultant_n_per_m
    scores = {}
    if include_capacity_scores:
        scores = _capacity_scores(resultant, in_band, continental, young_oceanic, belt, capacities)
        scores["spatial_column_capacity"] = _spatial_capacity_scores(
            resultant,
            in_band,
            continental,
            young_oceanic,
            belt,
            hc,
            hm,
            moho_lag,
            craton_hc,
        )
    event = {
        "kind": kind,
        "myr": float(world.elapsed_years / 1e6),
        "plate": int(plate.plate_id),
        "node_count": int(len(source_points)),
        "mesh_node_count": int(len(used)),
        "mesh_method": mesh_method,
        "mesh_repair": mesh_repair,
        "primary_mesh_triangle_count": int(primary_triangle_count),
        "mesh_band_coverage": band_coverage,
        "mesh_craton_band_coverage": craton_band_coverage,
        "mesh_triangle_count": int(len(triangles)),
        "candidate_band_node_count": int(band.sum()),
        "candidate_band_craton_node_count": int((band & craton).sum()),
        "band_node_count": int(in_band.sum()),
        "band_nodes_omitted_by_mesh": int(band.sum() - in_band.sum()),
        "band_craton_node_count": int((in_band & craton[used]).sum()),
        "band_craton_nodes_omitted_by_mesh": int((band & craton).sum() - (in_band & craton[used]).sum()),
        "band_max_resultant_n_per_m": float(np.max(resultant[in_band])),
        "band_p50_resultant_n_per_m": float(np.median(resultant[in_band])),
        "band_p95_resultant_n_per_m": float(np.quantile(resultant[in_band], 0.95)),
        "solver_stop_code": int(result.solver_stop_code),
        "equilibrium_residual_n": float(result.equilibrium_residual_n),
        "equilibrium_residual_over_load": float(
            result.equilibrium_residual_n / max(np.linalg.norm(result.equilibrated_force_n), 1.0)
        ),
        "topology_audit": topology_audit,
        "capacity_scores": scores,
    }
    if capture_node_samples:
        craton_fraction = np.divide(
            craton_hc, hc, out=np.zeros_like(hc, dtype=float), where=hc > 0.0
        )
        crust_labels = np.full(len(hc), "oceanic_mature", dtype=object)
        crust_labels[continental] = "continental"
        crust_labels[young_oceanic & ~continental] = "oceanic_young_ridge"
        sample_indices = np.flatnonzero(in_band)
        event["_node_samples"] = {
            "seed": None,
            "event_index": None,
            "kind": kind,
            "myr": float(world.elapsed_years / 1e6),
            "plate": int(plate.plate_id),
            "source_cell_index": used[sample_indices].tolist(),
            "n_ext_n_per_m": resultant[in_band].astype(float).tolist(),
            "h_crust_m": hc[in_band].astype(float).tolist(),
            "h_mantle_m": hm[in_band].astype(float).tolist(),
            "moho_lag_c": moho_lag[in_band].astype(float).tolist(),
            "crust_type": crust_labels[in_band].tolist(),
            "craton_crust_fraction": np.clip(craton_fraction[in_band], 0.0, 1.0).astype(float).tolist(),
            "inherited_belt": belt[in_band].astype(bool).tolist(),
            "solver_stop_code": int(result.solver_stop_code),
            "equilibrium_residual_over_load": float(
                result.equilibrium_residual_n / max(np.linalg.norm(result.equilibrated_force_n), 1.0)
            ),
            "mesh_band_coverage": band_coverage,
            "mesh_craton_band_coverage": craton_band_coverage,
        }
    return event

def _coverage_report(events: list[dict]) -> dict:
    """Report complete-candidate coverage and eligible connected-source solve coverage."""
    eligible = [
        event for event in events
        if event.get("topology_audit", {}).get("source_component_count") == 1
    ]
    scored = [
        event for event in eligible
        if "skipped" not in event
        and event.get("solver_stop_code") == 1
        and event.get("equilibrium_residual_over_load", 1.0) <= 0.01
        and event.get("mesh_band_coverage", 0.0) >= MIN_DIAGNOSTIC_BAND_COVERAGE
        and event.get("mesh_craton_band_coverage", 0.0) >= MIN_DIAGNOSTIC_BAND_COVERAGE
    ]

    def group_result(group: list[dict]) -> dict:
        group_eligible = [
            event for event in group
            if event.get("topology_audit", {}).get("source_component_count") == 1
        ]
        group_scored = [
            event for event in group_eligible
            if "skipped" not in event
            and event.get("solver_stop_code") == 1
            and event.get("equilibrium_residual_over_load", 1.0) <= 0.01
            and event.get("mesh_band_coverage", 0.0) >= MIN_DIAGNOSTIC_BAND_COVERAGE
            and event.get("mesh_craton_band_coverage", 0.0) >= MIN_DIAGNOSTIC_BAND_COVERAGE
        ]
        return {
            "candidate_events": len(group),
            "source_connected_candidates": len(group_eligible),
            "valid_shell_solves": len(group_scored),
            "coverage_fraction_all_candidates": len(group_scored) / len(group) if group else None,
            "coverage_fraction_source_connected": (
                len(group_scored) / len(group_eligible) if group_eligible else None
            ),
        }

    by_kind = {
        kind: group_result([event for event in events if event.get("kind") == kind])
        for kind in ("split", "failed_rift")
    }
    time_bins = (
        ("0-200_Myr", 0.0, 200.0),
        ("200-500_Myr", 200.0, 500.0),
        ("500-1000_Myr", 500.0, 1000.000001),
    )
    by_world_time = {
        name: group_result([event for event in events if low <= event.get("myr", -1.0) < high])
        for name, low, high in time_bins
    }
    by_kind_and_time = {
        kind: {
            name: group_result([
                event for event in events
                if event.get("kind") == kind and low <= event.get("myr", -1.0) < high
            ])
            for name, low, high in time_bins
        }
        for kind in ("split", "failed_rift")
    }
    minimum_stratum = 0.80
    stratum_fractions = [
        value["coverage_fraction_source_connected"]
        for value in [
            *by_kind.values(),
            *by_world_time.values(),
            *(stratum for groups in by_kind_and_time.values() for stratum in groups.values()),
        ]
        if value["coverage_fraction_source_connected"] is not None
    ]
    overall_fraction_all = len(scored) / len(events) if events else 0.0
    overall_fraction_connected = len(scored) / len(eligible) if eligible else 0.0
    return {
        "candidate_events": len(events),
        "source_connected_candidates": len(eligible),
        "valid_shell_solves": len(scored),
        "coverage_fraction_all_candidates": overall_fraction_all,
        "coverage_fraction_source_connected": overall_fraction_connected,
        "minimum_overall_source_connected_coverage": 0.90,
        "minimum_coverage_per_seed_kind_and_time_stratum_source_connected": minimum_stratum,
        "minimum_band_and_craton_band_node_coverage": MIN_DIAGNOSTIC_BAND_COVERAGE,
        "by_kind": by_kind,
        "by_world_time_myr": by_world_time,
        "by_kind_and_world_time_myr": by_kind_and_time,
        "meets_seed_minimum": overall_fraction_connected >= minimum_stratum,
        "meets_overall_and_stratified_minimums": (
            overall_fraction_connected >= 0.90
            and all(fraction >= minimum_stratum for fraction in stratum_fractions)
        ),
        "disconnected_source_candidates_excluded_from_solver": len(events) - len(eligible),
    }


def _aggregate_spatial_capacity(events: list[dict]) -> dict:
    """Combine node-weighted exceedance counts and event-level median ratios."""
    aggregate: dict[str, dict[str, dict]] = {}
    for event in events:
        for mode, mode_result in event.get("capacity_scores", {}).get("spatial_column_capacity", {}).items():
            mode_summary = aggregate.setdefault(mode, {})
            for stratum, result in mode_result["by_stratum"].items():
                bucket = mode_summary.setdefault(stratum, {
                    "node_count": 0,
                    "over_capacity_nodes": 0,
                    "event_median_ratios": [],
                    "event_capacity_medians_n_per_m": [],
                })
                bucket["node_count"] += result["nodes"]
                bucket["over_capacity_nodes"] += result["over_capacity_nodes"]
                bucket["event_median_ratios"].append(result["n_ext_over_capacity_median_p95"][0])
                bucket["event_capacity_medians_n_per_m"].append(result["capacity_n_per_m_p10_median_p90"][1])
    for mode_summary in aggregate.values():
        for bucket in mode_summary.values():
            bucket["node_weighted_over_capacity_fraction"] = (
                bucket["over_capacity_nodes"] / bucket["node_count"] if bucket["node_count"] else None
            )
            bucket["median_event_median_n_ext_over_capacity"] = float(np.median(bucket.pop("event_median_ratios")))
            bucket["median_event_capacity_median_n_per_m"] = float(np.median(bucket.pop("event_capacity_medians_n_per_m")))
    return {
        "interpretation": "Diagnostic exceedance comparison only; capacities are not fitted to split/failed-rift labels and do not feed back into replay.",
        "by_thermal_mode_and_stratum": aggregate,
    }

def _run(
    seed: int, years: float, out_path: Path,
    node_samples_path: Path | None = None, capture_only: bool = False,
) -> dict:
    for name, value in __import__("compare_cratons")._DEFAULTS.items():
        setattr(cratons, name, value)
    cratons.CRATON_FORMATION_MYR = float("inf")
    capacities = {} if capture_only else _physical_capacities()
    events: list[dict] = []
    if node_samples_path is not None:
        node_samples_path.parent.mkdir(parents=True, exist_ok=True)
        node_samples_path.write_text("", encoding="utf-8")
    original_cut = cratons.record_rift_cut
    original_lowland_relief = getattr(
        initial_terrain, "apply_lowland_spherical_relief", lambda _world: None
    )

    def tapped(world, plate, points, cut_normal, band_sin, kind):
        started = time.perf_counter()
        row = _score_event(
            world, plate, points, cut_normal, kind, capacities,
            capture_node_samples=node_samples_path is not None,
            include_capacity_scores=not capture_only,
        )
        if row is not None:
            row["event_index"] = len(events)
            row["diagnostic_wall_s"] = time.perf_counter() - started
            node_samples = row.pop("_node_samples", None)
            if node_samples is not None and node_samples_path is not None:
                node_samples["seed"] = seed
                node_samples["event_index"] = row["event_index"]
                with node_samples_path.open("a", encoding="utf-8") as samples_file:
                    samples_file.write(json.dumps(node_samples, allow_nan=False) + "\n")
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
        "capture_only": capture_only,
        "gate1_baseline_commit": GATE1_BASELINE_COMMIT,
        "events_scored": len(events),
        "diagnostic_coverage": _coverage_report(events),
        "spatial_capacity_aggregate": _aggregate_spatial_capacity(events),
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
    parser.add_argument(
        "--node-samples-out", type=Path,
        help="Optional JSONL sidecar with exact cut-band N_ext and per-cell state.",
    )
    parser.add_argument(
        "--capture-only", action="store_true",
        help="Skip capacity scoring while exporting node samples for post-processing.",
    )
    args = parser.parse_args()
    if args.capture_only and args.node_samples_out is None:
        parser.error("--capture-only requires --node-samples-out")
    report = _run(args.seed, args.years, args.out, args.node_samples_out, args.capture_only)
    args.out.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in ("seed", "years", "events_scored", "wall_seconds")} ), flush=True)


if __name__ == "__main__":
    main()
