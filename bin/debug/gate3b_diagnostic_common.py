"""Shared helpers for reproducible RFC 0010 Gate 3b offline diagnostics."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from app import elevation_lines, geometry, lithosphere, mantle, torque


def write_json_report(path: Path, report: dict) -> None:
    """Write a deterministic diagnostic report, creating its parent directory."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def edge_connected_continental_components(world):
    """Return edge-connected continental node sets ordered by represented area."""
    components = []
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
        areas = plate.node_areas_m2()[indices]
        for label in range(count):
            selected = indices[labels == label]
            components.append((float(areas[labels == label].sum()), plate, selected))
    components.sort(key=lambda row: row[0], reverse=True)
    return components


def production_force_field(world, plate):
    """Gather the production boundary force field without duplicating call setup."""
    spacing = elevation_lines.line_spacing_rad(world.node_density)
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


def local_failure_components(
    plate,
    cracked,
    direction,
    stress,
    yield_strength,
    residual_strength_fraction,
    opening_threshold=0.5,
):
    """Remove locally failed, opening-supported edges and label remaining plate components."""
    n = plate.node_count()
    adjacency = plate.adjacency()
    rows = np.repeat(np.arange(n, dtype=np.int64), np.diff(adjacency.offsets))
    cols = adjacency.neighbours.astype(np.int64, copy=False)
    pairs = np.unique(np.sort(np.column_stack((rows, cols)), axis=1), axis=0)
    points, _ = plate.all_points_and_elevation()
    edge = points[pairs[:, 1]] - (
        np.sum(points[pairs[:, 1]] * points[pairs[:, 0]], axis=1)[:, None]
        * points[pairs[:, 0]]
    )
    edge = geometry.normalize(edge)
    opens_a = np.abs(np.sum(edge * direction[pairs[:, 0]], axis=1))
    opens_b = np.abs(np.sum(edge * direction[pairs[:, 1]], axis=1))
    supported = (opens_a > opening_threshold) & (opens_b > opening_threshold)
    supported &= stress[pairs[:, 0]] > residual_strength_fraction * yield_strength[pairs[:, 0]]
    supported &= stress[pairs[:, 1]] > residual_strength_fraction * yield_strength[pairs[:, 1]]
    broken = pairs[supported & cracked[pairs[:, 0]] & cracked[pairs[:, 1]]]
    keep = np.ones(len(pairs), dtype=bool)
    if len(broken):
        broken_keys = broken[:, 0] * n + broken[:, 1]
        pair_keys = pairs[:, 0] * n + pairs[:, 1]
        keep = ~np.isin(pair_keys, broken_keys, assume_unique=True)
    survivors = pairs[keep]
    graph = csr_matrix(
        (np.ones(2 * len(survivors)),
         (np.r_[survivors[:, 0], survivors[:, 1]], np.r_[survivors[:, 1], survivors[:, 0]])),
        shape=(n, n),
    )
    count, labels = connected_components(graph, directed=False)
    sizes = np.bincount(labels, minlength=count)
    return broken, count, labels, sizes, pairs


def map_spherical_coordinates(points):
    """Project unit-sphere points to an antimeridian-centered longitude/latitude map."""
    xyz = geometry.normalize(points)
    longitude = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))
    latitude = np.degrees(np.arcsin(np.clip(xyz[:, 2], -1.0, 1.0)))
    center = np.degrees(np.angle(np.mean(np.exp(1j * np.radians(longitude)))))
    longitude = (longitude - center + 180.0) % 360.0 - 180.0 + center
    return xyz, longitude, latitude


def add_native_midline_source(
    world,
    plate,
    target_points,
    source_geometry_fn,
    scaled_centers_fn,
    force_field_fn,
    spacing_rad,
):
    """Install the native-p95 geometric midline source and return its owned force field."""
    native_centers = list(world.mantle_centers)
    native_p95 = float(mantle.rad_per_yr_to_cm_per_yr(np.percentile(
        np.linalg.norm(mantle.flow_at(target_points, native_centers), axis=1), 95.0
    )))
    source_geometry = source_geometry_fn(plate, spacing_rad)
    centers, scale = scaled_centers_fn(
        source_geometry["source_centers_unit_strength"], target_points, native_p95
    )
    world.mantle_centers = native_centers + centers
    field, inputs, subducting = force_field_fn(world, plate)
    source_velocity = (
        mantle.flow_at(geometry.normalize(field.points_xyz), centers)
        * lithosphere.PLANET_RADIUS_M / torque.SECONDS_PER_YEAR
    )
    midline_force = (
        torque.ASTHENOSPHERE_VISCOSITY_PA_S / torque.SHEAR_ZONE_THICKNESS_M
        * field.areas_m2[:, None] * source_velocity
    )
    return native_centers, native_p95, source_geometry, centers, scale, field, inputs, subducting, midline_force
