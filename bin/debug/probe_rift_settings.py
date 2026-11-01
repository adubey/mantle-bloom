#!/usr/bin/env python3
"""Census the four rift settings on saved worlds and route a craton-free rift from each (RFC 0010).

Read-only exploration on saved worlds. Active rifting is plume-heated continental lithosphere,
passive rifting is the production flow-fit gate, slab rollback is an upper plate moving away from
its trench, and transtension is an oblique, opening continental boundary. From each site two arms
grow to the plate boundary by Dijkstra, and the plate must truly separate. Intact craton interiors
are blocked; craton edges, weak cratons and plume-heated interiors with an eroded keel are open at a
cost. No route is a valid outcome: the site still weakens.
The great circle through the same site and axis is the baseline. Thresholds are design choices.

Usage (one process per save; RAM-heavy on mature worlds):
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_rift_settings.py \
      --saves <a.mbworld> <b.mbworld> --out analysis/rfc0010/rift-settings.json \
      --image-dir <private-dir>
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy.cluster.vq import kmeans2
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

from app import cratons, elevation_lines, geometry, lithosphere, mantle, merge_split, persistence, torque

MIN_FRONT_KM = 300.0
RATE_CM_YR = 1.0
OBLIQUITY_DEG = (20.0, 70.0)
BACK_ARC_KM = 300.0
HEAD_AREA_KM2 = 1e5
BELT_RESISTANCE = 0.6
PLUME_RESISTANCE = 0.6
ALIGN_PENALTY = 2.0
CRATON_EDGE_KM = 200.0
CRATON_BLOCK_STRENGTH = 0.5
CRATON_COST = 4.0
KEEL_BREACH_M = 60_000.0
INBOARD_KM = 300.0
MIN_DAUGHTER_SHARE = 0.05
# An arm may not end nearer the site than this, widening until the daughters are real plates.
ARM_REACH_KM = (300.0, 600.0, 1000.0, 1500.0, 2000.0, 3000.0)
CM_PER_RAD_YR = lithosphere.PLANET_RADIUS_M * 100.0


def tangent(vectors, points):
    return vectors - np.einsum("ij,ij->i", vectors, points)[:, None] * points


def components(plate, mask):
    """Edge-connected index sets within `mask`."""
    indices = np.flatnonzero(mask)
    if not len(indices):
        return []
    adjacency = plate.adjacency()
    rows = np.repeat(np.arange(plate.node_count()), np.diff(adjacency.offsets))
    keep = mask[rows] & mask[adjacency.neighbours]
    graph = csr_matrix((np.ones(int(keep.sum())), (rows[keep], adjacency.neighbours[keep])),
                       shape=(plate.node_count(),) * 2)
    _count, labels = connected_components(graph, directed=False)
    return [indices[labels[indices] == label] for label in np.unique(labels[indices])]


def boundary_state(world, plate):
    spacing = elevation_lines.line_spacing_rad(world.node_density)
    reach = torque.boundary_force_reach_rad(spacing)
    cache = getattr(world, "boundary_search_cache", None)
    neighbours = plate.get_neighbours(world.plates)
    inputs = torque.gather_boundary_force_inputs(plate, neighbours, spacing, reach, cache)
    convergent, divergent, transform, _contested = torque.classify_boundary_nodes(
        plate, neighbours, inputs, reach, cache)
    points = geometry.normalize(np.asarray(inputs.own_points, dtype=float))
    direction = tangent(np.nan_to_num(np.asarray(inputs.direction_to_neighbor, dtype=float)), points)
    direction = geometry.normalize(direction)
    relative = np.cross(np.asarray(plate.omega) - inputs.neighbor_omega, points) * CM_PER_RAD_YR
    closing = np.einsum("ij,ij->i", relative, direction)
    shear = np.linalg.norm(relative - closing[:, None] * direction, axis=1)
    return {
        "spacing": spacing, "reach": reach, "inputs": inputs, "points": points,
        "direction": direction, "relative": relative, "closing": closing, "shear": shear,
        "absolute": np.cross(np.asarray(plate.omega), points) * CM_PER_RAD_YR,
        "convergent": convergent, "divergent": divergent, "transform": transform,
        "continental": torque.own_is_continental(plate, inputs),
        "subducting": torque.subducting_boundary_mask(plate, inputs, reach),
    }


def _front(state, indices):
    """Site on a boundary front, moved inboard, with its trench-normal direction."""
    points, direction = state["points"], state["direction"]
    middle = indices[np.argmin(np.linalg.norm(points[indices] - points[indices].mean(axis=0), axis=1))]
    outward = geometry.normalize(direction[indices].mean(axis=0, keepdims=True))[0]
    return int(middle), outward


def _inboard(state, point, outward):
    """The plate's node nearest a point moved inboard from a boundary front."""
    target = geometry.normalize(point - outward * INBOARD_KM * 1e3 / lithosphere.PLANET_RADIUS_M)
    return int(np.argmax(state["points"] @ target))


def _sites_rollback(plate, state):
    band = state["convergent"] & ~state["subducting"] & state["inputs"].neighbor_is_oceanic
    away = np.einsum("ij,ij->i", state["absolute"], state["direction"]) < -RATE_CM_YR
    sites = []
    length_km = state["spacing"] * lithosphere.PLANET_RADIUS_M / 1e3
    for indices in components(plate, band & away):
        if len(indices) * length_km < MIN_FRONT_KM:
            continue
        middle, outward = _front(state, indices)
        inland = geometry.normalize(state["points"][middle] - outward * BACK_ARC_KM * 1e3 / lithosphere.PLANET_RADIUS_M)
        seed = int(np.argmax(state["points"] @ inland))
        sites.append({"setting": "rollback", "seed": seed, "extension": outward, "front_nodes": len(indices)})
    return sites


def _sites_transtension(plate, state):
    angle = np.degrees(np.arctan2(state["shear"], np.maximum(-state["closing"], 1e-9)))
    oblique = (-state["closing"] > RATE_CM_YR) & (angle >= OBLIQUITY_DEG[0]) & (angle <= OBLIQUITY_DEG[1])
    band = (state["divergent"] | state["transform"]) & oblique & state["continental"]
    sites = []
    length_km = state["spacing"] * lithosphere.PLANET_RADIUS_M / 1e3
    for indices in components(plate, band):
        if len(indices) * length_km < MIN_FRONT_KM:
            continue
        middle, outward = _front(state, indices)
        middle = _inboard(state, state["points"][middle], outward)
        extension = geometry.normalize(state["relative"][indices].mean(axis=0, keepdims=True))[0]
        sites.append({"setting": "transtension", "seed": middle, "extension": extension,
                      "front_nodes": len(indices)})
    return sites


def _sites_active(plate, state):
    heated = (plate.collect("plume_heat_myr") > 0.0) & state["continental"]
    areas = plate.node_areas_m2()
    sites = []
    for indices in components(plate, heated):
        if areas[indices].sum() / 1e6 < HEAD_AREA_KM2:
            continue
        points = state["points"]
        centre = geometry.normalize(np.average(points[indices], axis=0, weights=areas[indices]))
        seed = int(np.argmax(points @ centre))
        # A plume head has no preferred axis; take the plate's long axis as the line of least work.
        centred = points - points.mean(axis=0)
        _u, _s, vt = np.linalg.svd(centred, full_matrices=False)
        sites.append({"setting": "active", "seed": seed, "extension": vt[0], "front_nodes": len(indices)})
    return sites


def _sites_passive(world, plate):
    points, _ = plate.all_points_and_elevation()
    if plate.node_count() < max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density)):
        return []
    velocities = mantle.flow_at(points, world.mantle_centers)
    _, radius = geometry.bounding_sphere(points)
    relax = max(min(1.0, radius / merge_split.SPLIT_SIZE_CERTAIN_RIFT_RAD),
                min(1.0, plate.internal_stress / merge_split.STRESS_CERTAIN_RIFT))
    if merge_split._fit_residual_rms(points, velocities, plate.omega) < merge_split.SPLIT_RMS_RESIDUAL_THRESHOLD * (1.0 - relax):
        return []
    _, labels = kmeans2(velocities, k=2, minit="++", seed=plate.plate_id)
    if len(np.unique(labels)) < 2:
        return []
    refined = merge_split._refine_split_cut(points, labels == 0)
    if refined is None:
        return []
    mask_a, cut_normal = refined
    pole_a = mantle.fit_euler_pole(points[mask_a], velocities[mask_a])
    pole_b = mantle.fit_euler_pole(points[~mask_a], velocities[~mask_a])
    if np.linalg.norm(pole_a - pole_b) < merge_split.SPLIT_MIN_POLE_SEPARATION * (1.0 - relax):
        return []
    points = geometry.normalize(points)
    on_cut = np.flatnonzero(np.abs(points @ cut_normal) <= np.percentile(np.abs(points @ cut_normal), 5))
    centre = geometry.normalize(points.mean(axis=0, keepdims=True))[0]
    seed = int(on_cut[np.argmax(points[on_cut] @ centre)])
    return [{"setting": "passive", "seed": seed, "extension": cut_normal, "front_nodes": 0}]


def _graph(plate, points, extension, forbidden, resistance):
    adjacency = plate.adjacency()
    rows = np.repeat(np.arange(plate.node_count()), np.diff(adjacency.offsets))
    cols = adjacency.neighbours
    keep = ~forbidden[rows] & ~forbidden[cols]
    rows, cols = rows[keep], cols[keep]
    chord = points[cols] - points[rows]
    length = np.maximum(np.linalg.norm(chord, axis=1), 1e-12)
    along = chord / length[:, None]
    tangent_axis = extension - np.einsum("ij,j->i", 0.5 * (points[rows] + points[cols]), extension)[:, None] * (
        0.5 * (points[rows] + points[cols]))
    norm = np.maximum(np.linalg.norm(tangent_axis, axis=1), 1e-12)
    opening = np.abs(np.einsum("ij,ij->i", along, tangent_axis / norm[:, None]))
    weight = length * 0.5 * (resistance[rows] + resistance[cols]) * (1.0 + ALIGN_PENALTY * opening)
    return csr_matrix((weight, (rows, cols)), shape=(len(points),) * 2)


def _walk(predecessors, target):
    route = [int(target)]
    while predecessors[route[-1]] >= 0:
        route.append(int(predecessors[route[-1]]))
    return route[::-1]


def route_rift(plate, points, seed, extension, forbidden, resistance, reach_rad):
    """Two node-disjoint arms from `seed` to the plate boundary, either side of the strike."""
    degree = np.diff(plate.adjacency().offsets)
    edge = degree < np.median(degree)
    if forbidden[seed]:
        open_cells = np.flatnonzero(~forbidden)
        if not len(open_cells):
            return None, "no_open_cells"
        seed = int(open_cells[np.argmax(points[open_cells] @ points[seed])])
    strike = geometry.normalize(np.cross(points[seed], extension))
    side = (points - points[seed]) @ strike
    edge = edge & (np.arccos(np.clip(points @ points[seed], -1.0, 1.0)) >= reach_rad)
    distances, previous = dijkstra(_graph(plate, points, extension, forbidden, resistance), indices=seed,
                                   return_predecessors=True)
    targets = np.flatnonzero(edge & np.isfinite(distances) & (side > 0.0))
    if not len(targets):
        return None, "no_first_arm"
    first = _walk(previous, targets[np.argmin(distances[targets])])
    blocked = forbidden.copy()
    blocked[first[1:]] = True
    distances, previous = dijkstra(_graph(plate, points, extension, blocked, resistance), indices=seed,
                                   return_predecessors=True)
    targets = np.flatnonzero(edge & np.isfinite(distances) & (side < 0.0))
    if not len(targets):
        return None, "no_second_arm"
    second = _walk(previous, targets[np.argmin(distances[targets])])
    return second[::-1] + first[1:], "routed"


def split_by_route(plate, points, route, areas):
    """Daughter mask from removing route cells, or None if the plate does not separate."""
    keep = np.ones(plate.node_count(), dtype=bool)
    keep[route] = False
    pieces = sorted(components(plate, keep), key=lambda piece: -areas[piece].sum())
    if len(pieces) < 2:
        return None
    label = np.full(plate.node_count(), -1)
    label[pieces[0]] = 0
    label[pieces[1]] = 1
    assigned = label >= 0
    _d, nearest = cKDTree(points[assigned]).query(points[~assigned])
    label[~assigned] = label[assigned][nearest]
    return label == 0


def craton_barrier(plate, cratons_mask, spacing):
    """Intact interiors stop a rift. Edges and weak cratons don't, and neither does a heated
    interior whose keel has gone: plume heat alone never opens a craton."""
    adjacency = plate.adjacency()
    rows = np.repeat(np.arange(plate.node_count()), np.diff(adjacency.offsets))
    interior = cratons_mask.copy()
    for _ in range(max(1, round(CRATON_EDGE_KM * 1e3 / (spacing * lithosphere.PLANET_RADIUS_M)))):
        exposed = np.zeros(plate.node_count(), dtype=bool)
        exposed[rows[~interior[adjacency.neighbours]]] = True
        interior &= ~exposed
    strong = cratons.strength(plate.collect("craton_crust_m")) >= CRATON_BLOCK_STRENGTH
    breached = (plate.collect("plume_heat_myr") > 0.0) & (
        plate.collect("mantle_lithosphere_thickness_m") < KEEL_BREACH_M)
    return interior & strong & ~breached, interior


def route_resistance(plate, belts, forbidden):
    resistance = np.ones(plate.node_count())
    resistance[belts] = BELT_RESISTANCE
    resistance[plate.collect("plume_heat_myr") > 0.0] = PLUME_RESISTANCE
    craton_strength = cratons.strength(plate.collect("craton_crust_m"))
    return np.where(craton_strength > 0.0, 1.0 + CRATON_COST * craton_strength, resistance)


def score_site(world, plate, state, site, cratons_mask, belts):
    points, areas = state["points"], plate.node_areas_m2()
    spacing = state["spacing"]
    forbidden, interior = craton_barrier(plate, cratons_mask, spacing)
    resistance = route_resistance(plate, belts, forbidden)
    extension = geometry.normalize(np.asarray(site["extension"], dtype=float))
    seed = site["seed"]
    band = np.sin(merge_split.FAILED_RIFT_BAND_MULT * spacing)
    # The baseline great circle passes through the site, normal to the extension.
    normal = geometry.normalize(extension - (extension @ points[seed]) * points[seed])
    line_band = np.abs(points @ normal) < band
    row = {"setting": site["setting"], "plate_id": int(plate.plate_id),
           "plate_area_share": float(areas.sum() / (4.0 * np.pi * lithosphere.PLANET_RADIUS_M**2)),
           "front_km": site["front_nodes"] * spacing * lithosphere.PLANET_RADIUS_M / 1e3,
           "line_craton_cells": int((cratons_mask & line_band).sum()),
           "line_craton_area_km2": float(areas[cratons_mask & line_band].sum() / 1e6),
           "seed_blocked": bool(forbidden[seed])}
    route, mask, status = None, None, "no_first_arm"
    for reach_km in ARM_REACH_KM:
        candidate, status = route_rift(plate, points, seed, extension, forbidden, resistance,
                                       reach_km * 1e3 / lithosphere.PLANET_RADIUS_M)
        if candidate is None:
            continue
        pieces = split_by_route(plate, points, candidate, areas)
        if pieces is None:
            status = "does_not_separate"
            continue
        if min(areas[pieces].sum(), areas[~pieces].sum()) >= MIN_DAUGHTER_SHARE * areas.sum():
            route, mask, status = candidate, pieces, "separates"
            row["arm_reach_km"] = reach_km
            break
        status = "daughter_too_small"
    row["status"] = status
    if route is None:
        return row, None
    segments = np.linalg.norm(np.diff(points[route], axis=0), axis=1).sum()
    ends = float(np.arccos(np.clip(points[route[0]] @ points[route[-1]], -1.0, 1.0)))
    smaller = float(min(areas[mask].sum(), areas[~mask].sum()) / areas.sum())
    craton_area = float(areas[cratons_mask].sum())
    pieces = components(plate, cratons_mask)
    divided = [piece for piece in pieces if mask[piece].any() and (~mask[piece]).any()]
    cut = len(divided)
    # The craton a cut detaches: an edge sliver and a craton halved both count as divided.
    detached = [min(areas[piece][mask[piece]].sum(), areas[piece][~mask[piece]].sum()) for piece in divided]
    share = max((part / areas[piece].sum() for part, piece in zip(detached, divided)), default=0.0)
    min_nodes = max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density))
    row.update(route_km=float(segments * lithosphere.PLANET_RADIUS_M / 1e3),
               tortuosity=float(segments / max(ends, 1e-9)), smaller_daughter_fraction=smaller,
               smaller_daughter_nodes=int(min(mask.sum(), (~mask).sum())),
               meets_min_nodes=bool(min(mask.sum(), (~mask).sum()) >= min_nodes),
               route_craton_cells=int(cratons_mask[route].sum()), cratons_divided=cut,
               craton_detached_km2=float(sum(detached) / 1e6), craton_largest_detached_share=float(share),
               route_blocked_cells=int(forbidden[route].sum()), route_cells=len(route),
               route_interior_cells=int(interior[route].sum()),
               craton_area_km2=craton_area / 1e6, belt_share=float(belts[route].mean()))
    return row, {"route": route, "mask": mask, "line_band": line_band, "cratons": cratons_mask}


def draw(path, plate, state, row, snapshot):
    points = state["points"]
    centre = geometry.normalize(points.mean(axis=0, keepdims=True))[0]
    east = geometry.normalize(np.cross([0.0, 0.0, 1.0], centre))
    north = np.cross(centre, east)
    x, y = points @ east, points @ north
    scale = 720.0 / max(np.ptp(x), np.ptp(y))
    px = 40 + (x - x.min()) * scale
    py = 760 - (y - y.min()) * scale
    image = Image.new("RGB", (800, 800), "white")
    canvas = ImageDraw.Draw(image)
    for i in range(len(points)):
        fill = "#b9cfe4" if snapshot["mask"][i] else "#c9e2c4"
        canvas.rectangle((px[i] - 3, py[i] - 3, px[i] + 3, py[i] + 3), fill=fill)
    for i in np.flatnonzero(snapshot["line_band"]):
        canvas.point((px[i], py[i]), fill="#777777")
    for i in np.flatnonzero(snapshot["cratons"]):
        canvas.rectangle((px[i] - 3, py[i] - 3, px[i] + 3, py[i] + 3), fill="#d9a321")
    route = snapshot["route"]
    canvas.line([(px[i], py[i]) for i in route], fill="#c0202c", width=4)
    canvas.text((10, 8), f"{row['setting']} plate {row['plate_id']} tortuosity {row['tortuosity']:.2f}", fill="black")
    image.save(path)


def survey(path: Path, image_dir: Path | None):
    with path.open("rb") as handle:
        world = persistence.load_world_file(handle)
    rows = []
    best = {}
    for plate in world.plates:
        state = boundary_state(world, plate)
        cratons_mask = plate.collect("craton_crust_m") > 0.0
        belts = plate.collect("belt_formed_years") > 0.0
        sites = (_sites_active(plate, state) + _sites_passive(world, plate)
                 + _sites_rollback(plate, state) + _sites_transtension(plate, state))
        for site in sites:
            row, snapshot = score_site(world, plate, state, site, cratons_mask, belts)
            rows.append(row)
            if snapshot is not None and row["line_craton_cells"] and row["meets_min_nodes"]:
                if row["line_craton_area_km2"] > best.get(row["setting"], (0.0,))[0]:
                    best[row["setting"]] = (row["line_craton_area_km2"], plate, state, row, snapshot)
    if image_dir is not None:
        image_dir.mkdir(parents=True, exist_ok=True)
        for setting, (_area, plate, state, row, snapshot) in best.items():
            draw(image_dir / f"{path.stem}-{setting}.png", plate, state, row, snapshot)
    return {"save": path.name, "seed": int(world.seed), "elapsed_myr": float(world.elapsed_years / 1e6),
            "plates": len(world.plates), "sites": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--saves", type=Path, nargs="+", required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--image-dir", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    report = {"knobs": {"min_front_km": MIN_FRONT_KM, "rate_cm_yr": RATE_CM_YR, "obliquity_deg": OBLIQUITY_DEG,
                        "back_arc_km": BACK_ARC_KM, "head_area_km2": HEAD_AREA_KM2,
                        "belt_resistance": BELT_RESISTANCE, "plume_resistance": PLUME_RESISTANCE,
                        "align_penalty": ALIGN_PENALTY, "craton_edge_km": CRATON_EDGE_KM,
                        "craton_block_strength": CRATON_BLOCK_STRENGTH, "craton_cost": CRATON_COST,
                        "keel_breach_m": KEEL_BREACH_M},
              "worlds": [survey(path, args.image_dir) for path in args.saves]}
    report["runtime_seconds"] = time.perf_counter() - started
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True, default=float) + "\n")


if __name__ == "__main__":
    main()
