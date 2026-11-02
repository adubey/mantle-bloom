#!/usr/bin/env python3
"""Census the four rift settings and route a craton-free rift from each (RFC 0010, #737).

Two read-only modes. `--saves` surveys snapshots. `--seed` is Gate D1: it runs one fresh world
with production topology unchanged and no detected rift acting on it, detects at the point in the
step where the lifecycle will sit, and tracks each merged front until it persists for ONSET_MYR.
`--summarize` turns the per-seed reports into the gate's tables.

Active rifting is plume-heated continental lithosphere,
passive rifting is the production flow-fit gate, slab rollback is an upper plate moving away from
its trench, and transtension is an oblique, opening continental boundary. From each site two arms
grow to the plate boundary by Dijkstra, and the plate must truly separate. Intact craton interiors
are blocked; craton edges, weak cratons and plume-heated interiors with an eroded keel are open at a
cost. No route is a valid outcome: the site still weakens.
The great circle through the same site and axis is the baseline. Thresholds are design choices.

Usage (one process per save or seed; RAM-heavy on mature worlds):
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_rift_settings.py \
      --saves <a.mbworld> <b.mbworld> --out analysis/rfc0010/rift-settings.json \
      --image-dir <private-dir>
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_rift_settings.py \
      --seed <n> --myr 200 --step-myr 1 --out <seed.json> --save-dir <private-dir>
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_rift_settings.py \
      --summarize <seed-a.json> <seed-b.json> <seed-c.json> --out analysis/rfc0010/gate-d1.json
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import resource
import statistics
import time
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw
from scipy.cluster.vq import kmeans2
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

from app import (
    collision_polarity,
    cratons,
    elevation_lines,
    geometry,
    lithosphere,
    mantle,
    merge_split,
    persistence,
    plumes,
    torque,
)
from app import world as world_mod

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
# New divergent boundaries thin and stretch what they touch, so routes prefer to stand off a craton.
CRATON_MARGIN_KM = 150.0
CRATON_MARGIN_COST = 1.5
KEEL_BREACH_M = 60_000.0
INBOARD_KM = 300.0
MIN_DAUGHTER_SHARE = 0.05
# An arm may not end nearer the site than this, widening until the daughters are real plates.
ARM_REACH_KM = (300.0, 600.0, 1000.0, 1500.0, 2000.0, 3000.0)
CM_PER_RAD_YR = lithosphere.PLANET_RADIUS_M * 100.0
ONSET_MYR = 2.0
LAPSE_MYR = 5.0
PROPAGATION_KM_PER_MYR = 100.0
RETRY_MYR = 50.0
COOLDOWN_MYR = merge_split.SPLIT_MIN_AGE_YEARS / 1e6
# The RFC leaves "near an old site" open; this is the census's own design choice.
RETRY_NEAR_KM = 500.0
# A rift whose route could not be planned is assumed to hold its plate for the RFC's reference
# 3,000 km arm at the design speed.
UNROUTED_GROWTH_MYR = 15.0
PLATE_CEILING = 40
# Older ocean is the lower plate only past this age gap (#737): at small offsets the literature
# doesn't support age deciding polarity, so those fronts stay undecided. Design choice, frozen
# before the amendment ran and not tuned on the saves that motivated it.
AGE_POLARITY_MIN_GAP_MYR = 10.0
# The ocean-ocean polarity read counts as working above this decided share and below this flip
# share, both by qualifying front (#737); design choices frozen before the gate ran.
POLARITY_MIN_DECIDED = 0.8
POLARITY_MAX_FLIPPED = 0.1
TOPOLOGY_PHASE = "Merging and splitting plates"
SPLIT_MESSAGE = re.compile(r"Plate (\d+) split into plates")
MERGE_MESSAGE = re.compile(r"Plates (\d+) and (\d+) (?:collided and merged|had overlapped so long they fused)")
FAILED_MESSAGE = "began rifting but the rift failed"


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
    # No shared cache: the census runs mid-step and must not touch the step's own search state.
    cache = None
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


def merged_fronts(points, neighbour_ids, indices, link_rad):
    """Merge boundary fragments: nodes facing one neighbour and within `link_rad` of each other
    are one front, so a stretch that drops a node or two isn't counted as several (#737)."""
    fronts = []
    for neighbour in np.unique(neighbour_ids[indices]):
        facing = indices[neighbour_ids[indices] == neighbour]
        pairs = cKDTree(points[facing]).sparse_distance_matrix(
            cKDTree(points[facing]), link_rad, output_type="coo_matrix")
        _count, labels = connected_components(pairs, directed=False)
        fronts.extend((int(neighbour), facing[labels == label]) for label in np.unique(labels))
    return fronts


def front_km(state, indices):
    """Along-boundary length: distinct neighbour nodes faced, each one spacing long. Unlike a node
    count it doesn't grow with the width of the boundary band (#468)."""
    faced = np.unique(state["inputs"].neighbor_node_index[indices])
    return len(faced) * state["spacing"] * lithosphere.PLANET_RADIUS_M / 1e3


def _band_fronts(state, band):
    link = collision_polarity.FRONT_LINK_SPACINGS * state["spacing"]
    indices = np.flatnonzero(band & (state["inputs"].neighbor_plate_id >= 0))
    if not len(indices):
        return []
    fronts = merged_fronts(state["points"], state["inputs"].neighbor_plate_id, indices, link)
    return [(neighbour, nodes) for neighbour, nodes in fronts if front_km(state, nodes) >= MIN_FRONT_KM]


def front_polarity(world, plate, state, neighbour, nodes, plates_by_id):
    """"upper" or "lower" for `plate` at a front, from the consumption evidence alone, or
    "ambiguous" or "none". Slab evidence can't decide here: every oceanic node counts as
    down-going, so an ocean-ocean front has slab votes on both sides that cancel."""
    other = plates_by_id.get(neighbour)
    if other is None:
        return "none"
    other_nodes = state["inputs"].neighbor_node_index[nodes]
    fronts = {plate.plate_id: state["points"][nodes],
              neighbour: other.all_points_and_elevation()[0][other_nodes]}
    frames = {plate.plate_id: plate.frame, neighbour: other.frame}
    radius = collision_polarity.EVIDENCE_LOOKUP_SPACINGS * state["spacing"]

    def select(side, store):
        rows = (store.source == collision_polarity.SOURCE_CONSUMPTION) & (
            store.neighbour == (neighbour if side == plate.plate_id else plate.plate_id))
        near = np.zeros(len(store), dtype=bool)
        if np.any(rows):
            stored = geometry.to_world(frames[side], store.local_points()[rows])
            near[rows] = np.isfinite(cKDTree(fronts[side]).query(stored, distance_upper_bound=radius)[0])
        return near

    toward, against = collision_polarity._tally(world, (plate.plate_id, neighbour), select)
    source = collision_polarity.SOURCE_CONSUMPTION
    verdict = collision_polarity._verdict(float(toward[source]), float(against[source]))
    return {"a": "upper", "b": "lower", None: "none"}.get(verdict, verdict)


def age_verdict(own_created_years, other_created_years, now_myr):
    """"upper" or "lower" for the plate owning `own_created_years`, older ocean going under, or
    "tied" inside AGE_POLARITY_MIN_GAP_MYR. Cells from t = 0 and untracked cells read as the
    run's full age: all starting ocean is one age, so two such sides tie instead of guessing."""
    own = now_myr - max(float(np.median(own_created_years)) / 1e6, 0.0)
    other = now_myr - max(float(np.median(other_created_years)) / 1e6, 0.0)
    if abs(own - other) < AGE_POLARITY_MIN_GAP_MYR:
        return "tied"
    return "lower" if own > other else "upper"


def front_age_verdict(world, plate, state, neighbour, nodes, plates_by_id):
    other = plates_by_id.get(neighbour)
    if other is None:
        return "tied"
    other_nodes = state["inputs"].neighbor_node_index[nodes]
    return age_verdict(plate.collect("node_created_years")[nodes], other.collect("node_created_years")[other_nodes],
                       world.elapsed_years / 1e6)


def _sites_rollback(world, plate, state, plates_by_id, polarity_lookup=None):
    continental = state["continental"]
    away = np.einsum("ij,ij->i", state["absolute"], state["direction"]) < -RATE_CM_YR
    band = state["convergent"] & state["inputs"].neighbor_is_oceanic & away & ~(continental & state["subducting"])
    sites = []
    for neighbour, nodes in _band_fronts(state, band):
        upper_type = "continental" if continental[nodes].mean() >= 0.5 else "oceanic"
        polarity = (polarity_lookup(plate, state, neighbour, nodes) if polarity_lookup
                    else front_age_verdict(world, plate, state, neighbour, nodes, plates_by_id))
        # Against ocean, continent is the upper plate by buoyancy; only ocean-ocean needs evidence.
        if upper_type == "oceanic" and polarity != "upper":
            continue
        middle, outward = _front(state, nodes)
        sites.append({"setting": "rollback", "seed": _inboard(state, state["points"][middle], outward),
                      "extension": outward, "front_km": front_km(state, nodes), "neighbour": neighbour,
                      "nodes": nodes, "upper_type": upper_type, "polarity": polarity})
    return sites


def _sites_transtension(plate, state):
    angle = np.degrees(np.arctan2(state["shear"], np.maximum(-state["closing"], 1e-9)))
    oblique = (-state["closing"] > RATE_CM_YR) & (angle >= OBLIQUITY_DEG[0]) & (angle <= OBLIQUITY_DEG[1])
    band = (state["divergent"] | state["transform"]) & oblique & state["continental"]
    sites = []
    for neighbour, nodes in _band_fronts(state, band):
        middle, outward = _front(state, nodes)
        extension = geometry.normalize(state["relative"][nodes].mean(axis=0, keepdims=True))[0]
        sites.append({"setting": "transtension", "seed": _inboard(state, state["points"][middle], outward),
                      "extension": extension, "front_km": front_km(state, nodes), "neighbour": neighbour,
                      "nodes": nodes})
    return sites


def polarity_fronts(world, plate, state, plates_by_id):
    """Every qualifying convergent front against ocean, whatever the upper plate is doing, for
    Gate D1's check of polarity. Both sides are returned, since rollback looks fronts up from either
    plate; `counted` keeps each ocean-ocean front in the statistics once per pair."""
    rows = []
    band = state["convergent"] & state["inputs"].neighbor_is_oceanic
    for neighbour, nodes in _band_fronts(state, band):
        kind = "continent_ocean" if state["continental"][nodes].mean() >= 0.5 else "ocean_ocean"
        rows.append({"neighbour": neighbour, "nodes": nodes, "kind": kind,
                     "counted": kind == "continent_ocean" or plate.plate_id < neighbour,
                     "age": front_age_verdict(world, plate, state, neighbour, nodes, plates_by_id),
                     "verdict": front_polarity(world, plate, state, neighbour, nodes, plates_by_id)})
    return rows


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
        sites.append({"setting": "active", "seed": seed, "extension": vt[0],
                      "area_km2": float(areas[indices].sum() / 1e6), "nodes": indices})
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
    return [{"setting": "passive", "seed": seed, "extension": cut_normal, "relax": float(relax)}]


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


def longest_arm_km(points, route, seed):
    """Length of the longer arm, either side of `seed`, along a route through it."""
    steps = np.linalg.norm(np.diff(points[route], axis=0), axis=1)
    at = route.index(seed) if seed in route else len(route) // 2
    return float(max(steps[:at].sum(), steps[at:].sum()) * lithosphere.PLANET_RADIUS_M / 1e3)


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


def craton_margin(plate, cratons_mask, spacing):
    """Non-craton cells within CRATON_MARGIN_KM of a craton."""
    adjacency = plate.adjacency()
    rows = np.repeat(np.arange(plate.node_count()), np.diff(adjacency.offsets))
    near = cratons_mask.copy()
    for _ in range(max(1, round(CRATON_MARGIN_KM * 1e3 / (spacing * lithosphere.PLANET_RADIUS_M)))):
        grown = near.copy()
        grown[adjacency.neighbours[near[rows]]] = True
        near = grown
    return near & ~cratons_mask


def route_resistance(plate, belts, forbidden, margin):
    resistance = np.ones(plate.node_count())
    resistance[belts] = BELT_RESISTANCE
    resistance[plate.collect("plume_heat_myr") > 0.0] = PLUME_RESISTANCE
    resistance[margin] *= CRATON_MARGIN_COST
    craton_strength = cratons.strength(plate.collect("craton_crust_m"))
    # A weak craton must not undercut its own margin, or routes step into it to avoid the margin.
    return np.where(craton_strength > 0.0,
                    np.maximum(1.0 + CRATON_COST * craton_strength, CRATON_MARGIN_COST), resistance)


def score_site(world, plate, state, site, cratons_mask, belts):
    points, areas = state["points"], plate.node_areas_m2()
    spacing = state["spacing"]
    forbidden, interior = craton_barrier(plate, cratons_mask, spacing)
    margin = craton_margin(plate, cratons_mask, spacing)
    resistance = route_resistance(plate, belts, forbidden, margin)
    extension = geometry.normalize(np.asarray(site["extension"], dtype=float))
    seed = site["seed"]
    band = np.sin(merge_split.FAILED_RIFT_BAND_MULT * spacing)
    # The baseline great circle passes through the site, normal to the extension.
    normal = geometry.normalize(extension - (extension @ points[seed]) * points[seed])
    line_band = np.abs(points @ normal) < band
    row = {"setting": site["setting"], "plate_id": int(plate.plate_id),
           "plate_area_share": float(areas.sum() / (4.0 * np.pi * lithosphere.PLANET_RADIUS_M**2)),
           "front_km": site.get("front_km", 0.0),
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
               longest_arm_km=longest_arm_km(points, route, seed),
               tortuosity=float(segments / max(ends, 1e-9)), smaller_daughter_fraction=smaller,
               smaller_daughter_nodes=int(min(mask.sum(), (~mask).sum())),
               meets_min_nodes=bool(min(mask.sum(), (~mask).sum()) >= min_nodes),
               route_craton_cells=int(cratons_mask[route].sum()), cratons_divided=cut,
               craton_detached_km2=float(sum(detached) / 1e6), craton_largest_detached_share=float(share),
               route_blocked_cells=int(forbidden[route].sum()), route_margin_cells=int(margin[route].sum()), route_cells=len(route),
               route_interior_cells=int(interior[route].sum()),
               craton_area_km2=craton_area / 1e6, belt_share=float(belts[route].mean()))
    return row, {"route": route, "mask": mask, "line_band": line_band, "cratons": cratons_mask}


def detect_sites(world, plate, state, plates_by_id, polarity_lookup=None):
    return (_sites_active(plate, state) + _sites_passive(world, plate)
            + _sites_rollback(world, plate, state, plates_by_id, polarity_lookup) + _sites_transtension(plate, state))


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
    plates_by_id = {plate.plate_id: plate for plate in world.plates}
    for plate in world.plates:
        state = boundary_state(world, plate)
        cratons_mask = plate.collect("craton_crust_m") > 0.0
        belts = plate.collect("belt_formed_years") > 0.0
        for site in detect_sites(world, plate, state, plates_by_id):
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


@dataclass
class Track:
    """One detection followed across steps: a front in its plate's own frame, or for active and
    passive the plate and setting alone, since their sites move between steps."""

    track_id: int
    plate_id: int
    setting: str
    neighbour: int
    first_myr: float
    last_myr: float
    local: np.ndarray | None
    onset_myr: float | None = None
    row: dict | None = None
    verdict: str | None = None
    axes: list = field(default_factory=list)


class Tracker:
    """Persistence for detections. Before onset a missed step resets the clock; after onset a
    track bridges gaps of up to `lapse_myr`, as a running rift would."""

    def __init__(self, onset_myr, lapse_myr, match_rad):
        self.onset_myr, self.lapse_myr, self.match_rad = onset_myr, lapse_myr, match_rad
        self.open: list[Track] = []
        self.closed: list[tuple[Track, str]] = []
        self.next_id = 0

    def end_plate(self, plate_id, cause):
        """A rift ends if its parent plate splits or merges, even where the plate keeps its id."""
        self.closed.extend((track, cause) for track in self.open if track.plate_id == plate_id)
        self.open = [track for track in self.open if track.plate_id != plate_id]

    def _claim(self, plate_id, setting, neighbour, local, claimed):
        candidates = [track for track in self.open
                      if (track.plate_id, track.setting, track.neighbour) == (plate_id, setting, neighbour)
                      and track.track_id not in claimed]
        if local is not None:
            candidates = [track for track in candidates if track.local is None or np.any(
                cKDTree(track.local).query(local, distance_upper_bound=self.match_rad)[0] < np.inf)]
        return min(candidates, key=lambda track: track.first_myr, default=None)

    def update(self, now, detections, live_plates):
        """Advance one step. `detections` are (plate, setting, neighbour, local points) tuples.
        Returns (index into detections, track) for each detection now past onset for the first
        time, and every open track each detection landed on, by detection index."""
        claimed, landed = set(), {}
        order = sorted(range(len(detections)), key=lambda i: -(0 if detections[i][3] is None else len(detections[i][3])))
        for index in order:
            plate_id, setting, neighbour, local = detections[index]
            track = self._claim(plate_id, setting, neighbour, local, claimed)
            if track is None:
                track = Track(self.next_id, plate_id, setting, neighbour, now, now, local)
                self.next_id += 1
                self.open.append(track)
            track.last_myr, track.local = now, local
            claimed.add(track.track_id)
            landed[index] = track
        keep = []
        for track in self.open:
            if track.plate_id not in live_plates:
                self.closed.append((track, "plate_gone"))
            elif track.track_id in claimed:
                keep.append(track)
            elif track.onset_myr is None:
                self.closed.append((track, "lapsed"))
            elif now - track.last_myr > self.lapse_myr:
                self.closed.append((track, "lapsed"))
            else:
                keep.append(track)
        self.open = keep
        started = []
        for index, track in landed.items():
            if track.onset_myr is None and now - track.first_myr >= self.onset_myr - 1e-6:
                track.onset_myr = now
                started.append((index, track))
        return started, landed


def _detect_step(census, world, now):
    """Every plate's sites at the point in the step where the lifecycle will sit."""
    plates = [plate for plate in world.plates if plate.node_count() > 0]
    plates_by_id = {plate.plate_id: plate for plate in plates}
    sites, polarity = [], []
    states = [(plate, boundary_state(world, plate)) for plate in plates]
    for plate, state in states:
        polarity.extend((plate, state, row) for row in polarity_fronts(world, plate, state, plates_by_id))
    census.observe_polarity(polarity, now, plates_by_id)
    for plate, state in states:
        for site in detect_sites(world, plate, state, plates_by_id, census.frozen_polarity):
            sites.append((plate, state, site))
    detections = [(plate.plate_id, site["setting"], site.get("neighbour", -1),
                   geometry.to_local(plate.frame, state["points"][site["nodes"]])
                   if site["setting"] in ("rollback", "transtension") else None)
                  for plate, state, site in sites]
    started, landed = census.tracker.update(now, detections, set(plates_by_id))
    for index, track in landed.items():
        if track.setting == "passive":
            track.axes.append(np.asarray(sites[index][2]["extension"], dtype=float))
    for index, track in started:
        plate, state, site = sites[index]
        track.row = census.event_row(world, plate, state, site, track, now)


class Census:
    """One fresh world's detections, polarity reads and production topology."""

    def __init__(self, spacing):
        self.spacing = spacing
        match = collision_polarity.FRONT_MATCH_SPACINGS * spacing
        self.tracker = Tracker(ONSET_MYR, LAPSE_MYR, match)
        # A front keeps its polarity for life, surviving the gaps a collision front record does.
        self.polarity_tracker = Tracker(0.0, collision_polarity.FRONT_EXPIRY_YEARS / 1e6, match)
        self.polarity_steps = Counter()
        self.age_steps = Counter()
        self.polarity_flips = Counter()
        self.polarity_last: dict[int, str] = {}
        self.age_reversed: set[int] = set()
        self.polarity_fronts = Counter()
        self.head_credits: list[tuple[float, int]] = []
        self.production = Counter()
        self.plate_series: list[tuple[float, int]] = []
        self.consumption_step: dict[tuple[int, int], Counter] = {}
        self.deform_position: dict[int, int] = {}
        self.consumption_order = Counter()

    def note_consumption(self, plate, neighbour_ids, is_oceanic, neighbours):
        """Oceanic cells lost at each ocean-ocean pair this step, in deform order: a read on
        whether the evidence reflects a slab or just which plate deformed first."""
        by_id = {other.plate_id: other for other in neighbours}
        ids = np.asarray(neighbour_ids)[np.asarray(is_oceanic, dtype=bool)]
        for other, count in zip(*np.unique(ids[ids >= 0], return_counts=True)):
            neighbour = by_id.get(int(other))
            if plate.crust_type == "continental" or neighbour is None or neighbour.crust_type == "continental":
                continue
            pair = tuple(sorted((plate.plate_id, int(other))))
            self.consumption_step.setdefault(pair, Counter())[plate.plate_id] += int(count)

    def settle_consumption(self):
        for pair, cells in self.consumption_step.items():
            earlier, later = sorted(pair, key=lambda plate_id: self.deform_position[plate_id])
            self.consumption_order["pair_steps"] += 1
            self.consumption_order["earlier_plate_cells"] += cells[earlier]
            self.consumption_order["later_plate_cells"] += cells[later]
            self.consumption_order["only_earlier_lost"] += cells[later] == 0
            self.consumption_order["only_later_lost"] += cells[earlier] == 0
        self.consumption_step = {}
        self.deform_position = {}

    def observe_polarity(self, rows, now, plates_by_id):
        detections = [(plate.plate_id, row["kind"], row["neighbour"],
                       geometry.to_local(plate.frame, state["points"][row["nodes"]]))
                      for plate, state, row in rows]
        _started, landed = self.polarity_tracker.update(now, detections, set(plates_by_id))
        for index, track in landed.items():
            row = rows[index][2]
            kind, verdict, age = row["kind"], row["verdict"], row["age"]
            first_read = track.track_id not in self.polarity_last
            reversed_age = track.verdict is not None and age != "tied" and age != track.verdict
            if track.verdict is None and age != "tied":
                track.verdict = age
            if not row["counted"]:
                continue
            self.polarity_steps[kind, verdict] += 1
            last = self.polarity_last.get(track.track_id)
            self.polarity_fronts[kind] += first_read
            if verdict in ("upper", "lower"):
                self.polarity_flips[kind] += bool(last) and last != verdict
                self.polarity_last[track.track_id] = verdict
                self.age_steps[kind, "agrees_consumption"] += verdict == track.verdict
                self.age_steps[kind, "against_consumption"] += track.verdict is not None and verdict != track.verdict
            else:
                self.polarity_last.setdefault(track.track_id, "")
            self.age_steps[kind, track.verdict or "undecided"] += 1
            # Starting ocean is one age, so early fronts tie; reported apart, never a limit.
            self.age_steps[kind, "late_decided" if track.verdict else "late_undecided"] += now > 100.0
            self.age_steps[kind, "reversed_fronts"] += reversed_age and track.track_id not in self.age_reversed
            if reversed_age:
                self.age_reversed.add(track.track_id)

    def frozen_polarity(self, plate, state, neighbour, nodes):
        """The polarity the ocean-ocean front through `nodes` froze, from `plate`'s side."""
        local = geometry.to_local(plate.frame, state["points"][nodes])
        track = self.polarity_tracker._claim(plate.plate_id, "ocean_ocean", neighbour, local, set())
        return "tied" if track is None or track.verdict is None else track.verdict

    def event_row(self, world, plate, state, site, track, now):
        cratons_mask = plate.collect("craton_crust_m") > 0.0
        belts = plate.collect("belt_formed_years") > 0.0
        scored, _snapshot = score_site(world, plate, state, site, cratons_mask, belts)
        seed_point = state["points"][site["seed"]]
        head = [at for at, host in self.head_credits
                if host == plate.plate_id and 0.0 <= track.first_myr - at <= plumes.PLUME_HEAT_MYR]
        row = {"setting": track.setting, "plate_id": track.plate_id, "neighbour": track.neighbour,
               "first_myr": track.first_myr, "onset_myr": now, "route_status": scored["status"],
               "site_local": [round(float(v), 4) for v in geometry.to_local(plate.frame, seed_point)],
               "plate_area_share": round(scored["plate_area_share"], 4),
               "plate_oceanic": bool(plate.crust_type != "continental"),
               "front_km": round(scored["front_km"], 1), "after_head": bool(head)}
        for key in ("upper_type", "polarity"):
            if key in site:
                row[key] = site[key]
        if "relax" in site:
            row["relax"] = round(site["relax"], 3)
        if "area_km2" in site:
            row["head_area_km2"] = round(site["area_km2"], 0)
        if "longest_arm_km" in scored:
            row["longest_arm_km"] = round(scored["longest_arm_km"], 1)
            row["growth_myr"] = round(scored["longest_arm_km"] / PROPAGATION_KM_PER_MYR, 2)
        return row

    def note_topology(self, world, previous_years):
        for logged_at, message in world.events:
            if logged_at <= previous_years:
                continue
            split, merge = SPLIT_MESSAGE.search(message), MERGE_MESSAGE.search(message)
            if split:
                self.tracker.end_plate(int(split.group(1)), "production_split")
                self.production["splits"] += 1
            elif merge:
                self.tracker.end_plate(int(merge.group(1)), "plate_merged")
                self.tracker.end_plate(int(merge.group(2)), "plate_merged")
                self.production["merges"] += 1
            elif FAILED_MESSAGE in message:
                self.production["failed_rifts"] += 1

    def finish(self, now):
        rows, pre_onset = [], Counter()
        for track, ended in self.tracker.closed + [(track, "run_end") for track in self.tracker.open]:
            if track.row is None:
                pre_onset[f"{track.setting}:{ended}"] += 1
                continue
            axes = track.axes
            turn = [float(np.degrees(np.arccos(np.clip(abs(a @ b), 0.0, 1.0)))) for a, b in zip(axes, axes[1:])]
            rows.append({**track.row, "last_myr": track.last_myr, "ended": ended,
                         "duration_myr": round(track.last_myr - track.first_myr, 2),
                         **({"axis_turn_deg_median": round(statistics.median(turn), 1)} if turn else {})})
        rows.sort(key=lambda row: (row["onset_myr"], row["plate_id"], row["setting"]))
        return rows, pre_onset


def world_fingerprint(world):
    """A digest of the state the census must leave alone, for the read-only check."""
    digest = hashlib.sha256()
    for plate in sorted(world.plates, key=lambda plate: plate.plate_id):
        digest.update(np.asarray([plate.plate_id, plate.node_count(), plate.internal_stress, plate.age_years]).tobytes())
        digest.update(np.asarray(plate.omega, dtype=float).tobytes())
        digest.update(np.asarray(plate.frame, dtype=float).tobytes())
        for name in ("crustal_thickness_m", "mantle_lithosphere_thickness_m", "craton_crust_m", "plume_heat_myr"):
            digest.update(plate.collect(name).tobytes())
    return digest.hexdigest()


def run_seed(seed, myr, step_myr, save_dir=None, detect=True):
    """One fresh world to `myr`, returning its census. With `detect` off it is the control arm of
    the read-only check."""
    world = world_mod.generate_world(seed)
    census = Census(elevation_lines.line_spacing_rad(world.node_density))
    credit = plumes._credit_head_stress

    def observed_credit(plates_data, continental, host_plate):
        credited = credit(plates_data, continental, host_plate)
        if credited:
            census.head_credits.append((world.elapsed_years / 1e6, host_plate))
        return credited

    consumption = collision_polarity.record_consumption

    def observed_consumption(world_, plate, world_xyz, is_oceanic, neighbour_ids, neighbours):
        census.note_consumption(plate, neighbour_ids, is_oceanic, neighbours)
        consumption(world_, plate, world_xyz, is_oceanic, neighbour_ids, neighbours)

    plate_class = type(world.plates[0])
    deform = plate_class.deform

    def observed_deform(plate, *args, **kwargs):
        census.deform_position[plate.plate_id] = len(census.deform_position)
        return deform(plate, *args, **kwargs)

    plumes._credit_head_stress = observed_credit
    collision_polarity.record_consumption = observed_consumption
    plate_class.deform = observed_deform
    started = time.perf_counter()
    detect_seconds = 0.0
    try:
        while world.elapsed_years < myr * 1e6 - 1.0:
            previous = world.elapsed_years
            for message in world_mod.step_world_progress(world, step_myr * 1e6):
                if detect and message[0] == "progress" and message[2] == TOPOLOGY_PHASE:
                    begun = time.perf_counter()
                    _detect_step(census, world, world.elapsed_years / 1e6)
                    detect_seconds += time.perf_counter() - begun
            census.note_topology(world, previous)
            census.settle_consumption()
            census.plate_series.append((world.elapsed_years / 1e6, len(world.plates)))
            print(json.dumps({"seed": seed, "myr": world.elapsed_years / 1e6, "plates": len(world.plates),
                              "open": len(census.tracker.open),
                              "wall_s": round(time.perf_counter() - started)}), flush=True)
    finally:
        plumes._credit_head_stress = credit
        collision_polarity.record_consumption = consumption
        plate_class.deform = deform
    events, pre_onset = census.finish(world.elapsed_years / 1e6)
    report = {"seed": seed, "myr": myr, "step_myr": step_myr, "fingerprint": world_fingerprint(world),
              "plates": census.plate_series[::10] + census.plate_series[-1:],
              "plates_max": max((count for _at, count in census.plate_series), default=0),
              "production": dict(census.production), "plume_heads_emplaced": int(world.plume_heads_emplaced),
              "continental_heads_credited": len(census.head_credits), "events": events,
              "pre_onset_tracks": dict(pre_onset),
              "polarity_front_steps": {f"{kind}:{verdict}": count for (kind, verdict), count in census.polarity_steps.items()},
              "age_polarity_steps": {f"{kind}:{name}": count for (kind, name), count in census.age_steps.items()},
              "polarity_fronts": dict(census.polarity_fronts), "polarity_flips": dict(census.polarity_flips),
              "consumption_order": dict(census.consumption_order),
              "detect_seconds": round(detect_seconds, 1), "runtime_seconds": round(time.perf_counter() - started, 1),
              "peak_rss_gb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**30, 2)}
    if save_dir is not None:
        save_dir.mkdir(parents=True, exist_ok=True)
        (save_dir / f"seed{seed}-{round(world.elapsed_years / 1e6)}myr.mbworld").write_bytes(persistence.save_world_bytes(world))
    return report


def gated_events(events):
    """The events a running rift would let start: one rift per plate at a time, the per-plate
    cooldown after it, and no repeat near an old site of that setting for RETRY_MYR. A rift holds
    its plate for the shorter of its route's growth and its driver's life."""
    kept, busy_until, sites = [], {}, []
    for event in sorted(events, key=lambda row: row["onset_myr"]):
        onset, plate = event["onset_myr"], event["plate_id"]
        if onset < busy_until.get(plate, -1.0):
            continue
        site = np.asarray(event["site_local"])
        near = np.cos(RETRY_NEAR_KM * 1e3 / lithosphere.PLANET_RADIUS_M)
        if any(old_plate == plate and setting == event["setting"] and onset - at < RETRY_MYR
               and old @ site >= near * np.linalg.norm(site) * np.linalg.norm(old)
               for old_plate, setting, at, old in sites):
            continue
        kept.append(event)
        sites.append((plate, event["setting"], onset, site))
        hold = min(event.get("growth_myr", UNROUTED_GROWTH_MYR), event["last_myr"] - onset + LAPSE_MYR)
        busy_until[plate] = onset + hold + COOLDOWN_MYR
    return kept


def _setting_row(events, gated, reports):
    span = sum(report["myr"] for report in reports) / 100.0
    per_seed = {report["seed"]: sum(event["setting_seed"] == report["seed"] for event in events) for report in reports}
    routed = [event for event in events if "growth_myr" in event]
    closed = [event for event in routed if event["ended"] != "run_end"]
    outlast = [event for event in closed if event["duration_myr"] >= event["growth_myr"]]
    return {"events": len(events), "events_per_100myr": round(len(events) / span, 2),
            "gated_per_100myr": round(len(gated) / span, 2), "events_by_seed": per_seed,
            "seeds_with_events": sum(count > 0 for count in per_seed.values()),
            "median_duration_myr": statistics.median(e["duration_myr"] for e in closed) if closed else None,
            "median_growth_myr": statistics.median(e["growth_myr"] for e in routed) if routed else None,
            "routed": len(routed), "closed_routed": len(closed), "outlast_growth": len(outlast),
            "ended": dict(Counter(event["ended"] for event in events))}


def polarity_summary(reports):
    steps, fronts, flips = Counter(), Counter(), Counter()
    for report in reports:
        steps.update(report["polarity_front_steps"])
        fronts.update(report["polarity_fronts"])
        flips.update(report["polarity_flips"])
    result = {}
    for kind in ("ocean_ocean", "continent_ocean"):
        count = {verdict: steps[f"{kind}:{verdict}"] for verdict in ("upper", "lower", "ambiguous", "none")}
        total = sum(count.values())
        decided = count["upper"] + count["lower"]
        result[kind] = {"front_steps": total, "verdicts": count, "fronts": fronts[kind], "flips": flips[kind],
                        "decided_share": round(decided / total, 3) if total else None,
                        "flip_share": round(flips[kind] / fronts[kind], 3) if fronts[kind] else None}
    oceanic = result["ocean_ocean"]
    oceanic["works"] = bool(oceanic["front_steps"] and oceanic["decided_share"] >= POLARITY_MIN_DECIDED
                            and oceanic["flip_share"] <= POLARITY_MAX_FLIPPED)
    # Against a continent the ocean must be the lower plate, so this is the evidence's truth test.
    verdicts = result["continent_ocean"]["verdicts"]
    decided = verdicts["upper"] + verdicts["lower"]
    result["continent_ocean"]["ocean_lower_share"] = round(verdicts["upper"] / decided, 3) if decided else None
    return result


def age_polarity_summary(reports):
    """The age read on ocean-ocean fronts, judged by the limits the consumption read failed."""
    steps = sum((Counter(report.get("age_polarity_steps", {})) for report in reports), Counter())
    fronts = sum(report["polarity_fronts"].get("ocean_ocean", 0) for report in reports)
    count = {name: steps[f"ocean_ocean:{name}"] for name in ("upper", "lower", "undecided")}
    total = sum(count.values())
    if not total:
        return None
    decided = count["upper"] + count["lower"]
    reversed_fronts = steps["ocean_ocean:reversed_fronts"]
    late = steps["ocean_ocean:late_decided"]
    agree, against = steps["ocean_ocean:agrees_consumption"], steps["ocean_ocean:against_consumption"]
    return {"front_steps": total, "verdicts": count, "fronts": fronts, "reversed_fronts": reversed_fronts,
            "decided_share": round(decided / total, 3), "reversed_share": round(reversed_fronts / fronts, 3) if fronts else None,
            "decided_share_after_100myr": round(late / (late + steps["ocean_ocean:late_undecided"]), 3) if late else None,
            "consumption_agreement": round(agree / (agree + against), 3) if agree + against else None,
            "works": bool(decided / total >= POLARITY_MIN_DECIDED and reversed_fronts / max(fronts, 1) <= POLARITY_MAX_FLIPPED)}


def summarize(reports):
    """Gate D1's tables from one report per seed. A worst case for plate count: every gated event
    breaks up and replaces production's splits, while merges stay at production's rate."""
    events = [{**event, "setting_seed": report["seed"]} for report in reports for event in report["events"]]
    gated = [{**event, "setting_seed": report["seed"]} for report in reports
             for event in gated_events(report["events"])]
    settings = {}
    for setting in ("active", "passive", "rollback", "transtension"):
        chosen = [event for event in events if event["setting"] == setting]
        settings[setting] = _setting_row(chosen, [event for event in gated if event["setting"] == setting], reports)
    rollback = [event for event in events if event["setting"] == "rollback"]
    split_rollback = {kind: _setting_row([event for event in rollback if event["upper_type"] == kind],
                                         [event for event in gated if event["setting"] == "rollback" and event["upper_type"] == kind],
                                         reports)
                      for kind in ("oceanic", "continental")}
    passive = [event for event in events if event["setting"] == "passive"]
    span = sum(report["myr"] for report in reports) / 100.0
    plate_rows = []
    for report in reports:
        hundreds = report["myr"] / 100.0
        final = report["plates"][-1][1]
        own = sum(event["setting_seed"] == report["seed"] for event in gated) / hundreds
        splits = report["production"].get("splits", 0) / hundreds
        routed = sum(event["setting_seed"] == report["seed"] and event["route_status"] == "separates"
                     for event in gated) / hundreds
        plate_rows.append({"seed": report["seed"], "final_plates": final, "max_plates": report["plates_max"],
                           "gated_per_100myr": round(own, 2), "production_splits_per_100myr": round(splits, 2),
                           "merges_per_100myr": round(report["production"].get("merges", 0) / hundreds, 2),
                           "worst_case_plates": round(final + (own - splits) * hundreds, 1),
                           "worst_case_plates_routed_only": round(final + (routed - splits) * hundreds, 1),
                           "compatible": bool(final + (own - splits) * hundreds < PLATE_CEILING)})
    polarity = polarity_summary(reports)
    # Oceanic-upper rollback rests on the polarity read, so it can't carry the setting if that fails.
    counted = dict(settings)
    age_polarity = age_polarity_summary(reports)
    if not (age_polarity or polarity["ocean_ocean"])["works"]:
        counted["rollback"] = split_rollback["continental"]
    return {"settings": settings, "rollback_by_upper_plate": split_rollback, "polarity": polarity,
            "age_polarity": age_polarity, "plate_count": plate_rows, "passive_after_head": sum(event["after_head"] for event in passive),
            "passive_events": len(passive), "span_100myr": span,
            "pre_onset_tracks": dict(sum((Counter(report["pre_onset_tracks"]) for report in reports), Counter())),
            "gate_settings_passing": [setting for setting, row in counted.items() if row["seeds_with_events"] >= 2]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--saves", type=Path, nargs="+")
    parser.add_argument("--seed", type=int)
    parser.add_argument("--myr", type=float, default=200.0)
    parser.add_argument("--step-myr", type=float, default=1.0)
    parser.add_argument("--no-detect", action="store_true", help="control arm for the read-only check")
    parser.add_argument("--save-dir", type=Path)
    parser.add_argument("--summarize", type=Path, nargs="+")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--image-dir", type=Path)
    args = parser.parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if args.seed is not None:
        report = run_seed(args.seed, args.myr, args.step_myr, args.save_dir, detect=not args.no_detect)
        args.out.write_text(json.dumps(report, indent=1, sort_keys=True, default=float) + "\n")
        return
    if args.summarize:
        reports = [json.loads(path.read_text()) for path in args.summarize]
        summary = summarize(reports)
        args.out.write_text(json.dumps(summary, indent=1, sort_keys=True, default=float) + "\n")
        return
    started = time.perf_counter()
    report = {"knobs": {"min_front_km": MIN_FRONT_KM, "rate_cm_yr": RATE_CM_YR, "obliquity_deg": OBLIQUITY_DEG,
                        "back_arc_km": BACK_ARC_KM, "head_area_km2": HEAD_AREA_KM2,
                        "belt_resistance": BELT_RESISTANCE, "plume_resistance": PLUME_RESISTANCE,
                        "align_penalty": ALIGN_PENALTY, "craton_edge_km": CRATON_EDGE_KM,
                        "craton_block_strength": CRATON_BLOCK_STRENGTH, "craton_cost": CRATON_COST,
                        "keel_breach_m": KEEL_BREACH_M, "craton_margin_km": CRATON_MARGIN_KM,
                        "craton_margin_cost": CRATON_MARGIN_COST, "onset_myr": ONSET_MYR,
                        "lapse_myr": LAPSE_MYR, "propagation_km_per_myr": PROPAGATION_KM_PER_MYR,
                        "retry_myr": RETRY_MYR, "retry_near_km": RETRY_NEAR_KM,
                        "age_polarity_min_gap_myr": AGE_POLARITY_MIN_GAP_MYR},
              "worlds": [survey(path, args.image_dir) for path in args.saves]}
    report["runtime_seconds"] = time.perf_counter() - started
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True, default=float) + "\n")


if __name__ == "__main__":
    main()
