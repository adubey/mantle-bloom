#!/usr/bin/env python3
"""Prototype a strength-following split cut as a graph minimum cut (#389, RFC 0003).

The model cuts a splitting plate along a great circle placed by the mantle flow alone, so
the cut runs through cratons wherever the flow puts it. PR #418 tilted that great circle
and found a craton-free admissible alternative for only 12% of craton-crossing cuts: a
rigid line can't detour. Real rifts follow weak lithosphere and bend around strong blocks
(the East African Rift around the Tanzania craton). This probe asks whether a cut that is
free to bend can avoid the craton while still separating the flow's two regimes.

For every cut that goes ahead in an ordinary replay (split or failed rift, via a read-only
tap on cratons.record_rift_cut, so the trajectory is unchanged), it builds the plate's
k-nearest-neighbour node graph and finds the minimum cut between the nodes beyond a
corridor of half-width W on either side of the flow cut. An edge costs its length times
(floor + mean strength of its ends). Strength is RFC 0003's cut score, max(craton
strength, keel term), so the cut prefers thin, non-cratonic lithosphere, and the floor
makes it prefer a short boundary too. W bounds how far the rift may wander from the flow's
line, as the tilt cone did for #418.

Each min cut is admissible if both halves are at least SPLIT_MIN_NODES, their poles at
least the relaxed pole separation apart, their RMS fit residual within 10% of the flow
cut's, and each half's largest connected component at least 99% of it (RFC 0010's gates).
#418's area misfit gate is reported, not applied: it only caps how far the cut moves, and
the corridor does that job here. Craton exposure is measured two ways:
- "through craton": any craton node in the band the phase 1 cut counts and a failed rift
  use (FAILED_RIFT_BAND_MULT spacings), measured from the boundary's nodes so a bent cut and
  the great circle are judged alike; and the craton area in that band, pooled with the flow
  cut kept wherever no minimum cut is admissible (the operational exposure);
- "line through craton": a craton node at either end of a link the cut crosses.
Boundary length is the summed chord length of the crossed links.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/prototype_mincut_split.py \
        --out /tmp/mincut.json
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/prototype_mincut_split.py \
        --seeds 7 --years 2e8 --out /tmp/mincut-smoke.json
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/prototype_mincut_split.py \
        --coverage-only --out /tmp/coverage.json

`--coverage-only` skips the replay and reports the generated world's craton coverage: craton
area and its share of the surface, and the largest plates' area and craton area.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import compare_cratons
import numpy as np
from app import cratons, geometry, mantle, merge_split
from app import world as world_mod
from app.elevation_lines import PLANET_RADIUS_KM, line_spacing_rad
from app.lithosphere import REFERENCE_HM_CONTINENTAL_M
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import breadth_first_order, connected_components, maximum_flow
from scipy.spatial import cKDTree

# Graph: each node links to its K nearest neighbours, dropping links longer than
# EDGE_MAX_MULT spacings so the graph doesn't bridge gaps in the plate.
K_NEIGHBOURS = 8
EDGE_MAX_MULT = 2.5
# RFC 0003's keel term: 0 at the 100 km reference, 1 at a 200 km keel.
KEEL_FULL_RATIO = 2.0
# Max-flow needs integer capacities: a unit-strength edge one spacing long carries this.
CAPACITY_SCALE = 1000.0
CORRIDORS_KM = (500.0, 1000.0, 1500.0, 2500.0)
FLOORS = (0.05, 0.25)
RESIDUAL_LIMIT = 1.10


def _strength(plate) -> np.ndarray:
    craton = cratons.strength(plate.collect("craton_crust_m"))
    hm = plate.collect("mantle_lithosphere_thickness_m")
    keel = np.clip((hm / REFERENCE_HM_CONTINENTAL_M - 1.0) / (KEEL_FULL_RATIO - 1.0), 0.0, 1.0)
    return np.maximum(craton, keel)


def _graph(points: np.ndarray, spacing: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Undirected edges (i < j) of the plate's kNN graph, with chord lengths."""
    k = min(K_NEIGHBOURS + 1, len(points))
    dist, idx = cKDTree(points).query(points, k=k)
    rows = np.repeat(np.arange(len(points)), k - 1)
    cols, dist = idx[:, 1:].ravel(), dist[:, 1:].ravel()
    keep = dist < EDGE_MAX_MULT * spacing
    i, j = np.minimum(rows, cols)[keep], np.maximum(rows, cols)[keep]
    pairs, first = np.unique(i * len(points) + j, return_index=True)
    return pairs // len(points), pairs % len(points), dist[keep][first]


def _min_cut(n: int, i, j, length, strength, floor, spacing, source_set, sink_set) -> np.ndarray:
    """Source-side mask of the minimum cut separating source_set from sink_set."""
    cost = length / spacing * (floor + 0.5 * (strength[i] + strength[j]))
    cap = np.maximum(1, np.rint(cost * CAPACITY_SCALE)).astype(np.int32)
    src, snk = n, n + 1
    big = np.int32(2**28)
    s_nodes, t_nodes = np.flatnonzero(source_set), np.flatnonzero(sink_set)
    rows = np.concatenate([i, j, np.full(len(s_nodes), src), t_nodes])
    cols = np.concatenate([j, i, s_nodes, np.full(len(t_nodes), snk)])
    data = np.concatenate([cap, cap, np.full(len(s_nodes), big), np.full(len(t_nodes), big)]).astype(np.int32)
    graph = csr_matrix(coo_matrix((data, (rows, cols)), shape=(n + 2, n + 2)))
    graph.sum_duplicates()
    flow = maximum_flow(graph, src, snk).flow
    residual = (graph - flow).tocsr()
    residual.data = (residual.data > 0).astype(np.int8)
    residual.eliminate_zeros()
    reached = breadth_first_order(residual, src, directed=True, return_predecessors=False)
    mask = np.zeros(n + 2, dtype=bool)
    mask[reached] = True
    return mask[:n]


def _boundary_band(points, i, j, mask, band_chord) -> np.ndarray:
    """Nodes within band_chord of a node on the cut (an edge whose ends differ in side)."""
    crossing = mask[i] != mask[j]
    edge_nodes = np.unique(np.concatenate([i[crossing], j[crossing]]))
    if len(edge_nodes) == 0:
        return np.zeros(len(points), dtype=bool)
    dist, _ = cKDTree(points[edge_nodes]).query(points, distance_upper_bound=band_chord)
    return np.isfinite(dist)


def _halves(points, velocities, mask, crossing_length) -> dict:
    out = {"count_min": int(min(mask.sum(), (~mask).sum())), "crossing_length": float(crossing_length)}
    if out["count_min"] < 3:
        return {**out, "pole_sep": 0.0, "rms": np.inf}
    rss, poles = 0.0, []
    for half in (mask, ~mask):
        pole = mantle.fit_euler_pole(points[half], velocities[half])
        poles.append(pole)
        rss += merge_split._fit_residual_rms(points[half], velocities[half], pole) ** 2 * half.sum()
    return {**out, "pole_sep": float(np.linalg.norm(poles[0] - poles[1])), "rms": float(np.sqrt(rss / len(points)))}


def _components(n, i, j, mask) -> list[float]:
    """Largest connected component's share of each half."""
    shares = []
    for half in (mask, ~mask):
        keep = half[i] & half[j]
        idx = np.flatnonzero(half)
        remap = np.full(n, -1)
        remap[idx] = np.arange(len(idx))
        g = coo_matrix((np.ones(keep.sum()), (remap[i[keep]], remap[j[keep]])), shape=(len(idx), len(idx)))
        _, labels = connected_components(g, directed=False)
        shares.append(float(np.bincount(labels).max() / len(idx)) if len(idx) else 0.0)
    return shares


def _relaxed_gates(world, plate, points) -> tuple[float, int]:
    """maybe_split_plate's relaxed pole-separation threshold and min_nodes for this plate."""
    _, radius_rad = geometry.bounding_sphere(points)
    relax = max(min(1.0, radius_rad / merge_split.SPLIT_SIZE_CERTAIN_RIFT_RAD),
                min(1.0, plate.internal_stress / merge_split.STRESS_CERTAIN_RIFT))
    return merge_split.SPLIT_MIN_POLE_SEPARATION * (1.0 - relax), max(1, round(merge_split.SPLIT_MIN_NODES * world.node_density))


def _evaluate(world, plate, points, cut_normal, kind) -> dict:  # noqa: C901 - one flat report per cut
    started = time.perf_counter()
    spacing = line_spacing_rad(world.node_density)
    band_chord = 2.0 * np.sin(0.5 * merge_split.FAILED_RIFT_BAND_MULT * spacing)
    velocities = mantle.flow_at(points, world.mantle_centers)
    areas = plate.accounting_areas_m2(spacing)
    craton = plate.collect("craton_crust_m") > 0.0
    strength = _strength(plate)
    threshold, min_nodes = _relaxed_gates(world, plate, points)
    n = len(points)
    i, j, length = _graph(points, spacing)
    proj = points @ cut_normal
    flow_mask = proj > 0.0
    flow_band = _boundary_band(points, i, j, flow_mask, band_chord)
    flow_crossing = flow_mask[i] != flow_mask[j]
    flow_fit = _halves(points, velocities, flow_mask, length[flow_crossing].sum())
    row = {
        "kind": kind, "myr": world.elapsed_years / 1e6, "plate": plate.plate_id, "nodes": n,
        "plate_area_mkm2": float(areas.sum() / 1e12),
        "plate_craton_area_mkm2": float(areas[craton].sum() / 1e12),
        "flow": {"through_craton": bool(craton[flow_band].any()),
                 "line_through_craton": bool(craton[np.concatenate([i[flow_crossing], j[flow_crossing]])].any()),
                 "through_craton_greatcircle": bool(craton[np.abs(proj) < np.sin(merge_split.FAILED_RIFT_BAND_MULT * spacing)].any()),
                 "band_area_mkm2": float(areas[flow_band].sum() / 1e12),
                 "band_craton_area_mkm2": float(areas[flow_band & craton].sum() / 1e12),
                 "band_strength_area_mkm2": float((strength * areas)[flow_band].sum() / 1e12),
                 **flow_fit},
        "threshold": threshold, "min_nodes": min_nodes, "variants": {},
    }
    for corridor_km in CORRIDORS_KM:
        edge = np.sin(corridor_km / PLANET_RADIUS_KM)
        source_set, sink_set = proj >= edge, proj <= -edge
        for floor in FLOORS:
            tag = f"w{corridor_km:g}_f{floor:g}"
            if not source_set.any() or not sink_set.any():
                row["variants"][tag] = {"skipped": "corridor wider than plate"}
                continue
            t0 = time.perf_counter()
            mask = _min_cut(n, i, j, length, strength, floor, spacing, source_set, sink_set)
            cut_s = time.perf_counter() - t0
            band = _boundary_band(points, i, j, mask, band_chord)
            crossing = mask[i] != mask[j]
            fit = _halves(points, velocities, mask, length[crossing].sum())
            comp = _components(n, i, j, mask)
            v = {
                "through_craton": bool(craton[band].any()),
                "band_craton_area_mkm2": float(areas[band & craton].sum() / 1e12),
                "band_strength_area_mkm2": float((strength * areas)[band].sum() / 1e12),
                "misfit": float(areas[mask != flow_mask].sum() / areas.sum()),
                "residual_ratio": fit["rms"] / flow_fit["rms"] if flow_fit["rms"] > 0 else np.inf,
                # Summed chord length of the links each cut crosses: the quantity the cut cost
                # weights, so the ratio compares boundary lengths, not link counts (#421 review).
                "length_ratio": fit["crossing_length"] / flow_fit["crossing_length"] if flow_fit["crossing_length"] > 0 else np.inf,
                "line_through_craton": bool(craton[np.concatenate([i[crossing], j[crossing]])].any()),
                "largest_component_share": comp,
                "size_ok": fit["count_min"] >= min_nodes,
                "poles_ok": fit["pole_sep"] >= threshold,
                "cut_s": cut_s,
            }
            v["residual_ok"] = v["residual_ratio"] <= RESIDUAL_LIMIT
            v["connected"] = min(comp) >= 0.99
            v["admissible"] = v["size_ok"] and v["poles_ok"] and v["residual_ok"] and v["connected"]
            row["variants"][tag] = v
    row["eval_s"] = time.perf_counter() - started
    return row


def _summarize(rows: list[dict]) -> dict:
    through = [r for r in rows if r["flow"]["through_craton"]]
    out = {"cuts": len(rows), "flow_through_craton": len(through),
           "flow_line_through_craton": sum(r["flow"]["line_through_craton"] for r in rows),
           "flow_through_craton_greatcircle": sum(r["flow"]["through_craton_greatcircle"] for r in rows),
           "eval_s_mean": float(np.mean([r["eval_s"] for r in rows])) if rows else None, "variants": {}}
    flow_craton = sum(r["flow"]["band_craton_area_mkm2"] for r in rows)
    for tag in rows[0]["variants"] if rows else ():
        vs = [(r, r["variants"][tag]) for r in through if "skipped" not in r["variants"][tag]]
        alls = [r["variants"][tag] for r in rows if "skipped" not in r["variants"][tag]]
        out["variants"][tag] = {
            "evaluated_through": len(vs),
            "craton_free": sum(not v["through_craton"] for _, v in vs),
            "craton_free_admissible": sum(not v["through_craton"] and v["admissible"] for _, v in vs),
            "line_craton_free": sum(not v["line_through_craton"] for _, v in vs),
            "line_craton_free_admissible": sum(not v["line_through_craton"] and v["admissible"] for _, v in vs),
            "fails_connected": sum(not v["connected"] for _, v in vs),
            # Pooled over craton-crossing cuts: the craton a rift's operational band would
            # thin, minimum cut over flow cut, keeping the flow cut where none is admissible.
            "band_craton_area_vs_flow_through_admissible": (
                sum(v["band_craton_area_mkm2"] if v["admissible"] else r["flow"]["band_craton_area_mkm2"] for r, v in vs)
                / max(sum(r["flow"]["band_craton_area_mkm2"] for r, _ in vs), 1e-300)) if vs else None,
            "fails_size": sum(not v["size_ok"] for _, v in vs),
            "fails_poles": sum(not v["poles_ok"] for _, v in vs),
            "fails_residual": sum(not v["residual_ok"] for _, v in vs),
            "admissible_all_cuts": sum(v["admissible"] for v in alls),
            "band_craton_area_vs_flow": (sum(v["band_craton_area_mkm2"] for v in alls) / flow_craton) if flow_craton else None,
            "median_residual_ratio": float(np.median([v["residual_ratio"] for v in alls])) if alls else None,
            "median_misfit": float(np.median([v["misfit"] for v in alls])) if alls else None,
            "median_length_ratio": float(np.median([v["length_ratio"] for v in alls])) if alls else None,
            "max_cut_s": max((v["cut_s"] for v in alls), default=None),
        }
    return out


def _run(seed: int, years: float) -> dict:
    rows: list[dict] = []
    original = cratons.record_rift_cut

    def tapped(world, plate, points, cut_normal, band_sin, kind):
        rows.append(_evaluate(world, plate, np.asarray(points, dtype=float), np.asarray(cut_normal, dtype=float), kind))
        return original(world, plate, points, cut_normal, band_sin, kind)

    cratons.record_rift_cut = tapped
    started = time.perf_counter()
    try:
        compare_cratons.run_job(seed, "on", years, 100e6, 5e6, {})
    finally:
        cratons.record_rift_cut = original
    return {"seed": seed, "wall_s": time.perf_counter() - started, "summary": _summarize(rows), "cuts": rows}


def _coverage(seed: int) -> dict:
    world = world_mod.generate_world(seed)
    spacing = line_spacing_rad(world.node_density)
    plates = []
    for plate in world.plates:
        areas = plate.accounting_areas_m2(spacing)
        plates.append((float(areas.sum() / 1e12), float(areas[plate.collect("craton_crust_m") > 0.0].sum() / 1e12)))
    plates.sort(reverse=True)
    surface = sum(area for area, _ in plates)
    craton = sum(c for _, c in plates)
    return {"seed": seed, "plates": len(plates), "surface_mkm2": surface, "craton_mkm2": craton,
            "craton_share_of_surface": craton / surface, "largest_plates_area_craton_mkm2": plates[:5]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default="7,11,2024")
    parser.add_argument("--years", type=float, default=1e9)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--coverage-only", action="store_true")
    args = parser.parse_args()
    if args.coverage_only:
        coverage = [_coverage(seed) for seed in map(int, args.seeds.split(","))]
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(coverage, indent=1) + "\n")
        print(json.dumps(coverage, indent=1))
        return
    app_root = Path(merge_split.__file__).resolve().parents[2]
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=app_root, text=True).strip()
    result = {"app_path": str(Path(merge_split.__file__).resolve()), "model_commit": commit, "years": args.years,
              "corridors_km": CORRIDORS_KM, "floors": FLOORS, "runs": []}
    for seed in map(int, args.seeds.split(",")):
        result["runs"].append(_run(seed, args.years))
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=1, default=float) + "\n")
        print(json.dumps({"seed": seed, **result["runs"][-1]["summary"]}, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
