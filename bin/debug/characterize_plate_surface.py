#!/usr/bin/env python3
"""Issue #228 Phase 0a: deterministic characterization of a world's plate surface.

Writes one JSON document per input `.mbworld` describing how healthy the plate terrain
representation is, in terms a later `PlateWithSparseQuadPatch` can be compared against. The
sections are ordered from representation-specific to representation-neutral:

- `line_topology` -- legacy `ElevationLine` structure only: line-length histogram, one-node
  row stubs, arcs per row, along-row and across-row spacing. Has no quad equivalent; it
  exists to quantify the fragmentation the migration is meant to remove.
- `sample_cloud` -- needs only each plate's node positions and local frame: nearest-
  neighbour spacing, local refinement jumps, neighbourhood anisotropy and its alignment with
  the plate-local row direction, boundary fraction, thin tendrils, connected components.
- `implied_cells` -- the mesh-quality metrics defined in
  docs/plate-surface-baseline.md (aspect ratio, skew, folded/collapsed cells,
  non-conforming edges, cell area), measured on the structured quads *implied* by
  pairing adjacent rows. This is the bar a real quad mesh has to clear.
- `coverage` -- whole-sphere sampling: uncovered/overlapped sphere fraction, nodes inside
  another plate's outline, polygon area vs nominal node area, per-node represented area.
- `conservation` -- extensive-field totals (nominal node area), intensive-field
  distributions, categorical counts. See the field-policy table in
  docs/plate-surface-baseline.md for which
  aggregation each field gets.

Everything in the JSON is a deterministic function of the input world and the sample count
(Fibonacci sphere sampling, no RNG); floats are rounded to 6 significant figures so reruns
diff cleanly. Wall-clock timings of the derived-geometry rebuilds are written separately to
`<name>.timings.json` since they are not deterministic. `--render` also writes small Eckert IV
PNGs of the plates/platesDetail/elevation views.

`--advance-steps N` steps each loaded world N * 100 kyr before characterizing it (written as
`<name>+<N>steps.json`). Stepping is deterministic, so a behaviour-preserving refactor (issue
#228 Phase 1) must reproduce these files byte for byte.

Usage (from repo root, this repo's venv):
    backend/.venv/bin/python bin/debug/characterize_plate_surface.py WORLD.mbworld [...] \\
        --out analysis/issue228-phase0a/results [--render] [--samples 1000000] [--advance-steps 10]
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BACKEND_DIR = REPO_ROOT / "backend"
if str(BACKEND_DIR) not in sys.path:
    sys.path.insert(0, str(BACKEND_DIR))

import numpy as np  # noqa: E402
from PIL import Image  # noqa: E402
from scipy.sparse import coo_matrix  # noqa: E402
from scipy.sparse.csgraph import connected_components, dijkstra  # noqa: E402
from scipy.spatial import SphericalVoronoi, cKDTree  # noqa: E402

from app import elevation_lines, geometry, lithosphere, persistence, plates as plates_mod, render_image  # noqa: E402
from app.elevation_lines import (  # noqa: E402
    DEFRAG_CONNECT_RADIUS_MULT,
    IRREGULARITY_TOLERANCE,
    effective_is_continental_from_codes,
    line_spacing_rad,
)
from app.world import step_world  # noqa: E402

# Neighbourhood radius for the boundary/width/anisotropy metrics, in multiples of the
# world's line spacing. 1.5 takes the two along-row neighbours plus two or three nodes from
# each adjacent row on a healthy lattice (adjacent rows' thetas are not aligned, so a
# diagonal neighbour sits anywhere from 1.0 to ~1.4 spacings away).
NEIGHBOUR_RADIUS_MULT = 1.5
MAX_NEIGHBOURS = 16
# A node is on the boundary when its neighbours leave an empty angular sector wider than
# this. An interior lattice node's widest gap is ~90 degrees; a node on a straight edge sees
# 180; a convex corner more.
BOUNDARY_GAP_DEG = 135.0
# Anisotropy uses the k nearest same-plate neighbours regardless of radius, so a node on a
# stretched row still sees its (far) along-row neighbours.
ANISOTROPY_K = 8
ANISOTROPIC_RATIO = 2.0
ANISOTROPY_RATIO_CAP = 100.0
# Implied-cell edges longer than this bridge a gap rather than joining lattice neighbours --
# same multiple defragmentation uses to decide two nodes are one patch.
NONCONFORMING_EDGE_MULT = DEFRAG_CONNECT_RADIUS_MULT
# Adjacent rows further apart than this are not paired into implied cells (a missing row).
ROW_PAIR_MAX_GAP_MULT = 1.5
SMALL_COMPONENT_NODES = 10
# Two same-plate nodes closer than this (in spacings) stand in for the same patch of crust --
# "stacked" territory, e.g. two arcs of one row overlapping in theta. Reported on its own and
# excluded from the shape metrics below so they describe geometry, not duplication.
STACKED_MULT = 0.5

ELEVATION_PERCENTILES = (1, 5, 25, 50, 75, 95, 99)
SPACING_PERCENTILES = (1, 5, 50, 95, 99)


# --- small helpers ------------------------------------------------------------------------


def sig(x, digits: int = 6):
    """Round floats (recursively) to `digits` significant figures for stable diffs."""
    if isinstance(x, dict):
        return {k: sig(v, digits) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [sig(v, digits) for v in x]
    if isinstance(x, (np.floating, float)):
        x = float(x)
        if not np.isfinite(x) or x == 0.0:
            return x if np.isfinite(x) else None
        return float(f"{x:.{digits}g}")
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def pct(values: np.ndarray, qs=SPACING_PERCENTILES) -> dict:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {"n": 0}
    out = {f"p{q:02d}": float(np.percentile(values, q)) for q in qs}
    out.update(n=int(values.size), mean=float(values.mean()), min=float(values.min()), max=float(values.max()))
    return out


def fraction(mask: np.ndarray) -> float | None:
    return float(np.count_nonzero(mask)) / mask.size if mask.size else None


def fibonacci_sphere(n: int) -> np.ndarray:
    i = np.arange(n, dtype=float) + 0.5
    z = 1.0 - 2.0 * i / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    lon = np.pi * (1.0 + 5.0**0.5) * i
    return np.stack([r * np.cos(lon), r * np.sin(lon), z], axis=1)


def local_tangent(phi: np.ndarray, theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Plate-local east (increasing theta, i.e. along the row) and north unit vectors at
    (phi, theta), matching geometry.local_xyz's convention."""
    east = np.stack([-np.sin(theta), np.cos(theta), np.zeros_like(theta)], axis=-1)
    north = np.stack([-np.sin(phi) * np.cos(theta), -np.sin(phi) * np.sin(theta), np.cos(phi)], axis=-1)
    return east, north


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def plate_local_nodes(plate) -> tuple[np.ndarray, np.ndarray]:
    """(phi, theta) per node in `collect` order. Only this function and `line_topology`/
    `implied_cells` read `ElevationLine`s directly; a quad plate would supply its own."""
    lines = [line for line in plate.lines if len(line) > 0]
    if not lines:
        return np.zeros(0), np.zeros(0)
    phi = np.repeat(np.array([line.phi for line in lines], dtype=float), [len(line) for line in lines])
    return phi, np.concatenate([line.theta for line in lines])


def voronoi_areas(points_xyz: np.ndarray) -> np.ndarray:
    """Each node's spherical-Voronoi area (steradians) over *every* plate's nodes at once --
    the area a node actually stands in for, which sums to exactly 4*pi. Nodes within ~60 m
    of each other (co-located across an overlap) are merged into one generator, since
    SphericalVoronoi rejects near-duplicates, and share its area equally."""
    pts = points_xyz / np.linalg.norm(points_xyz, axis=1)[:, None]
    n = len(pts)
    pairs = cKDTree(pts).query_pairs(1e-5, output_type="ndarray")
    _, labels = connected_components(coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)), directed=False)
    first = np.unique(labels, return_index=True)[1]
    areas = SphericalVoronoi(pts[first]).calculate_areas()
    return areas[labels] / np.bincount(labels)[labels]


# --- legacy line topology -----------------------------------------------------------------


def _row_arc_intervals(lines) -> dict[float, list[tuple[float, float]]]:
    by_phi: dict[float, list[tuple[float, float]]] = defaultdict(list)
    for line in lines:
        by_phi[line.phi].append((float(line.theta[0]), float(line.theta[-1])))
    return by_phi


def _interval_stacking(ivs: list[tuple[float, float]]) -> tuple[int, float]:
    """(max number of one row's arcs covering a single theta, total theta length covered
    more than once counted with multiplicity -- sum of arc lengths minus their union).
    `_row_intervals` in plates.py merges such overlaps silently when building the outline."""
    events = sorted([(lo, 1) for lo, _ in ivs] + [(hi, -1) for _, hi in ivs], key=lambda e: (e[0], -e[1]))
    depth = best = 0
    for _, delta in events:
        depth += delta
        best = max(best, depth)
    union = 0.0
    cur_lo = cur_hi = None
    for lo, hi in sorted(ivs):
        if cur_hi is None or lo > cur_hi:
            if cur_hi is not None:
                union += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
        else:
            cur_hi = max(cur_hi, hi)
    union += cur_hi - cur_lo
    return best, sum(hi - lo for lo, hi in ivs) - union


LENGTH_BINS = ((1, 1), (2, 2), (3, 4), (5, 8), (9, 16), (17, 32), (33, 64), (65, 10**9))


def line_topology(world, spacing: float) -> dict:
    per_plate = []
    lengths_all = []
    along_all = []
    row_gap_all = []
    arcs_hist: dict[int, int] = defaultdict(int)
    overlapping_rows = max_depth = 0
    duplicated_length = 0.0
    for plate in world.plates:
        lines = [line for line in plate.lines if len(line) > 0]
        lengths = np.array([len(line) for line in lines], dtype=int)
        lengths_all.append(lengths)
        by_phi: dict[float, int] = defaultdict(int)
        for line in lines:
            by_phi[line.phi] += 1
            if len(line) >= 2:
                along_all.append(np.diff(line.theta) * np.cos(line.phi) / spacing)
        for arcs in by_phi.values():
            arcs_hist[arcs] += 1
        plate_dup = 0.0
        for phi, ivs in _row_arc_intervals(lines).items():
            depth, dup = _interval_stacking(ivs)
            if depth > 1:
                overlapping_rows += 1
                max_depth = max(max_depth, depth)
                plate_dup += dup * np.cos(phi) / spacing
        duplicated_length += plate_dup
        phis = np.sort(np.array(list(by_phi)))
        if len(phis) >= 2:
            row_gap_all.append(np.diff(phis) / spacing)
        per_plate.append(
            {
                "plate_id": plate.plate_id,
                "crust_type": plate.crust_type,
                "nodes": int(lengths.sum()),
                "lines": int(len(lengths)),
                "rows": len(by_phi),
                "one_node_lines": int(np.count_nonzero(lengths == 1)),
                "one_node_line_fraction": fraction(lengths == 1),
                "mean_nodes_per_line": float(lengths.mean()) if len(lengths) else None,
                "max_arcs_per_row": max(by_phi.values()) if by_phi else 0,
                "duplicated_row_length_spacings": plate_dup,
                "negligible_territory": bool(plate.has_negligible_territory()),
            }
        )
    lengths = np.concatenate(lengths_all) if lengths_all else np.zeros(0, dtype=int)
    along = np.concatenate(along_all) if along_all else np.zeros(0)
    row_gap = np.concatenate(row_gap_all) if row_gap_all else np.zeros(0)
    worst = sorted(per_plate, key=lambda p: (-p["one_node_lines"], p["plate_id"]))[:5]
    return {
        "lines": int(len(lengths)),
        "nodes": int(lengths.sum()),
        "one_node_lines": int(np.count_nonzero(lengths == 1)),
        "one_node_line_fraction": fraction(lengths == 1),
        "nodes_in_one_node_lines_fraction": float(np.count_nonzero(lengths == 1)) / max(1, int(lengths.sum())),
        "line_length_histogram": {
            (f"{lo}" if lo == hi else f"{lo}+" if hi > 10**8 else f"{lo}-{hi}"): int(np.count_nonzero((lengths >= lo) & (lengths <= hi)))
            for lo, hi in LENGTH_BINS
        },
        "arcs_per_row_histogram": {str(k): arcs_hist[k] for k in sorted(arcs_hist)},
        "rows": int(sum(arcs_hist.values())),
        "rows_with_overlapping_arcs": overlapping_rows,
        "max_arc_stacking_depth": max_depth,
        "duplicated_row_length_spacings": duplicated_length,
        "along_row_spacing_over_target": pct(along),
        "along_row_gaps_over_irregularity_tolerance_fraction": fraction(along > IRREGULARITY_TOLERANCE),
        "across_row_gap_over_target": pct(row_gap),
        "across_row_gaps_over_1_5_fraction": fraction(row_gap > 1.5),
        "negligible_territory_plates": [p["plate_id"] for p in per_plate if p["negligible_territory"]],
        "worst_plates_by_one_node_lines": worst,
        "per_plate": per_plate,
    }


# --- representation-neutral sample-cloud metrics ------------------------------------------


def sample_cloud(world, spacing: float) -> dict:
    """Metrics over each plate's node cloud. Positions are taken in plate-local xyz so the
    anisotropy orientation can be read against the plate-local row (east) direction."""
    nn_all, stacked_all, jump_all, ratio_all, align_all = [], [], [], [], []
    boundary_all, thin_all, depth_all = [], [], []
    per_plate = []
    components_total = isolated_total = small_component_nodes = outside_largest = 0
    for plate in world.plates:
        phi, theta = plate_local_nodes(plate)
        n = len(phi)
        if n == 0:
            continue
        xyz = geometry.local_xyz(phi, theta)
        tree = cKDTree(xyz)

        # Spacing, stacking, and refinement jumps. Neighbours closer than STACKED_MULT are
        # dropped (taking the next ones out) before the shape metrics.
        k_query = min(4 * ANISOTROPY_K + 1, n)
        dist, idx = tree.query(xyz, k=k_query)
        dist, idx = dist.reshape(n, -1)[:, 1:] / spacing, idx.reshape(n, -1)[:, 1:]
        nn = dist[:, 0] if dist.shape[1] else np.full(n, np.inf)
        nn_all.append(nn)
        stacked = nn < STACKED_MULT
        stacked_all.append(stacked)
        keep = np.isfinite(dist) & (dist >= STACKED_MULT)
        order = np.argsort(~keep, axis=1, kind="stable")[:, :ANISOTROPY_K]
        dist = np.take_along_axis(np.where(keep, dist, np.inf), order, axis=1)
        idx = np.take_along_axis(idx, order, axis=1)
        valid_k = np.isfinite(dist)
        nn_clean = np.where(valid_k[:, 0], dist[:, 0], np.inf) if dist.shape[1] else np.full(n, np.inf)
        if dist.shape[1] >= 4:
            neighbour_nn = nn_clean[np.where(valid_k[:, :4], idx[:, :4], 0)]
            with np.errstate(divide="ignore", invalid="ignore"):
                jump = np.maximum(neighbour_nn / nn_clean[:, None], nn_clean[:, None] / neighbour_nn)
            jump = np.where(valid_k[:, :4], jump, 1.0).max(1)
            jump_all.append(jump[np.isfinite(jump)])

        # Anisotropy of the k-nearest neighbourhood in the local tangent plane.
        east, north = local_tangent(phi, theta)
        if dist.shape[1] >= 4:
            offsets = xyz[np.where(valid_k, idx, np.arange(n)[:, None])] - xyz[:, None, :]
            e = np.einsum("nkj,nj->nk", offsets, east)
            nth = np.einsum("nkj,nj->nk", offsets, north)
            m = np.maximum(valid_k.sum(1), 1)
            cee, cnn, cen = (e * e).sum(1) / m, (nth * nth).sum(1) / m, (e * nth).sum(1) / m
            tr, det = cee + cnn, cee * cnn - cen * cen
            disc = np.sqrt(np.clip(tr * tr / 4 - det, 0.0, None))
            lmax, lmin = tr / 2 + disc, np.clip(tr / 2 - disc, 1e-30, None)
            # Clipped: a neighbourhood that is effectively one row wide is collinear, and
            # its unbounded ratio would swamp the mean.
            ratio = np.where(valid_k.sum(1) >= 3, np.minimum(np.sqrt(lmax / lmin), ANISOTROPY_RATIO_CAP), ANISOTROPY_RATIO_CAP)
            # cos(2*alpha) of the major axis vs local east: +1 row-aligned, -1 across rows.
            align = np.where(disc > 0, (cee - cnn) / (2 * np.where(disc > 0, disc, 1.0)), 0.0)
            ratio_all.append(ratio)
            align_all.append(align)
        else:
            ratio = np.full(n, np.inf)

        # Radius graph: boundary detection, depth, thin tendrils, components.
        r_dist, r_idx = tree.query(xyz, k=min(MAX_NEIGHBOURS + 1, n), distance_upper_bound=NEIGHBOUR_RADIUS_MULT * spacing)
        r_dist, r_idx = np.atleast_2d(r_dist), np.atleast_2d(r_idx)
        if r_dist.shape[0] != n:
            r_dist, r_idx = r_dist.T, r_idx.T
        valid = np.isfinite(r_dist) & (r_dist >= STACKED_MULT * spacing) & (r_idx < n)
        safe_idx = np.where(valid, r_idx, 0)
        offsets = xyz[safe_idx] - xyz[:, None, :]
        ang = np.arctan2(np.einsum("nkj,nj->nk", offsets, north), np.einsum("nkj,nj->nk", offsets, east))
        ang = np.where(valid, ang, np.inf)
        ang.sort(axis=1)
        counts = valid.sum(1)
        max_gap = np.full(n, 2 * np.pi)
        for c in np.unique(counts):
            rows = counts == c
            if c < 2:
                continue
            a = ang[rows, :c]
            gaps = np.diff(np.concatenate([a, a[:, :1] + 2 * np.pi], axis=1), axis=1)
            max_gap[rows] = gaps.max(1)
        boundary = max_gap > np.deg2rad(BOUNDARY_GAP_DEG)
        src = np.repeat(np.arange(n), valid.sum(1))
        dst = r_idx[valid]
        graph = coo_matrix((np.ones(len(src)), (src, dst)), shape=(n, n)).tocsr()
        if boundary.any():
            depth = dijkstra(graph, directed=False, unweighted=True, indices=np.flatnonzero(boundary), min_only=True)
        else:
            depth = np.full(n, np.inf)
        # Thin: a boundary node none of whose neighbours is interior (width <= 2 nodes).
        interior_neighbour = np.zeros(n, dtype=bool)
        np.logical_or.at(interior_neighbour, src, ~boundary[dst])
        thin = boundary & ~interior_neighbour
        boundary_all.append(boundary)
        thin_all.append(thin)
        depth_all.append(depth)

        labels = plates_mod.node_components(geometry.to_world(plate.frame, xyz), DEFRAG_CONNECT_RADIUS_MULT * spacing)
        sizes = np.bincount(labels)
        components_total += len(sizes)
        isolated_total += int(np.count_nonzero(sizes == 1))
        small_component_nodes += int(sizes[sizes < SMALL_COMPONENT_NODES].sum())
        outside_largest += int(n - sizes.max())

        # A disc of n lattice nodes has ~2*sqrt(pi*n) boundary nodes; >1 means a ragged or
        # elongated outline.
        per_plate.append(
            {
                "plate_id": plate.plate_id,
                "nodes": n,
                "stacked_fraction": fraction(stacked),
                "boundary_fraction": fraction(boundary),
                "boundary_roughness": float(np.count_nonzero(boundary)) / (2.0 * np.sqrt(np.pi * n)),
                "thin_fraction": fraction(thin),
                "anisotropic_fraction": fraction(ratio > ANISOTROPIC_RATIO),
                "components": int(len(sizes)),
                "nodes_outside_largest_component": int(n - sizes.max()),
            }
        )

    nn = np.concatenate(nn_all)
    stacked = np.concatenate(stacked_all)
    ratio = np.concatenate(ratio_all) if ratio_all else np.zeros(0)
    align = np.concatenate(align_all) if align_all else np.zeros(0)
    boundary = np.concatenate(boundary_all)
    thin = np.concatenate(thin_all)
    depth = np.concatenate(depth_all)
    aniso = ratio > ANISOTROPIC_RATIO
    finite_depth = depth[np.isfinite(depth)].astype(int)
    total_nodes = int(len(nn))
    return {
        "nodes": total_nodes,
        "nearest_neighbour_spacing_over_target": pct(nn),
        "stacked_fraction": fraction(stacked),
        "stacked_note": "nodes with a same-plate neighbour closer than STACKED_MULT spacings",
        "refinement_jump": pct(np.concatenate(jump_all)) if jump_all else {"n": 0},
        "refinement_jump_over_2_fraction": fraction(np.concatenate(jump_all) > 2.0) if jump_all else None,
        "anisotropy_ratio": pct(ratio),
        "anisotropic_fraction": fraction(aniso),
        "one_dimensional_fraction": fraction(ratio >= ANISOTROPY_RATIO_CAP),
        "anisotropic_row_alignment_mean_cos2": float(align[aniso].mean()) if aniso.any() else None,
        "anisotropic_across_row_fraction": fraction(align[aniso] < -0.5) if aniso.any() else None,
        "anisotropic_along_row_fraction": fraction(align[aniso] > 0.5) if aniso.any() else None,
        "boundary_fraction": fraction(boundary),
        "thin_fraction": fraction(thin),
        "boundary_depth_histogram": {str(d): int(c) for d, c in enumerate(np.bincount(np.clip(finite_depth, 0, 10)))},
        "boundary_depth_note": "graph hops to nearest boundary node; bucket 10 is >=10",
        "components": components_total,
        "isolated_nodes": isolated_total,
        "nodes_in_components_under_10": small_component_nodes,
        "nodes_outside_largest_component": outside_largest,
        "plates_with_multiple_components": sum(1 for p in per_plate if p["components"] > 1),
        "worst_plates_by_boundary_roughness": sorted(per_plate, key=lambda p: -p["boundary_roughness"])[:5],
        "per_plate": per_plate,
    }


# --- implied structured cells --------------------------------------------------------------


def _cell_metrics(a: np.ndarray, b: np.ndarray, d: np.ndarray, c: np.ndarray, spacing: float) -> dict[str, np.ndarray]:
    """Quad (a, b, d, c) -- a->b along the lower row, c->d along the upper -- projected to
    the tangent plane at its centroid. Counter-clockwise seen from outside the sphere."""
    centre = geometry.normalize(a + b + c + d)
    lat, lon = geometry.xyz_to_latlon(centre)
    east, north = local_tangent(lat, lon)
    corners = np.stack([a, b, d, c], axis=1) - centre[:, None, :]
    p = np.stack([np.einsum("nkj,nj->nk", corners, east), np.einsum("nkj,nj->nk", corners, north)], axis=-1) / spacing
    edges = np.roll(p, -1, axis=1) - p  # edge k: corner k -> k+1
    lengths = np.linalg.norm(edges, axis=-1)
    prev = np.roll(edges, 1, axis=1)
    cross = prev[..., 0] * edges[..., 1] - prev[..., 1] * edges[..., 0]
    dot = -(prev * edges).sum(-1)
    angle = np.degrees(np.arctan2(np.abs(cross), dot))
    area = 0.5 * np.sum(p[:, :, 0] * np.roll(p[:, :, 1], -1, axis=1) - np.roll(p[:, :, 0], -1, axis=1) * p[:, :, 1], axis=1)
    min_len = lengths.min(1)
    return {
        "aspect": np.where(min_len > 1e-9, lengths.max(1) / np.where(min_len > 1e-9, min_len, 1.0), np.inf),
        "skew_deg": np.abs(angle - 90.0).max(1),
        "folded": (cross <= 1e-12).any(1),
        "area": area,
        "max_edge": lengths.max(1),
    }


def implied_cells(world, spacing: float) -> dict:
    """Pair every lower-row segment (two consecutive nodes on one line) with the nearest
    nodes of the next row up. A collapsed cell (both lower nodes map to the same upper node)
    is the signature of a row stub or a row much sparser than its neighbour."""
    metrics = defaultdict(list)
    collapsed = unpaired_segments = segments = 0
    for plate in world.plates:
        lines = [line for line in plate.lines if len(line) > 0]
        by_phi: dict[float, list] = defaultdict(list)
        for line in lines:
            by_phi[line.phi].append(line)
        phis = sorted(by_phi)
        for lo_phi, hi_phi in zip(phis[:-1], phis[1:]):
            lower_segments = sum(max(0, len(line) - 1) for line in by_phi[lo_phi])
            segments += lower_segments
            if hi_phi - lo_phi > ROW_PAIR_MAX_GAP_MULT * spacing:
                unpaired_segments += lower_segments
                continue
            upper_theta = np.sort(np.concatenate([line.theta for line in by_phi[hi_phi]]))
            for line in by_phi[lo_phi]:
                if len(line) < 2:
                    continue
                pos = np.searchsorted(upper_theta, line.theta)
                left = np.clip(pos - 1, 0, len(upper_theta) - 1)
                right = np.clip(pos, 0, len(upper_theta) - 1)
                pick = np.where(np.abs(upper_theta[left] - line.theta) <= np.abs(upper_theta[right] - line.theta), left, right)
                up = upper_theta[pick]
                lower_xyz = geometry.local_xyz(np.full(len(line), lo_phi), line.theta)
                upper_xyz = geometry.local_xyz(np.full(len(line), hi_phi), up)
                m = _cell_metrics(lower_xyz[:-1], lower_xyz[1:], upper_xyz[1:], upper_xyz[:-1], spacing)
                is_collapsed = pick[:-1] == pick[1:]
                collapsed += int(np.count_nonzero(is_collapsed))
                for key, value in m.items():
                    metrics[key].append(value[~is_collapsed])
    if not metrics:
        return {"cells": 0}
    aspect = np.concatenate(metrics["aspect"])
    skew = np.concatenate(metrics["skew_deg"])
    folded = np.concatenate(metrics["folded"])
    area = np.concatenate(metrics["area"])
    max_edge = np.concatenate(metrics["max_edge"])
    cells = int(len(aspect))
    return {
        "lower_row_segments": segments,
        "segments_without_adjacent_row": unpaired_segments,
        "cells": cells,
        "collapsed_cells": collapsed,
        "collapsed_fraction": collapsed / max(1, cells + collapsed),
        "aspect_ratio": pct(aspect),
        "aspect_over_2_fraction": fraction(aspect > 2.0),
        "aspect_over_4_fraction": fraction(aspect > 4.0),
        "skew_deg": pct(skew),
        "skew_over_45_fraction": fraction(skew > 45.0),
        "folded_fraction": fraction(folded),
        "area_over_target": pct(area),
        "nonconforming_fraction": fraction(max_edge > NONCONFORMING_EDGE_MULT),
    }


# --- coverage and area ---------------------------------------------------------------------


def coverage(world, spacing: float, n_samples: int, timings: dict) -> dict:
    samples = fibonacci_sphere(n_samples)
    sample_sr = 4.0 * np.pi / n_samples
    nominal_sr = spacing**2
    live = [p for p in world.plates if p.node_count() > 0]
    containing = np.zeros(n_samples, dtype=np.int16)
    owner = np.full(n_samples, -1, dtype=np.int32)
    per_plate = []
    polygon_without_nodes = 0
    inside_foreign = 0
    total_nodes = 0
    t_outline = t_contains = 0.0
    spheres = {}
    for plate in live:
        pts, _ = plate.all_points_and_elevation()
        centre, radius = geometry.bounding_sphere(pts)
        spheres[plate.plate_id] = (centre, radius + 2 * spacing)
    for i, plate in enumerate(live):
        pts, _ = plate.all_points_and_elevation()
        total_nodes += len(pts)
        plate._invalidate_bounding_polygon()
        t0 = time.perf_counter()
        plate.get_bounding_polygon()
        t_outline += time.perf_counter() - t0
        centre, reach = spheres[plate.plate_id]
        near = np.flatnonzero(samples @ centre >= np.cos(min(np.pi, reach)))
        t0 = time.perf_counter()
        inside_idx = near[plate.contains_batch(samples[near])]
        t_contains += time.perf_counter() - t0
        containing[inside_idx] += 1
        owner[inside_idx] = np.where(owner[inside_idx] < 0, i, owner[inside_idx])

        # Area the outline claims vs the nominal area of the nodes it holds.
        tree = plate.get_node_kdtree()
        dist, _ = tree.query(samples[inside_idx])
        far = dist > 1.5 * spacing
        polygon_without_nodes += int(np.count_nonzero(far))
        polygon_sr = len(inside_idx) * sample_sr
        per_plate.append(
            {
                "plate_id": plate.plate_id,
                "nodes": len(pts),
                "polygon_area_over_nominal": polygon_sr / (len(pts) * nominal_sr),
                "loops": len(plates_mod._plate_outline_loops(plates_mod._row_intervals([line for line in plate.lines if len(line) > 0]))),
            }
        )

        # Nodes of this plate inside some other plate's outline (the stress-test invariant,
        # evaluated on every node rather than a 20-node sample).
        hit = np.zeros(len(pts), dtype=bool)
        for other in live:
            if other is plate:
                continue
            o_centre, o_reach = spheres[other.plate_id]
            cand = np.flatnonzero(pts @ o_centre >= np.cos(min(np.pi, o_reach)))
            if len(cand):
                hit[cand[other.contains_batch(pts[cand])]] = True
        inside_foreign += int(np.count_nonzero(hit))
        per_plate[-1]["nodes_inside_other_plate_fraction"] = fraction(hit)

    timings["outline_rebuild_s"] = t_outline
    timings["contains_batch_samples_s"] = t_contains

    all_points = np.concatenate([p.all_points_and_elevation()[0] for p in live])
    uncovered = containing == 0
    void_dist, _ = cKDTree(all_points).query(samples[uncovered]) if uncovered.any() else (np.zeros(0), None)
    overlap_tol = plates_mod.OVERLAP_TOLERANCE_MULT * spacing
    colocated = sum(int(v["overlap_mask"].sum()) for v in plates_mod.compute_node_overlap(live, overlap_tol).values())
    ratio = voronoi_areas(all_points) / nominal_sr
    loops = np.array([p["loops"] for p in per_plate])
    return {
        "samples": n_samples,
        "uncovered_sphere_fraction": fraction(uncovered),
        "uncovered_more_than_1_5_spacing_from_any_node_fraction": float(np.count_nonzero(void_dist > 1.5 * spacing)) / n_samples,
        "multiply_covered_sphere_fraction": fraction(containing >= 2),
        "max_plates_covering_one_point": int(containing.max()),
        "nodes_inside_another_plates_outline": inside_foreign,
        "nodes_inside_another_plates_outline_fraction": inside_foreign / max(1, total_nodes),
        "colocated_nodes_within_half_spacing": colocated,
        "nominal_node_area_over_sphere": total_nodes * nominal_sr / (4 * np.pi),
        "covered_area_over_nominal_node_area": (np.count_nonzero(~uncovered) * sample_sr) / (total_nodes * nominal_sr),
        "outline_area_farther_than_1_5_spacing_from_own_nodes_fraction": polygon_without_nodes / max(1, int(np.count_nonzero(~uncovered))),
        "voronoi_area_per_node_over_nominal": pct(ratio),
        "voronoi_area_under_0_5_fraction": fraction(ratio < 0.5),
        "voronoi_area_over_1_5_fraction": fraction(ratio > 1.5),
        "plates_with_holes": int(np.count_nonzero(loops > 1)),
        "total_boundary_loops": int(loops.sum()),
        "polygon_area_over_nominal_per_plate": pct(np.array([p["polygon_area_over_nominal"] for p in per_plate])),
        "per_plate": per_plate,
    }


# --- conservation / field totals ------------------------------------------------------------


EXTENSIVE_THICKNESS_FIELDS = (
    "crustal_thickness_m",
    "mantle_lithosphere_thickness_m",
    "lake_depth",
    "glacier_depth",
    "silt_depth",
    "soil_depth",
    "coal_deposit_m",
    "oil_gas_deposit_m",
    "mineral_deposit_m",
)
INTENSIVE_FIELDS = (
    "elevation",
    "crustal_thickness_m",
    "mantle_lithosphere_thickness_m",
    "channel_depth",
    "channel_width",
    "soil_mineral_content",
    "soil_organic_content",
    "divergent_age_myr",
)


def conservation(world, spacing: float) -> dict:
    area_m2 = lithosphere.node_area_m2(spacing)
    live = [p for p in world.plates if p.node_count() > 0]

    def cat(name):
        return np.concatenate([p.collect(name) for p in live])

    elevation = cat("elevation")
    codes = cat("crust_type_code")
    plate_cont = np.concatenate([np.full(p.node_count(), p.crust_type == "continental") for p in live])
    effective_cont = effective_is_continental_from_codes(codes, plate_cont)
    hc = cat("crustal_thickness_m")
    land = elevation > world.sea_level_m
    volcano = cat("is_volcano")
    active = cat("volcano_active_years_remaining") > 0
    reasons = cat("elev_change_reason").astype(int)
    created = cat("node_created_years")
    overlap_onset = cat("overlap_onset_years")
    totals = {f"{name}_km3": float(cat(name).sum()) * area_m2 / 1e9 for name in EXTENSIVE_THICKNESS_FIELDS}
    # The same totals weighted by each node's actual Voronoi area instead of the nominal
    # constant -- the gap between the two is how much the "every node is one nominal cell"
    # accounting the current conservation stats rely on misstates the real volume.
    voronoi_m2 = voronoi_areas(np.concatenate([p.all_points_and_elevation()[0] for p in live])) * lithosphere.PLANET_RADIUS_M**2
    voronoi_totals = {f"{name}_km3": float(np.dot(cat(name), voronoi_m2)) / 1e9 for name in EXTENSIVE_THICKNESS_FIELDS}
    return {
        "node_area_km2_nominal": area_m2 / 1e6,
        "sea_level_m": float(world.sea_level_m),
        "ocean_water_column_m": None if world.ocean_water_column_m is None else float(world.ocean_water_column_m),
        "plates": len(world.plates),
        "continental_plates": sum(1 for p in world.plates if p.crust_type == "continental"),
        "nodes": int(len(elevation)),
        "land_fraction_node": fraction(land),
        "land_area_km2": float(np.count_nonzero(land)) * area_m2 / 1e6,
        "extensive_totals": totals,
        "extensive_totals_voronoi_weighted": voronoi_totals,
        # stats.py's own conservation stat (plate-nominal continental), plus the node-effective split.
        "continental_plate_hc_volume_km3": float(hc[plate_cont].sum()) * area_m2 / 1e9,
        "effective_continental_hc_volume_km3": float(hc[effective_cont].sum()) * area_m2 / 1e9,
        "effective_continental_hc_volume_voronoi_km3": float(np.dot(hc[effective_cont], voronoi_m2[effective_cont])) / 1e9,
        "effective_oceanic_hc_volume_km3": float(hc[~effective_cont].sum()) * area_m2 / 1e9,
        "effective_continental_node_fraction": fraction(effective_cont),
        "elevation_land_m": pct(elevation[land], ELEVATION_PERCENTILES),
        "elevation_ocean_m": pct(elevation[~land], ELEVATION_PERCENTILES),
        "intensive": {name: pct(cat(name), ELEVATION_PERCENTILES) for name in INTENSIVE_FIELDS},
        "channel_nodes": int(np.count_nonzero(cat("channel_depth") > 0)),
        "lake_nodes": int(np.count_nonzero(cat("lake_depth") > 0)),
        "glacier_nodes": int(np.count_nonzero(cat("glacier_depth") > 0)),
        "volcano_nodes": int(np.count_nonzero(volcano)),
        "active_volcano_nodes": int(np.count_nonzero(active)),
        "crust_type_code_counts": {str(k): int(v) for k, v in enumerate(np.bincount(codes.astype(int), minlength=3))},
        "elev_change_reason_counts": {
            elevation_lines.ELEV_CHANGE_LABELS[k]: int(v) for k, v in enumerate(np.bincount(reasons, minlength=len(elevation_lines.ELEV_CHANGE_LABELS))) if v
        },
        "node_created_tracked_fraction": fraction(created >= 0),
        "node_created_years": pct(created[created >= 0], ELEVATION_PERCENTILES),
        "overlapping_nodes_now": int(np.count_nonzero(overlap_onset > 0)),
    }


# --- driver ---------------------------------------------------------------------------------


RENDER_VIEWS = ("plates", "platesDetail", "elevation")


ADVANCE_STEP_YEARS = 100_000


def characterize(path: Path, out_dir: Path, n_samples: int, render: bool, advance_steps: int = 0) -> dict:
    t0 = time.perf_counter()
    world = persistence.load_world_bytes(path.read_bytes())
    timings = {"load_s": time.perf_counter() - t0}
    t0 = time.perf_counter()
    for _ in range(advance_steps):
        step_world(world, ADVANCE_STEP_YEARS)
    timings["advance_s"] = time.perf_counter() - t0
    spacing = line_spacing_rad(world.node_density)
    name = path.stem + (f"+{advance_steps}steps" if advance_steps else "")
    result = {
        "world": {
            "file": path.name,
            "advanced_steps": advance_steps,
            "sha256": sha256(path),
            "seed": world.seed,
            "elapsed_myr": world.elapsed_years / 1e6,
            "steps_taken": world.steps_taken,
            "node_density": world.node_density,
            "climate_density": world.climate_density,
            "fluid_density": world.fluid_density,
            "line_spacing_km": spacing * elevation_lines.PLANET_RADIUS_KM,
            "plate_class": type(world.plates[0]).__name__ if world.plates else None,
        },
        "parameters": {
            "samples": n_samples,
            "neighbour_radius_mult": NEIGHBOUR_RADIUS_MULT,
            "boundary_gap_deg": BOUNDARY_GAP_DEG,
            "anisotropy_k": ANISOTROPY_K,
            "anisotropic_ratio": ANISOTROPIC_RATIO,
            "nonconforming_edge_mult": NONCONFORMING_EDGE_MULT,
            "row_pair_max_gap_mult": ROW_PAIR_MAX_GAP_MULT,
        },
    }
    for key, fn in (
        ("line_topology", lambda: line_topology(world, spacing)),
        ("sample_cloud", lambda: sample_cloud(world, spacing)),
        ("implied_cells", lambda: implied_cells(world, spacing)),
        ("coverage", lambda: coverage(world, spacing, n_samples, timings)),
        ("conservation", lambda: conservation(world, spacing)),
    ):
        t0 = time.perf_counter()
        result[key] = fn()
        timings[f"{key}_s"] = time.perf_counter() - t0
        print(f"  {key}: {timings[f'{key}_s']:.1f}s", flush=True)
    result = sig(result)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{name}.json").write_text(json.dumps(result, indent=1) + "\n")
    (out_dir / f"{name}.timings.json").write_text(json.dumps(sig(timings, 4), indent=1) + "\n")
    if render:
        for view in RENDER_VIEWS:
            # Palette-quantized: these are visual baselines for spotting row streaks, not
            # pixel references, and full-colour PNGs would bloat the checked-in results.
            png = Image.open(io.BytesIO(render_image.render_png(world, "eckert4", view, 1100, 611)))
            png.convert("RGB").quantize(colors=128, method=Image.Quantize.MEDIANCUT).save(out_dir / f"{name}.{view}.png", optimize=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("worlds", type=Path, nargs="+")
    parser.add_argument("--out", type=Path, default=REPO_ROOT / "analysis" / "issue228-phase0a" / "results")
    parser.add_argument("--samples", type=int, default=1_000_000)
    parser.add_argument("--render", action="store_true")
    parser.add_argument("--advance-steps", type=int, default=0)
    args = parser.parse_args()
    for path in args.worlds:
        print(path.name, flush=True)
        characterize(path, args.out, args.samples, args.render, args.advance_steps)


if __name__ == "__main__":
    main()
