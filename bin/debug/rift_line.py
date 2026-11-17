"""Rift lines for RFC 0018's D1 prototype (#888): the section pull, line placement, the strain integrator and the
two-daughter split. Pure numpy/scipy on arrays (no ``World``), so the integration and propagation checks can drive
it on manufactured plates. Not production code: #846 implements it.

RFC steps covered (Design, 2026-11-15 rift-line revision):

* step 1 (section pull): ``section_profile`` runs the line-cut running sums of proxy C and keeps, for every cell and
  each of six directions, the total pull ``P_n(s)`` across the straight section through the cell and the far-side
  area share;
* step 3 (rift line): ``candidate_corridors`` picks pull axes and corridors, ``find_line`` runs one Dijkstra search
  in a corridor, ``path_geometry`` gives the path normals and lengths, ``section_load`` the ratio ``Λ`` with its open
  branch, and ``advance_rift`` integrates strain and healing for one step in sub-steps;
* step 6 (breakthrough): ``split_masks`` removes the line, takes the two largest components, attaches strays and
  checks the daughters.

Readings of the RFC made here (also listed in ``rift_damage_mechanism``): ``Λ`` for an open line is reported as the
growth cap ``r_max Λ_grow`` (the RFC says it is always finite and ranks open lines by ``P``); a cell at its latch
strain keeps its ``N_eff`` in ``S`` for the step (it leaves ``S`` when the end-of-step latch marks it broken); edge
length between adjacent cells is counted as one per adjacency for the stray and ownership votes.
"""

from dataclasses import dataclass

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components, dijkstra
from scipy.spatial import cKDTree

import rift_damage_mechanism as mech
import rift_damage_prototype as proto
from rift_damage_mechanism import Params, RADIUS_M

from app import lithosphere, rheology, thin_shell

N_DIRECTIONS = len(proto.CUT_ANGLES)
PULL_NOISE = 1e-9  # relative to the total balanced force; keeps rounding noise from counting as a positive peak


# --- the section pull ----------------------------------------------------------------------------


@dataclass
class Sections:
    coords: np.ndarray  # (n, 2) plane coordinates (m) in the frame tangent at the plate's centroid
    normals: np.ndarray  # (6, 2) unit pull directions
    pull: np.ndarray  # (n, 6) total pull across the straight section through each cell (N), times k_cut
    far_share: np.ndarray  # (n, 6) area share of the far side of that section
    strike: np.ndarray  # (n, 6) coordinate along the section line (m)
    scale: float  # total balanced force magnitude (N): a peak below ``PULL_NOISE`` of it is not a positive pull


def plane_frame(points: np.ndarray):
    """Tangent frame at the plate's centroid and the cells' plane coordinates in metres."""
    centre = points.mean(axis=0)
    centre /= np.linalg.norm(centre)
    e_u, e_v = (axis[0] for axis in thin_shell._tangent_basis(centre[None, :]))
    return e_u, e_v, RADIUS_M * np.column_stack([points @ e_u, points @ e_v])


def section_profile(points, areas, force_n, k_cut: float) -> Sections:
    """Per direction, sort cells along ``n``; the pull across the section through cell ``x`` is the balanced force
    beyond ``x`` (plus half its own cell) projected on ``n``, times ``k_cut``; the far share is the area beyond."""
    balanced = proto.balanced_tangent_force(points, areas, force_n)
    e_u, e_v, coords = plane_frame(points)
    force_2d = np.column_stack([balanced @ e_u, balanced @ e_v])
    total = float(areas.sum())
    n = len(points)
    pull, share, strike = np.zeros((n, N_DIRECTIONS)), np.zeros((n, N_DIRECTIONS)), np.zeros((n, N_DIRECTIONS))
    normals = np.column_stack([np.cos(proto.CUT_ANGLES), np.sin(proto.CUT_ANGLES)])
    for d, normal in enumerate(normals):
        position, push = coords @ normal, force_2d @ normal
        order = np.argsort(position, kind="stable")
        beyond = np.empty(n)
        beyond[order] = np.cumsum(push[order][::-1])[::-1] - 0.5 * push[order]
        far_area = np.empty(n)
        far_area[order] = np.cumsum(areas[order][::-1])[::-1] - 0.5 * areas[order]
        pull[:, d], share[:, d] = k_cut * beyond, far_area / total
        strike[:, d] = coords @ np.array([-normal[1], normal[0]])
    return Sections(coords, normals, pull, share, strike, float(np.abs(force_2d).sum()))


@dataclass
class Corridor:
    direction: int
    cells: np.ndarray  # boolean mask over the plate's cells
    peak_pull: float  # N
    peak_share: float  # smaller-side area share at the peak section


def candidate_corridors(sections: Sections, p: Params) -> list[Corridor]:
    """Pull axes and corridors (step 3). For each direction, the peak of ``P`` over sections that leave at least
    ``a_min`` on each side; directions whose peak is positive and at least ``c_dir`` of the largest peak are
    candidates; the corridor is the contiguous run of sections (sorted along the pull) around the peak with
    ``P >= c_corr * peak`` and ``a_min`` on each side."""
    peaks = []
    for d in range(N_DIRECTIONS):
        smaller = np.minimum(sections.far_share[:, d], 1.0 - sections.far_share[:, d])
        valid = smaller >= p.a_min
        if not valid.any():
            peaks.append((-np.inf, -1, valid))
            continue
        masked = np.where(valid, sections.pull[:, d], -np.inf)
        best = int(np.argmax(masked))
        peaks.append((float(masked[best]), best, valid))
    top = max(peak for peak, _, _ in peaks)
    out = []
    noise = PULL_NOISE * sections.scale
    if top <= noise:
        return out
    for d, (peak, best, valid) in enumerate(peaks):
        if peak <= noise or peak < p.c_dir * top:
            continue
        position = sections.coords @ sections.normals[d]
        order = np.argsort(position, kind="stable")
        ok = valid[order] & (sections.pull[order, d] >= p.c_corr * peak)
        at = int(np.flatnonzero(order == best)[0])
        lo = hi = at
        while lo > 0 and ok[lo - 1]:
            lo -= 1
        while hi < len(order) - 1 and ok[hi + 1]:
            hi += 1
        cells = np.zeros(len(order), dtype=bool)
        cells[order[lo:hi + 1]] = True
        out.append(Corridor(d, cells, peak, float(min(sections.far_share[best, d], 1.0 - sections.far_share[best, d]))))
    return out


# --- the line ------------------------------------------------------------------------------------


def _trace(pred: np.ndarray, node: int) -> list[int]:
    chain = [int(node)]
    while pred[chain[-1]] >= 0:
        chain.append(int(pred[chain[-1]]))
    return chain


def find_line(graph: csr_matrix, corridor: np.ndarray, boundary: np.ndarray, strike: np.ndarray, cost: np.ndarray):
    """One Dijkstra search inside the corridor: the cheapest edge-connected chain from the corridor's plate-boundary
    cells in the lowest tenth of the strike coordinate to those in the highest tenth. ``cost`` is ``N_eff`` for
    intact cells and zero for broken cells; an edge costs the mean of its cells' costs times its length. Returns the
    chain of plate indices (ordered from the low end), or None when the corridor is not connected end to end."""
    idx = np.flatnonzero(corridor)
    if len(idx) < 2 or boundary[idx].sum() < 2:
        return None
    coo = graph[idx][:, idx].tocoo()
    local = cost[idx]
    floor = 1e-6 * max(float(local.max()), 1.0)
    weight = np.maximum(0.5 * (local[coo.row] + local[coo.col]), floor) * coo.data
    weighted = csr_matrix((weight, (coo.row, coo.col)), shape=(len(idx), len(idx)))
    edge, t = boundary[idx], strike[idx]
    low, high = np.percentile(t[edge], [10, 90])
    if low >= high:
        return None
    sources, targets = np.flatnonzero(edge & (t <= low)), np.flatnonzero(edge & (t >= high))
    if not len(sources) or not len(targets):
        return None
    dist, pred, _ = dijkstra(weighted, directed=False, indices=sources, min_only=True, return_predecessors=True)
    reached = dist[targets]
    if not np.isfinite(reached.min()):
        return None
    chain = _trace(pred, int(targets[int(np.argmin(reached))]))[::-1]
    return idx[np.asarray(chain)]


@dataclass
class Geometry:
    normal: np.ndarray  # (L, 2) in-plane unit normals (a, b) in each cell's tangent basis
    ds: np.ndarray  # (L,) line length through each cell (m)
    s: np.ndarray  # (L,) distance along the line from its first cell (m)
    length: float
    straight: float  # great-circle distance between the two end cells (m)


def path_geometry(points: np.ndarray, chain: np.ndarray, first: np.ndarray, second: np.ndarray) -> Geometry:
    """Normals and lengths of an ordered chain: the normal at a cell is perpendicular to the direction between its two
    line neighbours (in the cell's tangent basis); ``ds`` is half the centre-to-centre distances to its neighbours."""
    prev_i, next_i = np.r_[chain[0], chain[:-1]], np.r_[chain[1:], chain[-1]]
    tangent = points[next_i] - points[prev_i]
    ta, tb = np.sum(tangent * first[chain], axis=1), np.sum(tangent * second[chain], axis=1)
    norm = np.maximum(np.hypot(ta, tb), 1e-15)
    normal = np.column_stack([-tb / norm, ta / norm])
    step = RADIUS_M * np.arccos(np.clip(np.sum(points[chain[1:]] * points[chain[:-1]], axis=1), -1.0, 1.0))
    ds = np.zeros(len(chain))
    ds[:-1] += 0.5 * step
    ds[1:] += 0.5 * step
    s = np.concatenate([[0.0], np.cumsum(step)])
    straight = RADIUS_M * float(np.arccos(np.clip(points[chain[0]] @ points[chain[-1]], -1.0, 1.0)))
    return Geometry(normal, ds, s, float(step.sum()), straight)


def path_normal_load(tensor: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """``n̂ᵀ 𝐍 n̂`` (not clamped) for the line cells' tensors (L, 3) and normals (L, 2)."""
    return proto.normal_component(tensor, normal)


def section_load(pull: float, strength: float, p: Params):
    """(Λ, open) for a line with pull ``P`` and intact strength ``S``: ``P <= 0`` gives 0, ``S = 0`` with ``P > 0`` is
    open (reported at its cap ``r_max Λ_grow`` so the value is finite), otherwise ``P / S``."""
    if pull <= 0.0:
        return 0.0, False
    if strength <= 0.0:
        return p.r_max * p.lam_grow, True
    return pull / strength, False


def rank_key(pull: float, lam: float, is_open: bool):
    """Open lines rank above every other line (among themselves by ``P``); the rest by ``Λ``."""
    return (1, pull) if is_open else (0, lam)


def growth_ratio(lam: float, is_open: bool, p: Params) -> float:
    return p.r_max if is_open else min(lam / p.lam_grow, p.r_max)


def ligament_gives(pull: float, strength: float) -> bool:
    """The breakthrough test of step 6, ``P >= S`` with ``P > 0`` (an open line has ``S = 0``)."""
    return pull > 0.0 and pull >= strength


# --- band and strain integration -----------------------------------------------------------------


def band_reach(points: np.ndarray, chain: np.ndarray, eps: np.ndarray, p: Params):
    """The band at the step start: the line's cells plus every plate cell within ``w_i / 2`` of line cell ``i``, with
    ``w_i = min(3 w_0, w_0 (1 + ε_i))``. Returns the band's plate indices and the (band index, line index) pairs
    whose line cell reaches the band cell."""
    tree = cKDTree(points)
    rows, cols, members = [], [], []
    for k, cell in enumerate(chain):
        width = min(3.0 * p.w0_m, p.w0_m * (1.0 + float(eps[cell])))
        chord = 2.0 * np.sin(0.5 * (0.5 * width) / RADIUS_M)
        found = np.asarray(tree.query_ball_point(points[cell], chord * (1.0 + 1e-6)), dtype=np.int64)
        found = np.union1d(found, [cell])
        members.append(found)
        rows.append(found)
        cols.append(np.full(len(found), k))
    band = np.unique(np.concatenate(members))
    position = np.searchsorted(band, np.concatenate(rows))
    return band, position, np.concatenate(cols)


@dataclass
class LineProblem:
    """Everything ``advance_rift`` needs, frozen at the step start (RFC step 3)."""

    chain: np.ndarray  # plate indices, ordered along the line
    s: np.ndarray  # distance along the line (m)
    ds: np.ndarray  # line length through each cell (m)
    normal_load: np.ndarray  # n̂ᵀ 𝐍 n̂ at the line cells, without the plume term (N/m)
    base_quiet: np.ndarray  # k_cap N_cap M at the line cells, without plume weakening (N/m)
    base_plume: np.ndarray  # the same with the plume multiplier
    plume_left_yr: np.ndarray  # years of plume heat left in this step at each line cell (0 = none)
    intact: np.ndarray  # line cells that are not broken at the step start
    disk_rows: csr_matrix  # (len(chain), n_cells) area-weighted disk averages for the line cells
    band: np.ndarray  # plate indices of the band cells
    reach_band: np.ndarray  # band index of each (band, line) reach pair
    reach_line: np.ndarray  # line index of each reach pair
    latch: np.ndarray  # (n_cells,) strain at which each cell stops; held and broken cells hold their current strain
    capacity: "CapacityInputs | None" = None  # only for the ``live_ncap`` amendment


@dataclass
class CapacityInputs:
    """What the ``live_ncap`` amendment needs to recompute ``k_cap N_cap M`` of the line cells as they thin."""

    kind: np.ndarray
    hc: np.ndarray
    hm: np.ndarray
    lag: np.ndarray
    k_scale: np.ndarray  # k_cap (times the oceanic factor) per cell
    mult_quiet: np.ndarray
    mult_plume: np.ndarray
    eps_start: np.ndarray  # line cells' strain at the step start


@dataclass
class Advance:
    eps: np.ndarray
    warm: np.ndarray
    substeps: int
    lam_start: float
    open_start: bool
    ratio_start: float


def line_state(prob: LineProblem, eps: np.ndarray, warm: np.ndarray, t: float, p: Params):
    """(N_eff, path-normal load N⊥, plume_on) at the line cells at time ``t`` into the step, from the current
    strains: the plume load and weakening apply while ``t`` is before the cell's plume expiry."""
    plume_on = t < prob.plume_left_yr - 1e-9
    base_quiet, base_plume = prob.base_quiet, prob.base_plume
    if p.live_ncap and prob.capacity is not None:
        c = prob.capacity
        stretch = np.exp(-np.maximum(eps[prob.chain] - c.eps_start, 0.0))
        floor_c, floor_m = rheology.RIFT_CRITICAL_THICKNESS_M, lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M
        hc = np.maximum(c.hc * stretch, np.minimum(c.hc, floor_c))
        hm = np.maximum(c.hm * stretch, np.minimum(c.hm, floor_m))
        cap = mech.capacity(c.kind, hc, hm, c.lag) * c.k_scale
        base_quiet, base_plume = cap * c.mult_quiet, cap * c.mult_plume
    base = np.where(plume_on, base_plume, base_quiet)
    w = mech.softening(prob.disk_rows @ eps, prob.disk_rows @ warm, p)
    return np.maximum(base * w, 1.0), np.maximum(prob.normal_load + p.n_plume * plume_on, 0.0), plume_on


def along_line_distance(s: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Distance along the line (m) from each cell to the nearest reference cell (a boolean mask over the line)."""
    ref = np.sort(s[reference])
    at = np.searchsorted(ref, s)
    left = np.abs(s - ref[np.clip(at - 1, 0, len(ref) - 1)])
    right = np.abs(ref[np.clip(at, 0, len(ref) - 1)] - s)
    return np.minimum(left, right)


def advance_rift(prob: LineProblem, eps: np.ndarray, warm: np.ndarray, years: float, p: Params,
                 fixed_amplitude: float | None = None) -> Advance:
    """Strain, soften and heal the band for ``years`` in explicit sub-steps (RFC step 3).

    Each sub-step lasts ``min(t_left, Δε_max / max rate, time to the next plume expiry)``. At its start ``Λ``, the
    distances ``d_i`` and ``W`` are recomputed from the current strains; the line cell rate is
    ``A exp(-d_i / ℓ_tip)`` with ``A = ε̇_0 max(0, r - 1)^n``, a band cell takes the largest rate among the line cells
    that reach it, and a cell stops at its latch strain. Healing is the exact exponential over the sub-step on every
    cell of the plate. More than ``K_fail`` sub-steps in one call is an error. ``fixed_amplitude`` replaces ``A``
    (test hook for the propagation check at constant ``Λ``)."""
    eps, warm = eps.astype(float).copy(), warm.astype(float).copy()
    chain, band = prob.chain, prob.band
    expiries = np.unique(prob.plume_left_yr[(prob.plume_left_yr > 0.0) & (prob.plume_left_yr < years)])
    t, substeps = 0.0, 0
    first = None
    while t < years * (1.0 - 1e-12):
        substeps += 1
        if substeps > p.k_fail:
            raise RuntimeError(f"rift integration needed more than K_fail = {p.k_fail} sub-steps")
        n_eff, load, plume_on = line_state(prob, eps, warm, t, p)
        pull = float(np.sum(load * prob.ds))
        strength = float(np.sum((n_eff * prob.ds)[prob.intact]))
        lam, is_open = section_load(pull, strength, p)
        ratio = growth_ratio(lam, is_open, p) if pull > 0.0 else 0.0
        amplitude = fixed_amplitude if fixed_amplitude is not None else p.eps_dot0 * max(0.0, ratio - 1.0) ** p.n
        if first is None:
            first = (lam, is_open, ratio)
        marked = eps[chain] >= p.eps_cand * (1.0 - 1e-9)
        if not marked.any():
            marked = np.zeros(len(chain), dtype=bool)
            marked[int(np.argmax(load / n_eff))] = True  # the nucleus: the largest N⊥ / N_eff
        line_rate = amplitude * np.exp(-along_line_distance(prob.s, marked) / p.l_tip_m)
        band_rate = np.zeros(len(band))
        np.maximum.at(band_rate, prob.reach_band, line_rate[prob.reach_line])
        room = np.maximum(prob.latch[band] - eps[band], 0.0)
        band_rate = np.where(room > 1e-12, band_rate, 0.0)
        top = float(band_rate.max()) if len(band_rate) else 0.0
        dt = years - t
        if top > 0.0:
            dt = min(dt, p.d_eps_max / top)
        if p.crossing_stops:  # a line cell reaching eps_cand changes d_i for its neighbours: stop exactly there
            at = np.searchsorted(band, chain)
            below = (eps[chain] < p.eps_cand * (1.0 - 1e-9)) & (band_rate[at] > 0.0)
            if below.any():
                dt = min(dt, float(np.min((p.eps_cand - eps[chain][below]) / band_rate[at][below])))
        later = expiries[expiries > t + 1e-9]
        if len(later):
            dt = min(dt, float(later[0]) - t)
        d_eps = np.minimum(band_rate * dt, room)
        eps[band] += d_eps
        warm *= np.exp(-dt / p.tau_heal_yr)
        warm[band] += d_eps
        np.minimum(warm, eps, out=warm)
        t += dt
    return Advance(eps, warm, substeps, *(first if first else (0.0, False, 0.0)))


# --- breakthrough --------------------------------------------------------------------------------


def _components_without(graph: csr_matrix, removed: np.ndarray):
    keep = np.ones(graph.shape[0], dtype=bool)
    keep[removed] = False
    count, labels = connected_components(graph[keep][:, keep], directed=False)
    full = np.full(graph.shape[0], -1)
    full[keep] = labels
    return count, full


def side_votes(pattern: csr_matrix, side: np.ndarray, cells) -> np.ndarray:
    """Edges from ``cells`` to cells already on side 0 and on side 1 (adjacency counted once per edge)."""
    votes = np.zeros(2)
    for cell in np.atleast_1d(cells):
        row = pattern.indices[pattern.indptr[cell]:pattern.indptr[cell + 1]]
        votes += np.bincount(side[row][side[row] >= 0], minlength=2)
    return votes


def split_masks(graph: csr_matrix, areas: np.ndarray, chain: np.ndarray, p: Params, tie_break=None, points: np.ndarray | None = None):
    """Step 6's split. Remove the line, take the two largest remaining components as the daughters, give each line
    cell to the daughter it shares the most edges with (in rounds, so a cell with no labelled neighbour waits for one; ties through ``tie_break(side, cell)``, default side 0), give
    each other component of at most ``a_stray`` of the area to the daughter it shares the most edges with, and
    require both daughters connected and the smaller at least ``a_min`` of the area. Returns (mask of side A, None) or
    (None, reason).

    Reading of the RFC (flagged in #888): cells that were already cut off from the plate's main body before the line
    (isolated cells in the plate's own adjacency, 5 of 9 invalid lines on one development world) are islands, not
    pieces made by the line. They are left out of the stray and connectivity tests and go to the daughter of their
    nearest cell (by ``points``; side 0 without them)."""
    n = graph.shape[0]
    whole, whole_labels = connected_components(graph, directed=False)
    body = whole_labels == int(np.argmax(np.bincount(whole_labels))) if whole > 1 else np.ones(n, bool)
    islands = ~body
    count, labels = _components_without(graph, chain)
    labels = np.where(islands, -1, labels)
    live = np.unique(labels[labels >= 0])
    if len(live) < 2:
        return None, "no_separation"
    sizes = np.bincount(labels[labels >= 0], weights=areas[labels >= 0], minlength=count)
    order = [c for c in np.argsort(-sizes) if c in set(live.tolist())]
    main = order[:2]
    side = np.full(n, -1)
    for k, comp in enumerate(main):
        side[labels == comp] = k
    pattern = (graph > 0).astype(np.int8).tocsr()
    total = float(areas.sum())
    for comp in order[2:]:
        if sizes[comp] > p.a_stray * total:
            return None, "stray_too_large"
    waiting = list(chain)  # line cells in rounds: a cell with no labelled neighbour yet waits, so each side stays connected
    while waiting:
        assigned = {}
        for cell in waiting:
            votes = side_votes(pattern, side, cell)
            if votes.sum():
                assigned[cell] = int(np.argmax(votes)) if votes[0] != votes[1] else (tie_break(side, cell) if tie_break else 0)
        if not assigned:
            return None, "daughter_disconnected"
        for cell, label in assigned.items():
            side[cell] = label
        waiting = [cell for cell in waiting if cell not in assigned]
    for comp in order[2:]:
        members = np.flatnonzero(labels == comp)
        votes = side_votes(pattern, side, members)
        side[members] = int(np.argmax(votes)) if votes[0] != votes[1] else 0
    if islands.any():
        placed = np.flatnonzero(side >= 0)
        if points is not None and len(placed):
            nearest = placed[cKDTree(points[placed]).query(points[islands])[1]]
            side[islands] = side[nearest]
        else:
            side[islands] = 0
    mask = side == 0
    for part in (mask & ~islands, ~mask & ~islands):
        if not part.any() or connected_components(graph[part][:, part], directed=False)[0] != 1:
            return None, "daughter_disconnected"
    if min(areas[mask].sum(), areas[~mask].sum()) < p.a_min * total:
        return None, "daughter_small"
    return mask, None
