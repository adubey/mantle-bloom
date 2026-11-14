"""Offline prototype of RFC 0018's damage-band mechanism (D1, #888). Not production code: #846 implements it.

Usage: imported by ``check_rift_damage_prototype.py`` (unit-level checks) and ``run_rift_damage_d1.py`` (the
D1 driver). ``install(params)`` registers the three damage surface fields in this process, wraps
``merge_split.apply_topology_changes`` via ``maybe_split_plate`` and the margin-insertion path, and returns
a ``Mechanism`` whose ``report`` holds what the recorder needs. ``uninstall()`` restores the originals.

The mechanism follows the RFC's Design, steps 1-6:

* load: proxy C, ``rift_damage_prototype.cut_tensor`` scaled by ``k_cut`` (plume load added at expiry);
* strength: ``k_cap N_cap M W`` with non-local softening over a disk of diameter ``w_0``;
* strain: sub-stepped, with exact exponential healing, capped at ``K_sub``;
* thinning: flank transfer of the removed volume (``_allocate_stretch``'s rule, run in reverse);
* the broken latch, with no interior eruption;
* breakthrough: two Dijkstra legs through the band, the path-normal integrated ligament rule, and a
  mask-based split.

Readings of the RFC the prototype had to choose (each is a ``Params`` field or marked in the code, and is
listed in ``analysis/issue888/README.md`` for review):

* ``live_self_term``: the RFC freezes disk-averaged strains for the step but also sub-steps because "the rate
  depends nonlinearly on ε through W". Frozen throughout, W would be constant and sub-stepping pointless. So
  the neighbours' share of each disk average is frozen at the step start and the cell's own share follows its
  own strain. ``live_self_term=False`` freezes everything.
* ``N_cap`` is frozen at the step start (the RFC's operator split); thinning enters it at the next step.
* margin insertion inherits the area-weighted mean of the donors within ``K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION``
  hops times the stretch share, without ``_allocate_stretch``'s craton willingness weights.
* eruption resets are detected after the fact: a cell whose column became exactly the reference column
  and which was not before.
"""

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parent), str(HERE.parents[2] / "backend")]

import replay_gate2_capacity as gate2  # noqa: E402

from app import lithosphere, rheology  # noqa: E402
from app.surface_fields import SURFACE_FIELDS, RemapClass, _field  # noqa: E402

RADIUS_M = lithosphere.PLANET_RADIUS_M
YEAR_PER_MYR = 1.0e6
FIELDS = ("rift_strain", "rift_warm_strain", "rift_broken")
MULTIPLIER_FLOOR = 0.12
MULT_BELT, MULT_VOLCANO, MULT_PLUME, MULT_FAULT = 0.65, 0.65, 0.65, 0.60


@dataclass(frozen=True)
class Params:
    eps_dot0: float = 8e-8  # 1/yr
    n: float = 1.0
    eps_soft: float = 1.0
    s_fault: float = 0.40
    s_warm: float = 0.40
    s_max: float = 0.88
    tau_heal_yr: float = 62.8e6
    k_cap: float = 1.0
    k_cut: float = 1.33
    n_plume: float = 0.5e12  # N/m
    w0_m: float = 125e3
    eps_cand: float = 0.1
    eps_break: float = 2.0
    t_fail_yr: float = 5e6
    flank_rows: int = 2
    flank_rows_max: int = 4
    k_sub: int = 64
    d_eps_max: float = 0.05
    live_self_term: bool = True
    oceanic_latch_on_break_only: bool = False  # the predeclared oceanic fallback
    oceanic_k_cap_factor: float = 1.0


# --- fields --------------------------------------------------------------------------------------


def register_fields() -> None:
    """Add the damage fields to the surface registry for this process (idempotent). ``rift_strain`` is
    COUNTDOWN, which coarsening and merges take the maximum of, so a coarsened quartet keeps its most damaged
    cell; ``rift_warm_strain`` is intensive and clamped to ``rift_strain`` after every pass; ``rift_broken``
    is a provenance flag (any)."""
    SURFACE_FIELDS.setdefault("rift_strain", _field(float, 0.0, RemapClass.COUNTDOWN))
    SURFACE_FIELDS.setdefault("rift_warm_strain", _field(float, 0.0, RemapClass.INTENSIVE, coupled_to=("rift_strain",)))
    SURFACE_FIELDS.setdefault("rift_broken", _field(bool, False, RemapClass.BOOLEAN_PROVENANCE))


# --- strength ------------------------------------------------------------------------------------

HC_GRID = np.arange(2_000.0, 90_001.0, 2_000.0)
HM_GRID = np.arange(2_000.0, 200_001.0, 8_000.0)
LAG_GRID = np.array([0.0, 50.0, 100.0, 200.0, 400.0])
KINDS = ("continental", "oceanic_young_ridge", "oceanic_mature")
_TABLE_CACHE = Path.home() / "mantle-bloom-captures" / "issue888" / "ncap_table.npz"
_table: np.ndarray | None = None


def capacity_table() -> np.ndarray:
    """N_cap (N/m) on a (kind, Hc, Hm, Moho lag) grid from the #487 yield envelope, built once and cached."""
    global _table
    if _table is not None:
        return _table
    if _TABLE_CACHE.exists():
        data = np.load(_TABLE_CACHE)
        if data["hc"].shape == HC_GRID.shape and data["hm"].shape == HM_GRID.shape and data["lag"].shape == LAG_GRID.shape:
            _table = data["table"]
            return _table
    table = np.zeros((len(KINDS), len(HC_GRID), len(HM_GRID), len(LAG_GRID)))
    for ki, kind in enumerate(KINDS):
        for hi, hc in enumerate(HC_GRID):
            for mi, hm in enumerate(HM_GRID):
                for li, lag in enumerate(LAG_GRID):
                    table[ki, hi, mi, li] = gate2._column_capacity(kind, float(hc), float(hm), float(lag))
    _TABLE_CACHE.parent.mkdir(parents=True, exist_ok=True)
    np.savez(_TABLE_CACHE, table=table, hc=HC_GRID, hm=HM_GRID, lag=LAG_GRID)
    _table = table
    return table


def capacity(kind_index: np.ndarray, hc, hm, lag) -> np.ndarray:
    """Multilinear interpolation of ``capacity_table`` (clamped to the grid) for each cell."""
    table = capacity_table()
    out = np.empty(len(hc))
    for ki in range(len(KINDS)):
        sel = kind_index == ki
        if not np.any(sel):
            continue
        axes = []
        for value, grid in ((hc[sel], HC_GRID), (hm[sel], HM_GRID), (lag[sel], LAG_GRID)):
            x = np.clip(value, grid[0], grid[-1])
            lo = np.clip(np.searchsorted(grid, x, side="right") - 1, 0, len(grid) - 2)
            axes.append((lo, np.clip((x - grid[lo]) / (grid[lo + 1] - grid[lo]), 0.0, 1.0)))
        acc = np.zeros(int(sel.sum()))
        for corner in range(8):
            weight, index = np.ones(len(acc)), []
            for axis, (lo, frac) in enumerate(axes):
                bit = (corner >> axis) & 1
                weight = weight * (frac if bit else 1.0 - frac)
                index.append(lo + bit)
            acc += weight * table[ki][index[0], index[1], index[2]]
        out[sel] = acc
    return out


def softening(eps_bar, warm_bar, p: Params) -> np.ndarray:
    """W = 1 - min(s_max, s_fault φ(ε̄) + s_warm φ(ε̄_warm)), φ(ε) = min(1, ε / ε_soft); in [1 - s_max, 1]."""
    phi = lambda e: np.minimum(1.0, np.maximum(e, 0.0) / p.eps_soft)  # noqa: E731
    return 1.0 - np.minimum(p.s_max, p.s_fault * phi(eps_bar) + p.s_warm * phi(warm_bar))


def inherited_multiplier(belt, volcano, plume, fault) -> np.ndarray:
    """Product of the accepted inherited-weakness multipliers, clamped to [0.12, 1]."""
    m = np.ones(len(belt))
    for mask, factor in ((belt, MULT_BELT), (volcano, MULT_VOLCANO), (plume, MULT_PLUME), (fault, MULT_FAULT)):
        m = np.where(mask, m * factor, m)
    return np.clip(m, MULTIPLIER_FLOOR, 1.0)


# --- strain --------------------------------------------------------------------------------------


@dataclass
class StrainResult:
    eps: np.ndarray
    warm: np.ndarray
    capped: np.ndarray  # cells that hit K_sub
    substeps: int


def integrate_strain(load, base, eps, warm, nb_eps, nb_warm, self_weight, duration, latch_eps, p: Params) -> StrainResult:
    """Advance strain and healing for ``duration`` (yr, per cell) under a frozen load.

    ``base`` is ``k_cap N_cap M`` (N/m), so ``N_eff = base W``. The disk averages seen by a cell are
    ``nb + self_weight × its own value`` (neighbour part frozen, own part live unless ``live_self_term`` is
    off, in which case ``self_weight`` is taken as 0 and ``nb`` already holds the whole average). A cell stops at
    ``latch_eps`` (the strain at which it breaks). Healing is the exact exponential of each sub-step.

    Cells that need more than ``K_sub`` sub-steps are advanced straight to their latch strain and flagged:
    with ``Δε_max`` per sub-step ``K_sub`` sub-steps cover more than ``ε_break``, so the cap only binds on a
    runaway rate, which is a cell about to break anyway."""
    eps, warm = eps.astype(float).copy(), warm.astype(float).copy()
    remaining = np.maximum(np.broadcast_to(duration, eps.shape).astype(float), 0.0).copy()
    capped = np.zeros(len(eps), dtype=bool)
    active = (remaining > 0.0)
    substeps = 0
    for substeps in range(1, p.k_sub + 1):
        idx = np.flatnonzero(active)
        if not len(idx):
            substeps -= 1
            break
        w = softening(nb_eps[idx] + self_weight[idx] * eps[idx], nb_warm[idx] + self_weight[idx] * warm[idx], p)
        n_eff = np.maximum(base[idx] * w, 1.0)
        ratio = np.maximum(load[idx] / n_eff - 1.0, 0.0)
        rate = p.eps_dot0 * ratio ** p.n
        at_latch = eps[idx] >= latch_eps[idx] - 1e-12
        rate = np.where(at_latch, 0.0, rate)
        dt = np.where(rate > 0.0, np.minimum(remaining[idx], p.d_eps_max / np.maximum(rate, 1e-300)), remaining[idx])
        d_eps = np.minimum(rate * dt, np.maximum(latch_eps[idx] - eps[idx], 0.0))
        eps[idx] += d_eps
        warm[idx] = warm[idx] * np.exp(-dt / p.tau_heal_yr) + d_eps
        warm[idx] = np.minimum(warm[idx], eps[idx])
        remaining[idx] -= dt
        active[idx] = remaining[idx] > 1e-6 * np.maximum(duration if np.ndim(duration) == 0 else duration[idx], 1.0)
    over = np.flatnonzero(active)
    if len(over):
        capped[over] = True
        eps[over] = latch_eps[over]
        warm[over] = np.minimum(warm[over] + 0.0, eps[over])
    return StrainResult(eps, warm, capped, substeps)


# --- plate graph and disk averages ---------------------------------------------------------------


def plate_graph(plate, points):
    """Symmetric CSR adjacency with centre distances (m) as data, plus the plain pair list."""
    adjacency = plate.adjacency()
    n = len(points)
    counts = np.diff(adjacency.offsets)
    i = np.repeat(np.arange(n), counts)
    j = np.asarray(adjacency.neighbours)
    length = RADIUS_M * np.arccos(np.clip(np.sum(points[i] * points[j], axis=1), -1.0, 1.0))
    return csr_matrix((np.maximum(length, 1.0), (i, j)), shape=(n, n)), i, j


def disk_average_matrix(points, areas, radius_m):
    """Row-stochastic area-weighted average over the cells within ``radius_m`` (inclusive of itself)."""
    chord = 2.0 * np.sin(0.5 * radius_m / RADIUS_M)
    tree = cKDTree(points)
    neighbours = tree.query_ball_point(points, chord * (1.0 + 1e-6))
    rows = np.repeat(np.arange(len(points)), [len(v) for v in neighbours])
    cols = np.concatenate([np.asarray(v, dtype=np.int64) for v in neighbours])
    weights = areas[cols]
    totals = np.bincount(rows, weights=weights, minlength=len(points))
    return csr_matrix((weights / totals[rows], (rows, cols)), shape=(len(points), len(points)))


# --- flank transfer ------------------------------------------------------------------------------


def flank_reach(adjacency: csr_matrix, donors: np.ndarray, recipient_ok: np.ndarray, rows: int, rows_max: int):
    """For each donor the recipients within ``rows`` hops, widening one hop at a time to ``rows_max`` when it
    has none. Returns (csr donors x cells reach to recipients, boolean per donor: has a recipient)."""
    n = adjacency.shape[0]
    pattern = (adjacency > 0).astype(np.int8)
    visited = csr_matrix((np.ones(len(donors), dtype=np.int8), (np.arange(len(donors)), donors)), shape=(len(donors), n))
    frontier = visited.copy()
    result = csr_matrix((len(donors), n), dtype=np.int8)
    settled = np.zeros(len(donors), dtype=bool)
    ok = csr_matrix(np.where(recipient_ok, 1, 0).astype(np.int8)[None, :])
    for hop in range(1, rows_max + 1):
        frontier = (frontier @ pattern)
        frontier.data[:] = 1
        visited = (visited + frontier)
        visited.data[:] = 1
        if hop < rows:
            continue
        reach = visited.multiply(ok).tocsr()
        found = np.asarray(reach.sum(axis=1)).ravel() > 0
        take = found & ~settled
        if np.any(take):
            result = result + csr_matrix(reach.multiply(take[:, None].astype(np.int8)))
            settled |= take
        if settled.all():
            break
    return result.tocsr(), settled


@dataclass
class ThinningOutcome:
    hc: np.ndarray
    hm: np.ndarray
    material: np.ndarray
    craton: np.ndarray
    cover: np.ndarray
    cover_material: np.ndarray
    plateau: np.ndarray
    arc: np.ndarray
    moment: np.ndarray
    craton_lost_m3: float
    withheld: np.ndarray  # donors whose thinning was withheld


def thin_with_flank_transfer(areas, adjacency, d_eps, broken, columns: dict, p: Params) -> ThinningOutcome:
    """Pure-shear thinning of the straining cells by exp(-Δε), floored, with the removed volume moved to the
    non-straining, unbroken flank cells of the same plate. ``columns`` holds hc, hm, material, craton, cover,
    cover_material, plateau, arc, moment (per cell, thickness units). Volumes are computed from the post-clamp
    thicknesses, so a floored donor gives up only what it lost. A donor with no recipient within
    ``flank_rows_max`` hops is withheld: it neither thins nor gives volume."""
    out = {k: v.astype(float).copy() for k, v in columns.items()}
    hc, hm = out["hc"], out["hm"]
    straining = d_eps > 0.0
    donors = np.flatnonzero(straining)
    withheld = np.zeros(len(areas), dtype=bool)
    craton_lost = 0.0
    if not len(donors):
        return ThinningOutcome(**out, craton_lost_m3=0.0, withheld=withheld)
    ratio = np.exp(-d_eps[donors])
    ratio_c = np.where(hc[donors] > rheology.RIFT_CRITICAL_THICKNESS_M,
                       np.maximum(ratio, rheology.RIFT_CRITICAL_THICKNESS_M / np.maximum(hc[donors], 1.0)), 1.0)
    ratio_m = np.where(hm[donors] > lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M,
                       np.maximum(ratio, lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M / np.maximum(hm[donors], 1.0)), 1.0)
    recipient_ok = ~straining & ~broken
    reach, settled = flank_reach(adjacency, donors, recipient_ok, p.flank_rows, p.flank_rows_max)
    withheld[donors[~settled]] = True
    ratio_c = np.where(settled, ratio_c, 1.0)
    ratio_m = np.where(settled, ratio_m, 1.0)
    a_d = areas[donors]
    share = csr_matrix(reach.multiply(areas[None, :]))
    row_total = np.asarray(share.sum(axis=1)).ravel()
    share = csr_matrix(share.multiply((1.0 / np.where(row_total > 0.0, row_total, 1.0))[:, None]))
    received_area = np.maximum(areas, 1.0)

    def move(name: str, ratio_vec: np.ndarray, give_back: bool = True) -> None:
        values = out[name]
        lost = values[donors] * (1.0 - ratio_vec) * a_d  # volume
        values[donors] = values[donors] * ratio_vec
        if give_back:
            values += (share.T @ lost) / received_area

    move("hc", ratio_c)
    move("hm", ratio_m)
    for name in ("material", "cover", "cover_material", "plateau", "arc", "moment"):
        move(name, ratio_c)
    craton_lost = float(np.dot(out["craton"][donors] * (1.0 - ratio_c), a_d))
    out["craton"][donors] *= ratio_c
    return ThinningOutcome(**out, craton_lost_m3=craton_lost, withheld=withheld)
