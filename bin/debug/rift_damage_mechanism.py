"""Offline prototype of RFC 0018's damage-band mechanism (D1, #888). Not production code: #846 implements it.

Usage: imported by ``check_rift_damage_prototype.py`` (unit-level checks) and ``run_rift_damage_d1.py`` (the
D1 driver). ``install(params)`` registers the three damage surface fields in this process, wraps
``merge_split.apply_topology_changes`` via ``maybe_split_plate`` and the margin-insertion path, and returns
a ``Mechanism`` whose ``report`` holds what the recorder needs. ``uninstall()`` restores the originals.

The mechanism follows the RFC's Design (2026-11-15 rift-line revision), steps 1-6:

* load: proxy C, ``rift_damage_prototype.cut_tensor`` scaled by ``k_cut`` (plume load added until expiry), and the
  section pull across straight sections (``rift_line.section_profile``);
* strength: ``k_cap N_cap M W`` with non-local softening over a disk of diameter ``w_0``;
* rift line and strain (``rift_line``): one line per plate placed from the section pull, growing from a nucleus and
  its tips, integrated in sub-steps with exact exponential healing;
* thinning: flank transfer of the removed volume (``_allocate_stretch``'s rule, run in reverse), unlimited reach;
* the broken latch, with no interior eruption;
* breakthrough: the integrated ligament rule on the line, then a mask-based split.

Readings of the RFC the prototype had to choose (each is marked in the code and listed in
``analysis/issue888/README.md`` for review):

* ``N_cap`` is frozen at the step start (the RFC's operator split); thinning enters it at the next step.
* margin insertion inherits the area-weighted mean of the donors within ``K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION``
  hops times the stretch share, without ``_allocate_stretch``'s craton willingness weights.
* eruption resets are detected after the fact: a cell whose column became exactly the reference column
  and which was not before.
* the plume's expiry time differs per cell, so the sub-step grid includes every expiry in the step and each
  cell switches its plume load and weakening off at its own expiry.
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
    """RFC 0018 starting values (2026-11-15 rift-line revision). ``n`` and ``eps_soft`` are fixed by the RFC."""

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
    lam_grow: float = 0.3
    l_tip_m: float = 250e3
    r_max: float = 10.0
    a_min: float = 0.15
    a_stray: float = 0.01
    c_dir: float = 0.8
    c_corr: float = 0.8
    f_break: float = 0.5
    eps_cand: float = 0.1
    eps_break: float = 2.0
    t_fail_yr: float = 5e6
    flank_rows: int = 2
    k_fail: int = 4096
    d_eps_max: float = 0.02
    live_ncap: bool = True  # the line cells' N_cap follows their thinning inside the sub-steps (RFC revision 2026-11-15)
    crossing_stops: bool = True  # a sub-step ends when a line cell reaches eps_cand (RFC revision 2026-11-15)
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
HM_GRID = np.concatenate([np.arange(2_000.0, 24_000.0, 2_000.0), np.arange(24_000.0, 200_001.0, 8_000.0)])  # fine where thin
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
        if all(np.array_equal(data[k], grid) for k, grid in (("hc", HC_GRID), ("hm", HM_GRID), ("lag", LAG_GRID))):
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


@dataclass
class FlankPlan:
    """Where each band cell's removed volume goes. ``reach`` is (band cells x plate cells) with 1 at the recipients
    of that donor; ``held`` marks band cells with no recipient anywhere on the plate (they must not strain);
    ``hops`` is the reach each donor needed (0 for held cells)."""

    band: np.ndarray  # plate indices of the band cells
    reach: csr_matrix
    held: np.ndarray  # per band cell
    hops: np.ndarray  # per band cell


def flank_plan(adjacency: csr_matrix, band: np.ndarray, broken: np.ndarray, rows: int) -> FlankPlan:
    """Recipients are the unbroken plate cells outside the band. Each donor takes those within ``rows`` hops, and
    widens one hop at a time, with no upper limit, until it has some. A donor in a component with none is held."""
    n = adjacency.shape[0]
    donors = np.flatnonzero(band)
    recipient_ok = ~band & ~broken
    if not len(donors):
        return FlankPlan(donors, csr_matrix((0, n), dtype=np.int8), np.zeros(0, bool), np.zeros(0, int))
    pattern = (adjacency > 0).astype(np.int8)
    visited = csr_matrix((np.ones(len(donors), dtype=np.int8), (np.arange(len(donors)), donors)), shape=(len(donors), n))
    frontier = visited.copy()
    result = csr_matrix((len(donors), n), dtype=np.int8)
    ok = csr_matrix(np.where(recipient_ok, 1, 0).astype(np.int8)[None, :])
    settled = np.zeros(len(donors), dtype=bool)
    hops = np.zeros(len(donors), dtype=int)
    hop = 0
    while True:
        hop += 1
        frontier = frontier @ pattern
        frontier.data[:] = 1
        grown = visited + frontier
        grown.data[:] = 1
        stalled = grown.nnz == visited.nnz
        visited = grown
        if hop >= rows:
            reach = visited.multiply(ok).tocsr()
            take = (np.asarray(reach.sum(axis=1)).ravel() > 0) & ~settled
            if take.any():
                result = result + csr_matrix(reach.multiply(take[:, None].astype(np.int8)))
                settled |= take
                hops[take] = hop
        if settled.all() or stalled:
            break
    result.data[:] = 1
    return FlankPlan(donors, result.tocsr(), ~settled, hops)


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


def thin_with_flank_transfer(areas, plan: FlankPlan, d_eps, columns: dict) -> ThinningOutcome:
    """Pure-shear thinning of the straining band cells by exp(-Δε), floored, with the removed volume moved to the
    plan's flank recipients (same plate, outside the band, unbroken) in proportion to their area. ``columns`` holds
    hc, hm, material, craton, cover, cover_material, plateau, arc, moment (per cell, thickness units). Volumes are
    computed from the post-clamp thicknesses, so a floored donor gives up only what it lost. Craton crust thins in
    place and its lost volume is returned for the ``rifted_stretch_m3`` sink."""
    out = {k: v.astype(float).copy() for k, v in columns.items()}
    hc, hm = out["hc"], out["hm"]
    take = (d_eps[plan.band] > 0.0) & ~plan.held
    donors = plan.band[take]
    if not len(donors):
        return ThinningOutcome(**out, craton_lost_m3=0.0)
    ratio = np.exp(-d_eps[donors])
    ratio_c = np.where(hc[donors] > rheology.RIFT_CRITICAL_THICKNESS_M,
                       np.maximum(ratio, rheology.RIFT_CRITICAL_THICKNESS_M / np.maximum(hc[donors], 1.0)), 1.0)
    ratio_m = np.where(hm[donors] > lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M,
                       np.maximum(ratio, lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M / np.maximum(hm[donors], 1.0)), 1.0)
    a_d = areas[donors]
    share = csr_matrix(plan.reach[take].multiply(areas[None, :]))
    row_total = np.asarray(share.sum(axis=1)).ravel()
    share = csr_matrix(share.multiply((1.0 / np.where(row_total > 0.0, row_total, 1.0))[:, None]))
    received_area = np.maximum(areas, 1.0)

    def move(name: str, ratio_vec: np.ndarray) -> None:
        values = out[name]
        lost = values[donors] * (1.0 - ratio_vec) * a_d  # volume
        values[donors] = values[donors] * ratio_vec
        values += (share.T @ lost) / received_area

    move("hc", ratio_c)
    move("hm", ratio_m)
    for name in ("material", "cover", "cover_material", "plateau", "arc", "moment"):
        move(name, ratio_c)
    craton_lost = float(np.dot(out["craton"][donors] * (1.0 - ratio_c), a_d))
    out["craton"][donors] *= ratio_c
    return ThinningOutcome(**out, craton_lost_m3=craton_lost)
