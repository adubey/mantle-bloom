"""Load, strength and thinning kernels of RFC 0018's rift lines (#846). Pure numpy/scipy on arrays: no `World`.

`rift_lines.py` places a plate's rift line and integrates its strain with these; `rift_damage.py` runs both on a
plate each step. Everything here is behind `World.split_method == "rift_lines"`.

- **Parameters:** `RiftLineParams`, the method's knobs with RFC 0018's starting values. A world overrides them
  through `World.rift_line_params`.
- **Load:** proxy C of RFC 0018 D0 (#797), `cut_tensor`: the membrane tensor fitted to the pull across straight line
  cuts in six directions, from the torque engine's per-node forces with their rigid-rotation part removed.
- **Strength:** `capacity` interpolates the column's extensional capacity `N_cap` from a table of the RFC 0010 Gate 2
  yield envelope (`data/rift_capacity_table.npz`, built by `bin/debug/build_rift_capacity_table.py`);
  `inherited_multiplier` and `softening` are the weakness marks and the strain softening that scale it.
- **Thinning:** `flank_plan` and `thin_with_flank_transfer` thin the straining band and move the removed volume to
  the band's flanks, conserving every column total.
"""

from dataclasses import dataclass, fields
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.spatial import cKDTree

from . import lithosphere, rheology, thin_shell

RADIUS_M = lithosphere.PLANET_RADIUS_M
YEAR_PER_MYR = 1.0e6
# Inherited-weakness multipliers on N_cap (RFC 0018 step 2). They began as uncalibrated route-cost weights
# (RFC 0010's "belt and plume resistance 0.6"); `mark_power` deepens them and RFC 0021 tunes them.
MULTIPLIER_FLOOR = 0.12
MULT_BELT, MULT_VOLCANO, MULT_PLUME, MULT_FAULT = 0.65, 0.65, 0.65, 0.60
KINDS = ("continental", "oceanic_young_ridge", "oceanic_mature")
_TABLE_PATH = Path(__file__).resolve().parent / "data" / "rift_capacity_table.npz"
_table: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None = None


@dataclass(frozen=True)
class RiftLineParams:
    """The rift-line method's knobs, at RFC 0018's starting values. The six the RFC searches are `eps_dot0`, `s_warm`,
    `k_cap`, `n_plume`, `lam_grow` and `l_tip_m`; the strength-map knobs at the end came out of #888."""

    eps_dot0: float = 8e-8  # reference strain rate (1/yr) of a line at twice its growth threshold
    n: float = 1.0  # exponent of the rate law; fixed by the RFC
    eps_soft: float = 1.0  # strain over which softening saturates; fixed by the RFC
    s_fault: float = 0.40  # permanent softening of fully strained rock
    s_warm: float = 0.40  # softening of recently strained rock, which heals
    s_max: float = 0.88  # cap on the two softenings together
    tau_heal_yr: float = 62.8e6  # e-folding time of the warm strain
    k_cap: float = 1.0  # scale on every column's capacity
    k_cut: float = 1.33  # scale on the line-cut load (D0 calibration)
    n_plume: float = 0.5e12  # isotropic tension over plume heat (N/m)
    w0_m: float = 125e3  # width of a new rift band, and of the softening average
    lam_grow: float = 0.3  # pull over strength above which a line strains
    l_tip_m: float = 250e3  # distance over which strain rate falls off ahead of a rift's tips
    r_max: float = 10.0  # cap on pull over strength relative to `lam_grow`
    a_min: float = 0.15  # smallest daughter, as a share of the parent's area
    a_stray: float = 0.01  # largest stray piece a split may hand to a daughter
    c_dir: float = 0.8  # a direction is a candidate when its peak pull is this share of the best
    c_corr: float = 0.8  # a corridor is the run of sections with this share of the peak pull
    f_break: float = 0.5  # length share of the line that must be broken before it can split the plate
    eps_cand: float = 0.1  # strain at which a cell joins the rift
    eps_break: float = 2.0  # strain at which a cell breaks whatever its thickness
    flank_rows: int = 2  # hops a band cell's removed volume first looks for recipients
    k_fail: int = 4096  # sub-steps per step beyond which the integration is treated as a defect
    d_eps_max: float = 0.02  # largest strain increment of one sub-step
    oceanic_latch_on_break_only: bool = False  # RFC 0018's oceanic fallback: break only once both layers reach their floors
    oceanic_k_cap_factor: float = 1.0  # extra capacity scale on oceanic cells
    # Strength of ocean floor from its age (#888). Oceanic mantle lithosphere relaxes to 40 km whatever its age, which
    # makes old ocean floor the weakest rock on a plate, half as strong as marked continental crust. With this on,
    # N_cap of an oceanic cell uses a cooling plate's thickness instead. The plate's own thickness, and so its
    # isostasy, is untouched.
    oceanic_age_strength: bool = True
    ocean_cooling_m_per_sqrt_myr: float = 13_000.0  # 2.32 sqrt(kappa t), kappa = 1e-6 m2/s (half-space cooling); starting value
    ocean_plate_thickness_m: float = 95_000.0  # where old ocean floor stops thickening; starting value
    # Every continental cell has the same 100 km of mantle lithosphere, so a craton is no stronger than its
    # surroundings (#888). This adds a root to the thickness behind N_cap, in proportion to craton strength. 0 = off.
    craton_root_m: float = 0.0
    mark_power: float = 1.0  # the weakness multipliers raised to this power (2.8 takes 0.65 to 0.30)
    route_mark_power: float = 1.0  # the same, for the line's path cost only: the route follows marks more closely

    @classmethod
    def with_overrides(cls, overrides: dict | None) -> "RiftLineParams":
        """The starting values with `overrides` applied. An unknown name or a value of the wrong kind is an error, so
        a mistyped knob cannot silently leave the default in place."""
        # The kind comes from the default's type, not `f.type`, which is a string under postponed annotations.
        known = {f.name: type(f.default) for f in fields(cls)}
        values = {}
        for name, value in (overrides or {}).items():
            if name not in known:
                raise ValueError(f"unknown rift-line parameter {name!r}; known: {sorted(known)}")
            values[name] = _checked(name, known[name], value)
        return cls(**values)


# Knobs that divide, set a scale or bound a step: zero crashes or turns N_eff into NaN downstream (PR #950 review).
_POSITIVE = frozenset({"eps_dot0", "n", "eps_soft", "tau_heal_yr", "k_cap", "k_cut", "w0_m", "lam_grow", "l_tip_m", "r_max",
                       "eps_cand", "eps_break", "k_fail", "d_eps_max", "ocean_cooling_m_per_sqrt_myr", "ocean_plate_thickness_m"})
# Shares and softenings: a share above 1 is meaningless, and a softening of 1 removes all strength.
_AT_MOST_ONE = frozenset({"a_min", "a_stray", "c_dir", "c_corr", "f_break", "s_fault", "s_warm"})


def _checked(name: str, kind: type, value):
    """`value` as a `kind` knob, or ValueError. Integer knobs refuse fractions rather than truncating them, so a
    sweep cannot silently run a different setting from the one it records."""
    if kind is bool:
        if not isinstance(value, bool):
            raise ValueError(f"rift-line parameter {name!r} must be true or false, got {value!r}")
        return value
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value) or value < 0:
        raise ValueError(f"rift-line parameter {name!r} must be a finite number >= 0, got {value!r}")
    if kind is int and float(value) != int(value):
        raise ValueError(f"rift-line parameter {name!r} must be a whole number, got {value!r}")
    if name in _POSITIVE and value <= 0:
        raise ValueError(f"rift-line parameter {name!r} must be > 0, got {value!r}")
    if (name in _AT_MOST_ONE and value > 1) or (name == "s_max" and value >= 1):
        raise ValueError(f"rift-line parameter {name!r} must be at most 1 (s_max below 1), got {value!r}")
    return int(value) if kind is int else float(value)


# --- load ----------------------------------------------------------------------------------------

CUT_ANGLES = np.linspace(0.0, np.pi, 6, endpoint=False)
_CUT_DESIGN = np.column_stack([np.cos(CUT_ANGLES) ** 2, np.sin(CUT_ANGLES) ** 2, 2 * np.cos(CUT_ANGLES) * np.sin(CUT_ANGLES)])
_CUT_INVERSE = np.linalg.pinv(_CUT_DESIGN)


def balanced_tangent_force(points, areas, force_n) -> np.ndarray:
    """Tangent force with its area-weighted rigid-rotation part removed, as `thin_shell.solve_thin_shell` does. The
    remainder has zero net torque, so a line cut through the plate sees only the part that stresses it."""
    tangent = force_n - np.sum(force_n * points, axis=1)[:, None] * points
    xyz = points * RADIUS_M
    modes = np.stack([np.cross(axis, xyz) for axis in np.eye(3)], axis=2)  # (N, 3 xyz, 3 rotation axes)
    root = np.sqrt(np.maximum(areas, 1.0))
    design = (modes * root[:, None, None]).reshape(-1, 3)
    target = (tangent / root[:, None]).ravel()
    coefficients, *_ = np.linalg.lstsq(design, target, rcond=None)
    return tangent - areas[:, None] * (modes @ coefficients)


def plane_frame(points: np.ndarray):
    """Tangent frame at the plate's centroid and the cells' plane coordinates in metres."""
    centre = points.mean(axis=0)
    centre /= np.linalg.norm(centre)
    e_u, e_v = (axis[0] for axis in thin_shell._tangent_basis(centre[None, :]))
    return e_u, e_v, RADIUS_M * np.column_stack([points @ e_u, points @ e_v])


def _line_resultants(coords, force_2d, areas, spacing_m, angle) -> np.ndarray:
    """Pull across the straight line through each cell normal to `angle`: the balanced force on its far side along
    the normal, divided by the plate's width there (cell area within half a spacing of the line, over the spacing)."""
    normal = np.array([np.cos(angle), np.sin(angle)])
    position, push = coords @ normal, force_2d @ normal
    order = np.argsort(position, kind="stable")
    ranked = position[order]
    beyond = np.empty(len(position))
    beyond[order] = np.cumsum(push[order][::-1])[::-1] - 0.5 * push[order]  # strictly beyond, plus half its own cell
    cumulative = np.concatenate([[0.0], np.cumsum(areas[order])])
    upper = np.searchsorted(ranked, position + 0.5 * spacing_m, side="right")
    lower = np.searchsorted(ranked, position - 0.5 * spacing_m, side="left")
    width = np.maximum((cumulative[upper] - cumulative[lower]) / spacing_m, spacing_m)
    return beyond / width


def cut_tensor(points, areas, balanced, spacing_m) -> np.ndarray:
    """Proxy C: the membrane tensor (xx, yy, xy in each cell's `thin_shell._tangent_basis`, N/m) fitted to the
    line-cut resultants of the balanced force in six directions."""
    e_u, e_v, coords = plane_frame(points)
    force_2d = np.column_stack([balanced @ e_u, balanced @ e_v])
    resultants = np.column_stack([_line_resultants(coords, force_2d, areas, spacing_m, a) for a in CUT_ANGLES])
    frame = resultants @ _CUT_INVERSE.T  # (N, 3): xx, yy, xy in the (e_u, e_v) frame
    first, second = thin_shell._tangent_basis(points)
    r11, r12, r21, r22 = first @ e_u, first @ e_v, second @ e_u, second @ e_v  # rotation from frame to cell basis
    xx, yy, xy = frame[:, 0], frame[:, 1], frame[:, 2]
    return np.column_stack([
        r11 * r11 * xx + 2 * r11 * r12 * xy + r12 * r12 * yy,
        r21 * r21 * xx + 2 * r21 * r22 * xy + r22 * r22 * yy,
        r11 * r21 * xx + (r11 * r22 + r12 * r21) * xy + r12 * r22 * yy,
    ])


def normal_component(tensor: np.ndarray, normal_ab: np.ndarray) -> np.ndarray:
    """n̂ᵀ𝐍n̂ for in-plane unit normals given as (a, b) in the tangent basis."""
    a, b = normal_ab[:, 0], normal_ab[:, 1]
    return tensor[:, 0] * a * a + tensor[:, 1] * b * b + 2.0 * tensor[:, 2] * a * b


# --- strength ------------------------------------------------------------------------------------


def capacity_table():
    """(table, Hc grid, Hm grid, Moho-lag grid): N_cap (N/m) by crust kind, crust thickness, mantle-lithosphere
    thickness and Moho thermal lag, from the RFC 0010 Gate 2 yield envelope. Loaded once."""
    global _table
    if _table is None:
        data = np.load(_TABLE_PATH)
        _table = (data["table"], data["hc"], data["hm"], data["lag"])
    return _table


def capacity(kind_index: np.ndarray, hc, hm, lag) -> np.ndarray:
    """Multilinear interpolation of `capacity_table` (clamped to the grid) for each cell."""
    table, *grids = capacity_table()
    out = np.empty(len(hc))
    for ki in range(len(KINDS)):
        sel = kind_index == ki
        if not np.any(sel):
            continue
        axes = []
        for value, grid in zip((hc[sel], hm[sel], lag[sel]), grids):
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


def oceanic_strength_hm(hc, hm, age_yr, eps, p: RiftLineParams) -> np.ndarray:
    """Mantle-lithosphere thickness (m) that sets an oceanic cell's strength: a cooling plate of the cell's age, less
    its crust, thinned by the rift strain applied to it, and never below the model's own thickness. A cell without a
    known age (`age_yr` NaN) keeps the model's thickness."""
    total = np.minimum(p.ocean_cooling_m_per_sqrt_myr * np.sqrt(np.maximum(age_yr, 0.0) / YEAR_PER_MYR), p.ocean_plate_thickness_m)
    aged = np.maximum(total - hc, 0.0) * np.exp(-eps)
    return np.where(np.isnan(age_yr), hm, np.maximum(hm, aged))


def softening(eps_bar, warm_bar, p: RiftLineParams) -> np.ndarray:
    """W = 1 - min(s_max, s_fault φ(ε̄) + s_warm φ(ε̄_warm)), φ(ε) = min(1, ε / ε_soft); in [1 - s_max, 1]."""
    phi = lambda e: np.minimum(1.0, np.maximum(e, 0.0) / p.eps_soft)  # noqa: E731
    return 1.0 - np.minimum(p.s_max, p.s_fault * phi(eps_bar) + p.s_warm * phi(warm_bar))


def inherited_multiplier(belt, volcano, plume, fault, mark_power: float = 1.0) -> np.ndarray:
    """Product of the inherited-weakness multipliers, raised to `mark_power` and clamped to [0.12, 1]."""
    m = np.ones(len(belt))
    for mask, factor in ((belt, MULT_BELT), (volcano, MULT_VOLCANO), (plume, MULT_PLUME), (fault, MULT_FAULT)):
        m = np.where(mask, m * factor, m)
    return np.clip(m ** mark_power, MULTIPLIER_FLOOR, 1.0)


def disk_average_matrix(points, areas, radius_m) -> csr_matrix:
    """Row-stochastic area-weighted average over the cells within `radius_m` (inclusive of itself). Softening is
    averaged over it so a band cannot localize below the disk's width, whatever the cell spacing."""
    chord = 2.0 * np.sin(0.5 * radius_m / RADIUS_M)
    neighbours = cKDTree(points).query_ball_point(points, chord * (1.0 + 1e-6))
    rows = np.repeat(np.arange(len(points)), [len(v) for v in neighbours])
    cols = np.concatenate([np.asarray(v, dtype=np.int64) for v in neighbours])
    weights = areas[cols]
    totals = np.bincount(rows, weights=weights, minlength=len(points))
    return csr_matrix((weights / totals[rows], (rows, cols)), shape=(len(points), len(points)))


# --- thinning with flank transfer ----------------------------------------------------------------


@dataclass
class FlankPlan:
    """Where each band cell's removed volume goes. `reach` is (band cells x plate cells) with 1 at the recipients of
    that donor; `held` marks band cells with no recipient anywhere on the plate (they must not strain); `hops` is the
    reach each donor needed (0 for held cells)."""

    band: np.ndarray  # plate indices of the band cells
    reach: csr_matrix
    held: np.ndarray  # per band cell
    hops: np.ndarray  # per band cell


def flank_plan(adjacency: csr_matrix, band: np.ndarray, broken: np.ndarray, rows: int, group: np.ndarray | None = None) -> FlankPlan:
    """Recipients are the unbroken plate cells outside the band, of the donor's `group` when one is given (its crust
    type: continental crust moved onto ocean floor would be priced at oceanic density, PR #950 review). Each donor takes
    those within `rows` hops, and widens one hop at a time, with no upper limit, until it has some. A donor in a
    component with none is held."""
    n = adjacency.shape[0]
    donors = np.flatnonzero(band)
    recipient_ok = ~band & ~broken
    if not len(donors):
        return FlankPlan(donors, csr_matrix((0, n), dtype=np.int8), np.zeros(0, bool), np.zeros(0, int))
    pattern = (adjacency > 0).astype(np.int8)
    if group is None:
        reach, settled, hops = _flank_reach(pattern, donors, recipient_ok, rows)
        return FlankPlan(donors, reach, ~settled, hops)
    reach, settled, hops = csr_matrix((len(donors), n), dtype=np.int8), np.zeros(len(donors), bool), np.zeros(len(donors), int)
    for value in np.unique(group[donors]):
        rows_of = np.flatnonzero(group[donors] == value)
        r, s, h = _flank_reach(pattern, donors[rows_of], recipient_ok & (group == value), rows)
        place = csr_matrix((np.ones(len(rows_of), dtype=np.int8), (rows_of, np.arange(len(rows_of)))), shape=(len(donors), len(rows_of)))
        reach = reach + place @ r
        settled[rows_of], hops[rows_of] = s, h
    reach.data[:] = 1
    return FlankPlan(donors, reach.tocsr(), ~settled, hops)


def _flank_reach(pattern: csr_matrix, donors: np.ndarray, recipient_ok: np.ndarray, rows: int):
    """(reach, settled, hops) of `flank_plan` for one set of donors and the cells they may give to."""
    n = pattern.shape[0]
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
    return result.tocsr(), settled, hops


# The thickness-like columns thinning moves, in the order `ThinningOutcome` holds them.
COLUMNS = ("hc", "hm", "material", "craton", "cover", "cover_material", "plateau", "arc", "moment")


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
    cap_held_m3: float = 0.0  # crust and mantle a donor kept because its recipients were at the thickness cap


def thin_with_flank_transfer(areas, plan: FlankPlan, d_eps, columns: dict) -> ThinningOutcome:
    """Pure-shear thinning of the straining band cells by exp(-Δε), floored, with the removed volume moved to the
    plan's flank recipients (same plate, outside the band, unbroken) in proportion to their area. `columns` holds
    `COLUMNS` per cell, in thickness units. Volumes come from the post-clamp thicknesses, so a floored donor gives up
    only what it lost. Craton crust thins in place, as `_open_rift` does, and its lost volume is returned for the
    `rifted_stretch_m3` sink."""
    out = {k: v.astype(float).copy() for k, v in columns.items()}
    hc, hm = out["hc"], out["hm"]
    take = (d_eps[plan.band] > 0.0) & ~plan.held
    donors = plan.band[take]
    if not len(donors):
        return ThinningOutcome(**out, craton_lost_m3=0.0)
    ratio = np.exp(-d_eps[donors])
    floor_c, floor_m = rheology.RIFT_CRITICAL_THICKNESS_M, lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M
    ratio_c = np.where(hc[donors] > floor_c, np.maximum(ratio, floor_c / np.maximum(hc[donors], 1.0)), 1.0)
    ratio_m = np.where(hm[donors] > floor_m, np.maximum(ratio, floor_m / np.maximum(hm[donors], 1.0)), 1.0)
    a_d = areas[donors]
    share = csr_matrix(plan.reach[take].multiply(areas[None, :]))
    row_total = np.asarray(share.sum(axis=1)).ravel()
    share = csr_matrix(share.multiply((1.0 / np.where(row_total > 0.0, row_total, 1.0))[:, None]))
    received_area = np.maximum(areas, 1.0)
    lost_c, lost_m = hc[donors] * (1.0 - ratio_c) * a_d, hm[donors] * (1.0 - ratio_m) * a_d
    flow, leaves = _capped_flow(share, areas, ((hc, lost_c, lithosphere.MAX_CRUSTAL_THICKNESS_M),
                                               (hm, lost_m, lithosphere.MAX_MANTLE_LITHOSPHERE_THICKNESS_M)))
    cap_held = float(np.dot(lost_c + lost_m, 1.0 - leaves))
    # A donor keeps the share of its loss its recipients refused, the same share of every layer, so a capped donor
    # thins less but keeps its column's shape (thinning the dense mantle alone would lift it).
    ratio_c, ratio_m = 1.0 - (1.0 - ratio_c) * leaves, 1.0 - (1.0 - ratio_m) * leaves

    def move(name: str, ratio_vec: np.ndarray) -> None:
        values = out[name]
        lost = values[donors] * (1.0 - ratio_vec) * a_d  # volume
        values[donors] = values[donors] * ratio_vec
        values += (flow.T @ lost) / received_area

    move("hc", ratio_c)
    move("hm", ratio_m)
    for name in ("material", "cover", "cover_material", "plateau", "arc", "moment"):
        move(name, ratio_c)
    craton_lost = float(np.dot(out["craton"][donors] * (1.0 - ratio_c), a_d))
    out["craton"][donors] *= ratio_c
    return ThinningOutcome(**out, craton_lost_m3=craton_lost, cap_held_m3=cap_held)


def _capped_flow(share: csr_matrix, areas, layers):
    """Cap each recipient at each layer's ceiling, as every other path that adds crust does (orogeny, faults, merges;
    PR #950 review). `layers` is (thickness, each donor's loss in m³, ceiling) per layer. A recipient accepts the share
    of its incoming volume that fits under the tightest of its layers; each donor sends only what its recipients
    accept. Returns (flow with each donor's row summing to 1, the share of each donor's loss that leaves it)."""
    fits = np.ones(share.shape[1])
    for values, lost, cap in layers:
        incoming = share.T @ lost
        room = np.maximum(cap - values, 0.0) * areas
        over = incoming > room
        fits[over] = np.minimum(fits[over], room[over] / incoming[over])
    if np.all(fits >= 1.0):
        return share, np.ones(share.shape[0])
    accepted = csr_matrix(share.multiply(fits[None, :]))
    leaves = np.asarray(accepted.sum(axis=1)).ravel()
    flow = csr_matrix(accepted.multiply((1.0 / np.where(leaves > 0.0, leaves, 1.0))[:, None]))
    return flow, leaves
