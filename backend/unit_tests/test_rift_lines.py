from dataclasses import replace

import numpy as np
import pytest
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from app import rheology, thin_shell
from app import rift_lines as lines
from app import rift_strength as strength
from app.rift_strength import RADIUS_M, RiftLineParams

from .test_rift_strength import grid_edges, lattice

P = RiftLineParams()
MYR = 1e6
NUCLEUS_FACTOR = 0.65  # extra weakness of the middle cell, so the nucleus is there


class Toy:
    """A flat lattice on a sphere patch with one straight line along its middle row, under a uniform pull across it."""

    def __init__(self, p: RiftLineParams = P, pull_ratio: float = 2.0, nx: int = 31, ny: int = 9, h_m: float = 125e3):
        self.p = p
        u, v = (np.arange(nx) - (nx - 1) / 2) * h_m, (np.arange(ny) - (ny - 1) / 2) * h_m
        uu, vv = np.meshgrid(u, v)
        lon, lat = uu.ravel() / RADIUS_M, vv.ravel() / RADIUS_M
        self.points = np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
        n = len(self.points)
        r, c = grid_edges(nx, ny)
        length = RADIUS_M * np.arccos(np.clip(np.sum(self.points[r] * self.points[c], axis=1), -1.0, 1.0))
        self.graph = csr_matrix((length, (r, c)), shape=(n, n))
        self.areas = np.full(n, h_m * h_m)
        self.chain = np.arange(n).reshape(ny, nx)[ny // 2]
        self.disk = strength.disk_average_matrix(self.points, self.areas, p.w0_m / 2.0)
        self.hc, self.hm, self.lag = np.full(n, 35e3), np.full(n, 90e3), np.zeros(n)
        self.mult = np.full(n, strength.MULT_BELT)
        self.mult[self.chain[nx // 2]] *= NUCLEUS_FACTOR
        self.eps, self.warm, self.broken = np.zeros(n), np.zeros(n), np.zeros(n, dtype=bool)
        self.geometry = lines.path_geometry(self.points, self.chain, *thin_shell._tangent_basis(self.points))
        # The path-normal load that gives the requested initial Λ on the line (no softening yet).
        self.sigma = pull_ratio * p.lam_grow * float(np.mean(self.base()[self.chain]))

    def base(self) -> np.ndarray:
        return strength.capacity(np.zeros(len(self.hc), dtype=int), self.hc, self.hm, self.lag) * self.p.k_cap * self.mult

    def problem(self) -> lines.LineProblem:
        p, chain = self.p, self.chain
        band, reach_band, reach_line = lines.band_reach(self.points, chain, self.eps, p)
        floor = rheology.RIFT_CRITICAL_THICKNESS_M
        latch = np.minimum(p.eps_break, self.eps + np.log(np.maximum(self.hc, floor) / floor))
        latch = np.where(self.broken, self.eps, np.maximum(latch, self.eps))
        base = self.base()
        return lines.LineProblem(
            chain=chain, s=self.geometry.s, ds=self.geometry.ds, normal_load=np.full(len(chain), self.sigma),
            base_quiet=base[chain], base_plume=base[chain], plume_left_yr=np.zeros(len(chain)), intact=~self.broken[chain],
            disk_rows=self.disk[chain], band=band, reach_band=reach_band, reach_line=reach_line, latch=latch,
            capacity=lines.CapacityInputs(np.zeros(len(chain), dtype=int), self.hc[chain], self.hm[chain], self.lag[chain],
                                          np.full(len(chain), p.k_cap), self.mult[chain], self.mult[chain], self.eps[chain]))


def _connected(graph, mask) -> bool:
    return connected_components(graph[mask][:, mask], directed=False)[0] == 1


def test_section_load_is_always_finite_and_an_open_line_outranks_any_other():
    assert lines.section_load(0.0, 1.0, P) == (0.0, False) and lines.section_load(-1.0, 0.0, P) == (0.0, False)
    lam, is_open = lines.section_load(1.0, 0.0, P)
    assert is_open and lam == P.r_max * P.lam_grow
    assert lines.growth_ratio(lam, True, P) == P.r_max and lines.growth_ratio(10.0, False, P) == P.r_max
    assert lines.rank_key(1e9, 0.1, True) > lines.rank_key(1.0, 99.0, False) and lines.rank_key(2.0, 0, True) > lines.rank_key(1.0, 0, True)
    assert lines.ligament_gives(1.0, 0.0) and lines.ligament_gives(2.0, 1.0)
    assert not lines.ligament_gives(0.0, 0.0) and not lines.ligament_gives(-1.0, 0.0) and not lines.ligament_gives(1.0, 2.0)


def test_only_load_across_the_line_counts():
    normal = np.tile([0.0, 1.0], (10, 1))  # line along x, normal along y
    along, across = np.tile([2e11, 0.0, 0.0], (10, 1)), np.tile([0.0, 2e11, 0.0], (10, 1))
    assert np.allclose(lines.path_normal_load(along, normal), 0.0)
    assert np.allclose(lines.path_normal_load(across, normal), 2e11)


def test_path_length_through_a_cell_is_half_the_distance_to_each_neighbour():
    pts = np.column_stack([np.ones(3), np.array([0.0, 0.01, 0.02]), np.zeros(3)])
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    geometry = lines.path_geometry(pts, np.arange(3), *thin_shell._tangent_basis(pts))
    spacing = 0.01 * RADIUS_M
    assert np.allclose(geometry.ds / spacing, [0.5, 1.0, 0.5], rtol=1e-3)
    assert geometry.length == pytest.approx(2 * spacing, rel=1e-3) and geometry.straight == pytest.approx(geometry.length, rel=1e-3)


def _sections_for(force_sign: float):
    toy = Toy(nx=25, ny=15)
    side = np.where(toy.points[:, 1] > 0, 1.0, -1.0)  # each half pulled away from (or pushed toward) the middle
    force = force_sign * 1e11 * side[:, None] * np.array([0.0, 1.0, 0.0])[None, :]
    balanced = strength.balanced_tangent_force(toy.points, toy.areas, force)
    return toy, lines.section_profile(toy.points, toy.areas, balanced, P.k_cut)


def test_section_pull_peaks_where_two_halves_are_pulled_apart_and_respects_the_daughter_minimum():
    toy, sections = _sections_for(+1.0)
    peaks = [float(np.max(sections.pull[:, d])) for d in range(lines.N_DIRECTIONS)]
    expected = P.k_cut * 1e11 * float(np.sum(toy.points[:, 1] > 0))  # the force on one half
    assert max(peaks) == pytest.approx(expected, rel=0.1)
    corridors = lines.candidate_corridors(sections, P)
    assert corridors
    for corridor in corridors:
        smaller = np.minimum(sections.far_share[:, corridor.direction], 1 - sections.far_share[:, corridor.direction])
        assert corridor.peak_share >= P.a_min - 1e-9 and np.all(smaller[corridor.cells] >= P.a_min - 1e-9)


def test_a_plate_pushed_together_has_no_real_pull():
    _, pushed = _sections_for(-1.0)
    valid = lambda d: np.minimum(pushed.far_share[:, d], 1 - pushed.far_share[:, d]) >= P.a_min  # noqa: E731
    push = max(abs(float(np.min(np.where(valid(d), pushed.pull[:, d], np.inf)))) for d in range(lines.N_DIRECTIONS))
    for corridor in lines.candidate_corridors(pushed, P):  # the six-direction fit leaks a few per cent sideways
        assert corridor.peak_pull < 0.05 * push


def _corridor_lattice(nx=24, ny=13):
    index = np.arange(nx * ny).reshape(ny, nx)
    boundary = np.zeros((ny, nx), bool)
    boundary[0], boundary[-1], boundary[:, 0], boundary[:, -1] = True, True, True, True
    strike = np.broadcast_to(np.arange(nx, dtype=float), (ny, nx)).copy()  # the line runs along x
    corridor = np.zeros((ny, nx), bool)
    corridor[3:10] = True
    return lattice(nx, ny), index, boundary.ravel(), strike.ravel(), corridor.ravel(), nx, ny


def test_the_line_crosses_the_corridor_goes_around_strong_rock_and_follows_weak_rock():
    graph, index, boundary, strike, corridor, nx, ny = _corridor_lattice()
    straight = lines.find_line(graph, corridor, boundary, strike, np.ones(nx * ny))
    assert len(straight) >= nx - 3 and boundary[straight[0]] and boundary[straight[-1]]
    assert np.all(graph[straight[:-1], straight[1:]] > 0) and len(set(straight.tolist())) == len(straight)
    craton = np.zeros((ny, nx), bool)
    craton[4:9, 8:15] = True
    around = lines.find_line(graph, corridor, boundary, strike, np.where(craton.ravel(), 50.0, 1.0))
    assert not craton.ravel()[around].any()
    belt = np.full((ny, nx), 1.0)
    belt[8, :] = 0.2
    along_belt = lines.find_line(graph, corridor, boundary, strike, belt.ravel())
    assert np.isin(along_belt, index[8]).mean() > 0.8


def test_a_line_that_holds_a_rift_is_found_again_and_a_cut_corridor_gives_none():
    graph, index, boundary, strike, corridor, nx, ny = _corridor_lattice()
    belt = np.full((ny, nx), 1.0)
    belt[8, :] = 0.2
    first = lines.find_line(graph, corridor, boundary, strike, belt.ravel())
    broken_cost = np.where(np.isin(np.arange(nx * ny), first), 0.0, 1.0)  # broken cells cost nothing; the belt is gone
    assert np.array_equal(lines.find_line(graph, corridor, boundary, strike, broken_cost), first)
    cut = corridor.copy()
    cut.reshape(ny, nx)[:, 11:13] = False
    assert lines.find_line(graph, cut, boundary, strike, np.ones(nx * ny)) is None


def test_strain_starts_at_the_weakest_cell_spreads_along_the_line_and_stays_in_the_band():
    toy = Toy(pull_ratio=5.0)
    prob = toy.problem()
    early = lines.advance_rift(prob, toy.eps, toy.warm, 0.2 * MYR, P)
    nucleus = len(toy.chain) // 2
    assert int(np.argmax(early.eps[toy.chain])) == nucleus
    late = lines.advance_rift(prob, toy.eps, toy.warm, 5 * MYR, P)
    assert (late.eps[toy.chain] >= P.eps_cand).sum() > (early.eps[toy.chain] >= P.eps_cand).sum()
    outside = np.ones(len(toy.eps), dtype=bool)
    outside[prob.band] = False
    assert not late.eps[outside].any()
    assert np.all(late.warm <= late.eps + 1e-12)


def test_a_cell_stops_at_its_latch_strain():
    toy = Toy(pull_ratio=5.0)
    prob = toy.problem()
    eps0 = toy.eps.copy()
    eps0[toy.chain[3]] = 1.0
    prob.latch[toy.chain[3]] = 1.0  # already at its latch: held or broken
    out = lines.advance_rift(prob, eps0, toy.warm, 5 * MYR, P)
    assert out.eps[toy.chain[3]] == 1.0 and out.eps[toy.chain].max() > 0.2
    assert np.all(out.eps <= prob.latch + 1e-12)


def test_a_line_below_its_growth_threshold_does_not_strain():
    toy = Toy(pull_ratio=0.9)
    out = lines.advance_rift(toy.problem(), toy.eps, toy.warm, 20 * MYR, P)
    assert not out.eps.any() and out.ratio_start == pytest.approx(0.9, rel=1e-2)


def test_an_open_line_strains_at_the_capped_rate_and_one_without_pull_does_not():
    toy = Toy()
    prob = toy.problem()
    prob.intact[:] = False
    out = lines.advance_rift(prob, toy.eps, toy.warm, 5 * MYR, P)
    assert out.open_start and np.isfinite(out.eps).all() and out.eps[toy.chain].max() > 0.0
    assert np.all(out.eps <= prob.latch + 1e-12)
    prob.normal_load[:] = 0.0
    idle = lines.advance_rift(prob, toy.eps, toy.warm, 5 * MYR, P)
    assert not idle.eps.any() and not idle.open_start


def test_warm_strain_heals_the_same_amount_however_the_time_is_stepped():
    toy = Toy()
    toy.eps[:], toy.warm[:] = 1.0, 0.8
    prob = toy.problem()
    prob.normal_load[:] = 0.0
    for pieces in (1, 3):
        eps, warm = toy.eps.copy(), toy.warm.copy()
        for _ in range(pieces):
            out = lines.advance_rift(prob, eps, warm, P.tau_heal_yr / pieces, P)
            eps, warm = out.eps, out.warm
        assert np.allclose(eps, 1.0) and np.allclose(warm, 0.8 * np.exp(-1.0), rtol=1e-6)


def test_a_step_across_a_plume_expiry_equals_its_two_parts():
    q = replace(P, s_fault=0.0, s_warm=0.0, n_plume=0.5e12)
    toy = Toy(q, pull_ratio=1.2)
    heated = toy.problem()
    heated.plume_left_yr = np.full(len(heated.chain), 0.5 * MYR)
    crossed = lines.advance_rift(heated, toy.eps, toy.warm, 1 * MYR, q)
    first = lines.advance_rift(heated, toy.eps, toy.warm, 0.5 * MYR, q)
    cold = toy.problem()
    second = lines.advance_rift(cold, first.eps, first.warm, 0.5 * MYR, q)
    assert np.allclose(crossed.eps, second.eps, rtol=1e-6, atol=1e-9)
    never = lines.advance_rift(cold, toy.eps, toy.warm, 1 * MYR, q)
    assert never.eps[toy.chain].max() < crossed.eps[toy.chain].max()


def test_substeps_are_bounded_and_running_out_of_them_is_an_error():
    toy = Toy(pull_ratio=P.r_max * 1.2)
    prob = toy.problem()
    out = lines.advance_rift(prob, toy.eps, toy.warm, 1 * MYR, P)
    bound = int(np.ceil(1 * MYR * P.eps_dot0 * (P.r_max - 1) ** P.n / P.d_eps_max)) + 2 * len(prob.chain)
    assert out.substeps <= bound
    with pytest.raises(RuntimeError, match="K_fail"):
        lines.advance_rift(prob, toy.eps, toy.warm, 60 * MYR, replace(P, k_fail=3))


def test_the_tip_advances_at_the_speed_of_the_discrete_law():
    h, nx = 125e3, 31
    q = replace(P, l_tip_m=125e3, d_eps_max=0.05)
    toy = Toy(q, nx=nx, ny=3, h_m=h)
    toy.mult[:] = strength.MULT_BELT
    toy.mult[toy.chain[0]] *= NUCLEUS_FACTOR  # one front, starting at the line's first cell
    amplitude, x = q.eps_dot0, h / q.l_tip_m
    law_speed = (amplitude * q.l_tip_m / q.eps_cand) * x / (np.exp(x) - 1.0)  # m per year
    cross = np.full(nx, np.nan)
    eps, warm = toy.eps, toy.warm
    for k in range(1, int(1.3 * nx * h / law_speed / MYR) + 10):
        toy.eps = eps
        out = lines.advance_rift(toy.problem(), eps, warm, MYR, q, fixed_amplitude=amplitude)
        eps, warm = out.eps, out.warm
        cross[np.isnan(cross) & (eps[toy.chain] >= q.eps_cand)] = k
    cells = np.arange(nx)
    steady = ~np.isnan(cross) & (cells >= 6) & (cells <= nx - 4)
    assert steady.sum() >= 10
    speed = np.polyfit(cross[steady] * MYR, cells[steady] * h, 1)[0]
    assert speed == pytest.approx(law_speed, rel=0.2)


def test_a_line_across_a_plate_gives_two_connected_daughters():
    nx, ny = 20, 10
    graph, areas = lattice(nx, ny), np.ones(nx * ny)
    index = np.arange(nx * ny).reshape(ny, nx)
    mask, why = lines.split_masks(graph, areas, index[:, 10], P)
    assert why is None and 85 <= mask.sum() <= 115
    assert _connected(graph, mask) and _connected(graph, ~mask)


def test_a_split_keeps_small_strays_refuses_large_ones_and_refuses_a_small_daughter():
    nx, ny = 20, 10
    graph, areas = lattice(nx, ny), np.ones(nx * ny)
    index = np.arange(nx * ny).reshape(ny, nx)
    # A 3x3 pocket (4.5% of the plate) walled off by a ring beside the line.
    ring = np.unique(np.concatenate([index[2, 12:17], index[6, 12:17], index[2:7, 12], index[2:7, 16]]))
    cut, pocket = np.concatenate([index[:, 10], ring]), index[3:6, 13:16].ravel()
    assert lines.split_masks(graph, areas, cut, replace(P, a_stray=0.001)) == (None, "stray_too_large")
    mask, why = lines.split_masks(graph, areas, cut, replace(P, a_stray=0.05))
    assert why is None and len({bool(mask[c]) for c in pocket}) == 1 and not mask[pocket[0]]
    assert lines.split_masks(graph, areas, index[:, 2], P) == (None, "daughter_small")
    assert lines.split_masks(graph, areas, index[:, 2], replace(P, a_min=0.05))[0] is not None


def test_an_island_joins_the_nearest_daughter_and_a_thick_line_still_splits():
    nx, ny = 20, 10
    n = nx * ny
    index = np.arange(n).reshape(ny, nx)
    coo = lattice(nx, ny).tocoo()
    padded = csr_matrix((coo.data, (coo.row, coo.col)), shape=(n + 1, n + 1))  # one extra cell with no edges
    points = np.column_stack([np.arange(n + 1) % nx, np.arange(n + 1) // nx, np.zeros(n + 1)]).astype(float)
    points[n] = [2.0, 5.0, 0.0]  # beside the left half
    mask, why = lines.split_masks(padded, np.ones(n + 1), index[:, 9:12].ravel(), P, points=points)
    assert why is None and mask[n] == mask[index[5, 2]]
    assert _connected(padded[:n, :n], mask[:n]) and _connected(padded[:n, :n], ~mask[:n])
