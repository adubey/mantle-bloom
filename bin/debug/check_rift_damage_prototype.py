"""Unit-level checks of the RFC 0018 rift-line prototype pieces (D1, #888). Production tests come with #846.

Usage (repo root, this worktree's venv; the plate-level checks generate one small world, about a minute):
    backend/.venv/bin/python bin/debug/check_rift_damage_prototype.py [--seed 106533960] [--only NAME ...] [--synthetic]

Each ``check_*`` function asserts one of the RFC's acceptance tests (the matching bullet is in its docstring) and
the script prints PASS or the failure for each. Exit status is non-zero if any fails. Synthetic checks need no
world (``--synthetic`` runs only those); the plate-level checks share one generated development world. The
integration and propagation checks are in ``check_rift_line_d1.py`` and the line check on the 20 worlds in
``check_rift_lines_world_d1.py``.
"""

import argparse
import sys
import traceback
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

sys.path[:0] = [str(Path(__file__).resolve().parent)]

import check_rift_line_d1 as toy_mod  # noqa: E402
import rift_damage_mechanism as mech  # noqa: E402
import rift_damage_prototype as proto  # noqa: E402
import rift_damage_step as step_mod  # noqa: E402
import rift_line as rl  # noqa: E402
from rift_damage_mechanism import Params  # noqa: E402

from app import lithosphere, rheology  # noqa: E402

P = Params()
_WORLD = {}
MYR = 1e6


def lattice(nx: int, ny: int):
    """4-connected grid graph with centre distances 1000 m."""
    r, c = toy_mod.grid_edges(nx, ny)
    return csr_matrix((np.full(len(r), 1000.0), (r, c)), shape=(nx * ny, nx * ny))


def columns_for(n, rng):
    return {"hc": rng.uniform(20e3, 45e3, n), "hm": rng.uniform(60e3, 120e3, n), "material": rng.uniform(0, 30e3, n),
            "craton": rng.uniform(0, 10e3, n), "cover": rng.uniform(0, 500, n), "cover_material": rng.uniform(0, 200, n),
            "plateau": rng.uniform(0, 3e3, n), "arc": rng.uniform(0, 2e3, n), "moment": rng.uniform(0, 1e9, n)}


def toy(p: Params = P, pull_ratio: float = 2.0, nx: int = 31, ny: int = 9, h: float = 125e3) -> toy_mod.Toy:
    return toy_mod.Toy(h, nx, ny, p, pull_ratio)


# --- synthetic checks ----------------------------------------------------------------------------


def check_strength_bounds():
    """Bounds of the strength factors: W in [1 - s_max, 1], M in [0.12, 1]."""
    e = np.linspace(0, 6, 61)
    for sf, sw in ((0.4, 0.4), (0.9, 0.9), (0.0, 0.0)):
        w = mech.softening(e[:, None], e[None, :], Params(s_fault=sf, s_warm=sw))
        assert w.min() >= 1 - P.s_max - 1e-12 and w.max() <= 1.0 + 1e-12, (w.min(), w.max())
    flags = np.array([[a, b, c, d] for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)], bool)
    m = mech.inherited_multiplier(*flags.T)
    assert m.min() >= 0.12 and m.max() <= 1.0 and m[0] == 1.0, m
    assert abs(m[-1] - max(0.12, 0.65 ** 3 * 0.6)) < 1e-12


def check_softening_monotone():
    """Softening monotonicity: W never rises with strain or warm strain."""
    e = np.linspace(0, 5, 200)
    for fixed in (0.0, 0.5, 2.0):
        assert np.all(np.diff(mech.softening(e, fixed, P)) <= 1e-15)
        assert np.all(np.diff(mech.softening(fixed, e, P)) <= 1e-15)


def check_healing():
    """Healing: with no pull, warm strain decays by exp(-t/tau) on every cell, whatever the step; permanent strain never does."""
    t = toy()
    t.eps[:], t.warm[:] = 1.0, 0.8
    prob = t.problem()
    prob.normal_load[:] = 0.0
    for pieces in (1, 3):
        eps, warm = t.eps.copy(), t.warm.copy()
        for _ in range(pieces):
            out = rl.advance_rift(prob, eps, warm, 62.8e6 / pieces, P)
            eps, warm = out.eps, out.warm
        assert np.allclose(eps, 1.0) and np.allclose(warm, 0.8 * np.exp(-1.0), rtol=1e-6), (pieces, eps[:2], warm[:2])


def check_no_strain_off_band_and_latch():
    """No strain off the band; a cell at its latch strain does not strain further; a held cell (latch = strain) does not strain."""
    t = toy(pull_ratio=5.0)
    prob = t.problem()
    eps0 = t.eps.copy()
    eps0[t.chain[3]] = 1.0  # at the latch below
    prob.latch[t.chain[3]] = 1.0
    out = rl.advance_rift(prob, eps0, t.warm, 5 * MYR, P)
    outside = np.ones(len(eps0), dtype=bool)
    outside[prob.band] = False
    assert np.array_equal(out.eps[outside], eps0[outside]), "strain off the band"
    assert out.eps[t.chain[3]] == 1.0 and out.eps[t.chain].max() > 0.2, (out.eps[t.chain][:6],)
    assert np.all(out.eps <= prob.latch + 1e-12)


def check_plume_expiry():
    """Plume expiry: a step crossing the expiry equals its two parts run in sequence, and differs from no plume."""
    q = replace(P, s_fault=0.0, s_warm=0.0, n_plume=0.5e12)
    t = toy(q, pull_ratio=1.2)
    prob = t.problem()
    prob.plume_left_yr = np.full(len(prob.chain), 0.5 * MYR)
    crossed = rl.advance_rift(prob, t.eps, t.warm, 1 * MYR, q)
    first = rl.advance_rift(prob, t.eps, t.warm, 0.5 * MYR, replace(q))
    prob2 = t.problem()  # after the expiry: no plume
    second = rl.advance_rift(prob2, first.eps, first.warm, 0.5 * MYR, q)
    assert np.allclose(crossed.eps, second.eps, rtol=1e-6, atol=1e-9), (crossed.eps[t.chain][:5], second.eps[t.chain][:5])
    never = rl.advance_rift(prob2, t.eps, t.warm, 1 * MYR, q)
    always = rl.advance_rift(replace(prob, plume_left_yr=np.full(len(prob.chain), 2 * MYR)), t.eps, t.warm, 1 * MYR, q)
    assert never.eps[t.chain].max() < crossed.eps[t.chain].max() < always.eps[t.chain].max()


def check_flank_conservation_and_held():
    """Flank volume and tracer conservation, floors, craton booking, recipients never strain, far recipients with no reach limit,
    and a band with no recipient anywhere is held."""
    rng = np.random.default_rng(0)
    nx, ny = 16, 12
    g, n = lattice(nx, ny), nx * ny
    areas = rng.uniform(0.8, 1.2, n) * 4e9
    cols = columns_for(n, rng)
    band = (np.arange(n) % nx >= 6) & (np.arange(n) % nx <= 8)
    plan = mech.flank_plan(g, band, np.zeros(n, bool), P.flank_rows)
    assert not plan.held.any()
    d_eps = np.zeros(n)
    d_eps[band] = rng.uniform(0.05, 1.5, band.sum())
    out = mech.thin_with_flank_transfer(areas, plan, d_eps, cols)
    for name in ("hc", "hm", "material", "cover", "cover_material", "plateau", "arc", "moment"):
        before, after = float(np.dot(cols[name], areas)), float(np.dot(getattr(out, name), areas))
        assert abs(after - before) <= 1e-12 * abs(before), (name, before, after)
    lost = float(np.dot(cols["craton"] - out.craton, areas))
    assert abs(lost - out.craton_lost_m3) <= 1e-9 * max(lost, 1.0) and lost >= 0.0
    assert np.all(out.hc[band] <= cols["hc"][band] + 1e-9) and np.all(out.hc[~band] >= cols["hc"][~band] - 1e-9)
    # A band that is the whole plate has no recipient: every cell is held, and nothing moves.
    everything = np.ones(n, bool)
    plan_all = mech.flank_plan(g, everything, np.zeros(n, bool), P.flank_rows)
    assert plan_all.held.all()
    none = mech.thin_with_flank_transfer(areas, plan_all, np.full(n, 0.3), cols)
    assert np.array_equal(none.hc, cols["hc"]) and none.craton_lost_m3 == 0.0
    # The only recipients are far from the donors: the reach widens with no limit (here 9 hops) and volume is conserved.
    strip, n2 = lattice(14, 1), 14
    donors = np.zeros(n2, bool)
    donors[:10] = True
    plan2 = mech.flank_plan(strip, donors, np.zeros(n2, bool), P.flank_rows)
    assert not plan2.held.any() and plan2.hops[0] == 10 and plan2.hops[9] == 2, plan2.hops
    cols2 = columns_for(n2, rng)
    out2 = mech.thin_with_flank_transfer(np.full(n2, 4e9), plan2, np.where(donors, 0.5, 0.0), cols2)
    assert abs(np.dot(out2.hc - cols2["hc"], np.full(n2, 4e9))) <= 1e-9 * np.dot(cols2["hc"], np.full(n2, 4e9))
    # Floored donors give up only what they lost.
    cols3 = columns_for(n, rng)
    cols3["hc"][:], cols3["hm"][:] = 5_300.0, 2_100.0
    d3 = np.zeros(n)
    d3[band] = 2.0
    out3 = mech.thin_with_flank_transfer(areas, plan, d3, cols3)
    assert np.all(out3.hc[band] >= rheology.RIFT_CRITICAL_THICKNESS_M - 1e-6)
    assert np.all(out3.hm[band] >= lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M - 1e-6)
    assert abs(np.dot(out3.hc, areas) - np.dot(cols3["hc"], areas)) < 1e-9 * np.dot(cols3["hc"], areas)


def check_held_band_never_strains():
    """A band with no recipient, held over many steps: its strain, softening and latch never change."""
    t = toy(pull_ratio=6.0)
    prob = t.problem()
    in_band = np.ones(len(t.points), dtype=bool)
    plan = mech.flank_plan(t.graph, in_band, np.zeros(len(t.points), bool), P.flank_rows)
    assert plan.held.all()
    prob.latch = t.eps.copy()  # held: latch = current strain
    eps = t.eps.copy()
    for _ in range(40):
        eps = rl.advance_rift(prob, eps, t.warm, 5 * MYR, P).eps
    assert np.array_equal(eps, t.eps)


def check_section_load_branches():
    """The section load's branches: an open line with positive pull splits when the gates allow; an open or partly broken line with
    zero pull neither grows nor splits; open lines keep finite rates and strains and stop at their latch; no infinite value."""
    assert rl.section_load(0.0, 1.0, P) == (0.0, False)
    assert rl.section_load(-1.0, 0.0, P) == (0.0, False)
    lam, is_open = rl.section_load(1.0, 0.0, P)
    assert is_open and np.isfinite(lam) and lam == P.r_max * P.lam_grow
    assert rl.growth_ratio(lam, True, P) == P.r_max and rl.growth_ratio(10.0, False, P) == P.r_max
    assert rl.rank_key(1e9, 0.1, True) > rl.rank_key(1.0, 99.0, False) and rl.rank_key(2.0, 0, True) > rl.rank_key(1.0, 0, True)
    assert rl.ligament_gives(1.0, 0.0) and not rl.ligament_gives(0.0, 0.0) and not rl.ligament_gives(-1.0, 0.0)
    assert rl.ligament_gives(2.0, 1.0) and not rl.ligament_gives(1.0, 2.0)
    t = toy(pull_ratio=2.0)
    prob = t.problem()
    prob.intact[:] = False  # an open line
    out = rl.advance_rift(prob, t.eps, t.warm, 5 * MYR, P)
    assert out.open_start and np.isfinite(out.eps).all() and out.eps[t.chain].max() > 0.0
    assert np.all(out.eps <= prob.latch + 1e-12)
    prob.normal_load[:] = 0.0  # no pull: neither grows nor splits
    idle = rl.advance_rift(prob, t.eps, t.warm, 5 * MYR, P)
    assert np.array_equal(idle.eps, t.eps) and not idle.open_start
    assert not rl.ligament_gives(0.0, 0.0)


def check_substep_bound_and_kfail():
    """The sub-step bound (no band cell moves by more than d_eps_max) and the K_fail error."""
    t = toy(pull_ratio=P.r_max * 1.2)
    prob = t.problem()
    out = rl.advance_rift(prob, t.eps, t.warm, 1 * MYR, P)
    bound = int(np.ceil(1 * MYR * P.eps_dot0 * (P.r_max - 1) ** P.n / P.d_eps_max)) + 2 * len(prob.chain)
    assert out.substeps <= bound, (out.substeps, bound)
    try:
        rl.advance_rift(prob, t.eps, t.warm, 60 * MYR, replace(P, k_fail=3))
    except RuntimeError as error:
        assert "K_fail" in str(error)
    else:
        raise AssertionError("K_fail did not stop a run that needed more sub-steps")


def check_ligament_normal_load():
    """Path-normal loading: tension along the line never counts; across it does; compression across it does not."""
    n = 10
    normal = np.tile([0.0, 1.0], (n, 1))  # line along x, normal along y (tangent-basis a, b)
    along = np.tile([2e11, 0.0, 0.0], (n, 1))
    across = np.tile([0.0, 2e11, 0.0], (n, 1))
    assert np.allclose(rl.path_normal_load(along, normal), 0.0)
    assert np.allclose(rl.path_normal_load(across, normal), 2e11)
    assert np.all(np.maximum(rl.path_normal_load(-across, normal), 0.0) == 0.0)


def check_directional_sign():
    """Directional load signs: balanced forces pushing a plate's two halves together give no extension; pulling
    them apart does, and both agree with the sign of the cut resultants."""
    lat = np.linspace(-0.05, 0.05, 15)
    x, y = np.meshgrid(lat, lat, indexing="ij")
    points = np.column_stack([np.ones(x.size), x.ravel(), y.ravel()])
    points /= np.linalg.norm(points, axis=1)[:, None]
    areas = np.full(len(points), (0.0072 * mech.RADIUS_M) ** 2)
    sign = np.where(points[:, 1] > 0, 1.0, -1.0)
    for direction, expect_extension in ((+1.0, True), (-1.0, False)):  # + : each half pulled away from the middle
        force = direction * 1e11 * sign[:, None] * np.array([0.0, 1.0, 0.0])[None, :]
        tensor = P.k_cut * proto.cut_tensor(points, areas, force, 0.0072 * mech.RADIUS_M)
        lam, lam_min, _ = proto.principal(tensor)
        middle = np.abs(points[:, 1]) < 0.01
        top, bottom = float(np.median(lam[middle])), float(np.median(lam_min[middle]))
        if expect_extension:
            assert top > 0 and top > 5 * abs(bottom), (top, bottom)
        else:
            assert bottom < 0 and top < 0.15 * abs(bottom), ("compression registered as extension", top, bottom)


def _sections_for(force_sign: float):
    t = toy(nx=25, ny=15)
    pts = t.points
    side = np.where(pts[:, 1] > 0, 1.0, -1.0)  # +y half pulled +y, -y half pulled -y (tension across the middle)
    force = force_sign * 1e11 * side[:, None] * np.array([0.0, 1.0, 0.0])[None, :]
    return t, rl.section_profile(pts, t.areas, force, P.k_cut)


def check_section_pull_manufactured():
    """The section pull of manufactured force fields: two halves pulled apart peak at the middle section with the force on a half
    (times k_cut); pushed together gives no corridor; sections that leave less than a_min on a side are never candidates."""
    t, sections = _sections_for(+1.0)
    d_y = int(np.argmin(np.abs(sections.normals @ np.array([1.0, 0.0]))))  # the direction whose normal is along the lattice's y axis
    # the lattice's v axis maps to the plane frame's second axis; find the direction normal to the long axis by its peak size
    peaks = [float(np.max(sections.pull[:, d])) for d in range(rl.N_DIRECTIONS)]
    best = int(np.argmax(peaks))
    expected = P.k_cut * 1e11 * float(np.sum(t.points[:, 1] > 0))
    assert abs(peaks[best] - expected) <= 0.1 * expected, (peaks, expected, d_y)
    corridors = rl.candidate_corridors(sections, P)
    assert corridors and corridors[0].peak_share >= P.a_min - 1e-9
    for corridor in corridors:
        smaller = np.minimum(sections.far_share[:, corridor.direction], 1 - sections.far_share[:, corridor.direction])
        assert np.all(smaller[corridor.cells] >= P.a_min - 1e-9)
    _, pushed = _sections_for(-1.0)
    # The six-direction cut leaks a few per cent of a uniaxial load into the perpendicular direction (as in
    # check_directional_sign), so a compressed plate may keep a tiny positive peak: it must be far below the push.
    compressive = max(abs(float(np.min(np.where(np.minimum(pushed.far_share[:, d], 1 - pushed.far_share[:, d]) >= P.a_min, pushed.pull[:, d], np.inf))))
                      for d in range(rl.N_DIRECTIONS))
    for corridor in rl.candidate_corridors(pushed, P):
        assert corridor.peak_pull < 0.05 * compressive, ("compression registered as a pull", corridor.peak_pull, compressive)


def _line_lattice(nx=24, ny=13):
    """A lattice with a corridor across its middle band, for the line placement checks: returns (graph, boundary, strike, cost grid)."""
    graph = lattice(nx, ny)
    index = np.arange(nx * ny).reshape(ny, nx)
    boundary = np.zeros((ny, nx), bool)
    boundary[0], boundary[-1], boundary[:, 0], boundary[:, -1] = True, True, True, True
    strike = np.broadcast_to(np.arange(nx, dtype=float), (ny, nx)).copy()  # the line runs along x
    corridor = np.zeros((ny, nx), bool)
    corridor[3:10] = True
    return graph, index, boundary.ravel(), strike.ravel(), corridor.ravel(), nx, ny


def check_line_placement():
    """Line placement: around a craton, along a belt, once found it is found again after holding a rift; the line is a valid chain
    from end to end of the corridor; a corridor that is not connected end to end gives no line."""
    graph, index, boundary, strike, corridor, nx, ny = _line_lattice()
    cost = np.ones(nx * ny)
    straight = rl.find_line(graph, corridor, boundary, strike, cost)
    assert straight is not None and len(straight) >= nx - 3 and boundary[straight[0]] and boundary[straight[-1]]
    assert np.all(graph[straight[:-1], straight[1:]] > 0), "not an edge-connected chain"
    assert len(set(straight.tolist())) == len(straight)
    craton = np.zeros((ny, nx), bool)
    craton[4:9, 8:15] = True
    around = rl.find_line(graph, corridor, boundary, strike, np.where(craton.ravel(), 50.0, 1.0))
    assert not craton.ravel()[around].any(), "the line crossed the craton"
    belt = np.full((ny, nx), 1.0)
    belt[8, :] = 0.2
    along_belt = rl.find_line(graph, corridor, boundary, strike, belt.ravel())
    assert np.mean(index[8][1:-1][:, None] == along_belt[None, :]) > 0.8 or (np.isin(along_belt, index[8]).mean() > 0.8)
    cost_broken = np.where(np.isin(np.arange(nx * ny), along_belt), 0.0, belt.ravel())
    again = rl.find_line(graph, corridor, boundary, strike, cost_broken)
    assert np.array_equal(again, along_belt), "a line holding a rift was not found again"
    split = corridor.copy()
    split.reshape(ny, nx)[:, 11:13] = False
    assert rl.find_line(graph, split, boundary, strike, cost) is None, "an unconnected corridor gave a line"


def check_path_lengths():
    """ds of a path is half the centre distances to its neighbours (ends: half of one)."""
    from app import thin_shell

    pts = np.column_stack([np.ones(3), np.array([0.0, 0.01, 0.02]), np.zeros(3)])
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    first, second = thin_shell._tangent_basis(pts)
    geometry = rl.path_geometry(pts, np.arange(3), first, second)
    spacing = 0.01 * mech.RADIUS_M
    assert np.allclose(geometry.ds / spacing, [0.5, 1.0, 0.5], rtol=1e-3), geometry.ds / spacing
    assert np.isclose(geometry.length, 2 * spacing, rtol=1e-3) and np.isclose(geometry.straight, geometry.length, rtol=1e-3)


def check_split_masks():
    """Stray-piece attachment and its limit, the a_min daughter minimum, and a valid two-daughter split."""
    nx, ny = 20, 10
    graph = lattice(nx, ny)
    areas = np.ones(nx * ny)
    index = np.arange(nx * ny).reshape(ny, nx)
    line = index[:, 10]  # a vertical line splits the plate in two halves
    mask, why = rl.split_masks(graph, areas, line, P)
    assert mask is not None and 85 <= mask.sum() <= 115, why  # line cells tie on votes and go to side 0
    assert connected(graph, mask) and connected(graph, ~mask)
    # a stray piece (a 3x3 pocket, 4.5% of the plate) walled off by a ring: attached when the limit allows it, cancelling the split otherwise
    ring = np.unique(np.concatenate([index[2, 12:17], index[6, 12:17], index[2:7, 12], index[2:7, 16]]))
    cut = np.concatenate([line, ring])
    pocket = index[3:6, 13:16].ravel()
    mask3, why3 = rl.split_masks(graph, areas, cut, replace(P, a_stray=0.001))
    assert mask3 is None and why3 == "stray_too_large", why3
    mask3, why3 = rl.split_masks(graph, areas, cut, replace(P, a_stray=0.05))
    assert mask3 is not None and len({bool(mask3[c]) for c in pocket}) == 1 and mask3[pocket[0]] == mask3[index[5, 14] + 0] and not mask3[pocket[0]], why3
    # a_min: a line near the edge leaves a daughter under 15%
    near_edge = index[:, 2]
    mask4, why4 = rl.split_masks(graph, areas, near_edge, P)
    assert mask4 is None and why4 == "daughter_small", why4
    assert rl.split_masks(graph, areas, index[:, 2], replace(P, a_min=0.05))[0] is not None


def check_split_islands_and_deep_line():
    """Isolated cells of the plate (islands) are left out of the stray and connectivity tests and join the nearest daughter; a line
    three cells thick (its middle cells touch no daughter) still gives two connected daughters."""
    nx, ny = 20, 10
    graph = lattice(nx, ny)
    n = nx * ny
    index = np.arange(n).reshape(ny, nx)
    coo = graph.tocoo()
    padded = csr_matrix((coo.data, (coo.row, coo.col)), shape=(n + 1, n + 1))  # one extra cell with no edges
    areas = np.ones(n + 1)
    points = np.column_stack([np.arange(n + 1) % nx, np.arange(n + 1) // nx, np.zeros(n + 1)]).astype(float)
    points[n] = [2.0, 5.0, 0.0]  # next to the left half
    thick = index[:, 9:12].ravel()
    mask, why = rl.split_masks(padded, areas, thick, P, points=points)
    assert mask is not None, why
    assert mask[n] == mask[index[5, 2]] and connected(padded[:n, :n], mask[:n]) and connected(padded[:n, :n], ~mask[:n])


def connected(graph, mask) -> bool:
    from scipy.sparse.csgraph import connected_components

    return connected_components(graph[mask][:, mask], directed=False)[0] == 1


def check_tip_speed_against_law():
    """Tip start-up and steady speed against the discrete law under a constant Λ (the full sweep is in check_rift_line_d1.py)."""
    row = toy_mod.propagation(125e3, 125e3, 0.05, P)
    assert abs(row["nucleation_myr"] - row["law_nucleation_myr"]) <= 0.1 * row["law_nucleation_myr"], row
    assert row["within_20pct"], row


def check_capacity_and_disk_kernels():
    """The interpolated capacity matches a direct column evaluation, and disk averages are area-weighted rows of one."""
    import replay_gate2_capacity as gate2

    hc, hm, lag = np.array([35_000.0, 12_000.0, 8_000.0]), np.array([100_000.0, 40_000.0, 8_000.0]), np.array([0.0, 50.0, 0.0])
    kind = np.array([0, 2, 1])
    got = mech.capacity(kind, hc, hm, lag)
    for i, name in enumerate(mech.KINDS[k] for k in kind):
        direct = gate2._column_capacity(name, float(hc[i]), float(hm[i]), float(lag[i]))
        assert abs(got[i] - direct) <= 0.02 * direct, (name, got[i], direct)
    rng = np.random.default_rng(3)
    pts = rng.normal(size=(200, 3))
    pts /= np.linalg.norm(pts, axis=1)[:, None]
    pts = pts[pts[:, 2] > 0.8]
    areas = rng.uniform(1, 2, len(pts)) * 1e9
    disk = mech.disk_average_matrix(pts, areas, 400e3)
    assert np.allclose(np.asarray(disk.sum(axis=1)).ravel(), 1.0)
    assert np.allclose(disk @ np.full(len(pts), 3.0), 3.0)


# --- plate-level checks --------------------------------------------------------------------------


def world(seed):
    if "world" not in _WORLD:
        from app import world as world_module

        step_mod.install(P)
        w = world_module.generate_world(seed)
        w.simulate_climate_biomes = False
        _WORLD["world"] = w
    return _WORLD["world"]


def biggest(w):
    return max(w.plates, key=lambda p: p.node_count())


def check_remap_without_healing(seed):
    """Remap without healing: coarsening a quartet with one broken cell keeps it broken and its strain does not fall."""
    w = world(seed)
    plate = biggest(w)
    key = int(plate.cell_keys[len(plate.cell_keys) // 2])
    plate.refine_cells(np.array([key]))
    children = np.array([int(k) for k in __import__("app.sparse_quad_patch", fromlist=["x"]).child_cell_keys(np.array([key]))[0]])
    index = {int(k): i for i, k in enumerate(plate.cell_keys)}
    sel = np.array([index[int(c)] for c in children])
    eps, warm, broken = plate.collect("rift_strain"), plate.collect("rift_warm_strain"), plate.collect("rift_broken")
    eps[sel[0]], warm[sel[0]], broken[sel[0]] = 2.0, 2.0, True
    plate.set_fields_on_plate(rift_strain=eps, rift_warm_strain=warm, rift_broken=broken)
    plate.coarsen_cells(np.array([key]))
    at = int(np.flatnonzero(plate.cell_keys == key)[0])
    e, wm, b = plate.collect("rift_strain")[at], plate.collect("rift_warm_strain")[at], plate.collect("rift_broken")[at]
    assert e == 2.0 and b and 0.0 <= wm <= e, (e, wm, b)


def check_reset_on_eruption(seed):
    """Damage reset on eruption and fresh crust: a column that becomes exactly the reference column loses its damage."""
    w = world(seed)
    state = step_mod._STATE
    plate = biggest(w)
    n = plate.node_count()
    hc, hm = plate.collect("crustal_thickness_m"), plate.collect("mantle_lithosphere_thickness_m")
    hc[:5], hm[:5] = 12_345.0, 54_321.0
    plate.set_fields_on_plate(crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm)
    state.prev_reference[plate.plate_id] = (plate.cell_keys.copy(), step_mod._is_reference(plate.collect("crustal_thickness_m"), plate.collect("mantle_lithosphere_thickness_m")))
    hc[:3], hm[:3] = lithosphere.REFERENCE_HC_CONTINENTAL_M, lithosphere.REFERENCE_HM_CONTINENTAL_M  # erupted
    plate.set_fields_on_plate(crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm)
    eps, warm, broken = np.ones(n), np.ones(n) * 0.5, np.ones(n, bool)
    step_mod.reset_after_eruption(state, plate, eps, warm, broken)
    assert np.all(eps[:3] == 0) and np.all(~broken[:3]) and np.all(eps[3:] == 1.0) and np.all(broken[3:]), (eps[:6], broken[:6])


def check_margin_insertion_share(seed):
    """Margin insertion inherits only its stretched share of the donors' damage, never the broken flag."""
    w = world(seed)
    plate = biggest(w)
    n = plate.node_count()
    plate.set_fields_on_plate(rift_strain=np.full(n, 1.0), rift_warm_strain=np.full(n, 0.6), rift_broken=np.ones(n, bool))
    rifted = np.array([0, 1])
    original = step_mod._ORIGINAL.get("open_rift")
    step_mod._ORIGINAL["open_rift"] = lambda *a, **k: None
    try:
        step_mod.open_rift(plate, w, 0, rifted, rifted, np.array([0.25, 0.8]), 1e6)
    finally:
        step_mod._ORIGINAL["open_rift"] = original
    e, wm, b = plate.collect("rift_strain"), plate.collect("rift_warm_strain"), plate.collect("rift_broken")
    assert np.allclose(e[rifted], [0.25, 0.8]) and np.allclose(wm[rifted], [0.15, 0.48]) and not b[rifted].any(), (e[rifted], wm[rifted], b[rifted])
    assert np.all(e[2:] == 1.0)


def _fresh(plate):
    n = plate.node_count()
    plate.set_fields_on_plate(rift_strain=np.zeros(n), rift_warm_strain=np.zeros(n), rift_broken=np.zeros(n, bool))


def check_latch_no_eruption_and_applied_strain(seed):
    """The broken latch: overloaded line cells latch at the floor with no interior eruption, and ``rift_strain`` never exceeds the
    strain applied to the column."""
    w = world(seed)
    plate = biggest(w)
    state = step_mod._STATE
    old = state.params
    state.params = Params(eps_dot0=5e-6, k_cap=0.05, lam_grow=0.01)
    try:
        _fresh(plate)
        volcano_before = plate.collect("is_volcano").copy()
        hc0, hm0 = plate.collect("crustal_thickness_m").copy(), plate.collect("mantle_lithosphere_thickness_m").copy()
        step_mod.damage_step(w, plate, 1e6, state)
        one_eps, one_hc = plate.collect("rift_strain"), plate.collect("crustal_thickness_m")
        donors = one_eps > 0  # after one step no donor has received volume yet: its strain is what thinned its crust
        assert donors.any() and np.all(one_eps[donors] <= np.log(hc0[donors] / np.maximum(one_hc[donors], 1.0)) + 1e-6)
        for _ in range(5):
            step_mod.damage_step(w, plate, 1e6, state)
        hc, volcano, broken = plate.collect("crustal_thickness_m"), plate.collect("is_volcano"), plate.collect("rift_broken")
        assert broken.any(), "nothing latched under an overload"
        assert np.array_equal(volcano, volcano_before), "interior thinning ignited a volcano"
        assert hc[broken].min() >= rheology.RIFT_CRITICAL_THICKNESS_M - 1e-6
        assert not np.any(step_mod._is_reference(hc, plate.collect("mantle_lithosphere_thickness_m")) & ~step_mod._is_reference(hc0, hm0) & broken)
    finally:
        state.params = old
        _fresh(plate)


def check_damage_runs_for_every_plate(seed):
    """The damage update on a plate after another plate splits and on a daughter in cooldown."""
    w = world(seed)
    state = step_mod._STATE
    state.pass_key = None
    young = w.plates[-1]
    young.reset_age()
    step_mod.maybe_split_plate(w, w.plates[0], 1e6)
    assert set(state.steps) == {p.plate_id for p in w.plates}, "a plate was skipped by the damage pass"
    assert young.plate_id in state.steps and young.age_years < 1.0
    again = len(state.steps)
    step_mod.maybe_split_plate(w, w.plates[1], 1e6)  # same step: no second pass
    assert len(state.steps) == again and state.pass_key == w.elapsed_years


def check_broken_line_splits_and_f_break(seed):
    """The integrated ligament rule and the ``f_break`` rule on a real plate: a fully broken line with pull splits it into two
    connected plates that partition the parent and keep its ledgers; a line with less than ``f_break`` broken does not split."""
    w = world(seed)
    state = step_mod._STATE
    plate = max((p for p in w.plates if p.node_count() > 400), key=lambda p: p.node_count())
    plate._age_years = 1e9
    _fresh(plate)
    prep = step_mod.prepare_plate(w, plate, state.params, 1e6)
    assert prep.line is not None, "no rift line on the largest plate"
    chain = prep.line.chain
    eps = np.zeros(plate.node_count())
    broken = np.zeros(plate.node_count(), bool)
    state.pass_key = None
    # Few broken cells: P >= S can hold (a very strong pull) but the share is below f_break.
    few = chain[: max(1, len(chain) // 5)]
    eps[few], broken[few] = 1.9, True
    plate.set_fields_on_plate(rift_strain=eps, rift_broken=broken)
    old = state.params
    state.params = replace(old, k_cap=1e-3)  # the pull far exceeds the strength
    try:
        step_mod.damage_step(w, plate, 1e6, state)
        step = state.steps[plate.plate_id]
        assert step.line is not None and step.broken_share < old.f_break + 1e-9
        assert step_mod.breakthrough(w, plate, state) is None, "split with too little of the line broken"
        eps[chain], broken[chain] = 1.9, True
        plate.set_fields_on_plate(rift_strain=eps, rift_broken=broken)
        step_mod.damage_step(w, plate, 1e6, state)
        step = state.steps[plate.plate_id]
        assert step.broken_share >= old.f_break
        n_before = plate.node_count()
        crust_before = float(np.dot(plate.collect("crustal_thickness_m"), plate.node_areas_m2()))
        result = step_mod.breakthrough(w, plate, state)
        why = state.counters["rejected_because"]
        if result is None:  # a valid refusal (geometry); then the reason is recorded
            assert why, "no split and no recorded reason"
            return
        a, b = result
        assert a.node_count() + b.node_count() == n_before and a.plate_id != b.plate_id
        crust_after = sum(float(np.dot(d.collect("crustal_thickness_m"), d.node_areas_m2())) for d in (a, b))
        assert abs(crust_after - crust_before) <= 1e-9 * crust_before
        assert len(np.unique(np.concatenate([a.cell_keys, b.cell_keys]))) == n_before
        assert min(a.node_count(), b.node_count()) >= old.a_min * n_before - 1
    finally:
        state.params = old


def check_shadow_stats_skips_degenerate_plate():
    """A plate the thin-shell solver rejects (ValueError) is skipped and counted, not fatal (two stage 2 runs crashed on it, #888)."""
    import shadow_shell_check as shadow

    class World:
        plates = [object(), object()]

    original = shadow.d0.compute_record

    def reject(world, plate):
        raise ValueError("thin-shell solve requires at least four non-coplanar nodes")

    shadow.d0.compute_record = reject
    try:
        got = shadow.shadow_stats(World())
    finally:
        shadow.d0.compute_record = original
    assert got == {"plates": 0, "skipped_no_mesh": 2}, got


def check_oceanic_strength_follows_age():
    """Old ocean floor is stronger than young, is never thinner than the model's own lithosphere, and thins with rift strain."""
    p = Params()
    hc, hm = np.full(4, 11_000.0), np.full(4, 40_000.0)
    age = np.array([1e6, 25e6, 150e6, np.nan])
    got = mech.oceanic_strength_hm(hc, hm, age, np.zeros(4), p)
    assert got[0] == 40_000.0 and got[3] == 40_000.0, got  # young or undated floor keeps the model's thickness
    assert abs(got[1] - (13_000.0 * 5.0 - 11_000.0)) < 1e-6 and got[2] == p.ocean_plate_thickness_m - 11_000.0, got
    strained = mech.oceanic_strength_hm(hc, hm, age, np.full(4, 0.5), p)
    assert strained[2] < got[2] and np.all(strained >= hm), strained
    kind = np.full(4, 2)
    cap = mech.capacity(kind, hc, got, np.zeros(4))
    assert cap[2] > cap[1] > cap[0], cap


def check_strength_map_arms_default_off():
    """The strength-map arms change nothing at their defaults, and each moves strength the way it says."""
    p = Params()
    assert p.craton_root_m == 0.0 and p.mark_power == 1.0 and p.route_mark_power == 1.0 and not p.oceanic_age_strength
    kind, hc, lag = np.zeros(2, int), np.full(2, 40_000.0), np.zeros(2)
    cap = mech.capacity(kind, hc, np.array([100_000.0, 180_000.0]), lag)
    assert cap[1] > 2.0 * cap[0], cap  # a 180 km root at least doubles a craton's capacity
    mult = mech.inherited_multiplier(np.array([True, False]), np.zeros(2, bool), np.zeros(2, bool), np.zeros(2, bool))
    deep = np.clip(mult ** 2.8, mech.MULTIPLIER_FLOOR, 1.0)
    assert abs(deep[0] - 0.65 ** 2.8) < 1e-12 and deep[1] == 1.0, deep


CHECKS = [check_shadow_stats_skips_degenerate_plate, check_strength_map_arms_default_off, check_oceanic_strength_follows_age, check_strength_bounds, check_softening_monotone, check_healing, check_no_strain_off_band_and_latch, check_plume_expiry,
          check_flank_conservation_and_held, check_held_band_never_strains, check_section_load_branches, check_substep_bound_and_kfail,
          check_ligament_normal_load, check_directional_sign, check_section_pull_manufactured, check_line_placement, check_path_lengths,
          check_split_masks, check_split_islands_and_deep_line, check_tip_speed_against_law, check_capacity_and_disk_kernels]
PLATE_CHECKS = [check_remap_without_healing, check_reset_on_eruption, check_margin_insertion_share,
                check_latch_no_eruption_and_applied_strain, check_damage_runs_for_every_plate, check_broken_line_splits_and_f_break]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seed", type=int, default=106533960)
    parser.add_argument("--only", nargs="*")
    parser.add_argument("--synthetic", action="store_true")
    args = parser.parse_args()
    failed = 0
    for check in CHECKS + ([] if args.synthetic else PLATE_CHECKS):
        if args.only and check.__name__ not in args.only:
            continue
        try:
            check(args.seed) if check in PLATE_CHECKS else check()
            print(f"PASS  {check.__name__}")
        except Exception:  # noqa: BLE001
            failed += 1
            print(f"FAIL  {check.__name__}\n{traceback.format_exc()}")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
