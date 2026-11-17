import numpy as np
import pytest
from scipy.sparse import csr_matrix

from app import lithosphere, rheology
from app import rift_strength as strength
from app.rift_strength import RiftLineParams

P = RiftLineParams()


def grid_edges(nx: int, ny: int):
    """Both directions of every edge of the 4-connected nx x ny grid, as (rows, cols) index arrays."""
    index = np.arange(nx * ny).reshape(ny, nx)
    rows, cols = [], []
    for a, b in ((index[:, :-1], index[:, 1:]), (index[:-1, :], index[1:, :])):
        rows += [a.ravel(), b.ravel()]
        cols += [b.ravel(), a.ravel()]
    return np.concatenate(rows), np.concatenate(cols)


def lattice(nx: int, ny: int) -> csr_matrix:
    r, c = grid_edges(nx, ny)
    return csr_matrix((np.full(len(r), 1000.0), (r, c)), shape=(nx * ny, nx * ny))


def _columns(n, rng):
    return {"hc": rng.uniform(20e3, 45e3, n), "hm": rng.uniform(60e3, 120e3, n), "material": rng.uniform(0, 30e3, n),
            "craton": rng.uniform(0, 10e3, n), "cover": rng.uniform(0, 500, n), "cover_material": rng.uniform(0, 200, n),
            "plateau": rng.uniform(0, 3e3, n), "arc": rng.uniform(0, 2e3, n), "moment": rng.uniform(0, 1e9, n)}


def test_overrides_replace_starting_values_and_reject_mistakes():
    p = RiftLineParams.with_overrides({"k_cap": 0.14, "flank_rows": 3, "oceanic_age_strength": False})
    assert (p.k_cap, p.flank_rows, p.oceanic_age_strength) == (0.14, 3, False)
    assert RiftLineParams.with_overrides(None) == P and RiftLineParams.with_overrides({}) == P
    for bad in ({"k_cup": 1.0}, {"k_cap": -1.0}, {"k_cap": float("nan")}, {"k_cap": True}, {"oceanic_age_strength": 1.0}, {"k_cap": "1"}):
        with pytest.raises(ValueError):
            RiftLineParams.with_overrides(bad)


def test_capacity_reproduces_the_gate2_anchor_and_orders_rock_as_the_envelope_does():
    cap = lambda kind, hc, hm: float(strength.capacity(np.array([kind]), np.array([hc]), np.array([hm]), np.zeros(1))[0])  # noqa: E731
    # RFC 0018's anchor: the 100 km reference continental column (35 km crust over 65 km of mantle lithosphere).
    assert cap(0, 35e3, 65e3) == pytest.approx(6.28e12, rel=0.01)
    assert cap(0, 40e3, 180e3) > 2.0 * cap(0, 40e3, 100e3)  # a craton root at least doubles capacity
    assert cap(2, 11e3, 84e3) > cap(0, 35e3, 65e3) > cap(2, 11e3, 40e3) > cap(1, 11e3, 8e3)
    # Outside the table the nearest edge is used, so no column gets a negative or non-finite strength.
    assert np.isfinite(cap(0, 500.0, 500e3)) and cap(0, 500.0, 500e3) > 0.0


def test_ocean_floor_strength_follows_age_and_never_drops_below_the_model_thickness():
    hc, hm = np.full(4, 11_000.0), np.full(4, 40_000.0)
    age = np.array([1e6, 25e6, 150e6, np.nan])
    got = strength.oceanic_strength_hm(hc, hm, age, np.zeros(4), P)
    assert got[0] == 40_000.0 and got[3] == 40_000.0  # young or undated floor keeps the model's thickness
    assert got[1] == pytest.approx(13_000.0 * 5.0 - 11_000.0) and got[2] == P.ocean_plate_thickness_m - 11_000.0
    strained = strength.oceanic_strength_hm(hc, hm, age, np.full(4, 0.5), P)
    assert strained[2] < got[2] and np.all(strained >= hm)
    caps = strength.capacity(np.full(4, 2), hc, got, np.zeros(4))
    assert caps[2] > caps[1] > caps[0]


def test_softening_and_weakness_multipliers_stay_in_their_bounds():
    eps = np.array([0.0, 0.3, 1.0, 5.0])
    w = strength.softening(eps, eps, P)
    assert w[0] == 1.0 and np.all(np.diff(w) <= 0.0) and w[-1] == pytest.approx(1.0 - min(P.s_max, P.s_fault + P.s_warm))
    assert np.all(strength.softening(eps, eps, RiftLineParams(s_warm=0.8)) >= 1.0 - P.s_max - 1e-12)
    on, off = np.array([True, False]), np.array([False, False])
    assert strength.inherited_multiplier(on, off, off, off).tolist() == [strength.MULT_BELT, 1.0]
    assert strength.inherited_multiplier(on, on, on, on)[0] == pytest.approx(0.65 ** 3 * 0.60)
    deep = strength.inherited_multiplier(on, on, on, on, mark_power=2.8)
    assert deep[0] == strength.MULTIPLIER_FLOOR and deep[1] == 1.0
    assert strength.inherited_multiplier(on, off, off, off, mark_power=2.8)[0] == pytest.approx(0.65 ** 2.8)


def test_disk_average_rows_sum_to_one_and_include_the_cell():
    rng = np.random.default_rng(1)
    points = rng.normal(size=(200, 3)) * 0.05 + np.array([1.0, 0.0, 0.0])
    points /= np.linalg.norm(points, axis=1)[:, None]
    disk = strength.disk_average_matrix(points, rng.uniform(1.0, 2.0, 200), 150e3)
    assert np.allclose(np.asarray(disk.sum(axis=1)).ravel(), 1.0)
    assert np.all(disk.diagonal() > 0.0) and np.all(disk.diagonal() <= 1.0)


def test_line_cut_load_is_tension_when_halves_are_pulled_apart_and_not_when_pushed_together():
    lat = np.linspace(-0.05, 0.05, 15)
    x, y = np.meshgrid(lat, lat, indexing="ij")
    points = np.column_stack([np.ones(x.size), x.ravel(), y.ravel()])
    points /= np.linalg.norm(points, axis=1)[:, None]
    areas = np.full(len(points), (0.0072 * strength.RADIUS_M) ** 2)
    side = np.where(points[:, 1] > 0, 1.0, -1.0)
    middle = np.abs(points[:, 1]) < 0.01
    for direction in (+1.0, -1.0):  # +: each half pulled away from the middle
        force = direction * 1e11 * side[:, None] * np.array([0.0, 1.0, 0.0])[None, :]
        balanced = strength.balanced_tangent_force(points, areas, force)
        tensor = strength.cut_tensor(points, areas, balanced, 0.0072 * strength.RADIUS_M)
        mean, radius = 0.5 * (tensor[:, 0] + tensor[:, 1]), np.hypot(0.5 * (tensor[:, 0] - tensor[:, 1]), tensor[:, 2])
        top, bottom = float(np.median((mean + radius)[middle])), float(np.median((mean - radius)[middle]))
        if direction > 0:
            assert top > 0 and top > 5 * abs(bottom)
        else:  # the six-direction fit leaks a little of a uniaxial load sideways; it must stay far below the push
            assert bottom < 0 and top < 0.15 * abs(bottom)


def test_flank_transfer_conserves_every_column_and_books_craton_to_its_sink():
    rng = np.random.default_rng(0)
    nx, ny = 16, 12
    graph, n = lattice(nx, ny), nx * ny
    areas = rng.uniform(0.8, 1.2, n) * 4e9
    cols = _columns(n, rng)
    band = (np.arange(n) % nx >= 6) & (np.arange(n) % nx <= 8)
    plan = strength.flank_plan(graph, band, np.zeros(n, bool), P.flank_rows)
    assert not plan.held.any()
    d_eps = np.zeros(n)
    d_eps[band] = rng.uniform(0.05, 1.5, band.sum())
    out = strength.thin_with_flank_transfer(areas, plan, d_eps, cols)
    for name in ("hc", "hm", "material", "cover", "cover_material", "plateau", "arc", "moment"):
        assert float(np.dot(getattr(out, name), areas)) == pytest.approx(float(np.dot(cols[name], areas)), rel=1e-12)
    lost = float(np.dot(cols["craton"] - out.craton, areas))
    assert lost >= 0.0 and out.craton_lost_m3 == pytest.approx(lost, rel=1e-9)
    assert np.all(out.hc[band] <= cols["hc"][band] + 1e-9) and np.all(out.hc[~band] >= cols["hc"][~band] - 1e-9)


def test_a_band_with_no_recipient_is_held_and_a_far_recipient_is_still_reached():
    rng = np.random.default_rng(0)
    nx, ny = 16, 12
    graph, n = lattice(nx, ny), nx * ny
    areas, cols = np.full(n, 4e9), _columns(n, rng)
    plan_all = strength.flank_plan(graph, np.ones(n, bool), np.zeros(n, bool), P.flank_rows)
    assert plan_all.held.all()
    untouched = strength.thin_with_flank_transfer(areas, plan_all, np.full(n, 0.3), cols)
    assert np.array_equal(untouched.hc, cols["hc"]) and untouched.craton_lost_m3 == 0.0
    # The only recipients are far from the donors: the reach widens with no limit and volume is conserved.
    strip, n2 = lattice(14, 1), 14
    donors = np.zeros(n2, bool)
    donors[:10] = True
    plan = strength.flank_plan(strip, donors, np.zeros(n2, bool), P.flank_rows)
    assert not plan.held.any() and plan.hops[0] == 10 and plan.hops[9] == 2
    cols2 = _columns(n2, rng)
    out = strength.thin_with_flank_transfer(np.full(n2, 4e9), plan, np.where(donors, 0.5, 0.0), cols2)
    assert float(out.hc.sum()) == pytest.approx(float(cols2["hc"].sum()), rel=1e-12)


def test_a_floored_donor_gives_up_only_what_it_lost():
    rng = np.random.default_rng(0)
    nx, ny = 16, 12
    graph, n = lattice(nx, ny), nx * ny
    areas, cols = rng.uniform(0.8, 1.2, n) * 4e9, _columns(n, rng)
    cols["hc"][:], cols["hm"][:] = 5_300.0, 2_100.0
    band = (np.arange(n) % nx >= 6) & (np.arange(n) % nx <= 8)
    plan = strength.flank_plan(graph, band, np.zeros(n, bool), P.flank_rows)
    out = strength.thin_with_flank_transfer(areas, plan, np.where(band, 2.0, 0.0), cols)
    assert np.all(out.hc[band] >= rheology.RIFT_CRITICAL_THICKNESS_M - 1e-6)
    assert np.all(out.hm[band] >= lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M - 1e-6)
    assert float(np.dot(out.hc, areas)) == pytest.approx(float(np.dot(cols["hc"], areas)), rel=1e-12)
