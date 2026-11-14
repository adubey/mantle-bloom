"""Unit-level checks of the RFC 0018 damage-band prototype pieces (D1, #888). Production tests come with #846.

Usage (repo root, this worktree's venv; one small world is generated, about a minute):
    backend/.venv/bin/python bin/debug/check_rift_damage_prototype.py [--seed 106533960]

Each ``check_*`` function asserts one of the RFC's acceptance tests (the matching bullet is in its docstring) and
the script prints PASS or the failure for each. Exit status is non-zero if any fails. Synthetic checks need no
world; the plate-level checks share one generated development world.
"""

import argparse
import sys
import traceback
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

sys.path[:0] = [str(Path(__file__).resolve().parent)]

import rift_damage_mechanism as mech  # noqa: E402
import rift_damage_prototype as proto  # noqa: E402
import rift_damage_step as step_mod  # noqa: E402
from rift_damage_mechanism import Params  # noqa: E402

from app import lithosphere, rheology  # noqa: E402

P = Params()
_WORLD = {}


def lattice(nx: int, ny: int):
    """4-connected grid graph with unit spacing (m = 1000) and unit areas."""
    index = np.arange(nx * ny).reshape(ny, nx)
    rows, cols = [], []
    for a, b in ((index[:, :-1], index[:, 1:]), (index[:-1, :], index[1:, :])):
        rows += [a.ravel(), b.ravel()]
        cols += [b.ravel(), a.ravel()]
    r, c = np.concatenate(rows), np.concatenate(cols)
    return csr_matrix((np.full(len(r), 1000.0), (r, c)), shape=(nx * ny, nx * ny))


def columns_for(n, rng):
    return {"hc": rng.uniform(20e3, 45e3, n), "hm": rng.uniform(60e3, 120e3, n), "material": rng.uniform(0, 30e3, n),
            "craton": rng.uniform(0, 10e3, n), "cover": rng.uniform(0, 500, n), "cover_material": rng.uniform(0, 200, n),
            "plateau": rng.uniform(0, 3e3, n), "arc": rng.uniform(0, 2e3, n), "moment": rng.uniform(0, 1e9, n)}


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


def _strain(duration, pieces, *, load=2.0e12, base=1.0e12, self_weight=0.3):
    n = 5
    eps, warm = np.linspace(0, 0.6, n), np.linspace(0, 0.3, n)
    nb_e, nb_w = np.full(n, 0.2), np.full(n, 0.1)
    sw, latch = np.full(n, self_weight), np.full(n, 2.0)
    for _ in range(pieces):
        r = mech.integrate_strain(np.full(n, load), np.full(n, base), eps, warm, nb_e, nb_w, sw, np.full(n, duration / pieces), latch, P)
        eps, warm = r.eps, r.warm
    return eps, warm


def check_strain_timestep():
    """Timestep-robust strain integration: a step and its halves agree to a few percent; the quadrupled step too."""
    one, halves, quarters = _strain(2e6, 1), _strain(2e6, 2), _strain(2e6, 4)
    for a, b in ((one, halves), (halves, quarters)):
        assert np.max(np.abs(a[0] - b[0])) < 0.03, (a[0], b[0])
    # Without softening feedback the integral is exact: rate = eps_dot0 * (load/strength - 1).
    flat = Params(s_fault=0.0, s_warm=0.0)
    r = mech.integrate_strain(np.array([2.0]), np.array([1.0]), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1), np.zeros(1),
                              np.array([1e6]), np.array([2.0]), flat)
    assert abs(r.eps[0] - flat.eps_dot0 * 1e6) < 1e-12, r.eps


def check_healing():
    """Healing: with no load warm strain decays by exp(-t/tau); permanent strain never does."""
    n = 3
    r = mech.integrate_strain(np.zeros(n), np.full(n, 1e12), np.full(n, 1.0), np.full(n, 0.8), np.zeros(n), np.zeros(n), np.zeros(n),
                              np.full(n, 62.8e6), np.full(n, 2.0), P)
    assert np.allclose(r.eps, 1.0) and np.allclose(r.warm, 0.8 * np.exp(-1.0), rtol=1e-6), (r.eps, r.warm)


def check_latch_stops_strain():
    """The broken latch: a cell at its latch strain does not strain further, whatever the load."""
    n = 2
    r = mech.integrate_strain(np.full(n, 1e14), np.full(n, 1e12), np.array([2.0, 1.0]), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n),
                              np.full(n, 5e6), np.array([2.0, 1.2]), P)
    assert r.eps[0] == 2.0 and abs(r.eps[1] - 1.2) < 1e-9, r.eps


def check_plume_expiry():
    """Plume expiry: the split step differs from the averaged load, and equals its two parts run in sequence."""
    n = 1
    strength, years, half = 1.0e12, 1e6, 5e5
    args = dict(base=np.full(n, strength), nb_eps=np.zeros(n), nb_warm=np.zeros(n), self_weight=np.zeros(n), latch_eps=np.full(n, 2.0), p=Params(s_fault=0.0, s_warm=0.0))
    z = np.zeros(n)
    averaged = mech.integrate_strain(np.full(n, 0.8 * strength + 0.2 * strength), eps=z, warm=z, duration=np.full(n, years), **args)
    assert averaged.eps[0] == 0.0, "averaged load equals the strength: no strain"
    first = mech.integrate_strain(np.full(n, 1.2 * strength), eps=z, warm=z, duration=np.full(n, half), **args)
    second = mech.integrate_strain(np.full(n, 0.8 * strength), eps=first.eps, warm=first.warm, duration=np.full(n, half), **args)
    assert second.eps[0] > 0.0 and np.isclose(second.eps[0], args["p"].eps_dot0 * 0.2 * half)


def check_flank_conservation():
    """Flank volume and tracer conservation, floors, craton booking, and recipients never strain."""
    rng = np.random.default_rng(0)
    nx, ny = 16, 12
    g, n = lattice(nx, ny), nx * ny
    areas = rng.uniform(0.8, 1.2, n) * 4e9
    cols = columns_for(n, rng)
    d_eps = np.zeros(n)
    band = (np.arange(n) % nx >= 6) & (np.arange(n) % nx <= 8)
    d_eps[band] = rng.uniform(0.05, 1.5, band.sum())
    out = mech.thin_with_flank_transfer(areas, g, d_eps, np.zeros(n, bool), cols, P)
    for name in ("hc", "hm", "material", "cover", "cover_material", "plateau", "arc", "moment"):
        before, after = float(np.dot(cols[name], areas)), float(np.dot(getattr(out, name), areas))
        assert abs(after - before) <= 1e-12 * abs(before), (name, before, after)
    lost = float(np.dot(cols["craton"] - out.craton, areas))
    assert abs(lost - out.craton_lost_m3) <= 1e-9 * max(lost, 1.0) and lost >= 0.0
    assert out.hc.min() >= min(cols["hc"].min(), rheology.RIFT_CRITICAL_THICKNESS_M) - 1e-6
    assert np.all(out.hm >= np.minimum(cols["hm"], lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M) - 1e-6)
    assert np.all(out.hc[band] <= cols["hc"][band] + 1e-9) and np.all(out.hc[~band] >= cols["hc"][~band] - 1e-9)
    assert not out.withheld.any()


def check_flank_no_recipient_and_floor():
    """Flank transfer with no recipient (withheld, nothing moves) and with floored donors (volume and bounds hold)."""
    rng = np.random.default_rng(1)
    g, n = lattice(6, 6), 36
    areas = np.full(n, 4e9)
    cols = columns_for(n, rng)
    everything = np.full(n, 0.3)
    out = mech.thin_with_flank_transfer(areas, g, everything, np.zeros(n, bool), cols, P)
    assert out.withheld.all() and np.array_equal(out.hc, cols["hc"]) and out.craton_lost_m3 == 0.0
    # One donor beyond the widest reach of any recipient: a long strip where only the far end is free.
    strip, n2 = lattice(12, 1), 12
    d = np.zeros(n2)
    d[:10] = 0.5
    cols2 = columns_for(n2, rng)
    out2 = mech.thin_with_flank_transfer(np.full(n2, 4e9), strip, d, np.zeros(n2, bool), cols2, P)
    # Recipients are cells 10-11; a donor reaches one within 2 hops, widening to 4: cells 6-9 are served, 0-5 withheld.
    assert out2.withheld[:6].all() and not out2.withheld[6:].any(), out2.withheld
    # Floored donors give up only what they lost.
    cols3 = columns_for(n, rng)
    cols3["hc"][:] = 5_300.0
    cols3["hm"][:] = 2_100.0
    d3 = np.zeros(n)
    d3[14:16] = 2.0
    out3 = mech.thin_with_flank_transfer(areas, g, d3, np.zeros(n, bool), cols3, P)
    assert np.all(out3.hc[14:16] >= rheology.RIFT_CRITICAL_THICKNESS_M - 1e-6)
    assert np.all(out3.hm[14:16] >= lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M - 1e-6)
    assert abs(np.dot(out3.hc, areas) - np.dot(cols3["hc"], areas)) < 1e-6 * np.dot(cols3["hc"], areas)
    assert np.isclose(np.dot(out3.hc - cols3["hc"], areas), 0.0, atol=1e-3 * np.dot(cols3["hc"], areas) * 1e-6)


def check_ligament_rule():
    """The integrated ligament rule: tension along the path never splits; a fully broken band splits only while
    loaded across it; an intact ligament stronger than the path-normal load holds."""
    n = 10
    normal = np.tile([0.0, 1.0], (n, 1))  # path along x, normal along y (tangent-basis a, b)
    ds = np.full(n, 6e4)
    across = np.tile([0.0, 2e11, 0.0], (n, 1))
    along = np.tile([2e11, 0.0, 0.0], (n, 1))
    broken, strength = np.zeros(n, bool), np.full(n, 1e11)
    assert not proto.ligament_splits(along, normal, ds, ~broken, strength), "tension along the path split it"
    assert proto.ligament_splits(across, normal, ds, ~broken, strength)
    assert not proto.ligament_splits(across, normal, ds, ~broken, np.full(n, 3e11))
    full = np.ones(n, bool)
    assert proto.ligament_splits(across, normal, ds, ~full, strength)
    assert not proto.ligament_splits(np.zeros((n, 3)), normal, ds, ~full, strength), "a loadless broken band split"
    assert not proto.ligament_splits(-across, normal, ds, ~full, strength), "compression across a broken band split"
    # A band split into two broken stretches with a stronger gap: integrated, not per cell.
    mixed_broken = np.zeros(n, bool)
    mixed_broken[:7] = True
    assert proto.ligament_splits(across, normal, ds, ~mixed_broken, np.full(n, 5e11))  # 3 intact cells < the load
    assert not proto.ligament_splits(across, normal, ds, ~mixed_broken, np.full(n, 1.5e12))  # ... and above it


def check_directional_sign():
    """Directional load signs: balanced forces pushing a plate's two halves together give no extension; pulling
    them apart does, and both agree with the sign of the cut resultants."""
    lat = np.linspace(-0.05, 0.05, 15)
    x, y = np.meshgrid(lat, lat, indexing="ij")
    points = np.column_stack([np.ones(x.size), x.ravel(), y.ravel()])
    points /= np.linalg.norm(points, axis=1)[:, None]
    areas = np.full(len(points), (0.0072 * mech.RADIUS_M) ** 2)
    along_y = np.tile([0.0, 0.0, 1.0], (len(points), 1))
    sign = np.where(points[:, 1] > 0, 1.0, -1.0)
    for direction, expect_extension in ((+1.0, True), (-1.0, False)):  # + : each half pulled away from the middle
        force = direction * 1e11 * sign[:, None] * np.array([0.0, 1.0, 0.0])[None, :] + 0 * along_y
        tensor = P.k_cut * proto.cut_tensor(points, areas, force, 0.0072 * mech.RADIUS_M)
        lam, lam_min, _ = proto.principal(tensor)
        middle = np.abs(points[:, 1]) < 0.01
        top, bottom = float(np.median(lam[middle])), float(np.median(lam_min[middle]))
        if expect_extension:
            assert top > 0 and top > 5 * abs(bottom), (top, bottom)
        else:
            # The six-direction fit leaks about 10% of a uniaxial load into the other principal value at every scale (probed
            # at 0.05-0.2 rad); the leak is symmetric, so a compressive plate shows extension of about a tenth of it.
            assert bottom < 0 and top < 0.15 * abs(bottom), ("compression registered as extension", top, bottom)


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
    children = np.array([int(k) for k in mech.SURFACE_FIELDS and __import__("app.sparse_quad_patch", fromlist=["x"]).child_cell_keys(np.array([key]))[0]])
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
    eps = np.full(n, 1.0)
    plate.set_fields_on_plate(rift_strain=eps, rift_warm_strain=np.full(n, 0.6), rift_broken=np.ones(n, bool))
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


def check_latch_no_eruption(seed):
    """The broken latch: overloaded cells latch at the floor with no interior eruption (no reference column, no volcano)."""
    w = world(seed)
    plate = biggest(w)
    state = step_mod._STATE
    heavy = Params(eps_dot0=5e-6, k_cap=0.05)
    old = state.params
    state.params = heavy
    try:
        volcano_before = plate.collect("is_volcano").copy()
        hc0 = plate.collect("crustal_thickness_m").copy()
        for _ in range(6):
            step_mod.damage_step(w, plate, 1e6, state)
        hc, volcano = plate.collect("crustal_thickness_m"), plate.collect("is_volcano")
        broken = plate.collect("rift_broken")
        assert broken.any(), "nothing latched under an overload"
        assert np.array_equal(volcano, volcano_before), "interior thinning ignited a volcano"
        assert hc[broken].min() >= rheology.RIFT_CRITICAL_THICKNESS_M - 1e-6
        assert not np.any(step_mod._is_reference(hc, plate.collect("mantle_lithosphere_thickness_m")) & ~step_mod._is_reference(hc0, plate.collect("mantle_lithosphere_thickness_m")) & broken)
    finally:
        state.params = old


def check_damage_runs_for_every_plate(seed):
    """The damage update on a plate after another plate splits and on a daughter in cooldown."""
    w = world(seed)
    state = step_mod._STATE
    state.pass_key = None
    young = w.plates[-1]
    young.reset_age()
    before = {p.plate_id: p.collect("rift_strain").copy() for p in w.plates}
    step_mod.maybe_split_plate(w, w.plates[0], 1e6)
    assert set(state.steps) == {p.plate_id for p in w.plates}, "a plate was skipped by the damage pass"
    assert young.plate_id in state.steps and young.age_years < 1.0
    again = len(state.steps)
    step_mod.maybe_split_plate(w, w.plates[1], 1e6)  # same step: no second pass
    assert len(state.steps) == again and state.pass_key == w.elapsed_years
    assert any(not np.array_equal(before[p.plate_id], p.collect("rift_strain")) for p in w.plates)


def _force_band(w, plate, state):
    """A fully broken strip across the plate's widest axis, with tension across it, for the split checks."""
    step = state.steps[plate.plate_id]
    pts = step.points
    centre = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - centre, full_matrices=False)
    minor = vt[1]
    across = (pts - centre) @ minor
    width = 0.8 * mech.line_spacing_rad(w.node_density) if hasattr(mech, "line_spacing_rad") else 0.012
    band = np.abs(across) < width
    broken = band.copy()
    step.broken = broken
    step.eps = np.where(band, 1.5, 0.0)
    from app import thin_shell

    first, second = thin_shell._tangent_basis(pts)
    a, b = first @ minor, second @ minor
    tension = 1.0e12
    step.tensor = np.column_stack([tension * a * a, tension * b * b, tension * a * b])
    return band


def check_edge_to_edge_band_splits(seed):
    """A fully broken edge-to-edge band splits; the path is a valid chain and the daughters are two connected plates
    that partition the parent and keep its ledgers."""
    w = world(seed)
    state = step_mod._STATE
    plate = max((p for p in w.plates if p.node_count() > 400), key=lambda p: p.node_count())
    plate._age_years = 1e9
    step_mod.damage_step(w, plate, 1e6, state)
    band = _force_band(w, plate, state)
    assert band.sum() > 5
    n_before, crust_before = plate.node_count(), float(np.dot(plate.collect("crustal_thickness_m"), plate.node_areas_m2()))
    result = step_mod.breakthrough(w, plate, state)
    assert result is not None, ("no split", state.counters)
    a, b = result
    assert a.node_count() + b.node_count() == n_before and a.plate_id != b.plate_id
    crust_after = sum(float(np.dot(d.collect("crustal_thickness_m"), d.node_areas_m2())) for d in (a, b))
    assert abs(crust_after - crust_before) <= 1e-9 * crust_before
    keys = np.concatenate([a.cell_keys, b.cell_keys])
    assert len(np.unique(keys)) == n_before, "a cell belongs to both daughters"
    for d in (a, b):
        assert np.all(d.collect("rift_strain") >= 0.0)
    event = state.events[-1]
    assert event["kind"] == "split" and event["intact_share"] == 0.0 and event["path_cells"] >= 2


CHECKS = [check_strength_bounds, check_softening_monotone, check_strain_timestep, check_healing, check_latch_stops_strain,
          check_plume_expiry, check_flank_conservation, check_flank_no_recipient_and_floor, check_ligament_rule, check_directional_sign]
PLATE_CHECKS = [check_remap_without_healing, check_reset_on_eruption, check_margin_insertion_share, check_latch_no_eruption,
                check_damage_runs_for_every_plate, check_edge_to_edge_band_splits]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--seed", type=int, default=106533960)
    parser.add_argument("--only", nargs="*")
    args = parser.parse_args()
    failed = 0
    for check in CHECKS + PLATE_CHECKS:
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
