"""RFC 0018 D1's integration and propagation checks on manufactured plates (#888), before any parameter search.

Usage (repo root, this worktree's venv; no world is generated, a few minutes in all):
    backend/.venv/bin/python bin/debug/check_rift_line_d1.py [--out analysis/issue888/integration-propagation.json]

*Integration check.* A manufactured plate with a fixed pull and a fixed straight line (weak belt, with a weaker
volcano cell at its middle so two tips run) is integrated for 60 Myr with steps of 0.01 Myr (the reference), 0.5, 1,
2 and 5 Myr, once at the starting rate and once with the pull at the ``r_max`` cap (the case needing the most
sub-steps). After each step the thinning feeds back through ``N_cap`` (the operator split). Pass when, at every
common output time (every 10 Myr), accumulated ``rift_strain`` on the band is within 0.05 of the reference, the broken
sets differ by at most one cell at each tip, and the split time is within one coarse step of the reference.

*Propagation check.* A straight line at constant ``Λ`` (the amplitude ``A`` fixed at ``ε̇_0``, i.e. ``Λ = 2 Λ_grow``)
is run at the default and at half node density (spacing ``h`` and ``√2 h``) for ``ℓ_tip`` of 60, 125 and 250 km.
Nucleation time, first advance and steady tip speed are measured and set against the RFC's discrete law
(``ε_cand / A``, ``ε_cand e^x / A``, ``v = (A ℓ_tip / ε_cand) x / (e^x - 1)`` with ``x = h / ℓ_tip``). It is run
with ``Δε_max`` 0.05 (the RFC's value) and 0.005 to separate the law from the sub-step resolution. Pass when the
measured steady speed is within 20% of the law.
"""

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix

sys.path[:0] = [str(Path(__file__).resolve().parent)]

import rift_damage_mechanism as mech  # noqa: E402
import rift_line as rl  # noqa: E402
from rift_damage_mechanism import Params, RADIUS_M  # noqa: E402

from app import lithosphere, rheology  # noqa: E402

MYR = 1e6
NUCLEUS_FACTOR = 0.65  # extra weakness of the middle cell (a volcano) so the nucleus is there
OUTPUT_EVERY_MYR = 10.0  # common to every tested step length (0.5, 1, 2, 5)


# --- a manufactured plate ------------------------------------------------------------------------


def grid_edges(nx: int, ny: int):
    """Both directions of every edge of the 4-connected nx x ny grid, as (rows, cols) index arrays."""
    index = np.arange(nx * ny).reshape(ny, nx)
    rows, cols = [], []
    for a, b in ((index[:, :-1], index[:, 1:]), (index[:-1, :], index[1:, :])):
        rows += [a.ravel(), b.ravel()]
        cols += [b.ravel(), a.ravel()]
    return np.concatenate(rows), np.concatenate(cols)


class Toy:
    """A flat lattice on a sphere patch with one straight line along its middle row."""

    def __init__(self, h_m: float, nx: int, ny: int, p: Params, pull_ratio: float, n_nuclei: int = 1):
        self.p, self.h, self.nx, self.ny = p, h_m, nx, ny
        u = (np.arange(nx) - (nx - 1) / 2) * h_m
        v = (np.arange(ny) - (ny - 1) / 2) * h_m
        uu, vv = np.meshgrid(u, v)
        lon, lat = uu.ravel() / RADIUS_M, vv.ravel() / RADIUS_M  # equirectangular: exact spacing along the line
        self.points = np.column_stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
        n = len(self.points)
        index = np.arange(n).reshape(ny, nx)
        r, c = grid_edges(nx, ny)
        length = RADIUS_M * np.arccos(np.clip(np.sum(self.points[r] * self.points[c], axis=1), -1.0, 1.0))
        self.graph = csr_matrix((length, (r, c)), shape=(n, n))
        self.areas = np.full(n, h_m * h_m)
        self.chain = index[ny // 2]
        self.disk = mech.disk_average_matrix(self.points, self.areas, p.w0_m / 2.0)
        self.hc = np.full(n, 35e3)
        self.hm = np.full(n, 90e3)
        self.lag = np.zeros(n)
        mult = np.full(n, mech.MULT_BELT)
        self.nuclei = np.linspace(0, nx - 1, n_nuclei + 2)[1:-1].round().astype(int)
        mult[self.chain[self.nuclei]] *= NUCLEUS_FACTOR
        self.mult = mult
        self.eps, self.warm = np.zeros(n), np.zeros(n)
        self.broken = np.zeros(n, dtype=bool)
        base0 = self.base()[self.chain]
        # path-normal load giving the requested initial Λ on the line (W = 1): Λ0 = σ / mean(base)
        self.sigma = pull_ratio * p.lam_grow * float(np.sum(base0 * h_m) / (h_m * len(self.chain)))
        self.geometry = rl.path_geometry(self.points, self.chain, *_basis(self.points))
        self.t = 0.0
        self.max_substeps = 0

    def base(self) -> np.ndarray:
        cap = mech.capacity(np.zeros(len(self.hc), dtype=int), self.hc, self.hm, self.lag)
        return cap * self.p.k_cap * self.mult

    def problem(self) -> rl.LineProblem:
        p, chain = self.p, self.chain
        band, rb, rline = rl.band_reach(self.points, chain, self.eps, p)
        floor = rheology.RIFT_CRITICAL_THICKNESS_M
        latch = np.minimum(p.eps_break, self.eps + np.log(np.maximum(self.hc, floor) / floor))
        latch = np.where(self.broken, self.eps, np.maximum(latch, self.eps))
        base = self.base()
        return rl.LineProblem(
            chain=chain, s=self.geometry.s, ds=self.geometry.ds, normal_load=np.full(len(chain), self.sigma),
            base_quiet=base[chain], base_plume=base[chain], plume_left_yr=np.zeros(len(chain)), intact=~self.broken[chain],
            disk_rows=self.disk[chain], band=band, reach_band=rb, reach_line=rline, latch=latch,
            capacity=rl.CapacityInputs(np.zeros(len(chain), dtype=int), self.hc[chain], self.hm[chain], self.lag[chain],
                                       np.full(len(chain), p.k_cap), self.mult[chain], self.mult[chain], self.eps[chain]))

    def step(self, years: float, fixed_amplitude: float | None = None, thin: bool = True) -> rl.Advance:
        prob = self.problem()
        out = rl.advance_rift(prob, self.eps, self.warm, years, self.p, fixed_amplitude)
        self.max_substeps = max(self.max_substeps, out.substeps)
        d_eps = np.maximum(out.eps - self.eps, 0.0)
        self.eps, self.warm = out.eps, np.minimum(out.warm, out.eps)
        if thin:  # the operator split: thinning changes N_cap at the next step
            floor_c, floor_m = rheology.RIFT_CRITICAL_THICKNESS_M, lithosphere.MIN_MANTLE_LITHOSPHERE_THICKNESS_M
            self.hc = np.where(d_eps > 0, np.maximum(self.hc * np.exp(-d_eps), np.minimum(self.hc, floor_c)), self.hc)
            self.hm = np.where(d_eps > 0, np.maximum(self.hm * np.exp(-d_eps), np.minimum(self.hm, floor_m)), self.hm)
        newly = (d_eps > 0) & ~self.broken & ((self.eps >= self.p.eps_break - 1e-9) | (self.hc <= rheology.RIFT_CRITICAL_THICKNESS_M + 1e-6))
        self.broken |= newly
        self.t += years
        return out

    def splits(self) -> bool:
        """Step 6's test on the line: ``P >= S`` with ``P > 0`` and at least ``f_break`` of its length broken."""
        prob = self.problem()
        n_eff, load, _ = rl.line_state(prob, self.eps, self.warm, 0.0, self.p)
        pull = float(np.sum(load * prob.ds))
        strength = float(np.sum((n_eff * prob.ds)[prob.intact]))
        share = float(prob.ds[self.broken[self.chain]].sum() / prob.ds.sum())
        return rl.ligament_gives(pull, strength) and share >= self.p.f_break


def _basis(points):
    from app import thin_shell

    return thin_shell._tangent_basis(points)


# --- the integration check -----------------------------------------------------------------------


def integrate(step_myr: float, p: Params, pull_ratio: float, total_myr: float = 60.0, h_m: float = 125e3, nx: int = 41):
    toy = Toy(h_m, nx, 9, p, pull_ratio)
    outputs, split_at = {}, None
    steps = round(total_myr / step_myr)
    next_out = OUTPUT_EVERY_MYR
    for k in range(1, steps + 1):
        toy.step(step_myr * MYR)
        myr = k * step_myr
        if split_at is None and toy.splits():
            split_at = myr
        if abs(myr - next_out) < 1e-6:
            outputs[next_out] = (toy.eps[toy.chain].copy(), toy.broken[toy.chain].copy())
            next_out += OUTPUT_EVERY_MYR
    return outputs, split_at, toy


def integration_check(p: Params) -> dict:
    result = {}
    slow = replace(p, eps_dot0=p.eps_dot0 / 8.0)  # the capped pull with a rate low enough to leave several 10 Myr outputs to compare
    for name, pull_ratio, q in (("start_rate", 2.0, p), ("r_max_cap", p.r_max * 1.2, p), ("r_max_cap_slow", p.r_max * 1.2, slow)):
        reference, ref_split, toy = integrate(0.01, q, pull_ratio)
        centre = int(toy.nuclei[0])
        cases = {}
        for step in (0.5, 1.0, 2.0, 5.0):
            out, split, run = integrate(step, q, pull_ratio)
            worst, tips_ok, rows = 0.0, True, []
            for t, (eps, broken) in out.items():
                ref_eps, ref_broken = reference[t]
                err = float(np.max(np.abs(eps - ref_eps)))
                diff = np.flatnonzero(broken != ref_broken)
                per_tip = (int((diff < centre).sum()), int((diff > centre).sum()), int((diff == centre).sum()))
                worst = max(worst, err)
                tips_ok &= max(per_tip[0], per_tip[1]) <= 1 and per_tip[2] == 0
                rows.append({"t": t, "max_eps_error": round(err, 4), "broken_diff_cells_per_tip": per_tip[:2]})
            time_ok = (split is None) == (ref_split is None) and (split is None or abs(split - ref_split) <= step + 1e-9)
            cases[str(step)] = {"worst_eps_error": round(worst, 4), "eps_ok": worst <= 0.05, "tips_ok": bool(tips_ok), "split_myr": split,
                                "split_ok": bool(time_ok), "max_substeps_per_step": run.max_substeps, "pass": bool(worst <= 0.05 and tips_ok and time_ok), "rows": rows}
        result[name] = {"reference_split_myr": ref_split, "steps": cases}
    result["pass"] = all(case["pass"] for r in result.values() if isinstance(r, dict) for case in r["steps"].values())
    return result


# --- the propagation check -----------------------------------------------------------------------


def propagation(h_m: float, l_tip_m: float, d_eps_max: float, p: Params, nx: int = 31, step_myr: float = 1.0):
    q = replace(p, l_tip_m=l_tip_m, d_eps_max=d_eps_max)
    toy = Toy(h_m, nx, 3, q, 2.0)
    toy.mult[:] = mech.MULT_BELT  # uniform line: the nucleus is the first cell; measure one front
    toy.mult[toy.chain[0]] *= NUCLEUS_FACTOR
    toy.nuclei = np.array([0])
    amplitude = q.eps_dot0
    xs = h_m / l_tip_m
    law_v = (amplitude * l_tip_m * MYR / q.eps_cand) * xs / (np.exp(xs) - 1.0)  # m / Myr
    horizon_myr = 1.2 * nx * h_m / law_v + 10.0  # long enough for the front to cross the line at the law's speed
    cross = np.full(nx, np.nan)
    previous = toy.eps[toy.chain].copy()
    for k in range(1, round(horizon_myr / step_myr) + 1):
        toy.step(step_myr * MYR, fixed_amplitude=amplitude, thin=False)
        current = toy.eps[toy.chain]
        hit = np.isnan(cross) & (current >= q.eps_cand)
        for i in np.flatnonzero(hit):  # linear interpolation of the crossing inside the step
            frac = (q.eps_cand - previous[i]) / max(current[i] - previous[i], 1e-12)
            cross[i] = (k - 1 + frac) * step_myr
        previous = current.copy()
        if not np.isnan(cross).any():
            break
    x = h_m / l_tip_m
    law_start = q.eps_cand / amplitude / MYR
    law_first = q.eps_cand * np.exp(x) / amplitude / MYR
    law_speed = (amplitude * l_tip_m * MYR / q.eps_cand) * x / (np.exp(x) - 1.0) / 1e3  # km / Myr
    cells = np.arange(nx)
    ok = ~np.isnan(cross) & (cells >= 6) & (cells <= nx - 4)
    speed = None
    if ok.sum() >= 10:
        slope = np.polyfit(cross[ok], cells[ok] * h_m / 1e3, 1)[0]
        speed = float(slope)
    return {
        "h_km": h_m / 1e3, "l_tip_km": l_tip_m / 1e3, "d_eps_max": d_eps_max, "x": float(x),
        "nucleation_myr": None if np.isnan(cross[0]) else float(cross[0]), "law_nucleation_myr": float(law_start),
        "first_advance_myr": None if np.isnan(cross[1]) else float(cross[1]), "law_first_advance_myr": float(law_first),
        "speed_km_per_myr": speed, "law_speed_km_per_myr": float(law_speed),
        "speed_ratio": None if speed is None else float(speed / law_speed),
        "within_20pct": bool(speed is not None and abs(speed / law_speed - 1.0) <= 0.2),
    }


def propagation_check(p: Params) -> dict:
    rows = []
    for density, h in (("default", 125e3), ("half", 125e3 * np.sqrt(2.0))):
        for l_tip in (60e3, 125e3, 250e3):
            for d_max in (0.05, 0.005):
                row = propagation(h, l_tip, d_max, p)
                rows.append({"density": density, **row})
                print(json.dumps(rows[-1]), flush=True)
    by = {(r["density"], r["l_tip_km"], r["d_eps_max"]): r for r in rows}
    # the falsification test: half-density speed differs from the default's by more than 25% at every tested l_tip
    grid = {}
    for l_tip in (60.0, 125.0, 250.0):
        a, b = by[("default", l_tip, 0.05)]["speed_km_per_myr"], by[("half", l_tip, 0.05)]["speed_km_per_myr"]
        grid[str(l_tip)] = None if a is None or b is None else float(b / a)
    return {"rows": rows, "half_over_default_speed": grid,
            "grid_artifact": all(v is not None and abs(v - 1.0) > 0.25 for v in grid.values()),
            "pass": all(r["within_20pct"] for r in rows if r["d_eps_max"] == 0.05)}


def summarize(check: dict) -> str:
    parts = []
    for name in ("start_rate", "r_max_cap", "r_max_cap_slow"):
        for step, case in check[name]["steps"].items():
            parts.append(f"{name[:5]} dt={step}: eps_err={case['worst_eps_error']} tips={case['tips_ok']} split={case['split_myr']}/{check[name]['reference_split_myr']} sub<={case['max_substeps_per_step']} {'PASS' if case['pass'] else 'fail'}")
    return " | ".join(parts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", type=Path)
    parser.add_argument("--only", choices=("integration", "propagation"))
    args = parser.parse_args()
    p = Params()
    result = {"params": {k: v for k, v in p.__dict__.items()}, "cell_spacing_km_default": 125.0}
    if args.only != "propagation":
        result["integration"] = {}
        for name, q in (("rfc", p), ("crossing_stops", replace(p, crossing_stops=True)),
                        ("live_ncap+crossing_stops", replace(p, live_ncap=True, crossing_stops=True)),
                        ("live_ncap+crossing_stops+d_eps_max_0.02", replace(p, live_ncap=True, crossing_stops=True, d_eps_max=0.02)),
                        ("live_ncap+crossing_stops+d_eps_max_0.01", replace(p, live_ncap=True, crossing_stops=True, d_eps_max=0.01))):
            result["integration"][name] = integration_check(q)
            print(name, summarize(result["integration"][name]), flush=True)
    if args.only != "integration":
        result["propagation"] = propagation_check(p)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
