"""Bounds on the best stationary rainfall achievable with fixed transport and bounded sinks."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.optimize import brentq, linprog
from scipy.sparse.linalg import spsolve

from rfc0019_coarse_moisture import SECONDS_PER_DAY, net_export_operator, solve_budget

RAINOUT_DAYS = (5.0, 18.0)
LAND_SHARE_TOLERANCE = 0.05
RMS_TARGET = 0.50
TANGENT_LIMIT = 4
LP_SECONDS = 15
PRIMAL_TOLERANCE = 1e-8


def _enclosing_mass(export, supply, sink):
    operator = export + sink * sparse.eye(supply.size, format="csr")
    mass = spsolve(operator.tocsc(), supply)
    residual = supply - operator @ mass
    # Inverse is nonnegative and every column of the inverse sums to 1/sink.
    # This bounds the mass error by the l1 residual, with room for matvec rounding.
    roundoff = 64 * np.finfo(float).eps * (np.abs(supply).sum() + (abs(operator) @ abs(mass)).sum())
    return mass, float((np.abs(residual).sum() + roundoff) / sink)


def rainout_relaxation_lower_bound(export, supply, reference, area):
    """Bound RMS by keeping per-cell rain intervals and global mass conservation."""
    scale = float(supply.mean())
    supply = supply / scale
    reference = reference / scale
    fast, slow = (1 / day for day in RAINOUT_DAYS)
    fast_mass, fast_pad = _enclosing_mass(export, supply, fast)
    slow_mass, slow_pad = _enclosing_mass(export, supply, slow)
    mass_lo, mass_hi = np.maximum(fast_mass - fast_pad, 0), slow_mass + slow_pad
    # Slower removal increases column mass for a nonnegative transport operator.
    # Bound rain both by sink times mass and by the stationary ledger.
    outgoing = export.diagonal()
    incoming = sparse.diags(outgoing) - export
    net_lo = supply + incoming @ mass_lo - outgoing * mass_hi
    net_hi = supply + incoming @ mass_hi - outgoing * mass_lo
    rain_lo = np.maximum(slow * mass_lo, net_lo)
    rain_hi = np.minimum(fast * mass_hi, net_hi)
    total_supply = float(supply.sum())
    if np.any(rain_lo > rain_hi) or not rain_lo.sum() <= total_supply <= rain_hi.sum():
        raise AssertionError("the relaxed rain intervals do not contain a conserved budget")

    weights = 1.0 / (area / area.mean())
    def rain_for_shift(shift):
        return np.clip(reference + shift / weights, rain_lo, rain_hi)
    lo = float(np.min((rain_lo - reference) * weights)) - 1
    hi = float(np.max((rain_hi - reference) * weights)) + 1
    shift = brentq(lambda value: rain_for_shift(value).sum() - total_supply, lo, hi)
    rain = rain_for_shift(shift)
    # Evaluate the Lagrangian minimum to get a lower bound even if the shift
    # leaves a small error in matching total supply. Ignore land share here,
    # which enlarges the feasible set and makes this a conservative bound.
    cost = float(np.dot(weights, (rain - reference)**2))
    dual = cost - 2 * shift * float(rain.sum() - total_supply)
    norm = float(np.dot(weights, reference**2))
    return float(np.sqrt(max(dual / norm - 1e-8, 0)))


def _constraints(export, source, land, share_ref, fast, slow):
    count = source.size
    eye = sparse.eye(count, format="csr")
    land_export = sparse.csr_matrix(land.astype(float)[None, :]) @ export
    land_source = source[land].sum()
    lo = max(share_ref - LAND_SHARE_TOLERANCE, 0) * source.sum()
    hi = min(share_ref + LAND_SHARE_TOLERANCE, 1) * source.sum()
    matrix = sparse.vstack((export + slow * eye, -export - fast * eye, -land_export, land_export), format="csr")
    rhs = np.concatenate((source, -source, [hi - land_source, land_source - lo]))
    return matrix, rhs, lo, hi


def _tangent_bound(export, source, reference, weights, point, matrix, rhs, mass_limit):
    residual = source - export @ point - reference
    objective = float(np.dot(weights, residual**2))
    gradient = -2 * (export.T @ (weights * residual))
    linear = linprog(
        gradient, A_ub=matrix, b_ub=rhs, bounds=(0, None), method="highs-ipm",
        options={"time_limit": LP_SECONDS, "primal_feasibility_tolerance": PRIMAL_TOLERANCE,
                 "dual_feasibility_tolerance": PRIMAL_TOLERANCE},
    )
    if not linear.success:
        return None, None, str(linear.message)

    # Any nonpositive multipliers yield a bound once we charge for the
    # smallest reduced cost. Total mass is bounded by total source / slow.
    dual = np.minimum(linear.ineqlin.marginals, 0.0)
    reduced = gradient - matrix.T @ dual
    intercept = objective - float(np.dot(gradient, point))
    correction = min(float(reduced.min()), 0.0) * mass_limit
    dual_value = float(np.dot(rhs, dual))
    rounding_margin = 1e-8 * (1 + abs(intercept) + abs(dual_value) + abs(correction))
    lower = intercept + dual_value + correction - rounding_margin
    return lower, linear.x, "optimal"


def _upper_bound(grid, faces, source, scale, ref, weights, x, implied_rain, fast, slow, lo, hi):
    # Re-solve the physical ledger with a bounded sink so the witness does
    # not depend on the LP's equality or inequality feasibility tolerance.
    implied = np.divide(implied_rain, x, out=np.full_like(x, slow), where=x > 0)
    sink = np.clip(implied, slow, fast) / SECONDS_PER_DAY
    budget = solve_budget(grid.area, source, sink, faces)
    rain = budget.rain / scale
    if not lo <= rain[~grid.is_ocean].sum() <= hi:
        return None, None
    return float(np.dot(weights, (rain - ref)**2)), sink


def bounded_rainout(grid, faces, source, reference_rain_kg_s, baseline):
    """Return a conservative lower bound and any independently re-solved upper bound."""
    fastest, slowest = RAINOUT_DAYS
    fast, slow = 1 / fastest, 1 / slowest
    export = net_export_operator(grid.area, faces) * SECONDS_PER_DAY
    scale = float(source.mean())
    supply = source / scale
    ref = reference_rain_kg_s / scale
    weights = 1.0 / grid.area
    weights /= np.dot(weights, ref**2)
    land = ~grid.is_ocean
    share_ref = float(ref[land].sum() / ref.sum())
    matrix, rhs, lo, hi = _constraints(export, supply, land, share_ref, fast, slow)
    point = baseline.mass / (scale * SECONDS_PER_DAY)
    mass_limit = float(supply.sum() / slow)
    relaxation_lower = rainout_relaxation_lower_bound(export, source, reference_rain_kg_s, grid.area)
    lower = relaxation_lower**2
    upper = None
    witness = None
    statuses = []
    for _ in range(TANGENT_LIMIT):
        candidate_upper, sink = _upper_bound(
            grid, faces, source, scale, ref, weights, point, supply - export @ point, fast, slow, lo, hi,
        )
        if candidate_upper is not None and (upper is None or candidate_upper < upper):
            upper, witness = candidate_upper, sink
        if upper is not None and upper <= RMS_TARGET**2:
            break
        if lower > RMS_TARGET**2:
            break
        bound, next_point, status = _tangent_bound(export, supply, ref, weights, point, matrix, rhs, mass_limit)
        statuses.append(status)
        if bound is None:
            break
        lower = max(lower, bound)
        if lower > RMS_TARGET**2:
            break
        delta = export @ (next_point - point)
        residual = supply - export @ point - ref
        curvature = float(np.dot(weights, delta**2))
        # If the start misses land share, move to the LP's feasible domain first.
        land_rain = (supply - export @ point)[land].sum()
        alpha = np.clip(np.dot(weights, residual * delta) / curvature, 0, 1) if curvature > 0 else 1.0
        point = point + (alpha if lo <= land_rain <= hi else 1.0) * (next_point - point)
    if upper is not None and lower > upper + 1e-8:
        raise AssertionError("rainout bounds cross beyond numerical tolerance")
    if lower > RMS_TARGET**2:
        verdict = "excluded_at_fixed_sources_and_transport"
    elif upper is not None and upper <= RMS_TARGET**2:
        verdict = "feasible_witness"
    else:
        verdict = "inconclusive"
    return {
        "verdict": verdict,
        "rms_lower_bound": float(np.sqrt(lower)),
        "rms_relaxation_lower_bound": relaxation_lower,
        "rms_feasible_upper_bound": None if upper is None else float(np.sqrt(upper)),
        "rms_bound_gap": None if upper is None else float(np.sqrt(upper) - np.sqrt(lower)),
        "land_share_bounds": [float(lo / supply.sum()), float(hi / supply.sum())],
        "tangent_lp_statuses": statuses,
    }, witness


def _archived_control_bound(path, metrics):
    """Recover directed rates on the active support of a zero-diffusivity capture."""
    if metrics["clipped_negative_mass_kg"] != 0:
        raise ValueError("archive mass was clipped, so its active support cannot be recovered")
    with np.load(path) as arrays:
        mass, source, flux = (arrays[key] for key in ("mass_kg", "source_kg_s", "face_flux_kg_s"))
        first, second = arrays["face_first"], arrays["face_second"]
        active = mass > 0
        live = flux != 0
        positive = flux[live] > 0
        donor, receiver = np.where(positive, first[live], second[live]), np.where(positive, second[live], first[live])
        # With positive sinks, exact zero stationary mass has neither a source
        # nor a directed path from any source. Its outgoing rates do not enter
        # any admissible budget. Keep these cells in the score with zero rain.
        if np.any(source[~active]) or not np.all(active[donor] & active[receiver]):
            raise ValueError("archive has unexpected transport into zero-mass cells")
        rates = np.abs(flux[live]) / mass[donor]
        outgoing = np.bincount(donor, weights=rates, minlength=mass.size)
        incoming = sparse.csr_matrix((-rates, (receiver, donor)), shape=(mass.size, mass.size))
        export = sparse.diags(outgoing) + incoming
        residual = float(np.abs(source - arrays["rain_kg_s"] - export @ mass).sum() / source.sum())
        if residual > 1e-8:
            raise AssertionError("recovered rates do not reproduce the original control ledger")
        bound = rainout_relaxation_lower_bound(
            export * SECONDS_PER_DAY, source, arrays["ref_kg_s"], arrays["area_m2"],
        )
    return {"rms_lower_bound": bound, "control_reconstruction_residual": residual,
            "excludes_fixed_transport_and_source": bound > RMS_TARGET}


def main():
    parser = argparse.ArgumentParser(description="Bound rainout RMS from zero-mixing diagnostic captures.")
    parser.add_argument("captures", nargs="+", type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    report = {}
    for capture in args.captures:
        raw = (capture / "results.json").read_bytes()
        case = json.loads(raw)
        report[capture.name] = {
            "save": case["save"], "save_sha256": case["save_sha256"],
            "results_sha256": hashlib.sha256(raw).hexdigest(),
            "winds": {
                item["wind"]: _archived_control_bound(
                    capture / f"{item['wind']}-control.npz", item["arms"]["control"],
                )
                for item in case["results"]
            },
        }
    args.out.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
