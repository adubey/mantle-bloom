"""RFC 0015 gate G1 on the march (#766): response times of each ice body under arm (m3).

Gate branch only; `replay_g0.py --g1-myr T ...` installs it. On each named step (by closing
time) it keeps the real march, then re-marches the step from the same inputs:

- `ref`: the real opening ice, logged substep by substep (it must equal the real march);
- `cap<N>`: the real opening with the longest substep capped at N yr, for N in CONVERGENCE_CAPS.
  G1's convergence check compares each cap with the next smaller one;
- `fine x0.90` .. `fine x1.10`: the opening scaled by D2's perturbations, ±5% and ±10%, marched for
  FINE_WINDOW_YR at fixed FINE_SUBSTEP_YR substeps. A perturbation decays within the march's first
  1 kyr substep, so the full march can't resolve τ.

Bodies are the connected components of the opening ice (depth > 1 m) over the hydrology
neighbour graph; those over 10^4 km^2 are measured one by one and the rest pooled. For each
body and perturbation, the gap to the unperturbed fine march is followed through the window:

- `tau_efold_yr`: when the gap first falls to 1/e of its opening value (None if it never does
  within the window);
- `closing_ratio`, the gap at the window's end over its opening value.

Alongside: `turnover_yr`, the body's opening volume over its melt rate in `ref` (the RFC's
`H/|b|`), and frictional heat's melt equivalent on warm-based nodes beside the basal melt. One JSON line per step goes to
`<capture-dir>/g1.jsonl`.
"""

from __future__ import annotations

import json
import time

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

PERTURBATIONS = (0.90, 0.95, 1.05, 1.10)
CONVERGENCE_CAPS = (5000.0, 1000.0, 500.0, 250.0)
FINE_SUBSTEP_YR = 125.0
FINE_WINDOW_YR = 10_000.0
MIN_BODY_AREA_M2 = 1e10
ICE_NODE_DEPTH_M = 1.0
SECONDS_PER_YEAR = 3.15576e7
LATENT_HEAT_J_KG = 3.34e5


def body_labels(depth: np.ndarray, neighbor_idx: np.ndarray, areas: np.ndarray) -> tuple[np.ndarray, int]:
    """Per node: 0 for no ice, 1 for small bodies pooled, 2.. for each body over
    MIN_BODY_AREA_M2, largest first. Returns the labels and the number of measured bodies."""
    n, k = neighbor_idx.shape
    ice = depth > ICE_NODE_DEPTH_M
    rows = np.repeat(np.arange(n), k)
    cols = neighbor_idx.ravel()
    keep = (cols >= 0) & (cols < n)
    rows, cols = rows[keep], cols[keep]
    keep = ice[rows] & ice[cols]
    graph = coo_matrix((np.ones(int(keep.sum())), (rows[keep], cols[keep])), shape=(n, n))
    _, component = connected_components(graph, directed=False)
    labels = np.zeros(n, dtype=np.int64)
    component_area = np.bincount(component[ice], weights=areas[ice])
    big = [c for c in np.argsort(-component_area) if component_area[c] >= MIN_BODY_AREA_M2]
    labels[ice] = 1
    for rank, c in enumerate(big):
        labels[ice & (component == c)] = rank + 2
    return labels, len(big)


def _gap_times(times_ref, path_ref, times_cell, path_cell, opening_gap) -> dict:
    """How a body's gap to the unperturbed march decays through the window."""
    grid = np.union1d(times_ref, times_cell)
    gap = np.interp(grid, times_cell, path_cell) - np.interp(grid, times_ref, path_ref)
    if opening_gap == 0.0:
        return {"tau_efold_yr": None, "closing_ratio": None}
    ratio = gap / opening_gap
    below = np.flatnonzero(ratio <= np.exp(-1.0))
    return {"tau_efold_yr": float(grid[below[0]]) if len(below) else None, "closing_ratio": round(float(ratio[-1]), 5)}


def install(args, world) -> None:  # noqa: ARG001 - matches replay_g0's installers
    from app.hydroclimate import hydrology, ice_sheet  # noqa: E402

    real_march = hydrology._g0_march_ice
    wanted = [round(t * 1e6) for t in args.g1_myr]
    args.capture_dir.mkdir(parents=True, exist_ok=True)
    out = (args.capture_dir / "g1.jsonl").open("a")

    def logged_march(march_args, opening, longest, labels, areas, bins, fixed=None, years=None):
        """One march from `opening`, recording each body's volume after every substep. `fixed`
        holds every substep at that length; `years` shortens the march."""
        plain = ice_sheet.step_ice
        g0 = ice_sheet.G0
        kept = (g0.march_max_substep_years, g0.march_min_substep_years, g0.march_substep_years)
        times, volumes = [0.0], [np.bincount(labels, weights=opening * areas, minlength=bins)]

        def recording(*a, **k):
            new_depth, melt, ice = plain(*a, **k)
            times.append(times[-1] + a[-1])
            volumes.append(np.bincount(labels, weights=new_depth * areas, minlength=bins))
            return new_depth, melt, ice

        ice_sheet.step_ice = recording
        g0.march_max_substep_years = longest if fixed is None else fixed
        if fixed is not None:
            g0.march_min_substep_years = g0.march_substep_years = fixed
        try:
            depth, melt, ice = real_march(*march_args(opening, years))
        finally:
            ice_sheet.step_ice = plain
            g0.march_max_substep_years, g0.march_min_substep_years, g0.march_substep_years = kept
        return depth, melt, ice, np.array(times), np.array(volumes)

    def wrapped(world_, points, neighbor_idx, elevation, deflection, areas, prev_depth, precipitation,
                temperature, amplitude, heat, mobile_cover, is_ocean, allow_ocean_ice, plates_in_order, years):
        heat_snapshot = None if heat is None else heat.copy()
        depth, melt, ice = real_march(world_, points, neighbor_idx, elevation, deflection, areas, prev_depth, precipitation,
                                      temperature, amplitude, heat, mobile_cover, is_ocean, allow_ocean_ice,
                                      plates_in_order, years)
        target = next((t for t in wanted if abs(world_.elapsed_years - t) < 1000), None)
        if target is None:
            return depth, melt, ice
        started = time.perf_counter()
        kept_log = len(ice_sheet.G0_MARCH_LOG)
        kept_heat = ice_sheet.G0.frictional_heat_w_m2
        opening = np.clip(prev_depth, 0.0, None)
        labels, measured = body_labels(opening, neighbor_idx, areas)
        bins = measured + 2

        def march_args(start, span=None):
            return (world_, points, neighbor_idx, elevation, deflection, areas, start.copy(), precipitation, temperature,
                    amplitude, None if heat_snapshot is None else heat_snapshot.copy(), mobile_cover, is_ocean,
                    allow_ocean_ice, plates_in_order, years if span is None else span)

        longest = ice_sheet.G0.march_max_substep_years
        ref_depth, ref_melt, ref_ice, ref_times, ref_path = logged_march(march_args, opening, longest, labels, areas, bins)
        caps = {cap: logged_march(march_args, opening, cap, labels, areas, bins) for cap in CONVERGENCE_CAPS if cap < longest}
        fine = {f: logged_march(march_args, opening * f, longest, labels, areas, bins, fixed=FINE_SUBSTEP_YR, years=FINE_WINDOW_YR)
                for f in (1.0, *PERTURBATIONS)}
        del ice_sheet.G0_MARCH_LOG[kept_log:]
        ice_sheet.G0.frictional_heat_w_m2 = kept_heat
        ice_sheet.G0.reference_depth = None

        melt_volume = np.bincount(labels, weights=np.clip(ref_melt, 0.0, None) * areas, minlength=bins)
        bodies = []
        warm = getattr(ref_ice, "warm_based", None)
        for b in range(1, bins):
            in_body = labels == b
            if not in_body.any():
                continue
            v0 = float(ref_path[0][b])
            closes = {longest: float(ref_path[-1][b]), **{cap: float(cell[4][-1][b]) for cap, cell in caps.items()}}
            ladder = sorted(closes, reverse=True)
            row = {
                "body": "small" if b == 1 else b - 1,
                "area_km2": round(float(areas[in_body].sum()) / 1e6),
                "opening_km3": round(v0 / 1e9, 1),
                "closing_km3": round(float(ref_path[-1][b]) / 1e9, 1),
                "lat": round(float(np.degrees(np.arcsin(np.mean(points[in_body, 2]) / np.mean(np.linalg.norm(points[in_body], axis=1))))), 1),
                "turnover_yr": round(v0 / (melt_volume[b] / years)) if melt_volume[b] > 0 else None,
                "closing_by_cap_km3": {str(int(cap)): round(closes[cap] / 1e9, 1) for cap in ladder},
                "convergence_gap": {f"{int(a)}->{int(c)}": round(closes[c] / max(closes[a], 1.0) - 1.0, 5) for a, c in zip(ladder, ladder[1:])},
                "real_matches_ref": bool(abs(float(np.dot(depth, areas * in_body)) - float(ref_path[-1][b])) <= 1e-6 * max(v0, 1.0)),
                "fine": {f"x{f:.2f}": _gap_times(fine[1.0][3], fine[1.0][4][:, b], fine[f][3], fine[f][4][:, b], v0 * (f - 1.0))
                         for f in PERTURBATIONS},
            }
            bodies.append(row)
        frictional = None
        if heat_snapshot is not None:
            on_warm = areas if warm is None else areas * warm
            frictional = float(np.dot(np.clip(heat_snapshot, 0.0, None), on_warm)) * years * SECONDS_PER_YEAR / (ice_sheet.RHO_ICE * LATENT_HEAT_J_KG)
        ledger = ref_ice.ledger
        record = {
            "target_myr": target / 1e6, "step_years": years, "longest_substep_yr": longest,
            "ref_substeps": len(ref_times) - 1, "cap_substeps": {str(int(c)): len(cell[3]) - 1 for c, cell in caps.items()},
            "world_opening_km3": round(float(np.dot(opening, areas)) / 1e9, 1),
            "real_closing_km3": round(float(np.dot(depth, areas)) / 1e9, 1),
            "frictional_melt_km3": None if frictional is None else round(frictional / 1e9, 1),
            "melt_km3": round((ledger.surface_melt_m3 + ledger.basal_melt_m3) / 1e9, 1),
            "basal_melt_km3": round(ledger.basal_melt_m3 / 1e9, 1),
            "probe_seconds": round(time.perf_counter() - started, 1),
            "bodies": bodies,
        }
        out.write(json.dumps(record) + "\n")
        out.flush()
        print(json.dumps({k: v for k, v in record.items() if k != "bodies"}), flush=True)
        return depth, melt, ice

    hydrology._g0_march_ice = wrapped
