#!/usr/bin/env python3
"""How Earth-like is a world's relief? Initial-terrain metrics against Earth (#218).

For each generated seed (or a loaded save) it reports, per node (nodes are near-equal area):
- land hypsometry: percentiles, and the share of land above 1 km and 2 km;
- ocean: percentiles and IQR, the share shallower than 3 km and than 200 m, and the share
  of 0-200 m water within 300 km of land (shelves fringe coasts);
- the share of ocean not connected to the world ocean (enclosed seas), and of land under
  lakes;
- drainage (#218): the share of land whose routed water ends somewhere other than the ocean,
  split into ends at a frozen node (rivers stop there) and unfrozen ends, and the share of
  land in closed depressions deeper than 1-1,000 m (the rise water needs to reach the ocean);
- relief: RMS elevation difference between land points at 50-1,600 km, and the same for
  "lowland" pairs with both points below 1 km (how hilly the plains are);
- column backing: the largest |elevation - isostatic_elevation(Hc, Hm)| on any node with no
  ice, columns outside the Hc/Hm caps, and the share of land that starts collapse-mobile
  (Hc > orogeny.COLLAPSE_ONSET_HC_M) or with a Moho at or above orogeny.MELT_ONSET_MOHO_C.

Earth references (GMT earth_relief_10m, SRTM15+, area weighted; lowland relief box-averaged
to ~55-90 km to match ~60 km nodes) are in EARTH. `--step-myr` steps each world and
measures again at each listed time (in `--chunk-myr` steps), to follow how the generated
relief evolves. `--png`
also renders the elevation view of each measurement.

Usage:
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_initial_relief.py \
        --seeds 447276134,373936961 --cratons-first --synthetic-relief --out /tmp/relief.json
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_initial_relief.py \
        --seeds 447276134 --cratons-first --synthetic-relief --step-myr 1,5,20 --png /tmp/maps \
        --out /tmp/relief-steps.json
    PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_initial_relief.py \
        --seeds "" --load some.mbworld --out /tmp/relief-save.json
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

R_KM = 6371.0
PERCENTILES = (5, 25, 50, 75, 90, 95, 99)
LAGS_KM = (50, 100, 200, 400, 800, 1600)

EARTH = {
    "land_fraction": 0.291,
    "land_pct_m": {"5": 36, "25": 178, "50": 423, "75": 1027, "90": 2096, "95": 3002, "99": 4608},
    "ocean_pct_m": {"5": -5594, "25": -4871, "50": -4100, "75": -3053, "90": -618, "95": -100, "99": -18},
    "land_share_above_1km": 0.257,
    "land_share_above_2km": 0.107,
    "ocean_share_above_3km_depth": 0.243,
    "ocean_share_above_200m_depth": 0.069,
    "lake_share_of_land": 0.02,
    # Endorheic drainage, ~18-20% of land (literature). Closed depressions from a priority
    # flood of earth_relief_10m box-averaged to ~55 km (3x3 cells), ice-free land, straits
    # kept open by finding the world ocean at full resolution (bin/debug/measure_earth_depressions.py).
    "endorheic_share_of_land": 0.19,
    "depression_share_gt_m": {"1": 0.100, "10": 0.088, "50": 0.057, "100": 0.040, "250": 0.015, "500": 0.006, "1000": 0.0007},
    # Box-averaged to ~55 km (3x3 cells) and ~90 km (5x5): Earth's lowland relief at node scale.
    "lowland_rms_diff_m_by_km": {"50": 131, "100": 173, "200": 214, "400": 257},
    "lowland_rms_diff_m_by_km_box90km": {"50": 119, "100": 166, "200": 210, "400": 253},
}
MIN_LAND_NODES_FOR_RELIEF = 200


def _percentiles(values: np.ndarray) -> dict[str, float]:
    if not len(values):
        return {}
    return dict(zip(map(str, PERCENTILES), np.percentile(values, PERCENTILES).round(0).tolist()))


def _hypsometry(h: np.ndarray, land: np.ndarray) -> dict:
    ocean = ~land
    out = {
        "land_fraction": float(land.mean()),
        "land_pct_m": _percentiles(h[land]),
        "ocean_pct_m": _percentiles(h[ocean]),
    }
    if land.any():
        out["land_share_above_1km"] = float((h[land] > 1000).mean())
        out["land_share_above_2km"] = float((h[land] > 2000).mean())
    if ocean.any():
        out["ocean_iqr_m"] = float(np.subtract(*np.percentile(h[ocean], [75, 25])))
        out["ocean_share_above_3km_depth"] = float((h[ocean] > -3000).mean())
        out["ocean_share_above_200m_depth"] = float((h[ocean] > -200).mean())
    return out


def _shelf_and_seas(points: np.ndarray, h: np.ndarray, land: np.ndarray) -> dict:
    ocean = ~land
    out = {}
    shallow = ocean & (h > -200)
    if shallow.any() and land.any():
        d, _ = cKDTree(points[land]).query(points[shallow], workers=-1)
        out["shelf_share_within_300km_of_land"] = float((d < 2 * np.sin(300 / R_KM / 2)).mean())
    if ocean.sum() > 2:
        sea = points[ocean]
        tree = cKDTree(sea)
        spacing = float(np.median(tree.query(sea[:: max(1, len(sea) // 2000)], k=2)[0][:, 1]))
        pairs = tree.query_pairs(1.6 * spacing, output_type="ndarray")
        graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(sea), len(sea)))
        _, labels = connected_components(graph, directed=False)
        out["enclosed_sea_share_of_ocean"] = float(1 - np.bincount(labels).max() / len(sea))
    return out


def _relief_by_lag(points: np.ndarray, h: np.ndarray, land: np.ndarray) -> dict:
    """RMS height difference between land node pairs at each lag (+-10%), all pairs and pairs
    with both ends below 1 km, from up to 4,000 sampled nodes and up to four partners each,
    drawn at random from the ring so no direction is favoured."""
    lp, lh = points[land], h[land]
    if len(lp) < MIN_LAND_NODES_FOR_RELIEF:
        return {"relief_skipped": f"fewer than {MIN_LAND_NODES_FOR_RELIEF} land nodes"}
    rng = np.random.default_rng(0)
    sel = rng.choice(len(lp), min(4000, len(lp)), replace=False)
    tree = cKDTree(lp)
    all_pairs, lowland = {}, {}
    for lag in LAGS_KM:
        chord = 2 * np.sin(lag / R_KM / 2)
        a, b = [], []
        for i, js in zip(sel, tree.query_ball_point(lp[sel], chord * 1.1)):
            js = np.asarray(js, dtype=int)
            ring = js[np.linalg.norm(lp[js] - lp[i], axis=1) > chord * 0.9]
            far = rng.choice(ring, min(4, len(ring)), replace=False) if len(ring) else ring
            a.append(np.full(len(far), lh[i]))
            b.append(lh[far])
        a, b = np.concatenate(a), np.concatenate(b)
        if len(a):
            all_pairs[str(lag)] = float(np.sqrt(np.mean((a - b) ** 2)).round(0))
            low = (a < 1000) & (b < 1000)
            if low.any():
                lowland[str(lag)] = float(np.sqrt(np.mean((a[low] - b[low]) ** 2)).round(0))
    return {"land_rms_diff_m_by_km": all_pairs, "lowland_rms_diff_m_by_km": lowland}


def _drainage(world, land: np.ndarray) -> dict:
    """Where routed water ends, from the world's hydrology cache, and closed-depression depth."""
    import heapq

    from app.hydroclimate.erosion import climate_grid_indices

    fields, climate = world.hydrology_cache, world.climate_cache
    if fields is None or climate is None or len(fields.flow_target) != len(land):
        return {"drainage_skipped": "no hydrology cache"}
    end = np.where(fields.flow_target >= 0, fields.flow_target, np.arange(len(land)))
    for _ in range(64):
        end, previous = end[end], end
        if np.array_equal(end, previous):
            break
    rows, cols = climate_grid_indices(fields.points, *climate.air_temperature_c.shape)
    frozen_end = climate.air_temperature_c[rows, cols][end] < 0.0
    closed = land & ~fields.is_ocean[end]
    # Priority flood from the world ocean over the hydrology graph: the level water must rise
    # to before it can reach the sea. Water reaches node i from each node j listing i among
    # its neighbours, so the flood walks those edges in reverse.
    n = len(land)
    k = fields.neighbor_idx.shape[1]
    upstream = csr_matrix((np.ones(n * k), (fields.neighbor_idx.ravel(), np.repeat(np.arange(n), k))), shape=(n, n))
    upstream.sort_indices()
    indptr, indices = upstream.indptr.tolist(), upstream.indices.tolist()
    elevation = fields.elevation.tolist()
    filled = np.where(fields.is_ocean, fields.elevation, np.inf).tolist()
    heap = [(filled[i], int(i)) for i in np.flatnonzero(fields.is_ocean)]
    heapq.heapify(heap)
    while heap:
        level, i = heapq.heappop(heap)
        if level > filled[i]:
            continue
        for j in indices[indptr[i] : indptr[i + 1]]:
            rise = max(level, elevation[j])
            if rise < filled[j]:
                filled[j] = rise
                heapq.heappush(heap, (rise, j))
    depth = (np.asarray(filled) - fields.elevation)[land]
    n_land = max(int(land.sum()), 1)
    return {
        "closed_drainage_share_of_land": float(closed.sum() / n_land),
        "closed_drainage_frozen_end_share_of_land": float((closed & frozen_end).sum() / n_land),
        "closed_drainage_unfrozen_share_of_land": float((closed & ~frozen_end).sum() / n_land),
        "depression_share_gt_m": {str(t): float((depth > t).sum() / n_land) for t in (1, 10, 50, 100, 250, 500, 1000)},
    }


def _backing(world, land: np.ndarray) -> dict:
    from app import lithosphere, orogeny, plates

    hc = plates.collect_all_crustal_thickness(world.plates)
    hm = plates.collect_all_mantle_lithosphere_thickness(world.plates)
    elevation = plates.collect_all_elevation(world.plates)
    ice = plates.collect_all_glacier_depth(world.plates) > 0.0
    rho = np.concatenate(
        [lithosphere.node_crust_density(p.collect("crust_type_code"), p.crust_type) for p in world.plates]
    )
    residual = np.abs(elevation - lithosphere.isostatic_elevation(hc, hm, rho))
    lag = np.concatenate([p.collect("moho_thermal_lag_c") for p in world.plates])
    moho = orogeny.moho_temperature_c(hc, hm, lag)
    return {
        "max_unbacked_m_no_ice": float(residual[~ice].max()) if (~ice).any() else 0.0,
        "share_unbacked_over_1m_no_ice": float((residual[~ice] > 1.0).mean()) if (~ice).any() else 0.0,
        "columns_outside_hc_caps": int(((hc < lithosphere.MIN_CRUSTAL_THICKNESS_M - 1e-6) | (hc > lithosphere.MAX_CRUSTAL_THICKNESS_M + 1e-6)).sum()),
        "land_collapse_mobile_share": float((hc[land] > orogeny.COLLAPSE_ONSET_HC_M).mean()) if land.any() else 0.0,
        "land_moho_melt_eligible_share": float((moho[land] >= orogeny.MELT_ONSET_MOHO_C).mean()) if land.any() else 0.0,
        "land_hc_p50_p99_km": [round(float(np.percentile(hc[land], q)) / 1e3, 1) for q in (50, 99)] if land.any() else [],
    }


def measure(world) -> dict:
    from app import plates

    points = plates.collect_all_points(world.plates)[0]
    h = plates.collect_all_elevation(world.plates) - world.sea_level_m
    land = h > 0.0
    lakes = plates.collect_all_lake_depth(world.plates) > 0.0
    out = {"elapsed_years": float(world.elapsed_years), "nodes": int(len(h)), "sea_level_m": float(world.sea_level_m)}
    out.update(_hypsometry(h, land))
    out.update(_shelf_and_seas(points, h, land))
    out["lake_share_of_land"] = float((lakes & land).sum() / max(int(land.sum()), 1))
    out.update(_relief_by_lag(points, h, land))
    out["drainage"] = _drainage(world, land)
    out["backing"] = _backing(world, land)
    return out


def _render(world, path: Path) -> None:
    from app import render_image

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(render_image.render_png(world, "eckert4", "elevation", 1400, 700))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seeds", default="447276134,373936961,910211954,374464662")
    parser.add_argument("--load", action="append", default=[], help="a .mbworld save to measure (repeatable)")
    parser.add_argument("--cratons-first", action="store_true")
    parser.add_argument("--synthetic-relief", action="store_true")
    parser.add_argument("--node-density", type=float, default=4.0)
    parser.add_argument("--land-fraction", type=float, default=0.29)
    parser.add_argument("--continental-fraction", type=float, default=None)
    parser.add_argument("--step-myr", default="", help="comma-separated times (Myr) to step to and measure again")
    parser.add_argument("--chunk-myr", type=float, default=1.0, help="step size used to reach each --step-myr time")
    parser.add_argument("--png", type=Path, default=None, help="directory for elevation renders")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    from app import persistence
    from app import world as world_mod

    results = {"earth": EARTH, "runs": []}
    stops = sorted(float(t) for t in args.step_myr.split(",") if t)
    sources = [("seed", int(s)) for s in args.seeds.split(",") if s] + [("save", p) for p in args.load]
    for kind, source in sources:
        start = time.time()
        if kind == "seed":
            world = world_mod.generate_world(
                source, cratons_first=args.cratons_first, synthetic_relief=args.synthetic_relief,
                node_density=args.node_density, land_fraction=args.land_fraction,
                continental_fraction=args.continental_fraction,
            )
            label = f"seed{source}"
        else:
            world = persistence.load_world_bytes(Path(source).read_bytes())
            label = Path(source).stem
        run = {"label": label, "seconds_to_world": round(time.time() - start, 2), "measures": [measure(world)]}
        if args.png:
            _render(world, args.png / f"{label}-0myr.png")
        for stop in stops:
            while world.elapsed_years < stop * 1e6 - 1.0:
                world_mod.step_world(world, min(args.chunk_myr * 1e6, stop * 1e6 - world.elapsed_years))
            run["measures"].append(measure(world))
            if args.png:
                _render(world, args.png / f"{label}-{stop:g}myr.png")
        results["runs"].append(run)
        print(json.dumps({"label": label, "land_p50": run["measures"][0]["land_pct_m"].get("50")}), flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    main()
