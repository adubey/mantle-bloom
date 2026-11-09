#!/usr/bin/env python3
"""Issue #835: does coastal erosion breach a low sill between a sub-sea basin and open ocean?

Builds a controlled case in a generated world and steps the real `erosion.apply_erosion` with
the water-level datum held fixed (`apply_erosion` never moves sea level, and no other process
runs). On a continental patch the script carves: open ocean to the west, a flat sub-sea basin
to the east, land elsewhere, and a one-node-wide channel (the sill) from the ocean to the
basin at a chosen height above sea level. After each step it reports

* the sill's spill barrier -- the lowest over all paths from the basin to the open ocean of
  the highest node on the path (minimax over the hydrology neighbour graph),
* whether the basin is in the world sea, by the same rule the world uses (largest
  at-or-below-sea-level component, `hydrology._below_sea_components`),
* each leveling term at the sill nodes (wave/frost erosion before caps, the leveling grind,
  the datum, the fill received), and where the volume went (sill, basin, ocean, land).

Usage (from the repo root; one heavy run at a time, see the compute-budget skill):

    backend/.venv/bin/python bin/debug/diagnose_low_sill_breach.py \\
        --sill-heights=-3,0.5,1,2,3,6 --step-years 10000 --steps 30 \\
        --out ~/mantle-bloom-captures/issue835/sill-sweep.json
"""

from __future__ import annotations

import argparse
import copy
import heapq
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

from app import plates  # noqa: E402
from app.elevation_lines import PLANET_RADIUS_KM  # noqa: E402
from app.hydroclimate import erosion, hydrology  # noqa: E402
from app.world import generate_world  # noqa: E402

REGION_KM = 900.0
OCEAN_EDGE_KM = -120.0
BASIN_NEAR_KM, BASIN_FAR_KM, BASIN_HALF_WIDTH_KM = 120.0, 500.0, 300.0
OCEAN_FLOOR_M, BASIN_FLOOR_M, LAND_M = -1500.0, -35.0, 80.0


def _nodes(world):
    points, plates_in_order = plates.gather_node_positions(world.plates)
    return points, plates_in_order, plates.collect_all_elevation(plates_in_order)


def _write_elevation(plates_in_order, elevation):
    start = 0
    for plate in plates_in_order:
        count = plate.node_count()
        plate.set_fields_on_plate(
            elevation=np.ascontiguousarray(elevation[start : start + count])
        )
        start += count


def _pick_center(points, elevation):
    """A land node 350-600 km from the nearest real ocean node, with the most land nearby, so the
    patch sits on continental crust and its carved ocean side joins the world sea. Returns
    (center index, nearest ocean node index)."""
    from scipy.spatial import cKDTree

    chord = 2.0 * np.sin(REGION_KM / PLANET_RADIUS_KM / 2.0)
    ocean_idx = np.nonzero(elevation < -200.0)[0]
    land = np.nonzero(elevation > 0.0)[0]
    dist, nearest = cKDTree(points[ocean_idx]).query(points[land])
    km = dist * PLANET_RADIUS_KM
    ok = (km >= 350.0) & (km <= 600.0)
    if not np.any(ok):
        raise SystemExit(
            "no land node 350-600 km from the open ocean in this world; try another --seed"
        )
    candidates = land[ok]
    counts = cKDTree(points).query_ball_point(
        points[candidates], chord, return_length=True
    )
    best = int(np.argmax(counts))
    return int(candidates[best]), int(ocean_idx[nearest[ok][best]])


def _carve(
    points,
    center,
    ocean_node,
    sill_height_m,
    spacing_km,
    land_m=LAND_M,
    ocean_floor_m=OCEAN_FLOOR_M,
):
    """Return (elevation, masks) for the controlled layout around `center`, with the open ocean
    toward `ocean_node` (-x) and the basin on the far side (+x)."""
    c = points[center]
    toward_ocean = points[ocean_node] - c
    toward_ocean -= (toward_ocean @ c) * c
    east = -toward_ocean / np.linalg.norm(toward_ocean)
    north = np.cross(c, east)
    d = (points - c) * PLANET_RADIUS_KM
    x, y = d @ east, d @ north
    in_region = (points @ c > 0.0) & (np.hypot(x, y) < REGION_KM)
    # one-node-wide channel: the nodes nearest the x axis between ocean and basin
    along = in_region & (x >= OCEAN_EDGE_KM) & (x <= BASIN_NEAR_KM)
    channel = along & (np.abs(y) <= spacing_km * 0.5)
    ocean = in_region & (x < OCEAN_EDGE_KM)
    basin = (
        in_region
        & (x > BASIN_NEAR_KM)
        & (x < BASIN_FAR_KM)
        & (np.abs(y) < BASIN_HALF_WIDTH_KM)
    )
    land = in_region & ~channel & ~ocean & ~basin
    elevation = np.full(len(points), np.nan)
    elevation[land], elevation[ocean], elevation[basin], elevation[channel] = (
        land_m,
        ocean_floor_m,
        BASIN_FLOOR_M,
        sill_height_m,
    )
    return elevation, {
        "channel": channel,
        "ocean": ocean,
        "basin": basin,
        "land": land,
        "region": in_region,
    }


def _minimax_barrier(elevation, neighbor_idx, sources, targets):
    """Lowest possible highest-node over any path from `sources` to `targets`."""
    best = np.full(len(elevation), np.inf)
    heap = []
    for s in np.nonzero(sources)[0]:
        best[s] = elevation[s]
        heapq.heappush(heap, (best[s], int(s)))
    while heap:
        h, i = heapq.heappop(heap)
        if h > best[i]:
            continue
        if targets[i]:
            return float(h)
        for j in neighbor_idx[i]:
            nh = max(h, elevation[j])
            if nh < best[j]:
                best[j] = nh
                heapq.heappush(heap, (nh, int(j)))
    return float("inf")


def _connected_to_world_sea(world, points, elevation, area, basin):
    neighbor_idx = hydrology._build_neighbor_graph(points, world=world)
    below = elevation <= world.sea_level_m
    labels = hydrology._below_sea_components(below, neighbor_idx)
    sizes = np.bincount(labels[below], weights=area[below])
    sea = labels == int(np.argmax(sizes))
    return bool(np.any(sea & basin & below)), neighbor_idx


def _apply_ablation(mode):
    """Switch off one leveling term for an attribution run; returns the originals to restore.

    `no-fill`: `_spread_coastal_leveling` places nothing, so every source keeps its own rock
    (the function's own no-sink branch; mass stays closed). `datum-sea`: every node's leveling
    datum is the waterline, so the fill cannot lift a node above sea level and the grind
    planes anything above it."""
    originals = {}
    if mode == "no-fill":
        originals["_spread_coastal_leveling"] = erosion._spread_coastal_leveling

        def keep_all(
            points,
            elevation,
            openness,
            dist_to_land,
            sea_level_m,
            datum_m,
            source_amount,
            *a,
            **k,
        ):
            tagged = k.get("tagged_amount")
            return source_amount.copy(), (
                np.zeros_like(source_amount) if tagged is None else tagged.copy()
            )

        erosion._spread_coastal_leveling = keep_all
    elif mode == "datum-sea":
        originals["leveling_datum_m"] = erosion.leveling_datum_m

        def sea_datum(sea_level_m, openness, dist_to_land):
            return np.full(len(openness), float(sea_level_m))

        erosion.leveling_datum_m = sea_datum
    return originals


def _water_directions_two_sided(points, water, candidates, fetch_km):
    """For each candidate node: (any water within fetch_km, water on roughly opposite sides).
    Opposite = some pair of water neighbours whose tangent directions have cosine <= -0.3."""
    from scipy.spatial import cKDTree

    n = len(points)
    near, opposite = np.zeros(n, dtype=bool), np.zeros(n, dtype=bool)
    water_idx = np.nonzero(water)[0]
    if len(water_idx) == 0:
        return near, opposite
    tree = cKDTree(points[water_idx])
    chord = 2.0 * np.sin(fetch_km / PLANET_RADIUS_KM / 2.0)
    for i in np.nonzero(candidates)[0]:
        found = tree.query_ball_point(points[i], chord)
        if not found:
            continue
        near[i] = True
        v = points[water_idx[found]] - points[i]
        v -= (v @ points[i])[:, None] * points[i]
        norm = np.linalg.norm(v, axis=1)
        v = v[norm > 1e-12] / norm[norm > 1e-12, None]
        if len(v) >= 2 and (v @ v.T).min() <= -0.3:
            opposite[i] = True
    return near, opposite


def _apply_variants(args, points, get_elevation):
    """Experimental model variants for #835, applied as monkeypatches (not model changes).

    cap-variant: lets wave erosion cut a bench down to `reach_m` below the waterline at nodes
    near water (`any-water`: water within WAVE_FETCH_KM; `two-sided`: water on roughly opposite
    sides), by raising their drop-to-lowest-neighbour cap to that allowance.
    sea-relative: coastal and submarine erosion read height above sea level, not raw elevation.
    open-all-water: wave-exposure openness counts every below-sea node, enclosed seas included."""
    sea = args.sea_level
    restore = {}
    if args.infill_fraction is not None:
        restore["COASTAL_INFILL_MARINE_FRACTION"] = (
            erosion.COASTAL_INFILL_MARINE_FRACTION
        )
        erosion.COASTAL_INFILL_MARINE_FRACTION = args.infill_fraction
    if args.cap_variant != "none":
        restore["compute_slope"] = erosion.compute_slope
        original = erosion.compute_slope

        def capped(pts, elevation, world=None):
            slope, drop = original(pts, elevation, world=world)
            candidates = (elevation - sea <= 6.0) & (elevation - sea >= -args.reach_m)
            near, opposite = _water_directions_two_sided(
                pts, elevation < sea, candidates, args.fetch_km
            )
            flagged = opposite if args.cap_variant == "two-sided" else near
            allowance = np.clip(elevation - (sea - args.reach_m), 0.0, None)
            return slope, np.where(flagged, np.maximum(drop, allowance), drop)

        erosion.compute_slope = capped
    if args.drop_source_return:
        # Rock the spreads put back on a node that is itself a source this step (including
        # the keep-it-on-the-source fallback) is discarded instead. Breaks the budget.
        restore["_spread_marine_sediment"] = erosion._spread_marine_sediment
        restore["_spread_coastal_leveling"] = erosion._spread_coastal_leveling
        marine_fn, leveling_fn = (
            erosion._spread_marine_sediment,
            erosion._spread_coastal_leveling,
        )

        def marine_drop(pts, elevation, is_ocean, source):
            return np.where(
                source > 0, 0.0, marine_fn(pts, elevation, is_ocean, source)
            )

        def leveling_drop(*a, **k):
            fill, tagged = leveling_fn(*a, **k)
            source = a[6] > 0
            return np.where(source, 0.0, fill), np.where(source, 0.0, tagged)

        erosion._spread_marine_sediment = marine_drop
        erosion._spread_coastal_leveling = leveling_drop
    if args.no_refill_eroding:
        # Receivers of the marine and leveling spreads exclude nodes that lose sea-side rock
        # this step; unplaced rock falls back to its source, so mass is still conserved.
        restore["_spread_marine_sediment"] = erosion._spread_marine_sediment
        restore["_spread_coastal_leveling"] = erosion._spread_coastal_leveling
        marine, leveling = (
            erosion._spread_marine_sediment,
            erosion._spread_coastal_leveling,
        )

        def marine_wrapper(pts, elevation, is_ocean, source):
            return marine(pts, elevation, is_ocean & ~(source > 0), source)

        def leveling_wrapper(*a, **k):
            frame = sys._getframe(1)
            while "sea_side_erosion" not in frame.f_locals:  # tap wrappers sit between
                frame = frame.f_back
            eroding = frame.f_locals["sea_side_erosion"] > 0
            a = list(a)
            a[5] = np.where(eroding, -1e9, a[5])
            return leveling(*a, **k)

        erosion._spread_marine_sediment = marine_wrapper
        erosion._spread_coastal_leveling = leveling_wrapper
    if args.sea_relative:
        restore["coastal_erosion_amount"] = erosion.coastal_erosion_amount
        restore["submarine_erosion_amount"] = erosion.submarine_erosion_amount
        coastal, submarine = (
            erosion.coastal_erosion_amount,
            erosion.submarine_erosion_amount,
        )
        erosion.coastal_erosion_amount = lambda elevation, temperature_c, dt_myr: (
            coastal(elevation - sea, temperature_c, dt_myr)
        )
        erosion.submarine_erosion_amount = lambda elevation, slope, is_ocean, dt_myr: (
            submarine(elevation - sea, slope, is_ocean, dt_myr)
        )
    if args.open_all_water:
        restore["_coastal_openness"] = erosion._coastal_openness
        openness = erosion._coastal_openness
        erosion._coastal_openness = lambda pts, is_ocean: openness(
            pts, np.asarray(is_ocean, dtype=bool) | (get_elevation() < sea)
        )
    return restore


def census(world, points, elevation):
    """How much of the model's wave/frost band is far from any water? `coastal_erosion_amount`
    reads only elevation, so a low inland plain gets it too; the drop cap is what limits it."""
    from scipy.spatial import cKDTree

    sea = world.sea_level_m
    band = np.abs(elevation - sea) < erosion.COASTAL_EROSION_BAND_M
    weight = np.clip(
        1.0 - np.abs(elevation - sea) / erosion.COASTAL_EROSION_BAND_M, 0.0, 1.0
    )
    water = elevation < sea
    dist_km = cKDTree(points[water]).query(points)[0] * PLANET_RADIUS_KM
    out = {
        "sea_level_m": float(sea),
        "nodes": int(len(points)),
        "band_nodes": int(band.sum()),
    }
    for limit in (0.0, 150.0, 500.0):
        far = band & (dist_km > limit)
        out[f"band_nodes_beyond_{int(limit)}km_of_water"] = int(far.sum())
        out[f"band_weight_share_beyond_{int(limit)}km_of_water"] = float(
            weight[far].sum() / weight[band].sum()
        )
    out["land_nodes_in_band_inland_over_500km"] = int(
        (band & (elevation > sea) & (dist_km > 500.0)).sum()
    )
    return out


def run_case(base_world, center, ocean_node, sill_height_m, args, spacing_km):
    world = copy.deepcopy(base_world)
    world.sea_level_m = args.sea_level
    points, plates_in_order, _ = _nodes(world)
    elevation, masks = _carve(
        points,
        center,
        ocean_node,
        sill_height_m,
        spacing_km,
        args.land_m,
        args.ocean_floor_m,
    )
    elevation = (
        elevation + args.sea_level
    )  # layout heights are relative to the waterline
    keep = np.isfinite(elevation)
    base = plates.collect_all_elevation(plates_in_order)
    elevation = np.where(keep, elevation, base)
    _write_elevation(plates_in_order, elevation)
    area = erosion._gather_areas(world, plates_in_order)

    taps: dict[str, np.ndarray] = {}
    originals = {}
    ablated = _apply_ablation(args.ablate)
    variants = _apply_variants(
        args, points, lambda: plates.collect_all_elevation(plates_in_order)
    )

    ch = masks["channel"]

    def tap(name, pick):
        original = getattr(erosion, name)
        originals[name] = original

        def wrapper(*a, **k):
            out = original(*a, **k)
            taps[name] = pick(out, k)
            if name == "_spread_coastal_leveling":
                # what the pass was handed at the sill (m), against what it returned there
                taps["leveling_source_in"] = a[6][ch] / k["area_m2"][ch]
            return out

        setattr(erosion, name, wrapper)

    def fill_thickness(out, k):
        # Returns a volume because apply_erosion passes area_m2; divide back to metres.
        return out[0][ch] / k["area_m2"][ch]

    tap("coastal_erosion_amount", lambda out, k: out[ch].copy())
    tap("coastal_leveling_grind", lambda out, k: out[ch].copy())
    tap("leveling_datum_m", lambda out, k: out[ch].copy())
    tap("_spread_coastal_leveling", fill_thickness)
    tap("_coastal_openness", lambda out, k: out[ch].copy())
    tap("compute_slope", lambda out, k: out[1][ch].copy())

    # Per-term attribution at the sill (#835). apply_erosion builds `incoming` (deposit volume
    # by pathway) and its removal terms as locals just before _carry_overflow; read them from
    # the caller's frame instead of editing the model.
    original_carry = erosion._carry_overflow
    originals["_carry_overflow"] = original_carry

    def carry_wrapper(*a, **k):
        out = original_carry(*a, **k)
        frame = sys._getframe(1).f_locals
        area_ch = area[ch]
        incoming = frame["incoming"][ch][:, :-1] / area_ch[:, None]
        for col, name in enumerate(erosion.OVERFLOW_PATHWAYS):
            taps[f"deposit_{name}_m"] = incoming[:, col].copy()
        taps["deposit_overflow_net_m"] = (
            out.placed[ch][:, :-1] - frame["excess"][ch][:, :-1]
        )[:, erosion.OVERFLOW_SEDIMENT].sum(axis=1) / area_ch
        for name in (
            "erosion_amount",
            "sea_side_erosion",
            "ground_off",
            "flatten_removed",
            "removable_m",
            "drop_to_lowest_neighbor_m",
        ):
            taps[f"removal_{name}"] = np.asarray(frame[name])[ch].copy()
        for pathway in args.off_pathway:
            # Diagnostic only (breaks the erosion budget): drop one deposit pathway's column.
            col = erosion.OVERFLOW_PATHWAYS.index(pathway)
            frame["incoming"][:, col] = 0.0
            frame["excess"][:, col] = 0.0
            out.placed[:, col] = 0.0
        if args.no_deposit:
            # Diagnostic only (breaks the erosion budget): zero every pathway deposit in
            # place after the taps read it, so apply_erosion's elevation update sees none.
            frame["incoming"][:] = 0.0
            frame["excess"][:] = 0.0
            out.placed[:] = 0.0
        return out

    erosion._carry_overflow = carry_wrapper

    history = []
    breach_years = None
    try:
        for step in range(args.steps):
            before = plates.collect_all_elevation(plates_in_order).copy()
            hc_before = plates.collect_all_crustal_thickness(plates_in_order).copy()
            result = erosion.apply_erosion(world, years=args.step_years)
            after = plates.collect_all_elevation(plates_in_order)
            delta = after - before
            hc_change = (
                plates.collect_all_crustal_thickness(plates_in_order) - hc_before
            )
            volume = {
                name: float((delta * area)[m].sum())
                for name, m in masks.items()
                if name != "region"
            }
            volume["other"] = float((delta * area)[~masks["region"]].sum())
            connected, neighbor_idx = _connected_to_world_sea(
                world, points, after, area, masks["basin"]
            )
            barrier = _minimax_barrier(
                after, neighbor_idx, masks["basin"], masks["ocean"]
            )
            years = (step + 1) * args.step_years
            if breach_years is None and connected:
                breach_years = years
            history.append(
                {
                    "years": years,
                    "sill_max_m": float(after[ch].max()),
                    "sill_min_m": float(after[ch].min()),
                    "sill_elevation_change_m": [round(float(v), 3) for v in delta[ch]],
                    "sill_crust_thickness_change_m": [
                        round(float(v), 3) for v in hc_change[ch]
                    ],
                    "spill_barrier_m": barrier,
                    "basin_in_world_sea": connected,
                    "terms_at_sill": {
                        k: [round(float(v), 4) for v in taps[k]] for k in taps
                    },
                    "volume_change_m3": volume,
                    "budget_removed_m3": float(result.budget["removed_m3"])
                    if result
                    else None,
                    "budget_deposited_m3": float(result.budget["deposited_m3"])
                    if result
                    else None,
                }
            )
            if connected and args.stop_on_breach:
                break
    finally:
        for name, fn in originals.items():
            setattr(erosion, name, fn)
        for name, fn in {**ablated, **variants}.items():
            setattr(erosion, name, fn)
    return {
        "sill_height_m": sill_height_m,
        "channel_nodes": int(ch.sum()),
        "mask_nodes": {name: int(m.sum()) for name, m in masks.items()},
        "breach_years": breach_years,
        "history": history,
    }


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--seed", type=int, default=23)
    ap.add_argument("--num-plates", type=int, default=8)
    ap.add_argument("--node-density", type=float, default=1.0)
    ap.add_argument(
        "--sill-heights",
        default="-3,0.5,1,2,3,6",
        help="comma list, metres relative to sea level",
    )
    ap.add_argument("--step-years", type=float, default=10_000.0)
    ap.add_argument("--steps", type=int, default=30)
    ap.add_argument("--stop-on-breach", action="store_true")
    ap.add_argument(
        "--land-m",
        type=float,
        default=LAND_M,
        help="surrounding land elevation; set near the sill height as a control",
    )
    ap.add_argument(
        "--ocean-floor-m",
        type=float,
        default=OCEAN_FLOOR_M,
        help="open-ocean depth beside the sill, relative to sea level",
    )
    ap.add_argument(
        "--infill-fraction",
        type=float,
        default=None,
        help="override COASTAL_INFILL_MARINE_FRACTION (share of wave/submarine erosion sent to the leveling pool)",
    )
    ap.add_argument(
        "--ablate", choices=("none", "no-fill", "datum-sea"), default="none"
    )
    ap.add_argument(
        "--sea-level",
        type=float,
        default=0.0,
        help="fixed sea level (m); the layout is built relative to it",
    )
    ap.add_argument(
        "--cap-variant", choices=("none", "any-water", "two-sided"), default="none"
    )
    ap.add_argument(
        "--reach-m",
        type=float,
        default=2.6,
        help="wave-cut depth below the waterline for --cap-variant (2 m waves break near 2.6 m)",
    )
    ap.add_argument("--fetch-km", type=float, default=150.0)
    ap.add_argument(
        "--sea-relative",
        action="store_true",
        help="coastal/submarine erosion read height above sea level",
    )
    ap.add_argument(
        "--off-pathway",
        action="append",
        default=[],
        choices=list(erosion.OVERFLOW_PATHWAYS),
        help="drop one deposit pathway (repeatable; diagnostic, breaks the budget)",
    )
    ap.add_argument(
        "--drop-source-return",
        action="store_true",
        help="discard rock the spreads put back on nodes that are sources this step",
    )
    ap.add_argument(
        "--no-deposit",
        action="store_true",
        help="zero every pathway deposit in apply_erosion (diagnostic; breaks the budget)",
    )
    ap.add_argument(
        "--no-refill-eroding",
        action="store_true",
        help="marine and leveling spreads may not deposit on nodes losing sea-side rock",
    )
    ap.add_argument(
        "--open-all-water",
        action="store_true",
        help="openness counts enclosed seas as open water",
    )
    ap.add_argument(
        "--census",
        action="store_true",
        help="only report the wave band's distance from water, then exit",
    )
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    world = generate_world(
        seed=args.seed, num_plates=args.num_plates, node_density=args.node_density
    )
    points, _, elevation = _nodes(world)
    from scipy.spatial import cKDTree

    spacing_km = float(
        np.median(cKDTree(points).query(points, k=2)[0][:, 1]) * PLANET_RADIUS_KM
    )
    if args.census:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(census(world, points, elevation), indent=1))
        print(args.out.read_text())
        return
    center, ocean_node = _pick_center(points, elevation)
    cases = [
        run_case(world, center, ocean_node, float(h), args, spacing_km)
        for h in args.sill_heights.split(",")
    ]
    report = {
        "seed": args.seed,
        "node_density": args.node_density,
        "node_spacing_km": spacing_km,
        "step_years": args.step_years,
        "ablate": args.ablate,
        "land_m": args.land_m,
        "ocean_floor_m": args.ocean_floor_m,
        "sea_level_m": args.sea_level,
        "cap_variant": args.cap_variant,
        "reach_m": args.reach_m,
        "sea_relative": args.sea_relative,
        "open_all_water": args.open_all_water,
        "cases": cases,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=1))
    for case in cases:
        print(
            f"sill {case['sill_height_m']:+.1f} m: breach at {case['breach_years']} yr; final barrier {case['history'][-1]['spill_barrier_m']:.2f} m"
        )
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
