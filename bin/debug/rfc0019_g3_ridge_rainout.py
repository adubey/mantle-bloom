"""Test a conservative windward rainout rule on captured G3 tuning fields.

Usage: backend/.venv/bin/python bin/debug/rfc0019_g3_ridge_rainout.py CAPTURE [--maps DIR]
Compares the half-strength broad transport and overflow spread with two
interception strengths. Use --weak-sensitivity for the single weaker arm,
--mid-sensitivity for the intermediate 5% arm, or --dry-sensitivity to limit
the 10% arm by the control donor moisture. Use --strength-sensitivity to
compare all four strengths against the same control.
No saved world is loaded or stepped.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np

import rfc0019_g3_compare as g3
from rfc0019_coarse_moisture import (
    FaceInterception, TransportFaces, aggregate_inputs, build_transport_rates, distribute_rain,
    solve_bounded_budget,
)
from rfc0019_g3_archive_cap import (
    _belt_contrasts,
    _belt_row_shares,
    _broad_meridional_capacity,
    _convergence_correction,
    _distribute_overflow,
    _signed_faces,
)
from rfc0019_g3_bounded_column import COLD_COAST_DISTANCE_KM
from rfc0019_g3_mixing import GRID_SHAPE, RAIN_RATE, TUNING_SAVES
from rfc0019_g3_mixing import _regional_metrics

RIDGE_RISE_M = 250.0
RISE_SCALE_M = 1000.0
DRY_SCALE_KG_M2 = 10.0
ARMS = {"control": 0.0, "ridge_010": 0.10, "ridge_025": 0.25}
WEAK_ARMS = {"control": 0.0, "ridge_0025": 0.025}
MID_ARMS = {"control": 0.0, "ridge_005": 0.05}
DRY_ARMS = {"control": 0.0, "ridge_010_dry": 0.10}
STRENGTH_ARMS = {
    "control": 0.0, "ridge_0025": 0.025, "ridge_005": 0.05,
    "ridge_010": 0.10, "ridge_025": 0.25,
}


@dataclass(frozen=True)
class Crossings:
    face: np.ndarray
    toward_second: np.ndarray
    donor: np.ndarray
    crest: np.ndarray
    lee: np.ndarray
    rise_m: np.ndarray


def _ridge_crossings(grid, faces, elevation):
    """Follow each directed face across a crest and onward to the next lower land cell."""
    first, second = faces.first, faces.second
    nlat, nlon = grid.coarse_shape
    row, col = np.divmod(grid.coarse_cell, nlon)
    east_axis = row[first] == row[second]
    signed = faces.forward - faces.reverse
    toward_second = signed > 0.0
    step = np.sign(signed).astype(int)
    entry = np.where(toward_second, first, second)
    crest = np.where(toward_second, second, first)
    lee_row = row[crest] + step * ~east_axis
    lee_col = (col[crest] + step * east_axis) % nlon
    land_reservoir = np.full(nlat * nlon, -1, dtype=np.intp)
    land_reservoir[grid.coarse_cell[~grid.is_ocean]] = np.flatnonzero(~grid.is_ocean)
    lee_cell = np.clip(lee_row, 0, nlat - 1) * nlon + lee_col
    lee = land_reservoir[lee_cell]

    # Check the exact crest-to-lee face, summing its native boundary segments.
    face_pairs = np.column_stack((first, second))
    unique_pairs, pair_index = np.unique(face_pairs, axis=0, return_inverse=True)
    pair_flow = np.bincount(pair_index, weights=signed)
    flow_by_pair = {tuple(pair): value for pair, value in zip(unique_pairs, pair_flow, strict=True)}
    exit_first = np.where(step > 0, crest, lee)
    exit_second = np.where(step > 0, lee, crest)
    exit_flow = np.fromiter(
        (flow_by_pair.get((int(a), int(b)), 0.0) for a, b in zip(exit_first, exit_second, strict=True)),
        dtype=float, count=step.size,
    ) * step
    rise = elevation[crest] - elevation[entry]
    fall = elevation[crest] - elevation[np.maximum(lee, 0)]
    land = ~grid.is_ocean
    valid = ((step != 0) & (lee_row >= 0) & (lee_row < nlat) & (lee >= 0) & (exit_flow > 0)
             & land[entry] & land[crest] & land[np.maximum(lee, 0)]
             & (rise >= RIDGE_RISE_M) & (fall >= RIDGE_RISE_M))
    return Crossings(np.flatnonzero(valid), toward_second[valid], entry[valid], crest[valid], lee[valid], rise[valid])


def _interception(crossings, face_count, strength):
    forward = np.zeros(face_count)
    reverse = np.zeros(face_count)
    fraction = -strength * np.expm1(-crossings.rise_m / RISE_SCALE_M)
    east_or_north = crossings.toward_second
    forward[crossings.face[east_or_north]] = fraction[east_or_north]
    reverse[crossings.face[~east_or_north]] = fraction[~east_or_north]
    return FaceInterception(forward, reverse)


def _paired_rain(crossings, rain, control, grid):
    if not crossings.face.size:
        return {"crossings": 0}
    crest, lee = crossings.crest, crossings.lee
    native_area = grid.native_area.ravel()
    on_path = np.bincount(
        grid.native_ids.ravel(), weights=rain[::-1].ravel() * native_area, minlength=grid.area.size,
    ) / grid.area
    baseline = np.bincount(
        grid.native_ids.ravel(), weights=control[::-1].ravel() * native_area, minlength=grid.area.size,
    ) / grid.area
    crest_change = on_path[crest] - baseline[crest]
    lee_change = on_path[lee] - baseline[lee]
    return {
        "crossings": int(crest.size),
        "distinct_crests": int(np.unique(crest).size),
        "median_crest_rain_mm_yr": float(np.median(on_path[crest])),
        "median_lee_rain_mm_yr": float(np.median(on_path[lee])),
        "median_crest_minus_lee_mm_yr": float(np.median(on_path[crest] - on_path[lee])),
        "median_contrast_gain_mm_yr": float(np.median(crest_change - lee_change)),
        "median_lee_change_mm_yr": float(np.median(lee_change)),
        "crest_wetter_and_lee_drier_share": float(np.mean((crest_change > 0) & (lee_change < 0))),
    }


def _rain_metrics(fields, rain, area, distance, freeze):
    cold_land = ~fields.is_ocean & (fields.air_temperature_c < freeze)
    coastal = distance <= COLD_COAST_DISTANCE_KM
    metrics = _regional_metrics(fields, rain, area, distance, freeze)
    for name, mask in (("cold_coast", cold_land & coastal), ("cold_inland", cold_land & ~coastal)):
        metrics[name + "_rain_mm_yr"] = g3._weighted_mean(rain[mask], area[mask])
    return metrics | {
        "ocean_belt_contrasts": _belt_contrasts(fields, rain, area),
        "ocean_belt_wettest_row_shares": _belt_row_shares(fields, rain, area),
    }


def _distribution_metrics(fields, rain, area):
    """Summarize global and folded-latitude rain for the G3 density check."""
    mass_rate = rain * area / g3.SECONDS_PER_YEAR
    land = ~fields.is_ocean
    total = float(mass_rate.sum())
    land_total = float(mass_rate[land].sum())
    ocean_total = float(mass_rate[~land].sum())
    bands = []
    absolute_latitude = np.abs(fields.lat_deg)
    for lower in range(0, 90, 15):
        mask = (absolute_latitude >= lower) & (absolute_latitude < lower + 15)
        band_area = area[mask]
        band_land = land[mask]
        bands.append({
            "absolute_latitude_deg": [lower, lower + 15],
            "land_mean_mm_yr": g3._weighted_mean(rain[mask][band_land], band_area[band_land]),
            "ocean_mean_mm_yr": g3._weighted_mean(rain[mask][~band_land], band_area[~band_land]),
        })
    return {
        "total_rain_kg_s": total,
        "global_mean_rain_mm_yr": float(np.sum(rain * area) / area.sum()),
        "land_rain_share": land_total / max(total, 1e-300),
        "ocean_rain_share": ocean_total / max(total, 1e-300),
        "folded_latitude_bands": bands,
    }


def _place_surface_rain(grid, budget, overflow, area, fields, freeze):
    """Place reservoir rain on native cells, spreading only saturation overflow."""
    weights = np.ones(grid.native_ids.shape)
    ordinary_rain = budget.rain - overflow
    tolerance = 100 * np.finfo(float).eps * np.maximum(budget.rain, 1.0)
    if np.any(ordinary_rain < -tolerance):
        raise AssertionError("overflow exceeds the reservoir's total rain")
    ordinary = distribute_rain(grid, np.maximum(ordinary_rain, 0.0), weights)[::-1]
    excess = distribute_rain(grid, overflow, weights)[::-1]
    surface_class = np.where(fields.is_ocean, 2, fields.air_temperature_c < freeze)
    deposition = _distribute_overflow(excess, area, surface_class)
    before = np.bincount(surface_class.ravel(), weights=excess.ravel(), minlength=3)
    after = np.bincount(surface_class.ravel(), weights=deposition.ravel(), minlength=3)
    scale = max(float(budget.rain.sum()), 1e-300)
    class_residual = float(np.abs(after - before).sum() / scale)
    placed = float((ordinary + deposition).sum())
    mass_residual = abs(placed - scale) / scale
    residual = max(class_residual, mass_residual)
    if residual > 1e-12:
        raise AssertionError(f"overflow placement did not preserve surface-class budgets: {residual:g}")
    return (ordinary + deposition) / area * g3.SECONDS_PER_YEAR, residual


def _cold_inland_delivery(fields, distance, freeze, grid, faces, source, budget, surface_rain):
    """Separate atmospheric delivery from any spread across the coastal-distance mask."""
    cold_inland = ~fields.is_ocean & (fields.air_temperature_c < freeze) & (distance > COLD_COAST_DISTANCE_KM)
    inland = np.bincount(
        grid.native_ids.ravel(),
        weights=cold_inland[::-1].ravel() * grid.native_area.ravel(),
        minlength=grid.area.size,
    ) / grid.area
    inward = (inland[faces.second] - inland[faces.first]) * budget.face_flux
    atmospheric_rain = float(np.dot(budget.rain, inland))
    native_surface = surface_rain[::-1].ravel() * grid.native_area.ravel() / g3.SECONDS_PER_YEAR
    local_source = float(np.dot(source, inland))
    if abs(local_source + inward.sum() - atmospheric_rain) > 1e-12 * max(1.0, source.sum()):
        raise AssertionError("cold-inland regional budget did not close")
    inland_native = cold_inland[::-1].ravel()
    return {
        "inflow_kg_s": float(np.maximum(inward, 0).sum()),
        "outflow_kg_s": float(np.maximum(-inward, 0).sum()),
        "net_inflow_kg_s": float(inward.sum()),
        "local_source_kg_s": local_source,
        "atmospheric_rain_kg_s": atmospheric_rain,
        "spread_contribution_kg_s": float(native_surface[inland_native].sum()) - atmospheric_rain,
    }


def _wind_arrow_faces(grid, native_u, native_v):
    """Make coarse, locally directed neighbor transfers with shared coast exchange."""
    nlat, nlon = grid.coarse_shape
    res_area = grid.area
    ids = grid.native_ids.ravel()
    native_area = grid.native_area.ravel()
    u = np.bincount(ids, weights=native_u.ravel() * native_area, minlength=res_area.size) / res_area
    v = np.bincount(ids, weights=native_v.ravel() * native_area, minlength=res_area.size) / res_area
    row, col = np.divmod(grid.coarse_cell, nlon)
    class_lookup = np.full((nlat * nlon, 2), -1, dtype=np.intp)
    class_lookup[grid.coarse_cell, grid.is_ocean.astype(np.intp)] = np.arange(res_area.size)
    edge = np.linspace(-np.pi / 2, np.pi / 2, nlat + 1)
    latitude = 0.5 * (edge[:-1] + edge[1:])
    longitude_step = 2 * np.pi / nlon
    arcs = []

    def add_arc(origin, row_step, col_step, speed, distance):
        if speed == 0.0:
            return
        dest_row = row[origin] + row_step
        if not 0 <= dest_row < nlat:
            return
        dest_col = (col[origin] + col_step) % nlon
        cell = dest_row * nlon + dest_col
        surface = int(grid.is_ocean[origin])
        dest = class_lookup[cell, surface]
        if dest < 0:
            dest = class_lookup[cell, 1 - surface]
        if dest >= 0:
            arcs.append((int(origin), int(dest), abs(float(speed)) / distance))

    for origin in range(res_area.size):
        if u[origin] != 0.0:
            step = 1 if u[origin] > 0.0 else -1
            east_distance = 6_371_000.0 * np.cos(latitude[row[origin]]) * longitude_step
            add_arc(origin, 0, step, u[origin], east_distance)
        if v[origin] != 0.0:
            step = 1 if v[origin] > 0.0 else -1
            neighbor = row[origin] + step
            if 0 <= neighbor < nlat:
                north_distance = 6_371_000.0 * abs(latitude[neighbor] - latitude[row[origin]])
                add_arc(origin, step, 0, v[origin], north_distance)

    flux_by_pair = {}
    for origin, destination, rate in arcs:
        if row[origin] == row[destination]:
            eastward = (col[destination] - col[origin]) % nlon == 1
            first, second = (origin, destination) if eastward else (destination, origin)
        else:
            northward = row[destination] == row[origin] + 1
            first, second = (origin, destination) if northward else (destination, origin)
        pair = flux_by_pair.setdefault((first, second), [0.0, 0.0])
        direction = 0 if (first, second) == (origin, destination) else 1
        pair[direction] += rate * res_area[origin]

    if flux_by_pair:
        pairs = sorted(flux_by_pair)
        first = np.fromiter((pair[0] for pair in pairs), dtype=np.intp)
        second = np.fromiter((pair[1] for pair in pairs), dtype=np.intp)
        forward = np.fromiter((flux_by_pair[pair][0] for pair in pairs), dtype=float)
        reverse = np.fromiter((flux_by_pair[pair][1] for pair in pairs), dtype=float)
    else:
        first = second = np.array([], dtype=np.intp)
        forward = reverse = np.array([], dtype=float)
    coast = build_transport_rates(
        grid, np.zeros(grid.native_ids.shape), np.zeros(grid.native_ids.shape), 0.0,
    )
    has_exchange = coast.exchange > 0.0
    return TransportFaces(
        np.concatenate((first, coast.first[has_exchange])),
        np.concatenate((second, coast.second[has_exchange])),
        np.concatenate((forward, coast.forward[has_exchange])),
        np.concatenate((reverse, coast.reverse[has_exchange])),
        np.concatenate((np.zeros(forward.size), coast.exchange[has_exchange])),
    )


def _solve_arms(fields, archive, freeze, arms, *, solver="finite_volume", dry_scale=None):
    density_factor = fields.is_ocean.shape[0] // GRID_SHAPE[0]
    if fields.is_ocean.shape != (GRID_SHAPE[0] * density_factor, GRID_SHAPE[1] * density_factor):
        raise ValueError("capture grid must be density 1 or 4")
    grid = aggregate_inputs(fields.is_ocean[::-1], {}, factor=density_factor)
    if solver == "finite_volume":
        original = build_transport_rates(grid, fields.wind_u[::-1], fields.wind_v[::-1], 0.0)
        signed = original.forward - original.reverse
        correction = _convergence_correction(grid, original.first, original.second, signed)
        capacity = signed - correction + 0.5 * _broad_meridional_capacity(grid, original)
        faces = _signed_faces(original.first, original.second, capacity)
    elif solver == "wind_arrows":
        faces = _wind_arrow_faces(grid, fields.wind_u[::-1], fields.wind_v[::-1])
    else:
        raise ValueError(f"unsupported transport solver: {solver}")
    native_elevation = fields.elevation_m[::-1].ravel()
    elevation_by_reservoir = np.bincount(
        grid.native_ids.ravel(),
        weights=native_elevation * grid.native_area.ravel(),
        minlength=grid.area.size,
    ) / grid.area
    crossings = _ridge_crossings(grid, faces, elevation_by_reservoir)
    source = archive["source_kg_s"]
    cap = archive["cap_kg_m2"] * grid.area
    sink = np.full_like(source, RAIN_RATE)
    area = grid.native_area[::-1]
    distance = archive["coast_distance_km"]
    reference = archive["reference_precipitation_mm_yr"]
    maps, results = {}, {}
    moisture_factor = np.ones_like(crossings.rise_m)
    for name, strength in arms.items():
        interception = _interception(crossings, faces.first.size, strength * moisture_factor)
        budget, excess, _ = solve_bounded_budget(grid.area, source, sink, faces, cap, interception=interception)
        if dry_scale is not None and name == "control":
            donor = crossings.donor
            donor_column = budget.mass[donor] / grid.area[donor]
            moisture_factor = np.minimum(donor_column / dry_scale, 1.0)**2
        rain, placement_residual = _place_surface_rain(grid, budget, excess, area, fields, freeze)
        maps[name] = rain
        results[name] = _rain_metrics(fields, rain, area, distance, freeze) | _distribution_metrics(fields, rain, area) | {
            "ridge_transects": _paired_rain(crossings, rain, maps["control"], grid),
            "cold_inland_budget": _cold_inland_delivery(fields, distance, freeze, grid, faces, source, budget, rain),
            "raw_rms_to_legacy": g3._rms_pair(reference, rain, area)["raw_rel_rms"],
            "summed_cell_closure": budget.summed_cell_relative,
            "condensate_placement_relative_residual": placement_residual,
            "intercepted_rain_fraction": float(budget.intercepted_rain_kg_s.sum() / budget.rain.sum()),
        }
    return results, maps, area


def _compare_capture(capture, map_dir, latitude, freeze, arms, *, solver="finite_volume", dry_scale=None):
    provenance = json.loads((capture / "results.json").read_text())
    save = Path(provenance["save"])
    if TUNING_SAVES.get(save.name) != provenance["save_sha256"]:
        raise ValueError("expected a pinned tuning capture")
    density = int(provenance.get("density", provenance["grid_shape"][0] // GRID_SHAPE[0]))
    previous_path = capture / "overflow-placement" / "results.json"
    previous = json.loads(previous_path.read_text()) if previous_path.exists() else None
    result = {}
    for wind in ("diagnostic", "cfd"):
        with np.load(capture / f"{wind}-columns.npz", allow_pickle=False) as archive:
            fields = SimpleNamespace(
                is_ocean=archive["is_ocean"], lat_deg=latitude,
                air_temperature_c=archive["air_temperature_c"], elevation_m=archive["elevation_m"],
                wind_u=archive["wind_u_m_s"], wind_v=archive["wind_v_m_s"],
            )
            result[wind], maps, area = _solve_arms(
                fields, archive, freeze, arms, solver=solver, dry_scale=dry_scale,
            )
        if previous is not None:
            expected = previous[wind]["broad_transport_half_spread"]["raw_rms_to_legacy"]
            np.testing.assert_allclose(result[wind]["control"]["raw_rms_to_legacy"], expected, rtol=1e-10, atol=0)
        if map_dir is not None:
            for name in arms.keys() - {"control"}:
                if name != "ridge_025":
                    continue
                g3._save_common_scale_rain_map(
                    map_dir, save, density, wind, solver, name, maps["control"], maps[name], area,
                    labels=("Broad transport and spread", "With ridge rainout", "Change"),
                )
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture", type=Path)
    parser.add_argument("--maps", type=Path)
    parser.add_argument("--solver", choices=("finite_volume", "wind_arrows", "both"), default="both")
    sensitivity = parser.add_mutually_exclusive_group()
    sensitivity.add_argument("--weak-sensitivity", action="store_const", dest="profile", const="weak")
    sensitivity.add_argument("--mid-sensitivity", action="store_const", dest="profile", const="mid")
    sensitivity.add_argument("--dry-sensitivity", action="store_const", dest="profile", const="dry")
    sensitivity.add_argument("--strength-sensitivity", action="store_const", dest="profile", const="strengths")
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))
    from app.hydroclimate import climate, hydrology

    provenance = json.loads((args.capture / "results.json").read_text())
    height, width = provenance["grid_shape"]
    latitude, _, _ = climate._build_grid(height, width)
    arms = {"weak": WEAK_ARMS, "mid": MID_ARMS, "dry": DRY_ARMS,
            "strengths": STRENGTH_ARMS}.get(args.profile, ARMS)
    dry_scale = DRY_SCALE_KG_M2 if args.profile == "dry" else None
    solvers = ("finite_volume", "wind_arrows") if args.solver == "both" else (args.solver,)
    results = {
        solver: _compare_capture(
            args.capture, args.maps, latitude, hydrology.FREEZE_POINT_C, arms,
            solver=solver, dry_scale=dry_scale,
        )
        for solver in solvers
    }
    print(json.dumps(results, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
