"""Explore uniform column caps and remove transport convergence in archived G3 controls.

Usage: backend/.venv/bin/python bin/debug/rfc0019_g3_archive_cap.py CAPTURE_DIR [--cfd-save SAVE]
       backend/.venv/bin/python bin/debug/rfc0019_g3_archive_cap.py BOUNDED_CAPTURE --bounded-capture [--maps DIR]
Reads diagnostic-control.npz and cfd-control.npz from the targeted-belt captures.
With --cfd-save, read only the saved CFD wind and temperature arrays to test a local cap.
With --bounded-capture, use both wind modes' fields captured during the paired cap screen.
Add --belt-screen to test two broad rainout preferences with full convergence adjustment.
Add --belt-transport to test two amounts of broad meridional transport with uniform rainout.
Add --overflow-placement to test spreading excess precipitation inside the same surface class.
"""

from __future__ import annotations

import argparse
import json
import sys
import zipfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import sparse
from scipy.ndimage import gaussian_filter1d
from scipy.sparse.linalg import spsolve

from rfc0019_coarse_moisture import (
    TransportFaces,
    aggregate_inputs,
    build_transport_rates,
    solve_bounded_budget,
    solve_budget,
)
import rfc0019_g3_compare as g3
from rfc0019_g3_bounded_column import _downslope_pairs, _pattern_metrics
from rfc0019_g3_mixing import GRID_SHAPE, RAIN_RATE, SOURCE_DEPTH_M, TUNING_SAVES
from rfc0019_g3_mixing import _coast_distances, _regional_metrics

BELT_WIDTH_DEG = 5.0
BELT_STRENGTHS = {"belt_015": 0.15, "belt_045": 0.45}
TRANSPORT_STRENGTHS = {"broad_transport_half": 0.5, "broad_transport_full": 1.0}
MAX_BROAD_WIND_M_S = 2.0
OVERFLOW_WIDTH_DEG = 4.0
OVERFLOW_RANGE_DEG = 8.0


def _rainout_belt_score(grid, faces):
    """Bounded latitude-dependent rainout cue from the unmodified wind."""
    nlat, nlon = grid.coarse_shape
    row = grid.coarse_cell // nlon
    signed = faces.forward - faces.reverse
    inflow = np.bincount(row[faces.second], weights=signed, minlength=nlat)
    outflow = np.bincount(row[faces.first], weights=signed, minlength=nlat)
    area = np.bincount(row, weights=grid.area, minlength=nlat)
    sigma_rows = BELT_WIDTH_DEG * nlat / 180.0
    smoothed = gaussian_filter1d((inflow - outflow) / area, sigma=sigma_rows, mode="nearest")
    rms = np.sqrt(np.average(smoothed**2, weights=area))
    if rms == 0:
        return np.zeros_like(grid.area)
    return np.tanh(smoothed[row] / rms)


def _broad_meridional_capacity(grid, faces):
    """Signed smooth transport contribution, in m²/s, on density-1 north faces only."""
    nlat, nlon = grid.coarse_shape
    row = grid.coarse_cell // nlon
    northward = row[faces.second] == row[faces.first] + 1
    latitude_face = row[faces.first[northward]]
    signed = (faces.forward - faces.reverse)[northward]
    edge_lat = np.linspace(-np.pi / 2, np.pi / 2, nlat + 1)[1:-1]
    circumference = 2 * np.pi * g3.EARTH_RADIUS_M * np.cos(edge_lat)
    zonal_flow = np.bincount(latitude_face, weights=signed, minlength=nlat - 1)
    sigma_rows = BELT_WIDTH_DEG * nlat / 180.0
    speed = gaussian_filter1d(zonal_flow / circumference, sigma=sigma_rows, mode="nearest")
    speed = np.clip(speed, -MAX_BROAD_WIND_M_S, MAX_BROAD_WIND_M_S)
    broad = np.zeros_like(faces.forward)
    broad[northward] = (speed * circumference / nlon)[latitude_face]
    return broad


def _belt_contrasts(fields, rain, area):
    metrics = {}
    for center in (-60, 0, 60):
        offset = np.abs(fields.lat_deg[:, None] - center)
        center_mask = fields.is_ocean & (offset < 8)
        flank_mask = fields.is_ocean & (offset >= 8) & (offset < 20)
        center_rain = g3._weighted_mean(rain[center_mask], area[center_mask])
        flank_rain = g3._weighted_mean(rain[flank_mask], area[flank_mask])
        metrics[str(center)] = {
            "center_mm_yr": center_rain,
            "flank_mm_yr": flank_rain,
            "ratio": center_rain / flank_rain if center_rain is not None and flank_rain else None,
        }
    return metrics


def _belt_row_shares(fields, rain, area):
    ocean_rain = np.sum(np.where(fields.is_ocean, rain * area, 0), axis=1)
    shares = {}
    for center in (-60, 0, 60):
        rows = ocean_rain[np.abs(fields.lat_deg - center) < 8]
        total = rows.sum()
        shares[str(center)] = float(rows.max() / total) if total > 0 else None
    return shares


def _distribute_overflow(overflow, area, surface_class):
    """Place each donor's kg/s on nearby receiving cells at the same longitude."""
    nlat = overflow.shape[0]
    row_spacing = 180.0 / nlat
    radius = int(OVERFLOW_RANGE_DEG / row_spacing)
    normalizer = np.zeros_like(overflow)
    links = []
    for offset in range(-radius, radius + 1):
        donor = slice(max(0, -offset), min(nlat, nlat - offset))
        receiver = slice(max(0, offset), min(nlat, nlat + offset))
        same_class = surface_class[donor] == surface_class[receiver]
        distance = offset * row_spacing
        weight = np.exp(-0.5 * (distance / OVERFLOW_WIDTH_DEG)**2) * area[receiver] * same_class
        normalizer[donor] += weight
        links.append((donor, receiver, weight))
    available = overflow / normalizer
    disposed = np.zeros_like(overflow)
    for donor, receiver, weight in links:
        disposed[receiver] += available[donor] * weight
    return disposed


def _surface_rain(budget, overflow, area, fields, freeze_point, spread):
    rain = budget.rain.reshape(area.shape)[::-1]
    closure = 0.0
    if spread:
        excess = overflow.reshape(area.shape)[::-1]
        surface_class = np.where(fields.is_ocean, 2, fields.air_temperature_c < freeze_point)
        deposition = _distribute_overflow(excess, area, surface_class)
        original = np.bincount(surface_class.ravel(), weights=excess.ravel(), minlength=3)
        received = np.bincount(surface_class.ravel(), weights=deposition.ravel(), minlength=3)
        closure = float(np.sum(np.abs(received - original)) / budget.rain.sum())
        if closure > 1e-12:
            raise AssertionError(f"overflow placement did not preserve surface class budgets: {closure:g}")
        rain = rain - excess + deposition
    return rain / area * g3.SECONDS_PER_YEAR, closure


def _signed_faces(first, second, signed):
    return TransportFaces(first, second, np.maximum(signed, 0), np.maximum(-signed, 0), np.zeros_like(signed))


def _convergence_correction(grid, first, second, capacity):
    zeros = np.zeros(grid.native_ids.shape)
    geometry = build_transport_rates(grid, zeros, zeros, 1.0)
    np.testing.assert_array_equal(first, geometry.first)
    np.testing.assert_array_equal(second, geometry.second)
    # The unit-diffusivity conductances give face length divided by cell spacing.
    # Use their metric weights to remove the divergent part of the air capacity.
    edges = np.arange(first.size)
    incidence = sparse.csr_matrix((
        np.concatenate((np.ones(first.size), -np.ones(first.size))),
        (np.concatenate((first, second)), np.tile(edges, 2)),
    ), shape=(grid.area.size, first.size))
    laplacian = incidence @ sparse.diags(geometry.exchange) @ incidence.T
    divergence = incidence @ capacity
    potential = np.append(spsolve(laplacian[:-1, :-1].tocsc(), divergence[:-1]), 0.0)
    correction = geometry.exchange * (incidence.T @ potential)
    residual = incidence @ (capacity - correction)
    tolerance = 1e-8 * np.max(abs(capacity))
    # The pinned Laplacian omits the last reservoir equation; its residual is the
    # accumulated residual of all solved rows, plus sparse-sum roundoff. #913
    np.testing.assert_allclose(residual[:-1], 0.0, atol=tolerance)
    pinned_bound = np.abs(residual[:-1]).sum() + tolerance
    if abs(residual[-1]) > pinned_bound:
        raise AssertionError("pinned convergence residual exceeds accumulated free-row error")
    return correction


def _compare(archive):
    area, mass, source, flux, first, second, sea = (archive[key] for key in (
        "area_m2", "mass_kg", "source_kg_s", "face_flux_kg_s", "face_first", "face_second", "is_ocean",
    ))
    column = mass / area
    upwind_column = np.where(flux >= 0, column[first], column[second])
    capacity = np.divide(flux, upwind_column, out=np.zeros_like(flux), where=upwind_column > 0)
    grid = aggregate_inputs(sea.reshape(GRID_SHAPE), {})
    sink = np.full_like(area, RAIN_RATE)
    control = solve_budget(area, source, sink, _signed_faces(first, second, capacity))
    np.testing.assert_allclose(control.mass, mass, rtol=1e-10, atol=1e-10)
    correction = _convergence_correction(grid, first, second, capacity)
    lat = np.repeat(np.linspace(-89, 89, GRID_SHAPE[0]), GRID_SHAPE[1])
    polar_land = (abs(lat) >= 60) & ~sea
    reference_polar_rain = control.rain[polar_land].sum()
    result = {"zero_inferred_donor_faces": int(np.count_nonzero(upwind_column == 0)), "arms": []}
    for removal, limit in ((0.0, 25), (0.0, 50), (0.5, 50), (1.0, 50)):
        faces = _signed_faces(first, second, capacity - removal * correction)
        unbounded = solve_budget(area, source, sink, faces)
        bounded, overflow, _ = solve_bounded_budget(area, source, sink, faces, limit * area)
        result["arms"].append({
            "convergence_removed_fraction": removal,
            "uniform_cap_kg_m2": limit,
            "unbounded_peak_column_kg_m2": float(np.max(unbounded.mass / area)),
            "overflow_share": float(overflow.sum() / source.sum()),
            "polar_land_rain_relative_change": float(bounded.rain[polar_land].sum() / reference_polar_rain - 1),
            "polar_land_rain_mm_yr": float(bounded.rain[polar_land].sum() / area[polar_land].sum()
                                           * g3.SECONDS_PER_YEAR),
            "land_rain_share": float(bounded.rain[~sea].sum() / source.sum()),
            "summed_cell_closure": bounded.summed_cell_relative,
        })
    return result


def _saved_cfd_fields(save, healpix_grid):
    with zipfile.ZipFile(save) as zipped:
        manifest = json.loads(zipped.read("manifest.json"))
        metadata = manifest["world"]["atmosphere_cfd_state"]["fields"]
        grid = healpix_grid.build(metadata["grid"]["__healpix_grid__"])
        fields = []
        for name in ("u", "v", "temperature_c"):
            with zipped.open(f"a/{metadata[name]['__npy__']}.npy") as stream:
                values = np.load(stream, allow_pickle=False)
            fields.append(healpix_grid.resample_to_equirect(grid, values, *GRID_SHAPE))
        return fields


def _compare_projection(fields, grid, source, expected_mass, reference, climate, freeze_point, *, screen=None):
    faces = build_transport_rates(grid, fields.wind_u[::-1], fields.wind_v[::-1], 0.0)
    sink = np.full_like(source, RAIN_RATE)
    control = solve_budget(grid.area, source, sink, faces)
    np.testing.assert_allclose(control.mass, expected_mass, rtol=1e-10, atol=1e-10)
    signed = faces.forward - faces.reverse
    correction = _convergence_correction(grid, faces.first, faces.second, signed)
    cap = (climate.saturation_specific_humidity_g_kg(fields.air_temperature_c[::-1].ravel())
           * g3.AIR_DENSITY_KG_M3 * SOURCE_DEPTH_M / 1000.0)
    distance = _coast_distances(fields)
    area = grid.native_area[::-1]
    polar_land = ~fields.is_ocean & (abs(fields.lat_deg[:, None]) >= 60)
    pairs = None if fields.elevation_m is None else _downslope_pairs(fields, g3._windward_leeward_masks(fields, area))
    result, maps = {}, {}
    projected = signed - correction
    arms = [("legacy", None, sink, False), ("unbounded", None, sink, False), ("capped", signed, sink, False),
            ("half_projection_capped", signed - 0.5 * correction, sink, False),
            ("full_projection_capped", projected, sink, False)]
    if screen == "rainout":
        score = _rainout_belt_score(grid, faces)
        arms.extend((name, projected, sink * (1 + strength * score), False)
                    for name, strength in BELT_STRENGTHS.items())
    if screen in ("transport", "placement"):
        broad = _broad_meridional_capacity(grid, faces)
        strengths = TRANSPORT_STRENGTHS if screen == "transport" else {
            "broad_transport_quarter": 0.25, "broad_transport_half": 0.5,
        }
        for name, strength in strengths.items():
            capacity = projected + strength * broad
            arms.append((name, capacity, sink, False))
            if screen == "placement":
                arms.append((name + "_spread", capacity, sink, True))
    for name, capacity, rate, spread in arms:
        budget = control
        overflow = np.zeros_like(source)
        if capacity is not None:
            adjusted = _signed_faces(faces.first, faces.second, capacity)
            budget, overflow, _ = solve_bounded_budget(grid.area, source, rate, adjusted, cap * grid.area)
        rain, placement_closure = _surface_rain(budget, overflow, area, fields, freeze_point, spread)
        if name == "legacy":
            rain = reference
        maps[name] = rain
        result[name] = _regional_metrics(fields, rain, area, distance, freeze_point) | {
            "polar_land_rain_mm_yr": float(np.average(rain[polar_land], weights=area[polar_land])),
            "raw_rms_to_legacy": g3._rms_pair(reference, rain, area)["raw_rel_rms"],
            "hemisphere_bands": g3._hemisphere_bands(fields, rain, area),
            "land_rain_share": float(np.sum((rain * area)[~fields.is_ocean]) / np.sum(rain * area)),
        }
        if pairs is not None:
            result[name] |= _pattern_metrics(fields, rain, area, distance, freeze_point, pairs)
        if screen is not None:
            result[name]["ocean_belt_contrasts"] = _belt_contrasts(fields, rain, area)
        if screen == "placement":
            result[name]["ocean_belt_wettest_row_shares"] = _belt_row_shares(fields, rain, area)
            result[name]["condensate_placement_relative_residual"] = placement_closure
        if name != "legacy":
            result[name] |= {
                "max_column_kg_m2": float(np.max(budget.mass / grid.area)),
                "overflow_share": float(overflow.sum() / source.sum()),
                "summed_cell_closure": budget.summed_cell_relative,
            }
    return result, maps


def _compare_saved_cfd(archive, save, climate, healpix_grid, freeze_point):
    u, v, temperature = _saved_cfd_fields(save, healpix_grid)
    sea = archive["is_ocean"].reshape(GRID_SHAPE)[::-1]
    grid = aggregate_inputs(sea[::-1], {})
    latitude, _, xyz = climate._build_grid(*GRID_SHAPE)
    fields = SimpleNamespace(is_ocean=sea, world_xyz=xyz, air_temperature_c=temperature, lat_deg=latitude,
                             wind_u=u, wind_v=v, elevation_m=None)
    reference = (archive["ref_kg_s"] / grid.area).reshape(GRID_SHAPE)[::-1] * g3.SECONDS_PER_YEAR
    result, _ = _compare_projection(fields, grid, archive["source_kg_s"], archive["mass_kg"], reference,
                                    climate, freeze_point)
    return result


def _compare_bounded_capture(capture, map_dir, climate, freeze_point, *, screen=None):
    provenance = json.loads((capture / "results.json").read_text())
    save = Path(provenance["save"])
    if TUNING_SAVES.get(save.name) != provenance["save_sha256"]:
        raise ValueError("expected a pinned tuning capture")
    latitude, _, xyz = climate._build_grid(*GRID_SHAPE)
    report = {}
    for wind in ("diagnostic", "cfd"):
        with np.load(capture / f"{wind}-columns.npz", allow_pickle=False) as archive:
            fields = SimpleNamespace(
                is_ocean=archive["is_ocean"], world_xyz=xyz, lat_deg=latitude,
                air_temperature_c=archive["air_temperature_c"], elevation_m=archive["elevation_m"],
                wind_u=archive["wind_u_m_s"], wind_v=archive["wind_v_m_s"],
            )
            grid = aggregate_inputs(fields.is_ocean[::-1], {})
            reference = archive["reference_precipitation_mm_yr"]
            report[wind], maps = _compare_projection(fields, grid, archive["source_kg_s"], archive["control_mass_kg"],
                                                     reference, climate, freeze_point, screen=screen)
            captured = (archive["bounded_rain_kg_s"] / grid.area).reshape(GRID_SHAPE)[::-1] * g3.SECONDS_PER_YEAR
            np.testing.assert_allclose(maps["capped"], captured, rtol=1e-10, atol=1e-10)
        if map_dir is not None:
            g3._save_common_scale_rain_map(map_dir, save, 1, wind, "coarse", "full-projection",
                                          reference, maps["full_projection_capped"], grid.native_area[::-1])
            for name in BELT_STRENGTHS.keys() & maps.keys():
                g3._save_common_scale_rain_map(
                    map_dir, save, 1, wind, "coarse", name, maps["full_projection_capped"], maps[name],
                    grid.native_area[::-1], labels=("Uniform rainout", "Broad rainout preference", "Change"),
                )
            for name in TRANSPORT_STRENGTHS.keys() & maps.keys():
                g3._save_common_scale_rain_map(
                    map_dir, save, 1, wind, "coarse", name, maps["full_projection_capped"], maps[name],
                    grid.native_area[::-1], labels=("Adjusted transport", "With broad moisture transport", "Change"),
                )
            for name in maps:
                if name.endswith("_spread"):
                    g3._save_common_scale_rain_map(
                        map_dir, save, 1, wind, "coarse", name, maps[name.removesuffix("_spread")], maps[name],
                        grid.native_area[::-1], labels=("Local overflow", "Spread overflow", "Change"),
                    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("capture", type=Path)
    options = parser.add_mutually_exclusive_group()
    options.add_argument("--cfd-save", type=Path)
    options.add_argument("--bounded-capture", action="store_true")
    parser.add_argument("--maps", type=Path)
    screens = parser.add_mutually_exclusive_group()
    screens.add_argument("--belt-screen", dest="screen", action="store_const", const="rainout")
    screens.add_argument("--belt-transport", dest="screen", action="store_const", const="transport")
    screens.add_argument("--overflow-placement", dest="screen", action="store_const", const="placement")
    args = parser.parse_args()
    if args.maps is not None and not args.bounded_capture:
        parser.error("--maps requires --bounded-capture")
    if args.screen is not None and not args.bounded_capture:
        parser.error("the belt screens require --bounded-capture")
    if args.cfd_save is not None or args.bounded_capture:
        repo = Path(__file__).resolve().parents[2]
        sys.path.insert(0, str(repo / "backend"))
        from app import healpix_grid
        from app.hydroclimate import climate, hydrology

    if args.bounded_capture:
        report = _compare_bounded_capture(args.capture, args.maps, climate, hydrology.FREEZE_POINT_C,
                                          screen=args.screen)
    elif args.cfd_save is not None:
        provenance = json.loads((args.capture / "results.json").read_text())
        digest = g3._sha256(args.cfd_save)
        if not digest == provenance["save_sha256"] == TUNING_SAVES.get(args.cfd_save.name):
            parser.error("expected the pinned tuning save used by the capture")
        with np.load(args.capture / "cfd-control.npz", allow_pickle=False) as archive:
            report = _compare_saved_cfd(archive, args.cfd_save, climate, healpix_grid, hydrology.FREEZE_POINT_C)
    else:
        report = {}
        for wind in ("diagnostic", "cfd"):
            with np.load(args.capture / f"{wind}-control.npz", allow_pickle=False) as archive:
                report[wind] = _compare(archive)
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
