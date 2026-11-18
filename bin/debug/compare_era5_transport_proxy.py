#!/usr/bin/env python3
"""Compare a column-water × 10 m wind proxy with ERA5 pressure-level IVT (#923).

Uses the public WeatherBench 2 6-hourly ERA5 Zarr store. Data access is chunked by
calendar month; downloaded arrays are reduced immediately to the simulator's 90×180 grid.
See analysis/issue923/README.md for the preregistration, source limits and run command.
"""

from __future__ import annotations

import argparse
import calendar
from datetime import datetime, timedelta
import hashlib
import json
from pathlib import Path

import numpy as np
import xarray as xr


STORE = (
    "gs://weatherbench2/datasets/era5/"
    "1959-2023_01_10-wb13-6h-1440x721_with_derived_variables.zarr"
)
CHUNK_DAYS = 2
DATASET_URL = "https://weatherbench2.readthedocs.io/en/latest/data-guide.html"
G = 9.80665
EARTH_RADIUS_M = 6_371_000.0
LAT_BINS = np.arange(0.0, 91.0, 15.0)
SEASONS = {
    "DJF": (12, 1, 2),
    "MAM": (3, 4, 5),
    "JJA": (6, 7, 8),
    "SON": (9, 10, 11),
}
REQUIRED = (
    "specific_humidity",
    "u_component_of_wind",
    "v_component_of_wind",
    "surface_pressure",
    "total_column_water_vapour",
    "10m_u_component_of_wind",
    "10m_v_component_of_wind",
    "land_sea_mask",
)


def _lat_edges(centres: np.ndarray) -> np.ndarray:
    mids = (centres[:-1] + centres[1:]) / 2.0
    return np.concatenate(([-90.0], mids, [90.0]))


def _lon_edges(centres: np.ndarray) -> np.ndarray:
    mids = (centres[:-1] + centres[1:]) / 2.0
    step = float(np.median(np.diff(centres)))
    return np.concatenate(([centres[0] - step / 2.0], mids, [centres[-1] + step / 2.0]))


def _conservative_weights(
    source_lat: np.ndarray,
    source_lon: np.ndarray,
    target_lat: np.ndarray,
    target_lon: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return target-area-normalized latitude, longitude overlap matrices and sort order."""
    lat_order = np.argsort(source_lat)
    lat = source_lat[lat_order]
    se_lat = _lat_edges(lat)
    te_lat = _lat_edges(np.sort(target_lat))
    lat_weight = np.zeros((len(target_lat), len(lat)), dtype=np.float64)
    for ti, (south, north) in enumerate(zip(te_lat[:-1], te_lat[1:])):
        target_area = np.sin(np.deg2rad(north)) - np.sin(np.deg2rad(south))
        for si, (src_south, src_north) in enumerate(zip(se_lat[:-1], se_lat[1:])):
            lo, hi = max(south, src_south), min(north, src_north)
            if hi > lo:
                lat_weight[ti, si] = (
                    np.sin(np.deg2rad(hi)) - np.sin(np.deg2rad(lo))
                ) / target_area

    lon_order = np.argsort(source_lon)
    lon = source_lon[lon_order]
    se_lon = _lon_edges(lon)
    te_lon = _lon_edges(np.sort(target_lon))
    lon_weight = np.zeros((len(target_lon), len(lon)), dtype=np.float64)
    for ti, (west, east) in enumerate(zip(te_lon[:-1], te_lon[1:])):
        for si, (src_west, src_east) in enumerate(zip(se_lon[:-1], se_lon[1:])):
            overlap = 0.0
            for shift in (-360.0, 0.0, 360.0):
                lo = max(west, src_west + shift)
                hi = min(east, src_east + shift)
                overlap += max(0.0, hi - lo)
            lon_weight[ti, si] = overlap / (east - west)

    if not np.allclose(lat_weight.sum(axis=1), 1.0, atol=1e-10):
        raise ValueError("latitude remap weights do not cover every target cell")
    if not np.allclose(lon_weight.sum(axis=1), 1.0, atol=1e-10):
        raise ValueError("longitude remap weights do not cover every target cell")
    return lat_weight, lon_weight, (lat_order, lon_order)


def conservative_remap(
    field: np.ndarray,
    weights: tuple[np.ndarray, np.ndarray, tuple[np.ndarray, np.ndarray]],
) -> np.ndarray:
    lat_weight, lon_weight, (lat_order, lon_order) = weights
    ordered = np.take(np.take(field, lat_order, axis=-2), lon_order, axis=-1)
    # The weights use ascending latitude edges; return the simulator's north-to-south rows.
    return np.einsum("ay,...yx,bx->...ab", lat_weight, ordered, lon_weight, optimize=True)[..., ::-1, :]


def pressure_integral(
    integrand: np.ndarray,
    levels_hpa: np.ndarray,
    surface_pressure_pa: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Integrate q*u dp/g with below-ground masking and linear surface extension."""
    pressure_pa = np.asarray(levels_hpa, dtype=np.float64) * 100.0
    if np.any(np.diff(pressure_pa) <= 0):
        raise ValueError("pressure levels must be strictly increasing")
    valid = pressure_pa[None, :, None, None] <= surface_pressure_pa[:, None, :, :]
    count = valid.sum(axis=1)
    if np.any(count == 0):
        raise ValueError("a grid cell has no pressure levels above its surface")

    dp = np.diff(pressure_pa)[None, :, None, None]
    pairs = valid[:, :-1] & valid[:, 1:]
    result = np.sum(
        np.where(pairs, 0.5 * (integrand[:, :-1] + integrand[:, 1:]) * dp, 0.0),
        axis=1,
    )

    # Extend from the 50 hPa top level to zero using a zero integrand at 0 hPa.
    result += 0.5 * integrand[:, 0] * pressure_pa[0]

    last_idx = count - 1
    last_f = np.take_along_axis(integrand, last_idx[:, None], axis=1)[:, 0]
    last_p = pressure_pa[last_idx]
    above_idx = np.maximum(last_idx - 1, 0)
    prev_f = np.take_along_axis(integrand, above_idx[:, None], axis=1)[:, 0]
    prev_p = pressure_pa[above_idx]
    has_two = last_idx > 0
    slope = np.divide(
        last_f - prev_f,
        last_p - prev_p,
        out=np.zeros_like(last_f, dtype=np.float64),
        where=has_two,
    )
    bottom_dp = np.maximum(surface_pressure_pa - last_p, 0.0)
    surface_f = last_f + slope * bottom_dp
    result += 0.5 * bottom_dp * (last_f + surface_f)
    return result / G, int(np.count_nonzero(~has_two))


def _cell_area(latitude_desc: np.ndarray, longitude: np.ndarray) -> np.ndarray:
    lat_asc = np.sort(latitude_desc)
    edges = _lat_edges(lat_asc)
    d_sin = np.sin(np.deg2rad(edges[1:])) - np.sin(np.deg2rad(edges[:-1]))
    dlon = 2.0 * np.pi / len(longitude)
    return (EARTH_RADIUS_M**2 * d_sin[:, None] * dlon)[::-1, :]


def spherical_divergence(
    east: np.ndarray,
    north: np.ndarray,
    latitude_desc: np.ndarray,
) -> np.ndarray:
    """Finite-volume divergence on the periodic 2° sphere; cell values are center fluxes."""
    nlat, nlon = east.shape
    dlon = 2.0 * np.pi / nlon
    dlat = np.pi / nlat
    lat_edges = np.linspace(90.0, -90.0, nlat + 1)
    cos_edges = np.cos(np.deg2rad(lat_edges))

    east_face = 0.5 * (east + np.roll(east, -1, axis=1))
    west_face = 0.5 * (np.roll(east, 1, axis=1) + east)
    zonal_net = (east_face - west_face) * EARTH_RADIUS_M * dlat

    north_net = np.zeros_like(north, dtype=np.float64)
    for edge in range(1, nlat):
        face_north = 0.5 * (north[edge - 1] + north[edge])
        face_len = EARTH_RADIUS_M * max(cos_edges[edge], 0.0) * dlon
        north_net[edge - 1] -= face_north * face_len
        north_net[edge] += face_north * face_len
    area = _cell_area(latitude_desc, np.arange(nlon))
    return (zonal_net + north_net) / area


def _weighted_metrics(
    ref_e: np.ndarray,
    ref_n: np.ndarray,
    pred_e: np.ndarray,
    pred_n: np.ndarray,
    weights: np.ndarray,
) -> dict[str, float | None]:
    good = np.isfinite(ref_e) & np.isfinite(ref_n) & np.isfinite(pred_e) & np.isfinite(pred_n)
    w = np.where(good, weights, 0.0)
    total = float(w.sum())
    if total == 0.0:
        return {"nrmse": None, "opposing_area_fraction": None, "magnitude_bias": None, "mean_angle_deg": None, "signed_regression": None}
    re, rn = np.where(good, ref_e, 0.0), np.where(good, ref_n, 0.0)
    pe, pn = np.where(good, pred_e, 0.0), np.where(good, pred_n, 0.0)
    ref_sq = re * re + rn * rn
    pred_sq = pe * pe + pn * pn
    dot = re * pe + rn * pn
    denom = float(np.sum(w * ref_sq))
    rms = float(np.sqrt(np.sum(w * ((re - pe) ** 2 + (rn - pn) ** 2)) / total))
    nrmse = float(np.sqrt(np.sum(w * ((re - pe) ** 2 + (rn - pn) ** 2)) / denom)) if denom else None
    ref_mag, pred_mag = np.sqrt(ref_sq), np.sqrt(pred_sq)
    angle_valid = good & (ref_mag > 0.0) & (pred_mag > 0.0)
    angle_w = np.where(angle_valid, weights, 0.0)
    cos_angle = np.divide(
        dot,
        ref_mag * pred_mag,
        out=np.ones_like(dot, dtype=np.float64),
        where=angle_valid,
    )
    angle = np.degrees(np.arccos(np.clip(cos_angle, -1.0, 1.0)))
    angle_mean = float(np.sum(angle_w * angle) / angle_w.sum()) if angle_w.sum() else None
    mag_denom = float(np.sum(w * ref_mag))
    return {
        "nrmse": nrmse,
        "rmse_kg_m2_s": rms,
        "opposing_area_fraction": float(np.sum(w * (dot < 0.0)) / total),
        "magnitude_bias_fraction": float((np.sum(w * pred_mag) / mag_denom) - 1.0) if mag_denom else None,
        "mean_angle_deg": angle_mean,
        "signed_regression": float(np.sum(w * dot) / np.sum(w * pred_sq)) if np.sum(w * pred_sq) else None,
    }


def _aggregate_metrics(
    ref: np.ndarray,
    proxy: np.ndarray,
    weights: np.ndarray,
    months: np.ndarray,
    land_fraction: np.ndarray,
) -> dict[str, object]:
    # Arrays are month × component × latitude × longitude.
    result: dict[str, object] = {
        "overall": _weighted_metrics(ref[:, 0], ref[:, 1], proxy[:, 0], proxy[:, 1], weights[None, :, :]),
        "by_season": {},
        "by_surface_class": {},
        "by_absolute_latitude_band_deg": {},
    }
    for season, season_months in SEASONS.items():
        mask = np.isin(months, season_months)
        season_out: dict[str, object] = {}
        for surface, surface_mask in (
            ("land", land_fraction >= 0.5),
            ("ocean", land_fraction < 0.5),
        ):
            spatial_weight = weights * surface_mask
            season_out[surface] = _weighted_metrics(
                ref[mask, 0], ref[mask, 1], proxy[mask, 0], proxy[mask, 1], spatial_weight[None, :, :]
            )
        result["by_season"][season] = season_out

    for surface, surface_mask in (("land", land_fraction >= 0.5), ("ocean", land_fraction < 0.5)):
        result["by_surface_class"][surface] = _weighted_metrics(
            ref[:, 0], ref[:, 1], proxy[:, 0], proxy[:, 1], (weights * surface_mask)[None, :, :]
        )

    lat = np.arange(89.0, -90.0, -2.0)
    abs_lat = np.abs(lat)
    for low, high in zip(LAT_BINS[:-1], LAT_BINS[1:]):
        band_mask = (abs_lat[:, None] >= low) & (
            (abs_lat[:, None] < high) if high < 90.0 else (abs_lat[:, None] <= high)
        )
        result["by_absolute_latitude_band_deg"][f"{int(low)}-{int(high)}"] = _weighted_metrics(
            ref[:, 0], ref[:, 1], proxy[:, 0], proxy[:, 1], (weights * band_mask)[None, :, :]
        )

    group_metrics = []
    for season, season_months in SEASONS.items():
        time_mask = np.isin(months, season_months)
        for surface, surface_mask in (("land", land_fraction >= 0.5), ("ocean", land_fraction < 0.5)):
            metric = _weighted_metrics(
                ref[time_mask, 0], ref[time_mask, 1], proxy[time_mask, 0], proxy[time_mask, 1],
                (weights * surface_mask)[None, :, :],
            )
            group_metrics.append((season, surface, metric))
    result["broad_land_ocean_season_groups"] = {
        f"{season}_{surface}": metrics for season, surface, metrics in group_metrics
    }
    return result


def run(first_year: int, last_year: int, train_end_year: int, output: Path) -> None:
    fs = __import__("gcsfs").GCSFileSystem(token="anon")
    mapper = fs.get_mapper(STORE)
    metadata = mapper.get(".zmetadata")
    if metadata is None:
        raise ValueError("source store has no consolidated .zmetadata")
    metadata_sha256 = hashlib.sha256(metadata).hexdigest()
    with fs.open("weatherbench2/datasets/era5/LICENSE", "rb") as license_file:
        license_sha256 = hashlib.sha256(license_file.read()).hexdigest()
    ds = xr.open_zarr(mapper, consolidated=True, chunks={"time": 8})
    missing = [name for name in REQUIRED if name not in ds.variables]
    if missing:
        raise ValueError(f"source store lacks required variables: {missing}")
    selected = ds[list(REQUIRED)].sel(time=slice(f"{first_year}-01-01", f"{last_year}-12-31T23:59:59"))
    if selected.sizes.get("time", 0) == 0:
        raise ValueError("selected period has no timesteps")

    levels = np.asarray(ds["level"].values, dtype=np.float64)
    if levels.max() > 2_000.0:
        levels_hpa = levels / 100.0
    else:
        levels_hpa = levels
    src_lat = np.asarray(ds.latitude.values, dtype=np.float64)
    src_lon = np.asarray(ds.longitude.values, dtype=np.float64)
    target_lat = np.arange(89.0, -90.0, -2.0)
    target_lon = np.arange(1.0, 360.0, 2.0)
    weights = _conservative_weights(src_lat, src_lon, target_lat, target_lon)
    area = _cell_area(target_lat, target_lon)
    month_starts = [
        datetime(year, month, 1)
        for year in range(first_year, last_year + 1)
        for month in range(1, 13)
    ]

    references: list[np.ndarray] = []
    proxy_monthly: list[np.ndarray] = []
    proxy_submonthly: list[np.ndarray] = []
    monthly_div_ref: list[np.ndarray] = []
    monthly_div_proxy: list[np.ndarray] = []
    month_numbers: list[int] = []
    sample_counts: list[int] = []
    one_level_surface_extensions = 0
    remapped_land: np.ndarray | None = None

    for start in month_starts:
        if start.month == 12:
            next_start = datetime(start.year + 1, 1, 1)
        else:
            next_start = datetime(start.year, start.month + 1, 1)
        ref_e_sum = np.zeros((len(src_lat), len(src_lon)), dtype=np.float64)
        ref_n_sum = np.zeros_like(ref_e_sum)
        tcwv_sum = np.zeros_like(ref_e_sum)
        u10_sum = np.zeros_like(ref_e_sum)
        v10_sum = np.zeros_like(ref_e_sum)
        tcwv_u10_sum = np.zeros_like(ref_e_sum)
        tcwv_v10_sum = np.zeros_like(ref_e_sum)
        ntime = 0
        for offset in range(0, calendar.monthrange(start.year, start.month)[1], CHUNK_DAYS):
            chunk_start = start + timedelta(days=offset)
            chunk_end = min(chunk_start + timedelta(days=CHUNK_DAYS), next_start)
            chunk_ds = selected.sel(
                time=slice(chunk_start.isoformat(), (chunk_end - timedelta(seconds=1)).isoformat())
            ).load()
            chunk_count = int(chunk_ds.sizes.get("time", 0))
            if chunk_count == 0:
                continue
            q = np.asarray(
                chunk_ds["specific_humidity"].transpose("time", "level", "latitude", "longitude").values,
                dtype=np.float32,
            )
            u = np.asarray(
                chunk_ds["u_component_of_wind"].transpose("time", "level", "latitude", "longitude").values,
                dtype=np.float32,
            )
            v = np.asarray(
                chunk_ds["v_component_of_wind"].transpose("time", "level", "latitude", "longitude").values,
                dtype=np.float32,
            )
            ps = np.asarray(
                chunk_ds["surface_pressure"].transpose("time", "latitude", "longitude").values,
                dtype=np.float32,
            )
            if np.nanmedian(ps) < 2_000.0:
                ps *= 100.0
            ivt_e, n_single_e = pressure_integral(q * u, levels_hpa, ps)
            ivt_n, n_single_n = pressure_integral(q * v, levels_hpa, ps)
            one_level_surface_extensions += n_single_e + n_single_n
            ref_e_sum += np.nansum(ivt_e, axis=0)
            ref_n_sum += np.nansum(ivt_n, axis=0)

            tcwv = np.asarray(
                chunk_ds["total_column_water_vapour"].transpose("time", "latitude", "longitude").values,
                dtype=np.float32,
            )
            u10 = np.asarray(
                chunk_ds["10m_u_component_of_wind"].transpose("time", "latitude", "longitude").values,
                dtype=np.float32,
            )
            v10 = np.asarray(
                chunk_ds["10m_v_component_of_wind"].transpose("time", "latitude", "longitude").values,
                dtype=np.float32,
            )
            tcwv_sum += np.nansum(tcwv, axis=0)
            u10_sum += np.nansum(u10, axis=0)
            v10_sum += np.nansum(v10, axis=0)
            tcwv_u10_sum += np.nansum(tcwv * u10, axis=0)
            tcwv_v10_sum += np.nansum(tcwv * v10, axis=0)
            ntime += chunk_count
            del chunk_ds, q, u, v, ps, ivt_e, ivt_n, tcwv, u10, v10

        if ntime == 0:
            raise ValueError(f"no samples for {start:%Y-%m}")
        sample_counts.append(ntime)
        ref_e = conservative_remap(ref_e_sum / ntime, weights)
        ref_n = conservative_remap(ref_n_sum / ntime, weights)

        mean_tcwv = tcwv_sum / ntime
        mean_u10 = u10_sum / ntime
        mean_v10 = v10_sum / ntime
        p_a_e = mean_tcwv * mean_u10
        p_a_n = mean_tcwv * mean_v10
        p_b_e = tcwv_u10_sum / ntime
        p_b_n = tcwv_v10_sum / ntime
        proxy_a_e = conservative_remap(p_a_e, weights)
        proxy_a_n = conservative_remap(p_a_n, weights)
        proxy_b_e = conservative_remap(p_b_e, weights)
        proxy_b_n = conservative_remap(p_b_n, weights)

        references.append(np.stack((ref_e, ref_n)))
        proxy_monthly.append(np.stack((proxy_a_e, proxy_a_n)))
        proxy_submonthly.append(np.stack((proxy_b_e, proxy_b_n)))
        monthly_div_ref.append(spherical_divergence(ref_e, ref_n, target_lat))
        monthly_div_proxy.append(spherical_divergence(proxy_a_e, proxy_a_n, target_lat))
        month_numbers.append(start.month)

        if remapped_land is None:
            lsm = selected["land_sea_mask"]
            if "time" in lsm.dims:
                lsm = lsm.sel(time=f"{first_year}-01-01", method="nearest")
            remapped_land = conservative_remap(
                np.asarray(lsm.transpose("latitude", "longitude").values, dtype=np.float64), weights
            )

    ref = np.stack(references)
    proxy_a = np.stack(proxy_monthly)
    proxy_b = np.stack(proxy_submonthly)
    months = np.asarray(month_numbers)
    years = np.asarray([item.year for item in month_starts])
    if remapped_land is None:
        raise ValueError("land-sea mask was not loaded")
    train = years <= train_end_year
    train_weights = area[None, :, :]
    alpha_num = np.sum(train_weights * np.sum(ref[train] * proxy_a[train], axis=1))
    alpha_den = np.sum(train_weights * np.sum(proxy_a[train] ** 2, axis=1))
    alpha = float(alpha_num / alpha_den) if alpha_den else 0.0

    validation = years > train_end_year
    scaled_a = alpha * proxy_a
    scaled_b = alpha * proxy_b
    area_month = np.broadcast_to(area, (int(validation.sum()), *area.shape))
    metrics_a = _aggregate_metrics(ref[validation], scaled_a[validation], area, months[validation], remapped_land)
    metrics_b = _aggregate_metrics(ref[validation], scaled_b[validation], area, months[validation], remapped_land)

    div_ref = np.stack(monthly_div_ref)[validation]
    div_a = np.stack(monthly_div_proxy)[validation] * alpha
    div_error = np.sqrt(np.sum(area_month * (div_ref - div_a) ** 2) / np.sum(area_month * div_ref**2))
    integrated_div_ref = np.sum(area_month * div_ref)
    integrated_div_proxy = np.sum(area_month * div_a)
    result = {
        "issue": 923,
        "source": {
            "dataset": "WeatherBench 2 ERA5, public GCS Zarr, 6-hourly 13 pressure levels, 0.25 degree grid",
            "store": STORE,
            "dataset_guide": DATASET_URL,
            "period": f"{first_year}-{last_year}",
            "pressure_levels_hpa": levels_hpa.tolist(),
            "variables": list(REQUIRED),
            "consolidated_metadata_sha256": metadata_sha256,
            "license_sha256": license_sha256,
            "license_acceptance": "owner confirmed acceptance before data access",
            "frequency": "6-hourly",
            "mass_consistency": "standard ERA5; not the mass-consistent moisture-budget product",
        },
        "method": {
            "target_grid": "90x180, 2 degree regular latitude-longitude cell centers",
            "remapping": "separable first-order conservative area-overlap weights",
            "vertical_integration": "trapezoid in pressure; exclude levels below surface pressure; linearly extend lowest two valid levels to surface pressure; zero integrand at 0 hPa above 50 hPa top",
            "train_years": f"{first_year}-{train_end_year}",
            "validation_years": f"{train_end_year + 1}-{last_year}",
            "alpha_fit_on_monthly_mean_proxy": alpha,
            "monthly_sample_counts_6h": sample_counts,
            "single_valid_level_surface_extensions_per_component_sample": one_level_surface_extensions,
            "surface_class_threshold_land_fraction": 0.5,
            "divergence": "spherical finite-volume from remapped cell-centered vectors",
        },
        "validation": {
            "monthly_mean_column_water_x_monthly_mean_10m_wind": metrics_a,
            "mean_6_hourly_column_water_x_10m_wind": metrics_b,
            "proxy_a_scaled_by_training_alpha": True,
            "proxy_b_scaled_by_same_frozen_training_alpha": True,
            "divergence_nrmse_proxy_a": float(div_error),
            "area_integrated_divergence_reference_kg_s": float(integrated_div_ref),
            "area_integrated_divergence_proxy_a_kg_s": float(integrated_div_proxy),
            "divergence_source_limit": "No independent mass-consistent divergence or water-vapour tendency in this selected store.",
        },
        "validity": {
            "all_requested_months_present": len(sample_counts) == (last_year - first_year + 1) * 12,
            "minimum_samples_per_month": int(min(sample_counts)),
            "maximum_samples_per_month": int(max(sample_counts)),
            "global_divergence_closure_reference_abs_kg_s": float(abs(integrated_div_ref)),
            "global_divergence_closure_proxy_abs_kg_s": float(abs(integrated_div_proxy)),
            "finite_reference": bool(np.isfinite(ref).all()),
            "finite_monthly_mean_proxy": bool(np.isfinite(proxy_a).all()),
            "finite_submonthly_mean_proxy": bool(np.isfinite(proxy_b).all()),
            "land_fraction_min": float(np.nanmin(remapped_land)),
            "land_fraction_max": float(np.nanmax(remapped_land)),
            "global_divergence_closure_note": "Finite-volume periodic closed-sphere divergence should integrate to roundoff; this is an implementation check, not an independent atmospheric moisture-budget validation.",
        },
        "decision_thresholds": {
            "mismatch_prediction": "validation NRMSE > 0.25 or opposing vector area > 0.10 in at least two land/ocean-season groups",
            "evidence_against_mismatch": "validation NRMSE <= 0.15 and opposing vector area < 0.05 in every land/ocean-season group",
        },
    }
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    ds.close()
    print(json.dumps({"output": str(output), "alpha": alpha, "validation_nrmse": metrics_a["overall"]["nrmse"]}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--first-year", type=int, default=1979)
    parser.add_argument("--last-year", type=int, default=2022)
    parser.add_argument("--train-end-year", type=int, default=2000)
    parser.add_argument("--out", type=Path, required=True, help="compact JSON output path")
    args = parser.parse_args()
    if args.first_year >= args.last_year or not args.first_year <= args.train_end_year < args.last_year:
        parser.error("years must satisfy first-year <= train-end-year < last-year")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    run(args.first_year, args.last_year, args.train_end_year, args.out)


if __name__ == "__main__":
    main()
