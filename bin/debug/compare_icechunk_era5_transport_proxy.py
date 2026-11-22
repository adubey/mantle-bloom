#!/usr/bin/env python3
"""Compare a column-water × 10 m wind proxy with public Icechunk ERA5 IVT (#923).

ERA5 fields are streamed anonymously from the Earthmover S3 archive. A static NOAA land-sea
mask is used only for land/ocean stratification because the selected ERA5 archive omits LSM.
See analysis/issue923/README.md for source caveats and the exploratory protocol amendment.
"""

from __future__ import annotations

import argparse
import calendar
from datetime import datetime, timedelta
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time

import icechunk
import numpy as np
import xarray as xr

from compare_era5_transport_proxy import (
    _accumulate_surface_transport_sums,
    _transport_proxy_validity,
    _aggregate_metrics,
    _cell_area,
    _conservative_weights,
    add_comparison_arguments,
    conservative_remap,
    pressure_integral,
    spherical_divergence,
    validate_comparison_years,
)


BUCKET = "earthmover-icechunk-era5"
PREFIX = "icechunkV2"
BRANCH = "main"
NOAA_LAND_MASK_URL = (
    "https://psl.noaa.gov/thredds/dodsC/Datasets/ncep.reanalysis/surface/land.nc"
)
DATASET_URL = "https://registry.opendata.aws/earthmover-era5/"
CHUNK_DAYS = 1
SAMPLE_HOURS = 6


def _open_era5(snapshot_id: str | None = None) -> tuple[icechunk.Session, xr.Dataset, xr.Dataset]:
    storage = icechunk.s3_storage(
        bucket=BUCKET,
        prefix=PREFIX,
        region="us-east-1",
        anonymous=True,
    )
    repository = icechunk.Repository.open(storage)
    if snapshot_id:
        session = repository.readonly_session(snapshot_id=snapshot_id)
    else:
        session = repository.readonly_session(BRANCH)
    options = {"consolidated": False, "chunks": None, "zarr_format": 3}
    pressure = xr.open_zarr(session.store, group="pressure/spatial", **options)
    single = xr.open_zarr(session.store, group="single/spatial", **options)
    return session, pressure, single


def _retry_read(load, label: str, attempts: int = 3):
    """Retry transient source reads briefly; deterministic selection errors fail immediately."""
    for attempt in range(attempts):
        try:
            return load()
        except (KeyError, ValueError):
            raise
        except Exception as exc:
            if attempt + 1 == attempts:
                raise
            delay_seconds = 2**attempt
            print(
                f"{label} failed ({exc}); retrying in {delay_seconds}s "
                f"({attempt + 2}/{attempts})",
                flush=True,
            )
            time.sleep(delay_seconds)
    raise RuntimeError(f"unreachable retry state for {label}")


def _load_land_fraction(
    target_lat: np.ndarray,
    target_lon: np.ndarray,
    cache_path: Path,
) -> tuple[np.ndarray, str]:
    """Use a cached remap or fetch NOAA's static mask for the archive's missing LSM field."""
    grid_hash = hashlib.sha256(
        np.ascontiguousarray(target_lat, dtype=np.float64).tobytes()
        + np.ascontiguousarray(target_lon, dtype=np.float64).tobytes()
    ).hexdigest()
    if cache_path.exists():
        try:
            with np.load(cache_path, allow_pickle=False) as cache:
                if set(cache.files) != {"land_fraction", "source_mask_sha256", "grid_sha256"}:
                    raise ValueError("unexpected fields")
                land_fraction = cache["land_fraction"]
                source_hash = str(cache["source_mask_sha256"].item())
                cached_grid_hash = str(cache["grid_sha256"].item())
        except Exception as exc:
            raise ValueError(f"invalid land-mask cache {cache_path}: {exc}") from exc
        if land_fraction.shape != (len(target_lat), len(target_lon)):
            raise ValueError(f"wrong land-mask cache shape in {cache_path}")
        if cached_grid_hash != grid_hash or not np.isfinite(land_fraction).all():
            raise ValueError(f"land-mask cache does not match the requested target grid: {cache_path}")
        _record_checkpoint_checksum(cache_path, cache_path.parent)
        return land_fraction, source_hash

    def fetch_mask() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        with xr.open_dataset(
            NOAA_LAND_MASK_URL,
            engine="pydap",
            decode_times=False,
        ) as source:
            land = source["land"].isel(time=0).load()
            return (
                np.asarray(source.lat.values, dtype=np.float64),
                np.asarray(source.lon.values, dtype=np.float64),
                np.asarray(land.values, dtype=np.float64),
            )

    source_lat, source_lon, land_values = _retry_read(fetch_mask, "NOAA land mask")
    mask_sha256 = hashlib.sha256(np.ascontiguousarray(land_values).tobytes()).hexdigest()
    weights = _conservative_weights(source_lat, source_lon, target_lat, target_lon)
    land_fraction = conservative_remap(land_values, weights)
    _atomic_npz(
        cache_path,
        land_fraction=land_fraction,
        source_mask_sha256=np.asarray(mask_sha256),
        grid_sha256=np.asarray(grid_hash),
    )
    _record_checkpoint_checksum(cache_path, cache_path.parent)
    return land_fraction, mask_sha256


def _month_starts(first_year: int, last_year: int) -> list[datetime]:
    return [
        datetime(year, month, 1)
        for year in range(first_year, last_year + 1)
        for month in range(1, 13)
    ]


def _next_month(start: datetime) -> datetime:
    if start.month == 12:
        return datetime(start.year + 1, 1, 1)
    return datetime(start.year, start.month + 1, 1)


def _get_chunk(
    pressure: xr.Dataset,
    single: xr.Dataset,
    start: datetime,
    end: datetime,
) -> tuple[xr.Dataset, xr.Dataset]:
    duration_hours = int((end - start).total_seconds() // 3600)
    if duration_hours <= 0 or duration_hours % SAMPLE_HOURS:
        raise ValueError(f"invalid chunk interval: {start} to {end}")
    expected_times = np.asarray(
        [
            np.datetime64(start + timedelta(hours=hour), "ns")
            for hour in range(0, duration_hours, SAMPLE_HOURS)
        ],
        dtype="datetime64[ns]",
    )

    def load_pair() -> tuple[xr.Dataset, xr.Dataset]:
        try:
            pressure_chunk = pressure[["q", "u", "v"]].sel(valid_time=expected_times)
            single_chunk = single[["tcwv", "u10", "v10", "sp"]].sel(
                valid_time=expected_times
            )
        except KeyError as exc:
            raise ValueError(
                f"archive is missing a required 6-hourly timestamp for {start:%Y-%m-%d}"
            ) from exc
        pressure_times = np.asarray(pressure_chunk.valid_time.values)
        single_times = np.asarray(single_chunk.valid_time.values)
        if not np.array_equal(pressure_times, expected_times) or not np.array_equal(
            single_times, expected_times
        ):
            raise ValueError(f"pressure and surface timestamps do not match for {start:%Y-%m-%d}")
        return pressure_chunk.load(), single_chunk.load()

    return _retry_read(load_pair, f"ERA5 chunk {start:%Y-%m-%d}")


def _array(dataset: xr.Dataset, name: str, dims: tuple[str, ...]) -> np.ndarray:
    return np.asarray(dataset[name].transpose(*dims).values, dtype=np.float32)


def _accumulate_month(
    start: datetime,
    pressure: xr.Dataset,
    single: xr.Dataset,
    pressure_levels_hpa: np.ndarray,
) -> tuple[np.ndarray, ...]:
    next_start = _next_month(start)
    shape = (pressure.sizes["latitude"], pressure.sizes["longitude"])
    ref_e_sum = np.zeros(shape, dtype=np.float64)
    ref_n_sum = np.zeros(shape, dtype=np.float64)
    tcwv_sum = np.zeros(shape, dtype=np.float64)
    u10_sum = np.zeros(shape, dtype=np.float64)
    v10_sum = np.zeros(shape, dtype=np.float64)
    tcwv_u10_sum = np.zeros(shape, dtype=np.float64)
    tcwv_v10_sum = np.zeros(shape, dtype=np.float64)
    ntime = 0
    one_level_extensions = 0

    for offset in range(0, calendar.monthrange(start.year, start.month)[1], CHUNK_DAYS):
        chunk_start = start + timedelta(days=offset)
        chunk_end = min(chunk_start + timedelta(days=CHUNK_DAYS), next_start)
        p_chunk, s_chunk = _get_chunk(pressure, single, chunk_start, chunk_end)
        chunk_count = int(p_chunk.sizes.get("valid_time", 0))
        if chunk_count == 0:
            continue

        q = _array(p_chunk, "q", ("valid_time", "pressure_level", "latitude", "longitude"))
        u = _array(p_chunk, "u", ("valid_time", "pressure_level", "latitude", "longitude"))
        v = _array(p_chunk, "v", ("valid_time", "pressure_level", "latitude", "longitude"))
        ps = _array(s_chunk, "sp", ("valid_time", "latitude", "longitude"))
        if np.nanmedian(ps) < 2_000.0:
            ps *= 100.0

        ivt_e, n_single_e = pressure_integral(q * u, pressure_levels_hpa, ps)
        ivt_n, n_single_n = pressure_integral(q * v, pressure_levels_hpa, ps)
        one_level_extensions += n_single_e + n_single_n
        ref_e_sum += np.nansum(ivt_e, axis=0)
        ref_n_sum += np.nansum(ivt_n, axis=0)

        tcwv = _array(s_chunk, "tcwv", ("valid_time", "latitude", "longitude"))
        u10 = _array(s_chunk, "u10", ("valid_time", "latitude", "longitude"))
        v10 = _array(s_chunk, "v10", ("valid_time", "latitude", "longitude"))
        _accumulate_surface_transport_sums(
            tcwv,
            u10,
            v10,
            tcwv_sum,
            u10_sum,
            v10_sum,
            tcwv_u10_sum,
            tcwv_v10_sum,
        )
        ntime += chunk_count
        del p_chunk, s_chunk, q, u, v, ps, ivt_e, ivt_n, tcwv, u10, v10

    if ntime == 0:
        raise ValueError(f"no samples for {start:%Y-%m}")
    return (
        ref_e_sum / ntime,
        ref_n_sum / ntime,
        tcwv_sum / ntime,
        u10_sum / ntime,
        v10_sum / ntime,
        tcwv_u10_sum / ntime,
        tcwv_v10_sum / ntime,
        np.asarray([ntime, one_level_extensions], dtype=np.int64),
    )


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as temporary:
        temporary_path = Path(temporary.name)
        json.dump(payload, temporary, indent=2, sort_keys=True)
        temporary.write("\n")
        temporary.flush()
        os.fsync(temporary.fileno())
    os.replace(temporary_path, path)


def _atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
    ) as temporary:
        temporary_path = Path(temporary.name)
        temporary.write(value)
        temporary.flush()
        os.fsync(temporary.fileno())
    os.replace(temporary_path, path)


def _atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="w+b", dir=path.parent, prefix=f".{path.name}.", suffix=".npz", delete=False
    ) as temporary:
        temporary_path = Path(temporary.name)
        np.savez_compressed(temporary, **arrays)
        temporary.flush()
        os.fsync(temporary.fileno())
    os.replace(temporary_path, path)


def _record_checkpoint_checksum(path: Path, checkpoint_dir: Path) -> None:
    """Maintain the capture directory's SHA256SUMS without rescanning old checkpoints."""
    checksum_file = checkpoint_dir / "SHA256SUMS"
    records: dict[str, str] = {}
    if checksum_file.exists():
        for line in checksum_file.read_text(encoding="utf-8").splitlines():
            try:
                digest, relative_path = line.split("  ", 1)
            except ValueError as exc:
                raise ValueError(f"invalid checksum entry in {checksum_file}: {line!r}") from exc
            records[relative_path] = digest

    relative_path = path.relative_to(checkpoint_dir).as_posix()
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    previous_digest = records.get(relative_path)
    if previous_digest is not None and previous_digest != digest:
        raise ValueError(f"checkpoint checksum changed for {relative_path}")
    records[relative_path] = digest
    lines = [f"{records[name]}  {name}" for name in sorted(records)]
    _atomic_text(checksum_file, "\n".join(lines) + "\n")


def _month_checkpoint_path(checkpoint_dir: Path, start: datetime) -> Path:
    return checkpoint_dir / "months" / f"{start:%Y-%m}.npz"


def _load_month_checkpoint(
    path: Path,
    shape: tuple[int, int],
    expected_samples: int,
) -> dict[str, np.ndarray]:
    required = {"reference", "proxy_monthly", "proxy_submonthly", "sample_count", "extensions"}
    try:
        with np.load(path, allow_pickle=False) as checkpoint:
            if set(checkpoint.files) != required:
                raise ValueError(f"unexpected fields: {sorted(checkpoint.files)}")
            values = {name: checkpoint[name] for name in required}
    except Exception as exc:
        raise ValueError(f"invalid monthly checkpoint {path}: {exc}") from exc
    for name in ("reference", "proxy_monthly", "proxy_submonthly"):
        if values[name].shape != (2, *shape):
            raise ValueError(f"invalid {name} shape in {path}: {values[name].shape}")
        if not np.isfinite(values[name]).all():
            raise ValueError(f"non-finite values in {name} checkpoint {path}")
    if values["sample_count"].shape != () or values["extensions"].shape != ():
        raise ValueError(f"invalid count fields in {path}")
    if int(values["sample_count"]) != expected_samples:
        raise ValueError(
            f"unexpected sample count in {path}: {int(values['sample_count'])} "
            f"(expected {expected_samples})"
        )
    return values


def _summary(
    *,
    first_year: int,
    last_year: int,
    train_end_year: int,
    session: icechunk.Session,
    pressure_levels_hpa: np.ndarray,
    remapped_land: np.ndarray,
    land_mask_sha256: str,
    target_lat: np.ndarray,
    area: np.ndarray,
    month_starts: list[datetime],
    references: list[np.ndarray],
    proxy_monthly: list[np.ndarray],
    proxy_submonthly: list[np.ndarray],
    sample_counts: list[int],
    extension_counts: list[int],
) -> dict:
    ref = np.stack(references)
    proxy_a = np.stack(proxy_monthly)
    proxy_b = np.stack(proxy_submonthly)
    months = np.asarray([item.month for item in month_starts])
    years = np.asarray([item.year for item in month_starts])
    train = years <= train_end_year
    expected_train_months = (train_end_year - first_year + 1) * 12
    training_complete = int(train.sum()) == expected_train_months

    alpha_num = np.sum(area[None, :, :] * np.sum(ref[train] * proxy_a[train], axis=1))
    alpha_den = np.sum(area[None, :, :] * np.sum(proxy_a[train] ** 2, axis=1))
    fitted_alpha = float(alpha_num / alpha_den) if alpha_den else None
    alpha = fitted_alpha if training_complete else None

    validation = years > train_end_year
    validation_metrics = None
    divergence_metrics = None
    if training_complete and np.any(validation) and alpha is not None:
        scaled_a = alpha * proxy_a[validation]
        scaled_b = alpha * proxy_b[validation]
        ref_validation = ref[validation]
        validation_months = months[validation]
        metrics_a = _aggregate_metrics(
            ref_validation, scaled_a, area, validation_months, remapped_land
        )
        metrics_b = _aggregate_metrics(
            ref_validation, scaled_b, area, validation_months, remapped_land
        )
        div_ref = np.stack(
            [spherical_divergence(vector[0], vector[1], target_lat) for vector in ref_validation]
        )
        div_a = np.stack(
            [
                spherical_divergence(vector[0], vector[1], target_lat)
                for vector in proxy_a[validation]
            ]
        ) * alpha
        area_month = np.broadcast_to(area, div_ref.shape)
        div_numerator = np.sum(area_month * (div_ref - div_a) ** 2)
        div_denominator = np.sum(area_month * div_ref**2)
        integrated_div_ref = float(np.sum(area_month * div_ref))
        integrated_div_proxy = float(np.sum(area_month * div_a))
        validation_metrics = {
            "monthly_mean_column_water_x_monthly_mean_10m_wind": metrics_a,
            "mean_6_hourly_column_water_x_10m_wind": metrics_b,
            "proxy_a_scaled_by_training_alpha": True,
            "proxy_b_scaled_by_same_frozen_training_alpha": True,
        }
        divergence_metrics = {
            "divergence_nrmse_proxy_a": float(np.sqrt(div_numerator / div_denominator))
            if div_denominator
            else None,
            "area_integrated_divergence_reference_kg_s": integrated_div_ref,
            "area_integrated_divergence_proxy_a_kg_s": integrated_div_proxy,
            "global_divergence_closure_reference_abs_kg_s": abs(integrated_div_ref),
            "global_divergence_closure_proxy_abs_kg_s": abs(integrated_div_proxy),
            "divergence_source_limit": (
                "No independent mass-consistent divergence or water-vapour tendency in this archive."
            ),
        }

    completed_months = len(month_starts)
    expected_months = (last_year - first_year + 1) * 12
    complete_years = sorted({item.year for item in month_starts if item.month == 12})
    result_status = (
        "exploratory_source_substitution"
        if completed_months == expected_months
        else "exploratory_partial_source_substitution"
    )
    return {
        "issue": 923,
        "status": result_status,
        "progress": {
            "completed_months": completed_months,
            "expected_months": expected_months,
            "completed_years": complete_years,
            "last_complete_year": max(complete_years) if complete_years else None,
            "training_period_complete": training_complete,
            "validation_months_available": int(validation.sum()) if training_complete else 0,
        },
        "source": {
            "dataset": "Earthmover ERA5 hourly analysis, public Icechunk/Zarr v3 on anonymous AWS S3",
            "dataset_url": DATASET_URL,
            "bucket": BUCKET,
            "prefix": PREFIX,
            "branch": BRANCH,
            "snapshot_id": str(session.snapshot_id),
            "period": f"{first_year}-{last_year}",
            "frequency_used": "6-hourly samples at 00, 06, 12 and 18 UTC from the hourly archive",
            "grid": "0.25 degree, 721x1440",
            "pressure_levels_hpa": pressure_levels_hpa.tolist(),
            "pressure_variables": ["q", "u", "v"],
            "single_level_variables": ["tcwv", "u10", "v10", "sp"],
            "license": "CC BY 4.0 per AWS Open Data Registry entry",
            "land_sea_mask": {
                "source": NOAA_LAND_MASK_URL,
                "grid": "2.5 degree regular latitude-longitude",
                "sha256_uncompressed_values": land_mask_sha256,
                "remapping": "conservative area-overlap remap to 2 degree comparison grid",
            },
            "mass_consistency": "standard ERA5; not the mass-consistent moisture-budget product",
        },
        "method": {
            "target_grid": "90x180, 2 degree regular latitude-longitude cell centers",
            "remapping": "separable first-order conservative area-overlap weights",
            "source_chunk_days": CHUNK_DAYS,
            "sample_interval_hours": SAMPLE_HOURS,
            "vertical_integration": (
                "trapezoid in pressure; exclude levels below surface pressure; linearly extend "
                "lowest two valid levels to surface pressure; zero integrand at 0 hPa above 50 hPa top"
            ),
            "train_years": f"{first_year}-{train_end_year}",
            "validation_years": f"{train_end_year + 1}-{last_year}",
            "alpha_fit_on_monthly_mean_proxy": alpha,
            "partial_fit_alpha_not_for_acceptance": fitted_alpha if not training_complete else None,
            "monthly_sample_counts_6h": sample_counts,
            "single_valid_level_surface_extensions_per_component_sample": sum(extension_counts),
            "surface_class_threshold_land_fraction": 0.5,
            "divergence": "spherical finite-volume from remapped cell-centered vectors",
        },
        "validation": validation_metrics or {
            "status": "waiting for complete training period and validation data",
            "months_available": 0,
        },
        "divergence": divergence_metrics,
        "validity": {
            **_transport_proxy_validity(
                ref, proxy_a, proxy_b, sample_counts, expected_months, remapped_land
            ),
            "global_divergence_closure_note": (
                "Finite-volume periodic closed-sphere divergence should integrate to roundoff; "
                "this is an implementation check, not an independent atmospheric moisture-budget validation."
            ),
        },
        "interpretation_limits": [
            "This public source substitution was not part of the issue's preregistered source amendment; report results as exploratory.",
            "The archive omits a land-sea mask, so land/ocean strata use a separate 2.5 degree NOAA static mask.",
            "This comparison does not test RFC 0019 mass-consistency correction or atmospheric budget closure.",
        ],
    }


def run(
    first_year: int,
    last_year: int,
    train_end_year: int,
    output: Path,
    checkpoint_dir: Path,
    interim_output: Path,
) -> None:
    manifest_path = checkpoint_dir / "manifest.json"
    previous_manifest = None
    requested_snapshot_id = None
    if manifest_path.exists():
        try:
            previous_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            raise ValueError(f"invalid checkpoint manifest {manifest_path}: {exc}") from exc
        if not isinstance(previous_manifest, dict):
            raise ValueError(f"checkpoint manifest is not an object: {manifest_path}")
        requested_snapshot_id = previous_manifest.get("snapshot_id")
        if not isinstance(requested_snapshot_id, str) or not requested_snapshot_id:
            raise ValueError(f"checkpoint manifest has no snapshot ID: {manifest_path}")

    session, pressure, single = _open_era5(requested_snapshot_id)
    try:
        required_pressure = {"q", "u", "v"}
        required_single = {"tcwv", "u10", "v10", "sp"}
        missing_pressure = required_pressure - set(pressure.data_vars)
        missing_single = required_single - set(single.data_vars)
        if missing_pressure or missing_single:
            raise ValueError(
                f"archive lacks required fields: pressure={sorted(missing_pressure)}, "
                f"single={sorted(missing_single)}"
            )

        pressure_levels = np.asarray(pressure.pressure_level.values, dtype=np.float64)
        pressure_order = np.argsort(pressure_levels)
        pressure_levels_hpa = pressure_levels[pressure_order]
        pressure = pressure.assign_coords(pressure_level=pressure_levels).isel(
            pressure_level=pressure_order
        )
        src_lat = np.asarray(pressure.latitude.values, dtype=np.float64)
        src_lon = np.asarray(pressure.longitude.values, dtype=np.float64)
        target_lat = np.arange(89.0, -90.0, -2.0)
        target_lon = np.arange(1.0, 360.0, 2.0)
        weights = _conservative_weights(src_lat, src_lon, target_lat, target_lon)
        area = _cell_area(target_lat, target_lon)
        remapped_land, land_mask_sha256 = _load_land_fraction(
            target_lat,
            target_lon,
            checkpoint_dir / "land-mask.npz",
        )

        manifest = {
            "checkpoint_schema": 2,
            "algorithm": "icechunk-era5-transport-v2",
            "snapshot_id": str(session.snapshot_id),
            "source_branch": BRANCH,
            "period": [first_year, last_year],
            "train_end_year": train_end_year,
            "pressure_levels_hpa": pressure_levels_hpa.tolist(),
            "pressure_variables": ["q", "u", "v"],
            "single_variables": ["tcwv", "u10", "v10", "sp"],
            "sample_interval_hours": SAMPLE_HOURS,
            "timestamp_selection": "explicit expected UTC timestamps; pressure and surface coordinates must match",
            "source_chunk_days": CHUNK_DAYS,
            "target_latitudes": target_lat.tolist(),
            "target_longitudes": target_lon.tolist(),
            "land_mask_sha256": land_mask_sha256,
        }
        if previous_manifest is not None:
            if previous_manifest != manifest:
                raise ValueError(
                    f"checkpoint configuration differs from {manifest_path}; "
                    "use a new checkpoint directory to start a separate run"
                )
            _record_checkpoint_checksum(manifest_path, checkpoint_dir)
        else:
            existing_months = list((checkpoint_dir / "months").glob("*.npz"))
            if existing_months:
                raise ValueError(
                    f"found monthly checkpoints without a manifest in {checkpoint_dir}; "
                    "preserve them and choose a new checkpoint directory"
                )
            _atomic_json(manifest_path, manifest)
            _record_checkpoint_checksum(manifest_path, checkpoint_dir)

        month_starts = _month_starts(first_year, last_year)
        references: list[np.ndarray] = []
        proxy_monthly: list[np.ndarray] = []
        proxy_submonthly: list[np.ndarray] = []
        processed_months: list[datetime] = []
        sample_counts: list[int] = []
        extension_counts: list[int] = []

        for start in month_starts:
            checkpoint_path = _month_checkpoint_path(checkpoint_dir, start)
            if checkpoint_path.exists():
                expected_samples = (
                    calendar.monthrange(start.year, start.month)[1] * 24 // SAMPLE_HOURS
                )
                month = _load_month_checkpoint(
                    checkpoint_path,
                    area.shape,
                    expected_samples,
                )
                _record_checkpoint_checksum(checkpoint_path, checkpoint_dir)
            else:
                accum = _accumulate_month(start, pressure, single, pressure_levels_hpa)
                (
                    ref_e,
                    ref_n,
                    mean_tcwv,
                    mean_u10,
                    mean_v10,
                    mean_tcwv_u10,
                    mean_tcwv_v10,
                    counts,
                ) = accum
                ntime, extensions = map(int, counts)
                ref_e = conservative_remap(ref_e, weights)
                ref_n = conservative_remap(ref_n, weights)
                proxy_a_e = conservative_remap(mean_tcwv * mean_u10, weights)
                proxy_a_n = conservative_remap(mean_tcwv * mean_v10, weights)
                proxy_b_e = conservative_remap(mean_tcwv_u10, weights)
                proxy_b_n = conservative_remap(mean_tcwv_v10, weights)
                month = {
                    "reference": np.stack((ref_e, ref_n)),
                    "proxy_monthly": np.stack((proxy_a_e, proxy_a_n)),
                    "proxy_submonthly": np.stack((proxy_b_e, proxy_b_n)),
                    "sample_count": np.asarray(ntime, dtype=np.int64),
                    "extensions": np.asarray(extensions, dtype=np.int64),
                }
                _atomic_npz(checkpoint_path, **month)
                _record_checkpoint_checksum(checkpoint_path, checkpoint_dir)

            references.append(month["reference"])
            proxy_monthly.append(month["proxy_monthly"])
            proxy_submonthly.append(month["proxy_submonthly"])
            processed_months.append(start)
            sample_counts.append(int(month["sample_count"]))
            extension_counts.append(int(month["extensions"]))

            if start.month == 12:
                result = _summary(
                    first_year=first_year,
                    last_year=last_year,
                    train_end_year=train_end_year,
                    session=session,
                    pressure_levels_hpa=pressure_levels_hpa,
                    remapped_land=remapped_land,
                    land_mask_sha256=land_mask_sha256,
                    target_lat=target_lat,
                    area=area,
                    month_starts=processed_months,
                    references=references,
                    proxy_monthly=proxy_monthly,
                    proxy_submonthly=proxy_submonthly,
                    sample_counts=sample_counts,
                    extension_counts=extension_counts,
                )
                _atomic_json(interim_output, result)
                validation = result["validation"]
                nrmse = validation.get("monthly_mean_column_water_x_monthly_mean_10m_wind", {}).get(
                    "overall", {}
                ).get("nrmse")
                print(
                    f"Processed ERA5 year {start.year}; "
                    f"saved {interim_output}; validation NRMSE={nrmse}",
                    flush=True,
                )

        result = _summary(
            first_year=first_year,
            last_year=last_year,
            train_end_year=train_end_year,
            session=session,
            pressure_levels_hpa=pressure_levels_hpa,
            remapped_land=remapped_land,
            land_mask_sha256=land_mask_sha256,
            target_lat=target_lat,
            area=area,
            month_starts=processed_months,
            references=references,
            proxy_monthly=proxy_monthly,
            proxy_submonthly=proxy_submonthly,
            sample_counts=sample_counts,
            extension_counts=extension_counts,
        )
        if len(processed_months) != (last_year - first_year + 1) * 12:
            raise ValueError("not all requested months have checkpoints; refusing to write final result")
        _atomic_json(output, result)
        print(
            json.dumps(
                {
                    "output": str(output),
                    "interim_output": str(interim_output),
                    "checkpoint_dir": str(checkpoint_dir),
                    "snapshot_id": str(session.snapshot_id),
                    "alpha": result["method"]["alpha_fit_on_monthly_mean_proxy"],
                    "validation_nrmse": result["validation"]
                    ["monthly_mean_column_water_x_monthly_mean_10m_wind"]["overall"]["nrmse"],
                },
                indent=2,
            )
        )
    finally:
        pressure.close()
        single.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    add_comparison_arguments(parser)
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        help="monthly restart checkpoints (defaults outside the repository under ~/mantle-bloom-captures)",
    )
    parser.add_argument(
        "--interim-out",
        type=Path,
        help="atomic partial summary updated after each completed year",
    )
    args = parser.parse_args()
    validate_comparison_years(args, parser)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    checkpoint_dir = args.checkpoint_dir or (
        Path.home() / "mantle-bloom-captures" / "issue923" / "icechunk-era5-checkpoints"
    )
    interim_output = args.interim_out or args.out.with_name(
        f"{args.out.stem}.interim{args.out.suffix}"
    )
    run(
        args.first_year,
        args.last_year,
        args.train_end_year,
        args.out,
        checkpoint_dir,
        interim_output,
    )


if __name__ == "__main__":
    main()
