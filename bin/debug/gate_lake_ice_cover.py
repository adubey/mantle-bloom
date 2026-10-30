#!/usr/bin/env python3
"""RFC 0014 gate G3: compare the sinusoidal lake ice-free fraction with Benson observations (#540).

The observation window is 1990-2018, the years shared by the Benson freeze/thaw table and
WeatherBench2's 1990-2019 ERA5 hourly climatology. Lakes need at least 20 observed ice-duration
years in that window and a positive reported mean depth. The latter is required for the damping
fit. ERA5 values are sampled at the nearest 1.5-degree grid cell. No simulation save is needed.

Usage after `./bin/setup.sh` (downloads public source tables and necessary ERA5 chunks):

    backend/.venv/bin/python -m pip install numcodecs
    backend/.venv/bin/python bin/debug/gate_lake_ice_cover.py --out analysis/issue540/g3.json

For offline reruns, pass `--freeze-thaw` and `--physical` with local copies of the two NSIDC CSVs.
The complete source tables and ERA5 chunks are cached under `/tmp` by default and are not checked
in. Run from the repository root.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import urllib.request
from collections import defaultdict
from pathlib import Path

import numpy as np
from numcodecs import Blosc

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "backend"))

FREEZE_THAW_URL = "https://noaadata.apps.nsidc.org/NOAA/G01377/liag_freeze_thaw_table.csv"
PHYSICAL_URL = "https://noaadata.apps.nsidc.org/NOAA/G01377/liag_physical_character_table.csv"
CLIMATOLOGY = (
    "https://storage.googleapis.com/weatherbench2/datasets/era5-hourly-climatology/"
    "1990-2019_6h_240x121_equiangular_with_poles_conservative.zarr"
)
MASK_STORE = (
    "https://storage.googleapis.com/weatherbench2/datasets/era5/"
    "1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr"
)
FIRST_YEAR, LAST_YEAR = 1990, 2018
MIN_YEARS = 20
T_ICE_C = 0.0
SECONDS_PER_DAY = 86400.0
SECONDS_PER_YEAR = 365.25 * SECONDS_PER_DAY
REFERENCE_LAKES = (
    # HydroLAKES v1 lake-polygon centroids (EPSG:6933 centroid transformed to WGS84).
    # H5's mean temperatures were approximate and supplied neither coordinates nor amplitudes;
    # the check now samples both climate inputs from the shared ERA5 climatology at these points.
    {"name": "Superior", "hylak_id": 5, "latitude": 47.517370, "longitude": -87.757306,
     "observed_mm_yr": [557.0, 557.0]},
    {"name": "Caspian", "hylak_id": 1, "latitude": 41.514007, "longitude": 50.760228,
     "observed_mm_yr": [940.0, 1010.0]},
    {"name": "Great Salt Lake", "hylak_id": 67, "latitude": 41.197830, "longitude": -112.539890,
     "observed_mm_yr": [1000.0, 1140.0]},
    {"name": "Aral (Large Aral Sea)", "hylak_id": 13, "latitude": 44.980349, "longitude": 59.475963,
     "observed_mm_yr": [1000.0, 1200.0]},
    {"name": "Balkhash", "hylak_id": 12, "latitude": 46.281739, "longitude": 75.619560,
     "observed_mm_yr": [1000.0, 1000.0]},
    # Annual 1980-2020 budget from Lei et al. (2023), replacing Lei et al. (2021)'s
    # 975 +/- 142 mm over the ice-free May-December period.
    {"name": "Paiku Co", "hylak_id": 1474, "latitude": 28.889601, "longitude": 85.585580,
     "observed_mm_yr": [847.0, 893.0]},
    {"name": "Turkana", "hylak_id": 158, "latitude": 3.547667, "longitude": 36.184329,
     "observed_mm_yr": [2335.0, 2335.0]},
    {"name": "Chad", "hylak_id": 15, "latitude": 13.471631, "longitude": 13.983256,
     "observed_mm_yr": [2000.0, 2000.0]},
    {"name": "Dead Sea", "hylak_id": 1432, "latitude": 31.533605, "longitude": 35.490009,
     "observed_mm_yr": [1200.0, 1600.0]},
)


def _fetch(url: str, destination: Path | None = None) -> bytes:
    if destination is not None and destination.exists():
        return destination.read_bytes()
    with urllib.request.urlopen(url, timeout=120) as response:
        data = response.read()
    if destination is not None:
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(data)
    return data


def _metadata(store: str, name: str, cache: Path) -> dict:
    target = cache / (hashlib.sha256(f"{store}/{name}/.zarray".encode()).hexdigest() + ".json")
    return json.loads(_fetch(f"{store}/{name}/.zarray", target))


def _chunk(store: str, name: str, key: str, meta: dict, cache: Path) -> np.ndarray:
    target = cache / (hashlib.sha256(f"{store}/{name}/{key}".encode()).hexdigest() + ".bin")
    raw = _fetch(f"{store}/{name}/{key}", target)
    decoded = Blosc().decode(raw)
    values = np.frombuffer(decoded, dtype=np.dtype(meta["dtype"]))
    return values.reshape(meta["chunks"])


def _vector(store: str, name: str, cache: Path) -> np.ndarray:
    meta = _metadata(store, name, cache)
    return _chunk(store, name, "0", meta, cache).reshape(meta["shape"])


def _era5_grid(cache: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (lat, lon, t2m in C [lat, lon], ERA5 continentality [lat, lon])."""
    from era5_surface_wind_climatology import _array_meta, _fetch as wind_fetch
    from app.hydroclimate import biomes

    lat = _vector(CLIMATOLOGY, "latitude", cache).astype(float)
    lon = _vector(CLIMATOLOGY, "longitude", cache).astype(float)
    mask_meta = _array_meta(MASK_STORE, "land_sea_mask")
    land_raw = Blosc().decode(wind_fetch(f"{MASK_STORE}/land_sea_mask/0.0"))
    land = np.frombuffer(land_raw, dtype=np.dtype(mask_meta["dtype"])).reshape(mask_meta["chunks"])
    ocean = land.T < 0.5
    continentality = biomes.grid_continentality(ocean)

    name = "2m_temperature"
    meta = _metadata(CLIMATOLOGY, name, cache)
    # Preserve the 366-day climatology, average the four daily samples, and then average calendar
    # days into months. Chunk order follows WeatherBench2's [hour, day, lon, lat] layout.
    day_count = int(meta["shape"][1])
    day_chunk = int(meta["chunks"][1])
    lon_count = int(meta["shape"][2])
    lon_chunk = int(meta["chunks"][2])
    daily = np.empty((day_count, lon_count, int(meta["shape"][3])), dtype=np.float32)
    for day_i in range(math.ceil(day_count / day_chunk)):
        for lon_i in range(math.ceil(lon_count / lon_chunk)):
            block = _chunk(CLIMATOLOGY, name, f"0.{day_i}.{lon_i}.0", meta, cache)
            nd = min(day_chunk, day_count - day_i * day_chunk)
            nl = min(lon_chunk, lon_count - lon_i * lon_chunk)
            daily[day_i * day_chunk:day_i * day_chunk + nd,
                  lon_i * lon_chunk:lon_i * lon_chunk + nl, :] = block[:, :nd].mean(axis=0)[:, :nl, :]

    # 2m_temperature is Kelvin. WeatherBench2's day-of-year climatology has leap day at index 59;
    # omit it when forming the twelve 365-day calendar-month means.
    no_leap = np.concatenate((daily[:59], daily[60:365]), axis=0)
    month_days = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    monthly = []
    offset = 0
    for nd in month_days:
        monthly.append(no_leap[offset:offset + nd].mean(axis=0))
        offset += nd
    monthly_c = np.stack(monthly, axis=0) - 273.15
    mean_c = monthly_c.mean(axis=0).T
    amplitude_c = (monthly_c.max(axis=0) - monthly_c.min(axis=0)).T / 2.0
    return lat, lon, mean_c, amplitude_c, continentality


def _clean_number(value: str) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) and number > 0.0 else None


def _read_lakes(freeze_thaw: Path, physical: Path) -> tuple[list[dict], dict]:
    durations: dict[str, list[int]] = defaultdict(list)
    with freeze_thaw.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            try:
                year, duration = int(row["iceon_year"]), int(row["duration"])
            except (KeyError, TypeError, ValueError):
                continue
            if row.get("lakeorriver") == "L" and FIRST_YEAR <= year <= LAST_YEAR and duration >= 0:
                durations[row["lakecode"]].append(duration)

    characteristics = {}
    with physical.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if row.get("lakeorriver") == "L":
                characteristics[row["lakecode"]] = row

    lakes = []
    for code, years in sorted(durations.items()):
        row = characteristics.get(code)
        if row is None or len(years) < MIN_YEARS:
            continue
        depth = _clean_number(row.get("mean_depth", ""))
        try:
            lat = float(row["lat_decimal"])
            lon = float(row["lon_decimal"])
        except (KeyError, TypeError, ValueError):
            continue
        if depth is None or not math.isfinite(lat) or not math.isfinite(lon):
            continue
        lakes.append({
            "code": code,
            "name": row.get("lakename", code),
            "latitude": lat,
            "longitude": lon,
            "mean_depth_m": depth,
            "duration_days": float(np.mean(years)),
            "observation_years": len(years),
        })
    return lakes, {"lake_codes_with_durations": len(durations), "selected_lakes": len(lakes)}


def _nearest(values: np.ndarray, x: float, periodic: bool = False) -> int:
    if periodic:
        delta = np.abs((values - x + 180.0) % 360.0 - 180.0)
    else:
        delta = np.abs(values - x)
    return int(np.argmin(delta))


def _free_fraction(mean_c: np.ndarray, amplitude_c: np.ndarray) -> np.ndarray:
    amp = np.maximum(amplitude_c, 1e-9)
    x = np.clip((T_ICE_C - mean_c) / amp, -1.0, 1.0)
    return np.arccos(x) / np.pi


def _fit_depth_damping(depth: np.ndarray, mean: np.ndarray, amplitude: np.ndarray, observed: np.ndarray) -> tuple[float, float, float]:
    d_refs = np.geomspace(1.0, 2000.0, 2000)
    a_mins = np.linspace(0.0, 1.0, 201)
    best = (float("inf"), 0.0, 0.0)
    for d_ref in d_refs:
        factors = np.maximum(1.0 - depth[:, None] / d_ref, a_mins[None, :])
        predicted = 365.25 * (1.0 - _free_fraction(mean[:, None], amplitude[:, None] * factors))
        errors = np.median(np.abs(predicted - observed[:, None]), axis=0)
        j = int(np.argmin(errors))
        candidate = (float(errors[j]), float(d_ref), float(a_mins[j]))
        if candidate < best:
            best = candidate
    return best[1], best[2], best[0]


def _mae(observed: np.ndarray, predicted: np.ndarray) -> float:
    return float(np.median(np.abs(observed - predicted)))


def _evaporation_check(lat: np.ndarray, lon: np.ndarray, means: np.ndarray, amplitudes: np.ndarray) -> dict:
    rows = []
    for lake in REFERENCE_LAKES:
        i = _nearest(lat, lake["latitude"])
        j = _nearest(lon, lake["longitude"], periodic=True)
        mean_c = float(means[i, j])
        amplitude = float(amplitudes[i, j])
        phi = float(np.arccos(np.clip((T_ICE_C - mean_c) / max(amplitude, 1e-9), -1.0, 1.0)))
        free = phi / np.pi
        open_mean = mean_c + amplitude * math.sin(phi) / phi if phi > 0.0 else T_ICE_C
        baseline = max(50.0, 300.0 + 60.0 * mean_c)
        seasonal = max(50.0, 300.0 + 60.0 * open_mean) * free
        obs_lo, obs_hi = lake["observed_mm_yr"]
        observed = (obs_lo + obs_hi) / 2.0
        worsening_by_endpoint = [
            (abs(seasonal - value) - abs(baseline - value)) / value
            for value in (obs_lo, obs_hi)
        ]
        worsening = max(worsening_by_endpoint)
        observed_worsening = (abs(seasonal - observed) - abs(baseline - observed)) / observed
        rows.append({
            "lake": lake["name"],
            "hydrolakes_id": lake["hylak_id"],
            "centroid_latitude": lake["latitude"],
            "centroid_longitude": lake["longitude"],
            "era5_grid_latitude": float(lat[i]),
            "era5_grid_longitude": float(lon[j]),
            "T_mean_c_era5": mean_c,
            "A_c_era5": amplitude,
            "f_free": free,
            "T_open_c": open_mean,
            "baseline_mm_yr": baseline,
            "seasonal_mm_yr": seasonal,
            "observed_mm_yr_range": [obs_lo, obs_hi],
            "relative_error_change_fraction_at_midpoint": observed_worsening,
            "relative_error_change_fraction": worsening,
            "within_10_percent_worsening": worsening <= 0.10,
        })
    return {
        "status": "complete_with_published_observation_ranges",
        "source_note": "T_mean and A are sampled from the WeatherBench2 ERA5 1990-2019 climatology at equal-area centroids of the matching HydroLAKES v1 polygons. Observed rate ranges are transcribed from #519's H5 literature table, except Paiku Co, which uses the annual 870 +/- 23 mm/yr estimate of Lei et al. (2023); the 975 +/- 142 mm value in Lei et al. (2021) covers only May-December.",
        "metric": "maximum over each published observed-rate range endpoint of (seasonal relative absolute error - current relative absolute error); allowed increase <= 0.10",
        "observation_range_note": "The gate's reference observations include published ranges and point estimates of mixed precision; range endpoints are tested so a midpoint choice cannot mask a >10% worsening.",
        "all_within_10_percent": all(row["within_10_percent_worsening"] for row in rows),
        "lakes": rows,
    }


def _leave_one_out(lakes: list[dict]) -> dict:
    obs = np.array([lake["duration_days"] for lake in lakes])
    mean = np.array([lake["air_mean_c"] for lake in lakes])
    amplitude = np.array([lake["air_amplitude_c"] for lake in lakes])
    depth = np.array([lake["mean_depth_m"] for lake in lakes])
    fixed = 365.25 * (1.0 - _free_fraction(mean, amplitude))
    cv_fixed, cv_damped = [], []
    for i in range(len(lakes)):
        keep = np.arange(len(lakes)) != i
        d_ref, a_min, _ = _fit_depth_damping(depth[keep], mean[keep], amplitude[keep], obs[keep])
        cv_fixed.append(fixed[i])
        cv_damped.append(365.25 * (1.0 - _free_fraction(
            np.array([mean[i]]), np.array([amplitude[i] * max(1.0 - depth[i] / d_ref, a_min)])
        ))[0])
    return {"without_damping_median_absolute_error_days": _mae(obs, np.array(cv_fixed)),
            "with_damping_median_absolute_error_days": _mae(obs, np.array(cv_damped)),
            "folds": len(lakes), "fit_scope": "each fold fits D_ref and a_min on the other lakes"}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--freeze-thaw", type=Path)
    parser.add_argument("--physical", type=Path)
    parser.add_argument("--cache", type=Path, default=Path("/tmp/mantle-bloom-g3-cache"))
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    freeze_path = args.freeze_thaw or args.cache / "liag_freeze_thaw_table.csv"
    physical_path = args.physical or args.cache / "liag_physical_character_table.csv"
    if not freeze_path.exists():
        _fetch(FREEZE_THAW_URL, freeze_path)
    if not physical_path.exists():
        _fetch(PHYSICAL_URL, physical_path)
    lakes, counts = _read_lakes(freeze_path, physical_path)
    if len(lakes) < 50:
        raise SystemExit(f"Only {len(lakes)} eligible lakes; G3 requires at least 50.")

    latitudes, longitudes, air_mean, air_amplitude, continentality = _era5_grid(args.cache / "era5")
    from app.hydroclimate import biomes

    for lake in lakes:
        i = _nearest(latitudes, lake["latitude"])
        j = _nearest(longitudes, lake["longitude"], periodic=True)
        lake["air_mean_c"] = float(air_mean[i, j])
        lake["air_amplitude_c"] = float(air_amplitude[i, j])
        lake["continentality"] = float(continentality[i, j])
        lake["model_amplitude_c"] = float(biomes.seasonal_temp_amplitude(
            np.array([lake["latitude"]]), np.array([lake["continentality"]]), 23.5
        )[0])

    observed = np.array([lake["duration_days"] for lake in lakes])
    means = np.array([lake["air_mean_c"] for lake in lakes])
    amplitudes = np.array([lake["air_amplitude_c"] for lake in lakes])
    depths = np.array([lake["mean_depth_m"] for lake in lakes])
    free = _free_fraction(means, amplitudes)
    undamped = 365.25 * (1.0 - free)
    d_ref, a_min, fitted_mae = _fit_depth_damping(depths, means, amplitudes, observed)
    damped_amplitude = amplitudes * np.maximum(1.0 - depths / d_ref, a_min)
    damped = 365.25 * (1.0 - _free_fraction(means, damped_amplitude))
    amp_residual = np.array([lake["model_amplitude_c"] - lake["air_amplitude_c"] for lake in lakes])
    leave_one_out = _leave_one_out(lakes)

    output = {
        "gate": "RFC 0014 G3 (#540)",
        "source": {
            "ice_observations": "Benson et al. (2000), NSIDC G01377, Version 1",
            "freeze_thaw_url": FREEZE_THAW_URL,
            "freeze_thaw_sha256": _sha256(freeze_path),
            "physical_characteristics_url": PHYSICAL_URL,
            "physical_characteristics_sha256": _sha256(physical_path),
            "air_temperature": "ERA5 1990-2019 hourly climatology from WeatherBench2, 6-hourly x day-of-year, 1.5-degree conservative grid",
            "temperature_store": CLIMATOLOGY,
            "reference_lake_locations": {
                "dataset": "HydroLAKES v1 polygon layer",
                "url": "https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip",
                "sha256": "9ef63498569dd7ccd0b8d8852da026e009a75f42d293f00e483f3b9fcfe8e1f9",
                "method": "EPSG:6933 polygon centroid transformed to WGS84; nearest ERA5 1.5-degree grid cell",
            },
            "ice_duration_observation_window": [FIRST_YEAR, LAST_YEAR],
            "grid_sampling": "nearest 1.5-degree grid cell; 1990-2019 monthly means from daily climatology; February 29 omitted",
            "lake_sample": f"lake records with >= {MIN_YEARS} valid duration years in {FIRST_YEAR}-{LAST_YEAR}, positive reported mean depth, and available coordinates",
        },
        "sample": counts | {"eligible_lakes": len(lakes), "minimum_observation_years": MIN_YEARS,
                           "observation_years_median": float(np.median([lake["observation_years"] for lake in lakes]))},
        "proxy": {
            "formula": "f_free = arccos(clip((T_ice - T_mean) / A, -1, 1)) / pi; predicted ice days = 365.25 * (1 - f_free)",
            "T_ice_c": T_ICE_C,
            "without_damping_median_absolute_error_days": _mae(observed, undamped),
            "with_damping_full_sample_median_absolute_error_days": _mae(observed, damped),
            "depth_damping_fit": {"D_ref_m": d_ref, "a_min": a_min, "fitted_median_absolute_error_days": fitted_mae,
                                  "search": "D_ref geometric grid 1-2000 m (2000 values), a_min 0-1 (201 values); median absolute error"},
            "leave_one_lake_out": leave_one_out,
            "damping_reduces_full_sample_error": _mae(observed, damped) < _mae(observed, undamped),
            "damping_reduces_leave_one_out_error": (
                leave_one_out["with_damping_median_absolute_error_days"]
                < leave_one_out["without_damping_median_absolute_error_days"]
            ),
            "adopt_damping": (
                leave_one_out["with_damping_median_absolute_error_days"]
                < leave_one_out["without_damping_median_absolute_error_days"]
            ),
        },
        "model_amplitude": {
            "function": "app.hydroclimate.biomes.seasonal_temp_amplitude",
            "continentality": "ERA5 1.5-degree land-sea mask passed to biomes.grid_continentality",
            "axial_tilt_deg": 23.5,
            "median_absolute_error_c": float(np.median(np.abs(amp_residual))),
            "median_bias_c": float(np.median(amp_residual)),
            "observed_amplitude_median_c": float(np.median(amplitudes)),
            "model_amplitude_median_c": float(np.median([lake["model_amplitude_c"] for lake in lakes])),
        },
        "evaporation_calibration": {
            **_evaporation_check(latitudes, longitudes, air_mean, air_amplitude),
        },
        "lakes": [{k: round(v, 5) if isinstance(v, float) else v for k, v in lake.items()} for lake in lakes],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=1) + "\n", encoding="utf-8")
    print(json.dumps({k: output[k] for k in ("sample", "proxy", "model_amplitude", "evaporation_calibration")}, indent=1))


if __name__ == "__main__":
    main()
