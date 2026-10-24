#!/usr/bin/env python3
"""Zonal-mean 10 m wind speed over ocean and land from the ERA5 1990-2019 climatology, the
reference the diagnostic wind is calibrated against (#526).

Reads WeatherBench2's public ERA5 hourly climatology (1990-2019, 6-hourly x day-of-year,
1.5 deg conservative regrid; Rasp et al. 2024, doi:10.1029/2023MS004019; ERA5: Hersbach et al.
2020, doi:10.1002/qj.3803) straight from its public Zarr store over HTTPS, chunk by chunk, so
it needs no account and no xarray. `10m_wind_speed` is the climatological mean of the
instantaneous speed, which is what speed-linear consumers (evaporation, weathering) need; the
mean vector (u, v) is reported alongside it for comparison with the model's steady field.
The land-sea mask comes from the matching WeatherBench2 ERA5 store; a cell is ocean where the
mask is below 0.5.

Needs numpy and numcodecs (Blosc), which the backend venv does not pin. Run it from any
environment that has them, e.g.:

    uv run --with numpy --with numcodecs python bin/debug/era5_surface_wind_climatology.py \
        --out analysis/issue526/era5_wind_climatology.json

Downloads about 200 MB of compressed chunks and holds one 14 MB chunk at a time.
"""

from __future__ import annotations

import argparse
import json
import urllib.request

import numpy as np
from numcodecs import Blosc

CLIMATOLOGY = (
    "https://storage.googleapis.com/weatherbench2/datasets/era5-hourly-climatology/"
    "1990-2019_6h_240x121_equiangular_with_poles_conservative.zarr"
)
MASK_STORE = (
    "https://storage.googleapis.com/weatherbench2/datasets/era5/"
    "1959-2023_01_10-6h-240x121_equiangular_with_poles_conservative.zarr"
)
BAND_DEG = 10


def _fetch(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=120) as response:
        return response.read()


def _array_meta(store: str, name: str) -> dict:
    return json.loads(_fetch(f"{store}/{name}/.zarray"))


def _decode(store: str, name: str, key: str, meta: dict) -> np.ndarray:
    raw = Blosc().decode(_fetch(f"{store}/{name}/{key}"))
    return np.frombuffer(raw, dtype=np.dtype(meta["dtype"])).reshape(meta["chunks"])


def _time_mean(store: str, name: str) -> np.ndarray:
    """Mean over (hour, dayofyear) of a [4, 366, lon, lat] climatology variable."""
    meta = _array_meta(store, name)
    hours, days, n_lon, n_lat = meta["shape"]
    day_chunk = meta["chunks"][1]
    total = np.zeros((n_lon, n_lat))
    count = 0
    for i in range(-(-days // day_chunk)):
        block = _decode(store, name, f"0.{i}.0.0", meta).astype(np.float64)
        valid = min(day_chunk, days - i * day_chunk)
        total += block[:, :valid].sum(axis=(0, 1))
        count += hours * valid
    return total / count


def _area_mean(field: np.ndarray, weight: np.ndarray, mask: np.ndarray) -> float:
    return float((field * weight)[mask].sum() / weight[mask].sum())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out", help="write the compact JSON summary here")
    args = parser.parse_args()

    lat = np.frombuffer(Blosc().decode(_fetch(f"{CLIMATOLOGY}/latitude/0")), dtype="<f8")
    mask_meta = _array_meta(MASK_STORE, "land_sea_mask")
    land_fraction = _decode(MASK_STORE, "land_sea_mask", "0.0", mask_meta).astype(np.float64)
    speed = _time_mean(CLIMATOLOGY, "10m_wind_speed")
    u = _time_mean(CLIMATOLOGY, "10m_u_component_of_wind")
    v = _time_mean(CLIMATOLOGY, "10m_v_component_of_wind")

    lat_grid = np.broadcast_to(lat[None, :], speed.shape)
    weight = np.cos(np.radians(lat_grid))
    ocean = land_fraction < 0.5
    land = ~ocean
    everywhere = np.ones_like(ocean)

    bands = []
    for lo in range(-90, 90, BAND_DEG):
        in_band = (lat_grid >= lo) & (lat_grid < lo + BAND_DEG) if lo + BAND_DEG < 90 else (lat_grid >= lo)
        row = {"lat_lo": lo, "lat_hi": lo + BAND_DEG}
        for label, mask in (("ocean", ocean), ("land", land)):
            sel = in_band & mask
            if sel.any() and weight[sel].sum() > 0:
                row[f"{label}_speed_m_s"] = round(_area_mean(speed, weight, sel), 2)
                row[f"{label}_u_m_s"] = round(_area_mean(u, weight, sel), 2)
                row[f"{label}_v_m_s"] = round(_area_mean(v, weight, sel), 2)
        bands.append(row)

    summary = {
        "source": "ERA5 1990-2019 hourly climatology, WeatherBench2 240x121 conservative regrid",
        "global_speed_m_s": round(_area_mean(speed, weight, everywhere), 2),
        "ocean_speed_m_s": round(_area_mean(speed, weight, ocean), 2),
        "land_speed_m_s": round(_area_mean(speed, weight, land), 2),
        # Antarctica's katabatic plateau and Greenland dominate the high-latitude land mean;
        # the ice-free continents are the comparison for the model's land slowdown.
        "land_60s_to_70n_speed_m_s": round(_area_mean(speed, weight, land & (lat_grid > -60) & (lat_grid < 70)), 2),
        "bands": bands,
    }
    print(json.dumps(summary, indent=1))
    if args.out:
        with open(args.out, "w") as handle:
            json.dump(summary, handle, indent=1)
            handle.write("\n")


if __name__ == "__main__":
    main()
