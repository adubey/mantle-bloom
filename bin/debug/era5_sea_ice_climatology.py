#!/usr/bin/env python3
"""Observed perennial and seasonal sea-ice extent by latitude from the ERA5 1990-2019
climatology: the reference for RFC 0015 gate G5 (#548).

Reads WeatherBench2's public ERA5 hourly climatology (1990-2019, 6-hourly x day-of-year,
1.5 deg conservative regrid; Rasp et al. 2024, doi:10.1029/2023MS004019; ERA5: Hersbach et al.
2020, doi:10.1002/qj.3803) chunk by chunk with `era5_surface_wind_climatology.py`'s reader
(#526). `sea_ice_cover` is ERA5's sea-ice area fraction. A cell is ocean where the land-sea
mask from the matching WeatherBench2 ERA5 store is below 0.5. The conservative regrid leaves a
coastal cell without a value wherever any of its source cells is land, even when the mask
counts it as ocean (9 million km2 north of 40 N). A second version fills those ocean cells
with the mean of their valid neighbours (3 x 3, wrapping in longitude), repeated up to
`COAST_FILL_PASSES` times.

The day-of-year climatology's daily mean (over its four 6-hourly slots) is the
concentration. Following NSIDC's extent convention, a cell is ice-covered on a day when its
concentration is at least 15%. A cell is **perennial** when it is ice-covered on every day of
the climatological year (ice survives the summer), and **seasonal** when it is ice-covered on
some days but not all, matching D6's model states. The same statistics as the model driver
(`sea_ice_extent_stats.hemisphere_summary`) are reported, plus each hemisphere's extent on its
minimum and maximum day, to check against NSIDC's Sea Ice Index (G02135). Both are reported
with and without the coastal fill: unfilled undercounts coastal ice and filled overcounts it,
so the two bracket the observed extents.

Needs numpy and numcodecs (Blosc), which the backend venv does not pin. Run it from any
environment that has them, e.g.:

    uv run --with numpy --with numcodecs python bin/debug/era5_sea_ice_climatology.py \
        --out analysis/issue548/era5_sea_ice_climatology.json

Downloads about 100 MB of compressed chunks and holds one 14 MB chunk at a time.
"""

from __future__ import annotations

import argparse
import json
import warnings

import numpy as np
from era5_surface_wind_climatology import CLIMATOLOGY, MASK_STORE, _array_meta, _decode
from sea_ice_extent_stats import cell_areas_km2, hemisphere_summary

EXTENT_THRESHOLD = 0.15
COAST_FILL_PASSES = 3


def _daily_concentration(store: str, name: str) -> np.ndarray:
    """(dayofyear, lon, lat) mean over the four 6-hourly slots."""
    meta = _array_meta(store, name)
    days = meta["shape"][1]
    day_chunk = meta["chunks"][1]
    out = np.full(meta["shape"][1:], np.nan, dtype=np.float32)
    for i in range(-(-days // day_chunk)):
        chunk = _decode(store, name, f"0.{i}.0.0", meta)
        n = min(day_chunk, days - i * day_chunk)
        with warnings.catch_warnings():  # all-land cells have no value on any day
            warnings.simplefilter("ignore", RuntimeWarning)
            out[i * day_chunk:i * day_chunk + n] = np.nanmean(chunk[:, :n], axis=0)
    return out


def _fill_coast(conc: np.ndarray, ocean: np.ndarray) -> np.ndarray:
    """(day, lat, lon) concentration with missing ocean cells set to the mean of their valid
    3 x 3 neighbours, for up to COAST_FILL_PASSES passes."""
    conc = conc.copy()
    for _ in range(COAST_FILL_PASSES):
        missing = np.isnan(conc) & ocean[None]
        if not missing.any():
            break
        padded = np.pad(conc, ((0, 0), (1, 1), (0, 0)), constant_values=np.nan)
        shifts = [np.roll(padded[:, 1 + di:padded.shape[1] - 1 + di], dj, axis=2) for di in (-1, 0, 1) for dj in (-1, 0, 1)]
        stack = np.stack(shifts)
        valid = np.isfinite(stack)
        mean = np.where(valid.any(axis=0), np.nansum(stack, axis=0) / np.maximum(valid.sum(axis=0), 1), np.nan)
        conc = np.where(missing, mean, conc)
    return conc


def _climatology(lat: np.ndarray, area: np.ndarray, ocean: np.ndarray, conc: np.ndarray) -> dict:
    """Each hemisphere's extent on its minimum and maximum day, and the perennial and seasonal
    statistics, for (day, lat, lon) concentration `conc`."""
    covered = np.nan_to_num(conc, nan=0.0) >= EXTENT_THRESHOLD
    perennial = ocean & covered.all(axis=0)
    seasonal = ocean & covered.any(axis=0) & ~perennial
    extents = {}
    for name, side in (("north", lat > 0), ("south", lat < 0)):
        daily = (covered & ocean & side[None, :, None]).astype(float)
        extent_mkm2 = (daily * area[None]).sum(axis=(1, 2)) / 1e6
        extents[name] = {
            "min_day_of_year": int(np.argmin(extent_mkm2)) + 1,
            "min_extent_mkm2": round(float(extent_mkm2.min()), 3),
            "max_day_of_year": int(np.argmax(extent_mkm2)) + 1,
            "max_extent_mkm2": round(float(extent_mkm2.max()), 3),
        }
    # The observed counterpart of D6's Polar Sea Ice class (f_ice >= 0.5).
    half_year = ocean & (covered.mean(axis=0) >= 0.5)
    half_year_mkm2 = {name: round(float(area[half_year & side[:, None]].sum()) / 1e6, 3) for name, side in (("north", lat > 0), ("south", lat < 0))}
    return {"daily_extent": extents, "ice_half_the_year_mkm2": half_year_mkm2,
            "states": hemisphere_summary(lat, area, ocean, perennial, seasonal)}


def _vector(store: str, name: str) -> np.ndarray:
    meta = _array_meta(store, name)
    return _decode(store, name, "0", meta).reshape(meta["shape"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    lat = _vector(CLIMATOLOGY, "latitude")
    conc = _daily_concentration(CLIMATOLOGY, "sea_ice_cover")  # (day, lon, lat)
    mask_meta = _array_meta(MASK_STORE, "land_sea_mask")
    land = _decode(MASK_STORE, "land_sea_mask", "0.0", mask_meta).reshape(mask_meta["shape"])  # (lon, lat)

    # Model-side statistics use (lat, lon) rows.
    ocean = land.T < 0.5
    conc = np.transpose(conc, (0, 2, 1))
    area = cell_areas_km2(lat, ocean.shape[1])
    filled = _fill_coast(conc, ocean)
    result = {
        "source": CLIMATOLOGY,
        "mask_source": MASK_STORE,
        "extent_threshold": EXTENT_THRESHOLD,
        "ocean_cells_without_value": int((np.isnan(conc).any(axis=0) & ocean).sum()),
        "ocean_cells_without_value_after_fill": int((np.isnan(filled).any(axis=0) & ocean).sum()),
        "unfilled": _climatology(lat, area, ocean, conc),
        "coast_filled": _climatology(lat, area, ocean, filled),
    }
    with open(args.out, "w") as f:
        json.dump(result, f, indent=1)
    print(json.dumps(result, indent=1))


if __name__ == "__main__":
    main()
