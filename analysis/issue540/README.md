# Issue #540: RFC 0014 gate G3 — seasonal lake ice cover

## Question and preregistration

Test whether the sinusoidal ice-free fraction predicts observed complete lake-ice duration,
whether depth damping improves that prediction, and how the model's synthesized seasonal
temperature amplitude compares with ERA5. The related RFC acceptance check compares
`e_open(T_open) * f_free` with the existing `300 + 60*T` rate on #519's reference-lake table.

The observational window is 1990–2018, the common period covered by the G01377 records and the
ERA5 1990–2019 climatology. Include lakes with at least 20 valid ice-duration observations in
that window, valid coordinates, and positive reported mean depth. The depth requirement allows
the D2(a) damping fit. Fit `D_ref` and `a_min` only on this sample; additionally use leave-one-lake-
out fits to check whether any improvement generalizes. No simulation save or replay is used.

## Sources and reproduction

- Ice-on, ice-off and duration: Benson et al. (2000), [NSIDC G01377 Version 1](https://nsidc.org/data/g01377/versions/1), table [`liag_freeze_thaw_table.csv`](https://noaadata.apps.nsidc.org/NOAA/G01377/liag_freeze_thaw_table.csv).
- Lake depth: [G01377 physical-characteristics table](https://noaadata.apps.nsidc.org/NOAA/G01377/liag_physical_character_table.csv).
- Air temperature: [WeatherBench2 ERA5 hourly climatology](https://storage.googleapis.com/weatherbench2/datasets/era5-hourly-climatology/1990-2019_6h_240x121_equiangular_with_poles_conservative.zarr), 6-hourly by day-of-year at 1.5° resolution. For each lake, use the nearest grid cell; average the daily climatology into Gregorian monthly means after omitting February 29. The observed amplitude is half the warmest-minus-coldest monthly mean.
- Continentality: the matching WeatherBench2 ERA5 land-sea mask, passed through `biomes.grid_continentality`.

Input table SHA-256 values:

| Table | SHA-256 |
|---|---|
| G01377 freeze/thaw | `75bfd0fb1ad2f419e3487455c967d243d845cc4c14ef20836888b0869d4dae61` |
| G01377 physical characteristics | `0682de9103f35cb6f94e71ff495a45e2b1e68a1f884a42019916d1eb186b62ee` |

After `./bin/setup.sh`, install the driver's unpinned Zarr codec once in this worktree's venv.
The driver downloads the two source tables and caches ERA5 chunks under `/tmp`:

```bash
backend/.venv/bin/python -m pip install numcodecs
backend/.venv/bin/python bin/debug/gate_lake_ice_cover.py \
  --out analysis/issue540/g3.json
```

The gate ran on the issue branch based on `origin/main` `04d5c198b35d1951aa465c45a116fe066f3401b9`.
The output JSON contains the lake-level records, source checksums, fitted parameters and all
summary statistics. The source CSVs and downloaded ERA5 chunks remain outside the repository.

## Results

The sample has **54 lakes**, each with at least 20 observed-duration years in 1990–2018; the
median is 27 years. This meets the minimum of 50 lakes.

| Measure | Result |
|---|---:|
| Proxy median absolute error, no depth damping | **10.61 days** |
| Proxy median absolute error, damping fitted and measured on all lakes | 10.51 days |
| Leave-one-lake-out median absolute error, no damping | **10.61 days** |
| Leave-one-lake-out median absolute error, damping fit on the other 53 lakes | 11.28 days |
| Observed ERA5 seasonal amplitude, median | 12.08 °C |
| `biomes.seasonal_temp_amplitude`, median | 12.24 °C |
| Model amplitude median absolute error / bias | 2.08 °C / +0.75 °C |

The full-sample damping search reaches its `D_ref` upper bound (2,000 m) and `a_min = 0`; its
0.10-day improvement does not hold out of sample. **Keep amplitude damping off.** The undamped
ice-cover proxy meets the ≤30-day target by a wide margin, so G3(a) passes and the relation does
not need revision. G3(b) reports a small warm bias in the model-synthesized amplitude; the gate
sets no pass threshold for this comparison.

## Evaporation calibration check

The check now uses consistent climate inputs for the nine lakes in #519's H5 comparison:
monthly mean and half-range amplitude from the same 1990–2019 ERA5 climatology used by G3(a),
sampled at each lake's centroid. Centroids were calculated in EPSG:6933 from the matching
HydroLAKES v1 polygons and transformed to WGS84. The source archive is
[HydroLAKES polygons](https://www.hydrosheds.org/products/hydrolakes) at
`https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_polys_v10_shp.zip`, SHA-256
`9ef63498569dd7ccd0b8d8852da026e009a75f42d293f00e483f3b9fcfe8e1f9`. Coordinates, HydroLAKES
IDs, ERA5 grid cells, means, amplitudes, and per-lake results are in `g3.json`.

The observed rate ranges remain those transcribed in #519's H5 literature table. For Paiku Co,
the check uses the annual 870 ± 23 mm/yr estimate reported in
[Lei et al. (2023)](https://hess.copernicus.org/articles/27/4409/2023/hess-27-4409-2023.html).
This replaces the 975 ± 142 mm observation in
[Lei et al. (2021)](https://hess.copernicus.org/articles/25/3163/2021/), which covers only the
ice-free May–December period and is not comparable to the annual model rates.

For each lake, the reported change is the largest value over the observation-range endpoints of
`new relative absolute error − current relative absolute error`; a value above 0.10 would fail.
Checking both endpoints avoids relying on an arbitrary midpoint for a published range. All nine
remain within the allowance; the largest worsening is **3.91%** (Great Salt Lake). Thus the
evaporation calibration check passes against the published observations and stated uncertainty
ranges. The observations are heterogeneous in period and precision, as documented in #519; this
does not establish a universal evaporation law.

## Outcome

- G3(a): **pass** with observed ERA5 inputs, using the undamped amplitude; median absolute error
  is 10.61 days.
- Damping: **do not adopt**; its full-sample gain is negligible and its leave-one-lake-out error
  is worse.
- G3(b): model seasonal amplitudes are summarized above and lake by lake in `g3.json`.
- Evaporation check: **pass** across the nine published reference-lake observations; the worst
  endpoint increase in relative error is 3.91%.
- The G3 result is recorded in the RFC's G3 section; this report and `g3.json` provide its
  supporting provenance and lake-level results.
