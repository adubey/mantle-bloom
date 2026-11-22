# Issue #923: ERA5 transport proxy diagnostic

## Exploration summary

The preregistered Copernicus mass-consistent reference was replaced, at the project owner's
direction, with the public WeatherBench 2 ERA5 pressure-level subset. The owner confirmed
acceptance of the dataset's Copernicus licence before data access. The amendment is recorded
in [the issue discussion](https://github.com/adubey/mantle-bloom/issues/923#issuecomment-6907806350).

The preregistered WeatherBench 2 store contains six-hourly ERA5 on 13 pressure levels at 0.25°
resolution and ends on 2023-01-10. The comparison therefore uses complete years 1979–2022. It
tests the surface-wind proxy against standard ERA5 pressure-level IVT, not mass-consistent IVT.
It cannot validate RFC 0019's mass-consistency correction or atmospheric budget closure.

**Design principle 5** requires reporting column water, transport, and divergence in their
physical units and not labelling column water × 10 m wind as observed flux. **Design
principle 6** requires atmospheric water accounting and stationary closure ≤1e-8, but this
Earth-data comparison cannot establish that model invariant.

## Preregistered comparison

The owner-approved amendment, posted before data access, records the complete method and
thresholds on issue #923. In brief, compute 6-hourly IVT as the pressure integral of specific
humidity × pressure-level wind, mask levels below surface pressure, integrate by trapezoids,
extend the lowest two valid levels to surface pressure, and set the integrand to zero at 0 hPa
above the 50 hPa top. Conservatively remap vectors and proxy inputs to the simulator's 90×180
2° grid. Fit one global scalar to monthly-mean column water × monthly-mean 10 m wind on
1979–2000; freeze it for 2001–2022. The covariance control uses the monthly mean of
six-hourly column water × 10 m wind and the same frozen scalar.

The primary metric is area-weighted vector NRMSE. The preregistered mismatch screen is
validation NRMSE >0.25 or >10% opposing-vector area in at least two broad land/ocean-season
groups. Evidence against mismatch is NRMSE ≤0.15 and <5% opposing area in every such group.
Also report direction, magnitude bias, signed regression, seasons, land/ocean and 15°
absolute-latitude bands. Recompute spherical finite-volume divergence from the remapped
vectors; this is a proxy/reference comparison, not an independent budget validation.

## Reproduction

Install the analysis-only dependencies into the worktree's backend virtual environment:

```bash
backend/.venv/bin/pip install -r analysis/issue923/requirements.txt
```

Then run the chunked source-to-compact-summary diagnostic:

```bash
backend/.venv/bin/python bin/debug/compare_era5_transport_proxy.py \
  --out analysis/issue923/era5_transport_comparison.json
```

The script records the source store, relevant variable names, time coverage, pressure levels,
consolidated metadata checksum, method, metrics and validity checks. It streams one day of
six-hourly data at a time to bound the high-resolution pressure-level arrays' memory use; this
does not change the analysis period or any scientific calculation. It writes no downloaded
arrays to the repository. Raw reanalysis chunks are not checked in.

## Outcome

The first full pass accidentally used the coarser 1.5° WeatherBench 2 store (`240×121`) instead
of the preregistered 0.25° store (`1440×721`). Its output is preserved as an exploratory result
only; do not apply the preregistered decision thresholds to it. The intended 0.25° source metadata
was verified, but the full-resolution fetch was interrupted while xarray/dask was waiting on data
loads. The captured traceback ends in `KeyboardInterrupt`; it does not show a storage permission
exception, so the suspected access-control cause is not confirmed. No 0.25° result file was
written.

For transparency, the 1.5° exploratory pass gave an overall validation NRMSE of 0.580 and 12.7%
opposing-vector area for the monthly-mean proxy; its covariance control gave 0.569 and 11.8%.
At 1.5°, the primary proxy exceeded 10% opposing-vector area in all four land-season groups.
These values are not the preregistered result because source resolution is a material part of the
registered comparison.

The exploratory JSON and provenance are in
[`era5_transport_comparison_1p5_exploratory.json`](era5_transport_comparison_1p5_exploratory.json).
No model saves, G3 tuning candidates or held-out worlds were used.

## Open-source access probes

The NOAA/PSL NCEP/NCAR Reanalysis 1 one-day access probe succeeded over public THREDDS/OPeNDAP
for pressure humidity and winds, surface pressure, precipitable water, 10 m winds and land mask.
It is 4-times daily. The pressure-level humidity file has only eight levels (300–1000 hPa),
compared with 17 levels for wind, and the 10 m winds use a separate Gaussian grid. Using it would
therefore yield only a partial-column transport reference at a much coarser resolution than the
registered ERA5 comparison; the probe is not a scientific result, and no full NOAA comparison has
been started.

The next source probe succeeded against Earthmover's public Icechunk ERA5 archive on AWS S3 using
anonymous access. The pressure group exposes `q`, `u`, and `v` on 13 levels from 50 to 1000 hPa;
the single-level group exposes `tcwv`, `u10`, `v10`, and `sp`. A selected 2022-01-01 00 UTC sample
loaded for all seven variables at 0.25° resolution. The archive omits land-sea mask, so the new
diagnostic uses NOAA's static 2.5° NCEP/NCAR land-sea mask for the land/ocean strata. An initial
full pass reached 1991, then its execution session disappeared without writing a result file. It
had no intermediate checkpoints, so the completed years could not be recovered or summarized.

This source substitution happened after the issue's registered source amendment. Results from it
will be marked exploratory, and the original ERA5-specific decision thresholds will not be used as
a preregistered acceptance decision. The comparison still tests the vector proxy against standard
ERA5 pressure-level IVT; it does not test mass-consistent IVT or atmospheric budget closure.

## Icechunk ERA5 comparison

Install the updated analysis-only reader stack:

```bash
backend/.venv/bin/pip install -r analysis/issue923/requirements.txt
```

Run the source-substituted, 6-hourly comparison:

```bash
backend/.venv/bin/python bin/debug/compare_icechunk_era5_transport_proxy.py \
  --out analysis/issue923/era5_transport_comparison_icechunk_exploratory.json \
  --checkpoint-dir ~/mantle-bloom-captures/issue923/icechunk-era5-checkpoints
```

The script streams one day at a time from the `single/spatial` and `pressure/spatial` groups,
selects 00, 06, 12 and 18 UTC, conservatively remaps monthly reference/proxy fields to 90×180,
and records the Icechunk snapshot ID. The separate NOAA land mask is conservatively remapped to
that comparison grid. After each month it atomically saves compressed 2° checkpoint arrays outside
the repository. A manifest pins the checkpoint set to the Icechunk snapshot, date range, training
split, grid, variables, sample interval, and land-mask checksum. Restarting with the same arguments
skips completed months; incompatible checkpoints fail closed and are left untouched. After each
December it atomically refreshes the compact `*.interim.json` summary. Through 2000 it reports
training progress but no validation score; from 2001 onward it reports cumulative validation
metrics using the frozen 1979–2000 fit. The final JSON is written only after all requested months
are present. The checkpoints hold remapped monthly diagnostic fields, not raw ERA5 source arrays.

Source documentation:

- [NOAA/PSL NCEP/NCAR Reanalysis 1 overview](https://www.psl.noaa.gov/data/gridded/data.ncep.reanalysis.html)
- [NOAA/PSL pressure-level data catalog](https://psl.noaa.gov/thredds/catalog/Datasets/ncep.reanalysis/pressure/catalog.html)
- [NOAA/PSL surface data catalog](https://psl.noaa.gov/thredds/catalog/Datasets/ncep.reanalysis/surface/catalog.html)
- [Earthmover ERA5 on AWS Open Data Registry](https://registry.opendata.aws/earthmover-era5/)
