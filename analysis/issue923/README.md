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
consolidated metadata checksum, method, metrics and validity checks. It streams at most one
two-day block into memory and writes no downloaded arrays to the repository. Raw reanalysis
chunks are not checked in.

## Outcome

The first full pass accidentally used the coarser 1.5° WeatherBench 2 store (`240×121`) instead
of the preregistered 0.25° store (`1440×721`). Its output is preserved as an exploratory result
only; do not apply the preregistered decision thresholds to it. The intended 0.25° source metadata
has now been verified, and the correct run is pending.

For transparency, the 1.5° exploratory pass gave an overall validation NRMSE of 0.580 and 12.7%
opposing-vector area for the monthly-mean proxy; its covariance control gave 0.569 and 11.8%.
At 1.5°, the primary proxy exceeded 10% opposing-vector area in all four land-season groups.
These values are not the preregistered result because source resolution is a material part of the
registered comparison.

The exploratory JSON and provenance are in
[`era5_transport_comparison_1p5_exploratory.json`](era5_transport_comparison_1p5_exploratory.json).
No model saves, G3 tuning candidates or held-out worlds were used.
