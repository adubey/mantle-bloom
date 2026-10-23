# Issue #526: a physical 10 m diagnostic wind

The diagnostic wind (`climate.compute_wind`, the default `World.wind_model == "diagnostic"`)
had no physical speed. Its thermal term was measured per grid cell, so it changed with
`climate_density` (principles 4 and 5). This branch:

- calibrates the planetary wind to ERA5's 10 m mean speed;
- measures the thermal term per km across a fixed ground distance;
- defines the surface-wind contract (`climate.surface_wind`);
- recalibrates the two speed-linear consumers.

## ERA5 reference

[`era5_wind_climatology.json`](era5_wind_climatology.json) is produced by
[`bin/debug/era5_surface_wind_climatology.py`](../../bin/debug/era5_surface_wind_climatology.py).
Its input is the ERA5 1990–2019 hourly climatology (Hersbach et al. 2020,
doi:10.1002/qj.3803) as regridded to 1.5° by WeatherBench2 (Rasp et al. 2024,
doi:10.1029/2023MS004019). `10m_wind_speed` is the climatological mean of the instantaneous
speed. A cell counts as ocean where the ERA5 land-sea mask is below 0.5.

| \|lat\| band | ERA5 ocean S / N (m/s) | model before S / N | model after S / N |
|---|---|---|---|
| 0–10 | 5.64 / 5.30 | 6.07 / 6.11 | 5.63 / 5.66 |
| 10–20 | 6.72 / 6.64 | 6.53 / 6.64 | 6.44 / 6.56 |
| 20–30 | 6.68 / 6.42 | 7.44 / 7.54 | 6.62 / 6.71 |
| 30–40 | 7.50 / 7.13 | 8.52 / 8.57 | 7.29 / 7.33 |
| 40–50 | 9.66 / 8.09 | 9.49 / 9.59 | 8.54 / 8.63 |
| 50–60 | 10.57 / 8.16 | 10.39 / 10.61 | 8.90 / 9.09 |
| 60–70 | 8.62 / 7.34 | 11.33 / 11.47 | 7.98 / 8.07 |
| 70–80 | 6.85 / 6.28 | 11.98 / 12.06 | 6.77 / 6.79 |
| 80–90 | — / 5.89 | 12.53 / 12.47 | 6.35 / 6.07 |
| ocean mean | 7.25 | 8.16 | 7.02 |
| land mean | 3.46 (3.02 ice-free, 60°S–70°N) | 3.82 | 3.07 |

The model columns are the diagnostic wind at density 4 on save A (below). Its rows are 10°
bands of ocean cells. The model is symmetric about the equator, so its envelope
(`climate.OCEAN_WIND_SPEED_M_S`) averages ERA5's two hemispheres. The Southern Ocean's
9.7–10.6 m/s and the northern storm tracks' 8.1–8.2 m/s straddle it.

The mean scalar speed is several times the length of ERA5's mean vector in the mid-latitudes:
at 30–40° the mean vector is 1.2–1.3 m/s and the mean speed 7.1–7.5 m/s. The contract therefore
defines `hypot(wind_u, wind_v)` as the mean speed and the vector's direction as the prevailing
one. It is not a mean vector.

Direction is unchanged. The trades here run mostly meridional
(`CORIOLIS_DEFLECTION_GAIN = 1.8`). ERA5's trades have u/v ≈ 2.4 at 10–20°, which implies a
gain near 9. The gain is kept because the humidity sweeps' zonal/meridional weighting is
tuned against it.

## Resolution

These are land and ocean mean 10 m speeds (m/s) on save A, diagnostic wind, for the same
world at each `climate_density`:

| density | land before | land after | ocean after | thermal term on land, before | after |
|---|---|---|---|---|---|
| 0.5 | — | 3.08 | 7.02 | — | 0.68 |
| 1 | 3.80 | 2.98 | 7.01 | 2.08 | 0.86 |
| 2 | 3.81 | 3.03 | 7.02 | 1.35 | 0.86 |
| 4 | 3.82 | 3.07 | 7.02 | 0.86 | 0.86 |
| 8 | — | 3.09 | 7.02 | — | 0.86 |

Both "after" columns are measured at the final coefficient of 29 (m/s)/(°C/km). The before column's mean speeds barely moved, but its thermal
term fell by a factor of 2.4 from density 1 to 4. A per-km neighbour difference fixed that
only partly, because finer grids resolve steeper lapse-rate contrasts on mountain flanks.
The term is now differenced across a fixed ±2° span. From density 1 to 8 the land mean varies
by 3.8%; the rest comes from terrain smoothing that is still counted in cells. At density 0.5
the span is one 4° cell either side. `test_wind_speed_is_independent_of_climate_density` pins
the 5% tolerance on synthetic terrain, and `test_thermal_wind_term_is_per_km_of_ground_distance`
fails on `origin/main`.

## Downstream calibration

These are one-shot `compute_climate` comparisons against `origin/main` (`0fc4901`) on the
same world state. Save A is seed 628977350 at 226.7 Myr, SHA-256 `4e1e2ee3…b32189`. Save B is
seed 910211954 at 220.2 Myr, SHA-256 `bd487409…b060d2`.

| Consumer | How it reads the wind | Change |
|---|---|---|
| Weathering | linear in speed | `WEATHERING_COEFFICIENT` 3.0 → 3.6. A proxy (speed × humidity × relief over land) fell to 0.82× on A and 0.86× on B at density 4. |
| Moisture-flux convergence | linear in speed / reference speed | `MERIDIONAL_BASE_SPEED` (6.0) replaced as normaliser by `MFC_REFERENCE_WIND_SPEED_M_S` = 4.9. Land precipitation moves −1.2% (A) and +0.4% (B) at density 4, and +1.4% (A) and +2.9% (B) at density 2. |
| Ocean currents | linear (Ekman + swirl) | None. Mean current falls 13% with the wind, but nothing reads current magnitude: temperature advection uses direction, arrows scale by their maximum and swell markers are weighted samples. |
| Humidity | direction share (`abs(u) / (abs(u) + abs(v))`) | None. The direction is unchanged. |

With these, the one-shot land-biome map agrees with `origin/main` on 98–99% (A) and 97–98% (B)
of land cells.

## Agreement with the CFD

[`cfd_agreement.json`](cfd_agreement.json) comes from
[`bin/debug/measure_wind_cfd_agreement.py`](../../bin/debug/measure_wind_cfd_agreement.py).
It uses seed 7 at density 2 with 10 × 2 Myr steps, as in
`test_diagnostic_wind_model_tracks_the_cfd_biome_map`, but with the CFD world actually set to
`"cfd"`. The stress test doesn't set it, and `"diagnostic"` has been the default since
`dc8bb8b` (#529).

| | origin/main | branch |
|---|---|---|
| wind-speed rel-RMS | 0.731 | 0.700 |
| land-biome agreement | 0.606 | 0.620 |
| precipitation rel-RMS | 0.528 | 0.522 |
| CFD ocean mean (m/s) | 5.18 | 4.30 |
| diagnostic ocean mean (m/s) | 8.93 | 7.38 |

The CFD settles at 0.4–0.65× its own relaxation target (#530), and the stress test's floors
(0.78 and 0.25) fail on both arms once it compares the two models. Recomputing the CFD
world's climate with one input swapped at a time attributes the precipitation gap entirely to
the wind's convergence pattern. The diagnostic wind alone gives 0.786 land-biome agreement and
0.498 precipitation rel-RMS. The diagnostic air temperature alone gives 0.831 and 0.004.

## Replay A/B (pre-registered)

**Plan, recorded before the runs.** Both saves are replayed for 10 Myr in 100 × 100 ky steps,
with checkpoints every 2.5 Myr, under the saves' own `"diagnostic"` wind model. The arms are
`origin/main` `0fc4901` and this branch at the commit that adds this section. Runs go one at a
time, alone on the machine.

```bash
backend/.venv/bin/python <branch>/bin/debug/replay_surface_wind.py <local-save> --myr 10 \
  --checkpoint-myr 2.5 --out <scratch-output>/<arm> [--baseline <scratch-output>/main]
```

| Metric | Expected | Pass if |
|---|---|---|
| land weathering volume per 2.5 Myr | 0 to +5% | within ±5% at every checkpoint |
| global precipitation | −1 to −2% | within ±3% |
| land precipitation | −1 to +3% | within ±8% (one save jitters ~5% between neighbouring steps) |
| land-biome agreement vs main | ≥ 0.97 at 2.5 Myr | ≥ 0.90 at 10 Myr (state diverges over time) |
| land area fraction | ~0 | within ±0.005 |
| 10 m wind, ocean / land | −14% / −20% | reported (the intended change) |
| mean ocean current | −13% | reported (no consumer reads its magnitude) |

A smoke run on save A (3 steps, branch) gave land precipitation of 297, 288 and 281 mm in
consecutive 100 ky steps, so a single checkpoint can't resolve land-precipitation changes
below about 5%. Weathering moved less than 1% between steps.

## Replay results

All four runs ran alone, one at a time. The main arm is `0fc4901` and the branch arm
`ff373af`. Save A is seed 628977350 at 226.7 Myr and save B seed 910211954 at 220.2 Myr (SHAs
above). The compact per-checkpoint numbers are in [`replay-ab.json`](replay-ab.json).
Weathering is the land-node volume demand per 2.5 Myr interval; everything else is read at
the checkpoint. Changes are branch relative to main.

| Measure | A 2.5 | A 5 | A 7.5 | A 10 | B 2.5 | B 5 | B 7.5 | B 10 |
|---|---|---|---|---|---|---|---|---|
| land weathering volume | −5.0% | −3.8% | −2.2% | −1.8% | +4.0% | +2.4% | **+6.0%** | +3.8% |
| global precipitation | −1.6% | −1.5% | −1.5% | −1.5% | −1.8% | −1.8% | −1.6% | −1.2% |
| land precipitation | −2.1% | −1.4% | −0.5% | −1.2% | −3.8% | −3.4% | **−9.8%** | −0.6% |
| land-biome agreement | 0.949 | 0.948 | 0.937 | 0.939 | 0.936 | 0.927 | **0.673** | 0.913 |
| land area fraction (main → branch) | 0.162 → 0.162 | 0.162 → 0.162 | 0.161 → 0.161 | 0.161 → 0.161 | 0.322 → 0.322 | 0.323 → 0.323 | **0.352 → 0.322** | 0.321 → 0.321 |
| 10 m wind, land | −19.7% | −19.7% | −20.1% | −19.8% | −17.1% | −16.8% | −25.0% | −17.3% |
| 10 m wind, ocean | −13.9% | −13.9% | −13.9% | −13.9% | −14.4% | −14.4% | −12.9% | −14.3% |
| mean ocean current | −13.2% | −12.2% | −12.4% | −13.7% | −14.7% | −15.0% | −15.0% | −14.8% |

**Against the pre-registered thresholds**:

- Global precipitation, land area and the 10 Myr biome agreement pass on both saves.
- Weathering fails at one checkpoint (B at 7.5 Myr, +6.0% against ±5%), and A's 2.5 Myr value
  sits on the limit at −5.0%.
- Land precipitation fails at that same B checkpoint.
- Biome agreement at 2.5 Myr is 0.94–0.95, below the 0.97 I expected but not a pass criterion.

**Save B at 7.5 Myr.** The main arm's land fraction jumps from 0.323 to 0.352 and its mean land
elevation falls from about 2000 m to 1485 m, then both return by 10 Myr; the branch shows no
such excursion. The 0.67 biome agreement and the land-precipitation, land-wind and land-area
outliers all coincide with it. The extra land alone would cost about 0.08 of agreement (9% of
main's land cells are ocean in the branch), so the rest of that drop is climate on the
newly-low land. It is a transient in one arm, not a trend, though a single replicate per arm
can't say what triggers it. The weathering miss belongs to the whole 5–7.5 Myr interval, but at +6.0% it is in line
with B's other checkpoints.

**Weathering.** No single coefficient centres both saves. Land wind falls 19.7–20.1% on A
and about 17% on B (−25% at B's 7.5 Myr transient), because the envelope cuts high-latitude
and westerly speeds more than tropical ones, and the two worlds put their land at different
latitudes. At 3.6 the mean change is −3.2% on A and +4.0% on B, which brackets zero. A
coefficient tuned to either save alone would move the other by 7–8%.

**Intended changes.** The 10 m wind drops 14% over the ocean and 17–20% over land, toward
ERA5. Mean ocean current falls 12–15% with it; nothing reads current magnitude (see
*Downstream calibration*).

**Runtime and memory.** Branch runs took 734 s (A) and 1136 s (B), main runs 1006 s and
1226 s, with peak RSS of 1.8–2.0 GB (A) and 6.8–7.4 GB (B). The wind change adds two array
rolls per climate call. The runtime gap is larger than that and was not designed as a timing
test, so read it as noise.
