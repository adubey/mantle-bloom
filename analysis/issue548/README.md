# Issue 548: RFC 0015 seed sets, `DDF_ice` and gate G5 (sea ice)

Gate work for [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md) (tracking #544)
that needs neither a D1 prototype nor RFC 0014 D2(b): the seed sets, sourcing `DDF_ice` (G2
part 1), and diagnosis gate G5.

**Result: G5 passes at the starting knobs.** Without ice ages, D6's perennial sea-ice edge
falls within 10° of latitude of the observed edge in both hemispheres on all three diagnosis
saves, and the current anomaly moves the local edge by at least 2° in both hemispheres of each
save. The proxies stay. There are three findings to settle before D6 lands:
- the north sits 8–9° equatorward of the observed edge, near the limit;
- seasonal ice reaches about 50°, 10° further equatorward than observed, so the Polar Sea Ice
  class covers 1.5–3 times the observed area;
- at the full ice-age cooling, perennial ice covers 63–84% of the ocean poleward of 50°.

All runs are on `402f009` (no model changes from `origin/main`). They read each save
without stepping it, one at a time, with a peak RSS of 0.53–0.59 GB.

## Seed sets

Every candidate local save was surveyed with
[`survey_frozen_water_saves.py`](../../bin/debug/survey_frozen_water_saves.py). The survey
reads metadata, ice volume and polar ocean only, with no gate measure. All candidates are
at a tilt of 23.5° and a solar multiplier of 1, with ice ages off and the diagnostic wind.
Excluded:
- the development save (seed 308161711);
- the 19.2 Gyr ocean world (seed 456862899: no ice, 98% ocean);
- a density-1 save;
- a second save of seed 997271774.

| Set | Seed | Age | Density | SHA-256 | Ice volume | Visible ice | Ocean poleward of 60° N / S | Replayed with |
|---|---|---:|---:|---|---:|---:|---:|---|
| Diagnosis | 997271774 | 87.4 Myr | 4 | `08ad9461…ed762b` | 5.65 M km³ | 9.31 M km² | 8.5% / 5.6% | ice ages off |
| Diagnosis | 649658380 | 6.4 Myr | 4 | `ebbfb5da…73f7b3` | 3.56 M km³ | 4.24 M km² | 4.6% / 9.7% | ice ages off |
| Diagnosis | 374464662 | 132.6 Myr | 4 | `0f8c50c9…a4279c` | 1.86 M km³ | 4.04 M km² | 9.7% / 9.4% | ice ages, 8 Myr period (report-only arm: 1 Myr) |
| Held out | 447276134 | 116.4 Myr | 4 | `1ccae715…71c513` | 5.62 M km³ | 7.63 M km² | 6.7% / 8.1% | ice ages, 8 Myr period |
| Held out | 628977350 | 226.7 Myr | 4 | `707d35d3…e4acb7` | 18.82 M km³ | 15.40 M km² | 6.7% / 7.8% | ice ages off |
| Held out | 295569319 | 180.1 Myr | 4 | `1dd920a6…089840` | 7.75 M km³ | 8.46 M km² | 6.1% / 6.9% | ice ages off |

Seed 374464662 is the ice-age diagnosis save because it has the most ocean poleward of 60° in
both hemispheres, which G5 needs, and little ice, so the forcing will show. The held-out saves
are used for no gate or knob choice. The saves are local and not checked in.

## `DDF_ice` (G2 part 1)

These values were checked against the papers themselves.
- [Hock (2003)](https://doi.org/10.1016/S0022-1694(03)00257-9), Table 1: fits on bare ice over
  whole seasons or years run 5.5–8.2 mm w.e. d⁻¹ °C⁻¹. Norway, 1961–1990 best fits: 5.5–6.4.
  Satujökull: 7.7. Qamanârssûp sermia: 7.3 as a best fit and 8.2 over 512 days.
  Nordbogletscher: 7.5 over 415 days. Single-season values reach 20.0 (GIMEX, 1028 m).
- [Braithwaite (1995)](https://doi.org/10.3189/S0022143000017846): observed factors of 7.53 and
  8.19 at Nordbogletscher and Qamanârssûp sermia (simulated 7.67 and 8.12), "around
  8 mm d⁻¹ °C⁻¹". Snow factors are 36–44% of ice.

**Start 8.0, range 5.5–8.2**, with 20.0 run as an extreme. #368's unverified 7 lies inside
the range.

## G5: sea ice

### Observed reference

[`era5_sea_ice_climatology.py`](../../bin/debug/era5_sea_ice_climatology.py) reads ERA5's
1990–2019 sea-ice concentration climatology (WeatherBench2, 1.5°) and its land-sea mask. A
cell is ice-covered on a day at ≥ 15% concentration (NSIDC's extent convention). It is
*perennial* when ice-covered every day of the climatological year, and *seasonal* when
ice-covered on some days. The regrid leaves coastal ocean cells blank: 1,812 cells, 9 M km²
north of 40° N. So the reference is computed twice, with those cells left out and with them
filled from their neighbours. The two bracket NSIDC's
[Sea Ice Index](https://noaadata.apps.nsidc.org/NOAA/G02135/) monthly extents, averaged over
the same 1990–2019:

| M km² | NSIDC 1990–2019 | ERA5, coast left out | ERA5, coast filled |
|---|---:|---:|---:|
| Arctic minimum (September) | 5.66 | 5.84 | 7.00 |
| Arctic maximum (March) | 15.07 | 10.81 | 15.16 |
| Antarctic minimum (February) | 3.09 | 3.21 | 4.06 |
| Antarctic maximum (September) | 18.59 | 18.23 | 20.00 |

| Observed | North | South |
|---|---|---|
| Perennial equivalent-area edge | 75.3–76.7° | 69.5–70.5° |
| All-ice equivalent-area edge | 62.1–70.7° | 60.0–60.9° |
| Ice-covered at least half the year | 8.9–12.2 M km² | 12.8–14.5 M km² |

The *equivalent-area edge*
([`sea_ice_extent_stats.py`](../../bin/debug/sea_ice_extent_stats.py)) is the latitude
poleward of which the ocean area equals the hemisphere's ice area. The model side uses the
same code. A gap is measured to the farther end of the observed range.

### Model

[`gate_sea_ice.py`](../../bin/debug/gate_sea_ice.py) applies D6's balance on each save's own
climate grid (2°):
- `T_si` is the land construct at sea level plus `ocean_temperature_c` minus the zonal
  baseline, with the maritime amplitude;
- growth is Stefan's law; loss is `DDF_ice` × TDD plus `F_w`.

It runs each save as saved with ice ages off, and again at the glacial maximum. Ice-age
cooling is uniform in both temperature baselines, so the saved state at full intensity (8 °C)
is the ice-age arm, with no replay. The settings are `F_w` ∈ {1.3, 2, 2.7} W m⁻² and `DDF_ice`
∈ {5.5, 8, 20}. The current anomaly poleward of 50° spans −2.7 to +2.7 °C on all three saves.

**Without ice ages**, at the central setting (`F_w` = 2, `DDF_ice` = 8):

| Save | Hemisphere | Perennial / seasonal (M km²) | Perennial edge (no anomaly) | Gap to observed | Columns where the anomaly moves the edge ≥ 2° | Polar Sea Ice class (f_ice ≥ 0.5), M km² |
|---|---|---:|---:|---:|---:|---:|
| 997271774 | N | 17.9 / 29.4 | 68.2° (68.0°) | 8.5° | 121 of 180 | 30.4 |
| 997271774 | S | 8.1 / 30.7 | 68.5° (68.0°) | 2.0° | 62 of 103 | 19.6 |
| 649658380 | N | 6.6 / 23.3 | 68.3° (68.0°) | 8.4° | 15 of 114 | 15.5 |
| 649658380 | S | 17.9 / 40.4 | 68.4° (68.0°) | 2.1° | 84 of 180 | 34.2 |
| 374464662 | N | 18.1 / 32.8 | 68.3° (68.0°) | 8.4° | 108 of 180 | 32.6 |
| 374464662 | S | 17.0 / 36.8 | 68.6° (68.0°) | 1.9° | 84 of 180 | 31.1 |

**Sensitivity**, without ice ages: the perennial edge's gap to the observed edge, north / south.

| Setting | 997271774 | 649658380 | 374464662 |
|---|---:|---:|---:|
| `DDF_ice` 5.5 | 9.9° / 3.9° | **10.4°** / 3.8° | 9.7° / 3.7° |
| `DDF_ice` 8 (central) | 8.5° / 2.0° | 8.4° / 2.1° | 8.4° / 1.9° |
| `DDF_ice` 20 | 6.4° / 1.0° | 6.5° / 0.9° | 6.4° / 1.0° |
| `F_w` 1.3, `DDF_ice` 8 | 9.3° / 2.8° | 8.7° / 2.7° | 8.9° / 2.7° |
| `F_w` 2.7, `DDF_ice` 8 | 8.5° / 2.0° | 8.4° / 2.1° | 8.4° / 1.9° |

- **`F_w` hardly matters.** Bottom melt of 0.13–0.28 m/yr is small against 2–3 m of
  Stefan growth poleward of 66°.
- **`DDF_ice` moves the edge about 4° over its range.** Only the low end, 5.5, takes seed
  649658380's north past 10°, by 0.4° and only against the far end of the observed range.

**At the glacial maximum** (8 °C cooling, central setting):
- the perennial edge moves to 55.3–55.5° in both hemispheres, 15–21° equatorward of today's
  observed edge;
- perennial ice covers 73–81% of the ocean poleward of 50°;
- today's model draws 15–23 M km² of Polar Sea Ice in this state (the SST cutoff); D6's class
  would draw 30–59 M km² per hemisphere.

### Zonal bands by tilt

No anomaly, no ice ages, pure functions on a 0.25° grid
([`zonal_bands.json`](zonal_bands.json)): the latitude from which ice is perennial, and from
which it freezes at all, at `F_w` = 2.

| Tilt | Perennial, `DDF_ice` 5.5 / 8 / 20 | Freezing |
|---:|---|---:|
| 15° | 63.9° / 64.4° / 65.6° | 54.9° |
| 23.5° | 66.9° / 67.9° / 69.9° | 51.1° |
| 35° | 77.6° / 81.9° / none | 45.9° |

A larger tilt warms the polar summer and cools the winter, so perennial ice retreats
poleward while freezing reaches further equatorward. At 35° with `DDF_ice` = 20, no latitude
keeps perennial ice.

**Correction to RFC 0015 D6's scale check.** It says the bands are the same for `DDF_ice` = 4
and 7 mm. On a 0.25° grid they aren't: the perennial edge is 65.9° at 4 and 67.4° at 7 (both
at `F_w` = 2). The original check's 5° latitude rows (`measure_frozen_water.py`) couldn't
separate them. The 65.75° it quotes is the `DDF_ice` = 4 value.

### Outcome

- **Keep the proxies.** In each hemisphere of each diagnosis save, the equivalent-area edge
  without ice ages lies within 10° of the observed edge at the central setting (at most 8.5°),
  and the anomaly moves the local edge by ≥ 2° in each hemisphere of each save.
- **Not falsified.** Perennial ice covers ocean on every save, and without ice ages it covers
  16–41% of the ocean poleward of 50° across all settings, never most of it.

**For the owner, before D6's implementation issue:**
1. **The north is biased 6–10° equatorward** (66–70° against 75–77°) across the settings,
   leaving about 1.5° of margin at the central one. The south is within 1–4°. The knobs move
   the north by at most 4°, so the bias belongs to the proxy, not to a knob setting. It is
   consistent with the land construct at sea level being too cold in the Arctic summer for
   the sea surface.
2. **Seasonal ice reaches about 50.5° in both hemispheres** (observed all-ice edge: 62–71°N,
   60–61°S). D6's class (f_ice ≥ 0.5) then covers 15–34 M km² per hemisphere against the
   observed 9–15 M km² of ocean ice-covered at least half the year. D7 draws seasonal sea
   ice only through that class. Two options: raise the class threshold (G5 set it as a
   design choice), or give the freezing side the ocean's heat (an FDD proxy from the moderated
   air over open water). Both are revisions for `rfc-authoring`, not a G5 failure.
3. **At full ice-age cooling, perennial ice covers 63–84% of the ocean poleward of 50°.**
   The falsification criterion's "most of the ocean poleward of 50°" is read here as
   applying without ice ages, like its first clause. The glacial-maximum response is set by
   `ICE_AGE_MAX_COOLING_C` = 8 °C, which applies uniformly. The RFC should state which
   reading holds.

## Reproduce

```bash
cd backend
.venv/bin/python ../bin/debug/survey_frozen_water_saves.py <saves...> --out <scratch>/survey.jsonl
.venv/bin/python ../bin/debug/gate_sea_ice.py <save> --out ../analysis/issue548/g5-seed<seed>-<age>y.json
.venv/bin/python ../bin/debug/gate_sea_ice.py --zonal-only --out ../analysis/issue548/zonal_bands.json
cd ..
# Any Python with numpy and numcodecs:
python bin/debug/era5_sea_ice_climatology.py --out analysis/issue548/era5_sea_ice_climatology.json
```

Files: `g5-seed*.json` (per save: both arms, every setting; latitude bands at the central
setting only), `era5_sea_ice_climatology.json`, and `zonal_bands.json`.
