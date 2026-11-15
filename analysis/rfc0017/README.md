# RFC 0019: atmospheric moisture diagnosis

## G0: archived G1b evidence (partial)

The owner accepted [RFC 0019](../../rfc/0019-conservative-atmospheric-moisture.md) on
2026-11-14; [#898](https://github.com/adubey/mantle-bloom/issues/898) tracks the work.
This is an audit of existing outputs, not a new simulation or a completed G0.

The `rfc0016-g1b-outputs.tar.gz` archive from #771 has SHA-256
`90277c33c631239b480fe55060cb1a27deb531ea3f851d8467002131077828de`.
All 84 JSON files pass the archive's `SHA256SUMS` checks. Their recorded source commit is
`0bda513df0c9b0baeb3cd708be142f82eb332f1d`; the original run commands, checkpoint
hashes, controls and limitations are documented in [RFC 0016 G1b](../rfc0016/README.md#g1b-result).
The raw outputs remain private.

For each window, calculate each step's central-treatment / control ratio of
`land_rain_only.rain_m3`, then average the ratios after the first step, when the
humidity feedback has not yet acted. This reproduces #771's rounded values:

| Window | Step size | Mean ratio after first step | All steps, for context |
|---|---:|---:|---:|
| A | 100 kyr | 1.3807 | 1.2855 |
| D | 10 kyr | 1.4802 | 1.4682 |
| E | 10 kyr | 1.5078 | 1.4951 |

This confirms the arithmetic behind the reported 38–51% rise conditional on excluding the
initial step. It does not independently reproduce the simulation, establish a global water
budget, or establish that the rise persists on current source. The averaging convention
matters: including the initial, uncoupled step gives the lower ratios shown above.

Still required for G0: check exact input bytes and current-source compatibility, reproduce
`measure_climate_fields.py` for the named worlds, independently reproduce G1b where its
pinned code and inputs permit, trace the complete source/sink ownership, and establish
area-weighted totals. Do not mark the gate passed on this archive audit alone.

## Held-out registration before G0 replays

Registered on 2026-11-14 from the private bucket inventory, before any RFC 0019
model replay or tuning. These saves are reserved for G3 acceptance, not G0–G2
calibration or design selection. Selection used the save name, seed and hash only;
no climate outputs from either save were inspected for RFC 0019. SHA-1 identifies
the exact `.bloom` object in the bucket; verify its bytes before G3.

| Save | Bucket object SHA-1 |
|---|---|
| `mantle-bloom-seed507012855-479000000y.bloom` | `a249a559fac44e3f2b9d2dc76a1a26ca2711b659` |
| `mantle-bloom-seed762517225-135200000y.bloom` | `fe140dd3e8b458fe99a4aaefefc894b6647d3200` |

## G0: current-source writer and ownership inventory

The diagnostic backend is based on `9f8267c49b2da23c9be544f6de6606bf1f0e7135`
plus this RFC work; none of its changes alters the simulation. On a climate update,
`erosion.apply_erosion` computes `climate.compute_climate`, stores the result in
`world.climate_cache`, samples precipitation and humidity at terrain nodes, and passes
precipitation to `hydrology.compute_hydrology`. Climate uses the previous cache for its
vegetation input. The biome classifier reads precipitation and temperature; erosion and
geology also read them. Ice receives the partitioned precipitation and climate forcing.

| Process | Present writer or consumer | Ownership / unit limitation |
|---|---|---|
| Ocean moisture | `open_water_humidity_g_kg`, zonal and meridional sweeps | Resets humidity to a temperature-derived g/kg target; there is no ocean evaporation mass or ocean-water debit. |
| Land humidity sources | `_land_moisture_source`, `_vegetation_transpiration_source` | Lake/inland-sea and river depths supply g/kg proxies; frozen water contributes zero. Vegetation is prior precipitation divided by the fitted rain conversion, times 0.16 and a biome factor. No surface water account pays these inputs. |
| Atmospheric transport and rain | `_humidity_zonal_sweep`, `_humidity_meridional_sweep`, `compute_humidity`, `compute_precipitation` | Retention loss, sweeps' caps, noise and the humidity floor are unledgered; orographic loss is recorded separately. Rain combines humidity, orographic and convergence terms and a floor, without debiting humidity. The atmospheric field is not an extensive water mass. |
| Land rain and runoff | `glacier_rain_partition_mm`, `_partition_precipitation_sources`, `breaching.runoff_mm`, river routing | Wet/dry and snow/rain portions are partitioned; dry rain is split between Budyko runoff and aggregate ET, and routed water enters lakes or ocean. Aggregate dry ET is not a separately measured vegetation transpiration source. |
| River and lake losses | `river_transit_loss_fraction`, `river_evaporation_fraction`, lake phase balance | River losses and lake evaporation/sublimation are accounted in physical hydrology ledgers, but are not transferred as equal atmospheric inputs. Climate's lake/river proxies are separate functions of depth. |
| Snow and ice | `ice_sheet.SurfaceClimate.rates`, `glacier_rain_partition_mm`, `ice_sheet.step_ice` | Frozen precipitation feeds the glacier account; ice and lake losses and melt follow their respective ledgers. These do not close an atmospheric budget. |
| RFC 0016 treatment | `bin/debug/rfc0016_land_store.py` | The offline treatment debits a land store and injects its last-step transpiration divided by `PRECIP_MM_PER_G_KG` into climate. It is a diagnostic treatment, not the current production moisture source. |

The quantities therefore cannot all be expressed as physical atmosphere source/sink
volumes without introducing an unvalidated column or mixing-depth conversion. In
particular, the ocean reset is not evaporation, and a humidity-proxy integral is not a
water mass. The climate-grid inventory below uses exact spherical latitude-band areas
and explicitly keeps the two kinds of quantity separate. Inputs are the converted B2
objects; their current SHA-256 values are recorded in the output. Ocean and land rain
are diagnostic annual rates; they are not independent ledgered atmospheric sinks.

| Window | Land rain km³/yr | Ocean rain km³/yr | Ocean target | Lake/sea proxy | River proxy | Vegetation proxy |
|---|---:|---:|---:|---:|---:|---:|
| A | 82,680 | 495,142 | 4,541.54 | 49.39 | 4.83 | 64.52 |
| D | 80,593 | 489,436 | 4,512.76 | 32.16 | 3.97 | 63.31 |
| E | 79,173 | 490,865 | 4,521.43 | 34.34 | 4.30 | 60.21 |
| Source 220.2 Myr | 81,100 | 502,956 | 4,601.76 | 46.58 | 3.51 | 59.92 |
| Cold source 200.8 Myr | 116,975 | 573,761 | 5,183.55 | 24.48 | 1.83 | 91.48 |

The last four columns are area integrals in `10^12 m²·g/kg`, not fluxes. Ocean is
a target rather than an additive input; land values are pre-cap per-sweep proxies,
not effective source contributions after transport. They cannot be summed as a budget.
Reproduce with `backend/.venv/bin/python bin/debug/inventory_rfc0019_moisture.py
<save> --out <private-output.json>`; the script reports exact units and save hash.

As a separate existing physical land-water accounting example, the archived G1b control
runs' mean ledger terms (km³/yr, each per-step volume divided by that run's step size)
are:

| Window / step | Land rain | Dry ET | River evaporation | Lake evaporation | Lake sublimation | Ocean export |
|---|---:|---:|---:|---:|---:|---:|
| A / 100 kyr | 62,267 | 30,258 | 2,612 | 5,652 | 175 | 23,570 |
| D / 10 kyr | 67,299 | 32,348 | 3,154 | 4,465 | 2,101 | 25,229 |
| E / 10 kyr | 65,035 | 31,972 | 2,988 | 4,362 | 1,326 | 24,385 |

These are time averages of stepped, terrain-node land ledgers on the pinned G1b source,
not the climate-grid instantaneous snapshots above; they are not additive to those rows.
Storage changes, meltwater, and other ledger terms are omitted from this abbreviated
example, so its columns alone do not close the ledger.

### Input bytes and current climate snapshots

The five exact `.mbworld` byte sequences named by #771 were not returned by B2: lookup
resolved each legacy name to a converted `.bloom` object. Bucket metadata retains the
legacy SHA-256 values named by #771, but that metadata alone does not establish identity
of the loaded simulation state. These are the SHA-256 values of the downloaded bytes:

| Object | SHA-256 |
|---|---|
| `seed910211954-221.6myr.bloom` | `1f167a8d9e0edf06c3c52b23598dc2f776adc2927f843beb992b5193ad30fd67` |
| `seed910211954-225.1myr.bloom` | `c5b1f5d832e83e63f2a90a9f9f7bff03c885e831a0c6e585e22678a718168043` |
| `seed910211954-224myr.bloom` | `3e0fa48182b86b7d7e293cf73c180c6a4667f1b4c48dac7effaf57c9d2b87ca7` |
| `mantle-bloom-seed910211954-220200000y.bloom` | `8a51374c6b95335db7642be356b064d8237b2262909c3b4771fffa5cb8de43ed` |
| `mantle-bloom-seed317739385-200800000y.bloom` | `c7e8c4094376b2b57158e36b1066cf529eb475ac0bbf54a6b888b391f294d344` |

`measure_climate_fields.py <save> --myr 0 --out <private-output.json>` on current
source gave the following instantaneous snapshots (land and ocean mean precipitation in
mm/yr, land share of global precipitation in percent, land-mean humidity sources in
g/kg). `--myr 0` deliberately measures the loaded checkpoint; it does not step the world.

| World | Land rain | Ocean rain | Land rain share | Vegetation proxy | Lake/river proxy |
|---|---:|---:|---:|---:|---:|
| A | 468.393 | 1,484.478 | 14.31 | 0.366 | 0.307 |
| D | 451.589 | 1,475.991 | 14.14 | 0.355 | 0.202 |
| E | 445.221 | 1,477.459 | 13.89 | 0.339 | 0.217 |
| Cold source | 1,061.365 | 1,434.930 | 16.93 | 0.830 | 0.239 |
| Source 220.2 Myr, summary skipped | 493.922 | 1,454.181 | 13.89 | 0.365 | 0.305 |

On the 220.2 Myr source, the same script fails in its pre-climate ice/land summary:
the saved hydrology `is_ocean` array has 133,460 entries while the current plate area
array has 133,489. The separate climate-grid inventory above does produce a snapshot,
and retains the saved cache's own sampling semantics. A temporary copy of
`measure_climate_fields.py` with only its cached hydrology summary branch disabled
then produced the source climate snapshot shown above, including its latitude bands
and Köppen groups. This is a partial workaround: the original script still fails and
no ice or land percentage is reported for that snapshot. Reproduce by copying the
script and `replay_lib.py` to a scratch directory, changing `if hydro is not None:` to
`if False and hydro is not None:` in the copy, and running the copy with `--backend`
pointing to this checkout's backend, `--myr 0`, and the source save. Neither the
world nor its hydrology cache is modified by that skip.

### Independent G1b replay: window A

Using the diagnostic files from PR #764's head `ad3575850f02bcc2c593d5e80f7a7f94e71b0bc5`
and the pinned backend `0bda513df0c9b0baeb3cd708be142f82eb332f1d`, rerun A's
control and central treatment with the converted A save, `--myr 0.4 --step-years
100000 --checkpoint-myr 0.1`, adding `--treatment` for the central arm. Both runs'
`land_rain_only.rain_m3` match the corresponding archived values exactly at all four
steps. The control's dry ET, river and lake evaporation, and ocean export also match
exactly at all four steps. The mean of treatment/control rain ratios after the first
step is 1.380703848756227. This supports simulation parity for those measured outputs,
not byte identity or parity of every field.

Window D was independently replayed in both arms at the registered `--myr 0.4
--step-years 10000 --checkpoint-myr 0.1`. All 40 control and 40 central rain
values match the archived outputs exactly. Control dry ET, river and lake evaporation,
and ocean export also match at every step. The mean central/control rain ratio after
the first step is 1.480160784723084.

Window E was independently replayed with the same 40-step, 10 kyr schedule. Rain
matches the archived control and central runs exactly at every step; the control dry
ET, river and lake evaporation, and ocean export match exactly as well. The mean
central/control rain ratio after the first step is 1.5077775929811301. These three
paired replays verify the reported range on the converted inputs and pinned backend;
they do not establish full state parity or the physical atmospheric balance.

For all six independently reproduced A/D/E control and central runs, every step's
`land_rain_only` rain, dry ET, river evaporation, lake evaporation, and ocean export
matches the corresponding archived value exactly. This checks 84 reproduced steps
across the two schedules, not all 84 original G1b arms.

## G0 assessment

The 38–51% premise reproduces for A, D and E on the pinned backend and converted
checkpoints; the writer inventory and available source/sink and proxy totals are above.
The exact legacy input bytes were not available through B2, the original measurement
script fails on one source cache, and physical atmospheric evaporation/transpiration
fluxes and closure totals are not defined in the current model. The proxy totals do not replace the missing physical quantities. On 2026-11-14 the owner
accepted these documented limitations as satisfying G0 and authorized G1. This decision
does not establish atmospheric closure.

**Owner decision (accepted 2026-11-14):** Treat G0's input reproduction requirement as
satisfied by the recorded converted-file hashes plus the exact match of the five checked
ledger terms in both arms across A, D and E, and accept the source snapshot with its
cached land/ice summary omitted. Treat the area-weighted physical atmosphere source
and closure totals as unavailable in the current implementation, with the separately
reported proxy integrals and physical hydrology losses completing the inventory. This
would permit G1 and G2 to proceed without claiming atmospheric closure; the unmodified
script failure and missing original bytes remain open limitations. This interpretation is accepted for G0; the limitations remain.

## G1: literature and units (in progress)

Crossref resolves all three RFC references with the expected authors, titles and years.
The following checks use the sources themselves; they constrain scales, not synthetic
world calibration targets.

- [Trenberth (1998), Table I](https://staff.cgd.ucar.edu/trenbert/trenberth.papers/ClilmChg39_667-694.pdf),
  reports global-mean precipitable water 24.6 mm and precipitation 2.69 mm/day, so
  their ratio is 9.14 days. The same paper calculates 8.1 days from the inverse of
  the area-mean local depletion rate, and distinguishes it from the ratio of global
  means. Neither number implies a universal local rainout time. As liquid water
  equivalent, 1 mm of column water is 1 kg/m²; 1 mm/day is 1/86,400 kg/m²/s.
- [van der Ent et al. (2010)](https://doi.org/10.1029/2010WR009127), Table 1,
  reports a 40% continental-origin share of terrestrial precipitation for 1999–2008,
  and 57% of terrestrial evaporation returning to continental precipitation. Those
  are distinct source fractions, not local plant transpiration fractions.
- Oki and Kanae's 2006 water-cycle figure gives 436.5 and 65.5 thousand km³/yr for
  ocean evaporation and terrestrial evapotranspiration, respectively, and 391 and
  111 thousand km³/yr for ocean and land precipitation. These imply 86.95% and
  77.89% ocean shares, respectively; Antarctica is omitted in the land figure.
  See the figure reproduced by its author in [Oki (2016), slide 2](https://gewexevents.org/wp-content/uploads/Day3_Oki_IWMinLSMs.pdf).
- [ECMWF's ERA5 documentation](https://confluence.ecmwf.int/pages/viewpage.action?pageId=611555221)
  says to derive 2 m specific humidity from 2 m dew point and surface pressure,
  using saturation over water even for subfreezing dew points. Column vapour is
  kg/m²; vertically integrated q·wind·dp/g is kg/(m·s). Precipitation in the
  selected six-hour accumulation is metres of water equivalent per six hours, not
  an instantaneous rate; multiplying by 4,000 gives mm/day.

The WeatherBench2 1990–2019 climatology contains near-surface and column variables,
but not a native vertically integrated vapour flux. The planned same-time 2000 ERA5
sample instead integrates q·u and q·v over its 13 available pressure levels, excludes
below-ground levels and approximates the unobserved top and near-surface layers. It will
also integrate q alone and compare the result against ERA5's total-column vapour to
show the effect of that vertical approximation. A 12-month sample of 48-hour windows
is not a 30-year climatology and cannot by itself establish a precise universal
transport multiplier or mixing depth. The raw pressure-level probe is queued until
compute capacity is available.
