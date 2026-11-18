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

### G1: ERA5 sample and calibration screen

The WeatherBench2 1990–2019 climatology has near-surface and column variables but no
native vertically integrated vapour flux. The completed probe instead samples twelve
48-hour windows around the start of each month in 2000 at 1.5° resolution (24 days,
eight six-hourly samples per window). The January window begins on 1999-12-31. The
main fields are WeatherBench2 conservative regrids; its sampled 1.5° dewpoint is all
non-finite, so specific humidity is derived from valid native 0.25° hourly dewpoint
sampled at the exact 1.5° grid centers. This point sample does not match the area
averaging of the other fields.

Reproduce with `python bin/debug/era5_rfc0019_g1.py --out
analysis/rfc0017/g1-era5.json` in an environment with numpy and numcodecs. The script
records source URLs, the sampling schedule, method and SHA-256 of each compressed
source chunk. `g1-era5.json` is the compact checked-in output; it contains no saves or
full downloaded arrays. The integrated flux is an approximation from 13 pressure
levels (50–1000 hPa), omitting below-ground levels, holding the top value to 0 hPa,
and interpolating the surface layer from 2 m humidity and 10 m winds. A synthetic
check covers native-grid alignment and constant-humidity/wind integration. The selected
sample footprint is recorded per class because polar land has relatively few cells
passing the proxy-wind threshold.

| Class | Column water (kg/m²) | Pressure-integrated water (kg/m²) | Precipitation (mm/day) | Ratio of means (days) | Effective depth p10–p90 (km) | Signed IVT/proxy fit |
|---|---:|---:|---:|---:|---:|---:|
| Global | 23.69 | 23.98 | 2.93 | 8.10 | 1.69–2.45 | 1.06 |
| Land | 18.32 | 18.49 | 2.21 | 8.28 | 1.90–2.63 | 1.36 |
| Ocean | 25.86 | 26.19 | 3.21 | 8.05 | 1.65–2.38 | 1.05 |

The pressure-integrated water agrees with the supplied total-column water to within
1.3% in the global, land and ocean means, a useful check on the vertical approximation
for this sample. The signed fit regresses approximate integrated vapour transport on column water
times the 10 m wind vector. Its latitude-band values range 0.72–2.83; land and ocean
do not support one precise common mapping. Cell-level parallel ratios are broader and
can be negative where vectors oppose, so they are not interpreted as a positive local
constant. The 1.5–3.0 km effective-depth range and 0.5/1/2 wind-multiplier sensitivity
are provisional calibration ranges, not measured physical constants. Effective depth
is a column-to-near-surface conversion, not boundary-layer depth. The 24 days are not
a climatology and cannot establish interannual or seasonal uncertainty.

As a comparison, Trenberth's global means give 9.14 days; the ERA5 sample's 8.10-day
ratio is a same-grid sample ratio, not an annual estimate or a local rainout time. Its
column water and precipitation are within a plausible global scale. This is a unit and
scale sanity check only. Full per-band metrics and compressed-chunk hashes are in
[`g1-era5.json`](g1-era5.json). The earlier failed runs remain private diagnostic
artifacts; empty-bin and non-finite guards are retained in the script.

### G2: preliminary analytic prototype screen

`backend/.venv/bin/python bin/debug/rfc0019_transport_prototypes.py` compares two
diagnostic formulations on an equal-area lat/lon sphere: an upwind finite-volume
face-flux generator, and local wind-arrow transitions. Both solve the same absolute
stationary mass balance with a positive rain sink using a sparse direct solve. It
asserts nonnegative mass and per-cell and global closure below `1e-8` relative to
source. This is an isolated script; it does not change `World` or production climate.

On the 8×16 ocean-source case, both options close globally to below `1.3e-16`, with
maximum cell residual below `1.6e-17`; minimum column water is positive. Ocean-only
source tags produce nonzero rain over land, and ocean/land source superposition error
is below `2.3e-16`. An enhanced local sink increases rain at its target cell while
retaining global closure below `5.7e-16`. A separate 18-day single-cell availability
check caps requested rain of 1.98 kg to the 1.0 kg opening inventory. The script also
checks zero sources, calm cells, convergence/divergence, a closed wind loop with rain,
longitude wrap, equal-area pole boundaries, sinkless closed-class detection, and a
uniform-source resolution sweep. Uniform-source column water is 23.5872 kg/m² at
4×8, 8×16, and 16×32 for both schemes.

An initial combined 90×180 (climate density 1) run confirmed positive water and closure
but overlapped other simulations, so its timings were discarded. Separate sequential
runs after the budget showed no competing simulations measured finite-volume at 0.02865 s
and wind arrows at 0.02824 s for the sparse solve section. Both methods have effectively
the same cost at this scale; the 1.4% difference is below a useful tuning signal. Finite
volume closed to `1.25e-16` globally and `4.69e-18` per cell with a 22.48 kg/m² minimum;
wind arrows closed to zero globally and `5.95e-18` per cell with a 22.48 kg/m² minimum.
Peak RSS was not captured. The production default at density 4 is 360×720; direct-solver
memory and runtime there remain unmeasured. The cap check is a one-cell availability
identity, not a time-integrated transport update. The analytic G2 gate is complete, with
the one-cell cap scope and density-4 scaling limits recorded for owner review before G3.
Full outputs are in [`g2-prototype.json`](g2-prototype.json),
[`g2-prototype-90x180.json`](g2-prototype-90x180.json), and
[`g2-prototype-90x180-uncontended.json`](g2-prototype-90x180-uncontended.json).

## G3 preregistration and density-4 solver preflight

After reviewing G2, the owner approved moving to G3 on 2026-11-15 and approved the numeric
limits before any climate comparison. See the [G3 preregistration](g3-preregistration.md)
for the registered tuning/held-out split, limits, source conversion hypothesis, and metric
definitions. The tuning-only set is windows A, D, E and the cold 200.8 Myr source. The two
held-out saves remain acceptance-only and are not to be loaded before the tuning parameters
are frozen.

Before applying the density-4 solve to saved worlds, the synthetic finite-volume solver was
run once at 360×720 using the G2 prototype and its uniform source/circulation case. It closed
globally to `2.12e-15` relative and locally to `1.71e-17`; minimum column water was
`0.1113 kg/m²`; solve time was `0.6924 s`. This is only a density-4 linear-solver preflight,
not a climate result. Peak RSS was not captured. It does not establish climate-update runtime,
wind-mode behavior, or validation of the source mapping.

## G3 tuning screen: initial source mapping does not pass

The preregistered initial proxy-to-source mapping was screened on the four tuning saves
(windows A, D, E and the cold 200.8 Myr source), at densities 1 and 4, in diagnostic and
CFD wind modes, with both the finite-volume and wind-arrow solvers. The measured candidate
used effective depth 2.25 km, source timescale 9.1 days, rainout 5 days, and wind multiplier
0.5. The lagged climate cache was area-resampled to density 1 where needed. No held-out save
was opened.

**Disposition: stop before held-out testing.** The implementation conserves water, and the
grid-resolution check passes, but this tuning screen indicates the candidate cannot meet the
held-out climate-pattern target and the measured density-4 runtime exceeds its limit. Across
the four tuning worlds and both wind modes, the best-fit precipitation relative RMS after
allowing a single uniform scale is 0.669–0.695 at density 1 and 0.897–0.911 at density 4
for wind arrows. This diagnostic scale fit removes total-rain amplitude as a confound; it
still exceeds the 0.50 target by a wide margin, so the unscaled field cannot pass either.
The density-4 wind-arrow climate update is 2.97–3.65× the current-source update at the
selected settings.
The full density-1/4 tuning sweep on window A found no bracket setting meeting both the
0.50 best-fit RMS and 2× density-4 runtime limits: best density-4 best-fit RMS was 0.898,
and fastest density-4 update was 2.58×. The finite-volume density-4 path was much slower
(about 50–60× in the registered batch); timing component records are retained for diagnosis.

Other preregistered checks: candidate global rain changes by 0.10% or less between densities
1 and 4, and land/ocean shares shift by at most 0.63 percentage points. Maximum relative
water closure across the candidate batch is below `3e-15`; water and rain are nonnegative
after only roundoff-scale clipping. At the chosen settings, candidate global rain is 11–18%
below current-source rain on the tuning worlds, which is not a requirement violation by
itself but is consistent with the G0 direction. Climate-pattern error, rather than closure
or grid convergence, is the binding scientific failure.

The initial mapping was tested at the full preregistered 1.5/2.25/3 km depth, 5/9.1/18 day
rainout, and 0.5/1/2 wind-multiplier bracket on window A. Both diagnostic/CFD winds and both
solvers were included. Since even the best amplitude-adjusted spatial fit and fastest density-4
update miss their limits, the two acceptance-only saves remain sealed. These are tuning-screen
results, not held-out acceptance measurements. G3 is **not complete**. The follow-up source
mapping screen below did not resolve the spatial mismatch; the next design choice is whether
to revise the rainout/convection field or the transport approximation before any held-out
evaluation. Do not mark the RFC implemented or begin G4.

The full compact comparison is stored privately at
`~/mantle-bloom-captures/rfc0019/g3-registered-tuning-candidate.json` (SHA-256
`143e83b50e9c99c4d8545b2909a82a88c67c6cbc4a5be3b4e88b9f14e2624458`). The measurement
code at commit `25e4a3132fe06a3b23305aacaf01c5f08631ba91` is identified by SHA-256
`8f53951cd740905d6e9c48875bbe2a2995e0eeed810b1f69a50335da7c3b51f0` for
`bin/debug/rfc0019_g3_compare.py` and
`25853a43dabf03edc32891454f0e210252e81cc561d018f57a0f4519ed48b9cc` for
`bin/debug/rfc0019_transport_prototypes.py`. The run command was:

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_compare.py A D E cold \
  --density both --wind both --solver both \
  --effective-depth-km 2.25 --rain-days 5 --wind-multiplier 0.5 \
  --out ~/mantle-bloom-captures/rfc0019/g3-registered-tuning-candidate.json
```

### G3 source-mapping revision screen

At the owner's direction, two tuning-only source revisions were screened without opening
held-outs. First, one shared multiplier on the existing lake/river/vegetation proxies was
calibrated from the pooled source totals to bring aggregate land-origin input toward the RFC's
broad 13% Earth scale check (the complement of the cited 87% ocean contribution). A multiplier
of 5 raised the tuning-set land-source fraction from 2–4% to about 10–18%, depending on
world, wind mode, and resolution. It brought land rain closer on some cases but did not
improve density-4 best-fit spatial RMS: 0.896–0.922 versus 0.897–0.911 for the initial
mapping.

Second, the ocean target-humidity proxy was replaced by a bulk flux
`E_ocean = rho_air * C_E * |u| * max(q_sat(SST) - q_air, 0)`, using the legacy climate's
near-surface humidity as `q_air` for this offline screen. This follows the bulk-flux form
described in the source reference for the RFC. It is not a solved-candidate humidity and is
not ready as a production source law. With `C_E = 0.003` and the same land multiplier of 5,
the density-4 wind-arrow spatial RMS after a uniform amplitude fit was 0.892–0.925; closure
remained below `2.2e-16`, while runtime was 2.8–3.75× current-source time. The ocean share
of source was 81–89%, varying by wind mode and world. The coefficient bracket 0.001/0.003/0.005
and land multiplier bracket 1/5/10 did not approach the 0.50 spatial target. At `C_E=0.001`
global rain was 57–67% below the baseline; at `C_E=0.005` it was 30–46% above. Land multiplier
10 overproduced land rain and worsened the spatial fit; multiplier 1 had the lowest fitted
RMS (0.889 minimum) but left land sources at only a few percent of total.

**Conclusion:** source scaling and a bulk ocean humidity-deficit law are not sufficient to
recover the current spatial precipitation pattern on the tuning worlds. The search has not
produced a revised mapping that clears climate fidelity, and all tested density-4 updates
remain above the runtime limit. This is still tuning evidence, not a held-out result. The next
design question is whether to add a spatially varying rainout/convection law or revisit the
transport approximation. Keep the held-out saves sealed until a candidate is both plausible
and fast enough on the tuning set.

Compact outputs remain private and are checksummed in the RFC capture directory:

| Run | SHA-256 |
|---|---|
| Land-proxy multiplier 5 | `e1a914543e72400550ef910d008ae650912379fac72038b84d2d268af928d13f` |
| Bulk ocean flux, `C_E=0.003`, land multiplier 5 | `88c3a076ae8f42ffa741b6d8670f3c9bb6be211ff6904f12ef3da6cef97fb218` |
| Bulk coefficient bracket 0.001/0.003/0.005 | `25aac5b69084fd4e4af36bd1f10244fee60d3e2566a09a7c61250b51fad41c46` |
| Bulk ocean plus land-proxy multiplier bracket | `e35d97641b6183f245753548da7b9c2486ad92710d94d45576a3dcf7b53cb6c8` |

The land-proxy-only run used comparison-script SHA-256
`dc5cb14bf1244f3ae7a59f362408d197a2c1d304883417b365945f7c47cbc4ed`; the first bulk-flux
run used `f6a1f0b0311dbe62b3352f5d123a9db652c53bc0839a5aa357e5cd691f5ddaa9`; the bulk
coefficient sweep used `22be2ab75ccec8312124a7028821f18f41b8a5723b588d14d99c51ab757b8d03`;
and the bulk plus land-proxy sweep used
`7ab72dbd0383ea9dcd2b413d52ace16fb6df2ccd8d2c5eb2510f54de8743bd18`. All runs used code
commit `25e4a3132fe06a3b23305aacaf01c5f08631ba91` and the registered tuning-save hashes
listed above.

### G3 ERA5 spatial-rainout sensitivity

The G1 sample also supplies a diagnostic depletion time by dividing each latitude/surface
class's mean column water by its mean precipitation. The 24-day sample gives the following
ratios in days; these are not local rainout constants or a climatology.

| Absolute latitude | Land | Ocean |
|---|---:|---:|
| 0–15° | 8.12 | 9.08 |
| 15–30° | 11.56 | 12.65 |
| 30–45° | 8.79 | 6.22 |
| 45–60° | 5.61 | 3.84 |
| 60–75° | 4.84 | 3.46 |
| 75–90° | 3.79 | 5.88 |

I tested this piecewise 15-degree land/ocean rainout field with the bulk ocean source
(`C_E=0.003`), land-proxy multiplier 5, and the fixed G3 transport settings. Closure remained
below `4e-16`, but the amplitude-adjusted precipitation rel-RMS worsened to 0.766–0.805 at
density 1 and 0.934–0.946 at density 4. Density-4 update time was 3.08–3.93× baseline.
This coarse field derived from a short sample does not recover the tuning-set pattern; do not
carry it forward as the selected rain law.

The result narrows the next design choice: source reweighting and a first class-dependent
rainout field are insufficient. Before opening held-outs, the owner must decide whether to
revise moisture transport/convection more substantially or reopen the steady diagnostic
formulation. The private comparison is
`~/mantle-bloom-captures/rfc0019/g3-bulk-source-era5-spatial-rainout.json` (SHA-256
`19b2bbb71dbcb5098c817d2c42143499969a847818f38931899ec4b977b14360`); its measurement script
SHA-256 is `cc1cc30e594b87fcbd607077048f5e98cb4dd555a92b844a0b6ed6f2ad5a9cdc`.

### G3 steady-formulation screen: saturation-triggered rainout

After the source and coarse rainout screens failed, the owner directed the investigation to
revisit the steady formulation before opening held-out saves. I kept the bulk ocean source
(`C_E=0.003` with legacy humidity as ambient `q_air`), land proxy multiplier 5, 2.25 km
effective depth, 5-day background rainout, and 0.5 CFD-wind multiplier fixed. On the registered
tuning windows A, D, and E at density 4, I added a conservative nonlinear sink: background
rainout `W/τp` plus a one-day convective sink for column water above a specified fraction of
saturation. Saturation column water was estimated from the 2.25 km effective depth and the
local surface saturation humidity. This is a screening closure, not a calibrated cloud model.
The cold tuning save and both held-outs were not opened.

The sink closed the global budget exactly at displayed precision for every case. At an 80%
saturation threshold it activated in 3.3–3.6% of cells, with best-fit precipitation rel-RMS
0.906–0.907 and density-4 climate-update time 18.1–18.9× the current-source update. At 60%,
the active fraction rose to 7.5–8.1%, but best-fit rel-RMS was still 0.898–0.900 and runtime
was 17.6–21.2×. Both thresholds remain far outside the preregistered 0.50 spatial-error and
2× runtime limits. Raw precipitation rel-RMS was 2.06–2.18; the best-fit uniform amplitude
scale was only 0.175–0.191. The active-set method required 8–9 sparse solves per case, which
contributes to its poor measured runtime. This result rejects this particular thresholded
steady closure; it does not establish that every conservative steady formulation must fail.

**Disposition:** do not tune this sink further or open the held-out saves. The source field,
legacy ambient humidity, and the steady closure remain entangled, so these tuning-only
results cannot identify which physical assumption causes the spatial mismatch. A next
diagnostic should constrain the steady humidity/transport field against independent
vertically integrated water-vapor transport and column-water data before selecting another
rain law. If a steady closure still cannot meet the registered fidelity and runtime limits,
return to the owner with a concrete prognostic-atmosphere cost estimate; the existing
atmosphere advances one simulated day per tectonic step, so a multi-day humidity spin-up
would need an explicit cadence and runtime design. G3 is incomplete and G4 remains closed.

Private outputs are `~/mantle-bloom-captures/rfc0019/g3-steady-saturation-rh08.json` (SHA-256
`aadd401545ab2b0f436d3781049b3a589a3eb03f68bed9b5ecfdd8457fa4ecbf`) and
`g3-steady-saturation-rh06.json` (SHA-256
`8b8153751fafcbf1e321b2e613b4796273fdc2fa67a6fa8bdb532fc358bb3ec4`). Both used measurement
script SHA-256 `5aff7ef327ebd9ca0f1ce48ba580a037c1d465a3123ba711daca56caea8e3b4c`; measured
peak RSS was 1.86 GB and 1.71 GB, respectively. The first attempt to write the 80% capture
inside the sandbox was denied after its first case; the complete run was repeated with
authorized capture-directory access. Timings above are from the complete repeated run.

Reproduce each run on registered tuning saves A, D, and E with:

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_compare.py A D E \
  --density 4 --wind cfd --solver wind_arrows \
  --ocean-source-mapping bulk_flux --ocean-transfer-coefficient 0.003 \
  --land-proxy-multiplier 5 --effective-depth-km 2.25 --rain-days 5 \
  --wind-multiplier 0.5 --saturation-convection --convection-days 1 \
  --relative-humidity-threshold 0.8 \
  --out ~/mantle-bloom-captures/rfc0019/g3-steady-saturation-rh08.json
```

Repeat with `--relative-humidity-threshold 0.6` and output name `g3-steady-saturation-rh06.json`.

### G3 local orographic-lift screen (#913)

Issue #913 preregistered and tested one bounded local lift sink on the tuning saves A and cold
at density 1, using finite-volume transport and diagnostic and CFD winds. It uses 250 km
spherical terrain smoothing, positive windward ascent, a 2,250 m moisture scale, and a 5-day
maximum additional rainout time. Paired uniform controls use identical source fields, rainout,
masks and transport. Both saves close to machine precision. The lift sink strengthens the
paired windward/leeward rain ratio and modestly reduces best-fit precipitation RMS in both
wind modes; dry Köppen-B share also decreases. However, ≥60° land rain decreases in all four
save/wind combinations (9–53% on A and 13–38% on cold). The cold-land regression meets the
preregistered stop condition. Do not continue this local-lift law or open held-out saves; this
does not reject conservative moisture accounting or other orographic formulations. G3 remains
incomplete and the held-out acceptance gates and runtime limits are unchanged.

Private captures and maps are stored under `~/mantle-bloom-captures/issue913/`. The paired A
JSON SHA-256 is `e142bb03b043c5db71b532d920616ac8e8522d5646fcbb91657e41d596298751`; the cold
JSON SHA-256 is `ba60dd762d705c140f03cba38b4608a93b1a1fa6a5d80f6a837fa9d80bd0f68b`. The
reproducible preregistration, precise input identities, metrics, and limitations are in the
[issue #913 report](../issue913/README.md) and [lift preregistration](../issue913/lift-preregistration.md).

### G3 approved pattern checks: held-out result

On 2026-11-17 the owner approved replacing held-out precipitation RMS, latitude-band and
Köppen-share limits with the regional pattern checks in the
[G3 preregistration](g3-preregistration.md). The frozen `p = 0.25` ridge-rainout candidate
passed rain-belt, rain-shadow, cold-coast supply, closure and density-convergence checks on
the two held-out saves. It failed the dry-interior limit in 5 of 16 wind/solver/density
combinations; wind-arrow land/ocean shares also exceeded the retained 5-point held-out bound
at density 4. A density-4 offline candidate-update runtime proxy is below the 2× limit in all
wind/solver combinations (1.22–1.42×); production integration cost remains unmeasured. G3 has
not passed because of the pattern and held-out share failures.
See the [held-out pattern screen](g3-coarse-transport-plan.md#held-out-pattern-screen-for-frozen-ridge_025-candidate)
for metrics, exact save hashes and private result checksums, and the [runtime measurement](g3-coarse-transport-plan.md#density-4-runtime-proxy-for-frozen-ridge_025) for timing details. Do not tune against these acceptance-only saves.

### G3 final three-round tuning disposition: failed with owner exception

The owner authorized up to three tuning rounds to find the lowest-error eligible setting,
then an explicit G3 failure and exception if none passed. The rounds completed on 2026-11-18.
The best eligible setting, `alpha=1, L=1,500 km`, had a worst inland/coast rain ratio of
`1.515` on the cold CFD case, against the `0.6` limit. `alpha=2, L=1,500 km` improved this to
`1.300`, but retained only `77.5%` of baseline inland rain, below the `80%` safeguard, so it
was ineligible. No candidate was frozen or tested on new held-outs. The owner authorized RFC
continuation to G4 and implementation planning with this mismatch recorded as a modeling
simplification. **G3 is formally failed; the exception is not a pass.** The complete
preregistered outcomes and source fingerprint are in the [three-round tuning record](g3-three-round-plan.md#three-round-disposition-g3-failed-by-owner-exception).
