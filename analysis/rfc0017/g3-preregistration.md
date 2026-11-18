# RFC 0019 G3 preregistration

**Owner decision:** 2026-11-15, after G2 review. The owner authorized G3 and approved the
limits below before any G3 comparison. These limits are fixed for the registered held-out
cases; do not tune parameters or criteria against those cases.

On 2026-11-17 the owner approved replacing legacy climate-field agreement limits with the
regional acceptance checks below. Keep the legacy metrics as diagnostics. Freeze the
`p = 0.25` ridge-rainout candidate from the corrected tuning screen before evaluating the
held-out saves. Closure, density convergence and runtime limits remain gates.

## Registered cases

Use only the already registered tuning inputs for initial calibration: #771 windows A, D,
and E and the cold source at 200.8 Myr. Input names and SHA-256 digests are recorded in
[`README.md`](README.md#input-bytes-and-current-climate-snapshots). Compare both diagnostic
and CFD wind modes at climate densities 1 (90×180) and 4 (360×720). Compare the conservative
finite-volume and local wind-arrow stationary solvers at both resolutions so the density-4
runtime and memory cost can inform the solver choice.

The two saves in the [original held-out registry](README.md#held-out-registration-before-g0-replays)
were acceptance-only until their candidate evaluation; they have now been opened and are
spent. Do not use them for tuning or to accept a successor. Two replacement acceptance saves
and a bounded three-round tuning protocol were preregistered in
[`g3-three-round-plan.md`](g3-three-round-plan.md); verify those new save bytes before final
acceptance and do not inspect their climate fields while selecting parameters.

## Initial source mapping hypothesis

The current climate's surface inputs are humidity-equivalent proxies, not measured fluxes;
in particular, the ocean value is a target humidity rather than an evaporation ledger.
For this diagnostic comparison only, convert each proxy `q` (g/kg) into an effective column
water `W_proxy = (q / 1000) * rho_air * h_eff` (kg/m²), then an annual-mean source-rate
candidate `E_proxy = W_proxy / tau_source` (kg/m²/s). Use `rho_air = 1.2 kg/m³`, the RFC's
provisional `h_eff = 1.5–3.0 km` sensitivity, and `tau_source = 9.1 days` as the initial
case. Keep ocean target, lake/sea, river, and vegetation proxy contributions tagged and
reported separately. This conversion is a calibration hypothesis, not a measured physical
evaporation law or a production-ready surface-water debit. Any later source mapping that
changes this formula must be selected using tuning cases only and recorded before opening
held-outs.

Use `tau_rain = 9.1 days` initially for `P = W / tau_rain`. Retain the RFC's 5/9.1/18-day
rainout and 0.5/1/2 effective-wind sensitivity brackets on tuning cases. Treat surface wind
as a proxy for column transport: report `W*u` transport separately by mode, without calling
it observed IVT. Start with no humidity floor and no orographic enhancement; any enhancement
must be an additional rain sink debited once from available column water.

## Fixed acceptance limits

- Area-integrated source, transport, and rain closure per stationary climate update:
  relative residual no greater than `1e-8` of total source; column water and rainfall must be
  finite and nonnegative.
- Density-1 versus density-4 global rainfall changes by no more than 5% relative. Land and
  ocean rainfall shares change by no more than 5 percentage points.
- On each held-out save, at densities 1 and 4 and for diagnostic and CFD winds, land/ocean
  rainfall shares must remain within 5 percentage points of its current-source climate.
- For each pattern, evaluate eligible A/D/E tuning worlds and held-out worlds whose geography
  supports the measure, with the same masks and definitions at both densities and wind modes:
  - Equatorial ocean rain-belt center-to-flank mean ratio is at least `1.5`.
  - At least `60%` of eligible ridge crossings have wetter crests and drier lee sides, with
    positive median crest-minus-lee contrast gain.
  - Mean rain more than `1,500 km` inland is no greater than `0.6` times mean rain within
    `500 km` of a coast.
  - Cold land within `500 km` of a coast receives positive precipitation; there is no
    minimum cold-inland rate.
  Mark a pattern not evaluable when the registered geography lacks enough area or crossings
  for its measure; do not choose the eligible mask after seeing candidate rainfall.
- Preserve land/ocean rain shares, folded 15-degree bands, Köppen shares, and area-weighted
  precipitation relative RMS against the current-source field as diagnostic comparisons;
  they are no longer held-out acceptance limits.
- Median density-4 climate-update time is no greater than 2× the current-source climate
  update time, measured on the same world, wind mode, density, and process. Report each
  solver's peak RSS and timing; do not pool wind modes or resolutions.

## Measurement protocol

For each tuning world, wind mode, density, and solver: compute one current-source climate
snapshot and one conservative stationary moisture solution from the same terrain and wind
fields; use identical spherical cell areas and masks. If the saved lagged climate cache is at
another resolution, area-average its precipitation and temperature inputs and area-majority
its ocean mask onto the measurement grid before computing the snapshot. This preserves the
existing one-step vegetation input instead of silently disabling it at density 1. Record the
source-cache shape and whether it was resampled. Report land/ocean area-weighted rain,
folded 15-degree bands, Köppen group shares, precipitation relative RMS, each source-tag
amount, column transport, minimum column water, global and maximum-cell residuals, elapsed
solver time, total climate-update time, and peak RSS. Repeat timings enough times to report a
median after one warmup. Record exact save SHA-256, code commit, wind-state availability,
grid, solver, parameter values, and command for every measurement.

The candidate solve is stationary, so its per-update closure is the residual of the complete
source + horizontal transport − rain balance on that climate call. It does not establish
closure of a finite-duration atmosphere with changing storage; that limitation remains open
for later design review. First report the initial parameter case. Select any tuning changes
from the registered tuning cases, document why, then freeze the final parameter set before
loading either held-out save.
