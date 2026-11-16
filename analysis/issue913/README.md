# RFC 0019 G3 coarse-moisture follow-up

Issue [#913](https://github.com/adubey/mantle-bloom/issues/913), tracked under
[#898](https://github.com/adubey/mantle-bloom/issues/898). This is a diagnostic handoff;
no production climate behavior has changed.

## Principle and gate

The current diagnostic rain path violates design principle 6, **Account for material through
every process**: precipitation is not debited from carried atmospheric water. The invariant
is a finite, nonnegative stationary column-water field with source + transport − rain closure
no greater than `1e-8` of total source per update. G3 also retains its preregistered climate
structure and density-4 runtime limits; closure alone is not acceptance.

## Recovered inputs

The private G3 archive was verified against its SHA-256
(`606429f02fc93a7b966d2d1547953922005341110e75271d437ecead629a2d5d`); every member in
`SHA256SUMS` passed. The historical source snapshot is based on
`25e4a3132fe06a3b23305aacaf01c5f08631ba91`.

| Input | Name | SHA-1 | SHA-256 | Use |
|---|---|---|---|---|
| A | `seed910211954-221.6myr.bloom` | `0b67fd81b1be45fb6015baf62bcb520ce07de53f` | `1f167a8d9e0edf06c3c52b23598dc2f776adc2927f843beb992b5193ad30fd67` | Tuning; corrected-control screen complete |
| Cold | `mantle-bloom-seed317739385-200800000y.bloom` | `41075064222d1aa32a6aee9c5e032aa094dbaa71` | `c7e8c4094376b2b57158e36b1066cf529eb475ac0bbf54a6b888b391f294d344` | Tuning; pending first comparison |

Both pinned save identities are verified locally. Held-out saves remain sealed. Full JSON
captures, saves, and rain maps stay in the private issue capture directory, alongside its
verified `SHA256SUMS`; this report records only compact results.

## Corrected controls

Two diagnostic defects are fixed in `bin/debug/rfc0019_transport_prototypes.py` and
`bin/debug/rfc0019_g3_compare.py`:

1. The finite-volume builder flattened the negative-stride cell-area grid once per face. It
   now flattens once before building rates.
2. ERA5 class rainout was built in north-to-south map order while source, area, and wind were
   flipped into south-to-north solver order. The sink now uses the same flipped order.

On an 8×16 solver-order fixture, both solver schemes produced exactly the same rates, rain,
and residuals as the pre-fix uniform-sink control. An asymmetric land/ocean fixture confirmed
that the ERA5 class sink follows the matching hemisphere's surface mask.

## Window A initial control

The preregistered initial values were left unchanged (2.0 km effective depth, 9.1-day source
and rain timescales, unit wind and land-proxy multipliers). At 360×720, both solvers and wind
modes closed: global residuals were `0`–`5.4e-16`, and maximum cell residuals were `0`–`1.2e-17`.
Rain stayed finite and nonnegative after the solver's roundoff clipping. The water budget
closed, but this initial setting is not a candidate: density-4 best-fit precipitation rel-RMS
was `0.962`–`0.986`, far above the registered `0.50` limit. A single tuning setting cannot
decide whether the spatial error is caused by peaks or regional supply.

The instrumented diagnostic reported these median stage times (shared-machine screening only):

| Wind | Solver | Current update (s) | Diagnostic update (s) | Production-shaped estimate (s) | Duplicate sample (s) | Discarded zero-rain biome pass (s) | Rate build (s) | Solve + ledger (s) |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Diagnostic | Finite volume | 0.318 | 1.189 | 1.114 | 0.024 | 0.049 | 0.246 | 0.609 |
| Diagnostic | Wind arrows | 0.318 | 1.111 | 1.033 | 0.025 | 0.053 | 0.207 | 0.544 |
| CFD | Finite volume | 0.344 | 1.478 | 1.400 | 0.024 | 0.052 | 0.300 | 0.838 |
| CFD | Wind arrows | 0.344 | 1.239 | 1.162 | 0.023 | 0.054 | 0.247 | 0.661 |

The production-shaped estimate removes the repeated surface resample and the zero-rain biome
classification that the diagnostic later replaces with a classification from candidate rain.
It still exceeds the 2× budget in every A case. These timings overlapped other machine runs
and are not acceptance evidence. Peak process RSS was 1.25 GB.

The comparison used commit `25e4a3132fe06a3b23305aacaf01c5f08631ba91` and command shape
`backend/.venv/bin/python bin/debug/rfc0019_g3_compare.py <A-save> --density 4 --wind both
--solver both`. The initial capture is `g3-a-density4-corrected-controls.json`
(`fc71c33753f280cc33f3be5435429e8cab68eb567c1de994d9520c02b5e9d735`); the instrumented
rerun is `g3-a-density4-corrected-controls-instrumented.json`
(`772d20b38232f3c264e2561c493aee6b20d1f3505e326409289c97343d3b9911`). Numeric climate,
rain, source, and closure fields matched between those two captures.

The multiscale capture `g3-a-density4-multiscale.json` has SHA-256
`5ab94e116ea0bfdc0cd798fd4c5c45cbfadb46b15ae0fe8b86df63b2f12870b8`. Area-averaging both
density-4 rain maps to 90×180 reduces raw relative RMS from `3.67` to `1.68` for diagnostic
wind and `3.72` to `2.03` for CFD wind. Best-fit RMS improves from `0.962` to `0.849` and
from `0.965` to `0.893`, respectively. Aggregation removes some narrow-peak error, but the
remaining coarse-scale mismatch means peak width alone does not explain the initial case.

The scale and terrain-exposure summaries below use the wind-arrow solver; the finite-volume
solver has its own global and local closure results in the capture.

At this setting, land median rain is `2.3 mm/yr` (diagnostic) and `20.2 mm/yr` (CFD), versus
`156` and `357 mm/yr` in the corresponding current-source maps. High-latitude (≥60°) land
means are `116` versus `260 mm/yr` for diagnostic wind, and `289` versus `646 mm/yr` for CFD.
The Köppen dry-group share rises from `26.9%` to `46.5%` (diagnostic) and `53.3%` to `93.3%`
(CFD). These are initial-control results, not a selected or tuned candidate, but they point to
a regional land-water deficit in addition to unresolved peaks.

For an additional spatial check, land cells were ranked by the terrain-ascent projection of
the same wind used by each solve; the positive and negative extreme masks cover 4.4–5.6% of
land area. Current-to-candidate windward/leeward rain ratios fell from `2.98` to `1.78` in
diagnostic mode and `1.88` to `1.36` in CFD mode. This diagnostic mask is exploratory and is not
an acceptance metric; its exact cuts and areas are in the capture.

## Cold-window corrected control

The pinned cold save was compared at density 4 with the same initial parameters and script
commit. Capture `g3-cold-density4-corrected-controls.json` has SHA-256
`85e629b645c024a773eabbd64e9976836ef579d8cd181997d9d34f81c2535736`. Global closure was
`1.2e-16`–`3.0e-14`; peak process RSS was 1.65 GB. Timings shared the machine with other
simulations and are screening only.

| Tuning save | Wind | Native best-fit RMS | 90×180-equivalent RMS | Land mean rain, current → candidate (mm/yr) | ≥60° land rain, current → candidate (mm/yr) | Köppen B share, current → candidate | Windward/leeward ratio, current → candidate |
|---|---|---:|---:|---:|---:|---:|---:|
| A | Diagnostic | 0.962 | 0.849 | 468 → 296 | 260 → 116 | 26.9% → 46.5% | 2.98 → 1.78 |
| A | CFD | 0.965 | 0.893 | 650 → 341 | 646 → 289 | 53.3% → 93.3% | 1.88 → 1.36 |
| Cold | Diagnostic | 0.966 | 0.870 | 984 → 825 | 348 → 445 | 21.8% → 53.8% | 2.88 → 1.45 |
| Cold | CFD | 0.969 | 0.898 | 1,047 → 924 | 456 → 784 | 40.6% → 84.6% | 2.70 → 1.20 |

The cold case confirms the coarse-scale error does not disappear after aggregation, but its
high-latitude land rain increases under the initial candidate while its dry-biome share and
windward/leeward contrast worsen. A has less high-latitude rain and lower land rain. The sign
change in the high-latitude response means the aggregate mismatch is not a uniform land-water
shortage; sources, terrain uplift, temperature, and regional routing still need separation.
Neither initial case meets the precipitation-shape or density-4 runtime limits.

## Density-1 versus density-4 consistency

The density-1 corrected captures are `g3-a-density1-corrected-controls.json`
(`fdb7ab357f728449a17bf93a00818f8c4cbe84adf8db4ef0906b1bfa3a587315`) and
`g3-cold-density1-corrected-controls.json`
(`7d232d41fd4fbc97fdadb33dffabcc48670129237df491484b296b59a5d34e17`). The density-4 runs
used the same source mapping and parameters. All initial-control global-rain changes and
land/ocean share shifts were within the registered 5% / 5-point density consistency limits:

| Tuning save | Wind | Current global change | Current land-share change | Candidate global change | Candidate land-share change |
|---|---|---:|---:|---:|---:|
| A | Diagnostic | +1.39% | +2.36 pp | +0.10% | −0.54 pp |
| A | CFD | +2.87% | +3.54 pp | +0.02% | −1.09 pp |
| Cold | Diagnostic | +2.59% | +2.90 pp | +0.08% | −0.78 pp |
| Cold | CFD | +1.81% | +2.66 pp | +0.02% | +0.10 pp |

Candidate deltas shown are for wind arrows; finite-volume deltas also stayed within the
registered limits. The updated diagnostic records each scheme's scale metrics separately.
These initial-setting results establish resolution consistency for totals and surface shares,
not spatial acceptance; multiscale best-fit RMS remains `0.85`–`0.90` after area averaging.

## Window A paired lift screen

The preregistered lift law was compared with its uniform-sink control on the same A save at
90×180, using finite-volume transport and both wind modes. The capture is
`g3-a-density1-lift-comparison.json` (SHA-256
`e142bb03b043c5db71b532d920616ac8e8522d5646fcbb91657e41d596298751`). Four common-scale
rain maps are stored alongside the private capture. Peak RSS was 0.99 GB; runtime overlapped
two other simulations and is not treated as a performance measurement.

| A wind | Best-fit RMS, uniform → lift | Windward/leeward ratio, uniform → lift | ≥60° land rain, uniform → lift (mm/yr) | Köppen B share, uniform → lift | Lift share of rain | Global closure (uniform / lift) |
|---|---:|---:|---:|---:|---:|---:|
| Diagnostic | 0.861 → 0.833 | 1.40 → 3.31 | 114 → 53 | 43.2% → 40.6% | 25.2% | `1.3e-16` / `1.3e-16` |
| CFD | 0.866 → 0.844 | 1.39 → 2.89 | 89.9 → 81.5 | 89.6% → 87.2% | 26.7% | `0` / `0` |

The lift law produces its predicted stronger terrain contrast and a modest reduction in coarse
pattern RMS on A, while reducing dry-group share. It also lowers high-latitude land rain by
about 53% in diagnostic wind and 9% in CFD wind, a cold-land regression on this save. This
single tuning case is not enough to accept the local-lift path; the preregistered cold tuning
comparison is needed to decide whether the effect is consistent or save-dependent. No
held-out save was opened.

## Cold-window paired lift screen

The same paired run was made on the pinned cold tuning save with identical density, finite-
volume transport, wind modes, source settings, and lift law. Capture
`g3-cold-density1-lift-comparison.json` has SHA-256
`ba60dd762d705c140f03cba38b4608a93b1a1fa6a5d80f6a837fa9d80bd0f68b`. Four common-scale
maps are in the private captures directory. Peak RSS was 1.53 GB; runtime overlapped other
simulations and is not a performance measurement.

| Cold wind | Best-fit RMS, uniform → lift | Windward/leeward ratio, uniform → lift | ≥60° land rain, uniform → lift (mm/yr) | Köppen B share, uniform → lift | Lift share of rain | Global closure (uniform / lift) |
|---|---:|---:|---:|---:|---:|---:|
| Diagnostic | 0.874 → 0.847 | 3.84 → 9.49 | 527 → 326 | 51.1% → 44.3% | 29.4% | `1.2e-16` / `0` |
| CFD | 0.874 → 0.847 | 2.14 → 5.91 | 299 → 261 | 78.8% → 69.8% | 29.5% | `1.2e-16` / `0` |

This confirms the lift sink improves the predicted terrain contrast and modestly improves the
coarse-scale pattern on both tuning worlds. Dry-group share also falls in all four wind/save
pairs. But high-latitude land rain falls in all four: 53% and 9% on A, and 38% and 13% on
cold, for diagnostic and CFD wind respectively. The preregistered cold-land stop condition is
therefore met. Stop this local-lift formulation; it redistributes too much rain toward ascent
zones at the expense of high-latitude land delivery. This does not reject a conservative
moisture budget or every orographic process representation. No held-out save was opened, and
these density-1 screens do not satisfy G3 acceptance or runtime gates.

## Next

The comparison script now records native and density-1-equivalent rain RMS, area-weighted
intensity quantiles, hemisphere-specific bands, high-latitude land rain, Köppen dry share, and
terrain-ascent/descent rain pairs using shared masks. The single lift-aware arm is frozen in
[`lift-preregistration.md`](lift-preregistration.md): 250 km spherical terrain smoothing,
2,250 m moisture scale, and a 5-day cap on additional rainout. Initial controls on A and cold
are complete at both densities. The paired lift screens on A and cold both improved terrain
contrast and RMS, but reduced high-latitude land rain in diagnostic and CFD wind. The
preregistered stop condition is met; do not extend this local-lift law or tune against
held-outs. Any revised G3 design must protect cold high-latitude delivery while preserving the
moisture ledger, and requires owner review before held-out evaluation.
