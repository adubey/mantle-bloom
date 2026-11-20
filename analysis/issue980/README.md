# Issue 980: RFC 0021 G1 close-out, best-case check and 0.5 Myr repeat

Finishes G1 of [RFC 0021](../../rfc/0021-plate-force-calibration.md#diagnosis-gates) (tracking #957), after #975
(`analysis/issue975/`). Same harness and conditions: plates-only, rift lines at the rifting baseline (`mark_power` 2.8),
200 Myr, five stage 1 worlds (`analysis/issue975/stage1-seeds.json`), sampled every 10 Myr. Scores are in `scores/`;
`bin/debug/tabulate_g1_scores.py analysis/issue980/scores <names>` prints them.

## Best-case check (the RFC's falsification test)

The combined adopted corrections `b123` (B.1 + B.2 + B.3), brake off, `f_pull` 1.0, `ASTHENOSPHERE_VISCOSITY_PA_S`
2×10¹⁹, at the four slab-favourable corners of G2's ranges (module overrides applied after world generation).
`settings.json` is `analysis/issue975/settings.json`.

| | `b123` brake off (#975) | μ_slab 10²⁰, L_sink 150 km | 10²⁰, 660 km | 10²¹, 150 km | 10²¹, 660 km |
|---|---:|---:|---:|---:|---:|
| **T2 ratio (+10 Myr)** | 1.59 | **1.84** | 1.30 | 1.53 | 1.22 |
| T2 same time | 1.63 | 1.93 | 1.56 | 1.50 | 1.25 |
| Slab-attached share of plate area | 56% | 74% | 94% | 64% | 81% |
| T1 NNR speed, cm/yr | 4.25 | 9.24 | 21.63 | 5.74 | 7.88 |
| T3 pull share of drive | 60.0% | 58.1% | 83.7% | 57.1% | 85.1% |
| T4 MPa | 1.52 | 0.99 | 1.03 | 0.90 | 1.23 |
| T5 N/m | 1.90e+12 | 2.12e+12 | 2.50e+12 | 2.02e+12 | 2.16e+12 |
| T7 fastest/slowest world | 1.54 | 1.64 | 1.15 | 1.33 | 1.13 |
| Net slab work | +29.0% | +76.0% | +49.4% | +18.6% | +5.4% |
| Sunk-slab drag share of resistance | 37% | 10.0% | 38.0% | 41.9% | 78.8% |
| Splits per world per 100 Myr | 11.8 | 3.4 | 8.1 | 1.9 | 21.5 |
| Jammed collision-time drops (T8) | 3 | 6 | 14 | 0 | 9 |
| Plates (mean) | 17.5 | 14.7 | 25.4 | 13.1 | 27.5 |
| RFC 0018 targets met (of 10) | 6 | 6 | 6 | 6 | 6 |

Solver closure is at most 5.3×10⁻¹⁵ of the driving torque in every run and every sample converged.

**Verdict: no corner reaches 2.0; the best is 1.84.** By the RFC's falsification criterion ("Slab forces cannot drive
these plates") the premise fails, and the RFC says to revise toward alternative D (the column) or C (the flow). Per
the issue, no criterion was changed and G2 was not started: the verdict is the owner's to act on.

What the table shows, for that decision:

- **T2 does not rise with pull; it falls once most plates have a slab.** The slab-attached share of area goes from 56%
  to 94% as pull grows, so the "other" group (the denominator) shrinks to a few slow plates and the ratio compares
  nearly all plates with a remainder. At 10²⁰/660 km the plates move at 21.6 cm/yr (above the 15 cm/yr
  `MAX_PLATE_RATE`, which the quasi-static solve evidently does not clamp to; not investigated here) and pull is 84% of
  the drive, yet T2 is 1.30.
- So at these speeds T2 as defined stops measuring slab drive. Whether to change T2's definition is the owner's
  decision; this PR does not.
- The 10²¹ corners behave as the simple model said: slab drag still dominates (42% and 79% of the resistance) and
  speed rises little.
- T1 reaches its 3–5 band only at defaults (`b123`: 4.25); every corner overshoots it.

### Predictions against results

Posted on #980 before the runs. T2 was predicted 2.0–2.7 / 2.2–3.2 / 1.3–1.7 / 1.4–1.8 for the four corners in the
table's order and came out 1.84 / 1.30 / 1.53 / 1.22: both 10²⁰ corners missed (the 660 km one by a wide margin), the
10²¹ ones were close (1.53 inside, 1.22 under). T1 predictions were 6–9, 8–12, 4.5–6, 5–7 against 9.24, 21.6, 5.7,
7.9. The named risk (the slab-attached share growing and diluting the ratio) dominated. Splits were predicted to rise
50–150% at the 10²⁰ corners; they fell at 150 km (3.4) and fell by a third at 660 km (8.1).

## 0.5 Myr step repeat of `b123`, brake off

Worlds 106533960, 1097624328 and 1397817354, 200 Myr, `--step-myr 0.5`, against the same worlds at 1 Myr from the
`--world-checks` pass below (which reproduces #975's checked-in `b123-brake-off` score exactly).

| | 1 Myr | 0.5 Myr | Change | Criterion |
|---|---:|---:|---:|---|
| T1 NNR cm/yr | 4.14 | 4.27 | +3.3% | 15% |
| T2 ratio (+10 Myr) | 1.70 | 1.62 | −4.9% | 20% |
| T3 pull share | 57.4% | 56.8% | −1.0% | 20% |
| T4 MPa | 1.47 | 1.33 | −9.1% | 20% |
| T5 N/m | 1.84e+12 | 1.89e+12 | +3.1% | 20% |
| T7 | 1.54 | 1.53 | −0.3% | 20% |
| Splits per world per 100 Myr (not counted, #969) | 10.7 | 11.5 | +8% | |

**The step check passes.** (The prediction posted on #980 compared against an older capture of this setting left in
`issue975/b123-brake-off`; the checked-in score is the one reproduced here, so the posted 1 Myr figures
(T1 4.02, T2 1.67) are superseded by these.)

## Ledgers and alternative D (G1 items left open by #979)

`bin/debug/g1_world_checks.py`, run with `calibrate_plate_forces.py --world-checks`, reads the world's ledgers and
alternative D's effect at every sample. One pass of `b123-brake-off` on the five worlds
(`scores/b123-brake-off-worldchecks.json`). It reproduces #975's checked-in score for that setting to every digit
(T1 4.2479, 118 splits, all targets), so reading the ledgers and turning on `debug_diagnostics` for the Hm ledger leave
the run unchanged.

**Ledgers**, worst relative closure error over every sample of the five worlds:

| Ledger | Worst gap |
|---|---:|
| Continental material | 2.1×10⁻¹⁵ |
| Craton | 2.9×10⁻¹⁶ |
| Mobile cover | 0 |
| Hm (against the start-of-run inventory, signed residual of its identity) | 6.3×10⁻⁴ |

Three ledgers close to round-off. The Hm ledger does not: its residual is up to 0.06% of the initial Hm inventory, as
`hm_ledger.py` expects of a quantity that is thickness, not conserved mass, but it is not a numerical-precision
closure. RFC 0021 says the ledgers close before any correction becomes a default, and G2's preregistration should name
which Hm residual is acceptable. This probe does not decide that.

**Alternative D** (every oceanic cell's Hm set to B.3's age column, elevations re-isostasied, sea level re-solved for
the same water; geometric, so cut-off basins aren't separated), at 200 Myr, mean of five worlds:

| | Now | With D |
|---|---:|---:|
| Land fraction | 30.2% | 32.4% |
| Sea level change | | −570 m |
| Mean ocean depth | 3.93 km | 4.06 km |
| Mean oceanic Hm | 33.1 km | 62.0 km |

92% of oceanic area would deepen. Depth against age, pooled (m below sea level, GDH1 from Stein & Stein 1992):

| Age, Myr | Area, 10⁶ km² (5 worlds) | Now | With D | GDH1 |
|---|---:|---:|---:|---:|
| 0–10 | 285 | 3906 | 3415 | 3346 |
| 10–30 | 353 | 4252 | 4029 | 4176 |
| 30–60 | 304 | 4385 | 4641 | 4909 |
| 60–100 | 206 | 4372 | 4827 | 5361 |
| 100–200 | 288 | 4377 | 4804 | 5608 |
| over 200 | 122 | 3376 | 3673 | 5651 |

Today's floor is nearly flat with age (3.9–4.4 km). D follows GDH1 to within 0.15 km under 30 Myr, but stays shallower
beyond it by 0.27 km (30–60 Myr), 0.53 (60–100), 0.80 (100–200) and 2.0 km (over 200 Myr, floor that is undated or
anomalous and that D can't age). D moves the sea level by about 570 m and the land fraction by 2.2
points, so it would move every surface ledger the RFC says it would. These numbers are a first-order static estimate,
not a run with D on.

## Not done

- T6 (not in scope); G2; any change to targets, bars or defaults.
- T4 is not comparable across the B.4 settings (#975).
- The raw captures: `~/mantle-bloom-captures/issue980/` (corners, `step0.5/`, `worldchecks/`). The old
  `issue975/b123-brake-off` captures are not the run that produced the checked-in score; use `worldchecks/` for it.
