# #842 round 3: preregistration

Written 2026-11-13, before any round 3 run. The decision rule and thresholds below are
frozen from this commit. A change after a round 3 result is a threshold change, and #842
says a failed round 3 gets none: the result goes to the owner.

Question: can the trend-adaptive march (`trend-adapt`, longest substep **1 kyr**) replace
fixed-250 yr D2 as the production default, judged against RFC 0015's own criteria?

## Arms

| Arm | March |
|---|---|
| `d2` | Fixed 250 yr substeps, pair-mean settle test. Today's production march. |
| `trend-1k` | `MarchCost(trend=True, trend_max_substep_years=1000)`. Window 30 substeps, hold when the projected remaining change is < 0.1% (`trend_tol` 1e-3), double the substep while the change shrinks (`adaptive_growth` 2), back to 250 yr when the ratio exceeds 1.5. Every other knob as in `MARCH_COST_PRESETS["trend-adapt"]`. |
| `origin/main` | The pre-march solve, for the runtime baseline in E6 only. |

The `trend-1k` knobs are frozen before E5 runs. E5 uses the held-out set and may not be
used to choose any of them.

## Seed sets (RFC 0015 "Chosen (2026-10-23, #548)", unchanged)

- Diagnosis: seed 997271774 at 87.4 Myr, seed 649658380 at 6.4 Myr, seed 374464662 at
  132.6 Myr (replayed with ice ages at an 8 Myr period).
- Held out: seed 447276134 at 116.4 Myr (ice ages at an 8 Myr period), seed 628977350 at
  226.7 Myr, seed 295569319 at 180.1 Myr.

## Evaluations

| | What | Saves |
|---|---|---|
| E1 | 10-step replays at 0.25 Myr and 1 Myr steps, both arms | Diagnosis |
| E2 | End-of-run volume and area at 0.25 vs 1 Myr steps, both arms, 10 Myr | Diagnosis |
| E3 | 12-member opening-ice ensemble on both #741/#765 drop steps, `trend-1k` | Seed 997271774 |
| E4 | `d2` substep-convergence reference: longest substep 1000, 500, 250, 125 yr | Diagnosis |
| E5 | 10-step replays at 0.25 Myr and 1 Myr steps, both arms, plus the E2 step gap at 10 Myr | Held out |
| E6 | Median step time at density 4, plus #822's monster-test and CI-gate timings | Diagnosis |

## Thresholds: old and new

The owner relaxed five thresholds on 2026-11-13, before round 3 ran. The old values stay
here, and the final evaluation reports each result against both.

| # | Criterion | Old | New (gating) |
|---|---|---|---|
| 1 | Ledger residual (E1, E5, every step) | ≤ 0.1% of the step's accumulation | `trend-1k`'s residual ≤ 0.1% of accumulation **or** ≤ `d2`'s residual on the same step. The count of steps over 0.1% is reported for both arms. Meltwater ledgers to round-off is unchanged. |
| 2 | Step robustness, 0.25 vs 1 Myr, 10 Myr (E2, E5) | < 10% (held-out acceptance) and < 20% (G0, diagnosis) | < 20% on every save |
| 3 | Substep convergence (E4) | Within 1% (G1) | Reported only, no gate |
| 4 | Agreement with `d2` (E1, E5) | Within 0.87% (round 2 result, never a criterion) | Reported only, no gate |
| 5 | Runtime (E6) | Median step ≤ +15% vs `origin/main` | Median step ≤ +15% vs `d2` on the same saves and steps. The ratio to `origin/main` is reported. |

Unchanged and gating: no unforced swings (≤ 25% step-to-step after 2 Myr, ice ages
never, E1 and E5), E3's 0-of-12 invariant, held-out discipline.

Row 1 is a non-regression gate, because #852 shows the ice ledger missing 1e-3 on 183 of
200 steps at 100 kyr and longer for `d2` and `trend-adapt` alike. The absolute fix is #852's.
This round 3 shows only that the trend march isn't worse.

## Decision rule

Round 3 **passes** when all of these hold:

1. E1 and E5: gate 1 holds on every step of every replay, and no step-to-step change in
   global ice volume exceeds 25% after the first 2 Myr (ice ages off).
2. E2 and E5: the step gap is below 20% on every save for `trend-1k`.
3. E3: 0 of 12 members drop on both drop steps.
4. E6: `trend-1k`'s median step time is ≤ 1.15× `d2`'s.

If any condition fails, the result goes to the owner. No flip and no threshold change.
E4 and the agreement with `d2` are reported and don't gate. Neither does either arm's
absolute ledger count.

## Final evaluation must include

- The five old-to-new threshold changes above, each result stated against the old and the
  new threshold, including every case where the old one would have failed and the new one
  passes.
- Whether the result would change if the owner kept the old thresholds.
