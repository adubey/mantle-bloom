# #842 round 3: preregistration

Written 2026-11-13, before any round 3 run. The decision rule and thresholds below are
frozen from this commit. A change after a round 3 result is a threshold change, and #842
says a failed round 3 gets none: the result goes to the owner.

Question: can the trend-adaptive march (longest substep **1 kyr**) replace fixed-250 yr D2
as the production default, judged against RFC 0015's own criteria?

## Amendments before any run (2026-11-13)

The first draft (PR #870's first commit) was written before checking what `origin/main`
can run. Four things changed once the harness was built, all before any round 3 run:

1. **No `origin/main` arm.** `origin/main` now has the D2 march, so it *is* the `d2` arm.
   The pre-march solve the old runtime threshold measured against is gone. The old runtime
   check is estimated from #789's pre-march step times (row 5 below).
2. **E3 can't rerun #741's drop steps as they were.** Those came from the G0 prototype
   (`ea9a187`) on drift draws 1 and 2, and main has neither the prototype march nor the
   draw offset. E3 is redefined below: a drop-signature gate on every replay step, plus
   ensembles on main's own largest-loss step, with a positive control.
3. **The "meltwater ledgers to round-off" half of the ledger criterion can't be measured.**
   RFC 0015 D8's meltwater ledgers aren't implemented on main. The replays check #789's
   land-water ledger (≤ 0.1% of its external inputs) in its place.
4. **Step-gap denominator.** The gap is |V(0.25 Myr) − V(1 Myr)| / V(0.25 Myr) at the end
   of the run, for volume and for area (nodes with glacier depth > 0), and the larger of the
   two gates. #765's 9.71% used the 0.25 Myr run's volume.

## Arms

| Arm | March |
|---|---|
| `d2` | Fixed 250 yr substeps, pair-mean settle test. Today's production march (`MARCH_COST_PRESETS["d2"]`). |
| `trend-1k` | `MARCH_COST_PRESETS["trend-1k"]` = `MarchCost(trend=True, trend_max_substep_years=1000)`. Window 30 substeps, hold when the projected remaining change is < 0.1% (`trend_tol` 1e-3), double the substep while the change shrinks (`adaptive_growth` 2), back to 250 yr when the ratio exceeds 1.5. D2's pair test still runs beside it. |
| `d2-10k` | `d2` with `--substep-years 10000`. E3's positive control only: #741's failure form. |

The `trend-1k` knobs are frozen from this commit. E5 uses the held-out set and may not be
used to choose any of them.

## Seed sets (RFC 0015 "Chosen (2026-10-23, #548)", unchanged)

Checksums match `analysis/issue548/` (checked 2026-11-13).

- Diagnosis: seed 997271774 at 87.4 Myr, seed 649658380 at 6.4 Myr, seed 374464662 at
  132.6 Myr (replayed with ice ages at an 8 Myr period).
- Held out: seed 447276134 at 116.4 Myr (ice ages at an 8 Myr period), seed 628977350 at
  226.7 Myr, seed 295569319 at 180.1 Myr.

## Runs

Every run is one `bin/debug/probe_822_march_cost.py` process, one at a time under the
`compute-budget` skill. Runs that time something (all of E1, E2, E5, E6) start only when
`compute_budget.py check` lists no other run, since overlapping runs' times aren't
comparable. The ice-age saves get `--ice-age-period-years 8e6`.

| | What | Saves | Runs |
|---|---|---|---|
| E1 + E2 | 10 Myr replays: 40 × 0.25 Myr and 10 × 1 Myr, both arms | Diagnosis | 12 |
| E5 | The same | Held out | 12 |
| E3 | Opening-ice ensembles, one step, 3 arms × 13 members | Seed 997271774 | 78 (2 steps) |
| E4 | `d2` at substeps 1000, 500, 250, 125 yr, one step at 0.25 Myr and one at 1 Myr, from the save | Diagnosis | 24 |
| E6 | Median step time from E1 + E2; plus #822's gen30 1 Myr monster step and gen5 0.25 Myr step, both arms, report-only | Diagnosis, gen5, gen30 | 4 extra |

E1 and E5's ledger checks read the first 10 steps of each replay. Swings and step gaps read
all of it.

**E3, redefined.**

- **E3a (gating):** on every step of every `trend-1k` replay in E1, E2 and E5, no single
  substep moves global ice volume by more than 20% (`worst_substep_move` in the march log).
  That is #741's drop signature: one substep took the volume from 13.1 to 5.9 M km³.
- **E3b (ensembles):** two steps of seed 997271774, taken from its `d2` 0.25 Myr replay:
  the step ending at 97.4 Myr (#741's calendar step) and the step with the largest
  step-to-step ice loss after 2 Myr (the next largest, if that's the same step). Each step's
  opening state is an exact save from a rerun of that replay (`--write-save-after K`). Each
  step runs 13 members per arm: the control and 12 uniform opening-ice scalings of ±0.01,
  0.02, 0.05, 0.1, 0.2 and 0.5% (#741's ±0.01–0.5%).
  - A member **drops** when it closes below 88% of the unperturbed `d2` control. #741 used
    11.5 M km³ against about 13.
  - A step **counts as a drop step** when at least 8 of 12 `d2-10k` members drop (#741: 9 of
    12 at 10 kyr).
  - On a drop step, E3b passes when `trend-1k` drops no more members than `d2`, and none if
    `d2` drops none (the #741/#765 invariant: 0 of 12 at 250 yr and at 1 kyr).
  - A step that isn't a drop step is reported and doesn't gate, because the invariant can't
    be tested there. If neither step is a drop step, E3 rests on E3a, and the report says so.

## Thresholds: old and new

The owner relaxed five thresholds on 2026-11-13, before round 3 ran. The old values stay
here, and the final evaluation reports each result against both.

| # | Criterion | Old | New (gating) |
|---|---|---|---|
| 1 | Ledger residual (E1, E5, first 10 steps) | ≤ 0.1% of the step's accumulation | `trend-1k`'s residual ≤ 0.1% of accumulation **or** ≤ `d2`'s residual on the same step. The count of steps over 0.1% is reported for both arms. The land ledger stays ≤ 0.1% (amendment 3). |
| 2 | Step robustness, 0.25 vs 1 Myr, 10 Myr (E2, E5) | < 10% (held-out acceptance) and < 20% (G0, diagnosis) | < 20% on every save |
| 3 | Substep convergence (E4) | Within 1% (G1) | Reported only, no gate |
| 4 | Agreement with `d2` (E1, E5) | Within 0.87% (round 2 result, never a criterion) | Reported only, no gate |
| 5 | Runtime (E6) | Median step ≤ +15% vs `origin/main`, then the pre-march solve | Median step ≤ +15% vs `d2` on the same save and step size, on every diagnosis save at both step sizes |

Old row 5 is estimated, not measured: `trend-1k`'s 0.25 Myr median over #789's pre-march
medians (7, 6 and 8 s per step on the three diagnosis saves). Those ran on another day,
so the ratio is indicative only. #789 found `d2` itself at +76% to +213% on them.

Unchanged and gating: no unforced swings (≤ 25% step-to-step after 2 Myr, saves without ice
ages, all replay steps), E3, held-out discipline.

Row 1 is a non-regression gate, because #852 shows the ice ledger missing 1e-3 on 183 of
200 steps at 100 kyr and longer for `d2` and `trend-adapt` alike. The absolute fix is #852's.
This round 3 shows only that the trend march isn't worse.

## Decision rule

Round 3 **passes** when all of these hold:

1. E1 and E5: gate 1 and the land ledger hold on every checked step, and no step-to-step
   change in global ice volume exceeds 25% after the first 2 Myr (saves without ice ages).
2. E2 and E5: the step gap is below 20% on every save for `trend-1k`.
3. E3: E3a holds, and E3b passes on every step that counts as a drop step.
4. E6: `trend-1k`'s median step time is ≤ 1.15× `d2`'s on every diagnosis save at both step
   sizes.

If any condition fails, the result goes to the owner. No flip and no threshold change.
E4, the agreement with `d2`, either arm's absolute ledger count, and the E6 report-only
timings are reported and don't gate.

`bin/debug/score_842_round3.py` computes the verdict under both the new and the old
thresholds from the run JSONs. Under the old thresholds, row 5 uses the estimate above.

## Final evaluation must include

- The five old-to-new threshold changes above, each result stated against the old and the
  new threshold, including every case where the old one would have failed and the new one
  passes.
- Whether the result would change if the owner kept the old thresholds.
- The four amendments above, and for E3 whether either step was a drop step.
