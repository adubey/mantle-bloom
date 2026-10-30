# Issue #696: capped lake siblings recirculating surplus

The fix works: the capped-pair cycle is gone, the pass limit is no longer reached, and the first
step's final export, ocean overflow and lake volume match `origin/main` exactly. One preregistered
threshold failed because it was set wrongly; see [Results](#results). Final export reaches the
ocean in none of the 15 steps on either arm, and nearly all of it then leaves the model.

## Preregistration

**Change.** `_advance_prognostic_body` stops resending surplus around a closed pair: a basin held
at an artificial cap below the saddle it shares with a sibling that is full to that saddle. In
the step's last integration slice the return flow waits for the final export, at the same leaf
sink as before. No constant, cap or destination changes.

**Replays.** Windows A, B and C, five 100 kyr steps from each registered #680 checkpoint with
`ice_age_period_years = 0`, on a detached `origin/main` tree and on this branch. One replay at a
time. The solve is deterministic, so each arm runs once.

| Window | Checkpoint age (Myr) | SHA-256 |
|---|---:|---|
| A | 221.6 | `e4b0b4f09b87584aa7f77045f0f079f2309de0132049915a9aeb5589add27340` |
| B | 222.1 | `86fb5b7d2044846ac751906ff996184c7abe5454c93b39033bd82338e8749544` |
| C | 201.0 | `041b7430a6cf8867f2e2c76f18c3328e440add585ee4a2e4b3055e7766c2217e` |

**Primary metrics and decision thresholds.**

1. Steps reaching the pass limit: baseline about 11 of 15; pass if 0 of 15.
2. Internal transfer ÷ external supply: baseline up to about 500; pass if at most 2 in every step.
3. First-step final export and ocean overflow, before any silt difference feeds back: pass if
   each is within 0.1% of baseline in all three windows. A larger gap needs an explanation.
4. Normalized ledger residual stays at or below the baseline's 1.4e-12 order.

**Expected direction.** Passes fall to a few tens. Internal transfer falls by at least 90% in
steps with large cycles. Silt falls strongly in those steps, because silt reads supply that
counted every resend. That changes rims and water from step 2, so steps 2–5 are expected to
drift. This is not a fault in the fix; #693 owns the rest of the silt count.

**Controls and secondary measures.** Window B at 222.5 Myr ended below the limit with no final
export, so it should barely change. Also recorded: lake volume, anchor component `ΔV`, silt, step
time, where final export ends up, and each pair's parent outlet.

**Variability.** Deterministic replay of one world per window. A and B share a world, so they are
not independent. Steps 2–5 are a one-trajectory comparison, not a distribution.

## Results

Baseline `d709f2b5`, branch `6393180e`. Compact results are in
[`replay.json`](replay.json); the script is `bin/debug/measure_lake_cycle_export.py`.
The branch was then rebased over #697 and #695, which change no lake code or checkpoint-replay path.

```bash
backend/.venv/bin/python bin/debug/measure_lake_cycle_export.py <checkpoint.mbworld> \\
  --backend <tree>/backend --anchors analysis/issue696/anchor-window-a.json \\
  --myr 0.5 --out <scratch-output>
```

Each cell is baseline → branch. Volumes are m³; internal transfer is a multiple of external supply.

| Window, age (Myr) | Passes | Internal ÷ external | Final export | Ocean overflow | Lake volume | Silt |
|---|---:|---:|---:|---:|---:|---:|
| A, 221.7 | 991 → 61 | 148 → 22 | 7.35e+16 → 7.35e+16 | 2.47e+17 → 2.47e+17 | 1.82e+15 → 1.82e+15 | 3.87e+14 → 1.09e+14 |
| A, 221.8 | 992 → 74 | 155 → 16 | 8.35e+16 → 8.45e+16 | 2.58e+17 → 2.58e+17 | 2.46e+15 → 2.46e+15 | 4.7e+14 → 8.14e+13 |
| A, 221.9 | 1004 → 69 | 148 → 24 | 7.99e+16 → 7.82e+16 | 2.95e+17 → 2.95e+17 | 2.14e+15 → 2.13e+15 | 4.14e+14 → 1.27e+14 |
| A, 222.0 | 1033 → 76 | 150 → 23 | 7.7e+16 → 7.17e+16 | 2.63e+17 → 2.63e+17 | 2.07e+15 → 2.83e+15 | 4.21e+14 → 1.19e+14 |
| A, 222.1 | 989 → 82 | 415 → 23 | 2.83e+17 → 2.8e+17 | 2.61e+17 → 2.61e+17 | 2.24e+15 → 2.8e+15 | 1.24e+15 → 1.33e+14 |
| B, 222.2 | 1007 → 69 | 538 → 24 | 2.9e+17 → 2.9e+17 | 2.23e+17 → 2.23e+17 | 2.44e+15 → 2.44e+15 | 1.25e+15 → 1.08e+14 |
| B, 222.3 | 991 → 63 | 211 → 24 | 1.06e+17 → 1.07e+17 | 5.01e+17 → 5.01e+17 | 2.14e+15 → 2.14e+15 | 4.93e+14 → 1.17e+14 |
| B, 222.4 | 961 → 83 | 399 → 32 | 2.38e+17 → 3.1e+17 | 3.19e+17 → 2.57e+17 | 2.43e+15 → 2.47e+15 | 1.22e+15 → 1.55e+14 |
| B, 222.5 | 87 → 81 | 33 → 30 | 0 → 7.52e+16 | 5.74e+17 → 5.11e+17 | 2.16e+15 → 2.15e+15 | 1.59e+14 → 1.46e+14 |
| B, 222.6 | 972 → 84 | 379 → 30 | 2.4e+17 → 3.51e+17 | 4.06e+17 → 3.09e+17 | 2.43e+15 → 2.45e+15 | 1.03e+15 → 1.53e+14 |
| C, 201.1 | 26 → 26 | 10 → 10 | 0 → 0 | 4.4e+17 → 4.4e+17 | 4.82e+14 → 4.82e+14 | 2.9e+13 → 2.9e+13 |
| C, 201.2 | 33 → 33 | 15 → 15 | 0 → 0 | 4.46e+17 → 4.46e+17 | 1.55e+14 → 1.55e+14 | 4.12e+13 → 4.12e+13 |
| C, 201.3 | 31 → 31 | 13 → 13 | 0 → 0 | 4.76e+17 → 4.76e+17 | 5.04e+14 → 5.04e+14 | 3.69e+13 → 3.69e+13 |
| C, 201.4 | 489 → 27 | 22 → 13 | 6.05e+15 → 6.05e+15 | 4.8e+17 → 4.8e+17 | 4.4e+14 → 4.4e+14 | 5e+13 → 3.77e+13 |
| C, 201.5 | 504 → 29 | 329 → 13 | 2.08e+17 → 2.08e+17 | 1.77e+17 → 1.77e+17 | 1.02e+15 → 1.02e+15 | 4.33e+14 → 3.47e+13 |

### Preregistered checks

| # | Check | Threshold | Result | Verdict |
|---|---|---|---|---|
| 1 | Steps reaching the pass limit | 0 of 15 | 11 → 0 | Pass |
| 2 | Internal transfer ÷ external supply | ≤ 2 every step | 10–32 | **Fail** |
| 3 | First-step final export and ocean overflow | Within 0.1% | Relative difference 0 in all windows | Pass |
| 4 | Normalized ledger residual | Order of 1.4e-12 | Largest 1.4e-12 on both arms | Pass |

**Check 2 failed because the threshold was wrong, and it stays recorded as a failure.** It assumed
surplus moves about two hops, and ignored that one volume is counted again at every hop of a real
cascade. The control steps C 201.1–201.3 had no capped cycle and read 10–15 on both arms, so the
threshold could not have passed even there.

A follow-up step on window B at 222.2 Myr classified what remains. Of 2.55e19 m³ internal transfer,
1.54e19 m³ moved between reciprocal pairs, and 99.6% of that was between bodies both at a real rim.
That step had 114 merges, and its 70 reciprocal pairs match its 69 passes, which points to a merge
ladder climbing about one rung per pass. Capped pairs carried only 4.7e15 m³. The largest directional count, 3, was a one-way pair receiving three different volumes.

The ladder ends on its own and does not reach the pass limit, so it is not this issue's defect. It
still counts as silt supply, so #693 should measure it.

### What moved

**The issue's invariant holds.** Every first step matches baseline exactly. In steps with large
cycles, internal transfer falls 84–96%, and no pair sends more than three times against up to 1,032
on `origin/main`. On the fixture, internal transfer equals the surplus with 1 or 400 padding leaves, where
`origin/main` gives 8 and 806 times it. Those are the unit tests in `test_lakes.py`.

**Later steps drift through silt, as expected, and by more than the preregistration implied.**
Silt falls by up to 12.5 times, which moves rims from step 2. Window B's control step at
222.5 Myr therefore did not stay still: it picked up a 7.5e16 m³ final export, and its ocean overflow
fell 11%. Ocean overflow at 222.4 and 222.6 fell 20% and 24%. Most of that volume shifted into final
export, since ocean overflow plus final export differs by only about 2% in both steps.

**Window A's collapse did not happen on the branch.** The anchor component is unchanged through the
221.8 Myr event and then diverges:

| Age (Myr) | Anchor volume, baseline (km³) | Anchor volume, branch (km³) |
|---:|---:|---:|
| 221.7 | 249,389 | 249,389 |
| 221.8 | 879,186 | 879,186 |
| 221.9 | 660,973 | 668,129 |
| 222.0 | 227,732 | 930,602 |
| 222.1 | 144,267 | 727,038 |

This is one trajectory, and it has not been traced to a rim. It is consistent with #691's finding
that silt made up about 100 m of the 125 m laid on this rim, but it does not show that. Event `ΔV`
in the steps #680 registered is unchanged: A 221.8, B 222.3 and C 201.5 each match to four
significant figures.

**Runtime.** Step time falls where the cycle was largest. No other replay or test from this task ran alongside them.

| Window | Five steps (s) | Peak RSS (GB) |
|---|---:|---:|
| A | 109 → 65 | 4.0 → 4.5 |
| B | 93 → 60 | 5.4 → 5.2 |
| C | 97 → 93 | 5.7 → 6.0 |

## Where final export goes

Each step's final export was routed again alone through that step's flow targets and transit
losses. Over all 15 steps on the branch, it totals 1.93e18 m³, 12.8% of external supply (11.1% on
baseline).

| Destination | Share of final export |
|---|---:|
| Ocean | 0 |
| Lost in transit | 0.46% |
| Settled on a flow cycle | 98.7% |
| Settled at a dead end | 0.81% |

About 45% of it settles on dry cells. The flow cycle is the same pair seen by the river router: each
spilling sink points across the rim at the other, so `route_downstream` treats those nodes as stuck
and keeps the water there.

**Settled land water goes nowhere.** `settled_m3` is read only to report ocean-node discharge, and
`overflow_settled` only for the ocean figure. The ledger books final export as leaving the lake
system, which is why it closes, but nothing receives it: not a lake, not the ocean, and not a later
step's supply. It is not a silent ledger error, because the ledger's boundary stops at the export.
It is still the largest open question in this water path. Principle 6 applies directly.

## Recommendation: export through the parent's outlet

Use the parent's outlet, and deliver through the hierarchy instead of the river router.

- The tier cap is a numerical guard (#553), and the saddle is real. Without the cap, this water would
  fill the deep basin, merge the pair and leave over the parent's rim. The parent's outlet is
  therefore the physical destination, and the leaf sink is only where the loop happened to stop.
- Routing alone cannot show the gain. 84% of the volume has a parent outlet in another basin, and
  the river router strands that too. Only the 16% with a river outlet could be tested, and 94% of
  that reaches the ocean (15% of all final export).
- So the transfer must go to the body that owns the parent's outlet target, or the ocean, as an
  ordinary hierarchy transfer. That body may have room, in which case the water stays in lakes.

This changes where water goes, moves ocean overflow and lake storage, and needs its own replay
comparison. It also touches how a cap relates to a real saddle, so it may need an RFC 0014 revision.
