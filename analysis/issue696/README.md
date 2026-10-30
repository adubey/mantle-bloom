# Issue #696: capped lake siblings recirculating surplus

Status: preregistered. Results follow in this file after the replays.

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
