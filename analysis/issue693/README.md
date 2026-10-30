# Silt supply counts water once: replay comparison

Replay evidence for [#693](https://github.com/adubey/mantle-bloom/issues/693), measured after
[#696](https://github.com/adubey/mantle-bloom/issues/696).

**Principle:** principle 6, account for material through every process. RFC 0014 specifies silt
from a body's step-averaged supply, deposited once per body.

## Preregistration

Committed before any replay of this branch on the post-#696 base. An earlier run on the
pre-#696 base is recorded under [Superseded run](#superseded-run), so this comparison does not
hide what prompted it.

**Change.** A final body's silt supply is its external water plus overflow that entered from
another final body. Transfers between its own children are dropped, closed loops between final
bodies cancel, and supply is bounded by world external supply. Water routing is not touched.

**Replays.** Windows A, B and C, five 100 kyr steps from each registered #680 checkpoint with
`ice_age_period_years = 0`, on a detached `origin/main` tree (`23ad7d7`) and on this branch, one
at a time. Anchors are `analysis/issue696/anchor-window-*.json`.

| Window | Checkpoint age (Myr) | SHA-256 |
|---|---:|---|
| A | 221.6 | `e4b0b4f09b87584aa7f77045f0f079f2309de0132049915a9aeb5589add27340` |
| B | 222.1 | `86fb5b7d2044846ac751906ff996184c7abe5454c93b39033bd82338e8749544` |
| C | 201.0 | `041b7430a6cf8867f2e2c76f18c3328e440add585ee4a2e4b3055e7766c2217e` |

**Primary metric.** The largest supply any final body's silt implies in a step, divided by world
external supply.

**Decision rule.**

1. Baseline: if this ratio is at or below 1 on all 15 steps, the excess is gone, and #693 reduces
   to the invariant test. If any step exceeds 1, the excess survives #696 and this change stands.
2. Branch: the ratio must be at or below 1 on all 15 steps.
3. World bound: report every step where a body's supply before the bound exceeds world supply.
   If the bound binds, that is where it changes a through-flow lake's silt.

**Controls.** Normalized ledger residual stays at the 1.4e-12 order. Each window's first step
must match the baseline in lake volume, spill level, internal transfer and overflow, because no
silt from this change has yet fed back.

**Secondary metrics.** Anchor component silt, spill level, event `ΔV`, world lake silt and
internal transfer.

**Expected direction and size.** The baseline ratio stays above 1, because #696 left internal
transfer at 10–32× external supply in merge ladders between bodies at real rims. World lake silt
falls by more than 10× in those steps, since #696's silt of about 1e14 m³ per step implies about
20× world supply in total. Trajectories diverge from step 2 by amounts not predicted here. On the
pre-#696 base, window A's event silt moved from 1.878 to 0.0035 mm/yr.

**Variability.** The replay is deterministic, so each arm runs once. A and B share a world, so
they are not independent, and steps 2–5 are a one-trajectory comparison, not a distribution.

## Superseded run

Before #696 merged, this change was compared with `d709f2b` on window A and on the first two
steps of B's branch arm, then stopped when #693 became blocked. The ratio was 33–211 on the
baseline and 0.24–0.35 on the branch. Window C was not run. These numbers are not acceptance
evidence, because #696 changed the baseline.
