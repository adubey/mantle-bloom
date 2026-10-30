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

## Results

Baseline `23ad7d78`, branch `254729b7` (the code in this PR, before this write-up). Compact results
are in [`replay.json`](replay.json).

```bash
backend/.venv/bin/python bin/debug/measure_lake_connector_budget.py <checkpoint.mbworld> \
  --anchors analysis/issue696/anchor-window-a.json \
  --myr 0.5 --step-years 100000 --out <scratch-output>
```

**The excess survived #696, so this change stands.** No body on the branch exceeds the bound:

| Window, age (Myr) | Max body ÷ world supply | World lake silt (m³) | Anchor silt (mm/yr) | Anchor `ΔV` (km³) | Spill level (m) | Normalized residual |
|---|---:|---:|---:|---:|---:|---:|
| A, 221.7 | 13.0 → 0.31 | 1.09e14 → 7.07e12 | 0.0032 → 0.0032 | n/a | not spilling | 5.6e-13 → 5.6e-13 |
| **A, 221.8** | 15.6 → 0.35 | 8.14e13 → 4.77e12 | 0.0077 → 0.0024 | 629,794 → 629,794 | −805.6 → −805.6 | 5.7e-14 → 5.8e-14 |
| A, 221.9 | 21.8 → 0.31 | 1.27e14 → 8.37e12 | 0.0027 → 0.0027 | −211,054 → −211,027 | not spilling | 2.0e-13 → 2.0e-13 |
| A, 222.0 | 22.2 → 0.29 | 1.19e14 → 4.93e12 | 0.0080 → 0.0023 | 262,473 → 262,479 | −676.4 → −676.4 | 3.0e-13 → 2.9e-13 |
| A, 222.1 | 22.6 → 0.28 | 1.33e14 → 5.20e12 | 0.0028 → 0.0027 | −203,569 → −246,093 | not spilling | 6.1e-13 → 2.8e-14 |
| B, 222.2 | 23.1 → 0.25 | 1.08e14 → 4.47e12 | 0.0061 → 0.0061 | n/a | not spilling | 7.1e-14 → 7.1e-14 |
| **B, 222.3** | 23.8 → 0.30 | 1.17e14 → 4.96e12 | 0.0086 → 0.0052 | 131,625 → 131,625 | 664.8 → 664.8 | 1.2e-13 → 1.2e-13 |
| B, 222.4 | 30.2 → 0.32 | 1.55e14 → 4.90e12 | 0.0056 → 0.0056 | 4,594 → 4,748 | 668.5 → 668.5 | 9.2e-14 → 9.0e-14 |
| B, 222.5 | 29.9 → 0.33 | 1.46e14 → 4.87e12 | 0.0059 → 0.0059 | −302 → −303 | 670.1 → 670.1 | 1.4e-13 → 1.3e-13 |
| B, 222.6 | 29.1 → 0.25 | 1.53e14 → 5.27e12 | 0.0069 → 0.0068 | −88,261 → −88,527 | 671.7 → 671.8 | 7.5e-14 → 7.5e-14 |
| C, 201.1 | 10.5 → 0.33 | 2.90e13 → 2.99e12 | 0.0656 → 0.0021 | n/a | −1,471.2 → −1,471.2 | 1.4e-12 → 1.4e-12 |
| C, 201.2 | 15.4 → 0.38 | 4.12e13 → 2.64e12 | 0.0785 → 0.0019 | −261,592 → −261,587 | −1,469.6 → −1,471.2 | 3.6e-15 → 3.8e-15 |
| C, 201.3 | 12.7 → 0.33 | 3.69e13 → 2.96e12 | 0.0683 → 0.0018 | 369,382 → 369,373 | −1,454.3 → −1,457.7 | 1.7e-13 → 1.7e-13 |
| C, 201.4 | 12.9 → 0.38 | 3.77e13 → 2.80e12 | 0.0711 → 0.0021 | −369,616 → −369,622 | −1,465.4 → −1,470.5 | 3.8e-14 → 3.8e-14 |
| **C, 201.5** | 12.7 → 0.36 | 3.47e13 → 2.71e12 | 0.0792 → 0.0023 | 262,731 → 262,662 | −1,450.0 → −1,456.8 | 8.1e-14 → 8.1e-14 |

Each cell is baseline → branch. Bold rows are the #680 registered events. Spill level is the
narrowest listed lake around the anchor that is spilling; C's anchor component spans cells at
several levels, so its row is that lake's level, not the component's single surface.

### Preregistered checks

| Check | Threshold | Result | Verdict |
|---|---|---|---|
| Baseline ratio decides the scope | Any step above 1 means the change stands | 10.5–30.2, all 15 steps | Change stands |
| Branch ratio | At most 1 on all 15 steps | 0.25–0.38 | Pass |
| World bound | Report each step it binds | Largest supply before the bound equals the table's ratio on every step: 0 of 15 | Never bound |
| Ledger residual | 1.4e-12 order | Largest 1.4e-12 on both arms | Pass |
| First step of each window | Matches baseline | Lake volume, surface, rims, internal transfer and overflow are exactly equal | Pass |

The world-bound row settles the concern that the bound could reduce a through-flow lake's silt
where released storage exceeds world supply: that did not happen in these windows.

### What moved

**#696 removed the spike but not the count.** Merge ladders kept every body's supply 10–30× world
supply, on ordinary steps as well as events. World lake silt falls 10–32×, more than the expected
10×. The branch's total of about 5e12 m³ per step is close to `coefficient × external supply`
(about 4.5e12 m³), as it should be if water is counted about once.

**Water is moved only through silt feedback.** Internal transfer, which this change does not read,
is unchanged on every first step. From step 2 small differences follow from terrain: the
registered events keep their `ΔV` exactly in A and B and within 0.03% in C. Window A's 222.1 Myr
step is the largest difference, a `ΔV` of −246,093 against −203,569 km³.

**Spill levels are mostly unchanged in A and B.** Window C's listed lake sits up to 6.8 m lower,
because excess silt no longer builds its rim.

**Calibration consequence, not changed here.** `SILT_ACCUMULATION_COEFFICIENT` was set so total
silt matched the older rule on a reference replay (`analysis/issue363/`). That total is now
10–32× lower on these windows, so lakes silt in correspondingly slower. Recalibrating is out of
scope for #693, and is a decision for later work, not an automatic follow-up.

**Runtime.** Each window's five steps ran alone.

| Window | Time (s) | Peak RSS (GB) |
|---|---:|---:|
| A | 75 → 72 | 4.8 → 3.9 |
| B | 67 → 68 | 5.5 → 6.4 |
| C | 107 → 106 | 5.6 → 6.1 |

### Limits

Five steps on two worlds, and A and B share one. The windows were chosen for lake events, so they
do not sample all lake states. The bound on released storage was never exercised here, so it is
covered only by the unit test.

## Save compatibility

Unaffected. The change sits in the solver's per-step accounting, which is not pickled.
