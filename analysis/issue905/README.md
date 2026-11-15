# Issue #905: lake slice pass-bound failure at a saddle

A 0.1 Myr-step animation of seed 832424611 stopped at 216.1 Myr with
`RuntimeError: Lake slice failed to settle transfers and topology within its pass bound`.

## Mechanism

The pass loop in `_settle_prognostic_slice` exits only after a pass with no change to the
active body set and no deferred transfer. Its bound, `2 * len(leaves) + 3`, runs to thousands
on a real world, so hitting it means a cycle, not slow convergence.

The cycle is a period-2 merge/split at a saddle. `_integrate_body_volume` stopped Newton when
`|volume + loss - available| <= 1e-10 * available`. Over a long slice, `available` is mostly
the slice's losses: here about 200 times the storage, and in a mature world up to about 5,000
times. The solved storage was therefore loose by about 1e-8 of storage, against the 1e-10
saddle-storage test in `_active_bodies` (#768). On the `_MERGE_ELEVATION` pair with 1 m
sub-cell relief, a single 1,000 yr slice and supply at 1 - 1e-10 of the saddle budget:

| Pass state | Storage vs. saddle | `_active_bodies` verdict |
|---|---|---|
| children: one at its cap, one at the solver's landing | -5.1e-11 | merged |
| merged body, same water | -1.18e-10 | split |

The two states repeat bit for bit until the bound.

## Scan

There are 2,232 near-saddle cases: supply at 0 and ±1e-15 … ±1e-8 of the exact saddle budget
(31 values), sub-cell relief 0, 1 and 3 m, four slice/evaporation combinations, three
supply splits, and 1 or 4 substeps.

| Build | Failures |
|---|---|
| `origin/main` 849e30c1 | 3 (all with relief 1 m and one substep) |
| Branch | 0 |

## Replay

The run is labelled `seed832424611-216.1myr-torn`: the live world downloaded through
`GET /worlds/{id}/save` after the failed step, so it is partly stepped. It did not reproduce
the failure, which needs a body within solver slop of its saddle. It was stepped 20 times at
100 kyr with `bin/debug/capture_lake_solve_failure.py` in both trees:

| Metric (20 steps, 216.1 → 218.1 Myr) | `origin/main` | Branch |
|---|---|---|
| Mean lake ledger residual per step | 1.7e5 m³ | 4.4e2 m³ |
| Lake volume, steps 0–9 | 7.74e15 m³ at step 9 | identical to 10 digits |
| Lake volume, step 19 | 7.633e15 m³ | 7.655e15 m³ (+0.28%) |
| Mean step time | 7.06 s | 6.96 s |
| Level-solve budget evaluations, 3 steps | 46,509 | 42,244 |
| Repeat-state backstop firings | n/a | 0 |

Volumes match to 10 digits until step 10, then drift by at most 0.35% as different lake
depths feed back through climate and routing. In a one-step histogram, most level solves
converge in 2–6 iterations and none takes more than 17.

The convergence test is on the storage error (level error times wet area), not on the
budget residual. Losses dominate the slope of the budget against level over a long slice, so
the residual can't reach a storage-scale tolerance before the level reaches float precision.
A first version that tested the residual against storage took about 19 budget evaluations a
solve and ran about 15% slower.
