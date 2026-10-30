# Body-level lake sublimation: causal probe

Probe for [#694](https://github.com/adubey/mantle-bloom/issues/694), under
[RFC 0014](../../rfc/0014-lake-water-volume-and-freezing.md). It tests whether charging lake
sublimation inside the lake solve, per body, changes the G8 events and the step-length dependence
of the current scheme. D2(c) charges sublimation to the ice reservoir, so a design change would
need an RFC revision. This report authorizes no implementation.

**Principles:** principle 3 (model fast processes in aggregate), because today's amount depends
on stored water rather than a rate over the step, and principle 4 (treat each node as an area),
because cells in one lake lose different amounts.

## Preregistration

Committed before any run with a non-zero treatment. The only arm run earlier is the zero-rate
control below, which was needed to validate the script and patch.

**Source.** Both arms start from `origin/main` at `f4438e1c`, after #696 and #693.

**Arms.**

| Arm | Tree | Sublimation |
|---|---|---|
| Current | Detached `f4438e1c` | Per cell in `advance_lake_ice`, after the solve, limited by stored water |
| Body | Scratch worktree with [`body-sublimation.patch`](body-sublimation.patch) | `0.6 m/yr × ice cover × wet area`, inside each body's budget in the solve, and in the breach gate's demand |
| Body at 0.25× | Same tree, `--sublimation-scale 0.25` | As Body at 0.15 m/yr, the low end of the G5 range. 100 kyr only |

The patch stays out of `backend/`. In the Body arm the solve debits liquid, the ice keeps its
thickness trajectory, and ice and liquid exchange only the net change. Sublimation is booked once
at the solve's figure, and the ledger's evaporation excludes it. Ice on a cell the solve left dry
returns to liquid rather than to vapour. That simplification affects a small inventory and is a
limit of the probe.

**Windows.** The registered #680 checkpoints, `ice_age_period_years = 0`, with anchors
`analysis/issue696/anchor-window-*.json`.

| Window | Checkpoint age (Myr) | SHA-256 | Event step |
|---|---:|---|---:|
| A | 221.6 | `e4b0b4f09b87584aa7f77045f0f079f2309de0132049915a9aeb5589add27340` | 2 (221.8 Myr) |
| B | 222.1 | `86fb5b7d2044846ac751906ff996184c7abe5454c93b39033bd82338e8749544` | 2 (222.3 Myr) |
| C | 201.0 | `041b7430a6cf8867f2e2c76f18c3328e440add585ee4a2e4b3055e7766c2217e` | 2 (201.2 Myr) |

**Schedules.** Five 100 kyr steps per window, for the event and for sublimation over several
steps. A separate 0.4 Myr sidecar at 10 kyr steps, Current and Body only, for step-length
sensitivity and anchor connectivity. Runs go one at a time.

**Primary metric.** Relative change in event-component `|ΔV|`, Body against Current, at the event
step of the 100 kyr schedule. `ΔV` is the liquid volume change of the wet component nearest the
anchor, read as in #691 and #693. Liquid-plus-ice is reported beside it.

**Materiality (owner's G8 rule, #680).** Material if, in at least two of the three windows, either
`|ΔV|` changes by at least 10%, or the anchor's wet/dry class differs between arms at two
consecutive 10 kyr outputs. The 20% figure is descriptive. The 10% rule does not apply to the
10 kyr `ΔV`.

**Step-length metric.** Mean sublimation in km³/yr over the 10 kyr sidecar, divided by the mean
over the 100 kyr schedule, per arm and window.

**Accounting gate.** Normalized residual at most 0.1% of external supply on every step of every
arm, and the solve's sublimation equal to the ledger's in the Body arm. Failure stops
interpretation.

**Recommendation rule, applied in order.**

1. Accounting gate fails: no recommendation, return to accounting.
2. **RFC revision** if the G8 rule is met, or if in at least two windows Current's step-length
   ratio is at least 2 while Body's lies within 0.67–1.5.
3. **Rate review**, as well as or instead of the above, if Body's mean sublimation is at least 10%
   of external supply in at least two windows, or Body's cold-lake persistence falls at least 20
   percentage points below Current's in at least two windows.
4. **No change** if none of these applies.

**Secondary metrics.** Within-body level spread and area emptied after the solve (bodies above
10⁴ km², the Flicker size), global lake area, world liquid, overflow, cold-lake state (cells with
ice cover at least 0.5: area, liquid, component count, and step-to-step persistence of wet area),
and the ice-cap checks (maximum thickness at most 2 m, inventory, bedfast area split at 2 m of
water; G5's limit is at most 1% bedfast in deep water).

**Predictions.**

1. Current's sublimation in km³/yr is at least 3× higher at 10 kyr than at 100 kyr, because it
   is limited by stored water. Body's stays within 0.67–1.5.
2. Body sublimates 10–50× Current at 100 kyr, stays below the continuous figure, and takes
   10–50% of external supply.
3. Current emptying of solved cells largely disappears in Body, and within-body spread falls.
4. World liquid at step end rises in Body, despite the larger flux, because bodies keep their
   solved water. Overflow falls.
5. `|ΔV|` changes by at least 10% in at least two windows. This is the least certain prediction.
   The events are large merge steps, which may be insensitive to a loss of 0.6 m/yr.

**Variability.** The replay is deterministic, so each arm runs once. A and B share a world, so
they are not independent. Trajectories diverge after step 1. #687 showed partition-sensitive
topology, so each 10 kyr comparison is within its schedule.

## Zero-rate control

One 100 kyr step of window A at `--sublimation-scale 0` in both trees, run before registering. The
two reservoir updates should agree when no sublimation acts. They differ by 0.003% in world
liquid (4,504,995 → 4,504,878 km³) and 0.75% in ice inventory (15,510 → 15,627 km³), with the
anchor component's liquid and ice unchanged and both residuals at 5.6e-13. This validates the patch
to within the update's form.
