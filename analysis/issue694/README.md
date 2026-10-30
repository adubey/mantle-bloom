# Body-level lake sublimation: causal probe

Probe for [#694](https://github.com/adubey/mantle-bloom/issues/694), under
[RFC 0014](../../rfc/0014-lake-water-volume-and-freezing.md). It tests whether charging lake
sublimation inside the lake solve, per body, changes the G8 events and the step-length dependence
of the current scheme. D2(c) charges sublimation to the ice reservoir, so a design change would
need an RFC revision. This report authorizes no implementation.

**Principles:** principle 3 (model fast processes in aggregate), because today's amount depends
on stored water rather than a rate over the step, and principle 4 (treat each node as an area),
because cells in one lake lose different amounts.

## Conclusions

1. **The G8 rule is met, so the mechanism is material.** Event-component `|ΔV|` changes by at least
   10% in A and C. B moves 0.1%.
2. **Current sublimation depends on the step, and Body's does not.** Going from 100 to 10 kyr
   multiplies Current's km³/yr by 6–11 and Body's by 0.86–1.24.
3. **Current removes most of the water in cold lakes after the solve.** It empties 56–79% of the
   solved area of large bodies, against 0.3–3.2% in Body, so cold lakes barely exist at step end.
4. **Window C's anchor flicker is gone in Body.** Its anchor component reverses on every step in
   Current and on no step in Body, on both schedules. This reading is post hoc.
5. **Body at 0.6 m/yr takes 16–25% of external supply,** and its lakes had not settled after five
   steps. The rate needs review before any revision adopts it.

**Recommendation: an RFC revision and a rate review of `LAKE_ICE_SUBLIMATION_M_PER_YR`.** See
[Recommendation](#recommendation). The issue's "modest effect" expectation does not hold.

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

## Results

Both arms start from `f4438e1c`, and the preregistration is `acae885`. All 15 runs completed one at
a time, and checkpoint SHA-256 values match the table above. Compact results are in
[`replay.json`](replay.json).

```bash
backend/.venv/bin/python bin/debug/measure_lake_sublimation_probe.py <checkpoint.mbworld> \
  --anchors analysis/issue696/anchor-window-a.json \
  --backend <arm-tree>/backend --myr 0.5 --step-years 100000 --out <scratch-output>
```

### Decision checks

| Check | Threshold | Result | Verdict |
|---|---|---|---|
| Accounting: residual | At most 0.1% of supply | Largest 3.4e-11, every step of every arm | Pass |
| Accounting: booked once | Solve's sublimation equals the ledger's | Equal; evaporation + sublimation equals the solve's loss to 2e-16 | Pass |
| G8 limb (a), `|ΔV|` | At least 10% in two windows | A and C | **Met** |
| G8 limb (b), anchor class | Mismatch at two consecutive 10 kyr outputs in two windows | A only (longest run 9; B 0; C 1) | Not met alone |
| Step-length limb | Current's ratio at least 2, Body's in 0.67–1.5, in two windows | All three | **Met** |
| Rate limb: flux | Body at least 10% of supply in two windows | All three (16–25%) | **Met** |
| Rate limb: persistence | Body 20 points below Current in two windows | Largest gap in window means 8 points | Not met |

### Event `ΔV`

| Window | Current (km³) | Body (km³) | Change in `|ΔV|` | Body at 0.25× (km³) | Change |
|---|---:|---:|---:|---:|---:|
| A, 221.8 Myr | +629,794 | −80,955 | −87.1% | +763,489 | +21.2% |
| B, 222.3 Myr | +131,625 | +131,723 | +0.1% | +131,691 | +0.1% |
| C, 201.2 Myr | −261,587 | +1,397,267 | +434% | +2,285,992 | +774% |

Current reproduces the #693 branch values. Liquid-plus-ice changes these percentages by at most
0.2 points.

**These are not one event at two sizes in A and C.** In Body, A's anchor is dry at the event and
its nearest component lies 85 km away, so the comparison crosses different water. C's anchor sits in
a 4.8 M km² cold lake rather than Current's 163,000 km² remnant. The G8 rule is met as written, but
it is better read as "the event does not exist in this form" than as a change in magnitude. B's
target lake holds almost no ice, and its event is untouched.

### Sublimation and step length

| Window | Current, 100 kyr (km³/yr) | Current, 10 kyr | Ratio | Body, 100 kyr | Body, 10 kyr | Ratio | Body ÷ Current, 100 kyr |
|---|---:|---:|---:|---:|---:|---:|---:|
| A | 160 | 991 | 6.2 | 2,733 | 3,096 | 1.13 | 17 |
| B | 276 | 1,672 | 6.1 | 4,102 | 3,519 | 0.86 | 15 |
| C | 63 | 662 | 10.5 | 1,658 | 2,048 | 1.24 | 26 |

| Window | Share of external supply (%): Current, 100 kyr | Current, 10 kyr | Body, 100 kyr | Body, 10 kyr | Body at 0.25× |
|---|---:|---:|---:|---:|---:|
| A | 1.2 | 7.3 | 16.4 | 17.2 | 5.0 |
| B | 2.4 | 11.4 | 24.5 | 19.2 | 7.8 |
| C | 1.0 | 8.4 | 22.1 | 22.5 | 5.4 |

Current's flux scales nearly as 1/step, which is what a stock limit predicts. Even at 10 kyr it is
2–3× below Body, so the schedule does not converge on Body's figure within this range. Body never
exceeds the continuous figure.

### Lakes

Means over the five 100 kyr steps, Current → Body.

| Window | World liquid (M km³) | Lake area (M km²) | Cold liquid (km³) | Lake evaporation (10⁸ km³) | Overflow out (10⁸ km³) |
|---|---:|---:|---:|---:|---:|
| A | 2.5 → 9.4 | 5.4 → 14.8 | 39,600 → 6.4 M | 5.6 → 8.1 | 6.8 → 4.9 |
| B | 2.3 → 18.7 | 4.9 → 18.2 | 283 → 15.3 M | 5.1 → 7.5 | 6.1 → 4.1 |
| C | 0.5 → 6.9 | 1.6 → 6.4 | 59 → 6.0 M | 1.6 → 2.4 | 4.7 → 3.2 |

Sublimating more leaves more water. Current's solve produces this water and the phase step then
removes it, so cold lakes keep only a sliver of it at step end (B: 22 M km³ solved, about 300 km³
left). RFC 0014's "a cold lake keeps its water" does not hold at step end under Current.

**Within-body state, large bodies.**

| Window | Solved area emptied (%) | Largest spread (m) | Bodies mixed above 1 m | Largest median spread (m) |
|---|---:|---:|---:|---:|
| A | 55.9–71.1 → 0.6–3.2 | 3,211 → 44 | 11–12 → 13–28 of 32–46 | 0 → 4.5 |
| B | 72.8–73.9 → 0.5–1.0 | 3,105 → 44 | 9–10 → 11–27 of 29–40 | 0 → 6.6 |
| C | 68.8–78.5 → 0.3–1.1 | 2,622 → 37 | 4–7 → 4–7 of 18–27 | 0 → 0 |

Emptying and the extreme spreads fall by about two orders of magnitude. The count of bodies mixed
above 1 m does not fall, and it rises in A and B. Current's zero median is partly an artifact, since
emptied cells drop out of the spread. Body's small spreads come from ice exchange on partly flooded
cells, where a little depth moves the level a long way: the rate-zero control already shows them.

**Cold-lake persistence**, as window means, is 100% in Current and 92–98% in Body at 100 kyr, and
77–100% against 96–99% at 10 kyr. Current's 100 kyr figures rest on very little wet cold area
(0.04–0.2 M km² against 4.5–11.6 M km²), so they say little.

**Ice-cap checks.** Thickness stays at or below 2 m with no bound violations, and deep bedfast area
is 0.00% in every arm, within G5's 1%. Inventory rises from at most 1,507 km³ to at most 33,311 km³.
Shallow bedfast area rises from at most 0.2% to 18–27% of water-holding area at 100 kyr and 29–43%
at 10 kyr, because Body holds shallow cold fringes that Current empties. G5 reports this without a
limit, but it is well above G5's 0.5–11%, so a revision should look at it.

**Runtime** is unchanged: 63–111 s per 100 kyr window and 9–17 minutes per sidecar, with peak RSS
at most 7.1 GB.

### Predictions

| Prediction | Result |
|---|---|
| 1. Current at least 3× at 10 kyr; Body within 0.67–1.5 | Held: 6.1–10.5 and 0.86–1.24 |
| 2. Body 10–50× Current, below continuous, 10–50% of supply | Held: 15–26×, 16–25% |
| 3. Emptying largely disappears; spread falls | Partly held: emptying and largest spread fell, but mixed-body count did not |
| 4. World liquid rises and overflow falls | Held |
| 5. `|ΔV|` changes at least 10% in two windows | Held, though not as a change in event size |

### Post hoc: anchor reversals

Not preregistered, and read from the anchor component only, so this is not RFC 0014's Flicker
metric. It counts sign reversals where both changes exceed 50% of the prior volume.

| Window | Current, 100 kyr | Body, 100 kyr | Current, 10 kyr | Body, 10 kyr |
|---|---:|---:|---:|---:|
| A | 0 of 3 | 1 of 3 | 0 of 38 | 1 of 38 |
| B | 0 of 3 | 0 of 3 | 5 of 38 | 5 of 38 |
| C | 3 of 3 | 0 of 3 | 38 of 38 | 0 of 38 |

In C, Current's anchor alternates wet and dry at every 10 kyr output, and it stays wet throughout
in Body. Limb (b) does not count this, because a mismatch at every other output is never
consecutive. The mechanism behind the alternation has not been traced. A and B show no such
effect, so this is a lead for one window, not an explanation of the flicker.

### Limits

- Five 100 kyr steps on two worlds, and A and B share one. Every arm starts from a checkpoint shaped
  by Current, so step 1 and the event at step 2 include Body's spin-up.
- **Body has not settled.** At full rate, world liquid falls from 29.9 to 12.1 M km³ over B's five
  steps and from 12.4 to 7.9 M km³ over A's last three. At 0.25× it rises in A and B, and falls
  from 7.5 to 6.3 M km³ over C's last four steps.
- The windows were chosen for lake events, not for cold lakes.
- The patch is a probe. It returns ice on a dry cell to liquid rather than to vapour.
- The zero-rate control differs from Current by 0.75% in ice inventory.
- 0.6 m/yr and 0.15 m/yr are the only rates tested. No new source was reviewed for the rate.

## Recommendation

Applying the preregistered rule: the accounting gate passes, step 2 selects an RFC revision on both
limbs, and step 3 adds a rate review.

1. **Revise RFC 0014 D2(c).** Charge sublimation in the body's budget. The RFC already points both
   ways: D1's budget defines `Ē` as "including D2's sublimation", and the planned breach-gate demand
   uses "the solve's evaporation and sublimation", while D2(c) keeps "existing loss transfers ...
   as currently modeled". Current's step-end stripping also conflicts with D2(b)'s cold lake that
   "still holds liquid water". This restores principle 3, and principle 4 for the amount each cell
   loses.
2. **Review `LAKE_ICE_SUBLIMATION_M_PER_YR` first or alongside it.** At 0.6 m/yr, Body sends 16–25%
   of external supply through ice and lakes keep shrinking. The rate matches the roughly 60 cm/yr
   ablation the RFC cites for Lake Untersee, a perennially ice-covered Antarctic lake, and it
   applies to every ice-covered wet cell. Under
   Current that choice was hidden by the stock limit.
3. **Before adopting either,** run Body long enough to settle, from a checkpoint it has spun up,
   and report the Flicker metric, shallow bedfast area and ocean delivery. Validation needs the
   fresh held-out set G8 requires.

Whether to open that revision is the owner's decision. G8's materiality also does not make
sublimation the unique cause of the flicker.

## Save compatibility

Unaffected. Nothing in `backend/` changes.
