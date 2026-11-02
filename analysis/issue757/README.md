# Lake-ice sublimation inside the lake solve: replay comparison

Replay comparison for [#757](https://github.com/adubey/mantle-bloom/issues/757). It measures what
moving lake-ice sublimation into the lake solve's water balance does to #754's windows D and E
and to one other mature world.

**Principles:** 3, because the old after-solve loss stopped growing with step length, and 8,
because sublimation acted outside the solve's balance, so evaporation and sublimation never
shared one equilibrium against supply.

**The change.** Sublimation, `0.6 m/yr × (1 − f_free)`, is a loss in the solve's implicit
balance beside open-water evaporation, and the breach gate's rim-level demand includes it. The
solve opens with liquid and ice together, because floating ice displaces its own mass of water.
`advance_lake_ice` then splits the solved column between the capped ice layer and liquid, booking
only the net freeze or melt. Before, the solve ran without sublimation and `advance_lake_ice`
charged it afterwards. Below the 2 m cap that loss saturated at about 42 m per step. The capped
and `balanced` branches stripped the whole column.

## Preregistration

Committed before any replay.

**Arms.**

| Arm | Tree | Commit |
|---|---|---|
| main | Detached `origin/main` | `e8ef1e83` |
| branch | `claude/issue-757-sublimation-in-solver` | `0c5aac7f` |

**Windows.** 100 kyr steps with the ice-age cycle off, as in #754. D and E are #754's
registered checkpoints and anchors, and S is #754's registered seed-317739385 checkpoint. S is
#754's coldest case: about a fifth of its anchor catchment's precipitation falls frozen. Its
anchor lake, window C's in #694 and #696, swings by about 370,000 km³ with no routed supply.

| Window | Checkpoint | Age (Myr) | SHA-256 | Anchor | Horizon |
|---|---|---:|---|---|---|
| D | `seed910211954-225.1myr` | 225.1 | `6a8d6f90959b3cc4235a48ebbd6ddb0373041f79db131d9ac9f370dd6806bff7` | [`target-anchor-window-d.json`](../rfc0016/target-anchor-window-d.json) | 0.4 Myr, 4 steps |
| E | `seed910211954-224myr` | 224.0 | `4045b915a2d3d9e534a84089bba7277bbd488690e51d70f4e928f58c8cf674b1` | [`target-anchor-window-e.json`](../rfc0016/target-anchor-window-e.json) | 1.2 Myr, 12 steps |
| S | `seed317739385-201.6myr` | 201.6 | `ad8d34d60ac89db97efdfd722bc3db0287db22d8917288871502fe9a7747f94e` | [`anchor-seed317739385-201.6.json`](anchor-seed317739385-201.6.json), #754's rescan anchor | 1.2 Myr, 12 steps |

**Command.** One arm at a time, each after a `compute-budget` check (`--run sim --gb 5`; #754's
probes of D and E peaked at 3.3–4.7 GB):

```bash
backend/.venv/bin/python bin/debug/replay_lake_water_swing.py <checkpoint.mbworld> \
  --anchors <anchor.json> --myr <0.4 | 1.2> --step-years 100000 \
  --backend <arm tree>/backend --out <scratch>/issue757-<window>-<arm>.json
```

**Metrics.** All are read as #754 read them. Whole-world lake volume is liquid `lake_depth` times
cell area, with the step change taken between consecutive step outputs. The target is the wet
component (`lake_depth` above 1 m) under the anchor.

- **Primary:** E's lag-1 autocorrelation of the step change in whole-world lake volume, over its
  12 steps. About −1 is period-2 alternation. Main is expected to repeat #754's −1.000.
- **Secondary:** E's mean absolute step change from step 3, as a share of main's, and D's target
  ratio across the 225.3 Myr drop, step 1 against step 2. Main is expected to repeat #754's
  2,622,950 km³ and about 830×. S gets the same autocorrelation and step change.
- **Reported, no threshold:** whole-world lake liquid, ice and total per step; sea level per step
  and branch minus main at each window's end; the five largest lakes' areas, volumes, levels,
  floors and whether the sea-tier guard (3,500 m) holds them; the target's area and level; booked
  evaporation and sublimation; runtime and peak RSS. #754's knockout turned E's basin into a
  13.7 million km² sea at the 3,500 m guard and lowered sea level by about 85 m. Those quantities
  are reported explicitly.

**Predictions and thresholds.**

| # | Prediction | Supported if |
|---|---|---|
| 1 | E's period-2 alternation ends | Branch E lag-1 autocorrelation above −0.5 |
| 2 | E's step changes shrink | Branch E mean absolute step change from step 3 below 20% of main's |
| 3 | D's collapse goes | Branch D target ratio across the 225.3 Myr drop below 2× |
| 4 | Cold lakes keep their water, so lake volume rises and sea level falls | Branch whole-world lake volume at E's end at least 2× main's, and branch sea level lower than main's there |
| 5 | Sublimation booked rises by an order of magnitude | Branch sublimation per E step at least 10× main's, as #694's body arm found (15–26×) |

These follow from #754's K-sub knockout and #694's body-level arm. Both held more lake water and
lowered sea level. Rough sizes: E's basin holds a large inland sea, several million km² and up to
the 13.7 million km² of K-sub, smaller than K-sub's because sublimation now removes 0.6 m/yr over
its ice-covered share. Sea level falls by tens of metres, at most the 85 m of K-sub. S: its anchor
lake holds more water and stops reversing, as window C's lake did in #694's body arm.

**Validity.** Interpretation stops if any of these fails:

- **Control:** main must reproduce #754's controls on D and E. D's target is 2,293,616, 2,768,
  268,996 and 788,971 km³ over its four steps, and E's whole-world lake volume runs 1,449,874 to
  4,119,264 km³, both within 0.01%. The commits since #754's `0bda513d` change new-world
  generation, the viewer and validation counters, not stepping a save.
- **Accounting:** in every step of every arm, normalized lake and land residuals are at most
  1e-9.

**Reading a miss.** If the branch keeps E's autocorrelation at or below −0.5, the level lag or
another mechanism still alternates at world scale, and the issue's fix is incomplete there; the
unit tests' lagged case doesn't cover world-scale lag. A large lake that holds the 3,500 m guard
is the guard's work, not the basin's rim, and is reported as such.

**Limits.** D and E share seed 910211954, and S is a second seed. Replays are deterministic, so
each arm runs once. #758, which classifies enclosed ocean-floor basins as ocean, changes which
basin D and E measure; whichever of #757 and #758 lands second reruns these numbers.

## Results

All arms ran on 2026-11-02, one at a time, each after a `compute-budget` check. Another agent's
replay ran alongside most arms, so runtimes compare between arms but aren't clean. Compact
results are in [`replay-comparison.json`](replay-comparison.json), and the post hoc pass trace is
in [`posthoc-e-branch-passes.json`](posthoc-e-branch-passes.json).

**Runs.** Main ran once, on `e8ef1e83`. The branch ran three times; the tables report the last:

1. **Discarded:** all three branch arms failed before stepping. The saves refused to load,
   because the branch had removed `HydrologyFields.lake_ice_trace` and v5 saves carry it.
   `693cb2be` drops that retired field on load (`npz_save._RETIRED_FIELDS`) and changes no
   model code.
2. **Superseded:** the arms reran on `693cb2be`, with the model of the preregistered `0c5aac7f`.
3. **Reported:** review found that the readers of last step's lake before the solve still read
   liquid alone, while the solve opened with liquid and ice. `e5c697d2` gives them the same
   column, and the arms reran on it. D's and S's whole-world lake volumes moved by at most
   0.35% and 0.08% from run 2. E's drain moved slightly, and its refill came a step earlier and
   larger: 28.3 M km³ at step 10 instead of 13.0 M at step 11. That doubles E's step-change
   metric, so the timing of a refill is sensitive to small differences. All outputs are kept.

**Validity: both gates pass.** Main reproduces #754's D and E controls exactly: whole-world lake
volume, target volume and sea level differ by zero at every step. The worst normalized residuals
in any step of any reported arm are 5.8e-12 for lakes and 1.3e-12 for land, in E's branch arm.

### Decision checks

| # | Check | Threshold | Main | Branch | Verdict |
|---|---|---|---:|---:|---|
| 1 | E: lag-1 autocorrelation of step change | Above −0.5 | −1.000 | −0.286 | Supported, see below |
| 2 | E: mean absolute step change from step 3 | Below 20% of main's | 2,622,950 km³ | 6,686,041 km³ (255%) | **Not supported** |
| 3 | D: target ratio across the 225.3 Myr drop | Below 2× | 829× | 1.48× | Supported |
| 4 | E's end: whole-world lake volume; sea level | At least 2× main's; lower | 4.12 M km³; −748.1 m | 14.30 M km³ (3.5×); −773.2 m | Supported |
| 5 | E: sublimation per step | At least 10× main's | median 26 M km³ | median 345 M km³; median per-step ratio 15× | Supported |

The period-2 alternation is gone, but the swings aren't smaller. They are now a slower cycle of
fill, drain and refill, traced [below](#what-drives-the-cycle-post-hoc). Check 1 passes on
11 changes dominated by two fills, so read it as the end of the every-step flip, not as a stable
statistic. Prediction 2 was modelled on #754's K-sub knockout, where step changes fell to 15% of
control's. Without sublimation, K-sub's basin had its equilibrium at the 3,500 m guard and stayed
pinned there. With sublimation charged in the solve, the basin's equilibrium lies below the
guard, and nothing pins its level.

### By window

| Measure | D main | D branch | E main | E branch | S main | S branch |
|---|---:|---:|---:|---:|---:|---:|
| Lag-1 autocorrelation of step change | −1.000 | +1.000* | −1.000 | −0.286 | −0.888 | +0.626 |
| Mean absolute step change from step 3 (km³) | 186,542 | 3,825,385 | 2,622,950 | 6,686,041 | 257,591 | 436,558 |
| Whole-world lake volume, window mean (M km³) | 2.52 | 18.62 | 2.78 | 14.71 | 0.93 | 6.45 |
| Whole-world lake volume, end (M km³) | 2.32 | 10.02 | 4.12 | 14.30 | 0.52 | 6.34 |
| Lake ice, end (km³) | 1,257 | 20,197 | 1,468 | 25,229 | 591 | 9,861 |
| Lake area, end (M km²) | 4.09 | 13.83 | 4.44 | 16.57 | 1.65 | 6.59 |
| Sea level, end (m) | −744.0 | −765.1 | −748.1 | −773.2 | −1110.0 | −1120.7 |
| Largest branch-below-main sea-level gap (m) | | 78.3 | | 77.2 | | 18.8 |
| Target at end: area (M km²), volume (M km³) | 0.96, 0.79 | 9.53, 8.01 | 1.16, 2.29 | 11.09, 11.65 | 0.10, 0.02 | 4.57, 5.85 |
| Three largest lakes at end (M km²) | 1.51, 0.96, 0.53 | 9.53, 1.51, 0.68 | 1.50, 1.16, 0.68 | 11.09, 1.50, 1.43 | 0.64, 0.36, 0.18 | 4.57, 0.66, 0.33 |
| Sublimation per step, median (M km³) | 7.7 | 429 | 26.0 | 345 | 7.3 | 183 |
| Median step runtime (s) | 7.9 | 6.7 | 8.2 | 7.0 | 10.9 | 10.1 |
| Peak RSS (GB) | 2.9 | 4.3 | 4.4 | 4.5 | 7.1 | 6.4 |

\*D has four steps, so its autocorrelation rests on two pairs of changes and means little.

At the end of every branch run the target is the world's largest lake, holding several times
main's water. Mean sublimation per step is 17–34× main's (D 34×, E 17×, S 26×), comparable to
#694's 15–26× for its body-level arm. Every fourth step takes 20–46 s on both arms. S peaked near
7 GB, above the 5 GB the preregistration planned for. The replay changes the breach gate's demand
and the solve's loss together, so it doesn't separate their effects.

**E's basin and the 3,500 m guard.** At step 2 the branch fills E's basin to the sea-tier guard,
3,500 m above its floor: 27.87 M km³ over 13.55 M km², with sea level at −826.1 m, 77 m below
main's. #754's K-sub knockout had 27.88 M km³ over 13.54 M km² at −828.1 m at the same step. The
basin doesn't stay there. Over steps 3–9 it drains to 3.10 M km³ over 7.40 M km² and splits at
a saddle. At step 10 it refills to 28.32 M km³ over 13.70 M km², with sea level back at
−824.9 m, then drains again, to 11.65 M km³ over 11.09 M km² at step 12. Sea level follows, so
the cycle runs about eight steps. D fills at step 1 and drains over its four steps, too few to
refill. S's target fills to a 5.1 M km² lake, 2,303 m above its floor, drains to a 4.05 M km²
plateau for five steps, then rises again at step 12. On main, S's target alternates every step
between 260,000–390,000 km³ and 3,000–27,000 km³.

### What drives the cycle (post hoc)

Not preregistered. #754's pass tracer (`trace_lake_solve_passes.py`, not yet checked in) was run
on E's branch arm at `e5c697d2`, with its body-loss reading changed to count sublimation, which
this branch books apart from evaporation.

- **Refills come from flush passes, which charge no loss.** At steps 2 and 10 the anchor's
  sub-basin is at its saddle in the timed pass, spilling 222 M and 38 M km³. In the
  zero-duration flush passes that follow, the merged body receives 858 M and 205 M km³ and keeps
  24.6 M and 24.7 M km³, nearly all of each fill, with no evaporation or sublimation charged on
  it. This mechanism predates the change: #754 found the same flush-pass filling on main (3.9 to
  27.6 M km³ at 225.2 Myr), where after-solve sublimation stripped it within the step.
- **Drains are the basin relaxing toward an equilibrium that keeps moving.** On each drain step
  the timed pass is near balance, with loss exceeding supply by 0.5–1.4%. Over steps 3–9 the
  level falls about 1.5 km. The anchor's water-surface temperature, lagged through last step's
  level, rises from −4.0 °C to 10.2 °C, and its ice cover falls from 0.57 to 0.31. Loss per unit
  target area rises from 0.60 to 0.74 m/yr, while the basin's supply falls from 721 M to 542 M
  km³ per step. The equilibrium below the guard is this change's; the lag that moves it is
  #661's.

So the change removes the ice phase's period-2 artifact. What remains is a cycle of about eight
steps. Its fills come from flush passes that charge no loss over the step, a pre-existing breach
of principle 3. Its drains are the new equilibrium, kept moving by lagged feedbacks. This is one
traced window; S's step-12 rise wasn't traced.

### What this means

- **For #757:** sublimation now sets one balanced level with evaporation and supply, and its loss
  grows with step length; the unit tests show both. Cold lakes keep their water, so whole-world
  lake volume rises several times and sea level falls by 11–25 m at the windows' ends, and by up
  to 78 m while E's and D's basins are full.
- **Not fixed here:** the flush passes. When a split basin rejoins within a step, the step's
  surplus reaches the merged body with no loss charged over the step. The lagged surface
  temperature (#661) and wet-area rain partition also move a km-deep basin's equilibrium from
  step to step. Each needs its own issue and evidence. #758 would reclassify D and E's basin as
  ocean, which may remove this example, but not S's lake or the flush-pass mechanism.
- **For RFC 0016 G1b (#751):** D and E stay invalid as land-store tests, now for these reasons.
