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
