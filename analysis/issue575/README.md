# Issue 575: climate noise that persists across steps

Evidence for [#575](https://github.com/adubey/mantle-bloom/issues/575), a prerequisite of
[RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md)'s "No unforced swings" criterion.

**Principle:** #3, fast processes in aggregate over long steps. The textures stand for
geography, and an independent draw each step adds forcing that depends on step length.

**Invariant:** two climates 0.25 Myr apart on one surface differ by much less than two
independent draws (`test_climates_a_step_apart_differ_much_less_than_independent_draws`).

## Preregistration

Written before any replay.

**Change.** `_drifting_white_noise` blends independent keyframes with cos/sin weights over 40
Myr. The field stays standard normal at every instant, so nothing should move on average.

**Runs,** one at a time:

| Run | Tree | Save | Horizon |
|---|---|---|---|
| Climate A/B | detached `origin/main`, then the branch | seed 997271774 at 87.4 Myr, seed 649658380 at 6.4 Myr | 10 Myr at 0.25 Myr steps, `measure_climate_fields.py`, 1 Myr checkpoints |
| Arm (m3) | G0 prototype `ea9a187` plus this PR's `climate.py` change | seed 997271774 at 87.4 Myr | 10 Myr at 0.25 and 1 Myr steps, `replay_g0.py`, no `--freeze-climate-noise` |

**Primary metrics and thresholds:**

1. Arm (m3) passes G0 on seed 997271774 without the freeze: swing below 25%, step gap below
   10%. Expected near the frozen arm's 14.9% and 2.9%; it was 27.4% and 3.2% with the redraw.
2. Mean land precipitation and humidity, over the run's checkpoints, stay within 2% of
   `origin/main`. Expected about 0, since the noise has zero mean.
3. Land under visible ice shows no regression: the branch's end-of-run value lies within
   `origin/main`'s own checkpoint-to-checkpoint spread, or the gap is explained.

Each arm's checkpoint-to-checkpoint spread is the variability estimate. The branch's spread in
ice area should fall, since the redraw no longer moves it.

**Not measured:** a texture that holds for tens of Myr leaves its pattern on long-term
precipitation, where the redraw averaged out. These replays are too short to show any
landscape effect of that.

## Preregistration: is the land gap, or the step gap, one draw?

Written after the first runs, which missed two thresholds above: land humidity (+6.07% and
+2.81%) and land precipitation (+5.25% on seed 997271774); and arm (m3)'s step gap (12.2%
against 10%). The thresholds stay as written.

**Hypotheses.** The noise has the same distribution at every instant on both arms, so it can't
shift a mean directly. A shift is either one draw's overlap with this land, or persistence
acting through stored state (vegetation, lakes, rivers).

**Runs,** seed 997271774 at 87.4 Myr, one at a time:

| Run | Tree | Draws (`--noise-salt-offset`) |
|---|---|---|
| Climate, as above, recording moisture sources | `origin/main`, the branch | `origin/main` once; branch 0, 1, 2, 3 |
| Arm (m3), 0.25 and 1 Myr steps | G0 prototype with this PR's `climate.py` | 1, 2 |

**Decision rules:**

1. Land humidity and precipitation: if the four draws' mean gap to `origin/main` is within 2%,
   or the gaps change sign, the shift belongs to the draw and the "unchanged on average" check
   holds. If all four gaps exceed 2% with the same sign, the bias is real.
2. If the bias is real, the moisture sources show whether vegetation or water carries it.
3. Arm (m3): the change passes G0 only if the step gap is below 10% on all three draws. If it
   fails on some draws, the gap depends on the draw, and that goes to the owner with the numbers.

## Results so far (2026-10-31, interim)

Not final: the draw runs in the preregistration above haven't finished. Baseline `d106194`,
branch `0377e4e`, saves `08ad9461…` (seed 997271774) and `ebbfb5da…` (seed 649658380). Raw
outputs are in B2 as `issue575-replays-2026-10-31.tar.gz`.

| Branch vs `origin/main`, mean of checkpoints after the opening | Seed 997271774 | Seed 649658380 |
|---|---:|---:|
| Land humidity | +6.07% | +2.81% |
| Land precipitation | +5.25% | +1.24% |
| Ocean precipitation | −1.85% | −0.31% |
| Ice volume | +0.52% | −0.14% |
| Land under visible ice, end of run | 5.18% vs 5.60% | 2.43% vs 2.52% |

| Arm (m3), seed 997271774 | Swing, 0.25 / 1 Myr | Step gap | End volume, 0.25 / 1 Myr |
|---|---:|---:|---:|
| Per-step redraw (#554) | 27.4% / 10% | 3.2% | |
| Noise frozen (#554) | 14.9% / 2.4% | 2.9% | |
| Drifting noise | 3.1% / 2.5% | **12.2%** | 12.90 / 11.32 M km³ |

- Thresholds 1 (step gap) and 2 (land humidity and precipitation) are missed so far.
- The land gap appears within the first Myr, then fluctuates without growing.
- A run killed by the session's end at step 34 of 40 was discarded and rerun; its output is kept.
