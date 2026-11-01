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

## Results

Baseline `d106194`. Branch `0377e4e` for the first runs, then the same change rebased onto
`origin/main` for the draws; offset 0 reproduced the first run's land gaps to two decimals. Saves
`08ad9461…` (seed 997271774) and `ebbfb5da…` (seed 649658380). Raw outputs are in B2:
`issue575-replays-2026-10-31.tar.gz` and `issue575-replays-2026-11-01.tar.gz` (SHA-256
`e4cc13879a0e4dfe466cc586b28c067693503b1bd7053867fbeeff289fd708d0`), which also holds the
`--noise-salt-offset` patch for `replay_g0.py`. That patch applies to the G0 prototype only.

**Decision (owner, 2026-10-31):** G0's step-gap limit goes from 10% to 20%, in
[PR #735](https://github.com/adubey/mantle-bloom/pull/735). The first limit failed on every draw
of seed 997271774, so this report's second preregistration didn't decide this; the owner did.
Thresholds other than the step gap stay as written.

### Land climate, rule 1: holds

Branch vs `origin/main`, mean of checkpoints after the opening, seed 997271774:

| Metric | Offset 0 | 1 | 2 | 3 | Mean |
|---|---:|---:|---:|---:|---:|
| Land humidity | +6.07% | −3.73% | +1.62% | −1.85% | +0.53% |
| Land precipitation | +5.25% | −1.94% | +1.15% | −0.68% | +0.94% |
| Ocean precipitation | −1.85% | +0.58% | −0.74% | −1.08% | −0.77% |
| Ice volume | +0.52% | +0.74% | +0.69% | +0.39% | +0.59% |

The land gaps change sign and their means are within 2%, so the shift belongs to the draw. Single
draws still exceed 2%, so one seed's land can differ from `origin/main` by that much. Rule 2 isn't
needed. Seed 649658380 has one draw (+2.81% humidity, +1.24% precipitation), not four.

### Land under visible ice: within the range, all on one side

End of run on seed 997271774: 5.18%, 5.28%, 5.19% and 5.34%, against 5.60% for `origin/main`,
whose checkpoints after the opening span 5.07–5.64%. Every draw lies in that range, so read as a
range the threshold holds. Read as the largest change between checkpoints (0.34 points), offsets 0
and 2 exceed it (0.42 points). All four draws end below the single `origin/main` run, so a small
real drop in ice area isn't ruled out. Seed 649658380 has one draw: 2.43% vs 2.52%.

### Arm (m3), rule 3: passes both saves at the owner's 20%

| Save, drifting noise | Draw | Swing, 0.25 / 1 Myr | End volume, 0.25 / 1 Myr | Step gap |
|---|---:|---:|---:|---:|
| Seed 997271774 | 0 | 3.1% / 2.5% | 12.90 / 11.32 M km³ | 12.2% |
| | 1 | 19.2% / 4.1% | 10.53 / 12.31 M km³ | 17.0% |
| | 2 | 12.3% / 2.2% | 11.21 / 12.40 M km³ | 10.6% |
| Seed 649658380 | 0 | 5.9% / 5.2% | 2.11 / 2.21 M km³ | 5.1% |
| | 1 | 4.9% / 4.1% | 2.07 / 2.15 M km³ | 3.8% |
| | 2 | 6.0% / 4.6% | 2.25 / 2.35 M km³ | 4.5% |

| Earlier arms, #554 | Swing, 0.25 / 1 Myr | Step gap |
|---|---:|---:|
| Seed 997271774, per-step redraw | 27.4% / 10% | 3.2% |
| Seed 997271774, noise frozen | 14.9% / 2.4% | 2.9% |
| Seed 649658380, per-step redraw | 9.4% / 6.5% | 0.3% |

- Swing and step gap use #554's definitions, and offset 0 on seed 997271774 reproduces the first
  run. Offset 0's arm (m3) runs on that seed predate the rebase.
- Against the per-step redraw, the step gap rises on both saves. Its sign isn't consistent on seed
  997271774, so it isn't one step length gaining ice.
- **Why the gap grows isn't established.** The noise is a pure function of elapsed time, so both
  step lengths see one field, sampled at different times. No run here tests a cause.
- On seed 997271774, offset 0's 3.1% swing is the low end: the other draws reach 12–19%.
- **G0 isn't met on every diagnosis save.** The ice-age save (seed 374464662) wasn't rerun with
  drifting noise, and its runaway under the march stays with G4.
- The land gap appears within the first Myr, then fluctuates without growing.
- A run killed by the session's end at step 34 of 40 was discarded and rerun; its output is kept.
