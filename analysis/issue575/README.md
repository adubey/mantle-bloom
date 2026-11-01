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

## Preregistration: does the result hold on the merge base, and does draw 1 recover?

Written after review of the results below, before these runs. The review found that no reported
climate run used this PR's base, which adds persistent river water, and that draw 1 on seed
997271774 ends on a 19.2% drop with no recovery seen.

**Runs,** seed 997271774 at 87.4 Myr, one at a time:

| Run | Tree | Horizon |
|---|---|---|
| Climate A/B | `origin/main` at the PR's base, then the branch | 10 Myr at 0.25 Myr steps |
| Arm (m3), draw 1 | G0 prototype with this PR's `climate.py` | 14 Myr at 0.25 and 1 Myr steps |

**Decision rules:**

1. Base: the land humidity and precipitation gaps stay within 1 point of offset 0's +6.07% and
   +5.25%, and end-of-run ice area keeps its sign against `origin/main`. Otherwise the earlier
   climate tables don't carry over to this base, and the draws are rerun.
2. Draw 1: the first 10 Myr repeat the earlier run exactly, or the comparison is void. If every
   step to 14 Myr stays under the 25% swing limit and the end gap under 20%, it passes at that
   horizon. If either limit is broken, arm (m3) fails on this draw and that goes to the owner.
3. Arm (m3) stays on the G0 prototype's older base, so this measures noise, not river water.

## Results

Climate draws ran on `938180d` (base `0f4aac8`) against baseline `d106194`. The merge-base pair ran
at `8ebea73` and `d239b83` (that commit merged with this change). `main` moved on after that pair,
and it wasn't rerun. Arm (m3) ran on the G0 prototype's older base, so it measures noise only. Saves
`08ad9461…` (seed 997271774) and `ebbfb5da…` (seed 649658380). Raw outputs are in B2, with SHA-256:
- `issue575-replays-2026-10-31.tar.gz`: `a68d2d24ec2b3b832e3122973694fd2a4b413b63570950df4f1801b7a9056da5`
- `issue575-replays-2026-11-01.tar.gz`: `e4cc13879a0e4dfe466cc586b28c067693503b1bd7053867fbeeff289fd708d0`,
  with the `--noise-salt-offset` patch for `replay_g0.py`, which applies to the G0 prototype only
- `issue575-reruns-2026-11-01.tar.gz`: `046251c24465908dfb138f94d195ad6cbf316c57c2a3f3d577dbca848232c110`

**Decision (owner, 2026-10-31):** G0's step-gap limit went from 10% to 20% in
[PR #735](https://github.com/adubey/mantle-bloom/pull/735). The 10% limit failed on every draw of
seed 997271774, and the preregistered rules didn't cover that, so the owner decided it after seeing
these results. The other thresholds stay as written.

This section replaces the first version merged in
[PR #736](https://github.com/adubey/mantle-bloom/pull/736). Review found that it overstated the
rebase check and the ice-area result, and misread the step gap. The numbers in its tables were right.

### Summary against the issue's acceptance list

| Item | Status |
|---|---|
| Same-state renders agree | Met, and tested through `compute_climate` |
| Banding statistics unchanged | **Not shown.** No shift detected in land means, but a band moves; see below. Moved to #741 |
| Step-apart unit test | Met |
| Precipitation and land-ice area, no regression | Met by the owner's decision (2026-11-01); see below |
| Arm (m3) on seed 997271774 passes G0 | Met at the 20% limit, to 14 Myr on draw 1; G0 isn't met on every save |

### Land climate: no shift detected, which isn't the same as none

Branch vs `origin/main`, mean of checkpoints after the opening, seed 997271774:

| Metric | Offset 0 | 1 | 2 | 3 | Mean |
|---|---:|---:|---:|---:|---:|
| Land humidity | +6.07% | −3.73% | +1.62% | −1.85% | +0.53% |
| Land precipitation | +5.25% | −1.94% | +1.15% | −0.68% | +0.94% |
| Ocean precipitation | −1.85% | +0.58% | −0.74% | −1.08% | −0.77% |
| Ice volume | +0.52% | +0.74% | +0.69% | +0.39% | +0.59% |

The gaps change sign and their means are within 2%, so rule 1 holds. That rule is weak: it was
written after the misses, and sign changes are likely with four noisy draws. The means have standard
errors of about 2.2 points (humidity) and 1.6 points (precipitation), and there's one `origin/main`
run, so these runs can't resolve a shift of ±2%. Seed 649658380 has one draw (+2.81%, +1.24%).

Land precipitation by latitude, same runs:

| Band | `origin/main` (mm) | Offset 0 | 1 | 2 | 3 | Mean |
|---|---:|---:|---:|---:|---:|---:|
| 0–15° | 1,082 | +4.8% | −1.2% | −1.0% | −2.4% | +0.0% |
| 15–30° | 140 | +6.0% | +5.2% | +6.8% | +3.0% | +5.2% |
| 30–45° | 320 | +12.9% | −1.7% | −3.6% | +2.4% | +2.5% |
| 45–60° | 1,386 | +4.5% | −6.0% | +5.5% | +0.6% | +1.2% |
| 60–75° | 485 | +0.6% | +7.9% | +0.9% | −0.7% | +2.2% |
| 75–90° | 82 | +8.8% | −14.1% | −24.1% | −2.0% | −7.9% |

The 15–30° band gains on all four draws by more than 2%, which is rule 1's "real bias" condition,
though the rule was written for global land means. It's a dry band, and its base is one run.

### Land under visible ice: lower on every draw, not a regression (owner)

Seed 997271774 ends at 5.18%, 5.28%, 5.19% and 5.34%, against 5.60% for `origin/main`.

- The checkpoint range of `origin/main` is 5.07–5.64%, so every draw lies inside it, but the range
  includes an early transient. Its last three checkpoints are 5.56, 5.64 and 5.60, and no branch
  checkpoint on any draw exceeds 5.38.
- The mean over checkpoints is lower on all four draws: −5.37%, −3.58%, −4.81% and −1.32%, mean
  −3.77%. Rule 1 doesn't cover ice area, but by analogy this would fail both ways of holding.
- Seed 649658380, one draw, ends at 2.43% against 2.52%, below `origin/main`'s range of 2.51–2.66%.
  The gap is within `origin/main`'s largest change between checkpoints (0.10 points).
- The preregistered prediction that ice-area spread falls holds on seed 997271774 (0.56 to 0.25–0.36
  points) and reverses on seed 649658380 (0.15 to 0.27 points).

**Decision (owner, 2026-11-01):** this isn't a regression. The cause isn't known.

### Arm (m3): passes both saves at 20%

| Save, drifting noise | Draw | Swing, 0.25 / 1 Myr | End volume, 0.25 / 1 Myr | Step gap |
|---|---:|---:|---:|---:|
| Seed 997271774 | 0 | 3.1% / 2.5% | 12.90 / 11.32 M km³ | 12.2% |
| | 1 | 19.2% / 4.1% | 10.53 / 12.31 M km³ | 17.0% |
| | 1, to 14 Myr | 19.2% / 4.1% | 10.92 / 12.37 M km³ | 13.4% |
| | 2 | 12.3% / 2.2% | 11.21 / 12.40 M km³ | 10.6% |
| Seed 649658380 | 0 | 5.9% / 5.2% | 2.11 / 2.21 M km³ | 5.1% |
| | 1 | 4.9% / 4.1% | 2.07 / 2.15 M km³ | 3.8% |
| | 2 | 6.0% / 4.6% | 2.25 / 2.35 M km³ | 4.5% |

| Earlier arms, #554 | Swing, 0.25 / 1 Myr | Step gap |
|---|---:|---:|
| Seed 997271774, per-step redraw | 27.4% / 10% | 3.2% |
| Seed 997271774, noise frozen | 14.9% / 2.4% | 2.9% |
| Seed 649658380, per-step redraw | 9.4% / 6.5% | 0.3% |

Signed gap, (V(0.25 Myr) − V(1 Myr)) / V(0.25 Myr), at each 1 Myr checkpoint from +2 to +9 Myr and
then at the end, on seed 997271774:

| Draw | +2 to +9 Myr | End |
|---|---:|---:|
| 0 | +4.9% to +10.3% | +12.2% |
| 1 | +2.8% to +8.0% | −17.0% at +10 Myr, −13.4% at +14 Myr |
| 2 | +4.6% to +8.5% | −10.6% |

- **The 0.25 Myr run holds more ice than the 1 Myr run until a late drop.** The sign flips on draws 1
  and 2 only after one. So the end-of-run gap on these draws mostly records when a drop landed.
- **Draw 1's drop is a shift, not a blip.** Ice falls from 13.03 to 10.53 M km³ on one step at
  97.4 Myr and doesn't return to 13 by 14 Myr, peaking at 12.09. It moves +8.3% and −10.0% on
  successive steps at 100.4 and 100.65 Myr. After the drop the gap stays between −4.1% and −17.0%.
- Draw 1's 14 Myr runs repeat its earlier runs exactly for the first 10 Myr at both step lengths, and
  they meet rule 2 of the third preregistration: no swing reaches 25% and the end gap is under 20%.
- **Smooth noise doesn't remove the drops.** A 19.2% one-step drop under a field that changes about
  1% per step exceeds the frozen arm's 14.9%, and the gaps are several times its 2.9%. Fast ice
  following the redraw explains part of the swing at most. No run here finds the rest.
- "The swing falls" compares three draws with one redraw run, and seed 997271774's draws reach
  12–19%. Offset 0's 3.1% is the low end.
- Against the per-step redraw, the step gap rises on both saves. Seed 649658380 is the clearer case.
- **G0 isn't met on every diagnosis save.** The last arm (m3) run on the ice-age save failed its
  forcing conditions, and it hasn't been rerun with drifting noise.
- "Step robustness" stays at 10%. Arm (m3) exceeds that on every draw of seed 997271774.

### Merge-base check

Offset 0 at `8ebea73` against `origin/main` at the same commit, seed 997271774: land humidity +5.73%,
land precipitation +5.02%, mean land under visible ice −5.18%, ending 5.20% vs 5.55%. Rule 1 of the
third preregistration holds, so the climate tables carry over to that base. It's one draw. The
earlier check, that offset 0 reproduced after a rebase, showed nothing: that rebase changed no
simulation code.

### Not measured

- The 40 Myr period. A 10 Myr run covers a quarter of one keyframe interval, and a pattern now holds
  in the grid frame for tens of Myr. Offset 0 stays about 6% wetter on land for its whole run.
- Distributions beyond these bands and Köppen shares.
