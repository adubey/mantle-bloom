# Issue 741: climate noise over its 40 Myr period, and arm (m3)'s ice drops

Evidence for [#741](https://github.com/adubey/mantle-bloom/issues/741), which settles the banding
item of [#575](https://github.com/adubey/mantle-bloom/issues/575). Tracking:
[#544](https://github.com/adubey/mantle-bloom/issues/544).

**Principle:** #3, fast processes in aggregate over long steps. The noise stands for geography,
so swapping the redraw for a drift shouldn't move the climate on average.

**Invariant:** humidity and precipitation distributions by latitude band, averaged over a drift
period and over draws, match the per-step redraw's.

## Findings

1. **Banding: no shift on average, but land isn't shown within ±2%.** On 19 replays of 48 Myr, no
   metric is shifted on both saves. Ocean humidity is unchanged in every band, and global ocean
   precipitation is unchanged. Every land metric stays unresolved: its 95% interval contains 0
   but reaches past ±2% (global land precipitation +1.09% [−0.66, +2.84] and +0.20% [−3.38,
   +3.78]). The width comes from the drift's own spread: one world's land climate varies between
   drift draws about four times as much as between redraw draws.
2. **The 40 Myr period.** Global land offsets fade within a run: no draw holds a same-sign gap
   above 2% in all three windows. At 15–30°, 5 of 12 drift draws do, by up to 24%, for the whole
   48 Myr. That is short of "most draws", so by the preregistered rule it doesn't go further.
3. **Arm (m3)'s large ice drops are numerical.** On draw 1 the 19.2% drop happens in one 10 kyr
   substep of the march (13.10 to 5.93 M km³). The plateau ice there turns over about five times
   in 10 kyr. With the march's longest substep capped at 1 kyr, none of 12 perturbed runs drops,
   against 9 of 12 at 10 kyr. The noise isn't needed: holding it fixed leaves 8 of 12 dropping.
   Draw 2's 12.3% drop, on the same plateau body, goes too: 0 of 12 at the cap against 9 of 12.
4. **The +8.3% and −10.0% steps aren't changes in the ice.** Each closes inside its own march's
   scatter, and neither has a large substep.

**Owner's decisions, 2026-11-02:**
- **Finding 1:** #575's banding item is reworded to "no shift on average", and land within ±2%
  isn't required. The item is checked.
- **Finding 2:** the 15–30° persistence is left for now.
- **Finding 3:** it goes to a D2 revision, [#765](https://github.com/adubey/mantle-bloom/issues/765),
  and to G1 run on the march, [#766](https://github.com/adubey/mantle-bloom/issues/766).

## Preregistration

Written before any run in this report. The data already seen are #575's runs, including arm
(m3)'s 14 Myr series, whose per-step ledger motivated Part B.

### Part A: climate over 48 Myr

**The baseline isn't `origin/main`.** #736 merged the drifting noise, so `origin/main` is now the
treatment. The only backend change between `9d6fb5e` (before #736) and `bf58116` (`origin/main`
now) is `climate.py`, so `9d6fb5e` is the old `origin/main` exactly.

| Arm | Tree | Draws |
|---|---|---|
| Drift (D) | `bf58116` plus this branch's measurement script | offsets 0, 1, 2, 3 |
| Redraw (R) | `9d6fb5e` for draw 0; `bf58116` with `--noise-redraw` for the others | 0, 1, 2 |

**Replays,** one at a time: seed 997271774 from 87.4 Myr and seed 649658380 from 6.4 Myr, 48 Myr
at 0.25 Myr steps, a checkpoint each Myr. On seed 997271774 the run passes the keyframe at 120 Myr.

**Extension, fixed in advance.** If any primary metric's 95% half-width on seed 997271774 exceeds
2%, add D offsets 4–7 and R offset 3 there. This reads the width only, never the gap, so it can't
select for a result.

**Statistic.** Each run gives one number per metric: its mean over checkpoints after the opening.
Checkpoints within a run aren't independent, so the draw is the unit. The gap is
`(mean D − mean R) / mean R`, with a 95% Welch interval on the difference.

**Primary metrics:** global land and ocean precipitation and humidity, and the land precipitation
band at 15–30°. The rest are descriptive: every band's land and ocean means and p10, p50 and p90,
with about 40 metrics in each save, so about two 95% intervals should exclude 0 by chance.

**Decision rule,** per metric, for bands with at least 5% of land:

| Verdict | Condition |
|---|---|
| Unchanged | the interval lies inside ±2% |
| Shifted | the interval excludes 0 on both saves, with the same sign |
| Unresolved | anything else |

**Expected:** no shift. The noise has one distribution at every instant on both arms, so a real
shift would need a nonlinear path through stored state (vegetation, lakes, rivers), which the
redraw averages and a persistent pattern doesn't.

**The period itself.** Each D draw's gap to the mean of R, in windows 0–16, 16–32 and 32–48 Myr.
If the offset belongs to the texture's overlap with this land, it should weaken or change sign as
the keyframes hand over. If most draws keep a same-sign gap above 2% in all three windows, that's
persistence through stored state, and goes to the owner.

### Part B: what happens in the ice on the step of a large drop

**Tree.** G0 prototype `ea9a187` with #575's `climate.py`, arm (m3), seed 997271774, draw 1, 0.25
Myr steps, to 101.15 Myr. It must reproduce the earlier series; its ice volume at 97.4 Myr was
10.53 M km³.

**Method.** On each named step, re-march that step in-process from its real opening ice, with:

| Cell | Climate | Surface | Opening ice |
|---|---|---|---|
| `new_new` | this step's | this step's | real (must equal the real march) |
| `old_new` | last step's | this step's | real |
| `new_old` | this step's | last step's | real |
| `old_old` | last step's | last step's | real |
| `ice_plus`, `ice_minus` | this step's | this step's | real ±0.1% |

Climate is precipitation, temperature and amplitude at nodes. Surface is elevation, deflection,
mobile cover, and the ocean masks. `world` itself, for hypsometry and trapped water, stays this
step's in every cell.

**Steps:** the drop at 97.4, the step after at 97.65, the swing pair 100.4 and 100.65, and controls
96.4 and 99.15.

**Primary measures:** the cell's closing volume against `new_new`, and the march's within-step
scatter (standard deviation of substep volumes in the step's second half).

**Decision rules:**

1. Controls give the cells' noise: the largest |`old_old` − `new_new`| across control steps.
2. A group of inputs causes a drop if restoring last step's values alone removes at least half of
   it, by more than that noise.
3. If neither group does but `old_old` still drops, nothing changed in the inputs: the march is
   leaving a state of its own accord.
4. If `ice_plus` and `ice_minus` differ by as much as the drop, closing volume is sensitive to its
   start, and one step's change can't be attributed at all.
5. If a step's change is within its march's scatter, that step-end value is a sample of the
   scatter, not a change in the ice.

A climate result leads to one more run, holding the noise from 97.15 Myr, because the climate at
nodes also changes through the surface. Where the ice went, by body and region, comes from the
same captures.

### Part B amendment, after the first factorial

Written after the first factorial run and before any run below. It is post hoc: the design
changed because of a result.

Rule 4 fired on the drop step: ±0.1% opening ice gave 13.14 and 9.70 M km³, a gap larger than the
drop (13.03 to 10.53). One march can't be attributed, so the question becomes a distribution. The
capture also showed the moving ice on equatorial plateaus near 5 km, with nodes going between 0
and about 2 km in a step.

| Run | What | Reads |
|---|---|---|
| E1 | 96.4, 97.4 and 99.15: 12 uniform opening-ice scalings (±0.01–0.5%), this step's climate; on 97.4 also last step's climate | Does the drop step differ from ordinary steps, and do last step's inputs change the odds? |
| E2 | 97.4 with the noise held from 97.15 Myr, same 12 scalings | Is the noise's own change needed? |
| E3 | Scale the 97.4 opening by 1.001 and continue to 101.15 Myr | Does an avoided drop stay avoided, or does the low state arrive anyway? |
| E4 | Repeat the factorial on 97.4 and 100.4 | The node captures a filename bug lost; volumes must match the first run |

**Reading:** the drop counts as forced if at least 10 of 12 members drop (closing below 11.5 M km³)
under this step's inputs and at most 2 do at control steps. It counts as a sensitive switch if the
fraction lies between. E1 under last step's climate, and E2, are read as a change in that fraction;
with 12 members, a gap under 5 members isn't called a difference.

### Part B second amendment, after the ensembles

Written after E1–E4 and before any run below, so it is post hoc too.

E1–E3 found a change in the climate that isn't the noise, then E4's march path showed its form: on
the drop step one substep takes the volume from 13.1 to 5.9 M km³, and it rebounds to about 9.8.

| Run | What | Reads |
|---|---|---|
| E5 | Draw 1, 97.4: log every substep, capture the solve's node fields on any that moves volume over 15% | What the collapsing substep does, and where its ice goes |
| E6 | Draw 0, 97.4: the 12-member ensemble | Does a draw with no drop have the same switch? |
| E7 | Draw 2, 96.65, its largest step: the ensemble and substep log | Is draw 2's −12.3% the same event? |

**Reading.** The collapse is a defect of the solve if the ledger can't account for it as physical
melt, or if one substep removes ice no flux or melt rate at that substep's length could. It's the
same event on draw 2 if that step holds one substep over 30% on the same plateau body. E6 and E7 use
the same 5-member gap.

### Part B third amendment, after E5

Written after E5 and E6 and before E8, so it is post hoc too.

E5's collapsing substep is 10 kyr long, the march's longest, while the ice turns over about five
times in 10 kyr. The hypothesis is that a substep much longer than the ice's turnover lets the
march overshoot from the plateau's high state to its low one.

| Run | What | Reads |
|---|---|---|
| E8 | Draw 1, 97.4: E1's 12-member ensemble with the march's longest substep capped at 1 kyr | Does the drop need long substeps? |

**Reading.** Long substeps cause the drop if at most 2 of 12 members drop (against E1's 9) and no
member has a substep that loses over 20%. If 8 or more still drop with no such substep, the drop is
a real transition the long substeps only hasten. Anything between is reported as it stands.

### Part B fourth amendment, after E7 and E8

Written after E7 and E8 and before E9. E8 met the rule for "long substeps cause the drop" on draw
1. E7's draw-2 drop, on the same plateau body, has no substep over 15% but uses 2.5–10 kyr substeps.

| Run | What | Reads |
|---|---|---|
| E9 | Draw 2, 96.65: E7's ensemble with the longest substep capped at 1 kyr | Same rule as E8, against E7's 9 of 12 |

**Not tested:** other draws' drops beyond these. If time allows, the same factorial runs on the largest step of
draw 2. Any change to the drift period, the noise or the march is proposed separately.

## Results

### What was run

| Runs | Tree | Save (SHA-256) | Steps |
|---|---|---|---|
| Climate, drift draws 0–7 | `bf58116` (`origin/main`) with this branch's `measure_climate_fields.py` | seed 997271774 at 87.4 Myr (`08ad9461…`) | 192 × 0.25 Myr |
| Climate, redraw draw 0 | `9d6fb5e` (before #736) | same | same |
| Climate, redraw draws 1–3 | `bf58116`, `--noise-redraw` | same | same |
| Climate, drift 0–3, redraw 0–2 | as above | seed 649658380 at 6.4 Myr (`ebbfb5da…`) | same |
| Arm (m3) probes | G0 prototype `ea9a187` with #575's `climate.py` (`a94047d`), plus a probe patch | seed 997271774 at 87.4 Myr | 0.25 Myr |

```bash
# Climate, from backend/ in this branch's worktree:
.venv/bin/python ../bin/debug/measure_climate_fields.py <local-save> --myr 48 --step-years 250000 \
    --checkpoint-myr 1 [--noise-salt-offset N] [--noise-redraw] [--backend <before-736-worktree>/backend] \
    --out <scratch-output>/<save>-<arm><draw>.json
.venv/bin/python ../bin/debug/summarize_climate_noise_replays.py <scratch-output> --out ../analysis/issue741/climate.json
# Arm (m3), from backend/ in the prototype worktree with the probe patch:
.venv/bin/python ../bin/debug/replay_g0.py <local-save> --arm m3 --step-years 250000 --checkpoint-myr 0.25 \
    --noise-salt-offset 1 --myr 10 --factorial-myr 97.4 --ensemble-only [--ensemble-old-climate] \
    [--hold-noise-from 97.15] [--ensemble-max-substep-years 1000] [--substep-log-myr 97.4] --capture-dir <scratch-output>/e1
```

Climate replays took about 7 s a step (24 minutes a run) and peaked at 1.0–1.6 GB. One, seed
649658380's drift draw 3, was paused for an hour and resumed, so its step time reads 30 s; its
results aren't affected. Arm (m3) took 40–100 s a step at 1.1–1.6 GB, and each ensemble member
about 1 minute (2–2.5 minutes at the 1 kyr cap). Every run ran alone. Raw outputs, the probe patch
and the run scripts are in B2 as `issue741-runs-2026-11-02.tar.gz` (SHA-256
`672998797dc10ca9fd5ef1dbff6bac1b7cec2eb2d0e32418e8c2340a9918abc2`). `climate.json` and `ice.json`
here are their compact summaries.

### Part A: climate

Drift against redraw: mean of each run's checkpoints after the opening, then the gap between the
arms' means, with a 95% Welch interval. Seed 997271774 has 8 drift and 4 redraw draws, after the
preregistered extension; seed 649658380 has 4 and 3.

| Metric | Seed 997271774 | Seed 649658380 | Verdict |
|---|---|---|---|
| Land precipitation | +1.09% [−0.66, +2.84] | +0.20% [−3.38, +3.78] | Unresolved |
| Land humidity | +0.24% [−2.09, +2.57] | +1.04% [−2.13, +4.22] | Unresolved |
| Ocean precipitation | −0.29% [−0.94, +0.36] | −0.09% [−0.92, +0.74] | Unchanged |
| Ocean humidity | +0.01% [−0.08, +0.10] | +0.01% [−0.19, +0.21] | Unchanged |
| Land precipitation, 15–30° | +3.56% [−0.74, +7.87] | −2.97% [−18.49, +12.56] | Unresolved |

By band (land and ocean means; the p10, p50 and p90 of each are in `climate.json`):

| Band (share of land) | Land precipitation | Land humidity | Ocean precipitation | Ocean humidity |
|---|---|---|---|---|
| 0–15° (30%) | +0.42 / +2.44 | −0.39 / +3.36 | −0.18 / −0.04 | Unchanged |
| 15–30° (24–27%) | +3.56 / −2.97 | +0.78 / +0.62 | −0.68 / −1.23 | Unchanged |
| 30–45° (19–20%) | +0.72 / −3.03 | −0.23 / −1.49 | −0.41 / −0.23 | Unchanged |
| 45–60° (13–16%) | +1.90 / −0.71 | +1.34 / −0.20 | −0.19 / +0.38 | Unchanged |
| 60–75° (6–9%) | +1.84 / +0.24 | +0.21 / +0.31 | Unchanged | Unchanged |
| 75–90° (3–4%) | Under 5% of land | Under 5% of land | +0.31 / +0.28 | Unchanged |

Gaps in %, seed 997271774 / seed 649658380. Every cell not marked unchanged is unresolved; none is
shifted. Of 104 metrics, 35 are unchanged, 61 unresolved and 8 exempt; 7 intervals exclude 0 on one
save, about what 208 intervals at 95% give by chance, and none on both.

- **The spread is the drift's.** Land precipitation varies by ±2.04% between drift draws on seed
  997271774 against ±0.51% between redraw draws, and by ±5.11% against ±0.84% at 15–30°. A
  redraw run averages about 190 independent draws; a drift run sees about two.
- **#575's +5.2% at 15–30° was the draw.** Its first four draws gave +0.50% on this measure;
  eight give +3.56%, and seed 649658380 gives −2.97%.
- **Ice volume** isn't a primary metric: +3.44% [+2.47, +4.40] on seed 997271774, where every drift
  draw lies above every redraw draw, and −0.70% [−1.46, +0.05] on seed 649658380. The signs differ,
  so it isn't a shift. Land under visible ice: +0.96% and −1.46%, both unresolved.

**The period.** Each drift draw's land-precipitation gap to the redraw mean, by window:

| | 0–16 Myr | 16–32 Myr | 32–48 Myr |
|---|---|---|---|
| Global, seed 997271774 | −2.95 to +4.56 | −1.40 to +4.32 | −2.31 to +2.09 |
| 15–30°, seed 997271774, draw 5 | +12.04 | +13.95 | +12.84 |
| 15–30°, seed 997271774, draw 7 | +2.31 | +11.35 | +10.18 |
| 15–30°, seed 649658380, draw 2 | −14.41 | −24.16 | −13.45 |

The first and last windows share little noise: at 87.4 Myr the field is 0.96 of keyframe 2 and
0.28 of keyframe 3, and by 135.4 Myr it's 0.82 of keyframe 3 and 0.57 of keyframe 4, a
correlation of about 0.23. A band offset that holds through both is therefore either two
keyframes that happen to favour the same land or state carried by the surface, and these runs
can't tell which. Five of 12 draws is short of the preregistered "most", so this goes to the owner
as an observation, not a finding.

**Judgement for the owner.** Nothing shifts, and the ocean is unchanged. The land can't be shown
within ±2% at these draw counts, and the reason is the design: a slow texture leaves each world's
land more or less wet for tens of Myr. More draws would narrow the average but not that spread.
Checking #575's box on "no shift" is defensible; on "within ±2%" it isn't. The owner chose "no shift"
(2026-11-02).

### Part B: arm (m3)'s ice drops

All on seed 997271774, draw 1 unless stated. The first run reproduced #575's series exactly.

**Factorial on the drop step, 97.4 Myr** (opening 13.03 M km³; the real march closes at 10.53):

| Cell | Climate | Surface | Opening ice | Closing (M km³) |
|---|---|---|---|---|
| `new_new` | this step | this step | real | 10.53, equal to the real march |
| `old_new` | last step | this step | real | 12.74 |
| `new_old` | this step | last step | real | 10.54 |
| `old_old` | last step | last step | real | 12.96 |
| `ice_plus`, `ice_minus` | this step | this step | ±0.1% | 13.14, 9.70 |

Rule 4 applied: 0.1% of opening ice moved the result by more than the drop, so the amendment
replaced single cells with ensembles. A filename bug overwrote this run's node captures for 97.4
and 100.4; E4 repeated it and matched every volume.

**Ensembles,** 12 opening-ice scalings from ±0.01% to ±0.5%. "Drops" means closing below 11.5 M
km³; at 99.15 Myr, after the drop, the opening is already 10.94, so there it means a fall of over
10% from the opening. The first amendment didn't foresee that, and read literally it counts all 12.

| Step | Condition | Drop / 12 | Closing range (M km³) |
|---|---|---|---|
| 96.4 (control) | this step's climate | 0 | 13.01–13.18 |
| 99.15 (control) | this step's climate | 0 | 11.06–11.48 |
| 97.4 | this step's climate | 9 | 9.70–13.20 |
| 97.4 | last step's climate | 2 | 9.77–13.03 |
| 97.4 | noise held from 97.15 | 8 | 9.67–13.34 |
| 97.4, draw 0 | this step's climate | 2 | 10.97–12.98 |
| 96.65, draw 2 | this step's climate | 9 | 10.91–12.35 |
| **97.4** | **longest substep 1 kyr** | **0** | **13.14–13.25** |
| **96.65, draw 2** | **longest substep 1 kyr** | **0** | **12.27–13.21** |

E3 scaled the 97.4 opening by 1.001 and continued: the drop was avoided on that step and came on
the next (13.4 to 6.7 within one substep, closing at 9.81).

**What happens on the drop step** (E5, the real march):

| Substep | Length | Volume (M km³) | Accumulation | Surface melt | Visible nodes |
|---|---|---|---|---|---|
| 63 | 10 kyr | 13.12 → 13.10 | 76.1 | 74.2 | 3,796 |
| 64 | 10 kyr | 13.10 → 5.93 | 46.3 | 52.6 | 2,378 |
| 65 | 10 kyr | 5.93 → 12.73 | 63.9 | 55.7 | 3,292 |
| 66–112 | 0.6–5 kyr | settles at 9.4–11.0 | | | 3,200–3,500 |

Accumulation and melt in M km³ over the substep.

- **Where.** One body loses 2.54 M km³, more than the net 2.51, on an equatorial plateau (7.7° N,
  about 5.4 km, −11 °C mean), which held 10.47 of the 13.03 M km³; 504 of its 2,574 nodes went
  to bare ground. Its climate moved little on the step: +0.04 °C and +0.9% precipitation.
- **How.** The ledger closes: melt less accumulation, with basal melt, is the 7.17 M km³ lost. But
  10 kyr carries about five times the body's volume through it, and the march had grown its
  substep to 10 kyr because the previous changes were under 1%.
- **Read against the rules.** Rule 4 fired, so the cells alone attribute nothing. The ensembles'
  9 of 12 is one short of "forced". Under last step's climate the count falls by 7, so the
  climate's change on the step matters; with the noise held it falls by 1, so the noise's change
  doesn't. The capped ensemble meets the third amendment's rule for long substeps: 0 of 12, and
  no substep moves the ice by more than 1.5%.
- **Draw 2** drops 12.3% at 96.65 Myr on the same plateau body (1.65 M km³ of the loss), but
  gradually over the step, with no substep over 8.6%. Its march still runs 2.5–10 kyr substeps.
  At the 1 kyr cap, none of 12 drops and no substep moves over 1.7%, but the closes split: 8 near
  12.3 M km³ (6% below the 13.12 opening) and 4 near 13.2. The body still has two states close
  together; the cap removes the large fall between them, not the choice.

**The swing pair, 100.4 and 100.65 Myr.** Neither march has a substep over 13%. The march
oscillates through each step (standard deviation 0.46 M km³ over the second half of 100.65), and
the +8.3% step's close of 12.09 sits above every other cell (10.39–11.19), so by rule 5 both
values are samples of that scatter.

**What this means for RFC 0015.** In this prototype, how far the march's substep may grow sets
whether this body switches state. A cap tied to the ice's turnover, as G1's τ ≤ 0.1 Δt asks of the
whole step, would remove both draws' drops. Its cost isn't measured here beyond the 2–2.5× ensemble
time. This report proposes nothing: a change to D2 goes through `edit-rfc`.
