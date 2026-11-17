# Issue 888: RFC 0018 D1, rift-line prototype and its three checks

Interim report for [#888](https://github.com/adubey/mantle-bloom/issues/888). The prototype is the
rift-line version of RFC 0018 (PR #897 revision) in `bin/debug/rift_line.py`,
`rift_damage_mechanism.py` and `rift_damage_step.py` (draft PR #890). No parameter search has been
run: the RFC requires the three checks below to pass first. Development seeds only.

| Check | Command | Result |
|---|---|---|
| Propagation | `check_rift_line_d1.py --only propagation` | **Pass.** With the amended integrator the steady tip speed is within 2% of the discrete law at all 8 settings (default and half density, tip length 60/125/250/500 km) |
| Integration | `check_rift_line_d1.py --only integration` | **Fails as specified** (strain error 0.29–0.41 against 0.05). **Passes with the amended integrator** (0.035), see below |
| Line (20 worlds, 0 and 50 Myr) | `check_rift_lines_world_d1.py` | **Fails the validity part:** 66% of eligible plates give a valid split (need 90%). Craton part passes (89% of plates with craton) |

Compact results: `integration-propagation.json`, `line-check-summary.json`. Per-plate rows and
checksums are in the private captures directory.

## Cell spacing

The RFC says the default cell spacing is 62.5 km. It is 125 km (`line_spacing_rad(1.0)`), equal to
`w_0`. Consequences at default density: the band is one cell (the line); the non-local softening
disk (radius `w_0/2`) holds only the cell itself; the tip mesh factor is 0.58 (`x = 1` at
`ℓ_tip = 125 km`), not 0.77. With `ℓ_tip = 125 km` the steady speed is 58 km/Myr, not 77, and two
tips cross a 5,000 km line in about 45 Myr, not 34. With `ℓ_tip = 250 km` (now the starting value) the
speed is 154 km/Myr, 18 Myr to cross. Half-density tip speed relative to default: 0.56, 0.78, 0.89, 0.95
for `ℓ_tip` 60, 125, 250, 500 km.

## Integration check

Manufactured plate, fixed pull, straight line, steps 0.01 (reference), 0.5, 1, 2, 5 Myr, `ℓ_tip = 250 km`.
The integrator as written (frozen `N_cap` per step, sub-steps of `Δε_max = 0.05`) is first-order,
because a line cell reaching `ε_cand` changes its neighbours' rate only at the next sub-step, and
`N_cap` follows thinning only at the next step. Worst strain error by variant at the starting rate:

| Variant | 0.5 Myr | 1 | 2 | 5 |
|---|---:|---:|---:|---:|
| RFC as written | 0.29 | 0.34 | 0.36 | 0.41 |
| + stop sub-steps at `ε_cand` crossings | 0.10 | 0.14 | 0.20 | 0.30 |
| + `N_cap` of line cells follows thinning in the sub-steps | 0.068 | 0.073 | 0.076 | 0.078 |
| + `Δε_max = 0.02` (the amended integrator, now the default) | 0.035 | 0.035 | 0.035 | 0.036 |

The amended integrator passes every criterion (strain within 0.05, broken sets within one cell per
tip, split time within one coarse step), with at most 36 sub-steps per 1 Myr step at the starting
rate (103 at 5 Myr). The remaining error depends on `Δε_max`, not on the step. The amendments are
proposed in the RFC PR for #888.

## Line check

At 0 Myr (277 eligible plates) and 50 Myr (387), production splitter, 1 Myr steps; 0.5 Myr steps on
three seeds (115 plates) give the same picture. Two runs, both with the amended integrator unused
(no straining):

| | Valid share of lines | Valid share of eligible plates | Craton: line beats chord |
|---|---:|---:|---:|
| First run (vote rule, islands counted) | 75.6% | 66.1% | 88.8% |
| After the fixes below | **91.0%** | 79.7% | 88.8% |

The RFC's pass rule is 90% of lines (plates that have a line) and two thirds on craton, so the
second run passes both. 83 of 664 eligible plates (12.5%) have no line, because no direction has a
positive pull peak; the share of all eligible plates is reported alongside.

Fixes to the split (readings of the RFC, in `rift_line.split_masks`): cells isolated in the
plate's own adjacency (islands) are left out of the stray and connectivity tests and go to the
nearest daughter; line cells with no labelled neighbour wait for one, so each daughter stays
connected. The diagnosis (`diagnose_line_placement.py`, 8 worlds) found 5 of 9 invalid lines on one
world were islands and the breadth-first rule alone rescued 2 of 28 disconnections.

Remaining rejections: 50 of 52 are a stray piece above the 1% limit (median 3.4% of the plate),
1 no separation, 1 daughter under 15%. At 0.5 Myr steps 93–95% of lines are valid. Valid lines
split off a median 33% of the plate, are 1.30 times the straight distance between their ends, and
run 66% over oceanic crust; median `Λ` is 0.08 against `Λ_grow = 0.3` (4% of continental lines
reach it).

## Stage 2

The three stage 1 candidates (settings 0, 10 and 14 of `stage1-settings.json`, chosen by the
corrected scoring rules recorded on the issue) on all 20 development worlds: 200 Myr at 1 Myr
steps, plates-only, shadow shell solve every 10 Myr. 60 of 60 runs finished, largest ledger gap
7e-16. Two setting 0 runs first crashed in the shadow solve on a plate with fewer than four
non-coplanar cells; `shadow_stats` now skips and counts such a plate, and both were rerun from the
start. Scores are in `stage2-scores.json` (`score_stage1.py --runs <stage2> --baseline
<production-baseline>`); the production splitter makes 7.4 splits per world per 100 Myr on the
same worlds (`production-baseline-summary.json`).

| | Setting 0 | Setting 10 | Setting 14 | Target |
|---|---:|---:|---:|---|
| Targets met (of 10) | 6 | 6 | 6 | |
| Total violation | 2.50 | 2.21 | 2.48 | |
| Splits | 174 | 55 | 131 | |
| Worlds with a split | 20 | 17 | 20 | |
| Onset width (km) | 126 | 128 | 129 | 100–150 |
| Share of rifts wider at breakup | 0.92 | 0.94 | 0.93 | above 0.5 |
| Completed-rift median (Myr) | 25 | 35 | 25 | 10–50 |
| Transition ratio | **0.93** | **0.75** | **1.14** | at least 3 |
| Splits along inherited weakness | **16%** | **13%** | **17%** | at least 67% |
| Failed share | **82%** | **81%** | 80% | 20–80% |
| Localization violations | **9** | 0 | **8** | 0 |
| Smallest daughter | 15.3% | 16.7% | 15.6% | at least 15% |
| Line length ratio | 1.30 | 1.30 | 1.29 | at least 1.1 |
| Splits per world per 100 Myr | 4.35 | **1.38** | **3.28** | at least 3.7 |
| Oceanic share of successful rifts | 90% | 89% | 93% | fallback above 33% |
| Splits with a plume on the line | 30% | 29% | 58% | reported |

Bold values miss their target. By the predeclared ranking (targets met, then smaller total
violation) setting 10 is first, but it is the one that splits least: 1.4 per world per 100 Myr and
none at all on three worlds.

What the three settings share:

- **Rifts are oceanic.** 89–93% of successful rifts run mostly through oceanic crust (17, 6 and 9
  continental splits in settings 0, 10 and 14). This is above the one-third limit in every
  setting, so the RFC's oceanic fallback applies whichever is frozen.
- **No acceleration toward breakup.** Transition ratios are 0.75–1.14 against 3.
- **Lines do not follow inherited weakness.** Split lines are 25–28% weak rock against 23–25% for
  matched random paths.
- **About four fifths of candidates fail**, at the edge of the 20–80% range.

Localization fails only in the two settings that split often (largest damaged share 57% and 55%
of a plate, up to 5 damaged components at a split). Super-sea counts are 6,266, 1,614 and 4,380
checks, with the caveat recorded on the issue: they count checks, not basins.

Also reported: held band cells 81, 3 and 127; far flank transfers 8,387, 1,450 and 16,813; breakup
tests that did not give a valid split 513 of 687, 87 of 142 and 296 of 427; no strain off the
line; at most 30, 52 and 68 sub-steps per step.

### Shadow shell

D0's statistics for proxy C from 400 shadow solves per setting (20 worlds, every 10 Myr). Each
value is the median over the per-world, per-time checks; D0 pooled all plates of its 28 snapshot
worlds into one statistic, so the two columns are not the same estimator.

| | D0 threshold | D0, snapshots | Setting 0 | Setting 10 | Setting 14 |
|---|---|---:|---:|---:|---:|
| Sign agreement | at least 0.80 | 0.824 | 0.826 | 0.826 | 0.821 |
| Direction within 30° | at least 0.70 | 0.719 | 0.704 | 0.714 | 0.713 |
| Spearman | at least 0.5 | 0.692 | 0.682 | 0.678 | 0.682 |
| Top-decile ratio, continental | 0.5–2 | 0.57 | 0.51 | 0.52 | 0.52 |
| Top-decile ratio, oceanic | 0.5–2 | 0.54 | **0.40** | **0.46** | **0.41** |

The proxy agrees with the shell on evolving worlds about as well as it did at D0, except in
magnitude: it reads the most loaded oceanic cells at 0.40–0.46 of the shell value, below the 0.5
floor, and the continental ratio sits on the floor. Sign agreement drifts from 0.84 at 10 Myr to
0.80 at 200 Myr. Because the thresholds sit close to the medians, single checks often fall
below one (a single check passes sign in 63–65% of cases, direction in 52–56%, magnitude in
23–35%).

## Braking arms

Does the RFC 0008 contact brake slow continents enough to hold back rifting? Stage 1 setting 0 on
the five stage 1 worlds, 200 Myr at 1 Myr steps, plates-only, with the brake's two world tuners
(`run_rift_damage_d1.py --world-knobs`). 15 of 15 runs finished; summary in
`braking-summary.json`. Speed is the area-weighted mean of `|ω × r|` over crust of each kind,
sampled every 10 Myr. The brake-on arm reproduces the stage 1 runs split for split.

| | Brake on | `craton_drag_multiplier = 0` | `contact_drag_multiplier = 0` |
|---|---:|---:|---:|
| Continental crust speed (cm/yr) | 0.63 | 0.63 | 1.55 |
| Oceanic crust speed (cm/yr) | 0.77 | 0.77 | 2.12 |
| Share of lines above the growth threshold | 40% | 40% | 81% |
| Splits | 44 | 48 | 109 |
| Splits per world per 100 Myr (production 7.7) | 4.4 | 4.8 | 10.9 |
| Continent-majority splits | 4 | 4 | 10 |
| Oceanic share of successful rifts | 91% | 92% | 91% |
| Line length in oceanic crust | 66% | 66% | 68% |
| Splits mostly along inherited weakness | 20% | 31% | 61% |
| Weak share of split lines, against random paths | 29% / 23% | 35% / 27% | 56% / 44% |
| Transition ratio | 0.99 | 0.99 | 0.82 |
| Failed share | 78% | 78% | 89% |
| Localization violations | 0 | 0 | 8 |
| Targets met (of 10) | 8 | 8 | 6 |

- The whole brake costs a factor of about 2.5 in speed and in split rate. Its craton component
  alone changes neither.
- Plates are slow in every arm: 0.4–1.4 cm/yr braked and 1.3–2.9 cm/yr unbraked, by world,
  against about 4–5 cm/yr on Earth (#942).
- Braking does not decide where lines go. Two thirds of line length is in oceanic crust and nine
  tenths of successful rifts are oceanic in every arm, with crust about half continental by area.
  No split in any arm ran through continental crust alone.
- The inherited-weakness share rises to 61% without the brake, mostly because faster plates mark
  more rock as weak: random paths go from 23% to 44% weak. The lines' excess over random grows
  from 6 to 12 points.
- Without the brake more rifts fail and damage spreads again (8 localization violations).
