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
