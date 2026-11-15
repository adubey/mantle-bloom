# Issue 888: RFC 0018 D1, rift-line prototype and its three checks

Interim report for [#888](https://github.com/adubey/mantle-bloom/issues/888). The prototype is the
rift-line version of RFC 0018 (PR #897 revision) in `bin/debug/rift_line.py`,
`rift_damage_mechanism.py` and `rift_damage_step.py` (draft PR #890). No parameter search has been
run: the RFC requires the three checks below to pass first. Development seeds only.

| Check | Command | Result |
|---|---|---|
| Propagation | `check_rift_line_d1.py --only propagation` | **Pass.** Steady tip speed within 20% of the discrete law at all 6 settings (default and half density, tip length 60/125/250 km). Within 2% at `Δε_max` 0.005 |
| Integration | `check_rift_line_d1.py --only integration` | **Fails as specified** (strain error 0.40–0.82 against 0.05). Passes with three amendments, see below |
| Line (20 worlds, 0 and 50 Myr) | `check_rift_lines_world_d1.py` | **Fails the validity part:** 66% of eligible plates give a valid split (need 90%). Craton part passes (89% of plates with craton) |

Compact results: `integration-propagation.json`, `line-check-summary.json`. Per-plate rows and
checksums are in the private captures directory.

## Cell spacing

The RFC says the default cell spacing is 62.5 km. It is 125 km (`line_spacing_rad(1.0)`), equal to
`w_0`. Consequences at default density: the band is one cell (the line); the non-local softening
disk (radius `w_0/2`) holds only the cell itself; the tip mesh factor is 0.58 (`x = 1` at
`ℓ_tip = 125 km`), not 0.77. With `ℓ_tip = 125 km` the steady speed is 58 km/Myr, not 77, and two
tips cross a 5,000 km line in about 45 Myr, not 34. Half-density tip speed relative to default:
0.58, 0.83, 0.94 for `ℓ_tip` 60, 125, 250 km.

## Integration check

Manufactured plate, fixed pull, straight line, steps 0.01 (reference), 0.5, 1, 2, 5 Myr. The RFC
integrator (frozen `N_cap` per step, sub-steps of `Δε_max = 0.05`) is first-order in the
sub-step, because a line cell reaching `ε_cand` changes its neighbours' rate only at the next
sub-step, and `N_cap` follows thinning only at the next step. Worst strain error by variant at
the starting rate:

| Variant | 0.5 Myr | 1 | 2 | 5 |
|---|---:|---:|---:|---:|
| RFC as written | 0.40 | 0.44 | 0.55 | 0.82 |
| + stop sub-steps at `ε_cand` crossings | 0.11 | 0.17 | 0.26 | 0.42 |
| + `N_cap` of line cells follows thinning in the sub-steps | 0.079 | 0.084 | 0.088 | 0.089 |
| + `Δε_max = 0.02` | 0.035 | 0.036 | 0.038 | 0.038 |

The last variant passes every criterion (strain within 0.05, broken sets within one cell per tip,
split time within one coarse step) with at most 37 sub-steps per 1 Myr step at the starting rate
(139 at 5 Myr). The remaining error depends on `Δε_max`, not on the step.

## Line check

At 0 Myr (277 eligible plates) and 50 Myr (387), production splitter, 1 Myr steps; 0.5 Myr steps on
three seeds give the same picture. Of 664 plates, 83 have no line (no positive pull peak), 439 give
a valid split. Rejections among the 581 lines: 89 daughter disconnected after the line cells vote,
51 stray piece above 1%, 1 daughter under 15%, 1 no separation. Valid lines split off a median 33%
of the plate (quartiles 25–42%), are 1.30 times the straight distance between their ends, and run
65% over oceanic crust; median `Λ` is 0.08, against `Λ_grow = 0.3`.
