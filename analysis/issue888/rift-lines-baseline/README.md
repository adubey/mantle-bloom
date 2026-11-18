# Issue 888: D1's strength-map arms on the production rift-line method

RFC 0018's D1 search ran on the prototype in PR #890. The method that PR #950 puts in the model (`World.split_method =
"rift_lines"`, for #846) differs from it after two rounds of review:

- plume heat is kept at the split test, where the prototype dropped it;
- flank transfer respects the 84 km crust cap and moves crust only to cells of the same crust type;
- a donor whose recipients are all capped doesn't strain;
- the damage pass computes every plate before writing any.

So D1's tables don't describe the production method. This reruns D1's brake-off strength-map arms on it, as the
baseline for tuning under RFC 0021 (#942).

## Runs

- **Model:** PR #950 at `6bda6906`.
- **Scripts:** `bin/debug/run_rift_lines_arms.py` and `score_rift_lines.py` at `bd7b2c3b`.
- **Arms:** `arms.json`. Stage 1 setting 0 (`eps_dot0` 4.6e-8, `s_warm` 0.159, `n_plume` 6.2e11, `lam_grow` 0.544,
  `l_tip_m` 143 km) with ocean-floor strength from age, `k_cap` 0.14 and the contact brake off
  (`contact_drag_multiplier = 0`). Then:
  - craton root: `craton_root_m = 80e3`;
  - deeper marks: `mark_power = 2.8`;
  - route: `route_mark_power = 2.8`.

  The great-circle arm is the production splitter at its defaults (brake on), for the "splits happen" rate.
- **Worlds and stepping:** the five stage 1 worlds (106533960, 1097624328, 1397817354, 1629738962, 591453549),
  generated at node density 1, stepped plates-only at 1 Myr for 200 Myr.
- **Run status:** 25 of 25 runs finished.
- **Outputs:** `~/mantle-bloom-captures/issue888/rift-lines-baseline/<arm>/<seed>/`, with `SHA256SUMS` beside them.
  The scores are in `scores.json`.

The great-circle arm reproduces D1's production baseline world for world (17, 19, 20, 14 and 7 splits). So the rate
the rift-line arms are judged against is unchanged at 7.7 splits per world per 100 Myr, which puts the bar at 3.85.

## Scores

D1's prototype values for the same arm and worlds are in brackets.

| Brake off, `k_cap` 0.14 | Base | Craton root | Deeper marks | Route |
|---|---:|---:|---:|---:|
| Targets met (of 10) | 5 [5] | 4 [5] | **6** [7] | **6** [6] |
| Splits per world per 100 Myr (need 3.85) | 4.5 [4.8] | 3.5 [2.7] | 7.8 [9.4] | 3.1 [3.9] |
| Continent-majority splits | 16 of 45 [19 of 48] | 11 of 35 [5 of 27] | 25 of 78 [28 of 94] | 10 of 31 [13 of 39] |
| Splits mostly on weak rock (need 67%) | 53% [58%] | 49% [56%] | **79%** [82%] | **71%** [67%] |
| Weak share of split lines / random paths | 56% / 46% | 55% / 43% | 68% / 46% | 64% / 45% |
| Completed-rift median (Myr, 10–50) | 63 [68] | 70 [36] | 64 [46] | 52 [90] |
| Failed share (20–80%) | 93% [88%] | 90% [94%] | 88% [87%] | 93% [92%] |
| Transition ratio (need 3) | 2.75 [1.02] | 1.62 [1.53] | 1.50 [1.08] | 1.04 [1.33] |
| Localization violations | 2 [2] | 1 [1] | 9 [18] | **0** [2] |
| Split paths with plume heat | 29% [25%] | 23% [22%] | 35% [39%] | 29% [28%] |

Every arm meets the onset width (125–133 km), width growth, no-microplate (smallest daughter 17–20% of the parent)
and line-shape targets. Interior thinning closes every column's ledger to 8e-16 or better. Rejected
splits are almost all `stray_too_large` (1–72 per arm), and no step needed more than 30 sub-steps.

## Reading

- **No arm gets past 6 of 10.** Deeper marks lost the duration target (64 Myr against 46 in D1). Route marks now miss
  it narrowly (52 Myr against 90), but split too rarely (3.1 against the 3.85 bar).
- **Deeper marks still do the most:** the weakness target passes (79%), split lines cross 22 points more weak rock than
  random paths, and the split rate is the highest. Localization violations halved from D1 (9 against 18) but still fail.
- **Route marks are the cleanest split shapes:** no localization violations, 71% of splits on weak rock, and
  completed rifts take 52 Myr, just above the 50 Myr limit. They now split too rarely.
- **The base arm's transition ratio rose to 2.75,** the closest any arm has come to 3. D1's was 1.02. The other arms
  moved less. Which of the method changes did this hasn't been checked.
- **Failing on every arm:** the failed share (88–93% against at most 80%), the completed-rift duration (52–70 Myr
  against at most 50) and the transition ratio (1.0–2.75 against at least 3).
- **The share of split paths with plume heat barely moved** (23–35%), so fixing the plume read did not by itself
  change how often plumes sit on the lines that split.
- **Below-sea-level rift basins wider than 375 km** (the descriptive super-sea check) are common: 1,209–2,844 counted
  at the 5 Myr checks, as in D1 (1,810 on its base arm).

For RFC 0021, deeper marks and route marks are the two starting points. They pass the weakness target in different
ways, and between them they fail on rate, duration, localization and the failed share.

## Reproduce

```bash
backend/.venv/bin/python bin/debug/run_rift_lines_arms.py --arms analysis/issue888/rift-lines-baseline/arms.json \
    --out ~/mantle-bloom-captures/issue888/rift-lines-baseline --jobs 3
backend/.venv/bin/python bin/debug/score_rift_lines.py --runs ~/mantle-bloom-captures/issue888/rift-lines-baseline \
    --out analysis/issue888/rift-lines-baseline/scores.json
```

Each plates-only world took 5–6 minutes and used about 0.5 GB.
