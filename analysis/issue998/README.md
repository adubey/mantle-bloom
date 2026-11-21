# #998: model and Earth plate relative-speed distributions, rejection sampling

Exploration, no model change. Status: in progress; sections are added as they are measured.

## 1. Earth reference (PB2002)

`bin/debug/earth_speed_reference.py` -> `earth-reference.json`: 5,819 steps, 141 plate pairs, 106 contacts of at
least 500 km. Length-weighted; uncertainty by bootstrap over plate pairs. Classes by obliquity: sliding if
|closing| < 0.5 x speed, else converging or diverging by the sign of the closing rate.

| Class | Share of length | Median cm/yr | Under 2 | Over 8 | Lognormal mu, sigma (ln cm/yr) |
|---|---|---|---|---|---|
| Converging | 31% | 4.1 | 34% | 17% | 1.08, 1.10 |
| Diverging | 34% | 2.8 | 38% | 11% | 0.91, 1.05 |
| Sliding | 35% | 2.5 | 44% | 8% | 0.80, 1.06 |
| All | 100% | 3.0 | 39% | 12% | 0.92, 1.07 |

Area dependence is weak: over 106 contacts, ln speed on ln(larger area share) and ln(smaller area share) has
R^2 = 0.06; the smaller-plate coefficient is -0.13 (95% interval -0.24 to -0.02), the larger-plate coefficient 0.00
(-0.11 to 0.12).

## 2. Rejection criterion (fixed on 2026-11-20, before any model sample was scored)

A sample is one world at one time (plates-only, from `bin/debug/plate_speed_samples.py`). Every boundary node `i` in
the contact band has a length `w_i` (km), a relative speed `v_i` and a closing rate; its class `k_i` is by the same
obliquity rule as Earth's.

- **Score** `L` = sum_i w_i [ ln pi_(k_i) + ln f_(k_i)(v_i) ] / sum_i w_i: the length-weighted mean joint
  log-likelihood of class and speed under Earth. `pi_k` is Earth's share of boundary length in class `k`;
  `f_k` is the lognormal density of `ln max(v, 0.05)` with Earth's class (mu, sigma) from section 1. (The 0.05 cm/yr
  floor is the one PB2002 uses, whose speeds are quoted to 0.1.)
- **Threshold** `L*`: the 5th percentile of `L` over 2,000 random subsets of 50 of Earth's 141 plate pairs (all their
  steps; 50 is about the contact count of a 20-plate model world, so a sample's size noise is matched). A sample is
  **rejected** if `L < L*`.
- Earth's own scores come from a distribution fitted to itself, so about 5% of Earth subsets fall under `L*` by
  construction; that is the false-rejection rate of the test.
- **Reported, not used to reject:** the score split into its class-share and speed parts; per-class share under 2
  and over 8 cm/yr; the median contact speed against Earth's area regression (section 1); which plate-size pairs
  (larger and smaller area share, in terciles) carry the contacts with the lowest per-contact log-likelihood.
- Changing any number above after seeing results needs a dated note here saying why, and the old result stays.

## 3. Sampled worlds

Pending: 20 development seeds (`analysis/issue958/development-seeds.json`), `b123-brake-off` (the #995 setting) at
0, 1, 2, 5, 10, 20 and 50 Myr; `rifting-baseline` as the slow reference. Saves and raw samples are in
`~/mantle-bloom-captures/issue998/`, not in the repo.

First look (one seed, 591453549, `b123-brake-off`, not scored): median contact speed 1.6 cm/yr at t=0, 3.9 at
1 Myr, 7.9 at 5 Myr; `rifting-baseline` stays near 0.3 throughout the first 2 Myr.
