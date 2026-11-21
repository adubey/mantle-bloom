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
  `f_k` is the lognormal density in `v` (cm/yr) of `max(v, 0.05)`, Earth's class (mu, sigma) from section 1. (The 0.05 cm/yr
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

## 3. Sampled worlds and acceptance

20 development seeds (`analysis/issue958/development-seeds.json`), `b123-brake-off` (the #995 setting), plates-only,
1 Myr steps, sampled at 0, 1, 2, 5, 10, 20 and 50 Myr (`bin/debug/plate_speed_samples.py`; raw samples in
`~/mantle-bloom-captures/issue998/`, not in the repo). Scores: `scores-b123-brake-off.json` (`score_speed_samples.py`).
Threshold L* = -3.76 (Earth median -3.50).

| Myr | Accepted | Median L | Median contact speed, cm/yr |
|---|---|---|---|
| 0 | 20/20 | -3.12 | 2.3 |
| 1 | 2/20 | -4.16 | 6.1 |
| 2 | 0/20 | -4.55 | 7.3 |
| 5 | 0/20 | -4.80 | 8.1 |
| 10 | 0/20 | -4.80 | 8.2 |
| 20 | 0/20 | -4.69 | 8.2 |
| 50 | 0/20 | -4.76 | 8.4 |

Rejection is not a property of particular speed assignments: every seed passes at t=0 and fails from 2 Myr on. Every
(larger-plate, smaller-plate) area tercile of contacts is below L* too; large-large contacts score best
(-4.4 at 50 Myr) and small-small worst (-5.5), so size sorts the severity but no size pair is Earth-like.

## 4. Attribution (`attribute_speed_rejection.py` -> `attribution-b123-brake-off.json`)

Means over the 20 worlds. Shares are of the summed norms of the three driving terms; "/ driving" is the resisting
torque along the plates' motion axes over the driving torque along them.

| Myr | Slab pull share | Ridge push | Mantle flow | Pulled share of boundary | Basal drag / driving | Slab drag / driving | Area-weighted speed | Flow-only speed |
|---|---|---|---|---|---|---|---|---|
| 0 | 0.00 | 0.00 | 1.00 | 0.00 | 1.00 | 0.00 | 2.44 | 2.44 |
| 1 | 0.52 | 0.10 | 0.38 | 0.18 | 0.74 | 0.26 | 4.36 | 2.44 |
| 5 | 0.58 | 0.14 | 0.28 | 0.28 | 0.66 | 0.34 | 5.21 | 2.44 |
| 20 | 0.55 | 0.18 | 0.27 | 0.27 | 0.64 | 0.34 | 5.01 | 2.43 |
| 50 | 0.51 | 0.20 | 0.29 | 0.25 | 0.63 | 0.35 | 4.51 | 2.32 |

1. **The accepted state is not an Earth-like speed assignment.** At t=0 slab pull and ridge push are exactly zero in
   all 20 worlds and the plates move on mantle flow alone (2.44 cm/yr). The t=0 acceptance is therefore vacuous: the
   first step turns slab pull and ridge push on (paired change to 5 Myr: slab-pull share +0.58, 95% interval
   0.55 to 0.61; pulled boundary +0.28; speed +2.8 cm/yr, 2.5 to 3.0; all 20 seeds in the same direction). Why
   the t=0 solve has no boundary forces is not verified here; the boundary classification reads the previous
   step's rotation rates, which are not yet set at t=0 (`torque.classify_boundary_nodes`).
2. **Per plate, speed follows slab pull and ridge push.** Spearman rho of plate speed with slab pull per unit area is
   0.75 to 0.80 at every time (280-300 plate-samples each); with ridge push per area 0.35 at 1 Myr rising to 0.73 at
   50 Myr; with mantle flow per area 0.70 falling to 0.39. Speed against area is -0.14 at 1 Myr and -0.57 at 50
   (large plates slower; Earth's contact speeds are nearly independent of area, section 1).
   Correlations with the drag terms (-0.8 to -0.96) are responses to speed, not causes.
3. **Across seeds, at 1 Myr (the only time with both outcomes, 2 accepted, 18 rejected), the standardized
   differences (rejected minus accepted) are:** pulled boundary share d = 2.2, opening-past-gate share 1.4, plate
   count -1.1, largest-plate share 1.2, slab-bearing area share -0.5. With 20 worlds, 15 features and no multiplicity
   correction these are indicative only. At 5 to 50 Myr all are rejected, so only rank correlations with L remain:
   area-weighted speed (-0.68 to -0.75, nearly the definition), opening-past-gate share (-0.54 at 5 Myr), slab-attached
   area share (-0.47), largest-plate share (+0.57). Slab-pull share of driving itself does not separate worlds
   (rho -0.05 at 5 Myr, +0.41 at 20 Myr).
4. **Contacts are faster than the plates' speeds imply.** The median contact relative speed over the area-weighted plate
   speed is 1.6 to 1.9 in the model (8.1/5.2 at 5 Myr, 8.4/4.5 at 50) and 0.8 on Earth (3.0 median contact against
   3.76 cm/yr mean plate speed, no net rotation). Rescaling every speed to bring contacts to Earth's 3 cm/yr would
   put plates near 1.9 cm/yr, half of T1's target, so no single force scale meets both.

## 5. Answer (provisional)

- **Which feature predicts rejection?** The forces that appear at the first step, slab pull (52 to 58% of driving
  torque) and ridge push (10 to 20%), with plate speed tracking slab pull per area. Across worlds the pulled share of
  boundary and the opening share separate the two outcomes at 1 Myr; size and layout (plate count, largest plate)
  rank next. No driving-torque share separates worlds once all are rejected.
- **Tunable in RFC 0021 or a different driver?** The evidence points to a driver or structure problem more than a
  scale problem, by item 4: a global force or viscosity scale moves plate speeds and contact speeds together, and
  the two Earth targets sit at a different ratio (0.8) from the model's (1.6 to 1.9). That ratio would fall if
  neighbours moved more alike, as slab-pulled plates and their overriding neighbours do on Earth, which the model's
  both-plates-fast contacts do not. This is an inference from the ratio, **not tested**. The direct test is cheap:
  score the RFC 0021 corner settings and `rifting-baseline` with the same pipeline (about 35 minutes per setting). If
  no corner reaches L* at any time, tuning cannot fix it and a driver change is indicated (prescribed slab pull).
  Pending, with the `rifting-baseline` reference.
