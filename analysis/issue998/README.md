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

## 5. Amendment 1 (2026-11-21, made after scoring `rifting-baseline`; the one-sided result above stands)

`rifting-baseline` (brake on, old force forms, median contact speed 0.5 cm/yr, 39% of Earth's contacts are under
2 cm/yr and nearly all of its are) passed the section 2 test 20/20 with a median L of -2.85, better than Earth's
own -3.50. A mean log-likelihood rewards sitting at the density's mode, so it cannot say that a world is too
slow. The criterion as fixed above is therefore one-sided and incomplete; it is kept and reported as
"one-sided v1". Amended rule, both parts again from Earth's 2,000 random 50-pair subsets:

- reject if L is outside Earth's central 90% (below L* = -3.756 or above L_upper = -3.247), or
- reject if the length-weighted Kolmogorov-Smirnov distance of the sample's contact speeds from Earth's pooled
  contact speeds exceeds the 95th percentile of the same distance for the subsets (KS* = 0.199; the subsets are
  part of the pooled reference, so this threshold is a little low and the test a little strict).

By construction a true Earth subset is rejected up to about 15% of the time (three 5% tails, which overlap).

## 6. Other settings (same 20 seeds and 7 times)

Accepted worlds out of 20, amended rule (one-sided v1 in brackets), median contact speed in cm/yr at 10 Myr:

| Setting | 0 Myr | 1 Myr | 2-50 Myr | Median contact speed at 10 Myr | Plate speed (area-weighted) at 10 Myr |
|---|---|---|---|---|---|
| `rifting-baseline` (main) | 0 (20) | 0 (20) | 0 (20) | 0.5 | 0.75 |
| `b123-brake-off` | 2 (20) | 2 (2) | 0 (0) | 8.2 | 5.2 |
| 1e21 slab mantle, 150 km | 2 (20) | 0 (0) | 0 (0) | 11.2 | 6.9 |
| 1e21, 660 km | 2 (20) | 0 (0) | 0 (0) | 14.1 | 8.7 |
| 1e20, 150 km | 2 (20) | 0 (0) | 0 (0) | 19.9 | 11.2 |
| 1e20, 660 km | 2 (20) | 0 (0) | 0 (0) | 35.1 | 20.5 |

No setting is accepted at any time after t=0 under the amended rule (the 2 at t=0 are worlds whose flow-only
state falls inside the band by chance). Corner viscosity and sink length move plate speeds 4-fold and contact speeds
4-fold, never to Earth's.

Ratio of the median contact speed to the area-weighted plate speed (median contact speeds here are per contact, not
per node, which differs from Earth's 3.0 node-weighted by about 0.1 cm/yr):

| Setting | 5 Myr | 10 | 20 | 50 |
|---|---|---|---|---|
| `rifting-baseline` | 0.66 | 0.67 | 0.68 | 0.68 |
| `b123-brake-off` | 1.55 | 1.58 | 1.63 | 1.86 |
| 1e21, 150 km | 1.56 | 1.63 | 1.68 | 1.92 |
| 1e21, 660 km | 1.58 | 1.62 | 1.68 | 1.71 |
| 1e20, 150 km | 1.79 | 1.79 | 1.92 | 2.14 |
| 1e20, 660 km | 1.69 | 1.71 | 1.75 | 1.63 |
| Earth (PB2002, no net rotation) | 0.8 | | | |

The ratio is 1.55 to 2.1 in all five slab-pull-driven settings across a 4-fold range of plate speed. The
mantle-flow-driven baseline has 0.66, near Earth's 0.8. So slab-pull driving in the model puts both sides of a contact
into fast motion regardless of the force scale, and rescaling the force constants cannot reach Earth's ratio.

## 7. Are neighbouring plates coherent? (`contact_asymmetry.py` -> `contact-asymmetry.json`)

Hypothesis from the contact-to-plate speed ratio (section 6): at the model's converging contacts both plates are
fast, where Earth has a fast slab-pulled plate against a nearly still overriding plate. **Not supported.** The
plate-mean speed asymmetry min/max is alike on Earth and in the model (converging contacts: median 0.66 on Earth,
0.56 in `b123-brake-off` at 100 Myr; one-sided, under 0.25, in 10% of contact length in both).

What differs is direction. `opposition` = the contact's relative speed over the sum of the two plates' speeds:
about 1 when the plates move in opposite directions across the contact, near 0 when they move alike. Length-weighted
over contacts of 500 km or more, converging class (diverging and sliding give the same within 0.03):

| | Opposition, median | Contact length with opposition over 0.8 |
|---|---|---|
| Earth (PB2002, 75 contacts) | 0.41 | 22% |
| `rifting-baseline` at 50 Myr (flow-driven, slow) | 0.42 | 22% |
| `b123` brake on, 100 Myr | 0.56 | 31% |
| `b123-brake-off` at 10 Myr | 0.82 | 56% |
| `b123-brake-off` at 100 Myr | 0.84 | 64% |

Neighbouring plates in the slab-pull-driven, brake-off model move nearly against each other; on Earth, and in the
flow-driven baseline, they move mostly alike. The contact brake pulls the model part of the way back. This is a
plate-mean-speed approximation (it ignores the variation of speed along a plate and uses each plate's area-mean
no-net-rotation speed), and the opposition is nearly the same for converging, diverging and sliding contacts in both
Earth and the model, so it measures how coherent the plates' motions are, not how they relate to the boundary type.

## 8. Partial contact brake (`b123`, overall `contact_drag_multiplier` 0.25 and 0.35, craton 1.0; 20 seeds, 100 Myr)

Chosen from #995's one-seed sweep (section 9 there: 0.25 matched Earth's median contact speed), with 0.35 added
inside the sharp region between 0.25 and 0.5. Amended rule; one-sided v1 in brackets.

| Setting | Myr | Accepted | Median contact speed, cm/yr | Plate speed (area-weighted), cm/yr | Rejected by L low / L high / KS |
|---|---|---|---|---|---|
| brake 0.25 | 20 | 4/20 (8) | 4.7 | | |
| | 50 | 8/20 (14) | 4.1 | | |
| | 75 | 8/20 (16) | 3.6 | 2.47 | 4 / 1 / 12 |
| | 100 | 10/20 (16) | 3.6 | 2.36 | 4 / 1 / 9 |
| brake 0.35 | 20 | 8/20 (13) | 3.8 | | |
| | 50 | 9/20 (18) | 3.8 | | |
| | 75 | 9/20 (19) | 3.4 | 2.23 | 1 / 3 / 10 |
| | 100 | 12/20 (19) | 3.1 | 2.22 | 1 / 3 / 7 |

(The three counts are not exclusive: a sample can fail more than one.) These are the first settings with substantial
acceptance, 50% and 60% of worlds at 100 Myr; a world drawn from Earth passes about 85% of the time by construction
(and 10 or 12 of 20 has a wide interval, roughly 27-73% and 36-81%). The speed distribution is close to Earth's;
the remaining rejections are mostly the shape of the distribution (KS), not its level. Two costs:

- **Plate speed 2.2-2.4 cm/yr**, below T1's 3-5 cm/yr band (the worlds range 1.2-3.5).
- **Coherence is only partly restored**: opposition at converging contacts is 0.76 (brake 0.25) and 0.71 (0.35) at
  100 Myr against Earth's 0.41, brake off 0.84, brake on 0.56. Half the contact length still has the neighbours moving
  nearly against each other (opposition over 0.8: 0.50 and 0.48, Earth 0.22).

## 9. Static mantle-flow amplitude (`b123-brake-off`, centre strengths x2 and x4 after generation; 20 seeds)

| Flow scale | Myr | Accepted | Median contact speed, cm/yr | Plate speed, cm/yr | Flow share of driving | Opposition (converging) |
|---|---|---|---|---|---|---|
| x1 | 50 | 0/20 | 8.4 | 4.5 | 0.29 | 0.86 |
| x2 | 50 | 0/20 | 9.6 | 5.8 | 0.43 | 0.79 |
| x4 | 50 | 0/20 | 11.9 | 8.7 | 0.58 | 0.70 |

Stronger flow raises plate and contact speeds and moves the opposition only a little (0.86 to 0.70); it does not
make the contacts Earth-like. Mantle flow in the model is not the missing ingredient for contact speeds.

## 10. Answer

**Which force or flow feature predicts rejection?** Slab pull and ridge push switching on at the first step
(sections 3-4). Speeds are mantle-flow-only at t=0 (so t=0 acceptance is vacuous), and from the first step slab
pull is 52-58% of the driving torque and plate speed follows it per plate (Spearman 0.75-0.80). Across worlds the
pulled share of boundary and the opening share separate outcomes at 1 Myr; no torque share does later.

**Is it tunable in RFC 0021, or does it need a different driver?** Both parts of the answer matter:

1. *Contact-speed level: tunable, with a cost.* The existing contact brake multiplier at 0.25-0.35 (the multiplier
   already in RFC 0008/0021's reach) brings the contact-speed distribution to 50-60% acceptance at 100 Myr, where the
   brake at 0 gives 0% and at 1 gives about 5% (too slow). No force or viscosity corner reaches it (section 6), and
   neither does a stronger static mantle flow (section 9).
2. *Plate speed and coherence: not by these knobs.* The same setting leaves plate speed at 2.2-2.4 cm/yr (below
   T1) and neighbouring plates moving against each other (opposition 0.71-0.76 against Earth's 0.41), which the
   flow-driven baseline reproduces (0.38-0.42) but only at 0.4-0.8 cm/yr plate speed. Earth-like speed and
   Earth-like coherence together have not been reached by any setting tested. The structural candidate, not tested
   here, is that slab pull should act through a field that also drags the neighbouring plates (slab-induced flow),
   or be prescribed (the issue's alternative); this needs an RFC revision.

Limits: 20 development seeds, one force-form family (`b123`), 1 Myr steps, plates-only; acceptance rates have wide
intervals at n=20; the Earth reference is one planet's present state (141 plate pairs), and the amended rejection
rule was fixed after seeing the baseline pass (section 5).
