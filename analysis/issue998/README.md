# #998: model and Earth plate relative-speed distributions, rejection sampling

Exploration under RFC 0021 (tracking #957; related #995). Status: complete as of 2026-11-22. Sections 1-9 and 11-17 are
the evidence in the order it was gathered; this summary is the current answer. Section 10 is the answer as it stood
on 2026-11-21 and is kept, superseded by this summary. The one model-side change that came out of it (a mantle-flow
multiplier and the UI default below) is in the same PR; everything else is diagnostics in `bin/debug/`.

## Summary for the next agent (read this first)

**The question.** Do the model's plate relative speeds look like Earth's (PB2002), by rejection sampling, and which
driving force or flow feature predicts rejection; is it tunable in RFC 0021 or does it need a different driver?

**Answer in seven points** (all on the 20 development seeds of `analysis/issue958/development-seeds.json`, plates-only,
1 Myr steps; "plate-mean opposition" is the contact's relative speed over the sum of the two plates' area-mean speeds,
Earth 0.41 in the converging class; the node-velocity version in sections 16-17 gives Earth 0.52):

1. **t=0 acceptance is vacuous.** At t=0 speeds come from the static mantle flow alone (slab pull and ridge push are
   exactly zero until a slab exists). Every seed passes at t=0 and fails from 2 Myr on. Rejection is caused by slab pull
   and ridge push switching on at the first step (sections 3-4).
2. **With the brake off, contacts are too fast and too opposed:** median contact speed 8.4 cm/yr against Earth's 3.0,
   plate-mean opposition 0.83 against 0.41, 0/20 accepted (sections 3, 7).
3. **No single knob fixes both.** Force and viscosity corners (section 6), an overall partial brake (section 8), a
   brake extended to oceanic contacts (section 11), strong driving with a partial brake (section 12) and a stronger
   static flow alone (section 9) each move contact speed and coherence together: the brake buys coherence only by
   removing speed, driving buys speed but not coherence.
4. **Flow-dominated driving plus a partial brake does reach Earth's range.** Static flow x4 (flow about 60% of the
   driving torque) with the overall contact brake at 0.5: plate speed 3.9 cm/yr at 50 Myr (5.2 at 5 Myr), median
   contact speed 3.6, converging-class opposition 0.51, 10/20 worlds accepted at 50 Myr (section 13). On the
   node-velocity measure the model's mean opposition (0.56) and alike share (0.20) sit inside Earth's bootstrap
   interval (0.44 to 0.60 and 0.17 to 0.39, section 17). It is the only setting tested that does.
5. **Why: the flow's geometry sets neighbour coherence, and the flow places the trenches.** Contacts the flow carries
   alike end at opposition 0.34-0.40, contacts it carries oppositely at 0.71-0.81; where the flow converges contacts
   subduct 3-4 times as often as where it opens (section 16). Flipping every centre's sign moves subduction off the
   original downwellings (median percentile 0.77 to 0.37); turning the flow off leaves every plate still for 50 Myr
   (section 15). The static flow is the only starter of plate motion in this model.
6. **Tunable in RFC 0021, or a different driver?** Tunable only with a non-physical amplitude: x4 is a diagnostic, since
   flow alone then drives plates at about 18 cm/yr. The physical content of the result is that coherent, large-scale
   flow must carry a large share of the driving, and the model's flow is random, static, unrelated to its slabs, and
   has only 19-23% of contact length carried alike. A different driver is the real fix (candidates below); this
   needs an RFC revision or a new RFC.
7. **Earth is not as coherent as first assumed.** On node velocities Earth's contacts are 27% alike (under 0.3) and
   35% opposed (over 0.7); the model's flow x4 worlds are 19-20% and 33-36% (section 17). Earlier statements that no
   setting reaches Earth's coherence used the plate-mean measure (0.41); a target must name its measure.

**Settings tested at a glance** (converging-class plate-mean opposition; Earth 0.41; acceptance is the amended rule
of section 5; details in the named section):

| Setting | Plate speed (cm/yr) | Median contact speed (cm/yr) | Opposition | Accepted | Section |
|---|---|---|---|---|---|
| `b123-brake-off` (brake 0) | 4.5 | 8.4 | 0.83 | 0/20 at 50 Myr | 3, 7 |
| `b123` brake 0.25-0.35 | 2.2-2.4 | Earth-like | 0.71-0.76 | 50-60% at 100 Myr | 8 |
| oceanic-contact brake 0.5 (patch) | 2.9 | 4.2 | 0.70 | 2/20 at 50 Myr | 11 |
| slab-mantle 1e20 + brake 0.5 | 4.9 | 7.1 | 0.69 | 0/20 at 50 Myr | 12 |
| slab-mantle 1e21 + brake 0.5 | 2.8 | 4.7 | 0.75 | 0/20 at 50 Myr | 12 |
| flow x4, brake 0 | 8.7 | 11.9 | 0.70 | 0/20 at 50 Myr | 9 |
| flow x4 + brake 0.25 | 4.8 | 4.8 | 0.54 | 2/20 at 50 Myr | 13 |
| **flow x4 + brake 0.5** | **3.9** | **3.6** | **0.51** | **10/20 at 50 Myr** | 13 |

**What this PR changes in the product.** `World.mantle_flow_multiplier` (default 1.0, saves read 1.0, x1.0 bit-identical)
scales the static flow wherever it is read (`mantle.world_flow_at`). The `/worlds/{id}/controls` route now also takes
RFC 0021's force forms `ridge_push_form`, `slab_drag_form`, `age_force_column`. The UI's fresh-world defaults
(`frontend/src/uiDefaults.ts`) become: brake 0.5 (craton component 1.0), mantle flow x4, ridge push
`potential_energy`, slab drag `trench_normal`, age force column on. Backend model defaults, API and replay worlds are
unchanged. This replaces #1007's "braking 0x" default; a loaded world keeps its own stored values. The setting the
results were measured on also used RFC 0018's calibrated rift-line parameters (`analysis/issue975/settings.json`
`common_params`, `mark_power` 2.8); the UI does not push those, so UI worlds split plates with the model's
rift-line defaults and plate counts may differ from the 14-19 of these runs.

**Caveats that bound every number above.** 20 development seeds, so acceptance rates have wide intervals; the
rejection rule was amended after the first baseline passed (section 5, disclosed, one-sided results kept), so the
acceptance rates are not a held-out test: a held-out seed set has not been run. One force-form family (`b123`).
Plates-only. The x4 amplitude is unphysical. Earth is one planet's present state (106 contacts of 500 km or more).
Literature on upper-mantle flow speed is thin and model-dependent (no measured value; 20 cm/yr in a thin
asthenosphere channel in one convection model, more than 10x plate speed at a slab edge in Jadamec & Billen 2010);
no source gave the far-field speed. Pitfalls found on the way: `torque_norm_nm` omits the resisting terms (use
`torque_along_axis_nm`); the breakdown's `contact_brake` row does not see the in-process oceanic patch; opposition by
plate-mean speed and by node velocity differ (0.41 and 0.52 for Earth); scaling the flow by a constant leaves
opposition unchanged at t=0 (it is a ratio); WAIT from `compute_budget.py` can be the agent's own queue.

**Candidate next steps** (not authorized by #998; each is a model change and needs an RFC revision or a new RFC):

1. *Flow tied to slabs and downwellings.* Make downwelling strength follow existing slabs (a downwelling with a slab
   in it speeds up, the Jadamec & Billen result) and let slabs place downwellings that persist after the slab moves
   on; today the field is static and a trench drifts off its downwelling (section 14: the tendency fades 5 to 50 Myr).
2. *Coherent large-scale flow geometry.* A degree-2-like field (two broad upwellings and a downwelling girdle) in place
   of 8 random Gaussians, to raise the alike share of contact length from 19-23% toward Earth's 27%. Untested.
3. *Slab-induced suction.* A pull toward each slab on its neighbours, falling off with distance (the low-pressure zone
   behind a sinking slab). Untested; the in-process patch pattern of `plate_speed_samples.py --ocean-brake` is the
   template for a diagnostic.
4. *Prescribed slab pull*, the issue's alternative, if 1-3 do not hold up.
5. Evaluation: run the criterion on fresh seeds; compute Earth's targets on both opposition measures; use the
   calibrated rift-line parameters in the UI default if plate counts matter.

**Reproducing.** Each seed is about 1-1.5 minutes plates-only and about 1 GB (`compute-budget` first). Captures are in
`~/mantle-bloom-captures/issue998/` (not in the repo; regenerate with the commands in each script's docstring).

| Script (`bin/debug/`) | What it does | Output (`analysis/issue998/`) |
|---|---|---|
| `earth_speed_reference.py` | PB2002 per-class lognormal fits, bootstrap | `earth-reference.json` |
| `plate_speed_samples.py` | model contact samples; `--flow-scale`, `--ocean-brake` | captures only |
| `score_speed_samples.py` | the criterion (sections 2, 5) | `scores-<setting>.json` |
| `attribute_speed_rejection.py` | driving shares, resisting ratios, correlations | `attribution-<setting>.json` |
| `contact_asymmetry.py` | plate-mean asymmetry and opposition per class | `contact-asymmetry*.json` |
| `downwelling_check.py` | do trenches sit on the downwelling centres; flow off and sign-flip interventions | `downwelling-check*.json` |
| `contact_flow_coherence.py` | per-contact flow-only and actual coherence, subduction share | `contact-flow-coherence.json` |
| `earth_contact_coherence.py` | Earth's contact coherence from PB2002 poles | `earth-contact-coherence.json` |

Settings are in `analysis/issue975/settings.json` and `analysis/issue998/settings.json` (the flow-x4 settings are
`b123-brake-<0.25|0.5>-for-flow4`, run with `--flow-scale 4`).

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

## 10. Answer as of 2026-11-21 (superseded by the summary at the top; kept as written)

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

## 11. Exploration: a contact brake extended to oceanic contacts (in-process patch, not a model change)

`plate_speed_samples.py --ocean-brake G` adds, in the sampling process only, the same linear contact drag (G x 5e12 N/m)
at non-divergent oceanic band nodes (convergent, transform, contested); ridges stay free. Base: `b123-brake-off`, 20
development seeds, 50 Myr. Scores and attribution are in `scores-/attribution-b123-oceanbrake-<G>.json`, coherence in
`contact-asymmetry-oceanbrake.json`.

| G | plate speed at 50 Myr (cm/yr) | median contact speed (cm/yr) | opposition (Earth 0.41) | accepted at 50 Myr |
|---|---|---|---|---|
| 0 (brake off) | 4.25 (100 Myr) | 8.9 (100 Myr) | 0.83 | 0/20 |
| 0.1 | 4.0 | 6.8 | 0.79 | 0/20 |
| 0.25 | 3.4 | 5.8 | 0.75 | 0/20 |
| 0.5 | 2.9 | 4.2 | 0.70 | 2/20 |

Reading: the oceanic brake lowers contact speed and opposition together, as the continental brake did. It keeps a little
more plate speed at equal opposition (about 3.0 cm/yr at 0.70, against 2.3 at 0.73 for the overall 0.25 brake), but does
not reach Earth's coherence. Plate speeds are a single seed-set, and the torque-breakdown `contact_brake` row ignores the
patch (it builds its own mask), so speeds are the evidence. A real change would need an RFC revision.

## 12. Strong driving with an intermediate brake (slab-mantle viscosity 1e20 and 1e21, SUBDUCTION_LSINK 660 km)

Settings `corner-<1e20|1e21>-lsink660-brake-<0.25|0.5>` in `settings.json` (overall `contact_drag_multiplier` 0.25 or 0.5),
20 development seeds, 50 Myr. Scores, attribution and `contact-asymmetry-corner-brakes.json` are in this directory.

| Setting | peak plate speed (cm/yr) | plate speed at 50 Myr | median contact speed at 50 Myr | opposition (Earth 0.41) | accepted at 50 Myr |
|---|---|---|---|---|---|
| 1e20, brake 0.25 | 12.7 (5 Myr) | 8.1 | 9.7 | 0.72 | 0/20 |
| 1e20, brake 0.5 | 6.9 (5 Myr) | 4.9 | 7.1 | 0.69 | 0/20 |
| 1e21, brake 0.25 | 6.0 (10 Myr) | 4.2 | 6.6 | 0.80 | 0/20 |
| 1e21, brake 0.5 | 3.7 (10 Myr) | 2.8 | 4.7 | 0.74 | 0/20 |

Reading: no combination of driving strength and overall brake reaches Earth's opposition (0.41) with contact speed near 3
cm/yr. Plates overshoot early, then slow as the plate count rises from 14 to 19-21. Together with sections 9-11 this
supports the answer in section 10: the force-form family tested here cannot produce Earth's neighbour coherence by
retuning; a driver that couples adjacent plates' motion is needed (not tested).

## 13. Flow x4 with an intermediate brake (`--flow-scale 4`, settings `b123-brake-<0.25|0.5>-for-flow4`)

Static mantle-flow amplitude x4, slab pull and ridge push unchanged, overall `contact_drag_multiplier` 0.25 or 0.5.
20 development seeds, 50 Myr. Scores, attribution and `contact-asymmetry-flow4-brake.json` are in this directory.

| Brake | flow share of driving (5 / 50 Myr) | plate speed 5 / 50 Myr (cm/yr) | median contact speed at 50 Myr | opposition 5 / 50 Myr (Earth 0.41) | accepted at 50 Myr |
|---|---|---|---|---|---|
| 0.25 | 0.63 / 0.60 | 6.5 / 4.8 | 4.8 | 0.54 / 0.54 | 2/20 |
| 0.5 | 0.64 / 0.60 | 5.2 / 3.9 | 3.6 | 0.52 / 0.51 | 10/20 |

Reading: with flow carrying about 60% of the driving, the brake no longer trades coherence for speed. At brake 0.5 plate speed
is inside the 3-5 cm/yr band, the contact median is near Earth's, opposition is 0.51 (flow-only t=0 is 0.53, the unbraked
baseline 0.83) and 10/20 worlds are accepted. This is the first setting tested that approaches Earth's speed and coherence
together. Caveats: the amplitude x4 is a diagnostic (flow-only plate speed 18 cm/yr), the field is static and its centres are
unrelated to where slabs form, the amended criterion was set after seeing results, and n = 20.

## 14. Does subduction land on the model's downwelling centres? (`bin/debug/downwelling_check.py`)

Setting `b123-brake-0.5-for-flow4`, flow x4, 20 development seeds. At each time, the mean downwelling potential D (positive
over a downwelling centre) at subducting boundary nodes is compared with the same mean for 200 random redrawings of the
centres. Output: `downwelling-check.json`.

| Myr | mean D at subducting nodes | mean D, all band nodes | median percentile among random centres (0.5 = chance) | seeds above the 90th percentile (chance 10%) |
|---|---|---|---|---|
| 5 | 0.52 | 0.17 | 0.78 | 20% |
| 20 | 0.46 | 0.19 | 0.76 | 20% |
| 50 | 0.38 | 0.17 | 0.74 | 20% |

Reading: subducting nodes sit in more downwelling than chance, modestly. Nothing in the model ties trench position to the
centres, so the tendency is dynamical: the flow moves plates toward downwellings and they subduct there. It is weak and
fades from 5 to 50 Myr as plates move on and the static centres stay put. Ridges show a similar weak upwelling tendency
(median percentile 0.65-0.68).

## 15. Intervention on the flow (same 20 seeds, `b123-brake-0.5-for-flow4`)

Trench placement is always scored against the ORIGINAL centres (flow x4). Files `downwelling-check-flip.json` and `-off.json`.

| Arm | subducting nodes per seed (5 / 50 Myr) | median percentile vs random centres, at 5 / 20 / 50 Myr (0.5 = chance) | mean D at subducting nodes at 5 Myr |
|---|---|---|---|
| flow x4 (section 14) | 3,966 on seed 106533960 | 0.78 / 0.76 / 0.74 | +0.52 |
| every centre's sign flipped (-4) | 3,360 / 3,804 | 0.37 / 0.36 / 0.39 | -0.21 |
| flow off | 0 / 0 (20 of 20 seeds, to 50 Myr) | - | - |

Reading: with the flow flipped, subduction moves away from the original downwellings (onto the original upwellings), and
ridges move the other way (median ridge percentile 0.18-0.25). So in the model the flow places the trenches; the section 14
association is causal through plate motion. With the flow off nothing moves and nothing subducts in 50 Myr: the static
flow is the only starter of plate motion, because slab pull and ridge push need a slab to exist. Not a mirror image (0.37,
not 0.23): the world's geometry also sets trench positions.

## 16. Per-contact flow coherence (`bin/debug/contact_flow_coherence.py`; `contact-flow-coherence.json`)

`b123-brake-0.5-for-flow4`, 20 seeds, 5 / 20 / 50 Myr. Each contact is binned by what the flow alone does to its two plates.
Opposition here uses node velocities, not plate-mean speeds, so it is not comparable with sections 6-13.

| Bin (flow alone) | length share | subducting share (5 / 20 / 50 Myr) | actual opposition (5 / 20 / 50 Myr) |
|---|---|---|---|
| converging (closing > 0.5 cm/yr) | 0.48 / 0.52 / 0.51 | 0.65 / 0.58 / 0.52 | 0.58 / 0.56 / 0.55 |
| opening (closing < -0.5) | 0.46 / 0.42 / 0.42 | 0.15 / 0.15 / 0.19 | 0.59 / 0.58 / 0.58 |
| plates carried alike (flow opposition < 0.3) | 0.19 / 0.19 / 0.23 | 0.36 / 0.37 / 0.37 | 0.34 / 0.37 / 0.40 |
| plates carried oppositely (> 0.7) | 0.35 / 0.35 / 0.33 | 0.39 / 0.37 / 0.37 | 0.81 / 0.77 / 0.71 |

Reading: (1) where the flow converges, contacts subduct about 3-4 times as often as where it opens. (2) Actual coherence is
inherited from the flow's: contacts the flow carries alike end at opposition 0.34-0.40, contacts it carries oppositely at
0.71-0.81, and convergence alone makes no difference to opposition (0.55 against 0.58). (3) Only 19-23% of contact length is
carried alike by the model's random flow, and 33-36% oppositely; flow-only opposition overall is 0.57. (4) Subduction is not
rare at flow-alike contacts (37%), so it has other causes there. Earth's equivalent split is not yet computed.

## 17. Earth's contact coherence on the same measure (`bin/debug/earth_contact_coherence.py`; `earth-contact-coherence.json`)

PB2002 plate velocities at each boundary segment (Bird's Euler poles, no-net-rotation frame), per-contact opposition
sum(km |v_a - v_b|) / sum(km (|v_a| + |v_b|)), contacts of at least 500 km. Earth has 106 contacts; the interval is a
bootstrap over contacts. Model: `b123-brake-0.5-for-flow4`, 20 seeds (its net rotation is 0.3 cm/yr, so its frame is near NNR).

| | mean opposition | share alike (< 0.3) | share opposed (> 0.7) |
|---|---|---|---|
| Earth (PB2002) | 0.52 [0.44, 0.60] | 0.27 [0.17, 0.39] | 0.35 [0.21, 0.49] |
| model, flow x4 + brake 0.5, 5 / 20 / 50 Myr | 0.58 / 0.57 / 0.56 | 0.19 / 0.20 / 0.20 | 0.36 / 0.35 / 0.33 |
| model, flow alone (5 Myr geometry) | 0.57 | 0.19 | 0.35 |

Reading: on this measure Earth's neighbours are not much more alike than the model's flow x4 worlds: the model's mean and
shares sit inside or at the edge of Earth's sampling interval (alike share 0.20 against 0.27, interval 0.17-0.39).
Earth's plate-mean opposition (0.41, `contact-asymmetry.json`) and this local-velocity 0.52 differ, so a target must name
its measure. The unbraked baselines (plate-mean 0.83) are far outside either. The flow-alone and actual columns agree:
the model's coherence is inherited from its flow geometry.
