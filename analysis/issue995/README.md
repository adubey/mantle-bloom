# #995: RFC 0021's best-case check failed; alternatives D and C

## Summary and pause (2026-11-22)

RFC 0021 is paused and needs revision. This directory holds what #995 found; no model change ships from it.

**RFC 0021's own question (G1).** The best-case check failed T2, and the rescore shows why: T2's "other" group is not
slab-free (99-100% of its area has a subducting node, with slab pull torque as large as the attached group's), so no
force constant can separate the groups (section 1). D cannot move T2 by construction and C has no measured defect
(section 2). Two probe knobs (age-dependent force column, slab-history gate; sections 5-6) did not fix it. Earth's contrast is a
few plates with long trenches and most plates with almost none; the model's plates all converge somewhere at similar rates
(section 7). T2 needs a definition that measures slab drive, and G2 stays blocked until the RFC is revised.

**Relative speeds and slowdown (sections 8-10).** The contact brake, not slab pull, sets the model's neighbour relative
speeds (median 7.2 cm/yr at brake 0 against 3.2 on Earth; 3.2-3.9 at overall brake 0.25); the craton component has
little effect. Mean plate speed falls over time in most seeds and tracks the largest plate's growing share, not plate
count.

**Rifting (sections 11-20).** Splits are frequent (14-16 per 100 Myr) and about 20% of big-plate samples (10% of the
planet or more) break up within 20 Myr (about 1% per Myr), which is at the fast end of Earth's supercontinent durations
(0.3-1% per Myr), not below it. The hazard is flat in time spent big, where Pangea held about 75-125 Myr; early splits
are not a merger artifact. Craton-rich big plates split about 3 times less, independent of size. Persistent softening
does not raise the big-plate rate, `k_cap` 0.07 doubles all splits, and a time-dependent strength at the break-up test is
inert because the binding gates are the broken share and `split_masks`. Any time dependence has to act on damage
accumulation (the strain-rate probe was designed and not run).

**Not done.** Review fixes to `g1_world_checks.py` (wet-cohort depth table, per-correction ledger audit, signed Hm
residual, `--hm-ledger` docstring); the baseline Hm gap; the T1 NNR mean against the `MAX_PLATE_RATE` clamp; the RFC
revision itself; the strain-rate probe.

**Data.** Raw per-run outputs (`rift-state-*`, `rift-merge-*`, `slowdown-*`, `splits-*`, 110 files, 27 MB) are not in the
repo. They are in `~/mantle-bloom-captures/issue995/runs/` with a `SHA256SUMS`; the chain scripts here regenerate them
(7-10 min per run). The analysis scripts take `--dir` for that directory. Summaries (`*-summary.json`, `breakup-stat.json`,
`ramp-breakup-stat.json`, `craton-regression.json`, `rift-state-analysis.json`, `splits-analysis.json`,
`slowdown-analysis.json`) are checked in.

**Probe-only code (removed from this PR).** Sections 5 and 6 used two default-off knobs, `torque.FORCE_COLUMN_OWN_HM_CAP_M`
(a cap on the own-Hm floor of B.3's force column) and `torque.SLAB_HISTORY_GATE_SPACINGS` (pull only near floor the plate
consumed in the last 10 Myr), with unit tests. The owner dropped them from PR #1021, so the settings `b123-age-pull` and
`b123-slab-history` in `settings.json` no longer run on this tree; the code is in commit `8024d64a` of the PR's branch
history. The debug scripts in `bin/debug/` are the diagnostics for these sections.

Rescoring of the #980 captures (`~/mantle-bloom-captures/issue980/`, checksums verified), no new runs.
Script: `bin/debug/g1_t2_groups.py`; output `t2-groups.json`. T2 pairing is `calibrate_plate_forces.t2_ratio`'s
(speed 10 Myr after the sample; slab-attached = floor consumed in 10 Myr >= 10% of area).

## 1. Why the slab-free plates speed up with the slab-attached ones

T2's "other" group is not slab-free. 99-100% of its area (89% in the 10^20/660 km corner) has a subducting node, and
its slab pull torque is as large as the attached group's:

| Setting | Attached: share of area, speed | Other: share, speed | Mean plate area attached / other, Mkm2 | Slab pull torque, 1e25 N m, attached / other | Slab km per Mkm2 attached / other |
|---|---|---|---|---|---|
| b123 defaults | 57%, 5.16 | 43%, 3.24 | 21.8 / 60.0 | 17.3 / 22.2 | 305 / 164 |
| 1e20, 150 km | 75%, 10.92 | 25%, 5.94 | 34.3 / 47.9 | 17.8 / 17.7 | 308 / 157 |
| 1e20, 660 km | 94%, 23.35 | 6%, 17.90 | 24.6 / 8.5 | 71.2 / 26.8 | 422 / 319 |
| 1e21, 150 km | 66%, 6.68 | 34%, 4.38 | 32.9 / 72.5 | 14.4 / 17.1 | 282 / 142 |
| 1e21, 660 km | 81%, 8.34 | 19%, 6.84 | 18.2 / 30.5 | 56.5 / 49.0 | 385 / 224 |

(Plate-samples at 10^20/660 km: 1,846 attached, 345 other. Codex's counts were right.)

- Speed follows slab pull per unit area, not T2's grouping variable. Across plate-samples, correlation of ln(speed)
  with consumed share is -0.08 to 0.10; with pull torque per area it is 0.11 to 0.47; with ridge push per area 0.37 to 0.66.
- The "other" group differs from "attached" mainly in area (at 10^20/150 km: 48 against 34 Mkm2, equal pull, so
  pull per area is 0.37 against 0.52 of the attached). At 10^20/660 km the "other" plates are small (8.5 Mkm2)
  and slab-bearing, so their pull per area is as high as the attached group's and the ratio falls to 1.3.
- Strong pull compresses contrast. Ratio of speeds, top third over bottom third of plate-samples by slab pull per
  area: 2.63 (b123 defaults), 1.95, 1.72, 1.58, 1.27 (the four corners). Stronger pull makes more plates reach the
  slab-drag-limited speed, so they converge.
- Mantle flow is not what lifts the slow plates: flow-only speed is 2.2-2.6 cm/yr in every setting, and flow torque
  is 0.2-0.7 against 27-71 (x 1e25 N m) slab pull at 10^20/660 km.

Reading: T2 as defined cannot separate plates that have a slab from plates that don't, because in this model nearly
every plate has a trench with similar pull. No force constant fixes that.

## 2. D and C against T2

- **D cannot move T2.** B.3 already gives the forces the age column (`force_column_hm`), and D changes the plate's own
  Hm, so the forces that set speed are the same with D on. D is a bathymetry, sea-level and ledger question (G1's
  static estimate: land 30.2% to 32.4%, sea level -570 m, oceanic Hm 33 to 62 km). Ruled out as a fix for T2 by that
  argument (not run).
- **C is not indicated for T2 by a named number.** Flow torque is 1% of slab pull at the 660 km corners and 10-25% at the 150 km corners (flow-only speed 1.5-2.6 cm/yr throughout), and T7 (flow
  draw decides speed) passes at 1.5 against 3. C has no measured defect here to remove.

## 3. The Hm ledger residual in context

The gap is (live Hm inventory - start inventory - booked sources and sinks) / start inventory, worst over samples.
Five worlds, 200 Myr, final values: 6.3, 1.0, 2.6, 0.1, 4.7 (x 1e-4).

- It is a one-off step, not a drift. Each world's gap stays flat for tens of Myr and jumps at a few events (world
  106533960: 0 to 2.7e-4 at 20 Myr, to 6.2e-4 at 90 Myr, flat otherwise). Something writes Hm outside the booked
  accounts at those events. Not identified.
- Size, with an estimate: if the planet-mean Hm is about 40 km (not measured; the inventory isn't saved), 6.3e-4 is
  about 25 m of thickness averaged over the planet. Alternative D would change the oceanic Hm by +87% (33 to 62 km),
  about 1,000 times larger.
- What we can't say yet: whether the corrections cause it. The baseline's gap hasn't been measured. If the
  uncorrected model shows similar steps, it predates RFC 0021 and G2 only needs to show the arms don't increase it.

## 4. Why every plate has a slab: boundary geometry and slab age

Not initialization: the median plate's pulled-trench share of its boundary is 0.43-0.46 at every sample from 10 to
200 Myr (defaults and 10^20/150 km), through about 118 splits and a plate count that grows from 17 to 26.

`bin/debug/boundary_geometry_probe.py`, one world (591453549), defaults (`b123-brake-off`), samples at 50, 100, 150
and 200 Myr; `boundary-geometry-591453549.json`. Length-weighted over boundary-band nodes:

| Myr | Tangential speed > normal | Closing | Closing past 1 cm/yr (what pulls) | \|Normal\| < 1 cm/yr (the "transform" class) | Median relative speed, cm/yr |
|---|---|---|---|---|---|
| 50 | 49% | 48% | 43% | 9.9% | 8.1 |
| 100 | 50% | 50% | 41% | 17% | 7.4 |
| 150 | 48% | 47% | 43% | 9.4% | 8.8 |
| 200 | 47% | 49% | 44% | 9.1% | 8.0 |

- Half the boundary slides more than it closes, but `classify_boundary_nodes` calls a node a transform only when its
  normal rate is under 1 cm/yr, an absolute threshold, against a median relative speed of 8 cm/yr. So sliding
  boundary is labelled closing or opening by the sign of a small normal component, and 41-43% of boundary length
  is pulled as trench.
- Caveat: 50% is also what a boundary orientation independent of the relative motion gives. `direction_to_neighbor`
  is toward the nearest neighbour node, a noisy stand-in for the boundary normal. The finding is that the model's
  boundary orientation is uncorrelated with motion; Earth's isn't (ridges and trenches form roughly square to it).
  One world, one setting, four samples.

Age of the floor at the pulled-trench nodes (share of pulled length):

| Myr | 0-10 | 10-30 | 30-60 | 60-100 | over 100 | Force-column Hm over own Hm by more than 5% | Mean force Hm / own Hm, km |
|---|---|---|---|---|---|---|---|
| 50 | 20% | 24% | 25% | 6% | 25% | 59% | 62 / 47 |
| 100 | 15% | 22% | 17% | 13% | 33% | 49% | 72 / 59 |
| 150 | 20% | 24% | 17% | 12% | 27% | 56% | 68 / 52 |
| 200 | 24% | 30% | 13% | 7% | 26% | 56% | 68 / 55 |

About 40-55% of pulled length is floor under 30 Myr, and roughly half of it is floored at the cell's own Hm (47-59 km
at the pulled nodes), well above the age law (13 sqrt(age) - 11 km is 30 km at 10 Myr). Young floor pulls with an old
floor's column. This corrects section 2's "D cannot move T2": D would thin young floor and so lower its pull. How much
that moves T2 isn't measured.

### Angle test in place of the 1 cm/yr threshold (same world, same samples)

Transform where |normal| < sin(theta) |v_rel|. Shares of boundary-band length; "kept" is the share of today's pulled
length that the angle test would also pull (today's gate and down-going crust still apply).

| theta | Transform share of band | Closing share of band | Pulled share of band | Kept of today's pulled |
|---|---|---|---|---|
| today (1 cm/yr) | 9-17% | 41-44% | 41-42% | 100% |
| 15 deg | 13-17% | 39-42% | 36-38% | 84-90% |
| 30 deg | 31-33% | 31-33% | 30-31% | 70-74% |
| 45 deg | 47-50% | 23-25% | 22-24% | 52-57% |

Samples agree to within a few points at 50-200 Myr. These are band-wide totals on one world; how the removal falls
across plates (which is what T2 needs) is not measured here.

Per plate (pulled length over the plate's own band length; 20-27 plates per sample, same world):

| Myr | Test | p10 | Median | p90 | Plates under 5% |
|---|---|---|---|---|---|
| 50 | today | 0.32 | 0.45 | 0.55 | 0 of 21 |
| 50 | 30 deg | 0.22 | 0.32 | 0.45 | 0 |
| 50 | 45 deg | 0.15 | 0.25 | 0.38 | 0 |
| 100 | today | 0.00 | 0.46 | 0.52 | 2 of 20 |
| 100 | 30 deg | 0.00 | 0.32 | 0.44 | 2 |
| 100 | 45 deg | 0.00 | 0.23 | 0.34 | 2 |
| 150 | today | 0.17 | 0.40 | 0.48 | 2 of 25 |
| 150 | 30 deg | 0.13 | 0.28 | 0.38 | 2 |
| 150 | 45 deg | 0.10 | 0.20 | 0.31 | 2 |
| 200 | today | 0.27 | 0.44 | 0.53 | 1 of 27 |
| 200 | 30 deg | 0.18 | 0.34 | 0.44 | 1 |
| 200 | 45 deg | 0.13 | 0.24 | 0.32 | 1 |

An angle test scales every plate's pulled share down together (median 0.45 to 0.32 at 30 deg) and leaves the number of
near-zero plates where it is (1-2 of 20-27, as today). It does not create slab-free plates. The spread comes from
something else: the age and nature of the floor at each plate's trenches (section 4), or a slab history.

## 5. The age-dependence fix alone (cap on the force column's own-Hm floor)

`torque.FORCE_COLUMN_OWN_HM_CAP_M` (default infinity, the old form; unit-tested) caps the own Hm that floors B.3's
column. Setting `b123-age-pull` is `b123-brake-off` with the cap at 8 km (`YOUNG_RIDGE_HM_M`), so a dated oceanic cell's
column follows the cooling law. Same five stage 1 worlds, 200 Myr, 1 Myr steps; score `b123-age-pull-score.json`,
groups `t2-groups-age-pull.json`; raw runs in `~/mantle-bloom-captures/issue995/`.

| | b123-brake-off | with the cap | Target |
|---|---:|---:|---|
| T1 NNR cm/yr | 4.25 | 4.42 | 3-5 |
| T2 ratio (+10 Myr) | 1.59 | 1.66 | at least 2.5 (bar 2.0) |
| Slab-attached share of area | 56% | 62% | |
| T3 pull share | 60.0% | 59.7% | at least 50% |
| T4 MPa | 1.52 | 1.88 | at most 1 |
| T5 N/m | 1.91e12 | 1.97e12 | 1-4e12 |
| T7 | 1.54 | 1.57 | at most 3 |
| Group means attached / other, cm/yr | 5.16 / 3.24 | 5.28 / 3.17 | |

T2 moves from 1.59 to 1.66. The 0.5 Myr repeat of the baseline moved T2 by -4.9%, so this is inside the run-to-run
spread; it does not approach the 2.0 bar. T4 gets worse (it already failed). The fix is not enough on its own:
trenches on young floor were a large share of pulled length, but removing their excess column does not make plates
without slabs.

## 6. A slab-history gate on pull

`torque.SLAB_HISTORY_GATE_SPACINGS` (default infinity = off; unit-tested; the default-forms bit-identical test passes):
a node pulls only within 4 line spacings of oceanic floor its plate lost under a neighbour in the last 10 Myr
(`collision_polarity` consumption evidence, role lower). Setting `b123-slab-history`, same five worlds and settings as
section 5; score `b123-slab-history-score.json`, groups `t2-groups-slab-history.json`.

| | b123-brake-off | slab-history | Target |
|---|---:|---:|---|
| T1 NNR cm/yr | 4.25 | 4.43 | 3-5 |
| T2 ratio (+10 Myr) | 1.59 | 1.53 | at least 2.5 (bar 2.0) |
| Slab-attached share of area | 56% | 64% | |
| T3 | 60.0% | 59.4% | at least 50% |
| T4 MPa | 1.52 | 1.36 | at most 1 |
| T5 N/m | 1.91e12 | 1.94e12 | 1-4e12 |
| T7 | 1.54 | 1.49 | at most 3 |
| Median plate's pulled share of its band | 0.44 | 0.39 | |
| Plates with under 5% pulled | 0.1% | 0.7% | |

The gate removes about 11% of the pulled length and does not create slab-free plates. T2 does not move (and it
would be partly tautological if it did: the gate and T2's grouping read the same evidence). Reason: nearly every
closing trench has consumed some floor within 10 Myr, so consumption evidence is present on every plate that closes
anywhere. The "other" group still has a subducting node (100% of its area) and 16.6 against 18.2 (x 1e25 N m) of
attached plates' slab pull.

## 7. Earth's plate-boundary geometry (PB2002) against the model

`bin/debug/earth_boundary_shares.py` on Bird's PB2002 steps (`earth-boundary-shares.json`; 52 plates, 260,483 km of
boundary each counted once; the subducting plate is read from the pair's slash and checked on Kermadec, the Aleutians,
the Kuril and Cascadia). Per plate, pulled share = length of subduction zones where it is the subducting plate over
its perimeter, the same denominator idea as the probe's band. The model figures are section 4's and the per-plate table
(one world, defaults).

| | Earth (PB2002) | Model |
|---|---:|---:|
| Subduction share of boundary | 20% (counted once) | 41-44% of a plate's band (per plate) |
| Median plate's pulled share of its perimeter, area-weighted | 5.5% | 0.40-0.46 |
| p10 to p90, area-weighted | 0.9% to 28% | p10 0.00-0.32, p90 0.48-0.55 |
| Plates under 5% pulled | 35 of 52 (41% of area) | 0-2 of 20-27 |
| Boundary sliding more than closing or opening | 44% | 47-50% |
| Boundary with \|normal rate\| under 1 cm/yr (the transform class) | 45% | 9-17% |
| Boundary closing past 1 cm/yr | 27% | 43% (what pulls: 41-42%) |
| Median relative speed across boundaries, cm/yr | 3.0 | 7.4-8.8 |

Earth's pulled share is concentrated: Pacific 28%, Nazca-like and Philippine Sea plates 32-35%, Australia 19%; Africa,
Antarctica, South America, North America and Eurasia 1-6%; Somalia, Amur, Arabia and Okhotsk 0. About 5% of Earth's
area is on plates over 30% pulled.

- The sliding share matches (44% against 47-50%): the model's boundary orientation relative to motion is not the
  difference. The transform class does differ: Earth's 45% of boundary under 1 cm/yr is mostly because Earth's
  relative speeds are low (3 cm/yr median), so a 1 cm/yr threshold is a third of a typical rate. The model's median
  relative speed is 2.5-3 times Earth's, which empties the class. The earlier statement that Earth's trench share
  is about 35% of boundary was wrong; it is 20%.
- The model's plates all converge somewhere at about the same rate; Earth's contrast comes from a few plates with long
  trenches and most plates with almost none. That is the T2 contrast and the model has none of it.
- The T1 means agree (4.4 against 3.8 cm/yr), but relative speeds across boundaries do not (8 against 3), so the
  model's neighbouring plates move against each other more than Earth's.

## Not done

Review fixes to `g1_world_checks.py` (wet-cohort depth table, per-correction ledger audit, signed Hm residual,
`--hm-ledger` docstring); the baseline Hm gap; any run.

## 8. Relative speeds across contacts (contact graphs)

`bin/debug/plate_contact_graph.py` builds a plate graph (node = plate, edge = contact over 500 km, label =
length-weighted relative speed) for PB2002 and for generated worlds; `contact-graphs.png` / `.html` show Earth and three
seeds at 100 Myr with the b123 brake-off settings.

| Share of contact length | Earth | seeds 106533960 / 1097624328 / 591453549 |
|---|---:|---:|
| under 2 cm/yr | 39% | 0% / 0% / 7% |
| under 4 cm/yr | 61% | 5% / 0% / 22% |
| over 8 cm/yr | 9% | 49% / 42% / 45% |

Median relative speed by contact length: Earth 3.2 cm/yr, model 5.6-7.9. The model has no slow large-plate pairs.
Caveats: nodes are plates, not continents; one sample per world; one world (1097624328) had only 5 plates.

## 9. The contact brake explains the high relative speeds (sweep on seed 591453549, 100 Myr)

`contact_drag_multiplier` (overall) x `craton_drag_multiplier` (craton component), 13 runs
(`brake-sweep-summary.json`). b123 had the brake off (0.0); the default is 1.0.

| Overall | Craton | Median (cm/yr) | Under 2 | Under 4 | Over 8 |
|---|---:|---:|---:|---:|---:|
| Earth | | 3.2 | 39% | 61% | 9% |
| 0 | n/a | 7.2 | 7% | 22% | 45% |
| 0.25 | 0 / 0.5 / 1 | 3.3 / 3.9 / 3.2 | 29 / 24 / 24% | 64 / 53 / 61% | 7 / 3 / 1% |
| 0.5 | 0 / 0.5 / 1 | 1.2 / 0.9 / 1.8 | 69 / 70 / 54% | 86 / 83 / 68% | 0 / 0 / 5% |
| 0.75 | any | 0.5 | 77-83% | 88-92% | 0-1% |
| 1 | any | 0.4 | 78-81% | 89-91% | 0% |

The response is sharp between overall 0.25 and 0.5; the craton component matters little. At 0.25 the median and the
share under 4 cm/yr match Earth; the model still lacks Earth's few fast Pacific-rim contacts.

## 10. Mean plate speed falls over time (20 seeds, 200 Myr, plates-only, climate off)

Runs: brake off (`b123-brake-off`) and overall 0.25 / craton 1 (`brake-o0.25-c1`), 20 seeds each (3 of the first set
plus 17 drawn at random, none held out); `slowdown-*.json`, `slowdown-analysis.json`,
`bin/debug/slowdown_analysis.py`. Speed is the area-weighted NNR mean (T1, target 3-5 cm/yr).

- Last 50 Myr vs first 50 Myr: slower in 14 of 20 seeds (brake off, median -0.26 cm/yr, mean -0.56) and 19 of 20
  (brake 0.25, median -0.76, mean -0.64).
- Brake 0.25 puts mean speed at about 2-3 cm/yr, below T1's band for most seeds: the brake that fixes contact speeds
  costs some of T1.

Within-run correlations (each run's mean removed):

| | Brake off | Brake 0.25 |
|---|---:|---:|
| speed vs time | -0.42 | -0.47 |
| largest plate's area share vs time | +0.65 | +0.63 |
| speed vs largest share | -0.72 | -0.50 |
| speed vs time, given largest share | +0.09 | -0.23 |
| speed vs plate count, given time | +0.44 | -0.10 |

Plate count does not explain the slowdown. The largest plate's share does, fully with the brake off and partly with it
on. Larger plates are slower per plate: speed by area quartile (smallest to largest) is 9.1, 7.9, 6.4, 4.2 cm/yr brake
off and 3.0, 3.9, 3.4, 2.4 brake on. Samples within a run are not independent; read these as effect sizes.

## 11. Splits and merges

From the per-10-Myr log lines of the same 40 runs (merges = splits minus the change in plate count):

| 10-200 Myr | Brake off | Brake 0.25 |
|---|---:|---:|
| splits per run (mean) | 28.1 | 31.7 |
| merges per run (mean) | 13.7 | 6.9 |
| plates, 10 Myr -> 200 Myr (mean) | 15.2 -> 29.6 | 15.1 -> 39.9 |
| splits / merges, first 100 Myr | 8.4 / 3.5 | 9.7 / 1.9 |
| splits / merges, last 100 Myr | 18.2 / 9.4 | 20.0 / 5.0 |

- Splits are frequent (about 14-16 per 100 Myr) and rise with time; merges lag, so plate count roughly doubles. This
  does not support "too few rifts" in count.
- Splits are not slivers: in 15 earlier brake-off runs (`run.json`) the smaller daughter's median share is 0.27-0.33
  and none is under 10%.
- Parent and absorber attribution by 10-Myr plate-id matching (approximate, since several events per interval and
  ids can repeat): the then-largest plate is the parent of about 9-13% of new plates and absorbs about 20% of
  vanished ones. The largest plate's share grows in 59-61% of intervals.
- Not measured: whether the largest plates rift at the rate their area and stress imply (which needs the splitting
  plate's area at each event, which `step_all` does not record).

## 12. Which plates split (10 seeds x 2 settings, 200 Myr; `splits-*.json`, `splits-analysis.json`)

`run_rift_lines_d1.step_all` now records each split's parent area, size rank, plate count and the largest other plate.
273 splits with the brake off, 321 with overall 0.25 / craton 1.

| | Brake off | Brake 0.25 |
|---|---:|---:|
| parent is the largest plate | 7% | 6% |
| ... if any plate were equally likely | 5% | 4% |
| ... if splits were proportional to area | 28% | 23% |
| parent in top 3 | 25% | 21% |
| ... if equally likely | 15% | 12% |
| mean size percentile of parent (0 = largest) | 0.25 | 0.24 |
| mean parent share of planet | 8.1% | 6.9% |
| ... if proportional to area | 16.5% | 13.0% |

Splits favour larger plates over a uniform pick, but far less than area would: the largest plate is the parent of 6-7%
of splits against 23-28% under area-proportional splitting. The area-proportional rates are computed from the 10-Myr
samples of the same runs and are a reference, not an Earth claim. This is consistent with section 10: splits cut small
and mid-sized plates, and the largest plate's share grows. Not measured: why large plates resist rifting, and Earth's
large-plate rifting rate.

## 13. Why big plates rarely rift (rift-line state, 5 seeds x 2 settings, 200 Myr, code after #1002)

`bin/debug/big_plate_rift_probe.py` records, every 5 Myr, each plate's rift-line state; `bin/debug/rift_state_analysis.py`
summarises it (`rift-state-*.json`, `rift-state-analysis.json`). Big = at least 10% of the planet's area. Counts are
plate-samples, which repeat for a plate that stays in one state.

| Gate | Big (off / 0.25) | Smaller (off / 0.25) |
|---|---:|---:|
| line under half broken | 521 / 402 | 1995 / 2376 |
| no line | 49 / 65 | 1580 / 1826 |
| passes both tests | 8 / 5 | 240 / 997 |

- **Broken share of the line**: mean 0.115 / 0.101 on big plates, 0.162 / 0.315 on smaller ones (needs 0.5).
- **Pull is not the limit.** Pull on the line reaches its strength in 72-74% of big-plate samples. But the ratio P/S is
  small on big plates (median 1.4-2.1) and large on small plates (4.4 brake off, 14.9 at 0.25), because line strength
  grows with plate size. A small excess means slow strain, so the line stays under half broken. I have not checked the
  strain law itself.
- **Craton share is a weaker explanation than size.** Big plates are more cratonic (mean craton share 0.22 / 0.25 vs
  0.17 / 0.13, area vs craton Spearman 0.42-0.50). Broken share vs craton, controlling for log area: -0.13 (off) and
  -0.07 (0.25); broken share vs log area controlling for craton: -0.16 and -0.45. Pass rate by craton tercile is
  monotone with the brake off (10.8, 6.1, 2.6%) and not with it on (32, 13, 22%).
- **Rejected tests**: 76 of 236 (brake off) and 263 of 461 (0.25) break tests were rejected, almost all
  `stray_too_large` (a third damaged component larger than `a_stray` of the plate), so splits that clear both gates
  still fail often.

Caveats: this tests whether craton share predicts rift progress; it does not remove cratons, so it does not show cratons
are not the cause. An A/B with the craton resistance switched off would.

## 14. Literature: do failed rifts weaken or strengthen? (search snippets only, not full papers)

Correction: an earlier statement that Stein argued the Midcontinent Rift (MCR) is stronger was from memory and is not
supported. Stein et al.'s 2012 GSA abstract calls the MCR a weak intraplate zone (higher crustal temperature, faults),
so later reverse faulting on it is unsurprising. Stein, Stein and Elling 2018 (Tectonophysics 744) was not read beyond
its abstract.

| Process | Direction | Conditions and time scale | Source |
|---|---|---|---|
| Heat, faults, dikes | weakens | dikes through the lithosphere can cut a rift's resistance by about an order of magnitude | Nature Reviews Earth & Environment 2023 rift review |
| Thermal relaxation | strengthens, most in the thinnest part | during quiet periods the thinnest axis strengthens more than thicker surroundings, so later extension jumps sideways; failed rifts with little extension stay abandoned after short pauses, strongly extended ones need longer | Naliboff and Buiter, EGU2014-10482 |
| Critical relaxation time | strengthens | thinned zones strengthen until a neighbour is weakest; the time is shorter for hot lithosphere and low strain rates | Macquarie kinematic-model study |
| Cooling as a switch | strengthens | cooling can flip a rift between breakup and failure; a large drop in driving force halts it, a small one does not | Silva et al. 2025 Sci. Rep. |
| Mafic underplate | both | heats and weakens, but dense hard rock strengthens; depends on thermal state, timing, melt budget | Penn State record |
| Oslo Graben | weak | modelling: failure cannot be attributed to extra strength | Oxford record |
| Volatile depletion | strengthens (mantle) | failed Anza Rift mantle may be left refractory | EGU26-5886 |
| Dormant, not dead | weakness persists | failed rifts can reactivate after hundreds of millions of years if forces change; suture reuse favours sutures of large oceans | Nature Reviews 2023; EGU2015-1972 |

Model comparison: the model has one time scale, the 63 Myr healing of the warm-strain softening (`tau_heal_yr`) and no
strengthening on cooling. Numbers for the strengthening rate and abandonment time are still to be extracted (section 15).

## 15. Literature numbers for rift strengthening and abandonment

Read in full text: Silva and Choi 2025 (Sci. Rep. 15, 37023, "Continental rifts losing driving forces can still complete
breakup") and the Naliboff and Buiter EGU2014-10482 abstract (text extracted with `pdftotext`). Abstract only: Negredo,
Fernandez and Zeyen 1995 (EPSL 134, 87-98). Not readable: the Nature Reviews 2023 review (login redirect), so its
"order of magnitude" dike weakening and "hundreds of millions of years" reactivation remain search-snippet claims.

| Quantity | Value | Condition | Source |
|---|---|---|---|
| Critical relaxation time (minimum quiet time before deformation migrates) | a few Myr to over 130 Myr | stretching factor 1.65; shorter for hot lithosphere and low strain rate; no further breakdown in the abstract | Negredo et al. 1995 |
| Quiet time needed | small extension: stays abandoned even after short quiet phases; strong extension: longer | qualitative, no durations | Naliboff and Buiter 2014 |
| Their model | 800 x 240 km, 40 km crust, 120 km lithosphere, 1-2 cm/yr continuous rifting | cooling during pauses, strain-weakening healing included | Naliboff and Buiter 2014 |
| Rift-zone lithosphere thickness | 110 km at 8.5 Ma falls to 15 km by breakup at 10.8 Ma; in a force-reduced run it stayed near 90 km from 12.5 Ma | cooling thickening balances thinning by the reduced traction | Silva and Choi 2025 |
| Driving-force reduction | 25% can still break up under good timing; 50% generally fails | boundary traction 120-160 MPa; later start and slower reduction help | Silva and Choi 2025 |
| Rate at which a rift can still revive | 1 mm/yr for 8 Myr (6-14 Ma) then accelerated to breakup | the balance is current lithospheric strength against current driving force | Silva and Choi 2025 |
| Extension rates | 1.8 to 94 mm/yr in a constant-traction run; 0.2-0.7 mm/yr in failed runs | | Silva and Choi 2025 |

Silva and Choi's time scales are a few Myr to tens of Myr, with cooling thickening and driving force deciding the
outcome. Nobody in these sources gives a single strengthening time constant for a failed rift.

For scale (my arithmetic, not from the sources): a conductive e-folding time for a slab of thickness L is
L^2 / (pi^2 kappa); with kappa = 1e-6 m^2/s this is about 32 Myr for 100 km and 63 Myr for about 140 km. The model's
`tau_heal_yr` (62.8 Myr) is the diffusion time of lithosphere 140 km thick, so it is of the right order for cooling, but
the model applies it as decay of softening, not as strengthening of the line.

## 16. Share of big plates that break up (10 seeds per arm, 200 Myr, brake 0.25 x craton 1)

Statistic: of the plate-samples (every 5 Myr) where a plate holds at least 10% of the planet's area, the share where that
plate splits within the next 20 Myr. Bootstrap intervals resample seeds. Code: `bin/debug/big_plate_breakup_stat.py`;
numbers: `breakup-stat.json`. Arms: baseline; `heal250` (`tau_heal_yr` 250 Myr, so softening persists); `kcap07`
(`k_cap` 0.07).

| Arm | Big plate-samples | Break up within 20 Myr | 95% interval | Big plates that split | Splits per run (200 Myr) |
|---|---:|---:|---|---:|---:|
| Baseline | 966 | 20.3% | 17.5-23.5 | 30 of 50 | 33.8 |
| Persistent softening (heal250) | 940 | 14.6% | 12.0-19.2 | 26 of 48 | 31.2 |
| `k_cap` 0.07 | 741 | 30.5% | 25.6-37.6 | 35 of 49 | 67.2 |

- Persistent softening did not raise the big-plate break-up share (it is lower, with the intervals nearly touching). The
  hypothesis that weakening left by failed rifts helps big plates break is not supported at this setting.
- `k_cap` 0.07 raises it, but it roughly doubles splits of every size, so it is not a targeted fix.

Stratified, baseline:

| Stratum | Samples | Break up within 20 Myr | 95% interval | Implied hazard per Myr |
|---|---:|---:|---|---:|
| Craton share, lowest third | 322 | 30.7% | 20.9-43.1 | 1.8% |
| Craton share, middle third | 322 | 20.5% | 14.7-30.5 | 1.2% |
| Craton share, highest third | 322 | 9.6% | 4.1-18.4 | 0.5% |
| Big for under 25 Myr | 292 | 21.2% | 15.6-28.3 | 1.2% |
| Big for 25-50 Myr | 219 | 20.5% | 12.7-28.2 | 1.2% |
| Big for 50-100 Myr | 273 | 19.4% | 13.6-25.6 | 1.1% |
| Big for over 100 Myr | 182 | 19.8% | 8.6-30.7 | 1.1% |

Reading:
- Craton-rich big plates break up about three times less often than craton-poor ones. This agrees with the suspicion
  that too many cratons hold plates together, even though `craton_root_m` is 0 here, so the effect does not come from
  added line strength. Craton share is correlated with size and continental share, so it is not isolated.
- The hazard is flat in time spent big. On Earth, Pangea stayed together roughly 100-150 Myr before rifting, and
  Rodinia and Nuna longer (my recollection, to be checked against sources). The model breaks a plate about as often in its
  first 25 Myr as after 100 Myr, which is the clearest mismatch with the supercontinent cycle found so far.
- The brake-off arm (5 seeds) is the one that shows hazard rising with time big (16% under 25 Myr to 47% over 100 Myr),
  but its intervals are wide.

Limits: terciles and time bins pool plates from the same run, so samples are not independent; the bootstrap over seeds
accounts for that only partly. "Time big" starts at the first sample above 10%, not at assembly of a supercontinent. No
Earth reference rate has been measured here.

## 17. Craton share against size, and Earth supercontinent durations

### Craton share is not a size confound

`bin/debug/big_plate_craton_regression.py` fits a logistic regression of "splits within 20 Myr" on standardised area
share and craton share over the big plate-samples of each arm (`craton-regression.json`). Among big plates, area share and
craton share are barely correlated (0.11).

| Arm | Craton coefficient per SD | 95% interval (seed bootstrap) | Area coefficient per SD | 95% interval |
|---|---:|---|---:|---|
| Baseline | -0.57 | -1.11 to -0.15 | -0.13 | -0.33 to 0.25 |
| heal250 | -0.54 | -0.86 to -0.34 | -0.31 | -0.66 to 0.07 |
| `k_cap` 0.07 | -0.38 | -0.84 to -0.16 | -0.08 | -0.55 to 0.38 |

Craton-rich big plates split less in all three arms, the interval excludes zero each time, and area adds nothing once
craton share is in. Mechanism unknown: `craton_root_m` is 0, so cratons add no line strength in this code. Candidates to
test are the craton thinning resistance (`cratons.rift_resistance`), a craton's effect on the plate's pull and strength
along the line, and `stray_too_large` rejections. These are not checked.

### Earth durations (search snippets, mostly secondary sources; check primary papers before using)

| Supercontinent | Assembly to breakup | Source |
|---|---|---|
| Pangea | assembly about 320 Ma (fully assembled by 299-273 Ma), Central Atlantic opening 195-175 Ma, so roughly 75-125 Myr | [Torsvik, EGU abstract](https://meetingorganizer.copernicus.org/EGU2013/EGU2013-2515.pdf), [Heudes et al., EGU26-10868](https://meetingorganizer.copernicus.org/EGU26/EGU26-10868.html) |
| Rodinia | initial assembly about 1100 Ma, breakup 825-750 Ma, so about 300 Myr (Condie et al. give about 100 Myr for a narrower lifespan) | [Condie et al.](https://people.earth.yale.edu/sites/default/files/condie15.pdf) |
| Nuna | assembly 1700-1500 Ma, lifespan about 1600-1400 Ma, about 300 Myr | [Condie et al.](https://people.earth.yale.edu/sites/default/files/condie15.pdf), [EGU2018-5913](https://meetingorganizer.copernicus.org/EGU2018/EGU2018-5913.pdf) |
| Supercontinent cycle | 300-500 Myr (Wikipedia), about 600 +/- 100 Myr (EGU abstract) | [Supercontinent cycle](https://en.wikipedia.org/wiki/Supercontinent_cycle) |
| Farallon plate (oceanic) | split into Cocos and Nazca by about 23 Ma (rifting from about 25 Ma, complete by 19.5 Ma); slab-pull divergence or a fracture zone over the Galapagos hotspot proposed | [Lonsdale, Tectonophysics](https://recordcenter.sgc.gov.co/B6/21002100024423/Documento/Nativos/Apoyo/Anexo%201_ARTICULOS%20EN%20BASE%20DE%20DATOS%20REFERENCIA%20Y%20GEOQUIMICA/05%20Creation%20of%20the%20Cocos%20and%20Nazca%20plates%20by.pdf) |

My earlier recollection of Nuna breaking up at 1.8 to 1.3 Ga was not supported by these sources; the dates above replace it.

Reading against the model: a mean time as a supercontinent of about 100-300 Myr corresponds to an average hazard of
about 0.3-1% per Myr, against the model's flat 1.1% per Myr (mean about 90 Myr). The model is at the fast end of
Earth's range, not obviously below it. Earth's hazard is not flat either: Pangea survived its first tens of Myr, while a
newly big model plate splits at the same rate as an old one. The suspicion that big plates rift too rarely is not
supported by these numbers; the timing mismatch and the craton effect are the findings.

## 18. Time-dependent rifting: literature and what the model carries (search snippets, not full papers)

Earth: Pangea's break-up came roughly 75-125 Myr after assembly (section 17). Proposed reasons for the delay:

| Mechanism | What the sources say | Source |
|---|---|---|
| Thermal insulation under the supercontinent | Models give sub-continental heating of up to about 50 C (Yoshida 2013, cited in a search result); the size and existence are contested, and a 2012 abstract finds Earth-like models do not heat the sub-continental mantle above the sub-oceanic one on aggregation timescales. Post-breakup temperature anomalies of tens of degrees decay with e-folding about 50 Myr. | [GSA 2012 abstract](https://gsa.confex.com/gsa/2012AM/webprogram/Paper210141.html), [EGU2016-16611](https://meetingorganizer.copernicus.org/EGU2016/EGU2016-16611.pdf?pdf), [EGU2018-513](https://meetingorganizer.copernicus.org/EGU2018/EGU2018-513.pdf?pdf) |
| Upwelling or plume push beneath the supercontinent, with subduction retreat helping | A 3D spherical study concludes retreat of the surrounding subduction girdle helps, but sub-supercontinent upwelling is the essential driver; retreat stress is confined to about 600 km at the margin. | [Zhang et al.](https://par.nsf.gov/servlets/purl/10127684) |
| Slab pull from the girdle | Slabs pull plates quickly and tear continents apart; continental keels slow drift only when the continent is not attached to a subducting plate. | [PMC6821462](https://pmc.ncbi.nlm.nih.gov/articles/PMC6821462) |
| Thick-lithosphere control | Pangea's break-up phases followed a crescent of thick lithosphere around Tethys; the asthenosphere beneath Pangea was about 150 C warmer than below Panthalassa between 260 and 180 Ma. | [par.nsf.gov 10428142](https://par.nsf.gov/servlets/purl/10428142) |
| Immediate subduction trigger | Rifting began almost as soon as Pangea assembled, with subduction of the Pangean plate beneath its own margin proposed as the trigger. This conflicts with a long delay. | [Nature Geoscience](https://www.nature.com/articles/ngeo270) |

The sources disagree on which force dominates and on the delay itself. No source here gives a hazard or a relaxation time
that could be used as a model target.

What the model carries: baseline, big plates, 10 seeds.

| Time spent big | Samples | Median pull/strength | Median broken share | Break up within 20 Myr |
|---|---:|---:|---:|---:|
| under 25 Myr | 308 | 1.23 | 0.00 | 21.2% |
| 25-50 Myr | 231 | 1.22 | 0.02 | 20.5% |
| 50-100 Myr | 280 | 1.50 | 0.04 | 19.4% |
| over 100 Myr | 258 | 1.48 | 0.07 | 19.8% |

Break-up share by plate age (time since the plate formed): 22.8% under 10 Myr, 17.2% at 10-25, 21.6% at 25-50, 20.0% over
50 Myr. The hazard is flat in both clocks. Pull/strength and broken share do creep up with time big, so the line state has
some memory, but the splits do not follow it. The model has no mantle temperature or upwelling input to the rift line, so
neither insulation nor plume push can act on it as it stands; only slab pull enters the pull term.

Any time-dependent mechanism would first have to lower the early rate (about 21% in the first 25 Myr big), then raise
the later one. Adding weakening alone would raise the late rate and leave the early one where it is.

## 19. Are early big-plate splits a merger artifact? (baseline, 5 seeds, merges and daughter ids recorded)

The rerun (`merge-origin-chain.sh`, `rift-merge-<seed>.json`) reproduced the splits of the earlier runs exactly on all 5 seeds,
so recording merges changed nothing. Merges are not rare: 48 over the 5 seeds (6-16 per run) against 18-71 splits per run.

Big plate-samples (at least 10% of the planet's area) in the first 25 Myr after the plate became big, by how it became big
(a merge that kept its id, or a split daughter, in the 5 Myr before; otherwise it grew past 10% or was big at the start):

| How it became big | Samples | Break up within 20 Myr |
|---|---:|---:|
| Kept its id in a merge | 61 | 14.8% |
| Grew past 10% with no merge or split event | 40 | 35.0% |
| Big at the first sample (5 Myr) | 32 | 18.8% |
| Split daughter | 0 | none (daughters are never at least 10%) |

Reading: plates that become big through a merge split less often than average, so merging does not produce the early
splits. The highest rate belongs to plates that grew past 10% without an event, probably by gaining area at their margins.
The counts are small (30-60 samples per row, 5 seeds), so the 35% against 15% gap is suggestive, not established.
Failed rifts also call `reset_age`, so "plate age" in section 18 is not time since formation.

## 20. Strength ramp at the break-up test: no effect (3 arms x 5 seeds, 200 Myr, brake 0.25 x craton 1)

`big_plate_rift_probe.py --ramp EARLY LATE TAU_MYR` multiplies the line strength seen by the `P >= S` break-up test by
`late + (early - late) exp(-t/tau)`, where t is Myr spent at 10% of the planet or more. Plates under 10% and the damage
accumulation are untouched. Arms (`ramp-chain.sh`, `ramp-breakup-stat.json`):

| Arm (early, late, tau Myr) | Big plate-samples | Break up within 20 Myr | Splits per run |
|---|---:|---:|---:|
| Baseline, same 5 seeds | 424 | 20.0% | 39.6 |
| 3, 0.5, 50 | 424 | 20.0% | 39.0 |
| 3, 0.5, 100 | 429 | 18.4% | 38.4 |
| 1, 0.5, 50 | 424 | 20.3% | 39.6 |

On the three seeds I compared (591453549, 106533960, 390382636), the split times and plate ids are identical to the baseline,
with the same number of break tests, for both the 3x early boost and the weakening-only arm. So the ramp changed nothing there.

Why: the test is `P >= S` with `P > 0`, and an open line has `S = 0`. The line reaches the test only after at least half of
it is broken, when the remaining strength is small against the pull, so scaling it by 0.5 or 3 does not change the
outcome. The binding constraints are the broken-length share (`f_break`) and `split_masks`, not the strength comparison.
The earlier sample table (median pull/strength 1.2-1.5 for big plates) is read before the line has broken, not at the
test.

Consequence: a time-dependent strength cannot shape the break-up hazard through this test. It has to act on the damage
that accumulates in the line (strain rate, softening or `f_break`), or on which plates have a line to break.
This was a probe on baseline seeds only; the splitting seeds are not independent of the earlier runs, so the arms add
no new evidence about hazard beyond this null.
