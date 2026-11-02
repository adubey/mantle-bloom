# RFC 0010 probe: can a bent split cut avoid cratons?

[RFC 0003](../../rfc/abandoned/0003-rifting-around-cratons.md) tried to steer a plate split's cut off
cratons by tilting its great circle (PR #418). It rarely could: 12% of craton-crossing cuts
had an admissible craton-free alternative. A great circle moves as a whole, so this probe
asks whether a cut free to bend, a minimum cut through a strength graph
([RFC 0010](../../rfc/0010-rifts-follow-weak-lithosphere.md)), does better.

**Outcome: no, on the current craton layout.** At a 1,000 km corridor, 13 of 135
craton-crossing cuts (10%) get a craton-free admissible minimum cut, about the same as
#418's tilting. Where it is admissible it trims the craton a rift's band would thin by
about a quarter, at a 1% residual cost. The obstacle is the craton layout: each continent's interior is one craton, so there is no
route around it except a coastal sliver, which the size gate rejects. That is the
evidence for [RFC 0009](../../rfc/0009-craton-mosaic.md).

## Method

`bin/debug/prototype_mincut_split.py` replays fresh worlds with every default
(`compare_cratons.run_job`, cratons on, 5 Myr steps to 1 Gyr, 100 Myr checkpoints) on
seeds 7, 11 and 2024, one at a time, on a model identical to origin/main `1227c36`. A
read-only tap on `cratons.record_rift_cut` sees every cut that goes ahead (split or failed
rift) and computes, without changing the trajectory:

- the plate's 8-nearest-neighbour node graph, links over 2.5 spacings dropped;
- link cost = chord / spacing × (floor + mean strength), with strength = max(craton
  strength, keel term), RFC 0003's cut score;
- the minimum cut between nodes beyond a corridor of half-width W either side of the flow
  cut (`scipy.sparse.csgraph.maximum_flow`), for W = 500, 1,000, 1,500, 2,500 km and floor
  0.05 and 0.25;
- admissibility, RFC 0010's gates: both halves ≥ `SPLIT_MIN_NODES`, poles at least the
  relaxed separation apart, RMS residual ≤ 1.1 × the flow cut's, and each half's largest
  connected component ≥ 99% of it. Area misfit is reported, not gated: the corridor bounds
  the cut's movement instead;
- "through craton": any craton node within `FAILED_RIFT_BAND_MULT` (2.5) spacings of a node
  on a cut link, the band a failed rift thins. Measured the same way for the flow cut; it
  agrees with the great-circle band count to within two cuts (135 against 133);
- "line through craton": a craton node at either end of a cut link;
- **B**, operational band exposure: craton area in the band, summed over the 135
  craton-crossing cuts using the minimum cut where it is admissible and the flow cut where
  it isn't, over the same sum for the flow cuts;
- boundary length: the summed chord length of the links a cut crosses.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/prototype_mincut_split.py \
    --out $S/mincut_1gyr.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/prototype_mincut_split.py \
    --coverage-only --out $S/coverage.json
```

Runs took 794–1,075 s each, and reproduced the first version's trajectories cut for cut. Compact numbers: [mincut_probe.json](mincut_probe.json).

## Results (0–1 Gyr, three seeds pooled)

542 cuts went ahead (190, 180, 172); 135 crossed craton (69, 40, 26).

| Corridor W, floor 0.05 | Band craton-free | …and admissible | Line craton-free, admissible (L) | Failing size / poles / residual / connected | B (per seed; pooled) | Median residual ratio | Median misfit | Median length ratio |
|---|---|---|---|---|---|---|---|---|
| 500 km | 8 | 7 (5%) | 26 (19%) | 21 / 3 / 0 / 7 | 0.85, 0.84, 0.88; 0.85 | 1.001–1.004 | 4% | 0.87–0.89 |
| 1,000 km | 23 | 13 (10%) | 23 (17%) | 47 / 3 / 5 / 11 | 0.75, 0.66, 0.76; 0.73 | 1.007–1.016 | 9% | 0.80–0.83 |
| 1,500 km | 30 | 16 (12%) | 29 (21%) | 54 / 5 / 17 / 9 | 0.71, 0.64, 0.74; 0.69 | 1.015–1.028 | 13–16% | 0.75–0.78 |
| 2,500 km | 51 | 13 (10%) | 21 (16%) | 91 / 2 / 51 / 7 | 0.78, 0.78, 0.81; 0.78 | 1.065–1.093 | 29–30% | 0.62–0.64 |

"Failing" counts are over all 135 craton-crossing cuts, whether or not the minimum cut was
craton-free. Of the 135, 121 flow cuts also have craton on the line itself (67, 33, 21).
Floor 0.25 gives nearly the same cuts (8–12 band-free and 20–27 line-free admissible cuts
per corridor). The minimum cut takes at most 22 ms; the whole evaluation, about 33 ms per
cut.

- **The bent cut follows the flow.** Within 1,000 km it raises the fit residual by about 1%,
  and its boundary is 17–20% *shorter* than the great circle's, so the great circle wasn't
  the flow's best line either.
- **It rarely clears the craton.** Wider corridors find craton-free routes (51 of 135 at
  2,500 km), but they cut off a small piece of the plate: the size gate fails 91 times.
  Even the line alone clears craton in an admissible cut for at most 21% of cuts.
- **It trims the craton a rift would act on, modestly.** Counting only admissible minimum
  cuts, craton in the operational band falls by 27% at 1,000 km and 31% at 1,500 km. Wider
  corridors do worse, because more of their cuts fail the gates.
- **The control for RFC 0009's G1.** On today's layout L = 17% and B = 0.73 at 1,000 km,
  outside both of G1's building branches (L ≥ 20% and B ≤ 0.75, or L ≥ 50% and B ≤ 0.50).

The first version of this probe (reviewed in #421) counted `admissible` without the
connectivity gate, reported a link-count length ratio, and averaged band craton over all
minimum cuts, admissible or not. That made the cut look as if it halved band craton; the
operational figure above is the correct one. The craton-free admissible counts were
unchanged by the connectivity gate.

## Why: the craton is the interior

On craton-crossing cuts, the median plate is 34% craton (interquartile 18–51%), and craton
fills a median 44% of the flow cut's band (19–71%). At generation, craton covers 104, 82
and 46 M km² (20%, 16% and 9% of the surface) on seeds 7, 11 and 2024. The largest
craton-bearing plates carry 54 of 92 M km² (seed 7), 30 of 55 (seed 11) and 19 of 49
(seed 2024).

Seeding (875 km in) and formation (500 km in, 200 Myr quiet) make every continental
interior a single craton, and over 1 Gyr cratons rise to 88–97% of genuine continental area
on these runs (the compare_cratons checkpoints in this replay; #409 attributes the rise to
loss of non-cratonic crust). A cut that divides a continent into two viable halves has to
cross that interior. Steering the line, straight or bent, can't help until interiors have
weaker seams, which is RFC 0009's proposal.

## What this does not show

- Whether minimum cuts avoid cratons on a mosaic of cores and belts: RFC 0009's gate G1.
- How much rifting loss a quarter less craton in the band would save. Failed-rift
  loss scales with band craton directly, but most loss is stretch donation and divergent
  thinning along the boundary afterwards (#388), which this probe doesn't follow.
- "Craton-free" is strict: any craton node within ±312 km of the cut. A thinner band or a
  minimum craton thickness would count more cuts as clear.

## Rift settings and routes (2026-10-31)

Exploration for the detected-rift revision, not gate evidence. `bin/debug/probe_rift_settings.py`
looks for the four rift settings on saved worlds, then routes a rift from each site. Intact craton
interiors are blocked; edges within 200 km, cratons below half strength and heated interiors with an
eroded keel are open at a cost, and rock within 150 km of a craton costs 1.5 times as much. The great circle through the same site and direction is the
baseline. Thresholds are design choices, and no replay was run: each row is one snapshot, so these
are not rates. No save has plumes, so the keel rule itself was not exercised.

Saves, from the private B2 bucket and not available for download: seeds 317739385 (86 Myr),
447276134 (116), 308161711 (121), 628977350 (227) and 373936961 (244). Compact result:
[`rift-settings.json`](rift-settings.json).

| Setting | Sites | Separating route | Blocked cells entered | Routes touching craton | Routes dividing a craton | Line touches craton | Median tortuosity | Median smaller daughter |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Active | 0 (no plumes in any save) | — | — | — | — | — | — | — |
| Passive | 16 | 16 | 0 | 4 | 2 | 7 | 1.21 | 44% |
| Rollback | 0 | — | — | — | — | — | — | — |
| Transtension | 117 | 80 | 0 | 24 | 15 | 46 | 1.34 | 12.0% |

"Line touches craton" counts the separating routes whose same-site great circle has craton in its
band.

| Rule | Passive / transtension routes | Route cells inside craton | Beside craton | Routes dividing a craton | Largest craton share detached |
|---|---:|---:|---:|---:|---:|
| All craton plus a 1-cell buffer forbidden | 16 / 64 | 0 | not measured | 0 | 0 |
| Edges and weak cratons open | 16 / 84 | 421 | 690 | 19 | 35% |
| Margin cost added (current) | 16 / 80 | 421 | 185 | 17 | 25% |

A first margin trial, with weak cratons cheaper than the margin, raised route cells inside craton to
531 because routes stepped into them to avoid the margin. Weak cratons now cost at least the margin.

**Craton cost.** 17 of the 96 routes leave craton on both sides. The median detaches 9% of the craton
it divides, 8 detach more than 10%, and 2 more than 25%. These numbers were not tuned. Rifting craton
loss is mostly post-split divergent thinning and stretch donation, which this census cannot measure.

None of the 37 transtension failures was a true barrier: in all 37, open rock reaches the plate
boundary on both sides of the site, though that alone does not prove a larger route exists. In 34, a
route did separate the plate but its smaller daughter stayed under 5%, typically about 1.5%, because
the router takes the cheapest nearby boundary from a site only 300 km inboard. In 3, no separating
route was found at any reach. The probe's `status` records only the last reach tried, so
`no_first_arm` in the JSON does not mean no arm exists. 19 of the 117 sites began on a blocked cell,
including 11 of the failures. The probe's `meets_min_nodes` field fails on most rows, including
near-even splits, because `SPLIT_MIN_NODES` scales with density and not plate size.

![A routed rift on seed 373936961](rift-route-around-craton.png)

*Transtension route on plate 20, seed 373936961 at 244 Myr, tortuosity 1.79, drawn before the margin
cost was added. Gold cells are craton;
the red route divides the plate into blue and green daughters. The route follows the craton's edge
as designed, but at right the arms hook back along the plate margin and carve off a thin strip.
That is a route-quality defect the separation test does not catch, and tortuosity alone may not
either. The image was chosen by the probe as this world's largest great-circle craton overlap, so
it is not a typical route.*

**Why two settings found nothing.** All 16 local saves have zero plumes, so active rifting was not
tested. Rollback is a real gap: every oceanic node counts as down-going, so an ocean-ocean front has
no upper plate. Only 134 and 184 of 1,922 and 9,445 ocean-facing convergent nodes qualified on the
two worlds checked, and on the 244 Myr world 37.5% of those moved away from the trench faster
than 1 cm/yr, in runs too short for the 300 km front.

**What this does not show.** Transtension is over-counted, with 70 of 117 sites from one world,
because a front is broken into fragments. The sites are snapshots, not sustained settings. The
route was not run in a history, so this does not measure plate count, craton loss or daughter
motion after the split. Rollback's upper-plate signal is a proxy: upper-plate motion does not explain
every subduction zone [Heuret & Lallemand, 2005](https://doi.org/10.1016/j.pepi.2004.08.022).

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_rift_settings.py \
  --saves <saves> --out analysis/rfc0010/rift-settings.json --image-dir <private-dir>
```

## Gate D1: census (2026-11-01)

**Outcome: Gate D1 does not pass.** Active, passive and transtension produce events on all three
seeds. Rollback misses under the frozen rule, because ocean-ocean polarity failed its test, so
oceanic upper plates could not count and continental ones produced one event on one seed. The rate
is not compatible with fewer than 40 plates on any seed. Both are findings to diagnose, as the gate
expects. No detected rift acted on these worlds.

Settings were [frozen on #737](https://github.com/adubey/mantle-bloom/issues/737#issuecomment-6560714067)
before the runs: seeds 1183015207, 1486436815 and 1173383782, default worlds with plumes on, 1 Myr
steps to 200 Myr, base `8ebea73b`, probe sha256 `7d287ffb…`. Rates below are per world per 100 Myr.
Tables: [`gate-d1.json`](gate-d1.json).

| Setting | Events | By seed | Rate | Rate after gating | Median life, Myr | Median route growth, Myr | Closed routed events that outlast growth |
|---|---:|---:|---:|---:|---:|---:|---:|
| Active | 74 | 34 / 14 / 26 | 12.33 | 7.33 | 21.0 | 20.9 | 20 of 45 |
| Passive | 54 | 15 / 13 / 26 | 9.00 | 6.67 | 7.0 | 22.9 | 8 of 22 |
| Rollback, oceanic upper | 13 | 0 / 12 / 1 | 2.17 | 1.83 | 6.0 | 23.5 | 0 of 7 |
| Rollback, continental upper | 1 | 1 / 0 / 0 | 0.17 | 0.00 | — | — | — |
| Transtension | 145 | 50 / 46 / 49 | 24.17 | 8.67 | 11.5 | 22.8 | 29 of 88 |

Median life and "closed" exclude events still open at 200 Myr, and include events a production
split cut short. So `outlast` understates what a rift would get on passive and transtension, where
production splits ended 18 and 20 events. Median life is over closed routed events.

**Ocean-ocean polarity fails.** `collision_polarity` cannot supply a lower plate on an ocean-ocean
front. 21.3% of 7,562 qualifying front-steps were decided against the 80% limit, and 10.4% of 1,123
fronts flipped against the 10% limit; 78% were ambiguous. The read itself works: continent-ocean
fronts were 96.3% decided and 99.5% ocean-lower. The cause is in the model. Every contested oceanic
node on both plates may retreat, so nothing selects a lower plate, and at ocean-ocean pairs the plate
earlier in that step's shuffled deform order lost 80.8% of the cells (31,659 of 39,174). Only the
earlier plate lost cells in 3,337 pair-steps, only the later in 718. So the evidence mostly records
the shuffle. The per-cell `node_created_years` field is a candidate replacement, which these runs
did not test: a cell gets the creation time once, and remaps keep the earliest. On the 200 Myr saves,
the 51 ocean-ocean fronts have median side-to-side age gaps of 21–32 Myr by seed, and 43 of the 51
exceed 5 Myr. But 24–60% of ocean cells date from t = 0 and tie.
The 13 oceanic-upper events rest on this read and are not
reliable; 12 are on one seed.

**Plate count is incompatible on the frozen worst case.** Every seed fails, and fails even when only
separating routes count.

| Seed | Final plates | Gated events per 100 Myr | Production splits per 100 Myr | Worst case | Separating routes only |
|---|---:|---:|---:|---:|---:|
| 1183015207 | 26 | 27.0 | 5.0 | 70 | 48 |
| 1486436815 | 30 | 22.5 | 7.5 | 60 | 50 |
| 1173383782 | 26 | 24.0 | 9.0 | 56 | 45 |

This bound assumes every gated event breaks up, which D2 will lower. But production recorded no
merges in 600 Myr of history, and plate count climbs on every seed even at 5–9 splits per 100 Myr.
Nothing removes plates at the rate rifts could add them.

**Passive after a head.** 3 of 54 passive events follow a continental head on the same plate.
Production splits took another 25 passive tracks before onset; those are not events.

**After the freeze, reading the same reports (not frozen tests).**

- **Passive is permanently detected, on small plates.** The RFC feared this on very large plates.
  Instead, of 29 passive events open at 200 Myr, 14 have run over 100 Myr, the longest 196 Myr, on
  plates of 2–5% of the sphere. `relax` never exceeded 0.78 on any event, so relaxed gates
  are not the cause. Probably `plate.split` refuses these plates on size, so production never clears
  them, but this has not been traced.
- **Active outnumbers heads.** 74 events against 51 heads and 15 credited continental heads, and only
  16 events follow a credit on the same plate. The candidates are plume tails heating continents and
  heads that reach a plate other than the one credited; neither has been traced.
- **Transtension still repeats.** One plate gave 18 events on one seed, and gating removes about two
  thirds of the events.
- **Rollback fronts are brief.** Median 6 Myr, against about 24 Myr of route growth. None of the 7
  closed routed events outlast it.
- **Seed 1173383782 starts at 8 plates** at 1 Myr and has 18 by 11 Myr, on production alone. That is
  under D3's lower limit of more than 10 at every output.

Fingerprints are `cfe10acf0a71`, `f27b4ada46d4` and `30ecd0440c6f`. The fingerprint matched a
`--no-detect` control on a development seed, not on these seeds, to 12 Myr. Detection took 5–7% of
runtime; peak RSS was 0.6 GB. Step-size sensitivity was not tested. Reports and 200 Myr saves are
outside the repo, with checksums.

### Rollback amendment: cell age (2026-11-01)

**Outcome: the age read fails, and Gate D1 still does not pass.**

A [labeled rerun](https://github.com/adubey/mantle-bloom/issues/737#issuecomment-6561849550) of the same three seeds, with ocean-ocean polarity read from cell age. Probe hash unchanged during the runs. The model is untouched.

**Integrity:** all three fingerprints match D1, and active, passive and transtension events are identical on every seed. Only rollback changed.

| Age read | Result | Limit |
|---|---:|---:|
| Front-steps decided | 65.8% of 7,562 | at least 80% |
| Fronts reversed | 15.3% of 502 (77) | at most 10% |
| Decided after 100 Myr (reported only) | 82.2% | none |
| Agreement with consumption evidence | 63.1% | none |

Fronts here are tracked across gaps of up to 5 Myr, so the 502 is not comparable with D1's 1,123. 69 of the 77 reversals are on seed 1486436815.

| Rollback | Events by seed | Rate | After gating | Median life, Myr | Median route growth, Myr | Closed routed events that outlast growth |
|---|---:|---:|---:|---:|---:|---:|
| Oceanic upper | 0 / 10 / 0 | 1.67 | 1.33 | 12.0 | 19.6 | 2 of 5 |
| Continental upper | 1 / 0 / 0 | 0.17 | 0.00 | — | — | — |

**Reading:** rollback misses however the rule is taken. Under the frozen rule, oceanic events can't count. If they did, they are still on one seed of three. Age is a better read than consumption: decided share rises from 21.3% to 65.8%, and its events last twice as long (12 against 6 Myr). Its early failure is expected, since starting ocean is one age, but the reversals are not explained by that.

**Likely cause of the miss, read after the freeze:** ocean-ocean fronts are scarce on two worlds. Seed 1486436815 holds 6,409 of the 7,562 front-steps, against 595 and 558 on the others. Rollback may be limited by how few oceanic plates those worlds have, more than by polarity. This has not been traced further.

**Plate count** barely moves: worst case 70, 58 and 55, or 48, 48 and 44 on separating routes only. Every seed remains incompatible on the frozen bound.

Tables: [`gate-d1-rollback-amendment.json`](gate-d1-rollback-amendment.json). Sensitivity to the 10 Myr gap is untested.

**Why the age read missed, measured afterward** with `bin/debug/probe_age_polarity.py` on the same
seeds, fingerprints matching. Undecided: in 2,261 of 2,584 undecided front-steps (87.5%), at least half
the cells on both sides are starting ocean, so the median ages tie exactly; the other 323 are real
gaps under 10 Myr. A smaller gap would not help the 87.5%. These persist after 100 Myr (740 steps).
Reversed: the sampled cells change, because a fixed pair of cells keeps its gap. Between freeze and
reversal a median 25% of the front's nodes are still in place, 61% of reversals keep under half, and
only 4% keep them all; the sampled median age shifts by a median 16–33 Myr per side, and half of
frozen gaps were under 20 Myr. So a tracked front slides along the boundary onto other crust, and a
single verdict per front is too coarse or the track too loose. The 15.3% divides by all 502 fronts.
Among the 327 fronts that ever froze a verdict, 77 (23.5%) reversed. By seed that is 69 of 276, 8 of
39 and 0 of 12, so seed 1486436815 has more fronts, not a clearly worse rate.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_rift_settings.py \
  --seed <n> --myr 200 --step-myr 1 --out <seed.json> --save-dir <private-dir>
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_rift_settings.py \
  --summarize <a.json> <b.json> <c.json> --out analysis/rfc0010/gate-d1.json
```

### Rerun on dated rock (2026-11-01)

A [labeled rerun](https://github.com/adubey/mantle-bloom/issues/737#issuecomment-6569971922) of the
same three seeds on base `b3356fa6`, commit `8a082ebe`, after #745 and #752 dated starting rock.
Probe hash unchanged during the runs. It's the same frozen census, except that `age_verdict` reads
negative dates as given. Tables: [`gate-d1-rerun.json`](gate-d1-rerun.json).

**Outcome: the four-settings rule is met, narrowly, but the age read fails and the rate is still not
compatible with fewer than 40 plates.** Under the rerun's frozen rule, rollback pauses before D2.

| Age read on ocean-ocean fronts | Amendment | Rerun | Limit |
|---|---:|---:|---:|
| Front-steps decided | 65.8% | **95.9%** of 7,815 | at least 80% |
| Fronts reversed, of all fronts | 15.3% of 502 | **20.7%** of 487 (101) | at most 10% |
| Reversed, of fronts that froze a verdict | 23.5% of 327 | 24.3% of 416 | none |
| Agreement with consumption evidence | 63.1% | 46.7% | none |

Dating fixed the undecided share. It did not touch reversal: among fronts that froze a verdict, the
rate is unchanged. The headline rose because more fronts now freeze one.

| Events by seed | D1 and amendment | Rerun | Rerun rate, per world per 100 Myr |
|---|---:|---:|---:|
| Active | 34 / 14 / 26 | 38 / 13 / 27 | 13.00 |
| Passive | 15 / 13 / 26 | 11 / 14 / 25 | 8.33 |
| Rollback, oceanic upper | 0 / 10 / 0 | 0 / 17 / 1 | 3.00 |
| Rollback, continental upper | 1 / 0 / 0 | 1 / 1 / 0 | 0.33 |
| Transtension | 50 / 46 / 49 | 38 / 43 / 54 | 22.50 |

`main` also moved in hydroclimate between these columns, so a change can't be credited to ages alone.

**Rollback's pass is thin.** The failed read bars oceanic events from carrying the setting, so it passes
on two continental events. One starts at 3 Myr, lasts 2 Myr and ends in a production split; the other
starts at 155 Myr and lasts 25. Oceanic events would be on two seeds too, but 17 of 18 are on one seed,
and 1 of 10 closed routed events outlasts its route growth.

**Plate count.** Worst case 64, 61 and 58, or 46, 49 and 48 on separating routes only; every seed
remains incompatible. Production ends at 28, 27 and 23 plates with one merge in 600 Myr.

**Passive after a head:** 1 of 50. Fourteen passive events again last over 100 Myr, up to 191.

**Why fronts reversed, measured afterward** with `probe_age_polarity.py` on the same seeds, fingerprints
matching. Between freeze and reversal, a median 20% of the front's nodes are still in place, 67% of
reversals keep under half, and 8% keep them all. The median lag is 10 Myr. So the amendment's reading
holds: tracks slide onto other crust, and dated crust doesn't help. Only 32% of reversals started
from a gap under 20 Myr, so a wider gap would not clear them, and 53% end with a gap of at least 20 Myr.
By seed, 4 of 31, 89 of 353 and 8 of 32 frozen fronts reversed. Of 319 undecided steps, 115 have
both medians on starting rock and 204 have younger rock, all inside the 10 Myr gap. Whether a frozen
verdict or a per-step one is right can't be read from these runs; there is no independent truth for
ocean-ocean polarity, as consumption evidence mostly records the deform shuffle.

The dated pattern is a proxy, and the gap's sensitivity, step size and seed count remain untested.

**Owner decision (2026-11-02), after the rerun.** The project owner accepts the age read despite its
20.7% reversal, judging the 10% limit too strict against a geological estimate of about 15%. This
overrides a frozen limit after the result was seen; it is not a pass. Two cautions stand. Most
reversals come from tracks sliding onto other crust, not from plates changing roles, so they don't
measure the same thing as the 15% estimate. And 20.7% still exceeds that figure. With oceanic-upper
events counted, rollback has events on all three seeds (1 / 18 / 1), so all four settings pass the D1
rule. The plate-count bound still fails, and the owner chose to start D2 anyway, since D2's route
outcomes are what can tighten that worst case.
