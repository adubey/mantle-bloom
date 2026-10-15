# RFC 0010 probe: can a bent split cut avoid cratons?

[RFC 0003](../../rfc/0003-rifting-around-cratons.md) tried to steer a plate split's cut off
cratons by tilting its great circle (PR #418). It rarely could: 12% of craton-crossing cuts
had an admissible craton-free alternative. A great circle moves as a whole, so this probe
asks whether a cut free to bend, a minimum cut through a strength graph
([RFC 0010](../../rfc/0010-rifts-follow-weak-lithosphere.md)), does better.

**Outcome: no, on the current craton layout.** At a 1,000 km corridor, 13 of 135
craton-crossing cuts (10%) get a craton-free admissible minimum cut, about the same as
#418's tilting. The cut still halves the craton in its band, at a 1% residual cost. The
obstacle is the craton layout: each continent's interior is one craton, so there is no
route around it except a coastal sliver, which the size gate rejects. That is the
evidence for [RFC 0009](../../rfc/0009-craton-mosaic.md).

## Method

`bin/debug/prototype_mincut_split.py` replays fresh worlds with every default
(`compare_cratons.run_job`, cratons on, 5 Myr steps to 1 Gyr, 100 Myr checkpoints) on
seeds 7, 11 and 2024, one at a time, at origin/main `1227c36` with no model changes. A
read-only tap on `cratons.record_rift_cut` sees every cut that goes ahead (split or failed
rift) and computes, without changing the trajectory:

- the plate's 8-nearest-neighbour node graph, links over 2.5 spacings dropped;
- link cost = chord / spacing × (floor + mean strength), with strength = max(craton
  strength, keel term), RFC 0003's cut score;
- the minimum cut between nodes beyond a corridor of half-width W either side of the flow
  cut (`scipy.sparse.csgraph.maximum_flow`), for W = 500, 1,000, 1,500, 2,500 km and floor
  0.05 and 0.25;
- #418's gates that still apply: both halves ≥ `SPLIT_MIN_NODES`, poles at least the
  relaxed separation apart, RMS residual ≤ 1.1 × the flow cut's. Area misfit is reported,
  not gated: the corridor bounds the cut's movement instead;
- "through craton": any craton node within `FAILED_RIFT_BAND_MULT` (2.5) spacings of a node
  on a cut link. Measured the same way for the flow cut; it agrees with the great-circle
  band count to within two cuts (135 against 133).

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/prototype_mincut_split.py \
    --out $S/mincut_1gyr.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/prototype_mincut_split.py \
    --coverage-only --out $S/coverage.json
```

Runs took 730–813 s each. Compact numbers: [mincut_probe.json](mincut_probe.json).

## Results (0–1 Gyr, three seeds pooled)

542 cuts went ahead (190, 180, 172); 135 crossed craton (69, 40, 26).

| Corridor W, floor 0.05 | Craton-free cut | …and admissible | Failing size / poles / residual | Craton in band, vs flow cut (per seed) | Median residual ratio | Median misfit | Median length ratio |
|---|---|---|---|---|---|---|---|
| 500 km | 8 | 7 (5%) | 21 / 3 / 0 | 0.83, 0.77, 0.82 | 1.001–1.004 | 4% | 0.90 |
| 1,000 km | 23 | 13 (10%) | 47 / 3 / 5 | 0.66, 0.48, 0.60 | 1.007–1.016 | 9% | 0.84 |
| 1,500 km | 30 | 16 (12%) | 54 / 5 / 17 | 0.56, 0.32, 0.44 | 1.015–1.028 | 13–16% | 0.79 |
| 2,500 km | 51 | 13 (10%) | 91 / 2 / 51 | 0.34, 0.14, 0.25 | 1.065–1.093 | 29–30% | 0.65 |

"Failing" counts are over all 135 craton-crossing cuts, whether or not the minimum cut was
craton-free. Floor 0.25 gives nearly the same cuts (9–12 admissible craton-free cuts per
corridor). Both halves of every admissible craton-free cut are connected. The minimum cut
takes at most 22 ms; the whole evaluation, about 33 ms per cut.

- **The bent cut follows the flow.** Within 1,000 km it raises the fit residual by about 1%
  and is 16% *shorter* than the great circle, so the great circle wasn't the flow's best
  line either.
- **It rarely clears the craton.** Wider corridors find craton-free routes (51 of 135 at
  2,500 km), but they cut off a small piece of the plate: the size gate fails 91 times.
- **It does reduce the craton it crosses.** Craton area in the band falls by a third to a
  half at 1,000 km, and by two thirds or more at 2,500 km.

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
- How much rifting loss a cut with half the craton in its band would save. Failed-rift
  loss scales with band craton directly, but most loss is stretch donation and divergent
  thinning along the boundary afterwards (#388), which this probe doesn't follow.
- "Craton-free" is strict: any craton node within ±312 km of the cut. A thinner band or a
  minimum craton thickness would count more cuts as clear.
