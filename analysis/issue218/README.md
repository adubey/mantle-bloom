# Issue 218: initial terrain from a synthetic tectonic history

Freshly generated worlds don't look like Earth; worlds that have run for tens of Myr do. This
records the exploration in [#218](https://github.com/adubey/mantle-bloom/issues/218) and the
measurements behind `synthetic_relief.py`. The issue's Exploration summary has the
consolidated findings and the questions an RFC should settle. `results.json` holds every
number below.

## Setup

- origin/main `7f20853`, cratons-first layout, `node_density=4`, `land_fraction=0.29`.
- Seeds 447276134, 373936961, 910211954, 374464662, whose late saves (132–244 My) served as
  visual references.
- Earth: GMT `earth_relief_10m` (SRTM15+), area weighted. Lowland relief is box-averaged to
  ~55 km to match the ~60 km node spacing.
- Measures are per node. Cells vary up to ~1.4x in area, so per-node shares can differ
  slightly from area shares; generation keeps the land *node* count.
  - *Lowland relief* is the RMS elevation difference between land points 100 km apart, both
    below 1 km.
  - *Backed* means elevation equals `isostatic_elevation(Hc, Hm)`.

Regenerate the implementation rows with:

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_initial_relief.py \
    --seeds 447276134,373936961,910211954,374464662 --cratons-first --synthetic-relief \
    --out relief.json
```

Probes B1–B3 were scratch, elevation-only overwrites (crust thickness not rebalanced),
preregistered on #218 before running. Their designs are in the issue comments.
`synthetic_relief.py` implements B3 with every column backed by Hc.

## Results at t = 0

Ranges over the four seeds.

| world | land p50 (m) | land > 1 km | land > 2 km | ocean IQR (m) | ocean > −3 km | lowland relief at 100 km (m) | lakes, share of land |
|---|---|---|---|---|---|---|---|
| plates-first (1 seed) | 1665 | 0.65 | 0.38 | 303 | — | — | — |
| cratons-first | 556–1211 | 0.29–0.56 | 0.11–0.29 | 474–1688 | 0.19–0.25 | 62–73 | 0.006–0.033 |
| B1: layout and first motion | 351–471 | 0.09–0.19 | 0.008–0.05 | 2784–2893 | 0.30–0.32 | — | — |
| B2: + synthetic history, redrawn coasts | 405–591 | 0.17–0.29 | 0.07–0.14 | 2611–2924 | 0.35–0.37 | 81–104 | 0.016–0.032 |
| B3: + curved sutures, interior types | 308–499 | 0.13–0.26 | 0.05–0.14 | 2611–2924 | 0.35–0.37 | 143–151 | 0.059–0.104 |
| B3 at `continental_fraction` 0.4 | 338–531 | 0.08–0.28 | 0.014–0.12 | 1530–1674 | 0.21–0.23 | 142–154 | 0.055–0.094 |
| `synthetic_relief.py` (backed) | 330–473 | 0.15–0.24 | 0.053–0.092 | 2589–2647 | 0.32–0.34 | 150–157 | 0.070–0.101 |
| Earth | 423 | 0.26 | 0.107 | 1818 | 0.24 | 166–173 | ≈ 0.02 |

The implementation is fully backed: the largest |elevation − isostatic_elevation(Hc, Hm)| on
an ice-free node is 0 m, and no column is outside the Hc caps. It differs from B3 in its noise
streams, in leaving oceanic crust as generated (B3 drew unbacked ridges and trenches), and in
one fix. B3 measured the suture coordinate from each node's nearer craton, which is never
negative, so its warp lowered or split belts instead of moving them; the implementation
orients it by craton pair.

## Findings

- **Synthetic history supplies the missing relief.** Relief placed only from the first
  predicted motion (B1) leaves continents almost flat, because predicted convergence rarely
  touches continental crust. Sutures between neighbouring cratons and leading-edge
  cordilleras put land hypsometry on Earth's without any remap (B2).
- **Coastlines from distance to the crust edge** remove the noise-quantile inland seas and
  blob islands.
- **Interior terrain types** (shield, platform, hill country, plains) roughly double lowland
  relief, to within ~10% of Earth's.
- **The shallow-ocean surplus is continental crust.** At the 0.5 default, ~30% of the ocean
  floor is continental crust; at 0.4 the share of ocean shallower than 3 km matches Earth's.
- **Hc vs Hm.** Hm moves elevation about 11× less per km than Hc, so its main effect at
  generation is the thermal regime: collapse mobility, melting and delamination eligibility.
  Young belts start with burial lag, so they heat over ~25 Myr instead of starting hot.
- **The ocean can't be backed.** An 11 km oceanic column floats no shallower than about
  −4.4 km, so age-depth ridges and young floor need oceanic thermal isostasy, which the
  model lacks.
- **Lakes overshoot.** Scour pits and platform basins are sub-grid features that become
  node-sized closed basins: 6–10% of land against Earth's ~2%.

## Early evolution

Seed 910211954 stepped to 10 Myr in 1 Myr steps (`--step-myr 1,5,10`), with and without
`synthetic_relief`:

| world | t (Myr) | sea level (m) | land fraction | land p50 (m) | land > 2 km | ocean 0–200 m | lowland relief at 100 km (m) |
|---|---|---|---|---|---|---|---|
| synthetic relief | 0 | 0 | 0.283 | 379 | 0.092 | 0.070 | 156 |
| | 1 | −34 | 0.291 | 408 | 0.093 | 0.065 | 151 |
| | 5 | −166 | 0.331 | 494 | 0.098 | 0.021 | 143 |
| | 10 | −270 | 0.340 | 599 | 0.111 | 0.018 | 130 |
| cratons-first control | 0 | 0 | 0.283 | 556 | 0.127 | 0.027 | 68 |
| | 1 | −36 | 0.287 | 585 | 0.128 | 0.024 | 71 |
| | 5 | −167 | 0.298 | 704 | 0.148 | 0.014 | 92 |
| | 10 | −267 | 0.302 | 807 | 0.164 | 0.012 | 107 |

- The world steps without a jump: belts persist, rivers organize, and collapse-mobile land
  stays at 5–6%.
- Both worlds lose ~270 m of sea level in 10 Myr (the model's existing eustatic drift). With
  Earth-like shelves at 20–150 m, that drains the shelves: land grows by 5.7 points against
  1.9 for the control, and the 0–200 m share falls from 0.070 to 0.018. Realistic shelves make
  the existing sea-level drift visible quickly; an RFC should pair this relief with a
  sea-level budget that holds.
- Lowland relief erodes slowly (156 → 130 m), while the control roughens (68 → 107 m) as
  tectonics adds relief. The two converge, so the initial roughness is not lost in a step.
- Unbacked elevation grows equally in both runs (up to ~800 m on ~10% of nodes by 10 Myr),
  from stepping processes that move elevation without Hc (`back_elevation_gain_fields`
  documents this). It is not introduced by generation.

## Initial drainage

Generated worlds drained far less land to the sea at t = 0 than any stepped world.
`drainage.json` holds the numbers below. Measures are per node, as above.

- *Closed drainage* is land whose routed water ends somewhere other than the ocean. Rivers
  stop at a node below freezing, so ends at a frozen node are counted apart.
- *Depression depth* is how far water at a land node must rise before it can reach the
  ocean, from a priority flood over the hydrology graph. The Earth reference floods
  `earth_relief_10m` box-averaged to ~55 km (`bin/debug/measure_earth_depressions.py`).

**Cause 1: the river bootstrap had no time.** `initialize_river_channels` ran one hydrology
pass with `years=0`. That pass gave depression breaching no carving budget and let no basin
fill within the step (`_spill_within_step`), so every pit ended a river. The pass now gives
those two a drainage history of 1 Myr (`INITIAL_DRAINAGE_HISTORY_YEARS`); ice and the other
rates still see a zero-length step. Before is the same code with the history set to 0, which
matches origin/main `34fb86f`.

| world (4 seeds) | unfrozen closed drainage | frozen-end drainage | lakes, share of land |
|---|---|---|---|
| cratons-first, before | 0.31–0.48 | 0.12–0.48 | 0.006–0.033 |
| cratons-first, after | 0.24–0.42 | 0.14–0.49 | 0.004–0.030 |
| synthetic relief, before | 0.55–0.70 | 0.17–0.29 | 0.061–0.106 |
| synthetic relief, after | 0.48–0.64 | 0.18–0.29 | 0.060–0.086 |
| saves (6–220 My) after one 1 Myr step | 0.05–0.26 | 0.14–0.53 | 0.004–0.066 |
| Earth | ~0.19 (endorheic) | — | ~0.02 |

A second bootstrap pass with the climate held fixed changes little (0.483 → 0.466 on
910211954 with synthetic relief), so routing settles in one pass. A second pass with the
climate recomputed lowers closed drainage further (to 0.32–0.50 with synthetic relief),
because the new lakes wet the climate.

**What remains is climate and terrain.** At t = 0, median land precipitation is 55–226 mm/yr,
and 52–72% of land gets under 250 mm (saves: 155–354 mm and 44–59%). About 70% of the
remaining closed drainage lies in that arid land, where the breach gate's water balance
keeps basins closed, as it should. Synthetic relief also roughly doubles Earth's area of
closed depressions:

| share of land in a closed depression deeper than | 1 m | 100 m | 250 m | 500 m |
|---|---|---|---|---|
| Earth (~55 km) | 0.100 | 0.040 | 0.015 | 0.006 |
| synthetic relief | 0.20–0.24 | 0.09–0.11 | 0.032–0.035 | 0.004–0.010 |
| cratons-first | 0.07–0.12 | 0.012–0.038 | 0.001–0.016 | 0–0.008 |

Removing one synthetic-relief component at a time (seeds 910211954 and 447276134): interior
hills (shield, hill-country and plains noise) are the main source of depression area
(100 m share 0.087–0.105 → 0.032–0.065); orogens make the deep basins (250 m share
0.033–0.035 → 0.011); glacial pits and platform basins barely matter.

**Draining the synthetic surface.** The excess closed-depression area is not small pits. By
basin size, generated basins under 20 nodes cover 4.0–4.2% of land (Earth 4.3% under 20
cells), but basins of 20 or more nodes cover 15–20% (Earth 5.8%), and most spill below
500 m. High basins are, if anything, scarcer than on Earth (spill at 1 km or more: 0.5–1.8%,
Earth 3.2%). On Earth, 1.7% of land lies in 12 closed basins of 80,000 km² or more that
spill below 500 m (the largest ~413,000 km²).

`synthetic_relief.drain_land` therefore fills closed basins of at least 80,000 km² (~20
nodes at density 4) that spill below 500 m, by a priority flood from the world ocean with a
0.01 m/km drainage gradient (Barnes et al. 2014), except a seeded random set of them
covering about 1.7% of land (0.5–1.5 times that). It leaves high and small basins, and
basins around enclosed submerged holes, to the bootstrap's breach and water balance. The
fill is booked as mobile cover. Seeds and Earth as above; "before" is this branch with
`drain_land` stubbed out.

| density | surface | unfrozen closed drainage | lakes | depression > 1 m / > 100 m / > 250 m / > 500 m | lowland relief at 100 km (m) | land p50 (m) |
|---|---|---|---|---|---|---|
| 4x (4 seeds) | before | 0.48–0.64 | 0.060–0.086 | 0.19–0.24 / 0.09–0.11 / 0.032–0.035 / 0.004–0.010 | 148–159 | 350–482 |
| | after | 0.34–0.51 | 0.014–0.019 | 0.089–0.140 / 0.034–0.048 / 0.013–0.027 / 0.003–0.007 | 139–141 | 365–513 |
| 2x (2 seeds) | before | 0.67–0.70 | 0.073–0.076 | 0.17–0.21 / 0.06–0.11 / 0.010–0.042 / 0–0.006 | 139–148 | 394–396 |
| | after | 0.48–0.55 | 0.017–0.030 | 0.099–0.118 / 0.029–0.054 / 0.006–0.026 / 0–0.004 | 131–133 | 398–417 |
| 1x (2 seeds) | before | 0.62–0.70 | 0.066–0.067 | 0.14–0.18 / 0.06–0.09 / 0.011–0.028 / 0–0.005 | 160–171 | 389–440 |
| | after | 0.42–0.51 | 0.022 | 0.064–0.089 / 0.027–0.044 / 0.005–0.018 / 0–0.003 | 151–153 | 393–462 |
| Earth | | ~0.19 | ~0.02 | 0.100 / 0.040 / 0.015 / 0.006 | 166–173 | 423 |

- Depressions and lakes now sit near Earth at every density. The land fraction, the share
  above 2 km and the column backing are unchanged, and generation is no slower.
- Lowland relief falls by about 12%, because filled basins have flat floors.
- Unfrozen closed drainage stays above Earth's, mostly in the arid t = 0 climate's small dry
  pits and the kept large basins.
- The "before" depression rows were measured with an earlier `_drainage` that flooded along
  each node's own neighbour list rather than in reverse. On test spheres the two give the
  same shares.

After ten 100 kyr steps (the app's default step) on origin/main `34fb86f` and on this branch:

| seed | lakes | unfrozen closed drainage | depression > 1 m / > 100 m |
|---|---|---|---|
| 910211954, origin/main → branch | 0.093 → 0.039 | 0.26 → 0.23 | 0.197 / 0.082 → 0.141 / 0.045 |
| 447276134, origin/main → branch | 0.143 → 0.052 | 0.38 → 0.32 | 0.237 / 0.110 → 0.149 / 0.035 |

Lakes and shallow depressions grow over the first Myr of 100 kyr steps (from 1.4–1.9% and
9–14% at t = 0), as stepping roughens the surface.
