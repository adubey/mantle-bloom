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
- Measures are per node; nodes are near-equal area.
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
| `synthetic_relief.py` (backed) | 348–438 | 0.17–0.23 | 0.057–0.116 | 2589–2647 | 0.32–0.34 | 147–162 | 0.063–0.105 |
| Earth | 423 | 0.26 | 0.107 | 1818 | 0.24 | 166–173 | ≈ 0.02 |

The implementation is fully backed: the largest |elevation − isostatic_elevation(Hc, Hm)| on
an ice-free node is 0 m, and no column is outside the Hc caps. It differs from B3 only in its
noise streams and in leaving oceanic crust as generated (B3 drew unbacked ridges and trenches).

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
| synthetic relief | 0 | 0 | 0.283 | 395 | 0.116 | 0.069 | 147 |
| | 1 | −35 | 0.291 | 423 | 0.117 | 0.064 | 150 |
| | 5 | −166 | 0.331 | 505 | 0.120 | 0.020 | 142 |
| | 10 | −271 | 0.340 | 608 | 0.133 | 0.018 | 130 |
| cratons-first control | 0 | 0 | 0.283 | 556 | 0.127 | 0.027 | 68 |
| | 1 | −36 | 0.287 | 585 | 0.128 | 0.024 | 71 |
| | 5 | −167 | 0.298 | 704 | 0.148 | 0.014 | 92 |
| | 10 | −267 | 0.302 | 807 | 0.164 | 0.012 | 107 |

- The world steps without a jump: belts persist, rivers organize, and collapse-mobile land
  stays at ~7%.
- Both worlds lose ~270 m of sea level in 10 Myr (the model's existing eustatic drift). With
  Earth-like shelves at 20–150 m, that drains the shelves: land grows by 5.7 points against
  1.9 for the control, and the 0–200 m share falls from 0.069 to 0.018. Realistic shelves make
  the existing sea-level drift visible quickly; an RFC should pair this relief with a
  sea-level budget that holds.
- Lowland relief erodes slowly (147 → 130 m), while the control roughens (68 → 107 m) as
  tectonics adds relief. The two converge, so the initial roughness is not lost in a step.
- Unbacked elevation grows equally in both runs (up to ~800 m on ~10% of nodes by 10 Myr),
  from stepping processes that move elevation without Hc (`back_elevation_gain_fields`
  documents this). It is not introduced by generation.
