# Cratons-first generator: layout at t = 0 (#462)

[RFC 0009](../../rfc/0009-craton-mosaic.md)'s step "Cratons-first generator" (Design,
"Generation, cratons first (solution 1)"). This measures the layout the new generator builds
at t = 0 against the RFC's acceptance targets, over seeds 7, 11, 2024 and 17 more (1-6, 8-10,
12-19). The generator is behind a flag (`continent_layout.CRATONS_FIRST_DEFAULT = False`), so
`main`'s random worlds are unchanged: with the flag off, `generate_plates` on seeds 7, 11 and
2024 (with and without the UI's sliders) and `generate_world(7)`'s seeded cratons hash the
same as origin/main `dd3c620`.

The model description is in
[docs/simulation-model/foundations.md](../../docs/simulation-model/foundations.md#cratons-first-layout).

## Figures

![Seed 7 at t = 0, cratons-first layout](seed7_0myr_cratons_first.png)

*Seed 7 at t = 0, cratons-first layout, `land_fraction` 0.29, density 1. Brown: craton; tan:
continental crust (belt and margin); blue: oceanic crust; black: plate boundaries
(`bin/debug/render_craton_map.py --seed 7 --cratons-first --land-fraction 0.29`). Five
continents carry 47 cores. Boundaries reach the continents along belts between cores and
never cut one.*

![Seed 7 at t = 0, plates first](seed7_0myr.png)

*The same seed generated plates first (today's default), for comparison: every continental
edge is a plate boundary and each continental plate's craton sits in its middle.*

## Layout measures

`bin/debug/measure_plate_layout.py --cratons-first` on a quarter-degree grid; per-seed rows in
[layout_t0.json](layout_t0.json). Land uses the Generate World dialog's default
`land_fraction` 0.29; `generate_world(seed)` passes none, and that row is below.

| Measure (min / median / max over 20 seeds) | Target | Cratons first | Plates first | Seeds passing |
|---|---|---|---|---|
| Land edge within 300 km of a plate boundary | 15-50% (Earth 28%) | 16 / 26 / 33% | 45 / 73 / 100% | 20 / 20 |
| Plates mixed by land share | ≥ 50% (Earth 69%) | 50 / 72 / 92% | 0 / 23 / 56% | 20 / 20 |
| Boundary length with land on both sides | ≤ 25% (Earth 17%) | 21 / 27 / 32% | 8 / 15 / 22% | **4 / 20** |
| Boundary length with craton on both sides (a core crossed) | 0 | 0 | 0 | 20 / 20 |
| Continents with room for two cores but fewer than two | 0 | 0 | – | 20 / 20 |
| Belts, of continental area | 10-30% | 10.8 / 12.8 / 18.2% | – | 20 / 20 |
| Craton, of continental area | 35-50% | 36.3 / 40.0 / 40.3% | 19 / 27 / 49% | 20 / 20 |
| Continental area, of the sphere | 0.5 ± 0.05 | 0.499 / 0.500 / 0.501 | 0.19 / 0.37 / 0.63 | 20 / 20 |

Reported separately, as the RFC asks: the continental-crust edge within 300 km of a boundary
is 18 / 25 / 36% (plates first: 100% on every seed), and plates mixed by crust type are
50 / 77 / 90% (plates first: 0%). Boundary length with continental crust on both sides is
40 / 47 / 60%.

Craton area is below 40% on one seed (36.3%) because a continent smaller than one core disk
gets no core. The "plates first" craton share is the 875 km interior rule's.

Without a `land_fraction`, land is wherever the relief noise is high, independent of the
continents: land covers 7.5-47% of the sphere, and two of the land measures swing with it
(mixed by land 18 / 63 / 100%, 14 of 20 seeds passing; land on both sides 4 / 31 / 44%, 7 of
20 passing). The crust-based measures don't depend on land.

### The failing target: land on both sides of a boundary

Boundaries with land on both sides run 21-32% of boundary length, over the ≤ 25% target on
16 seeds. What was tried, each over the 20 seeds at `land_fraction` 0.29:

| Change | Land on both sides (min / median / max) | Effect on the others |
|---|---|---|
| `PLATE_COST_MARGIN` 3 → 6, 10 or 20 | 21 / 27-28 / 32% | none |
| Refinement 0 sweeps instead of 20 | 24 / 30 / 35% | land edge 14-27% |
| Refinement 200 sweeps | 21 / 27 / 32% | same as 20 |
| `CONTINENT_COUNT_RANGE` 3-8, 4-10 or 1-3 | 18-23 / 27-28 / 31-34% | mixed by land falls to 41-42% on some seeds |
| `CONTINENT_SIZE_CONCENTRATION` 5 | 25 / 27 / 36% | mixed by land min 42% |
| Dijkstra weighted by point cost, margin cheaper than ocean (0.33) | 15 / 23 / 29% | mixed by land min 43%, land edge 15-28% |

Refinement lowers J by 15-20% and moves land on both sides from 30% to 27%, but it converges
on a local minimum. A boundary across a continent is straight in its interior, where moving one
cell either way doesn't lower J, so no boundary cost gets it out. What sets the number is the
partition itself: half the sphere is continental, and Poisson-disc sites give plates of similar
size that the ±25% area band keeps similar, so a 100-190 M km² continent has to be divided
between several plates. The boundaries that do so run along belts, as the RFC intends, but
they still have land on both sides. Earth's plates are far less equal (the Pacific plate is
about a fifth of the sphere). Plates first scores better (15%) only because there every
continental edge is a plate boundary.

None of the RFC's knobs reaches the target without breaking the mixed-plate target. Whether
to relax the target, change the partition (unequal plate sizes, a non-local optimiser), or
carry this into G3 as it stands is the owner's call.

## Generation time

`generate_world` on seed 7 including the layout and plate build: 0.93-1.18 s at density 1
over the 20 seeds (plates first 0.84-1.13 s), and 3.1 s at density 4 (plates first 1.8 s). The
first generation in a process also builds the cached whole-sphere lattice; the slowest
density-1 run above (1.18 s, seed 7) is that one.

## Judgment calls

- **Drowned cores keep host thickness.** Land and sea come from the relief noise, which
  doesn't know where the continents are, and offshore subsidence then thins submerged
  continental crust. That left 26-31% of core area thinner than `CRATON_HOST_MIN_HC_M` on
  seeds 7, 11 and 2024 (up to 78% with no `land_fraction`), which seeding would have made into
  craton on thinned shelf. Core columns are raised to at least the host minimum before
  seeding, so drowned cores are submerged platforms (on the order of 1 km deep).
- **`BELT_CLOSURE_KM = 800`** (not an RFC knob): at 600 km belts were 9.7-16% of
  continental area.
- **A continent whose core target exceeds one disk always gets two nuclei.** Before, the
  noise let one core stretch over a target just above a disk (seeds 4 and 14).
- **`PLATE_COST_CORE` = `PLATE_COST_MARGIN`.** The RFC gives costs for ocean, belt and
  margin only; a boundary along a core's edge is costed like one through margin.
- **`continental_fraction` stays a plate share with the flag off.** `CONTINENTAL_FRACTION`
  (0.4) and the UI copy are unchanged so `main` doesn't move; the area share is
  `plates.CONTINENTAL_AREA_FRACTION = 0.5`, used only by the cratons-first path. The slider
  copy ("of the surface") changes when G3 makes the layout the default.

## Reproduce

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_plate_layout.py \
    --seeds 7,11,2024,1,2,3,4,5,6,8,9,10,12,13,14,15,16,17,18,19 \
    --cratons-first --land-fraction 0.29 --out /tmp/layout-cf.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_plate_layout.py \
    --seeds 7 --cratons-first --land-fraction 0.29 --density 4 --out /tmp/layout-cf-d4.json
```
