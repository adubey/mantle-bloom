# Continent generation: composite cratons and cratons at continental edges (#470)

A follow-up to the cratons-first generator (#462, PR #467) in
[RFC 0009](../../rfc/0009-craton-mosaic.md). #467's layout met most of the RFC's layout
targets, but its continents didn't look like Earth's:

1. **Every craton sat deep in its continent's interior.** Each continent grew outward from its
   cores by a uniform distance, so every craton was wrapped in an even rim of margin. On Earth
   many cratons sit at a continent's edge: West Africa, Australia's Pilbara and Yilgarn,
   Siberia, Baltica, and several in North America.
2. **Cratons looked like polka dots.** Every nucleus was its own craton, a noise-perturbed
   1,000 km disk at least 1,500 km from every other, so every world got 42-48 same-sized dots.

## Changes (`backend/app/continent_layout.py`)

- **Composite cratons, from a Chinese restaurant process.** Within a continent, each new
  nucleus starts a new craton with probability α / (α + n), where n is the continent's nuclei
  so far (`CRATON_CLUSTER_CONCENTRATION` α = 20). Otherwise it joins an existing craton, chosen
  in proportion to its nuclei, and is placed 0.8-1.3 core radii from one of that craton's
  nuclei (`NUCLEUS_JOIN_STEP_RANGE`), so the two disks overlap into one irregular core. The
  minimum spacing and the belt band now separate only different cratons, and a continent whose
  core target exceeds one disk still gets at least two cratons.
- **Cratons at continental edges.** A continent grows outward in order of growth time: the
  distance to its cratons plus the distance to its centre (`CONTINENT_CENTRE_PULL` = 1), over a
  growth speed exp(σ × noise) that varies smoothly at continent scale
  (`CONTINENT_GROWTH_SPEED_SIGMA` σ = 1), plus the outline noise. Cratons on a continent's far
  side, or where it grows slowly, end up near the coast.
- **Belt width 500 km** (`CRATON_BELT_WIDTH_KM`, the RFC's 200-600 km sweep range; it started
  at 300). With fewer, larger cratons there is less gap between them: at 300-400 km, belts fell
  to 7% of continental area on some seeds, under the RFC's 10-30%.
- **Belts close over the whole sphere.** With cratons at the coast, ocean beside a craton used
  to erode the closure and leave margin in the gap between two cratons. Ocean now dilates like
  any other point, and only the continent's own points become belt.
- **Shortfalls.** The centre pull can leave a continent's own region smaller than its area. The
  shortfall now comes from unclaimed points in the continent's growth order. Joined nuclei add
  less area each, so nucleus placement gets up to 1,000 tries (from 200).

## Figures

![Seed 7 before](before_seed7.png)
![Seed 7 after](after_seed7.png)

*Seed 7 at t = 0, `land_fraction` 0.29: #467's layout (top) and this one (bottom). Brown:
craton; tan: continental crust; teal: belt; blue: ocean; black: plate boundaries
(`bin/debug/render_layout_sheets.py --seeds 7`).*

All 40 seeds after the change:

![Seeds 7, 11, 2024 and 1-19 after](after_sheet1.png)

![Seeds 20-39 after](after_sheet2.png)

And before, as posted on #470: [seeds 7, 11, 2024 and 1-19](before_sheet1.png),
[seeds 20-39](before_sheet2.png).

## Measures

Over the same 40 seeds, `land_fraction` 0.29, density 1 (`measure_plate_layout.py
--cratons-first`; per-seed rows in [layout_t0.json](layout_t0.json)). The "before" craton
measures use the same code with α → ∞, no centre pull or speed field, and a 300 km belt.

| Measure (min / median / max) | Before | After |
|---|---|---|
| Cratons per world | 42 / 45 / 48 | 22 / 39 / 47 |
| Nuclei that joined an existing craton | 0% | 13 / 39 / 84% |
| Craton-area coefficient of variation | 0.14 / 0.16 / 0.19 | 0.39 / 0.59 / 1.29 |
| Cratons within 300 km of the continental-crust edge | 0 / 0 / 2% | 11 / 35 / 49% |
| Belts, of continental area (RFC 10-30%) | 11 / 13 / 18% | 13 / 17 / 20% |
| Craton, of continental area (RFC 35-50%) | 36 / 40 / 40% | 39 / 40 / 40% |
| Continental area (RFC 0.5 ± 0.05) | 0.50 | 0.50 |
| Continents with room for two cratons but fewer (RFC 0) | 0 | 0 |
| Cores crossed by a plate boundary (RFC 0) | 0 | 0 |

The plate-boundary targets move a little, because the boundaries now partition different
continents:

| Seeds passing, of 40 | Before (#467) | After |
|---|---|---|
| Land edge within 300 km of a boundary, 15-50% | 40 | 39 |
| Plates mixed by land ≥ 50% | 36 | 36 |
| Boundary with land on both sides ≤ 25% | 36 | 34 |
| All three | 33 | 29 |

Land on both sides is 15 / 22 / 27% (before 18 / 23 / 29%), and its interior part 10 / 17 / 25%;
which seeds lost a target and why isn't diagnosed. The craton count and edge share have no Earth target: Bleeker (2003) counts about
35 Archean cratonic pieces, and this layout makes 22-47 cratons. The knobs (α, the join step,
the centre pull, σ and the belt width) are design choices, set by this sweep over the 40 seeds:

| α, join step, centre pull, σ, belt width | Cratons | Joined | Near an edge | Belts |
|---|---|---|---|---|
| ∞, –, 0, 0, 300 km (before) | 40 / 44.5 / 48 | 0% | 0 / 0 / 2% | 10 / 14 / 17% |
| 2, 0.6-1.0, 0, 0, 300 km | 6 / 17 / 25 | 89% | 0% | 2 / 5 / 7% |
| 10, 0.6-1.0, 0, 0, 300 km | 18 / 30 / 39 | 60% | 0% | 5 / 9 / 13% |
| 20, 0.8-1.3, 0, 0, 300 km | 28 / 35 / 42 | 34% | 0% | 8 / 11 / 14% |
| ∞, –, 1, 0, 300 km | 40 / 44.5 / 48 | 0% | 0 / 11 / 29% | 10 / 14 / 17% |
| ∞, –, 1, 1, 300 km | 40 / 44.5 / 48 | 0% | 6 / 26 / 42% | 9 / 12 / 16% |
| 20, 0.8-1.3, 1, 1, 400 km | 20 / 37 / 45 | 35% | 7 / 31 / 59% | 7 / 13 / 16% |
| **20, 0.8-1.3, 1, 1, 500 km** | **22 / 38.5 / 47** | **39%** | **11 / 35 / 49%** | **13 / 17 / 20%** |

## Generation time

`generate_world` at density 1: 1.0-1.6 s over the 40 seeds (#467: 0.93-1.32 s), under the
RFC's 2 s.

## Reproduce

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_layout_sheets.py --out /tmp/sheets
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_plate_layout.py \
    --seeds 7,11,2024,1,2,3,4,5,6,8,9,10,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39 \
    --cratons-first --land-fraction 0.29 --out /tmp/layout.json
```
