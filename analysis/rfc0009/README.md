# RFC 0009 gate G1: does a craton mosaic give rifts a way around? (#426)

[RFC 0009](../../rfc/0009-craton-mosaic.md) proposes seeding each continent as a mosaic of
cratonic cores separated by belts that never become craton. Gate G1 asks whether the mosaic
lets RFC 0010's minimum cut route a rift between the cores. Its outcomes were fixed in the
RFC before this run.

**Outcome: neither building branch.** At the starting knobs (1,500 km core spacing, 300 km
belts), pooled over seeds 7, 11 and 2024, L = 24% and B = 0.91. The cheaper branch needs
L ≥ 20% *and* B ≤ 0.75. No seed passes it on its own. Band exposure is worse than on today's
layout with formation off (B = 0.77). Under the RFC's rule this is the falsification case:
don't build RFC 0010's cut. Belts also don't last. By 500 Myr they keep 5%, 19% and 34% of
their t = 0 area, under a third on two seeds, so under the RFC's other rule 0009 waits for
the sink RFCs (0001, 0002, 0005). Seed 2024's 34% means the stronger falsification, "every
seed under a third", is narrowly missed.

## Method

`bin/debug/probe_craton_mosaic.py` replays fresh worlds on origin/main `02d5145` (the probe
commit, `8ffdade`, adds only the script), with every default: 5 Myr steps to 1 Gyr and
100 Myr checkpoints. It runs one seed at a time and changes no saved state.

- **Mosaic seeding** (RFC 0009 solution 1). Per plate, the eligible columns are today's
  seeding host with the 500 km formation margin. Poisson-disc nuclei are placed
  `CORE_SPACING_KM` apart, from an RNG keyed on (seed, plate). Each column goes to its
  nearest nucleus by graph hops, within its connected piece of the interior. Columns within
  `BELT_WIDTH_KM / 2` of a partition boundary are belts. Everything else becomes craton, as
  seeding does today.
- **Belts** are carried by a probe-only HISTORY surface field, registered before the app
  loads. A merged or coarsened cell stays a belt if any contributor was one, which is the
  remap rule the RFC proposes.
- **Formation off** (`CRATON_FORMATION_MYR` = ∞), so the seeded layout is the only craton.
- **Every cut** that goes ahead (split or failed rift) is evaluated by
  `prototype_mincut_split._evaluate` at a 1,000 km corridor and floor 0.05, with all of RFC
  0010's admissibility gates. Strength is RFC 0010's: craton or keel, but
  `RIFT_BELT_STRENGTH` = 0 on a belt.
- **L**: the share of craton-crossing flow cuts that have an admissible minimum cut whose
  line touches no craton (no craton node at either end of a crossed link).
- **B**: craton area in the ±2.5-spacing band, using the minimum cut where it is
  admissible and the flow cut otherwise, divided by the same sum for the flow cuts.
- **Control** (`--layout blob`): today's 875 km seeding, also with formation off, so the
  two layouts differ only in where the seeded craton sits. RFC 0010's own control
  ([analysis/rfc0010](../rfc0010/README.md)) had formation on: L = 17% and B = 0.73.
- **Sensitivity** (`--belt-width-km 600`): the top of the RFC's belt-width sweep. This is
  not the gate.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py --out $S/g1_mosaic.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py --layout blob --out $S/g1_blob.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py --belt-width-km 600 --out $S/g1_mosaic600.json
```

Runs took 730–870 s each. Generation took 0.8–1.8 s with either layout. Compact numbers are
in [g1.json](g1.json).

## Results

### The gate (0–1 Gyr, 1,000 km corridor)

| Layout | Seed 7 L / B | Seed 11 L / B | Seed 2024 L / B | Pooled L / B | Craton-crossing cuts |
|---|---|---|---|---|---|
| **Mosaic, 300 km belts (gate)** | 9% / 0.94 | 37% / 0.83 | 33% / 0.98 | **24% / 0.91** | 53, 43, 27 |
| Today's layout, formation off (control) | 14% / 0.80 | 21% / 0.69 | 9% / 0.77 | 15% / 0.77 | 63, 33, 22 |
| Mosaic, 600 km belts (sensitivity) | 29% / 0.90 | 42% / 0.75 | 41% / 0.75 | 38% / 0.81 | 48, 59, 34 |

- **The mosaic raises L but not B.** Minimum cuts do run along belts. In one inspected seed 7
  cut, 37 of the 83 nodes on the cut were belt and 2 were craton. A clear line is also
  more common than on today's layout (24% against 15%). But the band a failed rift thins is
  ±312 km, wider than a 300 km belt. A cut down the middle of a belt still has cores on both
  sides of its band, so B stays near 1. Belts would have to be wider than the band, about
  625 km, before following one empties the band. Even at 600 km, B is 0.81 pooled, and only
  seeds 11 and 2024 reach 0.75.
- **The mosaic exposes more craton than today's layout.** It seeds with the 500 km
  formation margin, not 875 km, so it has more craton at t = 0 (125, 92 and 59 M km² against
  104, 82 and 46), and more of each plate is craton or belt.
- **Early cuts do better than late ones.** Over the first 300 Myr, the 300 km mosaic on
  seed 11 reaches L = 50% and B = 0.78. After 300 Myr the belts are mostly gone and both
  measures drift back toward the control (from the per-cut rows, which aren't checked in).

### Seeding (t = 0)

| | Seed 7 | Seed 11 | Seed 2024 |
|---|---|---|---|
| Seeded interior (M km²) | 150 | 111 | 71 |
| Belt share of it, 300 km / 600 km | 17% / 32% | 17% / 33% | 17% / 32% |
| Interior pieces with ≥ 2 nuclei, of which with ≥ 2 cores | 10 / 10 | 4 / 4 | 6 / 6 |
| Craton components (≥ 1 M km²), mosaic / today | 81 (61) / 23 (8) | 49 (37) / 10 (4) | 42 (28) / 12 (4) |

The seeding meets the RFC's t = 0 criteria: every piece large enough for two nuclei has at
least two cores, and belts are 10–30% of the seeded interior at the starting knobs.

### Belt survival (300 km belts)

Belt area counts belt-marked genuine continental columns without craton. "Marked on craton"
counts belt-marked cells that also carry craton, after a remap blended the two.

| Myr | Seed 7 belt (of t = 0) | Seed 11 | Seed 2024 | Seed 7 marked on craton (M km²) |
|---|---|---|---|---|
| 0 | 25.0 (100%) | 19.1 (100%) | 12.0 (100%) | 0.0 |
| 100 | 16.2 (65%) | 12.6 (66%) | 10.6 (89%) | 5.0 |
| 200 | 7.7 (31%) | 8.4 (44%) | 8.6 (72%) | 11.4 |
| 300 | 4.9 (20%) | 6.4 (34%) | 7.5 (63%) | 8.7 |
| 500 | 1.4 (5%) | 3.6 (19%) | 4.0 (34%) | 9.0 |
| 1,000 | 0.2 (1%) | 0.2 (1%) | 0.9 (8%) | 6.7 |

Genuine continental area falls too, by 69%, 43% and 33% at 500 Myr on the same runs. Belts
fall faster than that, by two routes:
- **Sinks.** Belts are non-cratonic host, so they are thinned out of the host window or
  destroyed like the rest of it (#409).
- **Merges onto craton.** When a quad merge or coarsening combines a belt cell with a
  craton cell, the result keeps both the mark (HISTORY) and part of the craton (extensive).
  On seed 7, 9 M km² of marked area sits on craton by 500 Myr. In the probe that cell still
  scores belt strength, but it counts as craton for L and B. The RFC's proposed remap rule
  ("keep a mark when any contributor has one") would make the same blend, and so would
  formation's exclusion of marked cells. The implementation would need to decide which
  state wins.

At 600 km, belts keep 7%, 15% and 26% of their area at 500 Myr: wider belts don't last
proportionally longer.

## What the outcome means

- **By the pre-registered rule**, G1 lands in the "Otherwise" branch: don't build RFC 0010's
  minimum cut. The belt-survival rule independently says the mosaic waits for the sink RFCs.
- **The falsification's stated reading** is "the flow puts regime boundaries through
  cratons". It fits L less well than B. L does rise with a mosaic (15% → 24%, or 38% with
  600 km belts), so the cut's geometry was part of the obstacle. B fails mainly because the
  operational band is wider than a belt. That is a property of how the model thins a failed
  rift (`FAILED_RIFT_BAND_MULT`), not of where the flow puts its boundaries.
- **Belts vanish within a few hundred Myr whatever their width.** No layout fix survives
  that. This matches #409's finding that non-cratonic host is what the model's sinks remove.

## What this does not show

- Whether a mosaic would pass G1 once the sinks are repaired (RFCs 0001, 0002, 0005), with
  belts added at collisions and rifts (solution 2's writers, not in this probe), or with
  formation on. The probe seeds belts only at generation.
- Whether B would fall with a failed-rift band scaled to the weak zone it follows. G1 fixed
  the band at today's `FAILED_RIFT_BAND_MULT`.
- How much rifting loss a mosaic saves after the split (#388). As in RFC 0010's probe, stretch
  donation and divergent thinning aren't followed.

# Layout: cratons follow plates (RFC 0009 revision)

G1's failure prompted mapping the worlds. Generation makes plates first, each wholly
continental or oceanic, and then puts a craton in the middle of each continental plate.
Earth works the other way round. This section is the evidence for RFC 0009's revision
(cratons first, plates placed around them).

## Figures

| | |
|---|---|
| ![Earth](earth_plates_cratons.jpg) | **Earth.** Red: the 16 major plates the model bundles (NNR-MORVEL56, `backend/app/data/plate_owner_grid.npz`). Orange and pink: shield and platform (cratons). USGS world geologic provinces, public domain. |
| ![Seed 7, t = 0](sim_seed7_0myr.png) | **Seed 7, t = 0, today's rules.** Brown craton, tan continental crust, blue oceanic crust, black plate boundaries. |
| ![Seed 7, t = 0, mosaic](sim_seed7_0myr_mosaic.png) | **Seed 7, t = 0, G1's mosaic** (teal belts). Cores are still framed by each plate's edges. |
| ![Seed 7, 500 Myr](sim_seed7_500myr.png) | **Seed 7, 500 Myr, today's rules.** Darker where plates overlap. |

```bash
# Earth: download the USGS map first (see the script's docstring).
PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_earth_cratons.py \
    --provinces provinces.jpg --out analysis/rfc0009/earth_plates_cratons.jpg
PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py --seed 7 --out sim_seed7_0myr.png
PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py --seed 7 --mosaic --out sim_seed7_0myr_mosaic.png
PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py --seed 7 --years 5e8 \
    --save-world seed7-500myr.mbworld --out sim_seed7_500myr.png
```

The simulation maps are on a plain lat/lon grid (origin/main `02d5145`). Where plates
overlap, the nearest cell alternates between them, so the script draws each pixel's
majority plate over a 9-pixel window and shades overlap. The Earth map's projection isn't
documented. It fits a Mercator projection to within 8 pixels at six coastline landmarks.

## Measures

`bin/debug/measure_plate_layout.py` on a quarter-degree grid ([plate_layout.json](plate_layout.json)).
For Earth it uses the 16 major plates and the real coastline: land stands in for continental
crust (leaving out the drowned shelves).

| | Earth | Seed 7 | Seed 11 | Seed 2024 | Seed 7 at 500 Myr |
|---|---|---|---|---|---|
| Continent edge within 300 km of a plate boundary | 28% | 100% | 100% | 100% | 64%* |
| Plates with mixed crust (5–95% continental) | 69% | 0% | 0% | 0% | 31% |
| Boundary length with continent on both sides | 17% | 33% | 20% | 27% | 28%* |
| Continental share of the surface | 29% (land) | 63% | 40% | 49% | 16% |
| Craton area within 500 km of a plate boundary | – | 0% | 0% | 0% | 89%* |

\* Inflated by the overlap zones, where the nearest cell's plate flips from pixel to pixel.
The 500 Myr row is from the saved world above, which isn't checked in.

- **At t = 0 every continental edge is a plate boundary, and no plate is mixed.** On Earth,
  most continental margins are passive, inside a plate, and most plates carry both continent
  and ocean. The premade "earth" world can't represent that either: it splits each real plate
  into land and ocean sub-plates (65 in all), because a generated plate has one crust type.
- **Cratons are placed by the tiling.** No craton lies within 500 km of a plate boundary,
  because seeding measures 875 km in from each plate's edge.
- **The runtime does make mixed plates.** By 500 Myr, 31% of seed 7's plates are mixed, as
  oceans open inside continental plates.
- **Earth's continents aren't free of plate boundaries.** 17% of Earth's boundary length has
  land on both sides: the Alpine–Himalayan belt, Siberia, the Middle East and East Africa.
  These run along young orogens and rifts, not through the old cratons.

# RFC 0009 gate G2: do mixed plates work? (#438)

The cratons-first generator makes most plates mixed (continental and oceanic crust) from
t = 0. Gate G2 asks whether the runtime handles that before the generator is built. Its
metrics and decision rule were fixed in the RFC before this run.

**Outcome: fails as written, on one metric, and the diagnosis points at the layout, not at
mixed plates.** Genuine continental area at 500 Myr is 3.9–4.5 points of the sphere higher
in the mixed runs than in the control on seeds 11 and 2024 (noise floors 1.8 and 1.2). Plate
count, overlap, merges and splits pass, and the ledgers behave the same as in the control.
The gap opens in the first 100 Myr, while the control's continents are still edged by plate
boundaries and the mixed runs' are not. After that the two lose continent at the same rate:
the same rule applied to the 100→500 Myr change passes on every seed. The extra crust the
control loses is mostly booked to deep subduction at trenches. Re-tiling moves continental
margins inside plates, the effect RFC 0009 wants, and the gate can't tell that apart from a
runtime defect. The project owner accepted this reading (see "Decision" below).

## Method

`bin/debug/probe_mixed_plates.py` runs fresh worlds on origin/main `652c2c0` (the probe adds
only the script), with every default: 5 Myr steps to 500 Myr and 100 Myr checkpoints. It runs
one world at a time and saves nothing. Per seed (7, 11 and 2024):

- **control**: `generate_world(seed)`, today's single-type plates.
- **mixed A and B**: the same finished world, re-tiled. Each draws a new
  `build_plate_tiling`, with the same plate count, from its own RNG stream, keyed on
  (seed, 438, run). Each new plate's cells copy every surface field from the nearest control
  node, so the crust, the seeded cratons and the ledger tracer are the same. A plate's nominal
  `crust_type` is the type covering the majority of its area. Its minority-type cells carry an
  explicit `crust_type_code`, as the RFC's generator will build them. The world then runs
  `finish_generation`, with the copied craton volume booked as `initial_m3`.
- **Metrics** (RFC, Design, G2):
  - genuine continental area (Hc ≥ `CRATON_HOST_MIN_HC_M`, `compare_cratons`' measure) and
    plate overlap (`stats`), as percentages of the sphere;
  - active plates;
  - cumulative merges and splits, from the logged events;
  - continental-material and craton ledger balance errors.
- **Diagnostics:** the share of plates mixed by continental area (5–95%), each
  continental-material ledger account, and `measure_plate_layout`'s share of continent edge
  within 300 km of a plate boundary.
- **Decision rule:** noise floor = max(|mixed A − mixed B|, floor), with floors of 1 plate,
  1 event or 1 percentage point. A metric fails when |mean(mixed) − control| > 2 × noise
  floor on at least two seeds. Plate count and overlap are judged at each checkpoint, the
  others at 500 Myr.

```bash
PYTHONPATH=backend:bin/debug backend/.venv/bin/python bin/debug/probe_mixed_plates.py --out $S/g2.json
PYTHONPATH=backend:bin/debug backend/.venv/bin/python bin/debug/probe_mixed_plates.py --summarize $S/g2.json
```

Each run took 290–440 s, and generation with re-tiling 1.4–4.2 s. Rerunning gave identical
numbers. Compact results are in [g2.json](g2.json).

## Results

At t = 0 the re-tiled runs have 60–100% of plates mixed (the control has none). Genuine
continental area matches the control to within 0.06 points, and craton area to within
0.3 M km².

### The gate (500 Myr)

| Metric | Seed 7: control / A / B | Seed 11 | Seed 2024 | Seeds failing |
|---|---|---|---|---|
| **Genuine continental area, % of sphere** | 17.6 / 17.3 / 21.1 | **15.8 / 19.5 / 21.2** | **14.6 / 19.1 / 18.0** | **2 (fails)** |
| Merges | 24 / 29 / 21 | 6 / 8 / 7 | 7 / 2 / 4 | 0 |
| Splits | 52 / 59 / 51 | 56 / 44 / 52 | 58 / 51 / 41 | 0 |
| Active plates at 500 Myr | 48 / 49 / 47 | 50 / 45 / 53 | 50 / 54 / 52 | 0 |
| Overlap at 500 Myr, % of sphere | 23.9 / 16.2 / 14.9 | 17.4 / 17.7 / 15.1 | 19.4 / 14.0 / 9.8 | 1 |

Active plates and overlap at 100–400 Myr fail on at most one seed each (see
`summary.metrics` in g2.json). Seed 11's failing gap is 4.5 points against a 1.8 floor, and
seed 2024's is 3.9 against 1.2. On seed 7 the two mixed runs differ by 3.7 points.

**Ledgers.** The craton ledger closes to 1e-15 in every run. The continental-material ledger
drifts by up to 1.9e-3 of its initial volume, in the control as much as in the mixed runs
(1.9e-3 on seed 2024's control, at most 1.0e-3 in a mixed run). The drift is on origin/main
already, most likely strike-slip shear (#299), so it doesn't separate the layouts.

### Diagnosis: the starting layout, not the runtime

| | Seed 7: control / A / B | Seed 11 | Seed 2024 |
|---|---|---|---|
| Continent edge within 300 km of a boundary, t = 0 | 100 / 26 / 20% | 100 / 18 / 14% | 100 / 15 / 14% |
| Same, 100 Myr | 71 / 59 / 65% | 91 / 49 / 67% | 80 / 49 / 52% |
| Genuine area change 0→100 Myr, points | −12.9 / −10.1 / −10.2 | −7.0 / −3.2 / −3.2 | −0.6 / +1.1 / +1.8 |
| Genuine area change 100→500 Myr, points | −20.6 / −23.5 / −19.8 | −8.9 / −9.0 / −7.3 | −9.0 / −6.1 / −8.0 |
| Deeply subducted by 500 Myr, M km³ | 3,890 / 2,904 / 4,386 | 5,611 / 3,538 / 3,015 | 2,537 / 1,180 / 2,040 |
| Juvenile additions by 500 Myr, M km³ | 2,536 / 1,478 / 2,005 | 3,350 / 1,902 / 2,051 | 1,926 / 977 / 1,395 |

- **Most of the gap opens in the first 100 Myr.** Re-tiling leaves only 14–26% of the
  continent's edge near a boundary, so most margins start passive. The runtime makes the
  control's plates mixed too (24–48% by 100 Myr), and by 500 Myr the two layouts expose
  their continents about equally (59–85%).
- **After 100 Myr the loss rates agree.** The RFC's rule applied to the 100→500 Myr change
  in genuine area passes on all three seeds: |mean mixed − control| is 1.0, 0.7 and 1.9
  points against noise floors of 3.7, 1.7 and 1.9.
- **The control's extra loss goes to the trench sink.** `deeply_subducted_m3`, the
  continental margin consumed at trenches, is the largest account in every run. It is
  2,000–2,600 M km³ larger in seed 11's control than in its mixed runs, and 500–1,400
  larger on seed 2024. That is offset in part by more arc crust (juvenile additions), which
  also needs an active margin.
- **The call sites can't be ruled out.** The audit (posted on #438) lists 22 runtime
  reads that apply a plate's majority type to all its nodes. Several of them suppress arcs or
  change the retreat rule on a mixed plate's minority margin (`torque.py:308`'s
  plate-level `neighbor_is_oceanic`, `lithosphere_plate.continental_arc_band`). They
  probably add to the lower juvenile additions in the mixed runs. The loss rate after
  100 Myr agrees within noise anyway, while 22–53% of the control's plates are mixed too.

## Decision

**Accepted by the project owner on 2026-10-15 (option 1 below):** the cratons-first
generator goes ahead on mixed plates, and the 22 call sites are fixed in #447 before G3.

The rule as written compares two layouts that differ in margin exposure by design, so a
failure on continental area can't separate "mixed plates are mishandled" from "passive
margins lose less crust". Rerunning without changes gives the same numbers. The options for
the project owner:

1. **Accept the diagnosis.** Record G2 as failing on continental area, attributed to the
   starting layout, with the 100→500 Myr comparison passing. Proceed with the cratons-first
   generator on mixed plates, and fix the per-node call sites (at least the boundary
   classification in `torque.py:308` and the arc band) as the tracking issue's "call-site
   fixes" step.
2. **Fix the call sites first and rerun G2.** This is the RFC's literal route. It is
   unlikely to close the gap, since the gap comes from exposure.
3. **Use the single-type fallback** (RFC, Design, G2).
