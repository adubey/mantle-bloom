# Craton layout scripts

These scripts measure and draw where plates, continents, cratons and belts sit relative to
each other. Most were written for [RFC 0009](../../rfc/0009-craton-mosaic.md), the
cratons-first layout
([simulation-model/foundations.md](../simulation-model/foundations.md#cratons-first-layout-continent_layoutpy-plate_partitionpy)). The run-time craton ledger is covered by the
[craton scripts](scripts-cratons.md).

Run them from the repository root. Scripts whose usage shows `PYTHONPATH=backend` need it.
See [Running the scripts](../debugging.md#running-the-scripts) for baselines.

## Measuring a layout

### `measure_plate_layout.py`

Measures where a world's plate boundaries run relative to its continents. On Earth, plate
boundaries mostly run through ocean. On a quarter-degree grid (each pixel takes its nearest
node) the script reports:

- boundary length and the share with land on both sides, on one side or on neither (for
  generated worlds, also continental crust on both sides);
- the share of the continental edge and of the coastline within `EDGE_NEAR_KM` of a boundary;
- mixed plates: the share of plates, and of plate area, that are 5–95% continental;
- craton area within 500 and 1,000 km of a boundary, and boundary length with craton on both
  sides;
- land-on-both-sides boundary length, split into continental interior and coast-adjacent
  (and, for a cratons-first world, along a belt).

`--cratons-first` generates with RFC 0009's layout (#462) and adds core, belt and margin
shares and generation time. `--land-fraction` passes the Generate World dialog's land slider.
`--earth` measures Earth's 16 real plates (`data/plate_owner_grid.npz`) against the real
coastline (`data/real_land_mask.npz`). `--load` measures a save.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_plate_layout.py \
    --seeds 7,11,2024 --earth --out /tmp/layout.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_plate_layout.py \
    --seeds 7,11,2024 --cratons-first --land-fraction 0.29 --out /tmp/layout-cf.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_plate_layout.py \
    --seeds "" --load seed7-500myr.mbworld --out /tmp/layout-500.json
```

### `probe_mixed_plates.py`

RFC 0009 gate G2 (#438): do mixed plates (continental and oceanic crust on one plate) behave
like single-type plates? Per seed it runs three worlds from the same starting crust, in 5 Myr
steps: `control` (`generate_world(seed)`), and `mixed A` and `mixed B`, the same world
re-tiled with a new plate tiling drawn from its own RNG stream. Each new plate copies its
surface fields from the nearest control node. At each checkpoint it records genuine
continental area, plate overlap, active plates, the mixed-plate share, merges and splits, the
ledger balance errors, each continental-material account, and the continental edge within
300 km of a boundary. `--summarize` applies the RFC's decision rule. A metric fails when
|mean(mixed) − control| > 2 × the noise floor on at least two seeds, where the noise floor is
max(|A − B|, floor).

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_mixed_plates.py --out /tmp/g2.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_mixed_plates.py \
    --seeds 7 --years 1e7 --every 5e6 --out /tmp/g2-smoke.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_mixed_plates.py --summarize /tmp/g2.json
```

### `compare_g2.py`

Compares two `probe_mixed_plates.py --out` files, usually `origin/main` against a branch
(#447). For each seed and run it writes, at the last checkpoint: genuine continental area,
juvenile additions, deep subduction, active plates, merges, splits and the worst ledger
balance error. Per seed and model it adds the control-versus-mixed gaps, with G2's noise
floor. Plain Python; no backend needed.

```bash
python3 bin/debug/compare_g2.py MAIN.json BRANCH.json --out analysis/issue447/g2_ab.json
```

### `probe_craton_mosaic.py`

RFC 0009 gate G3 (#472): does the cratons-first layout keep cratons whole? It compares the
production cratons-first generator with today's generator (`--layout blob`) on the same seeds
and step schedule, with craton formation off. Every split or failed rift is evaluated with
`prototype_mincut_split.py`'s minimum cut (1,000 km corridor, floor 0.05) under RFC 0010's
gates. Over the cuts whose band crosses craton, it reports L (line clearance) and B (band
exposure). Each 100 Myr checkpoint records belt area as a share of genuine continental area
and of its t = 0 value, craton cores per plate, the ledger balance and cut counters.
`--years 0` measures generation only. `--core-spacing-km` and `--belt-width-km` sweep the
layout's knobs.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py --out /tmp/mosaic.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py --layout blob --out /tmp/blob.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py --years 0 --out /tmp/mosaic-t0.json
```

### `measure_belt_marks.py`

Replay A/B for RFC 0009's belt marks (#461). `run` generates a default world per seed and
writes a row every `--every` years: belt area, craton area and volume, formed volume, genuine
continental area, both ledger balances, and the formation audit (craton formed on
belt-marked columns, which must stay 0). `bench` times `--steps` steps at `--density` after
`--warmup` steps, `--repeats` times, with the share spent in belt-mark code. `summarize`
compares two runs. Each row records `model_file`; the script imports its own tree's backend,
so a baseline needs a copy of the script inside the baseline worktree.

```bash
backend/.venv/bin/python bin/debug/measure_belt_marks.py run --seeds 7,11,2024 \
    --years 1e9 --every 100e6 --out /tmp/belts-branch.jsonl
backend/.venv/bin/python bin/debug/measure_belt_marks.py bench --seed 7 --density 4 --out /tmp/bench-branch.jsonl
backend/.venv/bin/python bin/debug/measure_belt_marks.py summarize /tmp/belts-main.jsonl /tmp/belts-branch.jsonl
```

## Maps

### `render_craton_map.py`

Draws a plain equirectangular map with four classes: oceanic crust, continental crust, craton
and plate boundaries (black, where the majority plate changes between pixels). Overlapping
plates are shaded darker (RFC 0009, #426). Elevation and climate are left out, because the
map only shows where cratons sit relative to plates. Options:

- `--years` steps the world first, in 5 Myr steps;
- `--load` draws a save instead of generating one, and `--save-world` saves the stepped
  world for a later `--load`;
- `--mosaic` seeds RFC 0009's mosaic (formation off) and draws belts;
- `--cratons-first` generates with the cratons-first layout;
- `--land-boundaries` draws coastlines and colours land-on-both-sides boundaries by their
  distance from a coast.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py --seed 7 --out /tmp/seed7.png
PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py \
    --seed 7 --cratons-first --land-fraction 0.29 --out /tmp/seed7-cratons-first.png
PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_craton_map.py --load world.mbworld --out /tmp/saved.png
```

### `render_layout_sheets.py`

Contact sheets of cratons-first layouts at t = 0, one map per seed (#470). Each seed is drawn
with `render_craton_map.render` and titled with its plate, continent and craton counts. The
per-seed PNGs and 4 × 5 sheets go to `--out`.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_layout_sheets.py \
    --seeds 7,11,2024 --land-fraction 0.29 --out /tmp/sheets
```

### `layout_explorer.py`

An interactive explorer for the cratons-first layout's tuning constants (#685). It serves one
page, `bin/debug/layout_explorer.html`, on localhost. Choose a port other than 8000 and 5173,
which belong to the user's own instance.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/layout_explorer.py --port 8765
```

**Knobs.** The sidebar lists the module constants in the script's `KNOBS` table:
- Continents: count range, Poisson-disc spacing, Dirichlet size concentration, outline and
  growth noise.
- Cratons: the craton-clustering concentration, craton share, nucleus steps, size scale,
  core radius and spacing.
- Belts: belt width, closure radius and lens.
- Plates: point sites, candidate partitions, refinement, boundary costs, auto plate count.

Each knob starts at its value on the checkout, so a branch's edits show as its defaults.
World settings (node density, plate count, continental and land fraction) are passed to
`lithosphere_plate.cratons_first_layout`.

**Jobs.** Each map is one job: a seed plus constant overrides. A job runs
`cratons_first_layout` in a process pool (`--workers`, default 4), with the overrides
applied to the modules for that job and reset before the next one. A job builds the layout
only, without plates or relief: about 1 s per layout at node density 1 and about 5 s at 4.
A layout that raises, such as the partition failure in #510, is shown as an error tile.

**Batches.** Each Generate adds a column of maps to the right of the earlier ones:
- **N seeds** at fixed knobs.
- **Sweep across columns**: one column per value of a knob, each with the same N seeds, so
  rows line up.
- **Sweep down one column**: one knob stepped over N values at a fixed seed.

A (min, max) knob such as `CONTINENT_COUNT_RANGE` is swept one of three ways: as a fixed
value (`k`–`k`, which sets the continent count itself), or by its min or max end with the
other end held at the sidebar's value. A value that would put min above max is rejected
before the batch runs.

A column's header lists its changed knobs, with the median and range of the layout
diagnostics (`ContinentLayout.diagnostics`). Stop cancels the jobs that haven't started.

**Display.** Maps can be coloured by geology, land and water, plate, craton or continent, with
overlays for plate boundaries, coastline, craton nuclei and plate sites. Click a map for an
enlarged view and all its diagnostics.

A knob changes the random draws that follow it. The same seed under different settings is
therefore a different world from that point on, not a slightly altered one. To compare
settings, use several seeds per column rather than single maps.

### `render_earth_cratons.py`

A reference figure for RFC 0009: Earth's 16 bundled plate outlines (NNR-MORVEL56, see
`real_plates.py`) drawn in red over the USGS world geologic provinces map, whose shield and
platform provinces are the cratons. The USGS map is public domain but not bundled; download it
first. The script fits a Mercator projection to it from six coastline landmarks.

```bash
curl -sSLo provinces.jpg https://upload.wikimedia.org/wikipedia/commons/a/a9/World_geologic_provinces.jpg
PYTHONPATH=backend backend/.venv/bin/python bin/debug/render_earth_cratons.py \
    --provinces provinces.jpg --out /tmp/earth_plates_cratons.jpg
```
