# Water and ice scripts

These scripts measure rivers, lakes, depression breaching, glaciers and the ice load.
Hydrology, lakes and ice are described in
[simulation-model/hydrology.md](../simulation-model/hydrology.md) and
[simulation-model/surface-processes.md](../simulation-model/surface-processes.md). For lake and
basin structure in a single save, the
[`stranded_basins` and `lake_hierarchy_diagnostics` CLIs](offline-cli.md) are quicker.

Every script here takes `--backend`, so the same script measures an `origin/main` worktree
and a branch (see [Running the scripts](../debugging.md#running-the-scripts)). They come in
two kinds:

- **Save replays** (`measure_ice_sheet.py`, `measure_lake_silt.py`, `measure_river_endings.py`)
  share `replay_lib.py`'s command line: a save, `--backend`, `--myr` (required),
  `--step-years` (default 100 kyr), `--checkpoint-myr` (default 1) and `--out` (required).
  Each records step runtime and peak RSS.
- **Seed runs** (`measure_breaching.py`, `measure_ice_loading.py`,
  `measure_meltwater_faults.py`, `measure_seasonal_melt.py`) share `seed_run_lib.py`'s command
  line: `--backend`, `--seeds`, `--steps`, `--years` (step length) and `--every`, plus `--out`.

## Rivers and lakes

### `measure_river_endings.py`

Measures where rivers end (#503). Each step, from that step's hydrology, it records:

- water reaching no outlet (m³/s): inflow into warm, ice-free land sinks
  (`flow_target == -1`) with no visible lake, plus inflow into `flow_target` cycles;
- the number of flow cycles and the discharge into the ocean;
- drawn river segments (`_rivers_to_draw`) whose downstream node is neither drawn, ocean nor
  a visible lake;
- long-enough river networks (`_river_network_is_long_enough`) by mouth type;
- warm single-node ice plugs, and the ice ledger's terms.

Values are averaged over each checkpoint window.

```bash
backend/.venv/bin/python bin/debug/measure_river_endings.py world.mbworld --myr 5 --out /tmp/branch.json
backend/.venv/bin/python bin/debug/measure_river_endings.py world.mbworld --myr 5 \
    --backend ../mantle-bloom-baseline/backend --out /tmp/main.json
```

### `measure_discharge_calibration.py`

Replays river and glacier calibration on the same save in two code trees (#355). First
create a starting save per density with the baseline backend (`--generate-seed`,
`--node-density`), then replay it on each tree. It reports raw stream-power and ice-abrasion
capacity, the erosion budget, area-weighted channel width and depth percentiles with their
land coverage, runtime and peak RSS. `--river-erosion-coefficient` and
`--width-growth-coefficient` override those constants. Commit only a compact comparison to
`analysis/issue355/`.

```bash
backend/.venv/bin/python bin/debug/measure_discharge_calibration.py /tmp/d1.mbworld \
    --backend ../mantle-bloom-baseline/backend --generate-seed 355004 --node-density 1
backend/.venv/bin/python bin/debug/measure_discharge_calibration.py /tmp/d1.mbworld \
    --backend ../mantle-bloom-baseline/backend --myr 10 --checkpoint-myr 1 --out /tmp/main-d1.json
backend/.venv/bin/python bin/debug/measure_discharge_calibration.py /tmp/d1.mbworld \
    --myr 10 --checkpoint-myr 1 --out /tmp/branch-d1.json
```

### `measure_lake_silt.py`

Measures lakes and lake silt (#363). At each checkpoint: land fraction, lake count, lake area
and water volume, and silt laid on lake beds (`HydrologyFields.silt_deposited`, after each
lake's room clamp). It also records the silt budget before that clamp, from a tap on
`lakes._water_balance`: `legacy` (the model-unit rule used before #363) and `real_per_unit`
(real inflow × step years, the silt at a coefficient of 1). Their ratio calibrates the
real-unit coefficient. Each is split by whether the lake was spilling.

```bash
backend/.venv/bin/python bin/debug/measure_lake_silt.py world.mbworld --myr 10 --out /tmp/branch.json
```

### `measure_breaching.py`

Measures least-cost depression breaching (`app/breaching.py`, #297). For each seed it
generates the world twice, with `hydrology.BREACH_DEPRESSIONS_ENABLED` off and on
(`--modes`), and steps both identically. Every hydrology call records pits, breached pits,
pits too costly to carve, and endorheic pits kept closed by their water balance; notch depth
and count; boundary passes lowered below the higher centre; standing lakes (flooded nodes,
area, spilling area, volume); and the breaching pass's wall time. The `off` run computes the
breaching statistics without applying them.

```bash
backend/.venv/bin/python bin/debug/measure_breaching.py --seeds 1 2 --steps 50 --years 100000
```

## Ice

### `measure_ice_sheet.py`

Measures glacier ice on a save (#353). At each checkpoint: land fraction; the visible ice cap
(the land share at or above `GLACIER_VISIBLE_DEPTH_M`), ice volume, and mean and maximum
depth; glacial erosion, both raw capacity and applied; lake water frozen into the ice; the
ice ledger and its worst residual; and the share of step time in `ice_sheet.step_ice`.

```bash
backend/.venv/bin/python bin/debug/measure_ice_sheet.py world.mbworld --myr 10 --out /tmp/branch.json
```

### `measure_seasonal_melt.py`

Measures glacier extent, to compare seasonal melt against annual-mean melt (#275 phase 2).
After every step: the fraction of land under visible ice, total ice volume, and glaciated
land by annual-mean temperature band. Ice that survives only because the annual mean is below
freezing shows up in the −15 to 0 °C bands.

```bash
backend/.venv/bin/python bin/debug/measure_seasonal_melt.py --seeds 1 2 --steps 30
```

### `measure_ice_loading.py`

Measures the isostatic ice-load deflection (#275 phase 3): land fraction, ice volume, the
stored `ice_load_deflection_m` (deepest, mean over glaciated land, land share depressed by
more than 10 m) and `drift_m`, the largest gap between the stored deflection and the one the
current ice implies.

```bash
backend/.venv/bin/python bin/debug/measure_ice_loading.py --seeds 1 2 3 --steps 60
```

### `measure_meltwater_faults.py`

Multi-seed long runs of earthquakes triggered by ice unloading (#275 phase 5). Every
`--every` steps: land fraction, `hc_at_max_fraction`, craton survival and ledger balance, the
continental-material reservoirs, earthquakes by trigger, the largest unloading-triggered
magnitude, and the largest one-step ice loss as a stress. Seeds run in parallel
(`--workers`). `--compare` compares two `--out` files. Long runs diverge once any field
differs, so `--controlled` deep-copies the world every `--every` steps and steps each copy
`--controlled-steps` more with and without `faults.trigger_unloading_earthquakes`.

```bash
backend/.venv/bin/python bin/debug/measure_meltwater_faults.py --seeds 1 2 3 --steps 100 --out /tmp/branch.jsonl
backend/.venv/bin/python bin/debug/measure_meltwater_faults.py --compare /tmp/main.jsonl /tmp/branch.jsonl
backend/.venv/bin/python bin/debug/measure_meltwater_faults.py --controlled --seeds 1 2 3 --steps 80 --every 20
```
