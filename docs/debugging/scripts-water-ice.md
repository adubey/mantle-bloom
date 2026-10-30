# Water and ice scripts

These scripts measure rivers, lakes, depression breaching, glaciers and the ice load.
Hydrology, lakes and ice are described in
[simulation-model/hydrology.md](../simulation-model/hydrology.md) and
[simulation-model/surface-processes.md](../simulation-model/surface-processes.md). For lake and
basin structure in a single save, the
[`stranded_basins` and `lake_hierarchy_diagnostics` CLIs](offline-cli.md) are quicker.

Most replay and measurement scripts take `--backend`, so the same script measures an
`origin/main` worktree and a branch (see [Running the scripts](../debugging.md#running-the-scripts)).
The RFC 0014 G5 cap prototype uses its own checkout's backend and restricts inputs to the
preregistered save checksums. Save replays and seed runs share these two command-line patterns:

- **Save replays** (`measure_ice_sheet.py`, `measure_lake_silt.py`, `measure_river_endings.py`,
  `measure_river_hops.py`, `trace_lake_drainage_phases.py`)
  share `replay_lib.py`'s command line: a save, `--backend`, `--myr` (required),
  `--step-years` (default 100 kyr), `--checkpoint-myr` (default 1) and `--out` (required).
  Each records step runtime and peak RSS.
- **Fixed-forcing solver controls** (`gate_lake_partition.py`) capture one registered #675
  event checkpoint and compare the lake solve under 10, 25 and 100 kyr partitions. Terrain,
  hierarchy, outlets and climate forcing remain fixed. The script requires the repository's
  backend and takes `deep` or `shallow`, a save, and `--out`.
- **Seed runs** (`measure_breaching.py`, `measure_ice_loading.py`,
  `measure_meltwater_faults.py`, `measure_seasonal_melt.py`) share `seed_run_lib.py`'s command
  line: `--backend`, `--seeds`, `--steps`, `--years` (step length) and `--every`, plus `--out`.

## Rivers and lakes

### `gate_lake_depth_caps.py`

Replays the RFC 0014 G4 diagnosis saves and records each lake- or sea-tier artificial depth-cap
binding by active body and step. It reports the body's leaf-sink identity, area, cap level,
direct external supply and positive overflow surplus; rim-limited spills and floor drains do
not count as depth-cap events. Commit only the compact diagnosis result in `analysis/issue553/`.

```bash
backend/.venv/bin/python bin/debug/gate_lake_depth_caps.py world.mbworld --myr 1 \
    --step-years 100000 --out /tmp/seed.json
```

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

### `measure_river_hops.py`

Counts river lines drawn across lakes (#520). At the loaded save and at each checkpoint, it
records spill hops (flow edges between non-neighbour nodes), the median flow-edge length, and,
for River Inspector segments and drawn map segments, how many are longer than 3× that median
and how many join non-neighbours. A backend without `hydrology.flow_segments` is measured from
each river node's own `flow_target` edge. `--render-dir` writes an elevation globe centred on
`--center LAT LON` at each checkpoint. `--timing N` adds the median of N cold (spill-path
cache cleared) and warm timings of an elevation map render, the `/worlds/{world_id}/rivers` payload and
`flow_segments` over every river node.

```bash
backend/.venv/bin/python bin/debug/measure_river_hops.py world.mbworld --myr 0.1 --checkpoint-myr 0.1 \
    --center 43 153 --render-dir /tmp/branch --out /tmp/branch.json
backend/.venv/bin/python bin/debug/measure_river_hops.py world.mbworld --myr 0.1 --checkpoint-myr 0.1 \
    --backend ../mantle-bloom-baseline/backend --center 43 153 --render-dir /tmp/main --out /tmp/main.json
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

### `trace_lake_drainage_phases.py`

Replays a saved lake-flicker event and maps the modeled drainage catchment in selected phases.
For each phase, it selects governing lake hierarchy bodies whose wet members overlap a
connected component containing the fixed geographic anchors. The output reports each matched
component's volume and wet area, the selected bodies' total member area, and the dry-land
catchment area. The optional PNG colors catchment land orange and connected lake water blue.
These are model hierarchy catchments, not independently measured watershed divides.

```bash
backend/.venv/bin/python bin/debug/trace_lake_drainage_phases.py world.mbworld \
    --anchors analysis/issue675/shallow-lake-anchor.json \
    --phase pre-drop=6 --phase low=7 --phase rebound=8 \
    --myr 0.8 --step-years 100000 --out /tmp/shallow-catchments.npz \
    --figure-out /tmp/shallow-catchments.png
```

### `gate_lake_partition.py`

Runs the RFC 0014 G7 fixed-forcing partition control on the registered deep step-5 or shallow
step-8 #675 event. Each case replays to a common checkpoint, captures the fixed terrain,
hierarchy/outlets, liquid and ice stores, and one 100 kyr interval of lake and ice climate, then
compares 10, 25 and 100 kyr integrations. It writes a compact JSON result and records endpoint
storage, anchor connectivity, component count, merge/split timing and ledger closure. The
registered held-out cases are for this attribution/robustness gate only; they must not be used
to select mechanism ranges for G8.

```bash
NUMBA_CACHE_DIR=/tmp/rfc0014-g7-numba-cache backend/.venv/bin/python \
  bin/debug/gate_lake_partition.py deep world.mbworld --out /tmp/deep.json
NUMBA_CACHE_DIR=/tmp/rfc0014-g7-numba-cache backend/.venv/bin/python \
  bin/debug/gate_lake_partition.py shallow world.mbworld --out /tmp/shallow.json
```

### `gate_lake_endpoint.py`

Refreshes RFC 0014 G1 (#589) on the three registered diagnosis saves under the current D2(c)
lake-ice sequence. It captures one production 10 kyr step's C liquid forcing and subsequent
ice-phase inputs, then compares endpoint equilibrium against an independently substepped direct
C integration for isolated leaf bodies (1,000 steps by default; `--transient-substeps 4000` for
a convergence check). Every body above 10,000 km² is measured, plus a seeded 5% sample of smaller
bodies. Merged bodies, leaves that could receive hierarchy transfers, spilling bodies, invalid
and unsampled bodies count as fallback. It records liquid endpoint error (the G1 criterion), a
secondary combined liquid-plus-ice comparison, relaxation and fill times, starting-volume
coverage and the capture step runtime. Direct offline measurements do not form a production
eligibility cache.

```bash
backend/.venv/bin/python bin/debug/gate_lake_endpoint.py \
    ~/Downloads/mantle-bloom-seed910211954-220200000y.mbworld \
    --out analysis/issue589/seed910211954.json
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

### `gate_lake_ice_cap.py`

Runs RFC 0014 gate G5 (#541) on only its three preregistered diagnosis saves. It applies a
closed-form Stefan conduction and loss balance for 10 Myr under fixed saved climate and lake
depth, reports cap thickness, frozen-to-bed area, volume and sublimation, and crosses the
0.25x/1x/4x sublimation and 0.5x/1x/2x conduction sensitivities. Sublimation scales with
`1 − f_free` for each cell. It rejects saves outside the
registered SHA-256 set. Commit only its compact JSON in `analysis/issue541/`.

```bash
backend/.venv/bin/python bin/debug/gate_lake_ice_cap.py \
    ~/Downloads/mantle-bloom-seed910211954-220200000y.mbworld \
    --out analysis/issue541/seed910211954.json
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

## Climate

### `measure_climate_fields.py`

Measures the climate fields on a save in either tree (#522). At each checkpoint it computes a
fresh climate snapshot and records, weighted by cos(latitude): land and ocean air temperature,
humidity (g/kg; an `origin/main` tree's unitless index is converted by `HUMIDITY_INDEX_G_KG`)
and precipitation; land's share of precipitation; the same by 15-degree band; Köppen group
shares of land; land fraction and the ice cap. At t = 0 both trees see the same loaded world,
so a difference there comes from the climate code alone.

```bash
backend/.venv/bin/python bin/debug/measure_climate_fields.py world.mbworld --myr 10 --step-years 250000 \
  --checkpoint-myr 2 --backend ../<baseline-worktree>/backend --out /tmp/main.json
```
