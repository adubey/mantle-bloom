# Elevation and relief scripts

These scripts check whether elevation is backed by its crustal column, how often land
saturates at `MAX_ELEVATION_M`, and how Earth-like a world's relief is. The isostasy model is
described in [simulation-model/foundations.md](../simulation-model/foundations.md).

Except for `measure_initial_relief.py`, these scripts come from the #171 sweep investigation.
They reuse [`sweep_lib.py`](scripts-infrastructure.md#sweep_libpy)'s seeds (`SWEEP_SEEDS`),
node density 1 and 10 Myr steps, so their numbers compare with the sweep's.

### `measure_initial_relief.py`

Compares a world's relief with Earth's (#218). For each generated seed or `--load`ed save it
reports land hypsometry (percentiles, share above 1 and 2 km); ocean depth (percentiles, IQR,
shares shallower than 3 km and 200 m, and shelf water within 300 km of land); enclosed seas
and lake-covered land; RMS relief between land points at 50–1,600 km, overall and for
lowland pairs; and column backing (the largest |elevation − isostatic elevation| on ice-free
nodes, columns outside the Hc/Hm caps, and land that starts collapse-mobile or with a hot
Moho). Earth references are in the script's `EARTH` table. `--step-myr` steps each world and
measures again at each listed time. `--png` renders the elevation view of each measurement.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_initial_relief.py \
    --seeds 447276134,373936961 --cratons-first --synthetic-relief --out /tmp/relief.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_initial_relief.py \
    --seeds 447276134 --cratons-first --synthetic-relief --step-myr 1,5,20 --png /tmp/maps \
    --out /tmp/relief-steps.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_initial_relief.py \
    --seeds "" --load world.mbworld --out /tmp/relief-save.json
```

### `diagnose_elevation_saturation.py`

How widespread land saturation at `MAX_ELEVATION_M` is, and for pinned nodes, whether Hc is
also at `MAX_CRUSTAL_THICKNESS_M` (#180). If Hc is at its cap, the column ran out of room. If
elevation clips before Hc or Hm reach their caps, the isostasy formula or its calibration is
the suspect. Takes `sweep_lib`'s `--seeds`, `--checkpoints` and `--node-density`, plus
`--collision-uplift-multiplier` and `--workers`. Writes to
`bin/debug/results/elevation_saturation.jsonl`.

```bash
backend/.venv/bin/python bin/debug/diagnose_elevation_saturation.py
backend/.venv/bin/python bin/debug/diagnose_elevation_saturation.py --seeds 829071382 --workers 1
```

### `drill_down_isostasy.py`

For every land node pinned at `MAX_ELEVATION_M`, recomputes
`lithosphere.isostatic_elevation(hc, hm, rho_c)` from the node's own column (#189). If the
recomputed value is also at the cap, the formula reads out the cap. If it is well under, the
stored elevation is "isostatic debt" written by a path that bypassed Hc/Hm. Runs the
`SWEEP_SEEDS` × collision-uplift grid {1.0, 0.25} to 120 Myr. No options.

```bash
backend/.venv/bin/python bin/debug/drill_down_isostasy.py
```

### `measure_transform_far_field_debt.py`

Measures the isostatic debt (elevation − `isostatic_elevation`) across all land at several
checkpoints of a long run (#189 follow-up). `transform_uplift` adds unbacked relief by
design, and it is the only remaining source of debt on columns with Hc/Hm tracking, so debt
growth reads how much that term injects. No options.

```bash
backend/.venv/bin/python bin/debug/measure_transform_far_field_debt.py
```

### `drill_down_volcanism.py`

Single-seed, step-by-step drill-down for #173 ("more volcanism → less land"). It runs seed
829071382 at `volcanism_multiplier` 1.0 and 4.0 side by side. Every step it records
continental crust volume, land volume above sea level, the elevation volcanism added, and net
erosion on `is_volcano` nodes against everywhere else. No options.

```bash
backend/.venv/bin/python bin/debug/drill_down_volcanism.py
```

### `measure_volcano_node_elevation.py`

Checks whether nodes that have erupted reach the elevation `ERUPTION_ELEVATION_M` and
`VOLCANIC_PLAIN_ELEVATION_M` were tuned for (#198, following #178). `is_volcano` marks every
node that was ever part of a volcanic field, mostly subsided oceanic crust, so the script
selects `mineral_deposit_m > 0` (the node has erupted at least once) and tracks active
volcanoes alongside. No options.

```bash
backend/.venv/bin/python bin/debug/measure_volcano_node_elevation.py
```
