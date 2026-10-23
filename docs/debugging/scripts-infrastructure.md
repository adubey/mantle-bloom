# Surface, save, sweep and profiling scripts

These scripts check the plate surface lattice, convert and replay saves, run parameter
sweeps and profile animation. This page also lists the shared libraries the other scripts
import.

## Surface lattice and saves

### `surface_parity.py`

Plate-surface audit runs and their quality gates (#247, #249). `audit` runs every seed (in
parallel with `--jobs`) and judges them: `comparison.json` holds every gate result and
`report.md` a summary, and the exit status is 1 when the verdict is `fail`. `run` does one
seed's world, optionally continuing a save (`--from-world`). `judge` re-judges an existing
output directory. Presets (`surface_parity.PRESETS`) are `smoke`, `standard`, `long` and
`issue147`, and any preset field can be overridden. `--render` writes elevation and
platesDetail PNGs at each checkpoint. Gates and presets are documented in
[surface-parity.md](../surface-parity.md).

```bash
backend/.venv/bin/python bin/debug/surface_parity.py audit --preset smoke --seeds 3 --out analysis/parity-smoke
```

### `measure_quad_area_churn.py`

Per-step churn in the area a quad surface represents, and the sea-level jitter it drives
(#259). It steps one world and after every step sums `SurfaceNodes.area_m2` over all plates
and reads `world.sea_level_m`. It then reports the standard deviation of each per-step change
over the second half of the run, and their correlation.

```bash
backend/.venv/bin/python bin/debug/measure_quad_area_churn.py --seed 1 --steps 40
```

### `convert_legacy_saves.py`

Converts line-backed `.mbworld` saves to sparse quads and judges the result (#248). For each
save it converts through `legacy_conversion`, writes the `ConversionReport`, round-trips the
converted world through the current save format, runs the hard-invariant audit
(`surface_parity.audit_world`), checks a second round trip reproduces its state hash, steps it
`--steps` times and checks the tolerances in
[save-compatibility.md](../save-compatibility.md). It writes `<save>.json` per save plus
`summary.json` and `summary.md`. `--write-converted` also writes each converted world as
`<save>-quad.mbworld`. `--recheck` re-judges existing output without converting.

```bash
backend/.venv/bin/python bin/debug/convert_legacy_saves.py SAVE.mbworld \
    --out analysis/issue248/conversions [--steps 3] [--render] [--jobs 4]
```

### `replay_logged_event.py`

Replays one step recorded by the server's replay log and checks that it reproduces (#442).
See [Replay logs](replay-logs.md) for how records are made and read.

```bash
backend/.venv/bin/python bin/debug/replay_logged_event.py \
    ~/.mantle-bloom/replay-logs/8442/<record> [--save-after /tmp/after.mbworld]
```

## Parameter sweeps

The sweep tools come from #171, which looked at land fraction trending low and ice-cap
fraction trending high on long runs. Raw and aggregated output goes to `bin/debug/results/`,
which is git-ignored.

### `run_sweep.py`

Runs `sweep_lib.run_one_job` for every (parameter, multiplier, seed) combination and appends
one JSON line per (job, checkpoint) to `--out`. It resumes: a triple that already has every
checkpoint in `--out` is skipped. `--params` defaults to `sweep_lib.DEFAULT_SWEEP_PARAMS`,
which leaves out four parameters with a repeated null result (#199); name one to recheck it.
A run with a different `--node-density` or `--checkpoints` should get its own `--out` file.

```bash
backend/.venv/bin/python bin/debug/run_sweep.py
backend/.venv/bin/python bin/debug/run_sweep.py --params rain_erosion,volcanism --workers 4
backend/.venv/bin/python bin/debug/run_sweep.py --node-density 2.0 --checkpoints 30,60,90,120,150,180 \
    --params avg_rotation_rate,volcanism --out bin/debug/results/sweep_results_density2.jsonl
```

### `aggregate_sweep.py`

Aggregates `run_sweep.py`'s rows into a mean and standard deviation per (parameter,
multiplier, checkpoint) across the sweep seeds, folding the shared baseline rows in as each
parameter's multiplier 1.0 point. It writes one small JSON file.

```bash
backend/.venv/bin/python bin/debug/aggregate_sweep.py \
    [--in bin/debug/results/sweep_results.jsonl] [--out bin/debug/results/sweep_summary.json]
```

## Profiling

### `profile_animation_60m.py`

Profiles a 60-frame, 60 Myr animation: seed 0 at node density 4, climate density 4, fluid
density 2, 10 steps of 100 kyr per frame, encoded to MP4 the way the animation route
does. It writes per-frame timings of each major step phase, render and encode to
`analysis/issue147-profile-20260922/` (`frames.csv`, `metadata.json`, `animation.mp4` and cProfile output).
It takes no options, and the output path is fixed. See [profiling.md](../profiling.md) for
profiling in general.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/profile_animation_60m.py
```

## Shared libraries

These modules have no command line of their own; scripts import them.

### `seed_run_lib.py`

The `--backend --seeds --steps --years --every` options and the generate-and-step loop for
`measure_breaching.py`, `measure_cap_overflow.py`, `measure_collision_denudation.py`,
`measure_ice_loading.py`, `measure_meltwater_faults.py`, `measure_orogenic_relief.py` and
`measure_seasonal_melt.py`. `import_backend` puts `--backend` first on `sys.path`, so a
script imports `app` only after calling it.

### `replay_lib.py`

The save-replay command line (`save`, `--backend`, `--myr`, `--step-years`,
`--checkpoint-myr`, `--out`, plus a script's own options through `add_arguments`) and the
step-and-checkpoint loop for `measure_ice_sheet.py`, `measure_lake_silt.py`,
`measure_river_endings.py` and `measure_river_hops.py`. `load_replay` puts `--backend` first
on `sys.path`.

### `sweep_lib.py`

The sweep's parameters (`PARAM_SPECS`), seeds (`SWEEP_SEEDS`, ten seeds), checkpoints
(30, 60, 90 and 120 Myr), step size (10 Myr) and node density (1.0), and `run_one_job`. Each
job starts from a fresh `generate_world`, never from a save. The elevation scripts on
[Elevation and relief scripts](scripts-elevation.md) and `rotation_ledger.py` reuse these
settings.

### `polarity_replay_lib.py`

Documented with the [collision polarity scripts](scripts-collisions.md#polarity_replay_libpy).
