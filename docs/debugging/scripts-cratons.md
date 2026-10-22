# Craton scripts

These scripts measure where craton volume goes: the craton ledger over a whole run, the site
that books each loss, and the loss caused by single merges and rift cuts. The ledger itself
and the Cratons map view are described in [Debug map views](map-views.md#cratons-craton),
and the model in
[simulation-model/tectonics.md](../simulation-model/tectonics.md#cratons-cratonspy). Losses
at plate contacts (jams, collision budgets) are measured by the
[collision scripts](scripts-collisions.md). The scripts for RFC 0009's cratons-first layout
are on [Craton layout scripts](scripts-craton-layout.md).

Run them from the repository root with `backend/.venv/bin/python`. Each one imports the
backend of its own checkout; see
[Running the scripts](../debugging.md#running-the-scripts) before using one for an
`origin/main` baseline.

## Long runs and ledger history

### `compare_cratons.py`

Long multi-seed comparison of worlds with cratons on and off (#274). Each (seed, variant)
job generates a default quad world and steps it to `--years`, with a checkpoint every
`--every` years. Each checkpoint records genuine continental area
(Hc >= `cratons.CRATON_HOST_MIN_HC_M`), node land fraction, the largest plate's genuine area,
craton extent and ledger, and tectonic activity: plate count, merges, splits, rift events and
mean rotation rate. The `off` variant disables craton seeding and formation, so every craton
resistance is inert. It also reports rifting by path per Gyr and the share of split and
failed-rift cuts that pass through craton ([RFC 0003](../../rfc/0003-rifting-around-cratons.md)).
`--set CRATON_NAME=VALUE` overrides a `cratons` constant. `--save-dir` writes each job's final
world as `seed<S>-<variant>-<Myr>myr.mbworld` for `ledger_history.py`.

```bash
backend/.venv/bin/python bin/debug/compare_cratons.py --seeds 7,11,2024 --years 500e6 \
    --out /tmp/cratons.jsonl --workers 6
backend/.venv/bin/python bin/debug/compare_cratons.py --summarize /tmp/cratons.jsonl
```

Several other scripts import `compare_cratons` for its job runner and measures:
`rebaseline_rfcs.py`, `craton_genuine_area.py` and `probe_mixed_plates.py`.

### `ledger_history.py`

Reads a save's `World.stats_history` (one row per stats sample, about 10 Myr apart on long
runs). It needs no replay. It prints three tables:

- checkpoints: craton volume and area, genuine continental area, land fraction and each
  craton sink's running total;
- craton destruction per Gyr by sink, as a share of the craton volume at the window start
  (rifting also split by path), with each sink's largest single-sample jumps. A sink made of
  a few large jumps comes from one-off events;
- continental-material fluxes per window in km³/yr: juvenile additions, remelted returns and
  the terminal sinks, against the mean live inventory.

```bash
backend/.venv/bin/python bin/debug/ledger_history.py world.mbworld [--json out.json]
```

### `rebaseline_rfcs.py`

Measures the RFC 0001–0007 baselines once per seed, serially (#403). It calls
`compare_cratons.run_job` unchanged (default quad worlds, cratons on, 5 Myr steps, 100 Myr
checkpoints), adding read-only taps for continental budgets, each collision merge's craton
loss, and trench attrition's share of `subducted_m3`. Saves go to `--save-dir` as scratch;
only the compact JSON belongs in `analysis/issue403/`.

```bash
backend/.venv/bin/python bin/debug/rebaseline_rfcs.py --out /tmp/tuning.json --save-dir /tmp/issue403
backend/.venv/bin/python bin/debug/rebaseline_rfcs.py --seeds 456862899 --years 3e9 \
    --out /tmp/long.json --save-dir /tmp/issue403
```

### `craton_genuine_area.py`

Replays #403's runs and adds `craton_genuine_area_km2`, the cratonic area inside the genuine
continental set, at each checkpoint (#409). `cratons.diagnostics` counts every cell with
craton crust, including cells thinned below `CRATON_HOST_MIN_HC_M`, so genuine minus craton
area is only a lower bound on non-cratonic host. The script wraps `compare_cratons._snapshot`
and runs `compare_cratons.run_job` as `rebaseline_rfcs.py` does, so its trajectory can be
checked against #403's JSON.

`--model-root` names the checkout whose model and `compare_cratons.py` are run. Both are
imported from there, never from this script's tree. Seeds run one at a time.

```bash
$BASE/backend/.venv/bin/python bin/debug/craton_genuine_area.py --model-root $BASE \
    --seeds 7,11,2024 --years 1e9 --out /tmp/genuine.jsonl
```

### `rfc0004_gates.py`

Computes [RFC 0004](../../rfc/0004-proto-cratons.md)'s formation and craton-share gates from
checked-in JSON (#409): `analysis/issue403/tuning.json`, `long.json` and
`analysis/issue409/craton_genuine_area.jsonl`. It needs no save or replay. Per run it reports
formation timing (the pulse at `CRATON_FORMATION_MYR` and formation after it), inventories
at each 100 Myr checkpoint, an upper bound on the host area formation converted against the
net loss of host, and each 500 Myr window's qualifiers. `--generation` also generates each
seed and measures the quiet host band that seeding leaves uncratonised, which needs the
backend.

```bash
python3 bin/debug/rfc0004_gates.py [--json analysis/issue409/gates.json]
PYTHONPATH=backend backend/.venv/bin/python bin/debug/rfc0004_gates.py --generation \
    --json analysis/issue409/gates.json
```

## Attributing a loss

### `split_craton_rework.py`

Replays a save and splits every craton-ledger sink by the call site that booked it (#332),
by tapping `cratons.record`. `collision_reworked_m3` is split into polarized transfer
(`crust_transfer._book`), own-plate suture accretion, anatexis relamination, orogen
relaxation, and phase-audit clips and residuals (`cratons.PhaseAudit.settle`). Checkpoints
also record the craton inventory and area, the mean relief of non-cratonic crust near craton
edges (all edges, and edges facing a convergent front), land and emergent continental area,
hypsometry, orogenic `ceiling_overflow_m3`, mean plate speed and merges so far.

With `--baseline` (this script's JSON from `main`, same save and span) it checks #340's
targets: the non-attrition craton loss rate, at most `NON_ATTRITION_LOSS_MAX_FRACTION` of
the baseline's, and craton-edge relief at convergent fronts, above the baseline's.
`--rescore` recomputes the targets from finished replays without stepping. Sites a branch
doesn't have don't appear.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/split_craton_rework.py \
    world.mbworld --myr 20.2 --out split.json [--baseline split-main.json]
PYTHONPATH=backend backend/.venv/bin/python bin/debug/split_craton_rework.py \
    --rescore split.json --baseline split-main.json
```

### `probe_craton_clips.py`

Finds the step function that leaves craton outside its column (#340). A phase audit clips
craton above its column's Hc or on a non-continental cell and books it to a guessed account,
so a clip in `split_craton_rework.py`'s output doesn't name its writer. This script replays a
save to `--from-myr`, then steps to `--to-myr` with taps on every audit's settle and on the
topology functions, printing each one that changes the craton outside columns (`excess`) or
the ledger by more than 1 km³. `--snapshot` saves the world at `--from-myr` and reuses it on
the next run. `--seed` generates a world instead of loading a save.

`--capture-dir` saves every collision merge with [RFC 0001](../../rfc/0001-merge-cap-craton-loss.md)'s
signature (more than `--min-share` of live craton booked, plus continental material clipped
above its column) as a pickle of the world just before the merge, for
`measure_merge_placement.py`. `--contacts` writes RFC 0008's per-step contact rows
(`probe_craton_contacts.Probe`) for the history before each captured merge.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_clips.py \
    world.mbworld --from-myr 238.4 --to-myr 240.45 --snapshot /tmp/world-238.mbworld
backend/.venv/bin/python bin/debug/probe_craton_clips.py --seed 7 --to-myr 400 \
    --step-years 5e6 --capture-dir /tmp/rfc0001 --contacts /tmp/rfc0001/contacts-7.jsonl
```

A [replay-log](replay-logs.md) record's `before.mbworld` is an ordinary save, and
`replay_logged_event.py` prints the matching `probe_craton_clips.py` command.

### `measure_merge_placement.py`

Replays one merge captured by `probe_craton_clips.py --capture-dir` (or
`probe_collision_budget.py --capture-dir`) and attributes its craton loss (RFC 0001 gate,
#440). It runs the real `merge_split._merge_plates_with_craton_audit` with read-only taps,
checks the replay books the same craton as the run, and reports each of RFC 0001's
hypotheses: room under the suture cap (H1), the eligible mask (H2), the volume each stage of
`_place_suture_crust` placed, craton carried with the overflow (H4), and excess left by later
writers (H3). `--contacts` adds the pair's RFC 0008 contact rows over `--history-myr` before
the merge.

```bash
backend/.venv/bin/python bin/debug/measure_merge_placement.py /tmp/rfc0001/seed7-*.pkl \
    --contacts /tmp/rfc0001/contacts-7.jsonl --json /tmp/rfc0001/placement-7.json
```

## Rift cuts

### `attribute_rift_boundary_loss.py`

Attributes divergent and stretch craton loss to the split cut that made the boundary
([RFC 0003](../../rfc/0003-rifting-around-cratons.md) phase 1b, #388). It replays a seed from
year 0 with probe-only surface fields. At each split, the daughters' nodes within the cut's
band are marked "through" if the band held craton, otherwise "clear". Later divergent
(`lithosphere_plate.deform_columns`) and stretch (`quad_tectonics._open_rift`) craton loss is
booked by those marks. Two sensitivity schemes run alongside: `ever` (a "through" mark is
never overwritten) and `wide` (twice the band). The summary reports F, the share of
`rifted_m3` that steering the cut could remove.

```bash
backend/.venv/bin/python bin/debug/attribute_rift_boundary_loss.py --seed 7 --years 1e9 \
    --out /tmp/issue388-seed7.jsonl
backend/.venv/bin/python bin/debug/attribute_rift_boundary_loss.py --summarize /tmp/issue388-seed*.jsonl
```

### `prototype_mincut_split.py`

Prototypes a split cut that bends around strong lithosphere as a graph minimum cut (#389,
RFC 0003; see [RFC 0010](../../rfc/0010-rifts-follow-weak-lithosphere.md)). For every cut in an
ordinary replay (a read-only tap on `cratons.record_rift_cut`, so the trajectory is
unchanged), it builds the plate's k-nearest-neighbour graph and finds the minimum cut between
the nodes beyond a corridor around the flow's great circle. Edge cost is length times (floor +
mean strength), with strength from RFC 0003's cut score. Each cut is checked against RFC
0010's admissibility gates. The output gives craton exposure in the failed-rift band and on
the cut line, and the boundary length. `--coverage-only` skips the replay and reports the
generated world's craton coverage.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/prototype_mincut_split.py --out /tmp/mincut.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/prototype_mincut_split.py \
    --seeds 7 --years 2e8 --out /tmp/mincut-smoke.json
```

`probe_craton_mosaic.py` reuses this minimum cut.
