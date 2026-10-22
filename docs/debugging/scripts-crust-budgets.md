# Crust and mantle budget scripts

These scripts account for continental crust (Hc), mantle lithosphere (Hm) and the
continental-material ledger: where collisions, retreat, topology clean-up, erosion and the
thickness caps lose or move material. For a per-phase Hc/Hm budget over a short interval,
start with the [`phase_budget_diagnostics` CLI](offline-cli.md#python--m-appphase_budget_diagnostics----per-phase-hchm-budget).
The ledgers are described in
[simulation-model/tectonics.md](../simulation-model/tectonics.md).

Scripts with a `--backend` option import that backend, so one script can measure both an
`origin/main` worktree and a branch. See
[Running the scripts](../debugging.md#running-the-scripts).

## Hm and crust transfer

### `attribute_hm_ratchet.py`

Attributes the climb of mantle lithosphere toward the Hm cap to step phases (#310). It
replays a save with `World.debug_diagnostics` on, so `phase_budget` books every instrumented
phase's Hm volume change, and adds what that budget can't see:

- `convergent_clip`: Hm that `rheology.apply_convergent_deformation` grew and then clipped at
  `MAX_MANTLE_LITHOSPHERE_THICKNESS_M`;
- `unattributed`: the whole-world Hm change minus every recorded phase's change;
- cap saturation at each checkpoint (global, cratonic and non-cratonic, per plate);
- land, emergent continental area, hypsometry and #314's collisional-shortening ledger
  (`phase_budget.SHORTENING_PHASE`).

`--compare` prints the same metrics for other saves, such as the checkpoint the issue's
numbers came from.

```bash
backend/.venv/bin/python bin/debug/attribute_hm_ratchet.py world.mbworld --myr 43.1 \
    --compare later.mbworld --out analysis/issue310/hm-budget.json
```

### `compare_crust_transfer.py`

Compares `attribute_hm_ratchet.py` replays of the same save and span, for example with and
without crust transfer (#320). For each run at the final checkpoint it reports land,
continental and emergent continental area, hypsometry, the continental-material inventory,
collision loss by account (total and per Myr), the continental Hm-cap fraction, runtime and,
where the engine has it, crust-transfer fronts, volumes and cost. `--baseline314` adds #314's
checked-in replay. Plain Python; no backend needed.

```bash
python3 bin/debug/compare_crust_transfer.py main=base.json transfer=transfer.json \
    --baseline314 analysis/issue314/cascade-seed997271774.json
```

### `summarize_issue317_suture_pairs.py`

Renders the suture Hm totals recorded for #317 (`hm_suture_budget_km3.by_pair` at the last
checkpoint) as a Markdown table by donor and neighbour pair. It reads
`analysis/issue317/cascade-after-316-seed997271774.json` and writes
`analysis/issue317/suture-pair-breakdown-after-316.md` by default (`--input`, `--output`).

```bash
python3 bin/debug/summarize_issue317_suture_pairs.py
```

## Collision sinks

### `calibrate_collision_sinks.py`

How much continental crust collisions lose, per step and compounded, under the current
constants or overrides (#276). It steps generated worlds (`--seeds`) or a save (`--save`) and
appends one JSON line per snapshot with two rates:

- **retreat-processed loss**: the share of the Hc that suture retreat donated
  (`suture_donated_m3`) lost to delamination and the no-outlet remainder;
- **whole-inventory loss**: the collision sinks as a share of the live inventory, compounded
  step by step (`collision_survival`), and the same for every sink (`all_sink_survival`), with
  a 500 Myr projection.

It also records land and continental-host area, the continental Hc distribution, area at the
Hc cap, the orogenic relief budget and the convergent ceiling overflow. `--set NAME=VALUE`
overrides a module constant in every `app.*` module that holds it. `--label` tags the rows.

```bash
backend/.venv/bin/python bin/debug/calibrate_collision_sinks.py --seeds 1 2 3 --steps 100 \
    --out bin/debug/results/276.jsonl
backend/.venv/bin/python bin/debug/calibrate_collision_sinks.py --save world.mbworld --steps 50 \
    --set SUTURE_ACCRETION_MAX_DELAMINATION_FRACTION=0.05 --label frac05 --out bin/debug/results/276.jsonl
```

### `attribute_suture_delamination.py`

Attributes every `delaminated_lower_crust_m3` increment and every suture no-outlet remainder
to the call site and front that produced it (#290). Generates `--seed` and steps it `--steps`
times of `--years`.

```bash
backend/.venv/bin/python bin/debug/attribute_suture_delamination.py --seed 3 --steps 100
```

### `attribute_topology_removed.py`

Classifies what is booked as `topology_removed_m3` (#305). It steps generated worlds and
classifies every continental-material booking by `merge_split.defragment_plates` or
`remove_defunct_plates` by the piece dropped: whether it carries continental crust, its node
count and its distance in node spacings to the nearest live plate. Prints one JSON summary
per seed.

```bash
backend/.venv/bin/python bin/debug/attribute_topology_removed.py --seeds 1 2 3 --steps 200
```

### `probe_quad_retreat_continental.py`

Splits the continental crust removed by quad boundary retreat (`_retreat`) by what happens
to it (#249): `accreted_kept` (put back on surviving cells), `accretion_cap_loss` (clipped by
`SUTURE_ACCRETION_MAX_HC_M`), `subducted_continental_plate` and `subducted_oceanic_plate`.

```bash
backend/.venv/bin/python bin/debug/probe_quad_retreat_continental.py --seed 1 --steps 240
```

## Thickness caps and erosion

These four scripts share `seed_run_lib.py`'s command line: `--backend`, `--seeds`,
`--steps`, `--years` (step length, default 1 Myr) and `--every` (snapshot interval in
steps), plus `--out`.

### `measure_cap_overflow.py`

Measures depositional Hc-cap overflow and where it ends up (#288). Every `--every` steps it
records Hc-cap saturation and thickening indicators (mean continental Hc, the area share at
or above `--belt-fraction` of the cap), erosion's overflow diagnostics by depositional
pathway (`--per-pathway` prints them), and every continental-material loss account, kept apart
so depositional overflow can be told from suture delamination and rift thinning.

```bash
backend/.venv/bin/python bin/debug/measure_cap_overflow.py --seeds 1 2 3 --steps 100
```

### `measure_collision_denudation.py`

Measures Hc-cap saturation and collision-belt denudation (#275 phase 4): `hc_at_max_fraction`
by node count and by area, land fraction, the ledger reservoirs, and for capped columns the
share that is land, their median slope and the share erosion didn't touch. Erosion totals per
window include collision-belt removal and `hc_cap_overflow_m3`.

```bash
backend/.venv/bin/python bin/debug/measure_collision_denudation.py --seeds 1 2 3 --steps 60
```

### `measure_orogenic_relief.py`

Where collision crust goes, and whether continents keep their area and volume (#290): land
fraction, continental-host area, live continental material, cumulative delamination, Hc-cap
saturation, continental Hc mean and p99, and continental Hc volume against the first
snapshot. On backends that have them it also records restite volume, mean Moho lag, the area
in each crust state and every `world.orogenic_relief_budget` account.

```bash
backend/.venv/bin/python bin/debug/measure_orogenic_relief.py --seeds 1 2 3 --steps 100
```

### `measure_erosion_budget.py`

Measures how closely one `erosion.apply_erosion` call conserves material (#275). It runs
erosion alone a few times on each generated world and, around every call, measures total Hc
volume, the continental tracer's live volume, the ledger's declared sinks (so a tracer drop
splits into declared and silent), tracer volume on oceanic hosts, and wall time. This script
has its own options (`--seeds`, `--steps`, `--years`, `--num-plates`, `--out`) and no
`--backend`.

```bash
backend/.venv/bin/python bin/debug/measure_erosion_budget.py --seeds 3 4 --steps 5
```
