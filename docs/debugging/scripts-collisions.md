# Collision and plate-motion scripts

These scripts measure what happens at plate contacts and in the plate-motion solve: craton
jams, collision overlap and the force that would stop it, collision polarity, the force-proxy
rules, and how rotation rate moves land. The model is in
[simulation-model/tectonics.md](../simulation-model/tectonics.md).

Each probe replays a full world in one process and only taps the step, so nothing it measures
feeds back into the run. Run one replay at a time: each holds a full world in memory
(0.6–2.7 GB). See [Running the scripts](../debugging.md#running-the-scripts).

## Contacts and collisions

### `probe_craton_contacts.py`

The measurement gate of [RFC 0008](../../rfc/0008-craton-collision-jamming.md) (#395): craton
jams at plate contacts. It replays one seed from generation and records, every step and for
every pair of plates in contact:

- **contact**: length, normal closing rate, and craton area within `--reach-km` (500 km);
- **jam**: contact nodes that are retreat candidates holding craton, or cratonic on both
  sides (`collision_polarity.retreat_candidates`, shared with the model's jam test). A pair is
  jammed when at least `MERGE_MIN_CONTACT_NODES` nodes are and it isn't diverging;
- **merge state**: `collision_progress` before and after the update, rolls, drops, and each
  fusion's kind, craton loss and craton-on-craton overlap;
- **craton loss near contacts**, cell by cell, booked to the nearest contact in reach;
- **buoyant pull**: the share of slab pull from continental or cratonic subducting nodes;
- **plates**: count, mean speed and the largest plate's share;
- **contact force** per metre of trench, as applied and as the gate estimates it.

The summary reports J (craton destruction near jammed contacts that didn't fuse within
50 Myr, as a share of all destruction) and the other gate measures.
`probe_craton_clips.py --contacts` and `probe_collision_budget.py` reuse its `Probe`.

```bash
backend/.venv/bin/python bin/debug/probe_craton_contacts.py --seed 7 --years 1e9 \
    --out /tmp/rfc0008-seed7.jsonl
backend/.venv/bin/python bin/debug/probe_craton_contacts.py \
    --summarize /tmp/rfc0008-seed*.jsonl --json analysis/rfc0008/gate.json
```

### `probe_collision_budget.py`

The collision-resolution gate (#452, exploration #449): where overlap builds at continental
contacts and what force would stop it. For every contact pair in scope (continental on both
sides, or craton within `--reach-km`), each step records:

- **G1, overlap**: overlapping area, penetration depth, and the double-counted volume D (both
  columns' Hc over the overlap ground). Each overlapping node is classed as craton-blocked,
  cap-blocked (within `CAP_TOL_M` of `SUTURE_ACCRETION_MAX_HC_M`) or deformable, and D's
  change is booked to that class;
- **G2/G3, forces**: the pair's motion systems from `torque.integrate_omega`, solved jointly.
  The *stall force* is the force per metre of trench that holds the edge's mean closing at
  zero. The *push* is the gravitational potential energy contrast of the contact band's
  columns against the reference continental column. The output also gives the slab-pull share
  of the stall force and H\*, the band thickness whose push equals it;
- **G4, budget**: crust driven in against D's change, retreat removal, band thickening and
  the deform pass's sinks;
- **G5**: at each merge, the absorbed plate's area and Hc over third plates, the overlap's
  age, and Hc and craton lost.

`--capture-dir` saves each merge with RFC 0001's signature for `measure_merge_placement.py`.
`--untapped` runs the same loop with no taps and writes only the end craton volume, to check
that the probe doesn't change the run. `--summarize --pool` pools pair-steps across runs.

```bash
backend/.venv/bin/python bin/debug/probe_collision_budget.py --seed 326494039 \
    --out /tmp/i452/seed326494039.jsonl --capture-dir /tmp/i452/captures
backend/.venv/bin/python bin/debug/probe_collision_budget.py --seed 326494039 --untapped \
    --out /tmp/i452/untapped-326494039.jsonl
backend/.venv/bin/python bin/debug/probe_collision_budget.py \
    --summarize /tmp/i452/seed*.jsonl --json analysis/issue452/gate.json
```

## Collision polarity

These three scripts examine `collision_polarity.py`, which decides which plate goes under at a
continental collision front (#318).

### `report_collision_polarity.py`

Loads a save (or generates `--seed`), steps it, and records at every checkpoint: which
evidence decided each new front (`consumption`, `slab`, `arc`), how often the fallback
decided and why, ambiguous and overruled evidence, record churn (fronts created, split,
merged, expired, dropped, fused), and cost (prepass seconds against the step, and boundary
searches run and saved by the shared cache). A loaded save has no contact history, so fronts
already in contact decide by fallback. `--warmup-myr` steps past that before the counters are
reset.

```bash
backend/.venv/bin/python bin/debug/report_collision_polarity.py \
    world.mbworld --myr 43.1 --out analysis/issue318/polarity-seed997271774.json
```

### `diagnose_polarity_fallbacks.py`

Replays a save and, for every front decided after the first step, records the pair's contact
history over the preceding steps: band nodes facing the other plate, how many converge, and
how many of them (or the node they face) are oceanic. It also records the nearest evidence on
each plate (distance in spacings, source, role, age).

```bash
backend/.venv/bin/python bin/debug/diagnose_polarity_fallbacks.py \
    world.mbworld --myr 43.1 --out analysis/issue318/fallback-diagnosis.json
```

### `diagnose_polarity_merge_conflicts.py`

Replays a save and, whenever front matching absorbs a record whose lower plate differs from
the surviving record's, records both: how each was decided, age, contact steps, split
lineage, stored size and the share of the merged front each covered.

```bash
backend/.venv/bin/python bin/debug/diagnose_polarity_merge_conflicts.py \
    world.mbworld --myr 43.1 --out analysis/issue318/merge-conflicts.json
```

### `polarity_replay_lib.py`

Library, not a script. The two `diagnose_polarity_*` scripts share its command line (`save`,
`--myr` default 43.1, `--step-years`, `--out`) and replay loop, which wraps a
`collision_polarity` function to record what it decides.

## Plate motion and forces

### `audit_force_rules.py`

[RFC 0012](../../rfc/0012-physics-first-simulation.md) gate G1 (#465): which force-proxy rules
bind. It replays one seed from generation and taps, without feeding back:

- **rate clamp**: every `mantle.clamp_rate` call by caller, and the share clipped at the
  0.5 cm/yr floor and the 15 cm/yr cap;
- **motion solve**: per `integrate_omega` call, the slowest relaxation time, the step over
  that time, and the distance from the quasi-static solution;
- **merges**: each fusion by path (timer and roll, or the deep-overlap backstop), with the
  contact's closing rate;
- **splits**: each `maybe_split_plate` call by the gate that decided it;
- **convergent yield**: nodes past yield or near the rate law's saturation, and
  `shortening.ring_cascade` demand absorbed, decayed and unrouted;
- **order dependence**: at `--order-myr` checkpoints, each plate's omega after shifting
  plates in `--order-shuffles` shuffled orders;
- **band rows**: boundary band nodes per distinct neighbour node.

`--profile` instead times one density's step phases, for the gate's cost row.

```bash
backend/.venv/bin/python bin/debug/audit_force_rules.py --seed 326494039 --myr 500 \
    --out analysis/rfc0012/rules_326494039.json
backend/.venv/bin/python bin/debug/audit_force_rules.py --seed 326494039 --profile \
    --density 4 --out analysis/rfc0012/profile_d4.json
```

### `probe_boundary_load_normalization.py`

Checks that slab pull, ridge push, slab drag and contact drag give the same force per metre
of boundary whatever the node density, band width (`BOUNDARY_FORCE_REACH_MULTIPLIER`) and
row sampling (RFC 0012 G1, #465 and #468). It builds two oceanic plates split along the
0/180 meridian: the west plate rotates east at `--rate-cm-yr`, making a trench at 0° and a
ridge at 180°. For the west plate's own band, as `shift_plate` reads it, each load per
represented front length is compared with its analytic value. Each density is built at
`--offsets` lattice alignments. A load passes if its spread across a sweep is at most 10%.
No save or seed needed.

```bash
backend/.venv/bin/python bin/debug/probe_boundary_load_normalization.py \
    --out analysis/rfc0012/load_normalization.json
```

### `probe_rotation_pileup.py`

Where the quad "pile-up" comes from (#289, found on #177). With `avg_rotation_rate` at 4×,
worlds keep less land than at 0.25× but much more land volume above sea level. This splits
the land-volume metric per (multiplier, seed, checkpoint) into the nominal-area sum, the
exact-area sum (`Plate.accounting_areas_m2`), and the exact sum with overlap counted once
(`stats.overlap_area_weights`). It also records represented and overlapped sphere fractions,
land area, mean land elevation and Hc, and the continental-material ledger.

```bash
backend/.venv/bin/python bin/debug/probe_rotation_pileup.py --seeds 350921662,57685824 \
    --out bin/debug/results/rotation_pileup_probe.jsonl
```

### `rotation_ledger.py`

Per-step land-area against crust-volume ledger for #177, the follow-up to #171's sweep, which
found rotation rate the strongest lever on land decline. It tracks translation-invariant
aggregates per multiplier (land, continental and total node counts and areas, continental
crust volume and mean continental Hc), so area lost to node attrition at retreating margins
shows as area falling faster than volume. Defaults to #171's seed 829071382 and writes to
`bin/debug/results/rotation_ledger.jsonl`.

```bash
backend/.venv/bin/python bin/debug/rotation_ledger.py --seed 829071382 --multipliers 0.25,4.0
```
