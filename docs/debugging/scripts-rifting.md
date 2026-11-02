# Rifting and thin-shell diagnostic scripts

These offline diagnostics support RFC 0010's Gate 2 assessment of local loading and lithosphere
capacity. Gate 2 is closed for the connected-source, single-plate diagnostic. Its capacity results
are conditional scenarios; they do not select a production rupture threshold or validate coupled
neighbor stresses. See the [RFC](../../rfc/0010-rifts-follow-weak-lithosphere.md) and
[issue 487 findings](../../analysis/issue487/README.md).

Run scripts from the repository root with the worktree's `backend/.venv/bin/python`. Several tools
replay a full world and use substantial memory, so follow the
[`compute-budget` skill](../../.agents/skills/compute-budget/SKILL.md) before starting one. Keep full saves and captures in scratch space; the compact outputs and findings belong
in `analysis/issue487/`. Each script's module docstring gives its complete usage and arguments.

### `audit_gate2_assumption_basis.py`

Calculates independent input-scale quantities for heat flow, strain rate, pore pressure and the
reference frictional K0 limit. It writes `gate2-assumption-basis.json`; it does not read event loads
or fit assumptions to old rift outcomes.

```bash
python3 bin/debug/audit_gate2_assumption_basis.py \
    --output analysis/issue487/gate2-assumption-basis.json
```

### `benchmark_gate2_coupled_shell.py`

Checks a manufactured two-region curved shell with a shared interface. The audit compares global
force and torque balance, interface reactions and recovered extensional resultants. It validates
the stated numerical benchmark, not a geological contact or finite-slip law.

### `calibrate_gate2_capacity.py`

Compares named column-capacity scenarios with published lithosphere strength scales and the matched
replay loads. The literature comparisons are order-of-magnitude checks with different spatial and
physical definitions; this is not a fit to legacy split locations or a calibrated global strength
field.

### `diagnose_gate2_column_envelope.py`

Computes brittle Mohr–Coulomb normal-fault strength and temperature/rate-dependent ductile flow
stress through representative lithosphere columns, then integrates the weaker additional strength
with depth. It reports a scenario envelope for missing local temperature, fluid, composition,
strain-rate, friction and background-stress state. It is a diagnostic rheology model, not production
rheology.

### `replay_gate2_assumption_screen.py`

Post-processes captured Gate 1 node loads without another world replay. It crosses the reference
stress, dry/wet mantle rheology and two pore-pressure profiles; its output records scenarios,
classification changes and input hashes. The near-critical K0 case is derived from the reference
friction angle and is a limiting sensitivity, not a measured stress estimate or probability bound.

### `replay_gate2_capacity.py`

Replays seeds 7, 11 and 2024 against the pinned Gate 1 source with diagnostic helpers overlaid. It
checks that candidate identities match Gate 1, builds a shell load for each craton-crossing event,
and compares physical column scenarios with a labeled phenomenological capacity alternative. The
capacity scores do not feed back into the replay. Follow the compute-budget skill before each seed, and
check the output's baseline count-match fields before interpreting results.

### `replay_gate2_coupled_event.py`

Replays one event for a source plate and its dominant force-contributing neighbor, then compares
pair-level stress under current contact loads and a bonded solve without the pair's mutual contact
drag. It also reports force/torque closure and rigid-mode projection sensitivity. These cases are
mechanical sensitivity probes, not a validated interface law or physical coupled-stress prediction.

### `sensitivity_gate2_capacity.py`

Re-scores captured node loads over one-factor-at-a-time capacity alternatives. It varies the
assumptions recorded in the script and compares each result against the same load samples. The
result is a sensitivity study, not a joint uncertainty distribution or fitted rupture label.

### `validate_thin_shell_plate_meshes.py`

Builds representative meshes from generated plates and checks clipping/coverage, manufactured
loads on the original curved mesh and a shallow-cap mapping, localized and paired clamp sensitivity,
and an independent spherical-shell solution. Its capacity response sweep checks a phenomenological
rule's behavior only; the chosen values are not physical strength estimates.

### `probe_rift_settings.py`

Looks for active, passive, rollback and transtension rift settings, routes a craton-free rift from
each site, and compares it with the great circle through that site. Read-only. `--saves` surveys
snapshots. `--seed` is RFC 0010 Gate D1: one fresh world with production topology unchanged, detecting
just before topology, merging boundary fragments into fronts and tracking them until they persist
for 2 Myr. `--no-detect` is its control arm, whose `fingerprint` must match. `--summarize` builds the
gate tables from the per-seed reports.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_rift_settings.py \
    --saves <a.mbworld> <b.mbworld> --out /tmp/rift-settings.json --image-dir /tmp/rift-images
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_rift_settings.py \
    --seed <n> --myr 200 --step-myr 1 --out <seed.json> --save-dir <private-dir>
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_rift_settings.py \
    --summarize <a.json> <b.json> <c.json> --out analysis/rfc0010/gate-d1.json
```

### `probe_age_polarity.py`

Runs the `probe_rift_settings.py` census on one seed and records why ocean-ocean age verdicts stay
undecided and why frozen verdicts reverse. Read-only; its fingerprint matches the census.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_age_polarity.py <seed> <out.json>
```
