# RFC 0010 diagnosis gate 1: event-stream reconciliation

## Result

The gate-one probe now records a paired shadow-calculator result for every production
rift-cut event, including the numeric size, pole-separation, residual and connectivity
margins. A clean three-seed replay completed on branch commit
`0118752d0ffe75060464204d5b94f81a60fe9a07` with cratons-first layout, craton formation
off, 2 Myr steps through 1 Gyr, a 1,000 km corridor, cost floor 0.05, and
`PYTHONHASHSEED=0`. Python and library versions, clean worktree status, and SHA256 hashes
of the probe, shadow calculator and production source are recorded in
[`gate1-summary.json`](gate1-summary.json).

| Seed | Tap rows | Production `all_cuts` | Flow-through candidates | Admissible | Line-clear admissible |
|---:|---:|---:|---:|---:|---:|
| 7 | 144 | 144 | 44 | 17 | 6 |
| 11 | 111 | 111 | 42 | 31 | 12 |
| 2024 | 125 | 125 | 44 | 28 | 5 |
| **Pooled** | **380** | **380** | **130** | **76** | **23** |

Each seed's sequential tap rows exactly match its final production ledger cut count. The
pooled through-craton rates are `L=0.177` and `B=0.837`. The gate margins are summarized
over those 130 flow-through candidates in `gate1-summary.json`; individual failure counts
can overlap when an event fails more than one gate.

Two independent 100 Myr runs with the same fixed hash seed produced identical event rows.
A clean 100 Myr run at archived model commit `68f8bf3` also produced the same rows as
`origin/main` under the same settings. More importantly, a second full 1 Gyr run on the
same source files and runtime reproduced every cut row from the first full run, for all
three seeds. Its model commit differs because the issue branch was rebased to include the
RFC tracker link, but its recorded source hashes are identical. The second run's compact
provenance and results are in the `repeat_rerun` section of `gate1-summary.json`; raw output
is `/tmp/issue487-gate1-rerun-2026-10-19.json`.

## Archived stream discrepancy

The archived production-probe and margin-instrumented series report the same git commit and
nominal seed/layout/step/duration settings, but their event signatures differ substantially
from each other and from the clean run. Signature identity is `(kind, elapsed Myr, plate ID,
node count)`. Only 3–6 events per seed match either archived series in the new full run.
The production series has 125/137/130 events for seeds 7/11/2024; the margin series has
130/116/115; the clean run has 144/111/125. Keep these trajectories and denominators
separate.

The archived JSONs omit working-tree fingerprints and runtime/library versions. Their
shared commit hash is not enough to determine whether uncommitted model or diagnostic code,
or a runtime difference, caused the divergence. The second clean run reproduces the first
clean run bit-for-bit, but the archived event streams diverge from it within the first 3–6
events. Counts and matching signature occurrences are tabulated in `gate1-summary.json`.

## Gate 1 disposition (2026-10-19): closed

Use the paired, provenance-recorded current runs as the repeatable baseline for subsequent
comparison. Explicitly retire the older production and margin archives from calibration because
their missing worktree/runtime provenance cannot be recovered and their event streams do not
match the repeatable runs. Keep the historical archive comparison as a limitation; do not pool
its counts or margins with the baseline. The paired baseline has 380 production cut rows
across seeds 7/11/2024 (144/111/125), with each seed's tap count matching its production
ledger count and all rows reproduced bit-for-bit in the second full run.

Signature identity `(kind, elapsed Myr, plate ID, node count)` is not enough to align all
events: among its 4, 6, and 4 matches to the archived margin streams, exact numeric/status
margin values match on only 3, 1, and 1 events respectively. Most discrepancies are small,
but one seed-2024 match has a 1,000 vs. 1,041 minimum-half-node count. The paired probe
avoids this ambiguity for the current baseline by evaluating shadow margins on the exact
production event before recording it. The archived mismatch still cannot be attributed to
a specific cause from saved metadata; keep diagnosis gate 1 open until the old runs are
recovered with provenance or explicitly excluded from calibration.

## Reproduction

Run from the repository root after `./bin/setup.sh`:

```bash
PYTHONHASHSEED=0 PYTHONPATH=backend backend/.venv/bin/python \
  bin/debug/probe_craton_mosaic.py \
  --layout cratons-first --seeds 7,11,2024 --years 1e9 --step 2e6 \
  --every 1e8 --out /tmp/issue487-gate1-current-main-instrumented.json
```

The full 446 KB output is intentionally left in `/tmp`; the compact statistics,
fingerprints, and archived-series comparison are checked in as `gate1-summary.json`.

## Gate 2 column-strength diagnostic

The standalone [`diagnose_gate2_column_envelope.py`](../../bin/debug/diagnose_gate2_column_envelope.py)
implements an illustrative depth-integrated normal-faulting strength envelope. Its input
assumptions, material-law sources, reference profiles, full scenario grid, and numerical
resolution comparison are documented in [`gate2-load-diagnosis.md`](gate2-load-diagnosis.md)
and [`gate2-column-envelope.json`](gate2-column-envelope.json). At 0.125 km depth spacing,
reference capacities are 23.79 TN/m (continental), 0.489 TN/m (young-ridge oceanic), and
5.50 TN/m (mature oceanic). The scenario-grid p10/median/p90 values are 5.45/12.70/27.55,
0.133/0.292/0.546, and 1.88/3.73/6.86 TN/m, respectively. These are not calibrated
predictions; the integrated capacity is especially sensitive to the assumed background
horizontal stress, pore pressure, mantle water state, and thermal profile. The shell diagnostic
now reports the largest positive principal stress increment against this additional capacity,
both in N/m, under an explicit normal-faulting `K0` assumption. An independent curved-sphere
benchmark recovers uniform isotropic tension to 1.18% mean error and 2.75% p95 error on a 17 × 33
mesh. Free in-plane displacement is selected as the isolated-plate diagnostic baseline; paired
clamps shift the mean by as much as 46.6%, and actual neighboring-plate coupling is outside the
solver. Gate 2 remains open for matched replay comparison of the physical and phenomenological
capacity choices. See the expanded evidence in
[`representative-mesh-summary.json`](representative-mesh-summary.json).

Reproduce with:

```bash
backend/.venv/bin/python bin/debug/diagnose_gate2_column_envelope.py \
  --dz-km 0.125 --output analysis/issue487/gate2-column-envelope.json
```

Reproduce the curved-shell and original-mesh boundary checks with:

```bash
PYTHONHASHSEED=0 PYTHONPATH=backend backend/.venv/bin/python \
  bin/debug/validate_thin_shell_plate_meshes.py \
  --seed 7 --num-plates 12 --node-density 0.5 \
  --out analysis/issue487/representative-mesh-summary.json
```
