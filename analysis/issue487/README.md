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
5.50 TN/m (mature oceanic). These are screening values, not calibrated predictions; integrated
capacity is sensitive to background horizontal stress, pore pressure, mantle water state, and
temperature profile.

The shell diagnostic reports the largest positive principal stress increment against this
additional capacity in N/m, under an explicit normal-faulting `K0` assumption. A curved-sphere
benchmark recovers uniform isotropic tension to 1.18% mean error and 2.75% p95 error on a
17 × 33 mesh. Free in-plane displacement is the isolated-plate baseline; paired clamps shift
the mean by as much as 46.6%, and actual neighboring-plate coupling remains outside the solve.

The first matched replay scored 130 candidates and reported 29 residual/load ratios above 1%.
That result is retained in `gate2-matched-replay.json` under
`pre_edge_connectivity_guard`; it included element meshes that touched only at vertices and
therefore were not connected for the shell solver's rigid-mode treatment. The nine flagged
cases were:

| Seed | World time (Myr) | Plate | Triangle component sizes | Old residual/load |
|---:|---:|---:|---:|---:|
| 7 | 38 | 0 | 3316, 1, 1 | 11.51% |
| 7 | 344 | 16 | 3861, 1 | 26.69% |
| 11 | 982 | 42 | 2586, 1 | 1.15% |
| 2024 | 90 | 1 | 11264, 14, 1 | 3.33% |
| 2024 | 306 | 9 | 5610, 2 | 5.85% |
| 2024 | 330 | 9 | 5560, 2 | 3.12% |
| 2024 | 348 | 9 | 6879, 2 | 1.89% |
| 2024 | 362 | 9 | 6804, 2 | 5.04% |
| 2024 | 452 | 31 | 2852, 3 | 6.43% |

The issue was unsupported mesh topology: separate triangle components introduced extra rigid
modes, but the solve projected only the plate-wide rigid modes. The residual therefore measured
an unsatisfied force balance on an invalid mesh, not evidence that the rock exceeded capacity
or that plate motion was miscalculated. A further edge-disconnected mesh at seed 7, 68 Myr,
plate 14 (components 5940 and 1) was rejected before solving; its residual is therefore not
part of the accepted-solve metrics.

### Connected-source mesh repair and replay coverage audit (2026-10-21)

The 130 Gate 1 craton-crossing candidates are complete production cut events (`split` or
`failed_rift`) whose cut band intersects at least one craton cell. A candidate is the event-level
load diagnostic target, not an independently nucleated local crack. Candidate identity and
production outcome are unchanged by this diagnostic.

The shell mesh now follows this order: require a source plate connected by exact cell-edge
adjacency; build the clipped curved mesh; if the source is connected but triangles are not edge
connected, try a spherical-hull mesh clipped with centroid and edge-midpoint probes; if that does
not reconnect it, retain the dominant triangle component only when it contains at least 99% of all
triangles. Every solve must then cover at least 95% of both candidate cut-band nodes and craton
cut-band nodes. Genuinely source-disconnected plates remain excluded as a whole. Keeping only their
largest component would invent a connected plate/load path, and solving their islands separately
would require an unsupported boundary/load condition.

The updated full 1 Gyr replay uses repeatable Gate 1 settings (seeds 7, 11, 2024; 2 Myr steps;
craton formation and post-PR-490 spherical relief disabled). Production cuts and ordered candidate
keys exactly match Gate 1 for every seed.

| Seed | Production cuts | Candidates | Source connected | Valid shell solves | Source-disconnected exclusions | Repaired meshes | All-candidate coverage |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 7 | 144 | 44 | 23 | 23 | 21 | 3 | 52.3% |
| 11 | 111 | 42 | 26 | 26 | 16 | 1 | 61.9% |
| 2024 | 125 | 44 | 25 | 25 | 19 | 6 | 56.8% |
| **Pooled** | **380** | **130** | **74** | **74** | **56** | **10** | **56.9%** |

All ten previously edge-disconnected meshes are now solved: six used the spherical-hull outline
probe repair (seed 2024 cases) and four dropped only one or two isolated triangles from the dominant
edge-connected component. The exact affected events were:

| Seed | Time (Myr) | Plate | Original triangle components | Repair |
|---:|---:|---:|---|---|
| 7 | 38 | 0 | 3316, 1, 1 | Keep dominant component (99.94% triangles) |
| 7 | 68 | 14 | 5940, 1 | Keep dominant component (99.98%) |
| 7 | 344 | 16 | 3861, 1 | Keep dominant component (99.97%) |
| 11 | 982 | 42 | 2586, 1 | Keep dominant component (99.96%) |
| 2024 | 90 | 1 | 11264, 14, 1 | Spherical hull, centroid + edge-midpoint outline probes |
| 2024 | 306 | 9 | 5610, 2 | Spherical hull, centroid + edge-midpoint outline probes |
| 2024 | 330 | 9 | 5560, 2 | Spherical hull, centroid + edge-midpoint outline probes |
| 2024 | 348 | 9 | 6879, 2 | Spherical hull, centroid + edge-midpoint outline probes |
| 2024 | 362 | 9 | 6804, 2 | Spherical hull, centroid + edge-midpoint outline probes |
| 2024 | 452 | 31 | 2852, 3 | Spherical hull, centroid + edge-midpoint outline probes |

All repaired solves retain at least 95% of the full cut band and craton subset; across all 74 solves
the minimum sample coverage is 99.67% and 99.44%, respectively. The solver stops with code 1 on
every solve. Residual/load is 0.0111% median, 0.0738% p90, 0.1362% p95, and 0.4514% maximum. The
old nine high-residual cases reflected unsupported mesh topology; the repaired replay does not
reproduce them as high-residual solves.

#### Minimum coverage and interpretation

For the explicitly scoped Gate 2 population, require valid shell solves for at least 90% of
craton-crossing candidates whose source plate is one exact edge-connected cell component overall,
and at least 80% within each seed, event kind (`split` / `failed_rift`), world-time bin (0–200,
200–500, 500–1000 Myr), and every nonempty event-kind × time-bin stratum. A scored event must also
meet the per-mesh 95% cut-band and craton-band sample coverage thresholds, stop code 1 and
residual/load ≤1%. The 90% and 80% thresholds are project completeness criteria, not physical
constants. Empty strata are not applicable. Source-disconnected candidates are outside this Gate 2
coverage denominator, but their number and distribution remain reported as a scope limitation.

The mesh and numerical criteria pass for all 74 source-connected candidates: coverage is **100%
(74/74)** overall and in every nonempty seed, kind, time and kind × time stratum. The 56
source-disconnected events (43.1% of the 130-event full replay) recur on 22 unique `(seed,
plate_id)` pairs; 46 events have no off-main source component larger than three cells, while 10
contain a detached component of at least six cells (one reaches 369 cells). Those are separate
source regions, not tiny mesh islands, so retaining just the dominant source component would change
the material domain and its load path. The exclusions are not evenly distributed:

| Candidate stratum | All candidates | Source connected and solved | Connected share of all |
|---|---:|---:|---:|
| Split | 69 | 33 | 47.8% |
| Failed rift | 61 | 41 | 67.2% |
| 0–200 Myr | 75 | 51 | 68.0% |
| 200–500 Myr | 34 | 14 | 41.2% |
| 500–1000 Myr | 21 | 9 | 42.9% |

The Gate 2 strength comparison is complete only for the explicitly scoped connected-source
population. Because connected-source eligibility is uneven (it covers 47.8% of splits versus 67.2%
of failed rifts, and 68.0%, 41.2%, and 42.9% across the three time bins), capacity fractions are
conditional and must not be generalized to disconnected-source events or the full Gate 1 event
population. Report all 56 disconnected-source events separately by seed/type/time; do not relabel
them as mesh failures. The connected-source coverage requirement passes. Gate 2 remains open for
physical source/capacity justification and broader validation, not for the excluded-population
coverage.

On the 22,896 accepted cut-band samples, median-capacity physical failed fractions are 0.274%
continental, 78.03% young-ridge oceanic proxy and 2.84% mature oceanic. At the phenomenological
1× reference, 2.82% of all samples fail (0.34% in marked belts, 3.50% outside). These are screening
outcomes only: saved state has no oceanic age field, and capacities remain conditional scenario
assumptions without spatial calibration.
Per-seed results, coverage strata, hashes and mesh repair details are in
[`gate2-matched-replay.json`](gate2-matched-replay.json); full outputs remain in `/private/tmp`.

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


## Gate 2 capacity scale calibration (2026-10-21)

A repeatable scale audit now checks the physical column envelope at 70, 90, 100, 135 and
150 km total continental lithosphere thickness, with the current 35 km crust held fixed. Its
provisional `K0=0.8` anchor gives 5.91 TN/m at 100 km, within the 3–8 TN/m regional EARS
integrated-strength estimate; the 100 km full factorial grid has median 4.43 TN/m and p10–p90
1.37–11.46 TN/m. A 150 km thick-craton scenario gives 19.02 TN/m, but the cited thickness model
is not a capacity observation. The matched replay loads are not rupture labels and were not used
to fit parameters. This calibrates a conditional scale bracket, not a production strength law.
See [`gate2-capacity-calibration.md`](gate2-capacity-calibration.md) and
[`gate2-capacity-calibration.json`](gate2-capacity-calibration.json).

### Thickness-aware spatial screen (2026-10-21)

The full replay was rerun on the pinned Gate 1 source, with the diagnostic force-field helpers
overlaid while leaving production torque routines unchanged. Production cut totals and ordered
candidate time/plate identities match Gate 1 for all three seeds. All 74 connected-source
candidates solve successfully; 56 disconnected-source candidates remain excluded. The local
capacity map uses each node's actual Hc/Hm and crust-type proxy, and compares steady conductive and
linearly interpolated Moho-lag profiles. It adds strata for craton share and inherited-belt marker.

Across 22,896 accepted cut-band nodes, 1,746 (7.63%) exceed local capacity in the steady-profile
reference; 1,738 (7.59%) do so with the Moho-lag interpolation. Median event-level median
`N_ext/N_cap` is 0.116; median event-level p95 ratio is 0.669. These are scenario-screening
statistics, not observed rupture rates or a fit to split/failed-rift outcomes. Local capacity
remains conditional on assumed thermal/rheology/stress state, so this step does not close Gate 2.
See [`gate2-capacity-state-audit.md`](gate2-capacity-state-audit.md) and the compact
[`gate2-spatial-capacity-replay.json`](gate2-spatial-capacity-replay.json); full event reports are
in `/private/tmp`.

### Capacity assumption sensitivity (2026-10-22)

The prior replay retained event summaries, not per-node load/state, so the matched Gate 1 events
were captured again with capacity scoring disabled. Cut totals, ordered candidate identities, and
all 74 connected-source solves match the accepted baseline. A 13-case one-factor-at-a-time
post-processing sweep then held the exact 22,896 node loads fixed while varying the documented
temperature, strain-rate, pore-pressure, `K0`, cohesion, friction and rheology endpoints; the
stored Moho-lag profile was compared separately at the reference scenario.

The reference exceedance fraction is 7.63%. The largest individual shift is `K0`: 4.44% of nodes
exceed capacity at `K0=1.0`, versus 12.47% at `K0=0.6`. Wet olivine gives 10.69%, and zero pore-
pressure ratio gives 4.90%. Relative node ranking stays highly correlated with the reference
(Spearman at least 0.9938), while 1,009 nodes (4.41%) exceed capacity across all tested one-factor
scenarios and 3,016 (13.17%) do so in at least one. This is not a joint uncertainty interval.

This identifies the unobserved stress and fluid assumptions as the main capacity sensitivities;
it does not select them or close Gate 2. Full assumptions, per-seed/stratum results and hashes are
in [`gate2-capacity-sensitivity.md`](gate2-capacity-sensitivity.md) and
[`gate2-capacity-sensitivity.json`](gate2-capacity-sensitivity.json). Exact captured per-node loads
remain outside the repository in `/private/tmp/`.
