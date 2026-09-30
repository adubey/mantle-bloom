# Line vs quad parity harness

Issue #247 (Phase 5a of #228). The harness runs line-backed and sparse-quad worlds from the same
seed and settings, records the same measurements for both, and judges them against the gates
below. It does not change production defaults. Diagnostic and test runs pick the surface with
`generate_world(surface="lines" | "quad")`; the harness exposes this as `--surfaces`.

Code: `backend/app/surface_parity.py` (runs, audits, metrics),
`backend/app/surface_parity_gates.py` (gates, report), `bin/debug/surface_parity.py` (CLI).

## Running it

```sh
cd backend
.venv/bin/python ../bin/debug/surface_parity.py paired --preset smoke --seeds 3 --out ../analysis/parity-smoke
```

`paired` runs every seed × surface combination, then writes `comparison.json` and `report.md`.
It exits with status 1 when the verdict is `fail`. Other subcommands:

- `run --seed S --surface quad|lines --out DIR` runs one world. Use it to spread a campaign
  across processes or machines. `--from-world SAVE.mbworld` continues a saved world instead of
  generating one; checkpoint ages then count from the save's age, and the save's own
  densities apply.
- `compare DIR` re-judges every run in a directory, for example after collecting runs made
  separately.

`--jobs N` runs worlds in parallel. Timings are only comparable between runs made at
`--jobs 1` on an otherwise quiet machine, so a performance comparison needs a sequential run.

| preset | density | step | checkpoints (Myr) | audit every | approx. cost per world |
|---|---:|---:|---|---:|---|
| `smoke` | 0.5 | 1 Myr | 0, 2, 4 | 1 step | ~10 s |
| `standard` | 1 | 1 Myr | 0, 30, 60, 120 | 5 steps | a few minutes |
| `long` | 1 | 1 Myr | 0, 30, 60, 120, 240, 400 | 10 steps | tens of minutes |
| `issue147` | 4 (climate 4, fluid 2) | 100 kyr | 0, 10, 30, 60 | 50 steps | hours; the #147 profile world |

Every preset field can be overridden: `--checkpoints-myr`, `--node-density`, `--step-years`,
`--audit-every`, `--samples`, `--climate-density`, `--fluid-density`. `--render` writes
`elevation` and `platesDetail` PNGs at each checkpoint. `--no-load-checks` skips the save/load
checks.

The smoke configuration is also the CI check:
`backend/stress_tests/test_surface_parity_smoke.py` (`./bin/stress_test.sh`, about 20 s). It
asserts that every hard invariant holds and that the metrics document reproduces byte for byte.

## Artifacts

| file | contents |
|---|---|
| `seed<S>-<surface>.json` | Deterministic metrics: config, checkpoints, audit violations, per-step climate series, load checks. Floats are rounded to 6 significant figures and keys sorted, so a rerun of the same code reproduces the file byte for byte. |
| `seed<S>-<surface>.timings.json` | Wall-clock only: generation, per-step phase buckets, derived-index builds/hits/build time, revision bumps, and index build/query time and memory at checkpoints. |
| `provenance.json` | Command line, git commit (and whether `backend/app` or `bin/debug` were dirty), Python/numpy/scipy versions, platform. |
| `comparison.json` | Every gate result, and one summary per gate that keeps its worst case. |
| `report.md` | The verdict, gate tables, side-by-side checkpoint metrics, and performance tables. |
| `renders/` | With `--render` only. |

### What is measured

**Audits** run every `audit_every` steps and at each checkpoint, so invariants are checked
throughout a run, not only at the end:

- fields are finite, have the plate's node count, and stay within their caps (elevation,
  Hc, Hm); node areas are positive;
- frames are proper rotations;
- quad leaf topology: no duplicate or overlapping leaves, and 2:1 balance;
- quad neighbours: in range, no self-loops, symmetric, and every pair shares an edge. The
  edge check uses exact integer cube-surface corners, so it also holds across 2:1 level jumps
  and cube-face seams;
- no folded cells (a corner that turns clockwise in the cell's tangent plane);
- **derived caches**: every *populated* plate cache (outline, outline KD-tree, node KD-tree,
  world points, and on quad plates centres, areas, adjacency, loops, row/column intervals and
  probes) equals a fresh rebuild from authoritative state. So do the world caches
  (`node_position_tree_cache`, `node_kdtree_cache`, `node_healpix_index_cache`);
- **revision tracking**: if a plate's local node set changed since the previous audit, its
  topology revision must have increased; if its frame changed, its geometry revision must
  have increased. Otherwise a cache keyed on revisions could be reused while stale;
- the atmospheric HEALPix grid keeps its resolution. It is Eulerian climate state and does
  not depend on node count.

**Checkpoints** record:

- area-weighted totals from each node's actual area (`SurfaceNodes.area_m2`): exact cell areas
  on quad plates, the line adapter's Voronoi-style estimate on line plates. These cover Hc/Hm
  volume, continental Hc, land fraction, elevation p05/p50/p95, plate counts and sea level;
- whole-sphere coverage from a Fibonacci sample tested against each plate's own
  `contains_batch`: uncovered, multiply covered, and void (uncovered and more than 1.5
  spacings from any node). Also the fraction of nodes inside another plate. Quad territory
  is exact; the line outline approximates its node cloud;
- neutral lattice quality (docs/plate-surface-baseline.md §1.1): stacked, anisotropic, row
  alignment, boundary, thin, and fragments;
- quad lattice quality: refined cells, one-cell-thin cells, hole loops, aspect ratio above
  4, skew above 45°;
- the legacy one-node-line fraction on line plates;
- climate stats, biome fractions, and hydrology summaries (rivers, lakes, glaciers, land sinks);
- **derived-index parity**. First, the cached world KD-tree must answer a query of the sample
  exactly like a freshly built one. Second, HEALPix nearest-node resampling is compared with
  the exact KD-tree: how often it picks the same node, its distance from the query in
  spacings, the elevation error, and how often it picks the wrong plate or the wrong side of
  the coast. These are reported overall and near coasts, plate boundaries, holes, poles
  (|lat| ≥ 80°) and the antimeridian (|lon| ≥ 177°).

**Load checks** run at each checkpoint. A save/load round trip must reproduce the
authoritative state hash and come back with every derived index dropped. At the final
checkpoint, one step from the loaded world must match one step from the in-memory world, which
still has warm caches. That shows no derived index acts as persistent authority.

**Timings** use phase buckets that match `step_world`'s sequence. `deform`, `topology`,
`gap_fill` and `overlap_tracking` add up to the `deform_topology` phase #228 calls out.
Derived-index builds are counted where the cache is filled (`Plate.get_node_kdtree`, the
outline and its tree, quad adjacency and loops, the line row lookup, and the world
KD-tree/HEALPix helpers). Rebuild frequency is reported as builds per step plus topology and
geometry revision bumps per step.

## Gates

A gate result is `pass`, `warn`, `fail`, `insufficient` (too little data, such as an
ensemble gate with fewer than 3 seeds) or `info`. `info` marks a failure on the *line* surface:
it is a baseline finding and never gates the quad surface. The verdict is `fail` if any gate
fails, otherwise `warn` if any gate warns, otherwise `pass`.

The tolerances come from docs/plate-surface-baseline.md §3. `l` is the paired line run at the
same age and seed.

| id | gate | applies to | fail | warn |
|---|---|---|---|---|
| H1 | quad leaf topology valid | quad, every audit | any violation | |
| H2 | quad neighbours valid | quad, every audit | any violation | |
| H3 | no folded cells | quad, every audit | any | |
| H4 | fields finite and sized; areas positive | both, every audit | any (quad) | |
| H5 | frames are proper rotations (1e-6) | both, every audit | any (quad) | |
| H6 | derived caches equal fresh rebuilds; exact KD-tree parity; atmosphere grid fixed | both | any (quad) | |
| H7 | revisions track node-set/frame changes | both, every audit | any (quad) | |
| H8 | save/load identical, derived indexes dropped | both, checkpoints | any (quad) | |
| H9 | loaded-world continuation identical | both, final checkpoint | any (quad) | |
| H10 | no stacked quad nodes (< 0.5 s) | quad, checkpoints | any | |
| H11 | elevation, Hc, Hm within their caps | both, every audit | | any (quad) |
| S1 | climate stats and hydrology finite | both, checkpoints | any (quad) | |
| C1 | uncovered sphere | checkpoints | > l + 0.5 pp | |
| C2 | void | checkpoints | > 0.05% | |
| C3 | multiply covered sphere | checkpoints | > l + 0.5 pp | > l |
| C4 | nodes inside another plate | checkpoints | > l + 0.5 pp | > l |
| K1 | Hc volume drift since 0 Myr, quad vs line | checkpoints > 0 | gap > 10 pp | gap > 5 pp |
| K2 | continental Hc drift, quad vs line | checkpoints > 0 | gap > 10 pp | gap > 5 pp |
| K3 | quad Σ cell area = sphere − uncovered + overlap | checkpoints | | off by > 1 pp |
| K4 | Hm volume drift, quad vs line | checkpoints > 0 | gap > 10 pp | gap > 5 pp |
| M1 | anisotropic fraction (ratio > 2) | quad, checkpoints | | > 1% |
| M2 | anisotropic row alignment \|cos 2α\| | quad, when ≥ 0.1% anisotropic | | > 0.3 |
| M3 | thin (tendril) fraction | quad, checkpoints | | > l or > 1% |
| M4 / M5 | cells with aspect > 4 / skew > 45° | quad, checkpoints | | > 0.1% |
| M6 | extra components + isolated nodes | quad, checkpoints | | > 0 |
| S2 | per-step jitter (std of first differences over the second half) of air/ocean temperature, precipitation, sea level, land fraction | per seed | | > 2 × l + slack |
| X1 | HEALPix vs KD-tree, per category | checkpoints | | same-node rate < l − 5 pp, or p95 distance > 1.25 × l |
| R1 | deform + topology s/step | per seed | | > 1.25 × l |
| R2 | total s/step | per seed | | > 1.25 × l |
| P1 | land fraction, Hc, continental Hc, plates, sea level, elevation p05/p50/p95 | ensemble ≥ 3 seeds, checkpoints > 0 | quad mean outside line mean ± 2σ (seed-to-seed) | |

The comparisons are statistical, never byte-for-byte. Three random streams are keyed on line
indices (docs/plate-surface-baseline.md §4.4), so line and quad trajectories diverge from the
first step. Coverage gates compare each quad run with its own paired line run. P1 compares
the two ensembles.

H11 warns instead of failing. The engine clamps the Hc/Hm caps only in specific code paths, and
at the time of writing both surfaces breach the Hm floor or ceiling within 15 Myr at density 1:
the line surface down to about 290 m, and quad up to 242 km against the 240 km cap. The Phase 0a
table lists the caps as exact, so a breach is an engine finding to fix or re-document, not a
quad regression.

Performance only warns. #228 asks that a dominant-phase regression be explicitly accepted or
fixed, so a `warn` on R1/R2 needs a decision in the campaign report (#249), not a silent pass.
For the `issue147` preset, the report adds a column with the #147 profile's per-step phase
means (`analysis/issue147-profile-20260922/frames.csv`). That run was under cProfile, so
compare it by ratios.

## Derived indexes

Neither index is authoritative plate state:

- exact `cKDTree`s: per plate (`Plate.get_node_kdtree`, outline tree) and whole world
  (`World.node_position_tree_cache`, `node_kdtree_cache`). They serve k-nearest, radius,
  distance, overlap-pair, fault and boundary queries;
- the node-cloud HEALPix index (`World.node_healpix_index_cache`): an optional approximate
  accelerator for nearest-value resampling under `node_cloud_resample_mode="healpix"`;
- the atmospheric HEALPix grid (`World.atmosphere_cfd_state`): independent Eulerian climate
  state at `fluid_density`'s resolution.

Plate caches are invalidated through `_invalidate_bounding_polygon` / `_reset_caches`. They
follow the geometry revision (rigid rotation) and the topology revision (any cell/line edit).
World caches are dropped at the start of every `step_world` and on load
(`persistence._drop_derived_caches`). `unit_tests/test_surface_parity.py` walks every quad
mutation (insert, remove, refine, coarsen, merge, split, failed rift, rotation, field-only
writes, pickling). For each one it checks that the right revision moved and that no populated
cache is left stale. Line mutations (`set_lines`, `replace_line`, rotation) are covered by
`test_plates.py` and `test_plate_surface_contract.py`. The run-time audit
repeats that check against real trajectories.
