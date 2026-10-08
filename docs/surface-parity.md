# Plate-surface audit harness

Issue #247 (Phase 5a of #228). The harness steps sparse-quad worlds, records their
measurements, and judges them against the gates below. It began as a line-vs-quad parity
harness; since the line surface was retired (#251), every run audits the quad surface on its
own and every gate has an absolute threshold.

Code: `backend/app/surface_parity.py` (runs, audits, metrics),
`backend/app/surface_parity_gates.py` (gates, report), `bin/debug/surface_parity.py` (CLI).

## Running it

```sh
cd backend
.venv/bin/python ../bin/debug/surface_parity.py audit --preset smoke --seeds 3 --out ../analysis/parity-smoke
```

`audit` runs every seed, then writes `comparison.json` and `report.md`. It exits with status 1
when the verdict is `fail`. Other subcommands:

- `run --seed S --out DIR` runs one world. Use it to spread a campaign
  across processes or machines. `--from-world SAVE.mbworld` continues a saved world instead of
  generating one; checkpoint ages then count from the save's age, and the save's own
  densities apply.
- `judge DIR` re-judges every run in a directory, for example after collecting runs made
  separately.

`--jobs N` runs worlds in parallel. Timings are only comparable between runs made at
`--jobs 1` on an otherwise quiet machine, so a performance measurement needs a sequential run.

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
| `seed<S>.json` | Deterministic metrics: config, checkpoints, audit violations, per-step climate series, load checks. Floats are rounded to 6 significant figures and keys sorted, so a rerun of the same code reproduces the file byte for byte. |
| `seed<S>.timings.json` | Wall-clock only: generation, per-step phase buckets, derived-index builds/hits/build time, revision bumps, and index build/query time and memory at checkpoints. |
| `provenance.json` | Command line, git commit (and whether `backend/app` or `bin/debug` were dirty), Python/numpy/scipy versions, platform. |
| `comparison.json` | Every gate result, and one summary per gate that keeps its worst case. |
| `report.md` | The verdict, gate tables, per-seed checkpoint metrics, and performance tables. |
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
  world points, cell centres, areas, adjacency, loops, row/column intervals and probes) equals a fresh rebuild from authoritative state. So do the world caches
  (`node_position_tree_cache`, `node_kdtree_cache`, `node_healpix_index_cache`);
- **revision tracking**: if a plate's local node set changed since the previous audit, its
  topology revision must have increased; if its frame changed, its geometry revision must
  have increased. Otherwise a cache keyed on revisions could be reused while stale;
- the atmospheric HEALPix grid keeps its resolution. It is Eulerian climate state and does
  not depend on node count.

**Checkpoints** record:

- area-weighted totals from each cell's exact area (`SurfaceNodes.area_m2`): Hc/Hm volume,
  continental Hc, land fraction, elevation p05/p50/p95, plate counts and sea level;
- whole-sphere coverage from a Fibonacci sample tested against each plate's own
  `contains_batch`: uncovered, multiply covered, and void (uncovered and more than 1.5
  spacings from any node). Also the fraction of nodes inside another plate;
- neutral lattice quality (docs/plate-surface-baseline.md §1.1): stacked, anisotropic, row
  alignment, boundary, thin, and fragments;
- quad lattice quality: refined cells, one-cell-thin cells, hole loops, aspect ratio above
  4, skew above 45°;
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
outline and its tree, quad adjacency and loops, and the world KD-tree/HEALPix helpers).
Rebuild frequency is reported as builds per step plus topology and geometry revision bumps per
step.

## Gates

A gate result is `pass`, `warn`, `fail` or `insufficient` (too little data to judge). The
verdict is `fail` if any gate fails, otherwise `warn` if any gate warns, otherwise `pass`.

The tolerances come from docs/plate-surface-baseline.md §3.

| id | gate | applies to | fail | warn |
|---|---|---|---|---|
| H1 | quad leaf topology valid | every audit | any violation | |
| H2 | quad neighbours valid | every audit | any violation | |
| H3 | no folded cells | every audit | any | |
| H4 | fields finite and sized; areas positive | every audit | any | |
| H5 | frames are proper rotations (1e-6) | every audit | any | |
| H6 | derived caches equal fresh rebuilds; exact KD-tree parity; atmosphere grid fixed | every audit | any | |
| H7 | revisions track node-set/frame changes | every audit | any | |
| H8 | save/load identical, derived indexes dropped | checkpoints | any | |
| H9 | loaded-world continuation identical | final checkpoint | any | |
| H10 | no stacked nodes (< 0.5 spacing) | checkpoints | any | |
| H11 | elevation, Hc, Hm within their caps | every audit | | any |
| S1 | climate stats and hydrology finite | checkpoints | any | |
| C2 | void | checkpoints | > 0.05% | |
| C3 | multiply covered sphere | checkpoints | > 2.0% | > 1.5% |
| K3 | Σ cell area = sphere − uncovered + overlap | checkpoints | | off by > 1 pp |
| M1 | anisotropic fraction (ratio > 2) | checkpoints | | > 1% |
| M2 | anisotropic row alignment \|cos 2α\| | when ≥ 0.1% anisotropic | | > 0.3 |
| M3 | thin (tendril) fraction | checkpoints | | > 1% |
| M4 / M5 | cells with aspect > 4 / skew > 45° | checkpoints | | > 0.1% |
| M6 | extra components + isolated nodes | checkpoints | | > 0 |

The gate ids keep their numbering from the line-vs-quad campaigns (#249), whose comparison
gates (C1, C4, K1, K2, K4, S2, X1, R1, R2, P1) were retired with the line surface. The report
still tabulates the per-step phase timings and the HEALPix-vs-KD-tree resampling metrics;
they are no longer gated.

H11 warns instead of failing. Before issue #256, the engine clamped the Hc/Hm caps only in
specific code paths, and worlds breached the Hm floor or ceiling within 15 Myr at density 1
(quad up to 242 km against the 240 km cap).
Every step now ends with a clamp on every plate (`lithosphere.clamp_column_caps`), so an H11
finding means the clamp itself has regressed. The Phase 0a table lists the caps as exact.

Every quad cell is authoritative territory, so any double claim C3 measures is real
duplicated crust. Quad plates also carry independently rotated lattices. At a
seam their cell edges generally cannot coincide, so claiming the last whole cell trades a
thin overlap for the uncovered sliver that leaving it empty would preserve. The #249 long and
stress campaigns found this overlap growing toward a stable 1.2--1.9% plateau by 400 Myr,
alongside much lower uncovered area, nodes inside another plate, and zero same-plate stacking.
The 1.5% warning keeps movement toward the envelope visible; 2.0% is the accepted ceiling.

Gap filling limits that tradeoff by sampling each candidate quad cell's interior and claiming
it only when more than half of the sampled footprint is still uncovered. A free cell centre
alone is insufficient: on independently rotated lattices the neighbour may already own most
of the cell around it. The density-0.5 issue #268 stress rerun reduced the 100--800 Myr C3
range from 1.5--2.2% to 1.1--1.8% without changing the fixed C3 bands.

Duplicated crust is not hidden from conservation totals. Hc and Hm totals sum every cell's
exact area, including cells in multiply covered territory, and K3 independently checks that
the summed quad cell area equals `sphere - uncovered + overlap` within its 1 pp sampling
tolerance. Thus C3 bounds the permitted duplication while K3 verifies that its area remains
explicitly accounted for.

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
follow the geometry revision (rigid rotation) and the topology revision (any cell edit).
World caches are dropped at the start of every `step_world` and on load
(`persistence._drop_derived_caches`). `unit_tests/test_surface_parity.py` walks every quad
mutation (insert, remove, refine, coarsen, merge, split, failed rift, rotation, field-only
writes, pickling). For each one it checks that the right revision moved and that no populated
cache is left stale. The run-time audit repeats that check against real trajectories.
