# Issue #1042 — shared Rust 3D spatial index (step 5 of #1046)

## Preregistration

Principle #2 (keep the planet moving at a useful pace) calls for compiled work where repeated
point queries limit simulation pace. Before this change, every nearest-neighbour and radius
query went through `scipy.spatial.cKDTree`: about 1.2 s of a 13 s young step and 1.6 s of a
14.6 s mature step, three-quarters of it in queries rather than builds (per-call-site profile
below). The candidate replaces those trees with `SpatialIndex` (`backend/app/spatial_index.py`),
backed by `mantle_bloom_rust.KdIndex` (`backend/rust/src/spatial/`).

The baseline is `origin/main` at `31094ae7`; the candidate is the issue branch. Both use
diagnostic wind, one warmup and three measured 100,000-year steps, run sequentially on the same
host. The young world is seed 0 at densities 4/4/2. The mature world is
`mantle-bloom-seed47235302-67900000y.bloom` (67.9 Myr, 28 plates), SHA-256
`eff1f6d7abdfe46c7af5b45ab23f136f3fac591f4d95f3457017c1e7d933216b`, the save #1046's earlier
steps used.

**Primary measure:** time spent in spatial queries and builds per step, summed over call sites
by `bin/debug/profile_kdtree_sites.py`. **Secondary:** mean step wall time and the cold and
warm `combined` render from `bin/debug/profile_phases.py`, and peak RSS from `/usr/bin/time -l`.
Expect at least a 30% cut in spatial time on both worlds, a step-time reduction of a few
hundred milliseconds, and no material RSS increase. As in step 3, a phase change below 10% or a
total-step change below 2% is inconclusive; replicate noisy results.

**Invariant:** model output is unchanged. `SpatialIndex` reproduces `cKDTree`'s distances bit
for bit (it probes the installed scipy's squared-distance rounding) and its threshold tests;
the two differ only on exact distance ties, where the index takes the lowest index and
`cKDTree` takes whichever tied point its traversal reaches first. Exact ties need coincident
points, so none is expected on these worlds. To check, replay the mature save for 1 Myr in ten
100,000-year steps on each arm with `attribute_hm_ratchet.py --myr 1 --checkpoint-myr 0.1` and
compare every checkpoint's metrics, inventories and ledgers. Any difference falsifies the
parity claim and must be explained before the speed result counts.

## Per-call-site profile on `origin/main`

`profile_kdtree_sites.py --steps 3 --renders 2` after one warmup step, on a quiet host. Each row
is milliseconds per step (or per render) charged to the call site that built or queried the
tree.

| Call site | Young ms/step | Mature ms/step | Shape |
|---|---:|---:|---|
| `torque.BoundarySearchCache.nearest` (per-plate trees) | 246.1 | 341.4 | ~200–350 queries of 4.6–7.5 k points |
| `quad_tectonics.grow_frontier` (build over concatenated neighbours) | 164.3 | 184.1 | 32–56 builds of 28–42 k points |
| `faults` shear/relief/unloading/stress/influence (all sites) | ~280 | ~520 | thousands of 5–20-point queries, per-step stress trees |
| `plates._plates_within` (outline trees) | 51.0 | 81.4 | ~670–1,170 bounded queries |
| `hydrology._nearest_hydro_node_idx` | 49.4 | 47.9 | 170 k points against the node cloud |
| erosion (`_coastal_openness`, `_spread_marine_sediment`, `compute_slope`, `_earthquake_erosion_multiplier`) | ~106 | ~142 | radius counts, k = 5–8 |
| **All spatial queries and builds** | **1,209.5** | **1,555.2** | 326 / 357 ms in builds |
| Render (`combined`, per render) | 20.2 | 20.1 | already Rust apart from the ocean-mask resample |

## Results

All runs were sequential on the same 16 GB, 10-core Apple Silicon host, with other agents'
runs coming and going (1-minute load 5–9). The baseline worktree is `31094ae7`. The first
candidate, used for the profiles and the first replay below, is `fffa09b3`. Its Rust path is
unchanged in the later commits, which change only the scipy fallback's tie order and accept
rank-list `k`.

### Spatial time per call site

`profile_kdtree_sites.py --steps 3 --renders 2`, baseline and candidate back to back.
Milliseconds per step (render: per `combined` render).

| Call sites | Young main | Young branch | Mature main | Mature branch |
|---|---:|---:|---:|---:|
| `torque` nearest boundary node (per-plate trees) | 288.9 | 81.1 | 337.6 | 111.4 |
| `faults` (all sites) | 269.0 | 155.0 | 422.0 | 236.1 |
| `grow_frontier` neighbour search | 160.1 | 11.2 | 179.9 | 15.1 |
| erosion (all sites) | 143.2 | 69.0 | 175.2 | 95.9 |
| `plates` (outline and node trees) | 94.7 | 48.5 | 106.1 | 46.9 |
| hydrology | 92.2 | 41.5 | 86.9 | 42.3 |
| climate | 64.9 | 31.2 | 60.7 | 30.2 |
| gaps, merge, world, other | 115.0 | 61.5 | 135.3 | 71.1 |
| **All spatial queries and builds** | **1,228.1** | **499.0** (−59%) | **1,503.7** | **648.8** (−57%) |
| of which builds | 311.5 | 145.2 | 342.0 | 155.3 |
| of which still bare `cKDTree` | 1,228.1 | 22.9 | 1,503.7 | 30.0 |
| Render, Python-side spatial | 17.9 | 8.3 | 18.4 | 9.4 |

The primary measure fell 59% and 57%, past the preregistered 30%. `grow_frontier`'s cut comes
from querying each neighbour's cached node tree instead of building a tree over their
concatenated nodes for every growth. The rest comes from the index itself: queries 2–6× faster
than `cKDTree` at the same answers, builds about 2× faster, and about 1–2 µs per small query
call against 5–7 µs.

### Step, render and memory

`profile_phases.py --scenarios render,step --steps 3 --warmup 1 --wind-model diagnostic`
under `/usr/bin/time -l`, three interleaved runs per arm. The third mature candidate run was
interrupted, so that row has two.

| Case | Mean step, main → branch | Cold render | Warm render | Peak footprint |
|---|---:|---:|---:|---:|
| Young, seed 0 (3 + 3 runs) | 13.07 → 11.88 s (14.46, 11.75, 13.01 → 12.41, 11.85, 11.37) | 1.41 → 0.87 s | 0.97 → 0.49 s | 1,191 → 1,257 MB |
| Mature, 67.9 Myr (3 + 2 runs) | 13.08 → 12.07 s (14.50, 12.37, 12.38 → 12.39, 11.75) | 1.29 → 0.84 s | 0.86 → 0.47 s | 1,224 → 1,242 MB |

The step reduction, about 1.0–1.2 s (8–9%), is close to the spatial time saved (0.7–0.9 s), but
run-to-run spread on this shared host is ±1 s, so read it as "about the spatial saving", not a
precise figure. Max RSS varied by up to 300 MB between runs of the same arm and is omitted;
peak memory footprint rose 1–6%. The render gains are consistent across runs and come from the renderer's own
k-d tree, which is now the shared one (`backend/rust/src/spatial/kdtree.rs`). Its traversal
rechecks each cell's bound when it is popped and prunes by each node's tight box; it answers the
same nearest-node queries (see the render check below).

### Parity: the preregistered invariant does not hold

The 1 Myr replay (`attribute_hm_ratchet.py <local-save> --myr 1 --checkpoint-myr 0.1`) moved
the model. The first checkpoint already differs (sea level −66.82439978 vs −66.82439977 m), and
the gap grows as small differences propagate through the coupled climate, hydrology and
erosion steps. Compact per-checkpoint values are in `replay-comparison.json`.

| At 68.9 Myr | Main | Branch | Relative |
|---|---:|---:|---:|
| Sea level (m) | −60.8954 | −60.9579 | 1.0e-3 |
| Land area (km²) | 171,862,737 | 171,941,408 | +4.6e-4 |
| Emergent continental area (km²) | 171,386,808 | 171,462,061 | +4.4e-4 |
| Mantle lithosphere volume (km³) | 40,616,446,488.168 | 40,616,446,488.161 | −1.7e-13 |
| Runtime (s) | 146.6 | 139.9 | |
| Peak memory (MB) | 1,102.8 | 1,065.6 | |

`compare_spatial_backends.py` found the cause. It answers every `SpatialIndex` call of one
mature step with a plain `cKDTree` too (`tie-probe-mature-1step.json`). Over about 2.2 million
query rows, no distance and no radius set differed. Every disagreement was an exact distance
tie between *distinct* points: 162 rows in `plates._plates_within`, 54 in
`climate._land_swirl_current`, 16 in `_spread_marine_sediment`, 13 in `faults.fault_influence`,
12 in `hydrology._build_neighbor_graph`, 9 in `compute_slope`, and a few elsewhere. The
preregistration assumed exact ties need coincident points. They don't: lattice and grid
symmetry makes distinct neighbours bit-identically equidistant somewhere in every step. Where
a caller uses only the distance (`_plates_within`, `fault_influence`,
`distance_from_land_approx`), a tie changes nothing. Where it uses the index (the land-swirl
currents, the hydrology neighbour graph, slopes, marine sediment), the tie picks a different but
equally near point. `cKDTree`'s pick depends on its build options, which already varied from one
call site to the next on main.

Reproducing `cKDTree`'s pick would mean reproducing its tree layout and traversal for each
build option. Instead, the lowest-index rule is now the contract on both backends: the scipy
fallback widens its k-NN query until no tie straddles the `k`-th slot and orders rows by
`(distance, index)`, so the Rust and fallback paths agree bit for bit, ties included
(`test_symmetric_ties_take_the_lowest_index`). The model output therefore moves once relative
to main, by a tie resolution rather than any distance or threshold change, in the land-area
range #1047's lake port showed on the same save. The project owner accepted this one-time model
change on review of PR #1058.

### Backend parity after review

The first fallback ordered tied rows by the *returned* distance. Distinct squared distances can
share a square root, so it could pick a farther point that Rust and `cKDTree` both reject (Codex
review of #1058). The fallback now re-orders such rows by the squared distance, rounded as scipy
rounds it, then index. With that fix (`942e56c3`):

- `compare_spatial_backends.py --reference fallback` over one mature step: every call agreed in
  values, order, dtype and shape, at every one of the 33 call sites.
- The 1 Myr replay with the Rust index and with only the spatial index forced onto its scipy
  fallback: bit-identical at all 11 checkpoints (runtime 122 s vs 136 s). The Rust replay is also
  identical to the earlier candidates' (`fffa09b3`, `08a63513`).

### Rendering is unchanged

`combined`, `biome`, `elevation` (plain, and with the mountain and plains toggles) in Eckert IV,
`combined` on the globe and `elevation` in Behrmann, rendered at 1100×611 from the loaded save
in each tree: all six PNGs are byte-identical between `origin/main` and the branch.

## Tectonic kernel comparison preregistration

Principle 2 calls for reducing repeated boundary searches and numerical passes while preserving
plate ownership and material budgets. Compare `origin/main` at `069d1d38` with the tectonic
kernel branch on the mature seed 47235302 save at 67.9 Myr (SHA-256
`eff1f6d7abdfe46c7af5b45ab23f136f3fac591f4d95f3457017c1e7d933216b`). Run ten
100,000-year steps, sequentially, using `attribute_hm_ratchet.py --myr 1 --checkpoint-myr 0.1`.
The primary invariant is matching land area, sea level, crust and mantle inventories, and
material ledgers at each checkpoint; the target is bit-identical results, with any deviation
investigated before accepting the port. The Python-fallback path on the candidate is an
additional control if the baseline and candidate diverge. Runtime and peak RSS are secondary
and a single run cannot establish small speed changes on a shared host. The kernels are
expected to reduce their own time without a material increase in total step time or memory;
changes under 2% of step time are inconclusive without replicates.
