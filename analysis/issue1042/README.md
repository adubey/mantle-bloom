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
