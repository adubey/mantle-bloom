# Profiling

A running record of where wall-clock time actually goes, measured against a real code
path rather than guessed at. Each section says what was measured, on what commit, with
what inputs, so a later run can be compared like-for-like.

## File > Make Animation (`POST /world/animate`), 60 frames

**Measured:** re-run 2026-09-01, commit `f68fa46` (`perf/vectorize-all-points-elevation`
-- `main`, which carries fixes 1-4 below, plus fix 5). 6-core Apple Silicon (2 performance + 4 efficiency
cores), Python 3.14.6, backend `.venv`. Drove the animate path directly (no vite/uvicorn
dev servers up this time), so less background contention than the earlier same-day re-run
at `75e8c95`. Supersedes both that run and the original at `c8f9990`, 2026-08-31; the
`c8f9990` numbers are kept in the tables for the before/after.

**Inputs** -- the frontend's own defaults (`frontend/src/App.tsx`, `frontend/src/FileModal.tsx`):
`eckert4` projection, `combined` view, 2200x1222 render (1100x611 display x RENDER_SCALE 2),
100 kyr/frame (1 step/frame x default 100 kyr/step), `node_density` = `climate_density` = 4.0,
`fluid_density` = 2.0, `wind_model` = `"diagnostic"` (the frontend default -- matters, see
the step breakdown), identity view rotation. `seed = 0`: 19 plates, ~131 K nodes; the
Biome/Combined render grid at `climate_density` 4.0 is **801 x 1601 = 1.28 M points**.

**Method:** drove `render_image.stream_animation_mp4` directly (same call `/world/animate`
makes) for two back-to-back 60-frame runs to get end-to-end timings, then a warmed
12-frame run under `cProfile` for the breakdown. The `combined` view path is
`_render_combined_view` -> `_biome_fields` + `_project_climate_grid` + `_fill_rects`.

### Headline (2026-09-01, `f68fa46`)

| | | vs `c8f9990` |
| --- | --- | --- |
| 60 frames, wall | **278 s and 302 s** on two back-to-back runs (**~4.6 and ~5.0 s/frame**) | 607 s (10.1 s/frame) -- **~2.0-2.2x faster** |
| First frame (render only, no step) | 1.8-2.0 s | 6.5 s |
| Steady state (step + render) | 4.4-7.5 s, mostly ~4.6-4.8 | 9.0-12.9 s |
| Split | **~40% rendering, ~60% simulation step** | was ~62% / ~38% -- rendering took the biggest cut |
| Output | 1.5-1.6 MB H.264/MP4 | 1.9 MB |

The old monotonic upward drift (9 s -> 13 s as the river/lake networks densified) is gone
-- run A was dead flat at 4.4-5.0 s/frame across all 60 frames; run B bounced 4.5-7.5 s
with a couple of mid-run contention spikes, not with simulated age. On an idle box expect
the low end.

### Rendering -- ~1.8-2.0 s/frame

Render-only cost is measured directly by frame 0 (no step): 1.80 s (run A), 2.04 s (run B).
Per-hotspot figures below are from the 12-frame `cProfile` run (profiler inflates absolute
time; ratios hold): `render_png` cumulative 1.97 s/frame there.

| Hotspot | Cost/frame | What it is |
| --- | --- | --- |
| `render_image._biome_fields` | **~1.3 s -> ~0.5 s** (fix 6) | a fresh `cKDTree(all_points)` + a `query` of all 1.28 M grid points, then `hydrology.sample_is_ocean` (a *second* fresh-tree full-grid query); 3x `_bilinear_resample` ~0.06 s; `climate.compute_climate_cached` is a cache hit. Fix 6 put `workers=query_workers(...)` on both full-grid queries -- microbench at the profiled grid size: `_biome_fields` end-to-end ~1.30 s -> ~0.53 s |
| `hydrology.sample_is_ocean` | **~0.38 s -> parallel** (fix 6) | was the single biggest line item -- a single-threaded fresh-tree full-grid `query`. Called once per render (full grid) and once per step; the render call now runs `workers=`-parallel (the per-step call's `query_xyz` is the coarse climate grid, also now parallel above its cutoff) |
| `biomes.smooth_biome_field` | ~0.2 s | Koppen classification + the `_neighbour_vote` boundary-cleanup pass |
| `_encode_image` + `GaussianBlur` | ~0.1 s | PNG encode ~0.07 s, blur ~0.03 s, both vectorized |
| `_project_climate_grid`, `_fill_rects`, `projections._eckert4_theta` | **~0** | fixes 1-3 below. Geometry served from the ring cache (0.000 s), `_fill_rects` njit'd (~0.014 s/call), the Eckert solve runs on ~801 unique latitudes (~0.0001 s). All three fell off the profile entirely. |

### Simulation step -- ~2.8-3.0 s/step

Steady-state frames run ~4.6-4.8 s (run A), so ~2.8-3.0 s of that is the step once the
~1.8 s render is subtracted. Per-hotspot figures from the same 12-frame `cProfile` run
(`step_world` cumulative 3.26 s/step there).

| Hotspot | Cost/step | Note |
| --- | --- | --- |
| `erosion.apply_erosion` | **~1.7 s -> ~1.5 s (fix 7) -> ~1.2 s (fix 8)** | `climate.compute_climate` fresh (~0.59 s), `hydrology.compute_hydrology` (~0.48 s), `_spread_coastal_leveling` (**~0.27 s -> ~0.02 s**, fix 8 -- was a k=48 k-d tree query over ~110 K "sources" that are really every ocean node; now pre-filtered to the few thousand within `INFILL_RANGE_RAD` of a sink), `hydrology._compute_basin_spill` (**~0.16 s -> ~0.02 s**, fix 7 -- was a pure-Python priority-flood, ~0.75 M `heapq.heappop` across the profiled run; now an njit heap kernel), `_coastal_openness` (**~0.15 s -> ~0.10 s**, fix 8 -- counts the smaller non-open neighbour set and subtracts), `lakes.build_lake_hierarchy` (~0.09 s) |
| plate movement (`shift` + `deform`) | ~1.4 s | `torque.gather_boundary_force_inputs` cumulative ~0.5 s but own time only ~0.05 s -- fix 4 landed, so its ~44 calls/step each hit a *cached* per-plate tree; what's left is `plates._plates_within` / `get_neighbours` (~0.2 s) and the `tree.query` itself. `lithosphere_plate._grow_or_shrink_line_for_deform` ~0.23 s. |
| `PlateWithLines.all_points_and_elevation` | **~0.05 s** | was ~0.26 s/step of per-node `geometry.latlon_to_xyz` (~234 K calls). Fix 5 landed: `PlateWithLines._get_world_points` builds the whole plate's world-space cloud in one `local_xyz` + one frame rotation and caches it; `latlon_to_xyz` is down to ~14 K calls / ~0.13 s across the 11 steps, and almost every `all_points_and_elevation` call in a step is a warm-cache read. |
| `world._advance_fluid_dynamics` | **0 (no-op)** | `World.wind_model` defaults to `"diagnostic"` (`frontend` `DEFAULT_WIND_MODEL`), and `_advance_fluid_dynamics` early-returns unless it is `"cfd"`. The old "runs in background threads" note only held when CFD wind was active. Diagnostic wind is rebuilt inside each `compute_climate`. |

The 12-frame `cProfile` still shows ~16 K short-lived thread joins (~27 s cumulative but
overlapped, off the critical path) -- these are scipy `cKDTree.query(workers=-1)` pools
from the many k-d tree queries across `erosion` / `climate` / `torque`, not fluid dynamics,
which no longer threads at all.

### Cheap -- measured, not worth touching

- PNG-encode then `Image.open`-decode round-trip inside `stream_animation_mp4`: **~0.04 s/frame**.
- libx264 encode + muxing: negligible at this frame count and 5/2 fps.

### Fixes

Landed since the original profile (each measured on the box it was written on):

1. ~~**Cache `_project_climate_grid` across frames + `lru_cache` `_biome_grid`.**~~ **Done**
   (`_PROJECT_GRID_CACHE` ring keyed on grid size + projection + rotation bytes + dimensions +
   padding; `_biome_grid` is `lru_cache`d). `_render_combined_view` 5.25 s (cold) -> 2.85 s
   (warm) per frame, ~2.4 s/frame saved.
2. ~~**Numba-jit `_fill_rects`.**~~ **Done** (`_fill_rects_kernel`, `@njit(cache=True)`,
   serial -- overlapping cells would race under `prange`). ~0.78 s -> ~0.014 s per call at
   the profiled grid size, run twice per render -- ~1.5 s/frame saved.
3. ~~**`_eckert4_theta`: fewer Newton iters + solve on unique latitudes.**~~ **Done**
   (`np.unique` on the flattened latitudes so ~801 distinct values are solved and scattered
   back; loop caps at 12 iters, breaks below 1e-14 correction -- ~5 in practice; bit-exact vs
   a 60-iter reference). Full-grid `eckert4` ~380 ms -> ~22 ms, ~8 passes/render -- roughly
   2.7 s/frame saved.
4. ~~**Build the neighbour `cKDTree` once per step in `torque`.**~~ **Done**
   (`Plate.get_node_kdtree`, a per-plate tree cached/invalidated in lockstep with the
   bounding-polygon caches; `gather_boundary_force_inputs` queries each neighbour's cached
   tree and keeps the elementwise-nearest via an argmin. Bit-exact vs the old combined-tree
   path.) `gather_boundary_force_inputs` cumulative 6.29 s -> 3.78 s over 6 steps; its own
   `tottime` is now near zero.
5. ~~**Vectorize `all_points_and_elevation` / `elevation_lines.world_xyz`.**~~ **Done**
   (`PlateWithLines._get_world_points` concatenates every non-empty line's `(phi, theta)` and
   runs one `local_xyz` + one frame rotation over the whole plate, cached in
   `_world_points_cache` and invalidated in lockstep with the bounding-polygon / node-kdtree
   caches -- on `rotate` or a node-set change, never an elevation-only edit.
   `all_points_and_elevation` returns that cached array paired with a fresh
   `collect("elevation")`. Matches the old per-line path to ~1 ULP.) ~0.26 s/step -> ~0.05
   s/step; microbench (9-plate world, `node_density` 4, 300 full-cloud sweeps) 2.94 s old ->
   1.24 s vectorized-cold -> 0.10 s warm cache.

6. ~~**`workers=` on the full-grid render k-d tree queries.**~~ **Done.** `_biome_fields`'
   `cKDTree(all_points).query(flat_xyz)` and `hydrology.sample_is_ocean`'s
   `cKDTree(hydro.points).query(...)` each built a fresh tree and ran a **single-threaded**
   query over all 1.28 M grid points every frame -- together ~1.0 s of the ~1.9 s render.
   Both now pass `workers=query_workers(len(flat_xyz))` (`-1` at this size), matching how
   `torque` / `climate` / the hydrology neighbour-index build already parallelize. Same
   one-liner applied to the sibling full-grid resamples on the same render path:
   `render_image._resource_fields`, `_render_geomorph_view`, `_render_elev_reason_view`, and
   `coastline._lake_mask_on_grid`. Microbench (seed 0, `node_density`/`climate_density` 4,
   801x1601 grid, 6-core Apple Silicon): `_biome_fields` end-to-end **1.30 s -> 0.53 s**
   (~0.77 s/frame), nearest-node output bit-identical (deterministic query). Speeds every
   `combined`/`biome`/`geomorph`/`resources` re-render, not just animation.

   ~~The tree itself is still rebuilt per render (it only changes on a step) -- caching it is
   a further, separate win left for later.~~ **Done** (`World.node_kdtree_cache`, populated by
   `render_image._node_cloud_and_tree` and reset at the top of `step_world` -- the node cloud
   is fixed between steps and no render runs mid-step; dropped on load like the other derived
   caches). Turned out smaller than the "~0.5-0.8 s" this note implied: at `node_density` 4
   (~131 K nodes) the `cKDTree` *build* is only ~20 ms -- the ~165 ms grid `query` is the real
   cost and fix 6 already parallelized that. Still worth it because a single
   `combined`/`elevation` render runs several separate resamples off one tree and an animation
   re-renders every frame without moving a node: measured ~35 ms/render on a 5-render
   view-switching sequence (4.61 s -> 4.44 s), render output bit-identical.

   Not cached: `hydrology.sample_is_ocean`'s `cKDTree(hydro.points)` (built once per
   `_biome_fields` / `_resource_fields` call). It keys off `world.hydrology_cache`, which is
   already per-step, so a `_ocean_tree` lazy attribute on `HydrologyFields` would be the same
   pattern for another ~20 ms/frame -- left as a small follow-up.

7. ~~**`hydrology._compute_basin_spill`'s `heapq.heappop` priority-flood** -- pure Python,
   ~0.21 s/step and grows with basin count. Candidate for a numba kernel or a
   `scipy.ndimage`-based watershed.~~ **Done** (`hydrology._basin_spill_kernel`, `@njit(cache=True)`).
   The k-NN graph isn't a grid, so `scipy.ndimage` watershed doesn't apply -- instead the
   multi-source minimax Dijkstra now runs over an array-backed binary min-heap, ordered
   lexicographically on `(key, node)` to match Python `heapq`'s tuple comparison. Every heap
   entry has a distinct `(cost, node)` key (`cost[j]` only changes on a strict decrease), so
   the min is unique and pop order is identical regardless of heap internals -- **bit-exact**
   with the old path (verified over every spill call across a 4-step seed-0 run, `cost` and
   `spill_target` both `array_equal`). Microbench (seed 0, `node_density` 4, ~131 K nodes):
   **~162 ms -> ~18 ms per call** (~9x), ~0.14 s/step off `apply_erosion`.

8. ~~**`erosion._spread_coastal_leveling` (~0.27 s/step) and `_coastal_openness` (~0.15 s/step)**
   -- the next tier down in `apply_erosion` after climate/hydrology.~~ **Done** -- neither was
   a Python neighbour sweep after all; both were dominated by one oversized scipy k-d tree
   query.
   - `_spread_coastal_leveling`: `source_amount` is `> 0` at ~110 K nodes (every ocean node
     carries a sliver of redirected submarine/coastal spoil), but the fill sinks are a thin
     coastal band (~800 nodes at the profiled size). The `cKDTree(points[sink_idx]).query(
     points[source_idx], k=48, workers=-1)` was therefore ~110 K x 48 -- **~85 ms**, and its
     `np.add.at` / `reduce` tails scaled with it too. Now a cheap `sink_tree.query_ball_point(
     sources, INFILL_RANGE_RAD, return_length=True)` pre-filter drops the source set to the
     ~2.7 K within reach of any sink before the k=48 query; the out-of-range sources keep
     their full amount, which is exactly the `any_reachable is False` branch that handled
     them before -- **bit-exact** (verified `array_equal` on both fields across a 4-step
     seed-0 run). ~197 ms -> ~19 ms per call.
   - `_coastal_openness`: two `query_ball_point(..., return_length=True)` radius counts, one
     over all nodes (`total`), one over the ~111 K open-ocean nodes. The second now runs over
     the ~20 K *non*-open nodes instead and recovers `ocean_count = total - non_open_count`
     (every node is in exactly one set, self included, so the subtraction is exact) -- same
     `ocean_count / max(total, 1)`, **bit-identical** output, ~6x smaller tree + query.
     ~0.15 s -> ~0.10 s per call (the `total` count is the irreducible half).

Still open, roughly in order:

9. ~~**`climate.compute_climate` (~0.59 s/step) and `hydrology.compute_hydrology` (~0.48 s/step)**
   -- with the coastal passes and the basin spill handled, these two are now what
   `apply_erosion` spends its time on. Not yet broken down.~~ **Broken down** (2026-09-01,
   `main` at `391cd71`). Re-measured directly -- 15 warmed `compute_climate` + `compute_hydrology`
   calls at the frontend defaults (seed 0, `node_density` 4, `climate_density` 4 -> **360 x 720
   grid**, ~131 K nodes), driving exactly the pair of calls `apply_erosion` makes, on 6-core
   Apple Silicon, unprofiled wall clock. Current cost is **climate ~0.50 s/call**, **hydrology
   ~0.23 s/call** -- hydrology already roughly halved vs the ~0.48 note (fix 7's njit basin-spill
   postdates that figure), climate a touch under ~0.59.

   Per-region wall (wrapper timers on the named sub-functions; nested regions indented, so a
   parent's number includes its children; `workers=-1` k-d tree pools spin up inside several
   of these and overlap, so the column over-counts the true critical path -- treat it as a
   ranking, not an additive budget):

   | `compute_climate` region | ms/call | what |
   | --- | --- | --- |
   | `compute_ocean_currents` | **~205** | `_land_swirl_current` ~120 (a fresh `cKDTree` over the land grid cells + a full-grid nearest query, then per-cell `np.cross`/`np.stack` over `(H, W, 3)` to build the tangent frame, and `_lon_grid` rebuilt twice); `_smooth_along_coast` ~59 (8 Jacobi passes, each 4x`np.roll` on `u` and `v` plus a re-deflect) |
   | `_sample_elevation_and_crust` | **~185** | a fresh `cKDTree` over the ~131 K node cloud + a full-grid query (259 K points), then `hydrology.sample_is_ocean` -- **a second fresh ~131 K-node tree** (last step's cached points) with its own full-grid query, ~86 of the 185 |
   | `biomes.smooth_biome_field` | ~80 | Koppen classification + `_neighbour_vote` boundary cleanup |
   | `compute_humidity` | ~69 | `_humidity_zonal_sweep` ~48 -- a **Python loop over all 720 columns x 2 laps**, each pass ~15 vectorized `(H,)` ops; a sequential scan along the sweep direction, not vectorizable as-is |
   | `compute_wind` | ~49 | `_sample_at_offset` x~11 (~30 total -- mountain-tangent pick + 5-step wake lookback, each a closed-form lat/lon->row/col gather with a `.tolist()`), `_smooth_field` x2, `_mountain_deflection` / `_mountain_wake_factor` |
   | `compute_moisture_flux_convergence` | ~39 | banding noise + flux-divergence sweep |
   | `_vegetation_transpiration_source` | ~28 | `biomes.classify_biomes` of last step's snapshot |
   | `advect_ocean_temperature`, `_build_grid` | ~7, ~6 | -- |

   | `compute_hydrology` region | ms/call | what |
   | --- | --- | --- |
   | `_build_neighbor_graph` | **~71** | a fresh `cKDTree(points)` (~131 K, `balanced_tree=False`) + one batched k=9 self-query. A **third** rebuild this step of essentially the same node-cloud tree `_sample_elevation_and_crust` already built |
   | `lakes.step_lakes` -> `build_lake_hierarchy` | **~58** | `_catchment_roots` ~22 (Python pointer-chase with path compression over all N nodes) + Python union-find (`find`, ~3.5 M calls across the 12-frame profile) over the catchment-boundary edges + a `_make_leaf` per catchment |
   | `_lake_component_sizes` -> `lake_components` | **~42** | a *separate* Python union-find over the whole k-NN graph (`is_lake.tolist()`, `neighbor_idx.tolist()`, nested loop), only to get each node's lake-component node count |
   | `_compute_basin_spill` | ~40 | njit heap kernel (fix 7) is ~20 of this; the rest is `neighbor_idx` marshalling into the kernel |
   | `route_downstream` (x2 paths) | ~28 | downstream flux accumulation |
   | `connected_ocean_mask` | ~24 | `coo_matrix` + `scipy` `connected_components` over the below-sea-level subgraph |
   | `_compute_flow_direction` x2 | ~22 | water + ice flow targets |

   Follow-ups, split out as their own items below (10-13).

   - Note on the sibling spreads: `_spread_marine_sediment` has the same ~110 K-source shape
     `_spread_coastal_leveling` did, but its targets are the *entire* ocean node set, not a
     thin band, so the source pre-filter doesn't apply -- nearly every source has an in-range
     lower-ocean target. `_spread_beach_sediment`'s sources are just river-mouth nodes (small).
     Neither showed up as a hotspot; left alone.

10. **Build the node-cloud `cKDTree` once per step and share it across climate + hydrology.**
    `_sample_elevation_and_crust` (climate) and `_build_neighbor_graph` (hydrology) each build a
    fresh `cKDTree` over the identical ~131 K-node `node_cloud` every step (~20-30 ms each), and
    the render path builds a third. `World.node_kdtree_cache` already exists for the render side
    (fix 6) and is reset by `step_world`; erosion's per-step `node_cloud` is that same cloud.
    Populate/read that cache from the two erosion-path callers too -> ~40-60 ms/step. (`sample_is_ocean`
    genuinely can't share -- it queries *last* step's cached points against a tree of *this*
    step's, on purpose -- but it can at least take `balanced_tree=False`.)

11. **numba the lake pipeline's Python union-find passes.** `lakes._catchment_roots`,
    `build_lake_hierarchy`'s boundary-merge loop, and `hydrology.lake_components` are ~100 ms/step
    of pure-Python pointer-chasing and `.tolist()` loops -- the same shape as the basin-spill
    priority-flood fix 7 replaced with an njit kernel. `_lake_component_sizes` -> `lake_components`
    in particular is a whole second union-find over the k-NN graph that runs every step purely
    for per-node component sizes.

12. **`climate._humidity_zonal_sweep` (~48 ms/step) is a 720-iteration Python column scan** (x2
    laps). It's a strict left-to-right recurrence along the sweep direction, so it can't be
    vectorized away, but it's a clean candidate for an njit prefix-scan kernel (fix 7 pattern).

13. **`climate._land_swirl_current` (~120 ms/step).** Fresh land-cell `cKDTree` + full-grid
    query, then a per-cell `np.cross` / `np.stack` tangent-frame build over `(H, W, 3)`. The
    `east`/`north` tangent frame is a pure function of the fixed grid -- precompute it once
    (module-level, keyed on `(height, width)`) instead of every call; `_lon_grid` is currently
    built twice per call on top of that.

## Stepping + rendering: cKDTree construction vs. query, and three alternative
## representations (2026-09-13, commit `88b4db6`)

**Measured:** commit `88b4db6` (`main`, post fault-based tectonics, frontier gap-fill,
lake/sea tiers, and the eustasy connectivity fix -- all landed after the profile above).
10-core Apple Silicon (M4) -- a different, faster box than the 6-core one the `f68fa46`
numbers above were measured on, so absolute wall-clock isn't directly comparable step-for-
step; the *proportions* (build vs. query, which call sites dominate) are what this section
leans on. Python 3.14.6, backend `.venv`.

**Inputs:** identical to the section above -- seed 0, `node_density=4.0`,
`climate_density=4.0`, `fluid_density=2.0`, `eckert4`/`combined`, 2200x1222, 100 kyr/step,
diagnostic wind. 19 plates, ~130,600 nodes; the combined-view render grid is 1,282,401 points.

**Method:** a `cKDTree` subclass (`InstrumentedKDTree`) swapped into every already-imported
`app.*` module's own `cKDTree` name (each does `from scipy.spatial import cKDTree`, binding
its own reference at import time, so this has to happen post-import, per-module) that records,
per call site (source file:line:function), construction count/time/point-count and
query/query_ball_point count/time/point-count separately. Driven against 9 real `step_world`
calls + 10 real `render_png` calls (same call `stream_animation_mp4` makes), then separately
under `cProfile` for the usual hotspot breakdown. This directly answers the question the three
alternative representations below all hinge on: is tree *construction* or tree *query* the
actual cost, and by how much.

### Headline: query dominates construction, roughly 2.6:1

Over the 9-step + 10-render run: **8,772 ms total in `.query`/`.query_ball_point`** vs.
**3,314 ms total in tree construction** (`cKDTree.__init__`), across every `cKDTree` built or
queried anywhere in plates/torque/faults/erosion/hydrology/climate/render. Per-tree
construction is cheap at this node count -- the priciest single builder
(`lithosphere_plate._claim_adjacent_territory`'s ~42K-point neighbour-cloud tree) is ~4.7 ms/
call; a full ~130K-node tree is ~15-20 ms to build. Nothing here is spending its time in
median-finding -- it's spending it walking the tree many, many times, and (see below) often
building the *same* tree redundantly before doing so.

Wall-clock: step ~2.9-3.1 s, render ~1.5 s at these settings on this box (render is down from
the ~1.8-2.0 s in the section above -- fix 6's `workers=` parallelism holding up; step is
roughly flat despite fix 7-9 landing, because fault generation, frontier gap-fill, and the
eustasy connectivity solve -- all added since -- backfilled the time those fixes freed).

### New since the last profile: redundant per-step tree construction

The single most actionable finding, and the direct, evidence-backed answer to "what if every
continent were a kd-tree, rebuilt every step": **that's already most of the architecture**
(`Plate.get_node_kdtree`, cached and invalidated only on `rotate()`/a node-set change -- see
fix 4 above), **but several call sites don't use the cache that already exists, and rebuild
an equivalent tree over the same ~130K-node cloud from scratch, once per step, independently
of each other:**

| Call site | Build cost | What it builds |
| --- | --- | --- |
| `climate._sample_elevation_and_crust` | ~15.6 ms/step | reads `World.node_position_tree_cache` if warm, else builds+populates it -- this one's already fixed (5ec022c) |
| `erosion.compute_slope` | ~17.3 ms/step | its own fresh `cKDTree(points)`, immediately `.query(points, k=5)` |
| `erosion._earthquake_erosion_multiplier` | ~15.5 ms/step | its own fresh `cKDTree(points)`, then a handful of tiny `query_ball_point` calls around each epicentre |
| `plates.compute_node_overlap` | ~15.5 ms/step | its own fresh `cKDTree(concatenated clouds)`, `.query_pairs(...)` |
| `hydrology._build_neighbor_graph` | ~5.9 ms/step | its own fresh `cKDTree(points, balanced_tree=False)`, `.query(points, k=9)` |
| `faults.generate_boundary_faults` | ~5.9 ms/step | its own fresh `global_tree` over all points |
| `erosion._route_wind_deposit` | ~2.7 ms/step | its own fresh tree over all points |
| `erosion._spread_marine_sediment` | ~5.3 ms/step (over the ocean subset) | its own fresh tree |
| `gaps._existing_node_tree` | ~15 ms on the (rare) steps it runs | its own fresh tree |

None of these are wrong -- every one gets the right answer -- but together they're
**~85-105 ms/step of pure `cKDTree.__init__` spent rebuilding what is, modulo ordering, the
same point cloud `World.node_position_tree_cache` already holds once climate.py has run this
step.** This is exactly item 10 above, except it turns out to be a wider problem than
originally scoped: fix 5ec022c only wired the cache through `climate.py`; `hydrology.py`,
`erosion.py`, `plates.py`, `faults.py`, and `gaps.py` each independently reinvented it. The
fix is mechanical and bit-exact (same points, same tree semantics, only the constant-factor
build cost changes) -- thread `world.node_position_tree_cache` through these call sites the
same way `climate.py` already does, falling back to a fresh build only where the query shape
genuinely differs (e.g. `_spread_marine_sediment`'s ocean-only subset can't reuse the
full-cloud tree, but could still build its subset tree with `balanced_tree=False,
compact_nodes=False` like its siblings already do).

### New since the last profile: two call sites reinvent torque.py's neighbour-tree fix, worse

`faults._plate_stress` (`~44K`-point concatenated-neighbour tree) and
`lithosphere_plate._claim_adjacent_territory` / `_fill_corner_notch_frontier` (`~42K`-point
concatenated-neighbour tree) each explicitly comment that the tree "is rebuilt from scratch
every call, a fresh combined-neighbour point set each time, so nothing to cache across
calls" -- which was exactly torque.py's problem before fix 4, and fix 4's actual solution
(query each neighbour's own cached `Plate.get_node_kdtree()` and keep the elementwise
nearest via `argmin`, rather than concatenating every neighbour into one fresh tree) applies
here unchanged, since it's the identical shape of computation ("nearest cross-plate boundary
node"). Measured cost of the current, uncached version:

| Call site | Calls/9 steps | Build | Query | Combined/step |
| --- | --- | --- | --- | --- |
| `faults._plate_stress` | 171 (~19/step) | 252 ms | 249 ms | ~56 ms/step |
| `lithosphere_plate._claim_adjacent_territory` | 153 (~17/step) | 717 ms | -- | ~80 ms/step |
| `lithosphere_plate._fill_corner_notch_frontier` | 153 (~17/step) | 714 ms | 399 ms | ~124 ms/step |

That's **~260 ms/step combined**, and unlike the previous section's fix, this one has an
exact, already-proven-correct blueprint sitting in the same codebase (`torque.
gather_boundary_force_inputs`) rather than needing new design work.

### Still the biggest single line item: `hydrology`'s ocean/sea nearest-node resample

`hydrology._nearest_hydro_node_idx` (backing `sample_is_ocean`/`sample_is_sea`) is the single
largest query-time consumer measured outside `torque`: 1,315 ms of query time across just 17
calls (some at the ~259K-point climate-grid resolution, some at the ~1.28M-point render-grid
resolution) plus 264 ms rebuilding `cKDTree(hydro.points)` from scratch each of those calls.
This is exactly the follow-up the previous section's fix 6 flagged and left undone ("`a
_ocean_tree` lazy attribute on `HydrologyFields` would be the same pattern... left as a small
follow-up") -- still true, still worth doing, and now clearly the biggest single win left in
this family: `HydrologyFields` already keys off `world.hydrology_cache` (per-step), so caching
`cKDTree(hydro.points)` there is a direct lazy-attribute add, no invalidation logic beyond what
already exists.

### Where the render-side query time actually goes

`render_image._biome_fields`'s `tree.query(flat_xyz)` (773-1,232 ms across 8-10 renders, no
stepping between) is the single biggest render cost, and its tree is already the cached,
shared one (fix 6/render_image._node_cloud_and_tree) -- there's no more construction cost left
to cut here. What's left is the query itself: 1.28M grid points, each asking "which of ~130K
moving nodes is nearest," `workers=-1` parallel. This is the same shape of cost as
`_sample_elevation_and_crust`'s climate-grid resample and `distance_from_land_approx`'s
coastline-distance resample -- all three are "fixed grid asks a moving point cloud who's
closest," and together they're roughly a third of all measured query time. See the HEALPix
discussion below for why this specific shape, not the other two, is the one worth
re-architecting rather than just caching harder.

### Answering the three alternative representations

**1. "Every continent as a kd-tree, rebuilt at every step."** Substantially already true --
`Plate.get_node_kdtree()` is exactly that, and it's cheap: building a ~7,000-node per-plate
tree or even the full ~130,600-node world tree costs single-digit-to-tens of milliseconds
(measured above). The actual waste this profile found isn't "we don't rebuild enough" or
"rebuilding is slow" -- it's that **several call sites rebuild an equivalent tree redundantly
within the same step** instead of sharing the one `World.node_position_tree_cache` (or a
per-plate `get_node_kdtree()`) already computed this step. The fix is "cache once per step
and thread it through more callers," which this codebase has already proven out twice (fix 4
for `torque.py`, fix 6/5ec022c for `climate.py`/render) -- it just hasn't reached
`hydrology.py`, `erosion.py`, `plates.py`, `faults.py`, or `gaps.py` yet. Rebuilding *across*
steps, by contrast, is unavoidable and correct: every plate rotates rigidly every step, so
every node's position genuinely changes and a stale tree would silently return wrong nearest-
neighbours.

**2. "Stored natively as an R-tree, allowing cheaper construction of kd-trees."** The evidence
here argues against it, on both halves of the premise:
- *Construction isn't the bottleneck to begin with* -- query time outweighs construction
  time 2.6:1 in this profile, and per-tree construction is already a few to tens of
  milliseconds even at ~130K points. There's no large "construction cost" left for a
  different data structure to undercut.
- *An R-tree's actual advantage -- amortized incremental insert/delete instead of a full
  rebuild -- doesn't apply to this workload.* That advantage pays off when most points are
  stable between rebuilds and only a few move (a particle simulation with local motion, a
  scene graph with occasional edits). Here, every plate's node cloud undergoes a **rigid
  rotation** every step (`Plate.rotate`) -- every single node's coordinates change, even if by
  a small angle -- so there is no stable majority to amortize an incremental update against;
  the correct operation is "throw away the old tree, build a new one from the rotated
  positions," which is exactly what happens today. Bulk-loading an R-tree (STR-packing) is
  also `O(n log n)`, the same order as `cKDTree`'s own median-split build -- nothing here is
  asymptotically cheaper.
- R-trees are optimized for bounding-box/range queries over extended geometries (their whole
  reason to exist); this workload is uniformly nearest-neighbour-over-point-clouds, which is
  `cKDTree`'s strong suit -- a single compiled call, batched, `workers=-1`-parallel. A Python
  `rtree`/libspatialindex binding has no batched multi-query equivalent and would almost
  certainly lose on the actual measured hot path (`tree.query(1.28M points)` in one call)
  even before counting construction.

Net: an R-tree would trade a data structure well-suited to this workload for one that isn't,
in exchange for solving a construction-cost problem that isn't actually the bottleneck.
Recommend against.

**3. "Elevation points fixed to HEALPix positions via a lookup table."** This is the one idea
of the three that changes the *asymptotic* shape of the actual bottleneck rather than just its
constant factor -- and this codebase has already built and validated the exact mechanism, just
for a different field. `healpix_grid.py` + `fluid_dynamics_healpix.py` (the opt-in `"cfd"`
wind model) replace `atmosphere_cfd.py`'s equirectangular grid with a fixed HEALPix pixel set
and an `ang2pix_nest_scalar` njit kernel -- an **O(1) per-point lookup with no tree at all**,
built specifically because (per that module's own docstring) an astropy-backed nearest-pixel
lookup "was `semi_lagrangian_advect`'s single costliest call." Two different places this idea
could go, with very different cost/benefit:
- *Snapping the plate/elevation node mesh itself to fixed HEALPix pixels* fights the physics
  rather than helping it: `ElevationLine` rows exist specifically to track continuous rigid
  rotation and deformation (row spacing, contested/open boundary claiming, stretch/shrink at a
  row's ends -- see `lithosphere_plate.py`), and every one of those operations assumes
  positions that move by arbitrary sub-pixel amounts each step. Quantizing to the nearest
  fixed pixel every step is itself a nearest-neighbour query (`ang2pix` per node is O(1), but
  *deciding what value belongs in each pixel as the underlying node moves through it* still
  needs the same resample this profile is trying to avoid) -- it relocates the cost rather
  than removing it, and would mean rebuilding plate tectonics as an Eulerian (fixed-grid)
  simulation instead of the current Lagrangian (moving-point) one. Not a data-structure swap;
  a different simulation.
- *Applying it to the render/climate/hydrology resample grids* is the actually promising
  version, and it's the direction this profile's own numbers point at: `render_image.
  _biome_fields`, `climate._sample_elevation_and_crust`, `hydrology._nearest_hydro_node_idx`,
  and `world.distance_from_land_approx` are all doing `cKDTree(moving_node_positions).query
  (fixed_grid_points)` -- the *gather* direction, `O(M log N)` for `M` grid cells and `N`
  nodes -- and together they're roughly a third of all measured query time. If that fixed
  grid were HEALPix (as `fluid_dynamics_healpix.py`'s already is) instead of equirectangular,
  the resample could run in the *scatter* direction instead: `ang2pix(node_positions)` is an
  `O(1)`-per-node, vectorized/njit lookup -- no tree, no `log N` -- that says which fixed pixel
  each of the ~130,600 moving nodes now falls into; splat elevation/lake/channel data into
  those pixels (keeping the same "nearest/owning node wins" tie-break the current resample
  uses) instead of asking each of 1.28M fixed cells "which moving node is closest." That's
  `O(N)` in the node count instead of `O(M log N)` in the grid resolution -- a real complexity
  win, not just a constant-factor one, and it's the only one of these three ideas that is. The
  cost is real too, though: it's a new grid representation for every render/climate/hydrology
  consumer (not just the CFD solver), a new resample/splat implementation, and a rethink of
  everywhere code currently assumes an equirectangular `(lat row, lon col)` grid --
  `projections.py`'s map projections, `elevation_lines`/hex export, and the frontend's own
  pixel mapping all lean on that shape today. Worth scoping as a real project, not a
  profiling-session patch -- but with a working precedent already in this repo to build it
  from.

### Suggested order

1. Thread `world.node_position_tree_cache` through `hydrology._build_neighbor_graph`,
   `erosion.compute_slope`/`_earthquake_erosion_multiplier`/`_route_wind_deposit`,
   `plates.compute_node_overlap`, and `faults.generate_boundary_faults` -- mechanical,
   bit-exact, ~85-105 ms/step.
2. Port `torque.gather_boundary_force_inputs`'s per-neighbour-cached-tree-plus-`argmin`
   pattern into `faults._plate_stress` and `lithosphere_plate._claim_adjacent_territory`/
   `_fill_corner_notch_frontier` -- ~260 ms/step, same proof of correctness fix 4 already did.
3. Cache `cKDTree(hydro.points)` on `HydrologyFields` the way item 10 originally asked for --
   ~150 ms/step-and-render combined, direct lazy-attribute add.
4. Scope HEALPix-ifying the render/climate/hydrology resample grid as its own project --
   the biggest number on this page (~a third of query time), the only complexity-class win of
   the three ideas, and the one with a working local precedent (`healpix_grid.py`).
5. Leave the R-tree idea aside -- no measured bottleneck it would address.

### Items 1-3, done (2026-09-13)

**Measured:** same box/inputs as the section above, before vs. after, via the same
instrumented-`cKDTree` harness (9 steps + 10 renders). Total `cKDTree` construction time
**3,314 ms -> 2,192 ms (-34%)**; total query time roughly flat (8,772 ms -> 8,998 ms, within
this box's own run-to-run noise -- see the caveat below). Wall-clock step/render means were
likewise flat (step 2.91 s -> 2.94 s, render 1.53 s -> 1.59 s): expected, since query -- which
these fixes deliberately didn't touch -- was already ~2.6x the size of construction, so cutting
a third of the smaller half moves the total by only a few percent. The win here is real but
modest; it doesn't change the headline number, only removes waste that was never load-bearing.

- **Item 1, landed as scoped, minus two call sites.** Added `plates.cached_node_position_tree
  (world, points)` (check `World.node_position_tree_cache`, build+populate on a miss, `world=
  None` always builds fresh for a direct-unit-test caller) and threaded it through
  `erosion.compute_slope`, `erosion._earthquake_erosion_multiplier`, `erosion._route_wind_
  deposit`, and `hydrology._build_neighbor_graph`'s `compute_hydrology` call site (all four
  already ran downstream of `world.py`'s per-step `gather_node_positions` call, strictly after
  shift/deform/topology-changes/gap-fill have settled the node set for the rest of the step --
  see `World.node_position_tree_cache`'s own docstring for why that window is safe to share
  across).

  **`faults.generate_boundary_faults` and `plates.compute_node_overlap` were deliberately left
  out**, on closer reading than the original write-up gave them: both run *before* that window
  closes -- `generate_boundary_faults` from `faults.update_faults`, called mid-`step_world`
  before topology changes/gap-fill can still add or remove nodes; `compute_node_overlap` from
  `merge_split.py` during the topology-change pass itself (and separately from `main.py`,
  outside `step_world` entirely). Sharing `world.node_position_tree_cache` from either would
  let a caller earlier in the step populate it with a node set gap-fill or topology changes
  are about to change, and nothing re-invalidates the cache in between -- climate.py's later
  read would then silently resample against a stale, too-small node cloud for the rest of the
  step. Not worth the correctness risk for ~15-20 ms/step apiece; left building their own
  fresh trees.

  Measured: the redundant full-node-cloud builds these four used to do independently
  (`erosion.py:628`/`596`/`701`, `hydrology.py:393`) have disappeared from the per-call-site
  breakdown entirely, replaced by one shared `plates.cached_node_position_tree` build/step
  (~155-160 ms across the 9-step run, i.e. ~17-18 ms/step for the *shared* tree -- the four
  call sites' own previous combined cost was ~370 ms across the same 9 steps).

- **Item 2, half landed, half reverted after measuring a regression.**
  `lithosphere_plate._claim_adjacent_territory` and `_fill_corner_notch_frontier` now share one
  concatenated-neighbour `cKDTree`, built once per `deform()` call (right before both, since
  neither touches any *neighbour's* own nodes -- only self's -- so the neighbour cloud can't
  have changed between the two calls) instead of each rebuilding an equivalent one from
  scratch. Measured: their combined construction cost **1,431 ms -> 745 ms over the 9-step
  run (-48%)**; query cost unchanged (same single combined tree, same single batched query
  each still runs -- the ~469-780 ms range measured for `_fill_corner_notch_frontier`'s own
  query across different runs is this box's own noise, not a code change, since nothing about
  its query path moved).

  `faults._plate_stress` got the *other* fix -- torque.gather_boundary_force_inputs' own
  per-neighbour-cached-tree-plus-`argmin` pattern -- and it measured as a clear **regression**,
  caught by re-profiling rather than assumed correct from the pattern's own pedigree: construction
  cost dropped to ~0 as expected (reusing each neighbour's already-cached `get_node_kdtree()`),
  but query cost grew **roughly 8x** (249 ms -> ~2,080 ms of query time across the 9-step run),
  because this function queries its *own* full `own_points` array once per neighbour in a
  Python loop -- the same shape torque.py's own fix uses, but torque calls this ~200+ times a
  step (so eliminating each call's fresh combined-tree build dominates), while `_plate_stress`
  is called only ~19 times a step (once per plate) -- there, the old combined-tree build was
  cheap enough (~250 ms total) that looping a full-size query per neighbour instead of running
  one query against one combined tree came out net negative. Reverted to the original
  concatenated-tree implementation, with a comment recording the measurement so this isn't
  re-attempted blind later. **The general lesson**: the per-neighbour-cached-tree pattern's
  win depends on how many times a step the *old* combined tree would have been rebuilt, not
  just on "does this look like the same computation as `gather_boundary_force_inputs`" --
  worth checking call frequency before porting it anywhere else.

- **Item 3, landed as scoped.** Added a module-level `_LAST_OCEAN_TREE` cache in
  `hydrology.py` (identity-keyed on the `HydrologyFields` instance, same convention
  `_LAST_NEAREST_HYDRO_NODE` right above it already uses) backing a new `_ocean_node_tree(
  hydro)` helper, and pointed `_nearest_hydro_node_idx` at it instead of building
  `cKDTree(hydro.points)` fresh on every call. Deliberately a plain module global, not a field
  on `HydrologyFields` itself: `world.hydrology_cache` is part of the pickled save-file object
  graph (persistence.py never drops it on load, unlike `node_kdtree_cache`/
  `node_position_tree_cache`), so a `cKDTree` living there would round-trip through every
  save/load file for no benefit. Measured: 17 calls across the 9-step+10-render run now
  trigger only 9 tree builds (one per distinct `HydrologyFields` this step created, reused by
  every resample against it -- climate-grid and render-grid alike -- instead of one build per
  call).

All three verified against the full unit + stress test suite (`pytest unit_tests
stress_tests`, 401 tests, all green including after the `_plate_stress` revert). One test
(`test_faults.py::test_fault_systems_spawn_with_long_master_traces_and_strand_families`)
failed identically on unmodified `main` before any of this section's changes too -- a
pre-existing issue, not caused by this work -- and turned out to already be fixed on
`origin/main` (`9937bf5`, merged as part of PR #116 while this section's changes were in
progress: the lone-fault length assertion's own tolerance was too tight for real
`BEND_MAX_FRACTION` wander on top of the length-clipped nominal trace). Rebasing this work
onto that `origin/main` picked the fix up for free.

## Item 4, scoped: HEALPix-ifying the render/climate/hydrology resample grid (2026-09-13)

Item 4 above was left as "scope this as its own project" rather than a patch, because it's the
one idea of the three that changes an asymptotic complexity class rather than a constant
factor -- which also makes it the one most worth getting the design right on paper before
writing code. This section is that scoping pass: no production code changed. It reads
`healpix_grid.py`/`fluid_dynamics_healpix.py` (the existing, validated precedent) and every
consumer named in the section above (`render_image.py`, `climate.py`, `hydrology.py`,
`world.py`) closely enough to turn "worth doing" into concrete phases, and surfaces the one
technical question (a fill step for scattered-but-sparse pixels) that has to be spiked with
real numbers before any of it is greenlit.

### What's actually in scope

The three-ideas section already separated "snap the physics mesh to HEALPix" (rejected -- fights
`ElevationLine`'s continuous rigid rotation) from "use HEALPix for the *resample* grids"
(promising). Reading each candidate consumer directly narrows that further:

| Consumer | Shape of the cost | In scope? |
| --- | --- | --- |
| `render_image._node_cloud_and_tree` gather -> `_biome_fields`/`_resource_fields`/`_render_grid_arrays` | `cKDTree(~130K moving nodes).query(fixed grid)`, nearest-value lookup only | **Yes** -- the core case |
| `climate._sample_elevation_and_crust` | Same shape, same node cloud, already shares `world.node_position_tree_cache` (item 1) | **Yes** -- do alongside render, same underlying tree |
| `hydrology._nearest_hydro_node_idx` (`sample_is_ocean`/`sample_is_sea`) | Same nearest-value shape, but deliberately queries *last* step's node positions against *this* step's grid (see item 10's own parenthetical) -- see caveat below | **Yes, with a semantics check** |
| `world.distance_from_land_approx` (`geology.py`'s caller) | Needs the **distance itself**, not just the nearest node's value -- `land_kdtree_cache.query(points)` returns `(dist, idx)` and only `dist` is used | **Different problem** -- see below, not a drop-in |
| `render_image._classify_terrain_relief` | `query_ball_point` (radius search, variable-count neighbours), not nearest-value | **Different problem** -- see below |
| `climate.py`'s own native simulation grid (`_build_grid`, `compute_climate`'s internal `(H,W)` physics arrays) | Not a moving-cloud resample at all -- it's climate's own fixed working grid | **Out of scope** -- see below |
| `fluid_dynamics_healpix.py` / CFD wind model | Already HEALPix | N/A, already done |
| `geodesic.py` hex export (`File > Export Hex Grid`) | A separate icosahedral tiling, separate output format, no render-grid resample in its critical path | **Out of scope**, unaffected either way |
| Frontend (`MapCanvas.tsx` `project`/`unproject`) | Operates on PNG bytes + its own lat/lon<->pixel math, never sees the backend's internal grid representation | **Unaffected** -- HEALPix would be purely a backend intermediate, invisible past `render_png`'s return value |

Climate's own simulation grid staying equirectangular is a deliberate exclusion, not an
oversight: `_humidity_zonal_sweep` (item 12) and the banding-noise/flux-divergence sweeps in
`compute_moisture_flux_convergence` are literal left-to-right/pole-to-pole row scans -- they
lean on "row `i` is a line of constant latitude, column `j+1` is column `j`'s immediate
eastward neighbour" being true by construction, which a HEALPix pixel index doesn't give for
free (its neighbour table is unordered/8-connected, not a scan direction). Porting *those*
kernels to HEALPix is a separate, unrelated project (arguably harder, since a prefix-scan
needs a consistent sweep order that nested-scheme HEALPix pixel numbering doesn't provide) --
this scope is only the moving-node-cloud -> fixed-grid *resample* step, which climate.py also
does (`_sample_elevation_and_crust`) but as input gathering, not as its own working grid.

### The one finding that reshapes the whole plan: `resample_to_equirect` is already O(1)/cell

`healpix_grid.resample_to_equirect` (used today by the CFD wind path to hand v2's HEALPix wind
state back to v1's equirectangular `climate.py`) does `pix = grid.ang2pix(lon_grid, lat_grid)`
then `field_pix[pix]` -- an `O(1)`-per-cell **lookup**, not a tree query, because it's going
*from* a HEALPix pixel array *to* arbitrary query points via the njit `ang2pix` kernel, same
direction as the scatter step below. That means the render/climate/hydrology migration doesn't
need to touch `_biome_grid`'s `(H, W)` array shape, `_bilinear_resample`, `biomes.
smooth_biome_field`, PNG encoding, the legend, or the frontend at all -- every one of those
stays exactly as it is today, fed by an `(H, W)` array. Only the step that currently *fills*
that array changes:

```
today:     cKDTree(node_xyz).query(grid_xyz)          -- O(M log N), M=grid cells, N=nodes
proposed:  ang2pix(node_xyz)  -> scatter onto HealpixGrid(npix ~ N)   -- O(N)
           ang2pix(grid_xyz)  -> resample_to_equirect(HealpixGrid)    -- O(M), O(1)/cell
```

This is a materially smaller and safer change than "migrate the render pipeline to HEALPix"
(the more invasive alternative -- painting HEALPix pixels directly via `_fill_rects`, the way
`_draw_climate_vectors_healpix` already does for the wind-vector overlay -- was considered and
set aside for this scope: it would remove the final `O(M)` hop too, but at the cost of porting
`smooth_biome_field`'s neighbour-vote, `_bilinear_resample`, and the legend/blur/encode chain to
an irregular flat pixel list, none of which have a HEALPix equivalent today). Worth revisiting
only if the `O(M)` resample-to-equirect hop itself ever shows up as a cost on its own merits --
nothing in this profile suggests it will, since it was already cheap enough to be the *existing*
production path for wind vectors.

### The open technical question: scatter is sparse, and sparse needs a fill step

The `O(N)` half (`ang2pix(node_xyz)`) is the genuinely new work, and it isn't a drop-in swap for
`cKDTree(...).query(...)` the way `resample_to_equirect` is, for one concrete reason: **nearest-
neighbour query never leaves a query point unanswered, but scattering leaves most pixels
empty.**

- `render_image.py`'s own comment on `grid_spacing_rad` (picking the finer of a fixed 100km and
  `plates.line_spacing_rad(world.node_density)`) establishes that the render grid is
  deliberately sized to match, not exceed, the physics node spacing -- so `HealpixGrid(nside)`
  should be picked so `npix` is the same order of magnitude as `N` (~130,600 nodes at the
  profiled density), not the render grid's own, much larger `M` (1,282,401 for the combined
  view). That's the right target for the `O(N)` scatter half.
- At `npix ~ N`, scattering `N` points into `npix` pixels via a single-valued last-write-wins
  `ang2pix` assignment leaves a real fraction of pixels with **zero** nodes (an occupancy-1
  scatter, even over a fairly even point distribution, doesn't fill every bin -- coastlines,
  recently-created gap-fill regions, and the poles are exactly where this is worst) and a
  smaller fraction with **more than one** (needing a deterministic tie-break -- "nearest to
  pixel center," matching today's exact-nearest-node semantics most closely, is the natural
  choice, not "last node in plate/index order," which would make output depend on plate
  iteration order).
- Every empty pixel then needs a value pulled from somewhere, and the only cheap source is the
  `HealpixGrid.neighbours` `(npix, 8)` table already built by `healpix_grid.build()`: a wavefront
  fill (repeatedly copy a filled pixel's value into its unfilled neighbours, like a multi-source
  BFS) is the same *shape* of algorithm as `hydrology._basin_spill_kernel` (fix 7) -- and should
  be written the same way, an `@njit` array-backed frontier, not a Python loop. Its cost is
  `O(N + npix)` only if the average empty-run length is a small constant (true if nodes are
  roughly evenly spread relative to `nside`); a pathological case (e.g. an early-game world with
  large unclaimed ocean regions, or a `node_density` chosen much coarser than the render's
  `nside`) could make the fill dominate the whole exercise instead of the tree query it's meant
  to replace.

This is the one number this scope doesn't have yet, and it's the one that decides whether the
rest of the plan is worth doing at all: **the fill step's real cost, on this codebase's real
node distributions, has to be measured before committing to phases 1+ below.**

### Two shapes that don't fit the scatter pattern at all

- **`world.distance_from_land_approx`** needs the distance value itself (`geology.py` uses it for
  distance-banded effects), not a nearest-node's attribute. A HEALPix analogue exists in
  principle -- a multi-source Dijkstra/BFS from every land pixel over `neighbours`/
  `neighbour_distance_m`, the exact same shape as the basin-spill kernel again -- but it's a
  genuinely different algorithm from "scatter + fill," not a variant of it, and its accuracy
  depends on how well `neighbour_distance_m`'s chord-projection approximates true geodesic
  distance over many hops (fine for CFD's own single-hop gradient use; unverified over the
  many-hop distances this function actually returns). Left as a separate follow-up, not part of
  this migration.
- **`render_image._classify_terrain_relief`** needs a variable-radius neighbour *set* (`elevation.max() - elevation.min()` over every node within 50km), not a nearest-value lookup. HEALPix's
  fixed neighbour table gives a fixed number of rings, not a fixed real-world radius, so
  "how many rings equal 50km at this `nside`" would need its own derivation and would only be
  approximately right (rings are equal-area, not equal-real-distance near the poles vs equator
  the same way `cKDTree.query_ball_point`'s radius already is exactly right by construction).
  Also left out of this migration; it's a small, occasional cost (only the Elevation view's two
  optional toggles) next to the render grid it's excluded from sharing.

### `hydrology`'s last-step/this-step semantics need a check, not an assumption

Item 10's own note ("`sample_is_ocean` genuinely can't share [the cached tree] -- it queries
*last* step's cached points against a tree of *this* step's, on purpose") still applies here in
a slightly different form: under scatter+fill, `hydro.is_ocean` (an array over *last* step's
node positions) would need to be scattered using *last* step's `ang2pix` assignment, then the
render/climate grid resampled from *this* step's filled HEALPix array using `resample_to_
equirect` against *this* step's `ang2pix`. That's two different `HealpixGrid` populations one
step apart, not one -- workable, but it means `sample_is_ocean`/`sample_is_sea`'s one-step
staleness tolerance (already documented as deliberate) needs to be re-verified empirically
against the HEALPix path rather than assumed to carry over unchanged, since "nearest node" and
"nearest filled HEALPix pixel after a wavefront fill" are not quite the same nearest-neighbour
operation at a coastline boundary, which is exactly where `is_ocean` correctness matters most.

### Suggested phases

0. ~~**Spike only, no production code.**~~ **Done, 2026-09-14 -- see below.** Build a
   throwaway `HealpixGrid(nside ~ node count)`, scatter a real world's node cloud (elevation,
   is_ocean) via `ang2pix` with the nearest-to-center tie-break, run the njit wavefront fill, and
   measure: fill-step wall clock across a few real saves (early-game sparse, late-game dense),
   empty-pixel fraction before/after fill, and end-to-end scatter+fill+resample_to_equirect wall
   clock vs. today's `cKDTree(...).query(...)` at the same effective resolution. This is the
   number that makes or breaks the rest of the plan -- don't move past this phase without it.
   **Result: ~2.9-3.8x faster at steady state, fill cost negligible -- go, with two caveats
   budgeted into phase 1 (see "Phase 0, done" below).**
1. **Done (2026-09-14).** Landed the shared `HealpixGrid` + scatter/fill/resample plumbing as
   a new, optional code path (`World.node_cloud_resample_mode`, default `"kdtree"`, opt into
   `"healpix"`) behind `render_image._node_cloud_and_tree` only -- see "Phase 1 landed" below
   for why `plates.cached_node_position_tree`, also named in this checklist item's original
   text, turned out not to belong in this phase after all.
2. Extend to `climate._sample_elevation_and_crust` (same node cloud, same tree) and `hydrology.
   _nearest_hydro_node_idx` (after the last-step/this-step check above is resolved) together,
   since they already share `world.node_position_tree_cache`/the item-3 ocean-tree cache and
   would share one `HealpixGrid` population per step the same way.
3. Leave `world.distance_from_land_approx` and `_classify_terrain_relief` on their current
   `cKDTree` paths (see above) -- not part of this migration, no shared infrastructure to gain
   from doing them at the same time.
4. Only after 1-2 are shipped and measured on a real animation run: consider flipping the
   default, and only then consider the more invasive "paint HEALPix pixels directly" alternative
   mentioned above, if the remaining `O(M)` `resample_to_equirect` hop turns out to matter after
   all (nothing in this profile suggests it will).

### Why this is scoped as a project, concretely

Each phase above is independently revertable and individually small, but the honest total is
still substantial: a new shared data structure touching five modules
(`render_image.py`/`climate.py`/`hydrology.py`/`plates.py`/`world.py`), one new numba kernel
(the wavefront fill) that needs the same bit-exactness discipline fix 7's basin-spill kernel
went through, a tie-break rule that has to be chosen and justified, and a semantics re-check for
hydrology's cross-step read. That's why item 4 stayed "scope it" rather than "do it" even after
items 1-3 landed in one session -- and why phase 0's spike, not this write-up, is the actual
next action if this gets picked up.

### Phase 0, done (2026-09-14): scatter+fill+resample beats `cKDTree.query` ~3x at steady state

**Measured:** throwaway spike script (not committed -- scatter/fill/resample reimplemented
standalone against `healpix_grid.py`'s real `HealpixGrid.build`/`ang2pix`/`resample_to_equirect`,
not wired into `render_image.py`), backend `.venv`, same box as the sections above. `seed=0`,
`node_density=climate_density=4.0` (this doc's own convention), grid `801x1601 = 1.28M` points
(the Biome/Combined render grid size from the animation-profile section above). `nside=128`
(`npix=196,608`, chosen as the smallest power-of-2 `nside` with `npix >= node count`, per the
write-up's own "same order of magnitude as N" guidance) against a node cloud of ~130.6K nodes.
Wavefront fill: a genuinely new njit kernel (this spike's own, plain multi-source BFS over
`HealpixGrid.neighbours`/`neighbour_valid`, array-backed frontier, no priority queue needed since
every hop costs the same) -- not adapted from `hydrology._basin_spill_kernel` as the write-up
above suggested, because that function no longer exists in this checkout (see caveat below).

**Headline:** at steady state (JIT warm, repeated calls against the same already-built node
cloud), scatter+fill+`resample_to_equirect` consistently beats a warm `cKDTree.query` by
**~2.9-3.8x** -- 5 back-to-back trials against the early-game world: 63.4-102.3 ms end-to-end vs.
216.7-240.6 ms for `cKDTree.query` alone (both distributions overlap run to run, but the ratio
never dropped below 2.35x once past the first call in the process). The fill step itself is
cheap and not the bottleneck the write-up worried it might be: **2.2-4.8 ms, 2-4 BFS rounds**,
0% empty pixels remaining every time, from a **~35% empty-before-fill** starting point (both the
early- and late-game worlds landed at 35.0-35.1% empty pre-fill, effectively identical). Scatter
(nearest-to-center tie-break, plain numpy + a Python loop over ~130K nodes, not yet its own njit
kernel) is the actual larger piece of the proposed path at ~54-82 ms, still well under
`cKDTree.query`'s own 217-240 ms. **This is the number the write-up said would make or break the
plan, and it clears the bar -- recommend proceeding to phase 1**, with the two caveats below
budgeted into that work rather than treated as settled.

**Caveat 1 -- single-sample timings on this box are noisy enough to invert the result.** The
very first measurement taken (early-game world, first call in a fresh process) showed
scatter+fill+resample at 231.8 ms against `cKDTree.query`'s 235.0 ms -- a statistical tie
(1.01x), which would have read as "phase 1 isn't worth it" if taken at face value. Five
repeated trials against the same node cloud (no world regeneration, so this isolates measurement
noise, not a real early-game-specific effect) settled to a consistent ~3x once past that first
call -- first-call cost is dominated by cold caches/page faults/thread-pool warm-up, not the
algorithm. Any phase-1 PR's own before/after numbers should report repeated trials the way the
"Items 1-3, done" section above already does, not a single sample.

**Caveat 2 -- accuracy tolerance vs. `cKDTree` needs to be a real gate, not a formality.**
Nearest-filled-HEALPix-pixel elevation disagreed with the `cKDTree` ground truth at the same
grid resolution by mean 33.8 m (early) / 81.2 m (late), p95 124.3 m / 389.6 m, but **max 4.16 km
/ 8.12 km** -- a handful of pixels disagree by more than almost any single real elevation swing
on this planet. Phase 1's own checklist item ("diffed pixel-for-pixel ... not just visually --
an exact `array_equal`/tolerance check") needs to actually inspect where those outlier pixels
land (a coastline/land-ocean BFS-fill boundary is the likely suspect, not confirmed here) before
picking a tolerance, rather than accepting a mean/p95-level number as sufficient.

**Node count and empty-fraction didn't grow between early- and late-game, contrary to the
write-up's own worry.** The write-up above flagged "an early-game world with large unclaimed
ocean regions" as a plausible pathological case for the fill step. Measured against 0 steps vs.
60 steps (100 kyr/step, ~6 My simulated) at `seed=0`: node count barely moved (130,587 ->
130,573) and empty-before-fill was effectively flat (35.0% -> 35.1%) -- plates tile the whole
globe by construction (Voronoi-based generation), so there's no growing "unclaimed" region over
time in this codebase's model, at least for this seed/density. Worth spot-checking a sparser
`node_density` and/or a sketch-based premade world with a deliberately uneven landmass before
treating the fill step's cost as settled for every world shape -- this run only ever saw one
node-density/seed combination.

**Aside -- two doc-drift findings, unrelated to the numbers above but worth flagging for
whoever picks up phase 1:**

- `hydrology._basin_spill_kernel`, cited both by the write-up above and by this doc's own
  "Why this is scoped as a project" paragraph as the precedent for "the same bit-exactness
  discipline" the wavefront-fill kernel should follow, no longer exists -- `hydrology.py`'s own
  current module docstring says it was removed because it could drift out of sync with
  `lakes.py`'s own, physically-correct spill hierarchy (`lakes.compute_spill_routing`). That's
  worth more than a shrug: it's a second instance of the exact failure mode phase 2's own
  hydrology cross-step semantics check is worried about (a derived resample diverging from its
  authoritative source), and whoever writes the real fill kernel should read that removal's
  history before using the old kernel as a model that isn't there anymore.
- This spike's first cut segfaulted (exit 139, no Python traceback) rather than raising, from a
  numba warm-up/precompile call that sliced `filled`/`values` to a small dummy length while
  still passing the full-size `neighbours`/`neighbour_valid` arrays -- an out-of-bounds read
  inside `nopython` code. Not a production bug (no production kernel exists yet), but a concrete
  footgun for phase 1: a kernel warm-up call needs an internally consistent dummy shape, never a
  slice of one array paired with another at full production size.

### Phase 1 landed (2026-09-14)

Implements `healpix_grid.scatter_node_indices` (HEALPix-pixel-center tie-break, via a single
`np.lexsort` rather than a per-pixel loop), a double-buffered `@njit` wavefront-fill kernel
(`_wavefront_fill_round`/`_wavefront_fill`, fixed lowest-neighbour-slot tie-break, raises rather
than ever returning an unfilled `-1`), and `NodePixelIndex` -- a `cKDTree.query()`-compatible
wrapper so every real `_node_cloud_and_tree` consumer (`_render_grid_arrays`, `_biome_fields`,
`_resource_fields`, and the elev-reason/crust-type/speckle views) needed zero changes. Gated by
a new `World.node_cloud_resample_mode` flag (`"kdtree"` default / `"healpix"`), modeled on
`gap_fill_algorithm` rather than `wind_model` -- backend/API-only for now
(`POST /world/controls`), no Controls-panel entry, since this is a proving-out flag, not yet a
user-facing tuning knob.

**Deviation from this checklist item's original text:** it named both
`render_image._node_cloud_and_tree` and `plates.cached_node_position_tree` as Phase 1's entry
points. Checking the latter's five real call shapes found three are flatly incompatible with a
HEALPix pixel-ownership index -- `hydrology._build_neighbor_graph`/`erosion.compute_slope`
both call `.query(..., k>1)` (self k-NN for a routing/slope graph) and
`erosion._earthquake_erosion_multiplier` calls `.query_ball_point` (a radius search) -- and the
remaining compatible, relevant caller (`climate._sample_elevation_and_crust`) was already this
document's own **Phase 2**, never Phase 1. So Phase 1 only touches
`render_image._node_cloud_and_tree`; `plates.cached_node_position_tree` is untouched, deferred
to Phase 2 alongside `climate._sample_elevation_and_crust`/`hydrology._nearest_hydro_node_idx`
as originally planned.

`_classify_terrain_relief`'s `query_ball_point` radius search (already out of scope per the
table above) still needs a real tree under `"healpix"` mode -- handled at its one call site in
`_render_grid_arrays`, which now builds a real `cKDTree` on demand (cached on the new
`World.node_kdtree_relief_cache`) only when the Elevation view's relief toggles are on, leaving
every other caller on the fast path.

**Numbers, measured on the real implementation** (not phase-0's throwaway script), same config
as phase-0's own spike (seed=0, node_density=climate_density=4.0, 4 steps, ~130.5K nodes,
nside=128/npix=196,608, 801x1601 render grid -- see `backend/stress_tests/test_healpix_
resample.py`):

- **Query speedup: ~8-9x** at steady state (already-built index/tree, 5 repeated trials) --
  notably better than phase-0's own ~2.9-3.8x, likely implementation/box-specific rather than a
  meaningful discrepancy (phase-0's own writeup already flagged single-box timing as noisy).
- **Scatter+fill rebuild cost: ~22-25ms** (5 repeated trials, includes the scatter/tie-break sort,
  not just the fill rounds), 3 rounds to close the gap to zero empty pixels -- comparable to,
  not cheaper than, a `cKDTree` build over the same node cloud (~20ms, per `_node_cloud_and_tree`'s
  own docstring), so the win is entirely on the query side, exactly as phase 0 predicted.
- **Accuracy vs. `cKDTree` ground truth, full 801x1601 render grid:** mean 53.5m, p95 256.9m,
  max 7048.2m -- squarely inside phase-0's own reported ranges (mean 34-81m, p95 124-390m, max
  4.2-8.1km), despite being a different (real, not throwaway) implementation.
- **Where the outliers actually land (phase-0 could only guess):** correlation between
  per-pixel disagreement and local elevation-gradient magnitude is **0.53**, and the worst
  disagreements concentrate in the highest-relief terrain, not randomly -- confirmed by direct
  measurement, not inferred. Coastlines are part of this (a dilated land/ocean boundary mask's
  mean disagreement there is ~10x the non-coastal mean: 435m vs. 44m in one measured run) since
  a land/ocean transition is by definition a steep one, but only explains a minority of large
  (>1km) disagreements on their own (~29% in one measured run) -- rugged interior terrain
  (mountain ridgelines, valley edges) drives the rest. So: high local relief generally, not
  "coastline specifically," is the real, confirmed explanation phase-0's writeup asked for.
- `hydrology`'s own one-step-behind `is_ocean`/`is_sea` semantics (this document's own
  "last-step/this-step" note above) are unaffected by Phase 1 -- `hydrology.py` isn't touched
  this phase at all (deferred to Phase 2 along with everything else hydrology-related).
