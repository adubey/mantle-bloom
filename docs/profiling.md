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

9. **`climate.compute_climate` (~0.59 s/step) and `hydrology.compute_hydrology` (~0.48 s/step)**
   -- with the coastal passes and the basin spill handled, these two are now what
   `apply_erosion` spends its time on. Both are called fresh once per step (climate is a cache
   hit at *render* time but not here). Not yet broken down.
   - Note on the sibling spreads: `_spread_marine_sediment` has the same ~110 K-source shape
     `_spread_coastal_leveling` did, but its targets are the *entire* ocean node set, not a
     thin band, so the source pre-filter doesn't apply -- nearly every source has an in-range
     lower-ocean target. `_spread_beach_sediment`'s sources are just river-mouth nodes (small).
     Neither showed up as a hotspot; left alone.
