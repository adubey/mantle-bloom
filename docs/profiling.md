# Profiling

A running record of where wall-clock time actually goes, measured against a real code
path rather than guessed at. Each section says what was measured, on what commit, with
what inputs, so a later run can be compared like-for-like.

## File > Make Animation (`POST /world/animate`), 60 frames

**Measured:** re-run 2026-09-01, commit `75e8c95` (`main` + `chore/biome-legend-reorder`;
that branch only reorders legend swatches, no effect on the animate path). 6-core Apple
Silicon (2 performance + 4 efficiency cores), Python 3.14.6, backend `.venv`, with VS Code
and the vite/uvicorn dev servers running -- there is real contention, see the variance in
the headline. Supersedes the original run at commit `c8f9990`, 2026-08-31, kept below for
the before/after.

**Inputs** -- the frontend's own defaults (`frontend/src/App.tsx`, `frontend/src/FileModal.tsx`):
`eckert4` projection, `combined` view, 2200x1222 render (1100x611 display x RENDER_SCALE 2),
100 kyr/frame (1 step/frame x default 100 kyr/step), `node_density` = `climate_density` = 4.0,
`fluid_density` = 2.0, `wind_model` = `"diagnostic"` (the frontend default -- matters, see
the step breakdown), identity view rotation.

**Method:** drove `render_image.stream_animation_mp4` directly (same call `/world/animate`
makes) for two full 60-frame runs to get end-to-end timings, then a warmed 9-frame run
under `cProfile` and a warmed 14-frame run with per-phase wall-time timers for the
breakdown. The `combined` view path is `_render_combined_view` -> `_biome_fields` +
`_project_climate_grid` + `_fill_rects`.

### Headline (2026-09-01)

| | | vs `c8f9990` |
| --- | --- | --- |
| 60 frames, wall | **288 s and 369 s** on two back-to-back runs (**~4.8 and ~6.1 s/frame**) | 607 s (10.1 s/frame) -- **~1.7-2.1x faster** |
| First frame (render only, no step) | 2.9-3.7 s | 6.5 s |
| Steady state (step + render) | 4.7-7.8 s | 9.0-12.9 s |
| Split | **~40% rendering, ~60% simulation step** | was ~62% / ~38% -- rendering took the biggest cut |
| Output | 1.6 MB H.264/MP4 | 1.9 MB |

The old monotonic upward drift (9 s -> 13 s as the river/lake networks densified) is gone
-- per-frame time now bounces 4.7-7.8 s roughly in step with machine load, not with
simulated age. The two 60-frame runs differing by 28% is that same contention; on an idle
box expect the low end.

### Rendering -- ~2.0-3.0 s/frame

Per-frame figures from the 9-frame `cProfile` run (profiler inflates absolute time; ratios
hold) cross-checked against the phase-timer run.

| Hotspot | Cost/frame | What it is |
| --- | --- | --- |
| `render_image._biome_fields` | **~1.3 s** | still the bulk. ~0.6 s of its own time is a fresh `cKDTree` over the whole node cloud + a **single-threaded** `query` of all 1.28 M grid points (no `workers=`); `hydrology.sample_is_ocean` ~0.4 s (a second full-grid fresh-tree single-threaded query); 3x `_bilinear_resample` ~0.15 s; `climate.compute_climate_cached` is now a cache hit |
| `biomes.smooth_biome_field` | ~0.2 s | Koppen classification + the boundary-cleanup neighbour vote |
| color composite + `GaussianBlur` + `_encode_image` | ~0.4 s | all vectorized; the blur itself ~0.05 s |
| `_project_climate_grid`, `_fill_rects`, `projections._eckert4_theta` | **~0** | fixes 1-3 below. Geometry served from the ring cache, `_fill_rects` njit'd (~0.015 s/call), the Eckert solve runs on ~801 unique latitudes. All three fell off the profile entirely. |

### Simulation step -- ~3.1-4.0 s/frame

| Hotspot | Cost/step | Note |
| --- | --- | --- |
| `erosion.apply_erosion` | **~1.7 s** | `climate.compute_climate` fresh (~0.6 s), `hydrology.compute_hydrology` (~0.5 s), `_spread_coastal_leveling` (~0.24 s), `hydrology._compute_basin_spill` (~0.19 s, ~1.05 M `heapq.heappop`), `_coastal_openness` (~0.14 s), `lakes.build_lake_hierarchy` (~0.09 s) |
| plate movement (`shift` + `deform`) | ~1.3 s | `torque.gather_boundary_force_inputs` ~0.6 s -- fix 4 landed, so its ~38 calls/step each hit a *cached* per-plate tree; what's left is the `tree.query` itself plus `all_points_and_elevation`'s ~0.26 s of per-node Python (fix 5, still open) |
| `world._advance_fluid_dynamics` | **0 (no-op)** | `World.wind_model` defaults to `"diagnostic"` (`frontend` `DEFAULT_WIND_MODEL`), and `_advance_fluid_dynamics` early-returns unless it is `"cfd"`. The old "runs in background threads" note only held when CFD wind was active. Diagnostic wind is rebuilt inside each `compute_climate`. |

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
   serial). ~0.78 s -> ~0.027 s per call at the profiled grid size, run twice per render --
   ~1.5 s/frame saved.
3. ~~**`_eckert4_theta`: fewer Newton iters + solve on unique latitudes.**~~ **Done**
   (`np.unique` on the flattened latitudes so ~801 distinct values are solved and scattered
   back; loop caps at 12 iters, breaks below 1e-14 correction -- ~5 in practice; bit-exact vs
   a 60-iter reference). Full-grid `eckert4` ~380 ms -> ~22 ms, ~8 passes/render -- roughly
   2.7 s/frame saved.
4. ~~**Build the neighbour `cKDTree` once per step in `torque`.**~~ **Done**
   (`Plate.get_node_kdtree`, a per-plate tree cached/invalidated in lockstep with the
   bounding-polygon caches; `gather_boundary_force_inputs` queries each neighbour's cached
   tree and keeps the elementwise-nearest via an argmin. Bit-exact vs the old combined-tree
   path.) `gather_boundary_force_inputs` cumulative 6.29 s -> 3.78 s over 6 steps.

Still open, roughly in order:

5. **`workers=` on the two full-grid render queries.** Both `_biome_fields`' own
   `cKDTree(all_points).query(flat_xyz)` and `hydrology.sample_is_ocean` build a fresh tree
   and run a **single-threaded** query over all 1.28 M grid points every frame. Pass
   `workers=-1` (the codebase already does this in `torque` / the land k-d tree via
   `plates.query_workers`) and/or cache the tree -- the queried node cloud only changes on a
   step, not between static-camera re-renders. ~0.6-1.0 s/frame, and it speeds every
   `combined`/`biome` re-render too, not just animation.
6. **Vectorize `all_points_and_elevation` / `elevation_lines.world_xyz`** to remove the
   per-node `latlon_to_xyz` calls (~0.26 s/step, ~180 K calls in the 9-frame run).
7. **`hydrology._compute_basin_spill`'s ~1.05 M `heapq.heappop`** -- a priority-flood in pure
   Python; ~0.19 s/step and grows with basin count. Candidate for a numba kernel or
   `scipy.ndimage`-based watershed.
