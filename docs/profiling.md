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
   padding; `_biome_grid` is `lru_cache`d). Measured on the same box: `_render_combined_view`
   dropped 5.25 s (cold) -> 2.85 s (warm cache) per frame, ~2.4 s/frame saved.
2. ~~**Numba-jit `_fill_rects`** (the codebase already JITs `atmosphere_cfd`). The
   1.28 M-iteration Python loop is ~1.5 s/frame.~~ **Done** (`_fill_rects_kernel`,
   `@njit(cache=True)`, serial -- overlapping cells would race under `prange`; the Python
   wrapper still does the vectorized clip/round and normalizes `colors` to the buffer dtype
   once). Microbenchmark on the same box at the profiled grid size (801x1601 points):
   ~0.78 s -> ~0.027 s per `_fill_rects` call, run twice per combined render -- ~1.5 s/frame
   saved.
3. ~~**`_eckert4_theta`: drop `_NEWTON_ITERS` 30 -> ~8** (Newton converges quadratically,
   double precision is reached well before 30) and/or solve on the ~801 unique latitudes
   rather than all 1.28 M points. Independent of fix 1.~~ **Done** (both: `np.unique` on the
   flattened latitudes so the solve runs on the ~801 distinct values and is scattered back;
   Newton loop now caps at 12 iters and breaks once the max correction drops below 1e-14 --
   ~5 iters in practice). Result is bit-exact against a 60-iteration reference. Microbenchmark
   on the same box at the profiled grid size (801x1601): `eckert4` full-grid call
   ~380 ms -> ~22 ms (~17x), ~8 passes/render -> roughly **2.7 s/frame saved**.
4. ~~**Build the neighbour `cKDTree` once per step in `torque`** and share it across every
   plate's `shift`/`deform` instead of rebuilding it ~24x.~~ **Done** (`Plate.get_node_kdtree`:
   a per-plate `cKDTree` over the plate's own node cloud, cached and invalidated in lockstep
   with the bounding-polygon caches -- i.e. on `rotate` / any node-set change, never an
   elevation-only edit -- and built with the same `balanced_tree=False, compact_nodes=False`
   fast-build flags v1's `deform` already uses. `gather_boundary_force_inputs` now queries
   each neighbour's cached tree and keeps the elementwise-nearest via an argmin over
   neighbours instead of concatenating every neighbour's nodes into a fresh tree per call;
   one plate is a neighbour of several others and is queried in both the shift and deform
   pass, so its cloud is treed ~once per step rather than ~24x.) Output is bit-exact against
   the old combined-tree path (`dist`/`direction`/`omega`/`is_oceanic`, verified over a
   10-plate world across the fresh state and three steps). Bench on the same box (10 plates,
   `node_density` 4, 6 steps): `gather_boundary_force_inputs` cumulative 6.29 s -> 3.78 s
   (its own `tottime` 4.49 s -> 0.15 s), ~0.42 s/step, ~3.6 s/step -> ~3.2 s/step wall.
5. ~~**Vectorize `all_points_and_elevation` / `elevation_lines.world_xyz` to remove the 237 K
   per-node calls.**~~ **Done** (both a vectorize and a cache). `PlateWithLines._get_world_points`
   concatenates every non-empty line's `(phi, theta)` and runs a single `local_xyz` + one
   frame rotation over the whole plate, instead of a small pair of numpy calls per line; the
   result is cached in `_world_points_cache` and invalidated in lockstep with the
   bounding-polygon / node-kdtree / row-lookup caches (i.e. on `rotate` or any node-set
   change, never an elevation-only edit -- same rule as fix 4). `all_points_and_elevation`
   now returns that cached array (read-only for callers, like `get_bounding_polygon()`)
   paired with a fresh `collect("elevation")`, since elevation mutates without a node-set
   change. Output matches the old per-line path to ~1 ULP (2e-16 on unit vectors, from BLAS
   matmul blocking at the larger size; full unit + stepping/plate/elevation-line stress
   suites stay green, including the "preserves spacing exactly" rigid-rotation checks).
   Microbench on a freshly generated 9-plate world (`node_density` 4, ~131 K nodes),
   300 full-plate-cloud sweeps: old per-line style 2.94 s -> vectorized-but-cold 1.24 s ->
   warm cache 0.10 s. In a real step almost every one of the ~24 `all_points_and_elevation`
   calls per plate hits the warm cache.
