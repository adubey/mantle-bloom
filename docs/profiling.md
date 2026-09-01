# Profiling

A running record of where wall-clock time actually goes, measured against a real code
path rather than guessed at. Each section says what was measured, on what commit, with
what inputs, so a later run can be compared like-for-like.

## File > Make Animation (`POST /world/animate`), 60 frames

**Measured:** commit `c8f9990`, 2026-08-31. Python 3.14, single machine, backend `.venv`.

**Inputs** -- the frontend's own defaults (`frontend/src/App.tsx`, `frontend/src/FileModal.tsx`):
`eckert4` projection, `combined` view, 2200x1222 render (1100x611 display x RENDER_SCALE 2),
100 kyr/frame (1 step/frame x default 100 kyr/step), `node_density` = `climate_density` = 4.0,
`fluid_density` = 2.0, identity view rotation.

**Method:** drove `render_image.stream_animation_mp4` directly (same call `/world/animate`
makes) for a full 60-frame run to get end-to-end timings, then a 9-frame run under
`cProfile` plus per-phase wall-time timers for the breakdown. Per-frame cost is stable
enough that the short run is representative. The `combined` view path is
`_render_combined_view` -> `_biome_fields` + `_project_climate_grid` + `_fill_rects`.

### Headline

| | |
| --- | --- |
| 60 frames, wall | **607 s (~10.1 s/frame)** |
| First frame (render only, no step) | 6.5 s |
| Steady state (step + render) | 9.0-12.9 s, mean 10.2 s |
| Split | **~62% rendering, ~38% simulation step** |
| Output | 1.9 MB H.264/MP4 |

Frame time drifts upward over the run (9 s -> 13 s) as plates deform and the river/lake
networks grow denser.

### Rendering -- ~5.9 s/frame

Costs below are totals over the 9-render profiled run (so ~divide by 9 for per-frame).

| Hotspot | Cost | What it is |
| --- | --- | --- |
| `projections._eckert4_theta` | **18.2 s (22% of total)** | 30-iteration Newton solve run over the full **801x1601 = 1.28 M-point** grid, ~8 projection passes per render |
| `render_image._fill_rects` | **14.1 s** | pure-Python `for i in range(1.28M)` slice-assign loop, run **twice** per render (colour buffer + alpha buffer) |
| `render_image._biome_fields` | 11.9 s | rebuilds a `cKDTree` over the whole node cloud and queries 1.28 M points every frame; plus `hydrology.sample_is_ocean` (5.9 s), `climate.compute_climate` (4.9 s), `biomes.smooth_biome_field` (2.6 s) |
| `_project_climate_grid` / `_corner_xy` | (subset of the projection cost above) | 4 extra full-grid projection passes per render purely to size each cell's rectangle |

Everything `_project_climate_grid` returns (`centers`, `half_w`, `half_h`, `scale`,
offsets) is a pure function of (projection, view rotation, grid dimensions, width, height).
**None of those change between animation frames**, yet it is fully recomputed every frame,
and `_biome_grid` itself (no `lru_cache`) rebuilds its lat/lon/xyz arrays each render too.
`_eckert4_theta` depends only on latitude, so it is solving 1.28 M points where ~801
distinct latitudes exist.

### Simulation step -- ~3.5 s/frame

Totals over the 8 profiled steps.

| Hotspot | Cost | Note |
| --- | --- | --- |
| `torque.gather_boundary_force_inputs` | **10.0 s** | called **192x** (~24/step) -- once per plate in `shift`, again per plate in `deform` -- each call builds a fresh `cKDTree` over all neighbour nodes |
| `erosion.apply_erosion` | 12.7 s (1.6 s/step) | `hydrology._compute_basin_spill` (~1 M `heapq.heappop`), `erosion._spread_coastal_leveling`, `lakes.build_lake_hierarchy` |
| `geometry.latlon_to_xyz` | 2.5 s | **237 K calls** via `elevation_lines.world_xyz` (234 K calls) from `PlateWithLines.all_points_and_elevation` -- per-node Python-level conversion |
| `world._advance_fluid_dynamics` | ~0 on the main thread | runs in background threads (540 thread joins, overlapped with the rest of the step) |

### Cheap -- measured, not worth touching

- The PNG-encode then `Image.open`-decode round-trip inside `stream_animation_mp4`:
  **~0.07 s/frame**.
- libx264 encode + muxing: negligible at this frame count and 5/2 fps.

### Highest-value fixes, roughly in order

1. ~~**Cache `_project_climate_grid`'s result across frames** (key on projection + rotation +
   dimensions), and `lru_cache` `_biome_grid`. Saves ~2.5 s/frame -- **~25% off the whole
   animation** -- and speeds every static-camera re-render, not just animation.~~ **Done**
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
4. **Build the neighbour `cKDTree` once per step in `torque`** and share it across every
   plate's `shift`/`deform` instead of rebuilding it ~24x.
5. Vectorize `all_points_and_elevation` / `elevation_lines.world_xyz` to remove the 237 K
   per-node calls.
