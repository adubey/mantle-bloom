# #801: proposed traces

Two layers. **User-visible traces** are the labels the run overlay shows, each with
Wikipedia and mantle-bloom doc links (`frontend/src/phaseInfo.ts`). **Profiling traces** are
finer spans recorded for developers and nested under a user-visible one, so every profiling
span can be rolled up into the bar the user sees. This replaces the issue's requirement that
the two match one to one. Shares are of one step or one generation on the young world (Y) and
the 67 Myr world (M); see [README.md](README.md). On the hosted server everything is slower, so
the shares matter more than the seconds.

Rule used: a phase is user-visible if it is at least ~5% of its bar in either world, or is a
concept worth a link. Everything else stays a profiling trace under the nearest user-visible
parent.

## Findings that change the existing labels

- `"Resolving glaciers"` is yielded before the wind solver and the whole of erosion, so it
  is on screen for 2.6 s (Y) to 6 s (M) of work that is mostly rivers, lakes and climate.
- `"Moving plates"` and `"Reshaping plates"` are yielded after each unit finishes. The motion
  solve, collision-polarity pass and fault update run with the previous label still showing.
- Labels should be emitted at the start of the work they name.

## User-visible traces (hierarchical)

The overlay shows one bar per top-level group and the current second-level label as text under
it, with its links. A phase gets its own second-level label if it is at least ~5% of its group
or is a concept worth a link, so sub-1% phases are leaves that cost no visual churn. Shares are
young world (Y) / 67 Myr world (M) of one step or generation.

### Building world

- **Laying out the planet**
  - Laying out continents and mountain belts (39% Y)
  - Dividing the sphere into plates (27% Y)
  - Placing sutures and cratons (existing label and links)
  - Building pseudo-history (32% when enabled; off by default)
- **Seeding the surface**
  - Seeding rivers and glaciers (8-10% Y)
  - Seeding climate and oceans (~10% Y; climate, bathymetry, rock ages)

### Loading a world

Headline: "Loading, rebuilding render cache & resolving climate".

- Reading the save (~1.6 s on the dev machine)
- Rebuilding render cache (k-d trees; ~1.2 s on the first render)
- Resolving climate (the first step also recomputes missing climate, hydrology and CFD state)

### Stepping

- **Moving plates** (51% of a young step, 35% with CFD wind, 59% of a 67 Myr step)
  - Solving mantle flow and plate motion (9% Y / 7% M of a step)
  - Detecting collisions (5% / 8%)
  - Reshaping plates (7% / 10%)
  - Updating faults (16% / 13%)
  - Filling seafloor gaps (8% / 10%)
  - Merging and splitting plates (finding rifts) (<1% / 2%), leaf
  - Moving magma (1% / 4%), leaf
  - Raising mantle plumes (mantle blooms) (<1%), leaf
- **Hydrometeorology (climate, water & ice)** (34% young, 54% with CFD wind, 34% mature)
  - Solving wind and ocean currents (31% with CFD wind, absent otherwise)
  - Resolving climate (8% / 4%)
  - Resolving rivers and lakes (20% / 31%)
  - Marching glaciers (7% / 10%)
- **Reshaping the surface** (13% young, 9% CFD, 5% mature)
  - Eroding and depositing sediment (erosion's own work, 4% / 2%)
  - Raising volcanoes, leaf
  - Forming natural resources, leaf
  - Balancing sea level (1-3%), leaf

Group shares are measured (the by-group rollup in `profile_phases.py`), and the groups follow the order the code
runs today. Rivers, lakes, glaciers and climate all run inside `apply_erosion`, so the
Hydrometeorology group is not contiguous with the others in time; the trace API needs to let a
span name its group explicitly rather than inferring it from call order.

## Profiling traces (nested under the above)

Measured with the expanded hook table (`hooks-v2-*.json`; per-step seconds, own time summed
over every parent a function runs under). Y = young world, M = 67 Myr world, both diagnostic
wind unless noted. Untraced time is 1.3% (Y), 0.9% (CFD) and 0.7% (M) of a step.

| Under | Profiling trace | Function | Y | M |
|---|---|---|---|---|
| Wind and ocean (CFD only) | stepping the CFD solver | `atmosphere_cfd.step_atmosphere_cfd` | 2.66 (CFD) | -- |
| Mantle flow and motion | coupled plate solve (own) | `torque.solve_plates_quasistatic` | 0.19 | 0.42 |
| | gathering boundary forces (3 callers) | `torque.gather_boundary_force_inputs` | 0.40 | 0.64 |
| | resolving blocked contacts | `torque._blocked_contact_response` | 0.15 | 0.30 |
| Detecting collisions | observe contacts (own) | `collision_polarity.observe_contacts` | 0.14 | 0.19 |
| Moving plates (shared) | outlining plate territory (4 callers) | `PlateWithSparseQuadPatch.outline_world` | not hooked | 2.92 |
| Reshaping plates | deform (own) | `PlateWithSparseQuadPatch.deform` | 0.28 | 0.30 |
| | growing plate frontiers (2 callers) | `quad_tectonics.grow_frontier` | 0.33 | 0.34+ |
| Updating faults | shear along faults | `faults._apply_plate_fault_shear` | 0.36 | 0.66 |
| | fault relief | `faults._apply_plate_fault_relief` | 0.30 | 0.66 |
| | boundary-fault generation | `faults.generate_boundary_faults` | 0.21 | 0.55 |
| | earthquakes | `faults._generate_earthquakes` | 0.12 | 0.27 |
| | ice-unloading earthquakes | `faults.trigger_unloading_earthquakes` | 0.11 | 0.28 |
| Filling gaps | finding gap points | `gaps._find_gap_points` | 0.16+0.05 | -- |
| | per-plate gap fill (own) | `quad_tectonics.fill_gap` | 0.03 | 0.00 |
| Moving magma | | `magma_transport.run_magma_transport` | 0.07 | 0.59 |
| Resolving rivers and lakes | lakes | `lakes.solve_lakes_prognostic` | 0.11 | 1.68 |
| | ordering the flow network (4 callers) | `hydrology._routing_order` | 0.33 | not hooked per caller |
| | routing flow downstream | `hydrology.route_downstream` | 0.39 | 0.43 |
| Marching glaciers | shaping ice sheets (250-700 calls) | `ice_sheet._shape` | 0.23 | 1.41 |
| Resolving climate | ocean currents | `climate.compute_ocean_currents` | 0.14 | -- |
| | humidity | `climate.compute_humidity` | 0.05 | -- |
| Balancing sea level | connected-ocean solve | `eustasy._solve_sea_level_connected` | 0.17 | -- |
| Erosion (own work) | sediment passes: marine, mass wasting, others | `erosion._spread_*`, `erosion._route_*` | 0.09 + 0.07 | -- |

Building world (synthetic history, 16.8 s): raising orogenic belts 2.1 s and shaping continents
2.0 s (`continent_layout`), 100 Voronoi partition attempts 3.0 s (`plate_partition`), pseudo-history
4.6 s own plus a coupled plate solve inside it.

Render cache (young world, cold): nearest-node lookup per pixel 0.66 s (0.29 s warm), climate
resample 0.49 s (absent warm), biome classification 1.55 s own on every render.

### What the finer hooks showed

- **The mature world's cost is spread differently.** Lakes (1.70 s) and ice shaping (1.43 s) are
  ~21% of a step, and plate outlining (below) is another 20%.
- **Shared functions hide their cost.** Boundary-force gathering runs under three phases and
  frontier growth under two; the per-label total in the profiler output is the number to use.
- **Erosion's own sediment passes are small** (about 0.3 s total on both worlds), so splitting
  "Eroding and depositing sediment" further would buy little.
- **CFD:** the solver step is the whole 2.7 s; forcing refresh is under 1.5% of a step.
- **Plate outlines are the largest single label on the mature world.** `hooks-v3` shows
  "Outlining plate territory" (`PlateWithSparseQuadPatch.outline_world`) at 2.92 s of a 14.9 s
  step (19.6%), recomputed ~29 times per step (about 10 plates, 3 times each). It is reached through
  `get_neighbours` from four phases: per-plate gap fill (1.00 s), collision polarity (0.76 s),
  plate deformation (0.72 s) and the plate-motion solve (0.2 s). A cProfile drill-down found the
  same thing, so it is not a profiler artifact. `get_bounding_polygon` is cached (1,608 calls,
  88 recomputes over 3 steps), but each plate's cache is dropped whenever its geometry changes,
  and each recompute costs ~100 ms (2.9 s over ~29 recomputes); tracked in #808. This explains most of the "own time" that collision polarity,
  deformation and gap fill showed at ~1 s each. It is a profiling-only trace, not a user-visible one.
- Still not split: the per-caller breakdown of the flow-network ordering on the mature world, and
  the remaining own time of `deform` (0.3 s) and `fill_gap`.

## Implementation

Built in `backend/app/trace.py` (see [architecture.md](../../docs/architecture.md#trace-spans)).

- **Spans.** `@trace.traced(label, group, user=True)` on 21 coarse functions, plus one
  `trace.span` around the plate-deform loop (a per-plate span would credit its share of the bar
  after the first plate), gives the 22 user-visible labels, with the group stated at the span so Hydrometeorology is not inferred from
  call order. A span is a no-op unless a `Recorder` is active.
- **Delivery.** The step and generate endpoints run the progress generator on a worker thread
  (`trace.run_traced`) so a label from inside erosion reaches the client when its work starts. A
  playing world publishes the phase on its runner; `GET /worlds/{id}/run` returns it as `phase`
  and the app reads it every 1.5 s. The generators, their tests, the replay log and
  `_drain_progress` are unchanged.
- **Fractions.** The share of each phase's typical time (measured on the young world) already
  finished, never decreasing, below 1 until the stream ends. On a mature world lakes and ice
  dominate, so the bar sits still through them.
- **Profiling layer.** `bin/debug/profile_phases.py` wraps the finer functions in plain spans at
  runtime, so production carries none of them, and prints the same tree with a by-group rollup.
- **Frontend.** The overlay shows the group above the phase and learn links for each new label;
  a test checks every user-visible label has links.

Not done: the Loading flow (`handleWorldReplaced`) is not streamed, so its overlay shows the
headline and no phases; the `Hydrometeorology` bar sits still through long lake solves on mature
worlds; weights were measured on one young seed.
