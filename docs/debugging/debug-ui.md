# Debug UI

## The `?deb` URL flag

Load the frontend with `?deb` anywhere in the URL query string (e.g.
`http://localhost:5173/?deb`) and a cluster of developer-only UI appears that is otherwise
**absent entirely, not just disabled** -- an ordinarily-generated world's sidebar and dialogs
stay uncluttered by default. It's a one-time flag read from `window.location.search` at module
load (`frontend/src/App.tsx`'s `DEBUG_UI` constant), not something that changes while the app
is open, and it works identically against the packaged prod server and the dev server -- it
only ever looks at the URL the browser loaded, never which backend is serving it.

`?deb` gates three things:

1. The Map View select's **"Debug >"** option group ([Debug map views](map-views.md)).
2. The **Event Console** side panel ([below](#event-console)).
3. The Generate World dialog's **"Debugging Worlds"** tab ([below](#debugging-worlds-tab----scripted-plate-scenarios)).

Everything on this page and on [Debug map views](map-views.md) requires `?deb` in the URL unless noted otherwise.

## Event Console

A side panel rendered only when `?deb` is set (`frontend/src/App.tsx`'s `EventConsole`).
Without `?deb` the underlying data is still reachable through the API -- `GET /world/summary`'s
`events` field -- only the in-app viewer is hidden. It also lists the records the server's
replay log wrote (the summary's `replay_logs`, see [Replay logs](replay-logs.md)), merged in
order of elapsed time.

<a id="event-log"></a>
### Event log

The `events` list on `GET /world/summary`, shown in the Event Console, logs lake
formation/splits and other discrete events.

**Lake-churn aggregation -- `lakes.summarize_lake_events`.** On a long run over a dithering
low-relief coast the lake solver produces hundreds of near-sea-level transient merge/split
transitions per My -- one pair per puddle per step (see
[GitHub issue #122](https://github.com/adubey/mantle-bloom/issues/122)). Left raw, these flood
the console and bury real basin/tectonic events.

`lakes.step_lakes` returns structured `lakes.LakeEvent`s (`kind` / `node_count` /
`elevation_m` / `basin_count`, plus a `.message` property with the same wording as before)
instead of pre-formatted strings. [`erosion.py`](../../backend/app/hydroclimate/erosion.py) runs a step's
events through `lakes.summarize_lake_events(events, world.sea_level_m)` before logging:

- A transition whose water surface is more than `lakes.NEAR_SEA_LEVEL_EVENT_BAND_M`
  (**15 m**) from the current sea level is a genuine basin event -- logged individually,
  unchanged.
- Transitions **within** that band are the coastal-pond churn. A lone one still logs
  verbatim; **two or more in one step collapse to a single line** --
  `"38 transient coastal ponds churned near sea level this step (22 merged, 16 split)."`

So a persistent deep endorheic basin (e.g. a real ~435-node lake oscillating near -1770 m)
stays visible in the console while the checkerboard shelf contributes at most one aggregate
line per step. The band is measured against `world.sea_level_m`, so it tracks a sea-level
control change. This aggregation always runs (it isn't gated by `?deb`) -- only the console
that displays the result is. Tests:
[`unit_tests/test_lakes.py`](../../backend/unit_tests/test_lakes.py)
(`test_summarize_lake_events_*`).

## "Debugging Worlds" tab -- scripted plate scenarios

The Generate World dialog's fourth tab (`POST /world/generate_debug`,
`backend/app/debug_worlds.py`) builds tiny, low-resolution (`node_density=0.5`, the coarsest
real choice), hand-placed plate configurations for fast iteration on the gap-filling problem,
independent of any real save's history. `GET /world/debug_scenarios` lists the current
options; picking one and pressing Generate replaces the current world exactly like the
Random/Human-made/Premade tabs, with `debug_diagnostics` already on (so the phase budget
records -- see [`python -m app.phase_budget_diagnostics`](offline-cli.md#python--m-appphase_budget_diagnostics----per-phase-hchm-budget)).

Each plate's motion is **pinned**, not torque-driven: `World.pinned_omegas` (`plate_id ->` a
fixed angular velocity) fixes that plate's degrees of freedom in the coupled solve every step.
Other plates still solve against the pinned plate's drag and contact loads. The batch then
uses `torque.apply_omega_and_rotate` to rotate all plates and report displacement. This is what
makes a scenario reliable to iterate on: "these two plates
diverge" stays true step after step, rather than however the force balance against an
otherwise-irrelevant `mantle_centers` field happens to settle it. `pinned_omegas` is empty for
every ordinarily-generated world -- this is purely a debug-world mechanism.

A scenario is a set of seed lat/lons, each plate's crust type, and a list of pairwise
`(plate_i, plate_j, "divergent" | "convergent")` boundary relationships;
`debug_worlds._omegas_from_relationships` converts that into each plate's own single rotation
(the sum of its own boundaries' tangential "push" contributions, converted to a rotation axis
via `omega = cross(seed_xyz, tangent_direction)`) so a plate sitting at a junction of several
boundaries gets one rotation reflecting all of them at once, the same way a real plate's
motion is a net result of every boundary force acting on it -- just picked directly here
instead of emerging from `torque.py`'s force balance.

Current scenarios (`debug_worlds.DEBUG_SCENARIOS`):

| Scenario | Layout | Purpose |
|---|---|---|
| `two_plate_divergent` | 2 plates, straight rift | Baseline: ordinary boundary advance and rift opening should close this alone, never leaving a real gap. |
| `two_plate_convergent` | 2 plates, closing boundary | Subduction/collision baseline. |
| `triple_junction_mixed` | 3 plates, one junction, 2 divergent legs + 1 convergent | Three plates receding from a shared point, so boundary growth has to keep pace with three-way divergence (confirmed as a gap source on a real save, seed 430031492). Start here. |
| `four_plate_grid` | 4 plates, 2x2, all edges divergent | Four simultaneous triple-junction-like corners at once, around one shared center point. |
| `five_plate_irregular` | 5 plates, irregular ring, mixed relationships | Closest single scenario to a real save's messiness while staying small enough to iterate on quickly. |

Workflow for chasing a junction gap: enable `?deb`, generate `triple_junction_mixed`, step it
a handful of times (Play/Stop, small "Years per step" for fine-grained observation), and watch
the "Added/removed points" and "Plate overlap age" views together.

These tiny scripted scenarios don't naturally exercise
`gaps.fill_gaps_by_growing_neighbours` (every plate starts adjacent to every gap that can open
here) -- checking the whole-sphere sweep needs a real or synthetically-vacated-plate world
instead, see `unit_tests/test_gaps.py`'s `fill_gaps_by_growing_neighbours` tests for the shape
of that setup.
