# Debugging & Diagnostic Tools

Tracking down a problem such as plate geometry degrading on long runs, craton lost in a merge
or a coastline dithering pixel by pixel usually needs a number the ordinary map views don't
show. This page indexes the tools for that. The details are on the pages under
[`debugging/`](debugging/).

| Page | What it covers |
|---|---|
| [Offline diagnostic CLIs](debugging/offline-cli.md) | `python -m app.*` dumps of a saved world: plate geometry, stranded basins, lake hierarchy, per-phase Hc/Hm budget |
| [Replay logs](debugging/replay-logs.md) | `--replay-log`: the server saves the world before a step whose health metrics go wrong, for an offline rerun |
| [Debug UI](debugging/debug-ui.md) | the `?deb` URL flag, the Event Console and the scripted "Debugging Worlds" scenarios |
| [Debug map views](debugging/map-views.md) | the Map View "Debug >" group: speckle, platesDetail, geomorph, elevReason, overlapAge, nodeAge, cratons, the Plate Inspector |
| [`bin/debug` scripts](#bin-debug-scripts) | about 70 checked-in measurement, probe and replay scripts, indexed below |

For map images of a saved world, use `GET /world/render?view=...`
([api-reference.md](api-reference.md)) or `bin/debug/render_craton_map.py --load` for a plain
plates-and-cratons map. Record the save or seed, elapsed time and view with the image.

## Running the scripts

The scripts in `bin/debug/` run offline: they generate a world or load a `.mbworld` save,
step it in the same process, and print or write JSON. Run them from the repository root with
the worktree's own `backend/.venv/bin/python`. Where a script's usage line shows
`PYTHONPATH=backend`, include it. Each script's module docstring is its full reference; the
pages below summarise it.

- **Baselines.** Scripts with a `--backend` option import the backend it names, so one copy
  of the script measures both an `origin/main` worktree and a branch. All other scripts put
  their own checkout's `backend/` first on `sys.path`, and `PYTHONPATH` does not override
  that. For an `origin/main` baseline, run a copy of the script inside a detached
  `origin/main` worktree. `craton_genuine_area.py` instead takes `--model-root`.
- **Memory.** A replay or probe holds a full world (0.6–2.7 GB). Run one at a time, never
  alongside another replay or a test run.
- **Output.** Write full outputs to scratch space. Commit only compact JSON and findings, in
  `analysis/issue<N>/` (see [AGENTS.md](../AGENTS.md#conventions)). Saves stay out of the
  repository.
- **New scripts.** Give a new script a usage docstring and add it to the index below and to
  its topic page.

<a id="bin-debug-scripts"></a>
## `bin/debug` scripts

### Cratons

[Craton scripts](debugging/scripts-cratons.md): the craton ledger and where its losses come
from.

| Script | Answers |
|---|---|
| [`compare_cratons.py`](debugging/scripts-cratons.md#compare_cratonspy) | multi-seed long runs with cratons on and off |
| [`ledger_history.py`](debugging/scripts-cratons.md#ledger_historypy) | craton and continental-material budgets over a save's whole stats history |
| [`rebaseline_rfcs.py`](debugging/scripts-cratons.md#rebaseline_rfcspy) | RFC 0001–0007 baselines with per-merge and attrition taps |
| [`craton_genuine_area.py`](debugging/scripts-cratons.md#craton_genuine_areapy) | cratonic area inside the genuine continental set, on #403's runs |
| [`rfc0004_gates.py`](debugging/scripts-cratons.md#rfc0004_gatespy) | RFC 0004's formation and craton-share gates from checked-in JSON |
| [`split_craton_rework.py`](debugging/scripts-cratons.md#split_craton_reworkpy) | each craton sink split by the call site that booked it |
| [`probe_craton_clips.py`](debugging/scripts-cratons.md#probe_craton_clipspy) | the step function behind a phase-audit clip; captures lossy merges |
| [`measure_merge_placement.py`](debugging/scripts-cratons.md#measure_merge_placementpy) | why one captured merge lost craton |
| [`attribute_rift_boundary_loss.py`](debugging/scripts-cratons.md#attribute_rift_boundary_losspy) | divergent and stretch loss by whether the split cut crossed craton |
| [`prototype_mincut_split.py`](debugging/scripts-cratons.md#prototype_mincut_splitpy) | whether a minimum-cut split could bend around cratons |

### Craton layout

[Craton layout scripts](debugging/scripts-craton-layout.md): where plates, continents and
cratons sit (RFC 0009).

| Script | Answers |
|---|---|
| [`measure_plate_layout.py`](debugging/scripts-craton-layout.md#measure_plate_layoutpy) | where plate boundaries run relative to continents and cratons, also for Earth |
| [`probe_mixed_plates.py`](debugging/scripts-craton-layout.md#probe_mixed_platespy) | whether mixed continental and oceanic plates behave like single-type ones (gate G2) |
| [`compare_g2.py`](debugging/scripts-craton-layout.md#compare_g2py) | two `probe_mixed_plates.py` results, `origin/main` against a branch |
| [`probe_craton_mosaic.py`](debugging/scripts-craton-layout.md#probe_craton_mosaicpy) | whether the cratons-first layout keeps cratons whole (gate G3) |
| [`measure_belt_marks.py`](debugging/scripts-craton-layout.md#measure_belt_markspy) | belt area, formation audit and step-time cost of belt marks |
| [`render_craton_map.py`](debugging/scripts-craton-layout.md#render_craton_mappy) | a plain map of plates, crust and cratons from a seed or save |
| [`render_layout_sheets.py`](debugging/scripts-craton-layout.md#render_layout_sheetspy) | contact sheets of cratons-first layouts at t = 0 |
| [`render_earth_cratons.py`](debugging/scripts-craton-layout.md#render_earth_cratonspy) | Earth's plates over its cratons, as a reference figure |

### Collisions and plate motion

[Collision and plate-motion scripts](debugging/scripts-collisions.md).

| Script | Answers |
|---|---|
| [`probe_craton_contacts.py`](debugging/scripts-collisions.md#probe_craton_contactspy) | craton jams, closing rates and craton loss at each plate contact (RFC 0008) |
| [`probe_collision_budget.py`](debugging/scripts-collisions.md#probe_collision_budgetpy) | collision overlap, the force that would stall it, and the contact's volume budget |
| [`report_collision_polarity.py`](debugging/scripts-collisions.md#report_collision_polaritypy) | which evidence decides collision polarity, fallback rates, churn and cost |
| [`diagnose_polarity_fallbacks.py`](debugging/scripts-collisions.md#diagnose_polarity_fallbackspy) | contact history behind each polarity fallback |
| [`diagnose_polarity_merge_conflicts.py`](debugging/scripts-collisions.md#diagnose_polarity_merge_conflictspy) | front merges whose records disagree on polarity |
| [`polarity_replay_lib.py`](debugging/scripts-collisions.md#polarity_replay_libpy) | library for the two polarity diagnostics |
| [`audit_force_rules.py`](debugging/scripts-collisions.md#audit_force_rulespy) | which force-proxy rules (rate clamp, merge roll, split gates) bind (RFC 0012) |
| [`probe_boundary_load_normalization.py`](debugging/scripts-collisions.md#probe_boundary_load_normalizationpy) | whether boundary forces per metre are independent of resolution |
| [`probe_rotation_pileup.py`](debugging/scripts-collisions.md#probe_rotation_pileuppy) | how much of fast rotation's land-volume rise is overlap and cell-area bias |
| [`rotation_ledger.py`](debugging/scripts-collisions.md#rotation_ledgerpy) | land area against crust volume as rotation rate changes |

### Rifting and thin-shell diagnostics

[Rifting and thin-shell scripts](debugging/scripts-rifting.md): RFC 0010's connected-source
Gate 2 load, mesh and capacity diagnostics.

| Script | Answers |
|---|---|
| [`audit_gate2_assumption_basis.py`](debugging/scripts-rifting.md#audit_gate2_assumption_basispy) | independent heat-flow, strain-rate and stress-assumption scales |
| [`benchmark_gate2_coupled_shell.py`](debugging/scripts-rifting.md#benchmark_gate2_coupled_shellpy) | manufactured force, torque and resultant checks across a shared two-plate mesh |
| [`calibrate_gate2_capacity.py`](debugging/scripts-rifting.md#calibrate_gate2_capacitypy) | column-capacity scales against published benchmarks and matched loads |
| [`diagnose_gate2_column_envelope.py`](debugging/scripts-rifting.md#diagnose_gate2_column_envelopepy) | depth-integrated normal-fault and ductile lithosphere capacity scenarios |
| [`replay_gate2_assumption_screen.py`](debugging/scripts-rifting.md#replay_gate2_assumption_screenpy) | crossed K0, mantle-rheology and pore-pressure scenarios on fixed loads |
| [`replay_gate2_capacity.py`](debugging/scripts-rifting.md#replay_gate2_capacitypy) | thin-shell loads and capacity alternatives on the pinned Gate 1 events |
| [`replay_gate2_coupled_event.py`](debugging/scripts-rifting.md#replay_gate2_coupled_eventpy) | pair-level stress sensitivity to contact loads and rigid-mode projections |
| [`sensitivity_gate2_capacity.py`](debugging/scripts-rifting.md#sensitivity_gate2_capacitypy) | one-factor capacity sensitivity on captured node loads |
| [`validate_thin_shell_plate_meshes.py`](debugging/scripts-rifting.md#validate_thin_shell_plate_meshespy) | generated-plate mesh coverage, curved-surface and edge-condition checks |

### Crust and mantle budgets

[Crust and mantle budget scripts](debugging/scripts-crust-budgets.md).

| Script | Answers |
|---|---|
| [`attribute_hm_ratchet.py`](debugging/scripts-crust-budgets.md#attribute_hm_ratchetpy) | which phases push Hm toward its cap |
| [`compare_crust_transfer.py`](debugging/scripts-crust-budgets.md#compare_crust_transferpy) | `attribute_hm_ratchet.py` replays side by side |
| [`summarize_issue317_suture_pairs.py`](debugging/scripts-crust-budgets.md#summarize_issue317_suture_pairspy) | #317's suture Hm totals by donor and neighbour pair |
| [`calibrate_collision_sinks.py`](debugging/scripts-crust-budgets.md#calibrate_collision_sinkspy) | continental crust lost to collisions, per step and compounded, with constant overrides |
| [`attribute_suture_delamination.py`](debugging/scripts-crust-budgets.md#attribute_suture_delaminationpy) | the call site behind each delamination and no-outlet booking |
| [`attribute_topology_removed.py`](debugging/scripts-crust-budgets.md#attribute_topology_removedpy) | what topology clean-up drops as `topology_removed_m3` |
| [`probe_quad_retreat_continental.py`](debugging/scripts-crust-budgets.md#probe_quad_retreat_continentalpy) | where continental crust removed by boundary retreat goes |
| [`measure_cap_overflow.py`](debugging/scripts-crust-budgets.md#measure_cap_overflowpy) | depositional Hc-cap overflow and where it ends up |
| [`measure_collision_denudation.py`](debugging/scripts-crust-budgets.md#measure_collision_denudationpy) | Hc-cap saturation and erosion of collision belts |
| [`measure_orogenic_relief.py`](debugging/scripts-crust-budgets.md#measure_orogenic_reliefpy) | whether continents keep area and volume through collisions |
| [`measure_erosion_budget.py`](debugging/scripts-crust-budgets.md#measure_erosion_budgetpy) | how closely one erosion call conserves material |

### Elevation and relief

[Elevation and relief scripts](debugging/scripts-elevation.md).

| Script | Answers |
|---|---|
| [`measure_initial_relief.py`](debugging/scripts-elevation.md#measure_initial_reliefpy) | how Earth-like a world's hypsometry, bathymetry and relief are |
| [`diagnose_elevation_saturation.py`](debugging/scripts-elevation.md#diagnose_elevation_saturationpy) | how much land sits at `MAX_ELEVATION_M`, and whether Hc is capped there too |
| [`drill_down_isostasy.py`](debugging/scripts-elevation.md#drill_down_isostasypy) | whether capped land is backed by its own column |
| [`measure_transform_far_field_debt.py`](debugging/scripts-elevation.md#measure_transform_far_field_debtpy) | isostatic debt over a long run |
| [`drill_down_volcanism.py`](debugging/scripts-elevation.md#drill_down_volcanismpy) | why more volcanism gave less land (#173) |
| [`measure_volcano_node_elevation.py`](debugging/scripts-elevation.md#measure_volcano_node_elevationpy) | the elevation of nodes that have erupted |

### Water and ice

[Water and ice scripts](debugging/scripts-water-ice.md).

| Script | Answers |
|---|---|
| [`measure_river_endings.py`](debugging/scripts-water-ice.md#measure_river_endingspy) | rivers that reach neither sea nor lake |
| [`measure_discharge_calibration.py`](debugging/scripts-water-ice.md#measure_discharge_calibrationpy) | river and glacier calibration on one save in two trees |
| [`measure_lake_silt.py`](debugging/scripts-water-ice.md#measure_lake_siltpy) | lake area, volume and silt |
| [`measure_breaching.py`](debugging/scripts-water-ice.md#measure_breachingpy) | depression breaching on and off |
| [`measure_ice_sheet.py`](debugging/scripts-water-ice.md#measure_ice_sheetpy) | glacier extent, volume, erosion and the ice ledger |
| [`measure_seasonal_melt.py`](debugging/scripts-water-ice.md#measure_seasonal_meltpy) | glaciated land by annual-mean temperature |
| [`measure_ice_loading.py`](debugging/scripts-water-ice.md#measure_ice_loadingpy) | isostatic ice-load deflection |
| [`measure_meltwater_faults.py`](debugging/scripts-water-ice.md#measure_meltwater_faultspy) | earthquakes triggered by ice unloading |
| [`measure_climate_fields.py`](debugging/scripts-water-ice.md#measure_climate_fieldspy) | temperature, humidity and precipitation by latitude band in two trees |

### Surface, saves, sweeps and profiling

[Surface, save, sweep and profiling scripts](debugging/scripts-infrastructure.md).

| Script | Answers |
|---|---|
| [`surface_parity.py`](debugging/scripts-infrastructure.md#surface_paritypy) | plate-surface audit runs and quality gates |
| [`measure_quad_area_churn.py`](debugging/scripts-infrastructure.md#measure_quad_area_churnpy) | per-step churn in represented area and the sea-level jitter it drives |
| [`convert_legacy_saves.py`](debugging/scripts-infrastructure.md#convert_legacy_savespy) | converting line-backed saves to quads, with audits |
| [`replay_logged_event.py`](debugging/scripts-infrastructure.md#replay_logged_eventpy) | whether a replay-log record reproduces |
| [`run_sweep.py`](debugging/scripts-infrastructure.md#run_sweeppy) | the multi-parameter sensitivity sweep |
| [`aggregate_sweep.py`](debugging/scripts-infrastructure.md#aggregate_sweeppy) | sweep rows summarised per parameter and checkpoint |
| [`profile_animation_60m.py`](debugging/scripts-infrastructure.md#profile_animation_60mpy) | per-frame timings of a 60 Myr animation |
| [`seed_run_lib.py`](debugging/scripts-infrastructure.md#seed_run_libpy) | library: generate-and-step runs with `--backend` |
| [`replay_lib.py`](debugging/scripts-infrastructure.md#replay_libpy) | library: save replays with `--backend` |
| [`sweep_lib.py`](debugging/scripts-infrastructure.md#sweep_libpy) | library: sweep parameters, seeds and job runner |
