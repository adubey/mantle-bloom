# Line vs quad parity report

**Verdict: PASS** — configs issue147; seeds 804913535.

Reproduce:

    bin/debug/surface_parity.py run --preset issue147 --seed 804913535 --surface lines --from-world /Users/adubey/Downloads/mantle-bloom-seed804913535-352400000y.mbworld --step-years 100000 --checkpoints-myr 0.5,1 --audit-every 1 --render --out ../analysis/issue249-campaign/repro-lines

Tolerances and gate definitions: docs/surface-parity.md. `info` rows are line-baseline
findings; they never gate the quad surface.

## Hard invariants

| gate | status | worst case |
|---|---|---|
| H4:fields | pass | violations=0, audits=11, surface=lines, seed=804913535 |
| H5:frames | pass | violations=0, audits=11, surface=lines, seed=804913535 |
| H6:derived_caches | pass | violations=0, audits=11, surface=lines, seed=804913535 |
| H7:revisions | pass | violations=0, audits=11, surface=lines, seed=804913535 |
| H8:load_round_trip | pass | checks=3, surface=lines, seed=804913535 |
| H9:load_continuation | pass | surface=lines, seed=804913535 |
| H11:field_caps | info (line baseline: 1 finding(s)) | violations=205, surface=lines, seed=804913535; first: {'detail': 'mantle_lithosphere_thickness_m: 1 values outside [2000, 240000] (range 1623.76..39276.7)', 'kind': 'field_bounds', 'plate_id': 1, 'step': 0, 'surface': 'lines'} |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=804913535 |

## Checkpoints, seed 804913535

| seed 804913535 | 352.4 Myr lines | 352.9 Myr lines | 353.4 Myr lines |
|---|---|---|---|
| plates | 35 | 35 | 35 |
| nodes | 158589 | 158593 | 158756 |
| land fraction | 0.3892 | 0.3924 | 0.3983 |
| sea level (m) | -826.2 | -831.5 | -837.7 |
| Hc volume (km³) | 1.31e+10 | 1.31e+10 | 1.3e+10 |
| continental Hc (km³) | 9.08e+09 | 9.07e+09 | 9.06e+09 |
| elevation p05/p50/p95 (m) | -5340/-2491/6064 | -5333/-2462/6081 | -5323/-2388/6115 |
| uncovered | 0.0232 | 0.02341 | 0.02261 |
| multiply covered | 0.004735 | 0.00392 | 0.004505 |
| void | 0.00014 | 0.0002 | 0.00025 |
| nodes inside other plate | 0.01379 | 0.01257 | 0.01383 |
| stacked | 0.275 | 0.2761 | 0.2767 |
| anisotropic | 0.04662 | 0.04726 | 0.04814 |
| row alignment | 0.8406 | 0.8422 | 0.8303 |
| thin | 0.05789 | 0.05775 | 0.05702 |
| one-node lines | 0.293 | 0.3036 | 0.3092 |
| quad one-cell-thin | — | — | — |
| air temperature mean (°C) | -7.643 | -7.687 | -7.675 |
| precipitation mean (mm) | 1021 | 1027 | 1027 |
| river fraction of land | 0.1134 | 0.1131 | 0.1132 |
| HEALPix same node | 0.5346 | 0.5358 | 0.534 |

## Performance (mean seconds per step)

| s/step | seed 804913535 lines | #147 profile |
|---|---|---|
| step_total | 7.167 | 10.26 |
| deform_topology | 2.744 | — |
| climate_erosion_hydrology | 1.009 | 1.073 |
| deform | 2.473 | 2.58 |
| faults | 0.711 | 1.157 |
| fluid_dynamics | 3.29e-06 | — |
| gap_fill | 0.03358 | — |
| magma_transport | 1.688 | — |
| overlap_tracking | 0.05423 | — |
| record_stats | 0.009285 | — |
| resource_formation | 0.04062 | — |
| sea_level | 0.2404 | 0.2404 |
| shift | 0.6346 | 0.8096 |
| stranded_basins | 0.002115 | — |
| topology | 0.1828 | 0.2014 |
| volcanism | 0.08588 | 0.0845 |

Index builds per step (mean):

- seed 804913535 lines: line_row_lookup 89.5, plate_node_kdtree 96.0, plate_outline 97.6, plate_outline_kdtree 95.8, world_node_kdtree 1.0

The #147 profile ran under cProfile, so its absolute times are inflated.

