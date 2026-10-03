# Line vs quad parity report

**Verdict: WARN** — configs issue147; seeds 804913535.

Reproduce:


Tolerances and gate definitions: docs/surface-parity.md. `info` rows are line-baseline
findings; they never gate the quad surface.

## Hard invariants

| gate | status | worst case |
|---|---|---|
| H1:quad_topology | pass | violations=0, audits=11, surface=quad, seed=804913535 |
| H2:quad_neighbours | pass | violations=0, audits=11, surface=quad, seed=804913535 |
| H3:quad_folded | pass | violations=0, audits=11, surface=quad, seed=804913535 |
| H4:fields | pass | violations=0, audits=11, surface=lines, seed=804913535 |
| H5:frames | pass | violations=0, audits=11, surface=lines, seed=804913535 |
| H6:derived_caches | pass | violations=0, audits=11, surface=lines, seed=804913535 |
| H7:revisions | pass | violations=0, audits=11, surface=lines, seed=804913535 |
| H8:load_round_trip | pass | checks=3, surface=lines, seed=804913535 |
| H9:load_continuation | pass | surface=lines, seed=804913535 |
| H10:quad_no_stacked_nodes | pass | surface=quad, seed=804913535 |
| H11:field_caps | info (line baseline: 1 finding(s)) | violations=19, surface=lines, seed=804913535; first: {'detail': 'mantle_lithosphere_thickness_m: 1 values outside [2000, 240000] (range 1623.76..39276.7)', 'kind': 'field_bounds', 'plate_id': 1, 'step': 0, 'surface': 'lines'} |

## Coverage and overlap

| gate | status | worst case |
|---|---|---|
| C1:uncovered | pass | age_myr=352.4, quad=0.00829, lines=0.0232, fail_above=0.02821, seed=804913535 |
| C2:void | pass | age_myr=352.4, quad=0.00012, lines=0.00014, fail_above=0.0005, seed=804913535 |
| C3:multiply_covered | pass | age_myr=352.4, quad=0.00713, warn_above=0.015, fail_above=0.02, surface=quad, seed=804913535 |
| C4:nodes_inside_other_plate | pass | age_myr=352.4, quad=0.01165, lines=0.01379, fail_above=0.01879, seed=804913535 |

## Conservation (actual node/cell areas)

| gate | status | worst case |
|---|---|---|
| K1:hc_volume_km3_drift | pass | age_myr=352.4, quad=0, lines=0, gap=0, warn_above=0.05, fail_above=0.1, seed=804913535 |
| K2:continental_hc_volume_km3_drift | pass | age_myr=352.4, quad=0, lines=0, gap=0, warn_above=0.05, fail_above=0.1, seed=804913535 |
| K3:quad_area_accounting | pass | age_myr=352.4, quad=0.9991, expected=0.9988, warn_above=0.01, seed=804913535 |
| K4:hm_volume_km3_drift | pass | age_myr=352.4, quad=0, lines=0, gap=0, warn_above=0.05, fail_above=0.1, seed=804913535 |

## Mesh quality (quad)

| gate | status | worst case |
|---|---|---|
| M1:anisotropic_fraction | pass | age_myr=352.4, quad=0.001235, lines=0.04662, warn_above=0.01, seed=804913535 |
| M2:row_alignment | pass | age_myr=352.4, quad=0.1901, lines=0.8406, warn_above=0.3, seed=804913535 |
| M3:thin_fraction | pass | age_myr=352.4, quad=0.001389, lines=0.05789, warn_above=0.01, seed=804913535 |
| M4:aspect_gt_4 | pass | age_myr=352.4, quad=0, warn_above=0.001, seed=804913535 |
| M5:skew_gt_45 | pass | age_myr=352.4, quad=0, warn_above=0.001, seed=804913535 |
| M6:fragments | warn | age_myr=352.4, quad=39, lines=0, seed=804913535 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=804913535 |
| S2:stability_air_temperature_mean_c | pass | quad=0.01764, lines=0.01924, warn_above=0.08848, seed=804913535 |
| S2:stability_land_fraction | pass | quad=0.000214, lines=0.000109, warn_above=0.002218, seed=804913535 |
| S2:stability_ocean_temperature_mean_c | pass | quad=0.001828, lines=0.002744, warn_above=0.05549, seed=804913535 |
| S2:stability_precipitation_mean_mm | pass | quad=9.607, lines=9.819, warn_above=24.64, seed=804913535 |
| S2:stability_sea_level_m | pass | quad=0.6746, lines=1.671, warn_above=4.343, seed=804913535 |

## Derived indexes: HEALPix vs KD-tree

| gate | status | worst case |
|---|---|---|
| X1:healpix_all | pass | age_myr=352.4, samples=200000, quad_same_node=0.5544, lines_same_node=0.5346, quad_distance_p95=1.179, lines_distance_p95=1.174, seed=804913535 |
| X1:healpix_antimeridian | pass | age_myr=352.4, samples=3333, quad_same_node=0.5464, lines_same_node=0.4599, quad_distance_p95=1.199, lines_distance_p95=1.18, seed=804913535 |
| X1:healpix_coast | pass | age_myr=352.4, samples=23978, quad_same_node=0.5527, lines_same_node=0.5465, quad_distance_p95=1.185, lines_distance_p95=1.178, seed=804913535 |
| X1:healpix_hole | pass | age_myr=352.4, samples=1413, quad_same_node=0.5216, lines_same_node=0.4092, quad_distance_p95=1.306, lines_distance_p95=1.311, seed=804913535 |
| X1:healpix_plate_boundary | pass | age_myr=352.4, samples=17107, quad_same_node=0.5266, lines_same_node=0.5088, quad_distance_p95=1.283, lines_distance_p95=1.309, seed=804913535 |
| X1:healpix_pole | pass | age_myr=352.4, samples=3038, quad_same_node=0.5402, lines_same_node=0.4908, quad_distance_p95=1.231, lines_distance_p95=1.23, seed=804913535 |

## Performance

| gate | status | worst case |
|---|---|---|
| R1:deform_topology_s_per_step | pass | quad=1.737, lines=2.953, ratio=0.5884, warn_above=1.25, seed=804913535 |
| R2:step_total_s_per_step | pass | quad=5.824, lines=7.958, ratio=0.7318, warn_above=1.25, seed=804913535 |

## Ensemble parity

| gate | status | worst case |
|---|---|---|
| P1:continental_hc_volume_km3 | insufficient | age_myr=352.4, seeds=1, lines_mean=9.08e+09, quad_mean=9.08e+09, note=needs >= 3 seeds |
| P1:elevation_p05 | insufficient | age_myr=352.4, seeds=1, lines_mean=-5340, quad_mean=-5426, note=needs >= 3 seeds |
| P1:elevation_p50 | insufficient | age_myr=352.4, seeds=1, lines_mean=-2491, quad_mean=-3420, note=needs >= 3 seeds |
| P1:elevation_p95 | insufficient | age_myr=352.4, seeds=1, lines_mean=6064, quad_mean=5820, note=needs >= 3 seeds |
| P1:hc_volume_km3 | insufficient | age_myr=352.4, seeds=1, lines_mean=1.31e+10, quad_mean=1.41e+10, note=needs >= 3 seeds |
| P1:land_fraction | insufficient | age_myr=352.4, seeds=1, lines_mean=0.3892, quad_mean=0.3313, note=needs >= 3 seeds |
| P1:plates | insufficient | age_myr=352.4, seeds=1, lines_mean=35, quad_mean=35, note=needs >= 3 seeds |
| P1:sea_level_m | insufficient | age_myr=352.4, seeds=1, lines_mean=-826.2, quad_mean=-826.2, note=needs >= 3 seeds |

## Checkpoints, seed 804913535

| seed 804913535 | 352.4 Myr lines | 352.4 Myr quad | 352.9 Myr lines | 352.9 Myr quad | 353.4 Myr lines | 353.4 Myr quad |
|---|---|---|---|---|---|---|
| plates | 35 | 35 | 35 | 35 | 35 | 36 |
| nodes | 158589 | 129602 | 158591 | 129570 | 158761 | 129506 |
| land fraction | 0.3892 | 0.3313 | 0.3924 | 0.3288 | 0.3982 | 0.3291 |
| sea level (m) | -826.2 | -826.2 | -831.5 | -810.5 | -838 | -815.9 |
| Hc volume (km³) | 1.31e+10 | 1.41e+10 | 1.31e+10 | 1.41e+10 | 1.3e+10 | 1.41e+10 |
| continental Hc (km³) | 9.08e+09 | 9.08e+09 | 9.07e+09 | 9.08e+09 | 9.06e+09 | 9.08e+09 |
| elevation p05/p50/p95 (m) | -5340/-2491/6064 | -5426/-3420/5820 | -5333/-2462/6081 | -5417/-3410/5858 | -5324/-2391/6114 | -5417/-3411/5983 |
| uncovered | 0.0232 | 0.00829 | 0.02343 | 0.00591 | 0.0226 | 0.00612 |
| multiply covered | 0.004735 | 0.00713 | 0.00392 | 0.004805 | 0.00448 | 0.00465 |
| void | 0.00014 | 0.00012 | 0.0002 | 0 | 0.00025 | 0 |
| nodes inside other plate | 0.01379 | 0.01165 | 0.01257 | 0.002717 | 0.01383 | 0.002417 |
| stacked | 0.275 | 0 | 0.2761 | 0 | 0.2767 | 0 |
| anisotropic | 0.04662 | 0.001235 | 0.04726 | 0.000641 | 0.04814 | 0.000672 |
| row alignment | 0.8406 | 0.1901 | 0.8422 | 0.2855 | 0.8303 | 0.3068 |
| thin | 0.05789 | 0.001389 | 0.05774 | 0.000841 | 0.05703 | 0.00088 |
| one-node lines | 0.293 | — | 0.3036 | — | 0.3096 | — |
| quad one-cell-thin | — | 0.002971 | — | 0.002987 | — | 0.003228 |
| air temperature mean (°C) | -7.643 | -7.643 | -7.687 | -8.07 | -7.674 | -8.052 |
| precipitation mean (mm) | 1021 | 1021 | 1027 | 1024 | 1027 | 1025 |
| river fraction of land | 0.1134 | — | 0.1131 | 0.1243 | 0.1132 | 0.1232 |
| HEALPix same node | 0.5346 | 0.5544 | 0.5358 | 0.5557 | 0.534 | 0.5567 |

## Performance (mean seconds per step)

| s/step | seed 804913535 lines | seed 804913535 quad | #147 profile |
|---|---|---|---|
| step_total | 7.958 | 5.824 | 10.26 |
| deform_topology | 2.953 | 1.737 | — |
| climate_erosion_hydrology | 1.123 | 0.949 | 1.073 |
| deform | 2.675 | 0.7376 | 2.58 |
| faults | 0.7586 | 0.7029 | 1.157 |
| fluid_dynamics | 5.8e-06 | 7.44e-06 | — |
| gap_fill | 0.03521 | 0.81 | — |
| magma_transport | 2.003 | 1.824 | — |
| overlap_tracking | 0.0569 | 0.09329 | — |
| record_stats | 0.01316 | 0.003896 | — |
| resource_formation | 0.04724 | 0.01791 | — |
| sea_level | 0.252 | 0.1518 | 0.2404 |
| shift | 0.7078 | 0.4167 | 0.8096 |
| stranded_basins | 0.002237 | 7.32e-05 | — |
| topology | 0.1853 | 0.09648 | 0.2014 |
| volcanism | 0.09087 | 0.01863 | 0.0845 |

Index builds per step (mean):

- seed 804913535 lines: line_row_lookup 89.5, plate_node_kdtree 96.0, plate_outline 97.6, plate_outline_kdtree 95.8, world_node_kdtree 1.0
- seed 804913535 quad: plate_node_kdtree 92.4, plate_outline 98.0, plate_outline_kdtree 95.0, quad_adjacency 39.1, quad_boundary_loops 62.7, world_node_kdtree 1.0

The #147 profile ran under cProfile, so its absolute times are inflated.

