# Line vs quad parity report

**Verdict: FAIL** — configs long; seeds 804913535.

Reproduce:

    bin/debug/surface_parity.py paired --preset long --checkpoints-myr 30,120,240,352 --seeds 804913535 --jobs 2 --render --out ../analysis/issue249-final/repro-equivalent

Tolerances and gate definitions: docs/surface-parity.md. `info` rows are line-baseline
findings; they never gate the quad surface.

## Hard invariants

| gate | status | worst case |
|---|---|---|
| H1:quad_topology | pass | violations=0, audits=37, surface=quad, seed=804913535 |
| H2:quad_neighbours | pass | violations=0, audits=37, surface=quad, seed=804913535 |
| H3:quad_folded | pass | violations=0, audits=37, surface=quad, seed=804913535 |
| H4:fields | pass | violations=0, audits=37, surface=lines, seed=804913535 |
| H5:frames | pass | violations=0, audits=37, surface=lines, seed=804913535 |
| H6:derived_caches | pass | violations=0, audits=37, surface=lines, seed=804913535 |
| H7:revisions | pass | violations=0, audits=37, surface=lines, seed=804913535 |
| H8:load_round_trip | pass | checks=5, surface=lines, seed=804913535 |
| H9:load_continuation | pass | surface=lines, seed=804913535 |
| H10:quad_no_stacked_nodes | pass | surface=quad, seed=804913535 |
| H11:field_caps | pass | violations=0, surface=lines, seed=804913535 |

## Coverage and overlap

| gate | status | worst case |
|---|---|---|
| C1:uncovered | pass | age_myr=0, quad=0.00533, lines=0.0202, fail_above=0.0252, seed=804913535 |
| C2:void | pass | age_myr=0, quad=0, lines=0, fail_above=0.0005, seed=804913535 |
| C3:multiply_covered | pass | age_myr=0, quad=0.00524, warn_above=0.015, fail_above=0.02, surface=quad, seed=804913535 |
| C4:nodes_inside_other_plate | pass | age_myr=0, quad=0.007697, lines=0.01452, fail_above=0.01952, seed=804913535 |

## Conservation (actual node/cell areas)

| gate | status | worst case |
|---|---|---|
| K1:hc_volume_km3_drift | fail | age_myr=120, quad=-0.02366, lines=-0.1319, gap=0.1082, warn_above=0.05, fail_above=0.1, seed=804913535 |
| K2:continental_hc_volume_km3_drift | warn | age_myr=120, quad=-0.06602, lines=-0.1161, gap=0.05006, warn_above=0.05, fail_above=0.1, seed=804913535 |
| K3:quad_area_accounting | pass | age_myr=30, quad=1.002, expected=1.002, warn_above=0.01, seed=804913535 |
| K4:hm_volume_km3_drift | warn | age_myr=352, quad=-0.323, lines=-0.4189, gap=0.09598, warn_above=0.05, fail_above=0.1, seed=804913535 |

## Mesh quality (quad)

| gate | status | worst case |
|---|---|---|
| M1:anisotropic_fraction | pass | age_myr=0, quad=9.39e-05, lines=0.000306, warn_above=0.01, seed=804913535 |
| M2:row_alignment | pass | age_myr=0, quad=-0.2187, note=too few anisotropic nodes to have a direction, seed=804913535 |
| M3:thin_fraction | pass | age_myr=0, quad=0, lines=6.13e-05, warn_above=6.13e-05, seed=804913535 |
| M4:aspect_gt_4 | pass | age_myr=0, quad=0, warn_above=0.001, seed=804913535 |
| M5:skew_gt_45 | pass | age_myr=0, quad=0, warn_above=0.001, seed=804913535 |
| M6:fragments | warn | age_myr=30, quad=1, lines=0, seed=804913535 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=804913535 |
| S2:stability_air_temperature_mean_c | pass | quad=0.2001, lines=0.1266, warn_above=0.3032, seed=804913535 |
| S2:stability_land_fraction | pass | quad=0.001149, lines=0.001299, warn_above=0.004597, seed=804913535 |
| S2:stability_ocean_temperature_mean_c | pass | quad=0.01784, lines=0.01348, warn_above=0.07697, seed=804913535 |
| S2:stability_precipitation_mean_mm | pass | quad=5.775, lines=5.786, warn_above=16.57, seed=804913535 |
| S2:stability_sea_level_m | pass | quad=5.024, lines=6.652, warn_above=14.3, seed=804913535 |

## Derived indexes: HEALPix vs KD-tree

| gate | status | worst case |
|---|---|---|
| X1:healpix_all | pass | age_myr=0, samples=200000, quad_same_node=0.5569, lines_same_node=0.5548, quad_distance_p95=1.177, lines_distance_p95=1.176, seed=804913535 |
| X1:healpix_antimeridian | pass | age_myr=0, samples=3333, quad_same_node=0.5497, lines_same_node=0.5344, quad_distance_p95=1.226, lines_distance_p95=1.238, seed=804913535 |
| X1:healpix_coast | pass | age_myr=0, samples=13046, quad_same_node=0.5544, lines_same_node=0.5476, quad_distance_p95=1.201, lines_distance_p95=1.211, seed=804913535 |
| X1:healpix_hole | insufficient | age_myr=0, samples=0, seed=804913535 |
| X1:healpix_plate_boundary | pass | age_myr=0, samples=14669, quad_same_node=0.536, lines_same_node=0.5374, quad_distance_p95=1.244, lines_distance_p95=1.244, seed=804913535 |
| X1:healpix_pole | pass | age_myr=0, samples=3038, quad_same_node=0.5405, lines_same_node=0.5596, quad_distance_p95=1.214, lines_distance_p95=1.199, seed=804913535 |

## Performance

| gate | status | worst case |
|---|---|---|
| R1:deform_topology_s_per_step | pass | quad=1.016, lines=2.023, ratio=0.5025, warn_above=1.25, seed=804913535 |
| R2:step_total_s_per_step | pass | quad=3.028, lines=4.148, ratio=0.7301, warn_above=1.25, seed=804913535 |

## Ensemble parity

| gate | status | worst case |
|---|---|---|
| P1:continental_hc_volume_km3 | insufficient | age_myr=30, seeds=1, lines_mean=8.07e+09, quad_mean=8.14e+09, note=needs >= 3 seeds |
| P1:elevation_p05 | insufficient | age_myr=30, seeds=1, lines_mean=-5486, quad_mean=-5457, note=needs >= 3 seeds |
| P1:elevation_p50 | insufficient | age_myr=30, seeds=1, lines_mean=-4941, quad_mean=-4906, note=needs >= 3 seeds |
| P1:elevation_p95 | insufficient | age_myr=30, seeds=1, lines_mean=2116, quad_mean=2222, note=needs >= 3 seeds |
| P1:hc_volume_km3 | insufficient | age_myr=30, seeds=1, lines_mean=1.1e+10, quad_mean=1.11e+10, note=needs >= 3 seeds |
| P1:land_fraction | insufficient | age_myr=30, seeds=1, lines_mean=0.3143, quad_mean=0.3078, note=needs >= 3 seeds |
| P1:plates | insufficient | age_myr=30, seeds=1, lines_mean=11, quad_mean=11, note=needs >= 3 seeds |
| P1:sea_level_m | insufficient | age_myr=30, seeds=1, lines_mean=-590.6, quad_mean=-545, note=needs >= 3 seeds |

## Checkpoints, seed 804913535

| seed 804913535 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 352 Myr lines | 352 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|
| plates | 9 | 9 | 11 | 11 | 21 | 25 | 32 | 34 | 45 | 42 |
| nodes | 32643 | 31961 | 32774 | 32231 | 34170 | 32927 | 37010 | 32787 | 37482 | 33140 |
| land fraction | 0.2737 | 0.2741 | 0.3143 | 0.3078 | 0.2703 | 0.2299 | 0.1804 | 0.1481 | 0.1211 | 0.09691 |
| sea level (m) | 0 | 0 | -590.6 | -545 | -1379 | -1163 | -2429 | -1658 | -2760 | -2013 |
| Hc volume (km³) | 1.13e+10 | 1.13e+10 | 1.1e+10 | 1.11e+10 | 9.79e+09 | 1.1e+10 | 6.55e+09 | 8.41e+09 | 5.05e+09 | 7.32e+09 |
| continental Hc (km³) | 8.04e+09 | 8.04e+09 | 8.07e+09 | 8.14e+09 | 7.11e+09 | 7.51e+09 | 3.72e+09 | 3.95e+09 | 1.99e+09 | 2.6e+09 |
| elevation p05/p50/p95 (m) | -5222/-4318/1951 | -5222/-4320/1953 | -5486/-4941/2116 | -5457/-4906/2222 | -6229/-4946/2631 | -5825/-4869/6338 | -6426/-5285/1650 | -5886/-5126/5210 | -6462/-5518/519.9 | -5906/-5317/3562 |
| uncovered | 0.0202 | 0.00533 | 0.02459 | 0.005855 | 0.04796 | 0.01003 | 0.06039 | 0.01336 | 0.06531 | 0.01549 |
| multiply covered | 0.01241 | 0.00524 | 0.004545 | 0.00793 | 0.01961 | 0.01206 | 0.00424 | 0.01349 | 0.005185 | 0.01342 |
| void | 0 | 0 | 0 | 0 | 2.5e-05 | 0 | 0.00199 | 0 | 0.00022 | 0 |
| nodes inside other plate | 0.01452 | 0.007697 | 0.01404 | 0.00726 | 0.0427 | 0.009476 | 0.01467 | 0.006405 | 0.02196 | 0.005069 |
| stacked | 0 | 0 | 0.003967 | 0 | 0.08112 | 0 | 0.2461 | 0 | 0.2436 | 0 |
| anisotropic | 0.000306 | 9.39e-05 | 0.002044 | 0.001582 | 0.004858 | 0.002369 | 0.01594 | 0.00613 | 0.03052 | 0.006307 |
| row alignment | 0.8335 | -0.2187 | -0.1449 | -0.1173 | 0.2515 | 0.253 | 0.5551 | 0.1737 | 0.5107 | 0.1075 |
| thin | 6.13e-05 | 0 | 0.004516 | 0.001831 | 0.007053 | 0.003614 | 0.01021 | 0.008265 | 0.01318 | 0.007996 |
| one-node lines | 0 | — | 0.06045 | — | 0.08988 | — | 0.1823 | — | 0.2408 | — |
| quad one-cell-thin | — | 0.000219 | — | 0.004995 | — | 0.008959 | — | 0.0133 | — | 0.01629 |
| air temperature mean (°C) | 8.022 | 8.051 | 5.354 | 5.126 | 0.7672 | -8.147 | -1.416 | -6.035 | -12.73 | -14.43 |
| precipitation mean (mm) | 1004 | 1004 | 1031 | 1046 | 1093 | 1091 | 1199 | 1170 | 1215 | 1174 |
| river fraction of land | — | — | 0.1078 | 0.1119 | 0.1358 | 0.1992 | 0.1926 | 0.2466 | 0.2232 | 0.3145 |
| HEALPix same node | 0.5548 | 0.5569 | 0.5588 | 0.5545 | 0.5447 | 0.5529 | 0.53 | 0.5537 | 0.5327 | 0.5514 |

## Performance (mean seconds per step)

| s/step | seed 804913535 lines | seed 804913535 quad |
|---|---|---|
| step_total | 4.148 | 3.028 |
| deform_topology | 2.023 | 1.016 |
| climate_erosion_hydrology | 0.7756 | 0.6961 |
| deform | 1.906 | 0.4902 |
| faults | 0.6769 | 0.8688 |
| fluid_dynamics | 4.86e-06 | 3.82e-06 |
| gap_fill | 0.01505 | 0.438 |
| magma_transport | 0.1299 | 0.1009 |
| overlap_tracking | 0.02255 | 0.04584 |
| record_stats | 0.008349 | 0.004442 |
| resource_formation | 0.02683 | 0.01301 |
| sea_level | 0.07549 | 0.05704 |
| shift | 0.3782 | 0.2585 |
| stranded_basins | 2.88e-05 | 6.48e-05 |
| topology | 0.07971 | 0.04238 |
| volcanism | 0.04921 | 0.01158 |

Index builds per step (mean):

- seed 804913535 lines: line_row_lookup 70.0, plate_node_kdtree 70.1, plate_outline 76.1, plate_outline_kdtree 74.7, world_node_kdtree 1.0
- seed 804913535 quad: plate_node_kdtree 77.7, plate_outline 90.1, plate_outline_kdtree 81.4, quad_adjacency 49.2, quad_boundary_loops 62.1, world_node_kdtree 1.0

