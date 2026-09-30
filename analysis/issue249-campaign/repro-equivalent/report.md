# Line vs quad parity report

**Verdict: FAIL** — configs long; seeds 804913535.

Reproduce:

    bin/debug/surface_parity.py paired --preset long --checkpoints-myr 30,120,240,352 --seeds 804913535 --jobs 2 --render --out ../analysis/issue249-campaign/repro-equivalent

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
| H11:field_caps | warn (line baseline: 1 finding(s)) | violations=162, surface=quad, seed=804913535; first: {'detail': 'mantle_lithosphere_thickness_m: 1 values outside [2000, 240000] (range 1935.3..172025)', 'kind': 'field_bounds', 'plate_id': 6, 'step': 50, 'surface': 'quad'} |

## Coverage and overlap (quad vs paired line run)

| gate | status | worst case |
|---|---|---|
| C1:uncovered | pass | age_myr=0, quad=0.00533, lines=0.0202, fail_above=0.0252, seed=804913535 |
| C2:void | pass | age_myr=0, quad=0, lines=0, fail_above=0.0005, seed=804913535 |
| C3:multiply_covered | fail | age_myr=240, quad=0.01209, lines=0.003755, fail_above=0.008755, seed=804913535 |
| C4:nodes_inside_other_plate | pass | age_myr=0, quad=0.007697, lines=0.01452, fail_above=0.01952, seed=804913535 |

## Conservation (actual node/cell areas)

| gate | status | worst case |
|---|---|---|
| K1:hc_volume_km3_drift | warn | age_myr=120, quad=-0.1864, lines=-0.1334, gap=0.05298, warn_above=0.05, fail_above=0.1, seed=804913535 |
| K2:continental_hc_volume_km3_drift | fail | age_myr=120, quad=-0.2677, lines=-0.1009, gap=0.1668, warn_above=0.05, fail_above=0.1, seed=804913535 |
| K3:quad_area_accounting | pass | age_myr=30, quad=0.9931, expected=0.993, warn_above=0.01, seed=804913535 |
| K4:hm_volume_km3_drift | fail | age_myr=120, quad=-0.02599, lines=0.09112, gap=0.1171, warn_above=0.05, fail_above=0.1, seed=804913535 |

## Mesh quality (quad)

| gate | status | worst case |
|---|---|---|
| M1:anisotropic_fraction | pass | age_myr=0, quad=9.39e-05, lines=0.000306, warn_above=0.01, seed=804913535 |
| M2:row_alignment | pass | age_myr=0, quad=-0.2187, note=too few anisotropic nodes to have a direction, seed=804913535 |
| M3:thin_fraction | pass | age_myr=0, quad=0, lines=6.13e-05, warn_above=6.13e-05, seed=804913535 |
| M4:aspect_gt_4 | pass | age_myr=0, quad=0, warn_above=0.001, seed=804913535 |
| M5:skew_gt_45 | pass | age_myr=0, quad=0, warn_above=0.001, seed=804913535 |
| M6:fragments | warn | age_myr=240, quad=1, lines=0, seed=804913535 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=804913535 |
| S2:stability_air_temperature_mean_c | warn | quad=0.3826, lines=0.1216, warn_above=0.2932, seed=804913535 |
| S2:stability_land_fraction | pass | quad=0.002614, lines=0.00208, warn_above=0.006161, seed=804913535 |
| S2:stability_ocean_temperature_mean_c | pass | quad=0.01443, lines=0.01288, warn_above=0.07576, seed=804913535 |
| S2:stability_precipitation_mean_mm | pass | quad=6.29, lines=6.215, warn_above=17.43, seed=804913535 |
| S2:stability_sea_level_m | warn | quad=49.74, lines=8.265, warn_above=17.53, seed=804913535 |

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
| R1:deform_topology_s_per_step | pass | quad=0.7623, lines=2.133, ratio=0.3573, warn_above=1.25, seed=804913535 |
| R2:step_total_s_per_step | pass | quad=2.871, lines=4.361, ratio=0.6583, warn_above=1.25, seed=804913535 |

## Ensemble parity

| gate | status | worst case |
|---|---|---|
| P1:continental_hc_volume_km3 | insufficient | age_myr=30, seeds=1, lines_mean=8.05e+09, quad_mean=8.1e+09, note=needs >= 3 seeds |
| P1:elevation_p05 | insufficient | age_myr=30, seeds=1, lines_mean=-5481, quad_mean=-5458, note=needs >= 3 seeds |
| P1:elevation_p50 | insufficient | age_myr=30, seeds=1, lines_mean=-4943, quad_mean=-4907, note=needs >= 3 seeds |
| P1:elevation_p95 | insufficient | age_myr=30, seeds=1, lines_mean=2122, quad_mean=2254, note=needs >= 3 seeds |
| P1:hc_volume_km3 | insufficient | age_myr=30, seeds=1, lines_mean=1.1e+10, quad_mean=1.1e+10, note=needs >= 3 seeds |
| P1:land_fraction | insufficient | age_myr=30, seeds=1, lines_mean=0.3132, quad_mean=0.3095, note=needs >= 3 seeds |
| P1:plates | insufficient | age_myr=30, seeds=1, lines_mean=11, quad_mean=11, note=needs >= 3 seeds |
| P1:sea_level_m | insufficient | age_myr=30, seeds=1, lines_mean=-616.5, quad_mean=-526.1, note=needs >= 3 seeds |

## Checkpoints, seed 804913535

| seed 804913535 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 352 Myr lines | 352 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|
| plates | 9 | 9 | 11 | 11 | 21 | 21 | 32 | 28 | 37 | 31 |
| nodes | 32643 | 31961 | 32831 | 31896 | 33475 | 32859 | 35799 | 32672 | 39091 | 32963 |
| land fraction | 0.2737 | 0.2741 | 0.3132 | 0.3095 | 0.2676 | 0.2173 | 0.1965 | 0.129 | 0.1946 | 0.0998 |
| sea level (m) | 0 | 0 | -616.5 | -526.1 | -1289 | -1352 | -2156 | -1884 | -2632 | -2135 |
| Hc volume (km³) | 1.13e+10 | 1.13e+10 | 1.1e+10 | 1.1e+10 | 9.77e+09 | 9.17e+09 | 6.92e+09 | 6.75e+09 | 5.63e+09 | 6.12e+09 |
| continental Hc (km³) | 8.04e+09 | 8.04e+09 | 8.05e+09 | 8.1e+09 | 7.23e+09 | 5.89e+09 | 4.39e+09 | 2.83e+09 | 3.19e+09 | 1.86e+09 |
| elevation p05/p50/p95 (m) | -5222/-4318/1951 | -5222/-4320/1953 | -5481/-4943/2122 | -5458/-4907/2254 | -6260/-4968/2712 | -5864/-4907/2957 | -6442/-5180/1920 | -5946/-5207/463.1 | -6375/-5256/1966 | -5954/-5307/-272.6 |
| uncovered | 0.0202 | 0.00533 | 0.0233 | 0.01206 | 0.04532 | 0.008345 | 0.07395 | 0.0108 | 0.06213 | 0.0126 |
| multiply covered | 0.01241 | 0.00524 | 0.004405 | 0.00508 | 0.007875 | 0.01162 | 0.003755 | 0.01209 | 0.00247 | 0.01197 |
| void | 0 | 0 | 0 | 0 | 0.000465 | 0 | 0.01132 | 0 | 0.000565 | 0 |
| nodes inside other plate | 0.01452 | 0.007697 | 0.0127 | 0.005894 | 0.02094 | 0.01038 | 0.01751 | 0.005662 | 0.01213 | 0.005218 |
| stacked | 0 | 0 | 0.004203 | 0 | 0.06922 | 0 | 0.2102 | 0 | 0.3153 | 0 |
| anisotropic | 0.000306 | 9.39e-05 | 0.001767 | 0.000784 | 0.004421 | 0.002617 | 0.02327 | 0.003734 | 0.04686 | 0.003125 |
| row alignment | 0.8335 | -0.2187 | -0.1712 | -0.2685 | 0.03578 | -0.0197 | 0.3696 | 0.03376 | 0.721 | -0.08226 |
| thin | 6.13e-05 | 0 | 0.001828 | 0.001035 | 0.0046 | 0.002769 | 0.01497 | 0.005111 | 0.02425 | 0.004247 |
| one-node lines | 0 | — | 0.05845 | — | 0.1028 | — | 0.2794 | — | 0.2224 | — |
| quad one-cell-thin | — | 0.000219 | — | 0.003386 | — | 0.007852 | — | 0.01108 | — | 0.01077 |
| air temperature mean (°C) | 8.022 | 8.051 | 5.29 | 4.991 | 0.376 | 0.5143 | 1.538 | 2.254 | -10.65 | -7.998 |
| precipitation mean (mm) | 1004 | 1004 | 1035 | 1042 | 1068 | 1128 | 1184 | 1213 | 1151 | 1224 |
| river fraction of land | — | — | 0.1103 | 0.1068 | 0.1735 | 0.1837 | 0.1582 | 0.2295 | 0.1807 | 0.1868 |
| HEALPix same node | 0.5548 | 0.5569 | 0.5571 | 0.5524 | 0.5491 | 0.5531 | 0.5327 | 0.5494 | 0.5272 | 0.5521 |

## Performance (mean seconds per step)

| s/step | seed 804913535 lines | seed 804913535 quad |
|---|---|---|
| step_total | 4.361 | 2.871 |
| deform_topology | 2.133 | 0.7623 |
| climate_erosion_hydrology | 0.838 | 0.79 |
| deform | 1.996 | 0.5243 |
| faults | 0.7016 | 0.8799 |
| fluid_dynamics | 7.04e-06 | 3.64e-06 |
| gap_fill | 0.01577 | 0.1473 |
| magma_transport | 0.1383 | 0.08812 |
| overlap_tracking | 0.02248 | 0.04146 |
| record_stats | 0.008463 | 0.004677 |
| resource_formation | 0.02965 | 0.0137 |
| sea_level | 0.07462 | 0.05918 |
| shift | 0.3933 | 0.2633 |
| stranded_basins | 1.88e-05 | 1.61e-05 |
| topology | 0.09874 | 0.04925 |
| volcanism | 0.04157 | 0.008829 |

Index builds per step (mean):

- seed 804913535 lines: line_row_lookup 62.2, plate_node_kdtree 62.3, plate_outline 68.2, plate_outline_kdtree 66.7, world_node_kdtree 1.0
- seed 804913535 quad: plate_node_kdtree 48.4, plate_outline 62.4, plate_outline_kdtree 49.6, quad_adjacency 27.9, quad_boundary_loops 39.5, world_node_kdtree 1.0

