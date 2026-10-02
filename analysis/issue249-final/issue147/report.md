# Line vs quad parity report

**Verdict: FAIL** — configs issue147; seeds 0.

Reproduce:

    bin/debug/surface_parity.py paired --preset issue147 --seeds 0 --jobs 1 --out ../analysis/issue249-final/issue147

Tolerances and gate definitions: docs/surface-parity.md. `info` rows are line-baseline
findings; they never gate the quad surface.

## Hard invariants

| gate | status | worst case |
|---|---|---|
| H1:quad_topology | pass | violations=0, audits=13, surface=quad, seed=0 |
| H2:quad_neighbours | pass | violations=0, audits=13, surface=quad, seed=0 |
| H3:quad_folded | pass | violations=0, audits=13, surface=quad, seed=0 |
| H4:fields | pass | violations=0, audits=13, surface=lines, seed=0 |
| H5:frames | pass | violations=0, audits=13, surface=lines, seed=0 |
| H6:derived_caches | pass | violations=0, audits=13, surface=lines, seed=0 |
| H7:revisions | pass | violations=0, audits=13, surface=lines, seed=0 |
| H8:load_round_trip | pass | checks=4, surface=lines, seed=0 |
| H9:load_continuation | pass | surface=lines, seed=0 |
| H10:quad_no_stacked_nodes | pass | surface=quad, seed=0 |
| H11:field_caps | pass | violations=0, surface=lines, seed=0 |

## Coverage and overlap

| gate | status | worst case |
|---|---|---|
| C1:uncovered | pass | age_myr=0, quad=0.003935, lines=0.01328, fail_above=0.01828, seed=0 |
| C2:void | pass | age_myr=0, quad=0, lines=0, fail_above=0.0005, seed=0 |
| C3:multiply_covered | pass | age_myr=0, quad=0.0038, warn_above=0.015, fail_above=0.02, surface=quad, seed=0 |
| C4:nodes_inside_other_plate | pass | age_myr=0, quad=0.005561, lines=0.01771, fail_above=0.02271, seed=0 |

## Conservation (actual node/cell areas)

| gate | status | worst case |
|---|---|---|
| K1:hc_volume_km3_drift | pass | age_myr=10, quad=-0.0272, lines=-0.0102, gap=0.017, warn_above=0.05, fail_above=0.1, seed=0 |
| K2:continental_hc_volume_km3_drift | pass | age_myr=10, quad=0.009266, lines=0.01438, gap=0.005116, warn_above=0.05, fail_above=0.1, seed=0 |
| K3:quad_area_accounting | pass | age_myr=10, quad=0.9994, expected=0.9994, warn_above=0.01, seed=0 |
| K4:hm_volume_km3_drift | fail | age_myr=60, quad=0.2295, lines=0.4192, gap=0.1897, warn_above=0.05, fail_above=0.1, seed=0 |

## Mesh quality (quad)

| gate | status | worst case |
|---|---|---|
| M1:anisotropic_fraction | pass | age_myr=0, quad=3.19e-05, lines=0.000145, warn_above=0.01, seed=0 |
| M2:row_alignment | pass | age_myr=0, quad=-0.01708, note=too few anisotropic nodes to have a direction, seed=0 |
| M3:thin_fraction | warn | age_myr=10, quad=0.002476, lines=0.001226, warn_above=0.001226, seed=0 |
| M4:aspect_gt_4 | pass | age_myr=0, quad=0, warn_above=0.001, seed=0 |
| M5:skew_gt_45 | pass | age_myr=0, quad=0, warn_above=0.001, seed=0 |
| M6:fragments | warn | age_myr=60, quad=1, lines=4, seed=0 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=0 |
| S2:stability_air_temperature_mean_c | pass | quad=0.03914, lines=0.04244, warn_above=0.1349, seed=0 |
| S2:stability_land_fraction | pass | quad=0.000269, lines=0.000342, warn_above=0.002683, seed=0 |
| S2:stability_ocean_temperature_mean_c | pass | quad=0.002085, lines=0.003097, warn_above=0.05619, seed=0 |
| S2:stability_precipitation_mean_mm | pass | quad=5.796, lines=5.676, warn_above=16.35, seed=0 |
| S2:stability_sea_level_m | pass | quad=0.9975, lines=4.117, warn_above=9.233, seed=0 |

## Derived indexes: HEALPix vs KD-tree

| gate | status | worst case |
|---|---|---|
| X1:healpix_all | pass | age_myr=0, samples=200000, quad_same_node=0.5458, lines_same_node=0.5559, quad_distance_p95=1.211, lines_distance_p95=1.183, seed=0 |
| X1:healpix_antimeridian | pass | age_myr=0, samples=3333, quad_same_node=0.5488, lines_same_node=0.5437, quad_distance_p95=1.257, lines_distance_p95=1.17, seed=0 |
| X1:healpix_coast | pass | age_myr=0, samples=4538, quad_same_node=0.5357, lines_same_node=0.5367, quad_distance_p95=1.274, lines_distance_p95=1.221, seed=0 |
| X1:healpix_hole | insufficient | age_myr=0, samples=0, seed=0 |
| X1:healpix_plate_boundary | pass | age_myr=0, samples=11096, quad_same_node=0.5333, lines_same_node=0.543, quad_distance_p95=1.266, lines_distance_p95=1.243, seed=0 |
| X1:healpix_pole | pass | age_myr=0, samples=3038, quad_same_node=0.549, lines_same_node=0.5346, quad_distance_p95=1.297, lines_distance_p95=1.267, seed=0 |

## Performance

| gate | status | worst case |
|---|---|---|
| R1:deform_topology_s_per_step | pass | quad=1.364, lines=1.562, ratio=0.8735, warn_above=1.25, seed=0 |
| R2:step_total_s_per_step | pass | quad=3.767, lines=8.451, ratio=0.4457, warn_above=1.25, seed=0 |

## Ensemble parity

| gate | status | worst case |
|---|---|---|
| P1:continental_hc_volume_km3 | insufficient | age_myr=10, seeds=1, lines_mean=4.45e+09, quad_mean=4.43e+09, note=needs >= 3 seeds |
| P1:elevation_p05 | insufficient | age_myr=10, seeds=1, lines_mean=-5199, quad_mean=-5226, note=needs >= 3 seeds |
| P1:elevation_p50 | insufficient | age_myr=10, seeds=1, lines_mean=-4864, quad_mean=-4867, note=needs >= 3 seeds |
| P1:elevation_p95 | insufficient | age_myr=10, seeds=1, lines_mean=962.9, quad_mean=780.3, note=needs >= 3 seeds |
| P1:hc_volume_km3 | insufficient | age_myr=10, seeds=1, lines_mean=8.47e+09, quad_mean=8.32e+09, note=needs >= 3 seeds |
| P1:land_fraction | insufficient | age_myr=10, seeds=1, lines_mean=0.1737, quad_mean=0.1627, note=needs >= 3 seeds |
| P1:plates | insufficient | age_myr=10, seeds=1, lines_mean=28, quad_mean=27, note=needs >= 3 seeds |
| P1:sea_level_m | insufficient | age_myr=10, seeds=1, lines_mean=-373.9, quad_mean=-384.2, note=needs >= 3 seeds |

## Checkpoints, seed 0

| seed 0 | 0 Myr lines | 0 Myr quad | 10 Myr lines | 10 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad |
|---|---|---|---|---|---|---|---|---|
| plates | 19 | 19 | 28 | 27 | 31 | 29 | 38 | 30 |
| nodes | 130587 | 125335 | 130475 | 125198 | 135191 | 125330 | 140244 | 125936 |
| land fraction | 0.153 | 0.153 | 0.1737 | 0.1627 | 0.1677 | 0.1543 | 0.156 | 0.1456 |
| sea level (m) | 0 | 0 | -373.9 | -384.2 | -695.5 | -542 | -865.9 | -547.3 |
| Hc volume (km³) | 8.56e+09 | 8.56e+09 | 8.47e+09 | 8.32e+09 | 8.69e+09 | 8.41e+09 | 8.79e+09 | 8.57e+09 |
| continental Hc (km³) | 4.39e+09 | 4.39e+09 | 4.45e+09 | 4.43e+09 | 4.51e+09 | 4.47e+09 | 4.4e+09 | 4.34e+09 |
| elevation p05/p50/p95 (m) | -5222/-4434/780.9 | -5222/-4434/781.1 | -5199/-4864/962.9 | -5226/-4867/780.3 | -5408/-5016/1399 | -5490/-5037/856.2 | -5454/-4953/1855 | -5582/-5026/1249 |
| uncovered | 0.01328 | 0.003935 | 0.02469 | 0.0058 | 0.02671 | 0.0064 | 0.03074 | 0.00723 |
| multiply covered | 0.01601 | 0.0038 | 0.00614 | 0.005195 | 0.001435 | 0.00658 | 0.001725 | 0.00611 |
| void | 0 | 0 | 2e-05 | 0 | 7.5e-05 | 0 | 0.001265 | 0 |
| nodes inside other plate | 0.01771 | 0.005561 | 0.007473 | 0.002548 | 0.004076 | 0.003599 | 0.005954 | 0.002176 |
| stacked | 0 | 0 | 0.008937 | 0 | 0.07282 | 0 | 0.1344 | 0 |
| anisotropic | 0.000145 | 3.19e-05 | 0.001019 | 0.001006 | 0.004386 | 0.001867 | 0.0137 | 0.002263 |
| row alignment | 0.6786 | -0.01708 | 0.06369 | 0.05686 | 0.4836 | 0.2201 | 0.7258 | -0.02588 |
| thin | 3.06e-05 | 2.39e-05 | 0.001226 | 0.002476 | 0.002352 | 0.003503 | 0.003608 | 0.004058 |
| one-node lines | 0.002784 | — | 0.0404 | — | 0.08951 | — | 0.152 | — |
| quad one-cell-thin | — | 0.00016 | — | 0.005144 | — | 0.006535 | — | 0.007679 |
| air temperature mean (°C) | 5.78 | 5.798 | 3.999 | 4.211 | 0.3401 | -0.2379 | -2.967 | -5.288 |
| precipitation mean (mm) | 1096 | 1096 | 1129 | 1137 | 1126 | 1140 | 1153 | 1160 |
| river fraction of land | — | — | 0.08387 | 0.09501 | 0.1018 | 0.1178 | 0.1135 | 0.1325 |
| HEALPix same node | 0.5559 | 0.5458 | 0.5557 | 0.5471 | 0.5477 | 0.5443 | 0.5414 | 0.5463 |

## Performance (mean seconds per step)

| s/step | seed 0 lines | seed 0 quad | #147 profile |
|---|---|---|---|
| step_total | 8.451 | 3.767 | 10.26 |
| deform_topology | 1.562 | 1.364 | — |
| climate_erosion_hydrology | 0.8317 | 0.742 | 1.073 |
| deform | 1.339 | 0.5454 | 2.58 |
| faults | 0.628 | 0.6846 | 1.157 |
| fluid_dynamics | 3.57e-06 | 1.35e-06 | — |
| gap_fill | 0.03938 | 0.6647 | — |
| magma_transport | 4.721 | 0.4461 | — |
| overlap_tracking | 0.04148 | 0.06877 | — |
| record_stats | 0.00758 | 0.00301 | — |
| resource_formation | 0.03229 | 0.01679 | — |
| sea_level | 0.2128 | 0.1789 | 0.2404 |
| shift | 0.4021 | 0.3149 | 0.8096 |
| stranded_basins | 1.32e-05 | 2.93e-05 | — |
| topology | 0.1416 | 0.08523 | 0.2014 |
| volcanism | 0.04971 | 0.01528 | 0.0845 |

Index builds per step (mean):

- seed 0 lines: line_row_lookup 80.4, plate_node_kdtree 80.5, plate_outline 90.6, plate_outline_kdtree 87.2, world_node_kdtree 1.0
- seed 0 quad: plate_node_kdtree 75.7, plate_outline 82.1, plate_outline_kdtree 79.9, quad_adjacency 33.5, quad_boundary_loops 53.7, world_node_kdtree 1.0

The #147 profile ran under cProfile, so its absolute times are inflated.

