# Line vs quad parity report

**Verdict: FAIL** — configs long; seeds 1, 2.

Reproduce:

    bin/debug/surface_parity.py paired --preset long --node-density 0.5 --checkpoints-myr 100,200,400,600,800 --audit-every 10 --samples 100000 --seeds 1,2 --jobs 2 --out ../analysis/issue249-campaign/stress

Tolerances and gate definitions: docs/surface-parity.md. `info` rows are line-baseline
findings; they never gate the quad surface.

## Hard invariants

| gate | status | worst case |
|---|---|---|
| H1:quad_topology | pass | violations=0, audits=81, surface=quad, seed=1 |
| H2:quad_neighbours | pass | violations=0, audits=81, surface=quad, seed=1 |
| H3:quad_folded | pass | violations=0, audits=81, surface=quad, seed=1 |
| H4:fields | pass | violations=0, audits=81, surface=lines, seed=1 |
| H5:frames | pass | violations=0, audits=81, surface=lines, seed=1 |
| H6:derived_caches | pass | violations=0, audits=81, surface=lines, seed=1 |
| H7:revisions | pass | violations=0, audits=81, surface=lines, seed=1 |
| H8:load_round_trip | pass | checks=6, surface=lines, seed=1 |
| H9:load_continuation | pass | surface=lines, seed=1 |
| H10:quad_no_stacked_nodes | pass | surface=quad, seed=1 |
| H11:field_caps | warn (line baseline: 2 finding(s)) | violations=312, surface=quad, seed=1; first: {'detail': 'mantle_lithosphere_thickness_m: 3 values outside [2000, 240000] (range 1696.08..240000)', 'kind': 'field_bounds', 'plate_id': 2, 'step': 60, 'surface': 'quad'} |

## Coverage and overlap (quad vs paired line run)

| gate | status | worst case |
|---|---|---|
| C1:uncovered | pass | age_myr=0, quad=0.00836, lines=0.03196, fail_above=0.03696, seed=1 |
| C2:void | pass | age_myr=0, quad=0, lines=0, fail_above=0.0005, seed=1 |
| C3:multiply_covered | fail | age_myr=100, quad=0.01215, lines=0.00663, fail_above=0.01163, seed=1 |
| C4:nodes_inside_other_plate | pass | age_myr=0, quad=0.01637, lines=0.0187, fail_above=0.0237, seed=1 |

## Conservation (actual node/cell areas)

| gate | status | worst case |
|---|---|---|
| K1:hc_volume_km3_drift | warn | age_myr=400, quad=-0.365, lines=-0.2938, gap=0.07125, warn_above=0.05, fail_above=0.1, seed=1 |
| K2:continental_hc_volume_km3_drift | fail | age_myr=200, quad=-0.1862, lines=-0.004331, gap=0.1818, warn_above=0.05, fail_above=0.1, seed=1 |
| K3:quad_area_accounting | pass | age_myr=100, quad=1.001, expected=1.001, warn_above=0.01, seed=1 |
| K4:hm_volume_km3_drift | fail | age_myr=600, quad=-0.2487, lines=-0.367, gap=0.1183, warn_above=0.05, fail_above=0.1, seed=1 |

## Mesh quality (quad)

| gate | status | worst case |
|---|---|---|
| M1:anisotropic_fraction | warn | age_myr=800, quad=0.01094, lines=0.09463, warn_above=0.01, seed=1 |
| M2:row_alignment | pass | age_myr=0, quad=0.1263, note=too few anisotropic nodes to have a direction, seed=1 |
| M3:thin_fraction | warn | age_myr=0, quad=0.000385, lines=0.000368, warn_above=0.000368, seed=1 |
| M4:aspect_gt_4 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M5:skew_gt_45 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M6:fragments | warn | age_myr=800, quad=3, lines=2, seed=1 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=1 |
| S2:stability_air_temperature_mean_c | warn | quad=2.358, lines=0.317, warn_above=0.6841, seed=2 |
| S2:stability_land_fraction | pass | quad=0.001486, lines=0.001218, warn_above=0.004436, seed=1 |
| S2:stability_ocean_temperature_mean_c | pass | quad=0.01168, lines=0.01004, warn_above=0.07009, seed=1 |
| S2:stability_precipitation_mean_mm | pass | quad=5.407, lines=5.327, warn_above=15.65, seed=1 |
| S2:stability_sea_level_m | warn | quad=68.6, lines=10.71, warn_above=22.42, seed=1 |

## Derived indexes: HEALPix vs KD-tree

| gate | status | worst case |
|---|---|---|
| X1:healpix_all | pass | age_myr=0, samples=100000, quad_same_node=0.5058, lines_same_node=0.5007, quad_distance_p95=1.112, lines_distance_p95=1.074, seed=1 |
| X1:healpix_antimeridian | pass | age_myr=0, samples=1668, quad_same_node=0.509, lines_same_node=0.473, quad_distance_p95=1.227, lines_distance_p95=1.141, seed=1 |
| X1:healpix_coast | warn | age_myr=600, samples=47, quad_same_node=0.4043, lines_same_node=0.4707, quad_distance_p95=0.9915, lines_distance_p95=1.277, seed=2 |
| X1:healpix_hole | insufficient | age_myr=0, samples=0, seed=1 |
| X1:healpix_plate_boundary | pass | age_myr=0, samples=12435, quad_same_node=0.4913, lines_same_node=0.485, quad_distance_p95=1.193, lines_distance_p95=1.163, seed=1 |
| X1:healpix_pole | pass | age_myr=0, samples=1520, quad_same_node=0.5441, lines_same_node=0.4783, quad_distance_p95=1.209, lines_distance_p95=1.128, seed=1 |

## Performance

| gate | status | worst case |
|---|---|---|
| R1:deform_topology_s_per_step | pass | quad=0.5196, lines=1.705, ratio=0.3048, warn_above=1.25, seed=1 |
| R2:step_total_s_per_step | pass | quad=2.206, lines=3.601, ratio=0.6127, warn_above=1.25, seed=1 |

## Ensemble parity

| gate | status | worst case |
|---|---|---|
| P1:continental_hc_volume_km3 | insufficient | age_myr=100, seeds=2, lines_mean=5.77e+09, quad_mean=5.83e+09, note=needs >= 3 seeds |
| P1:elevation_p05 | insufficient | age_myr=100, seeds=2, lines_mean=-6172, quad_mean=-5853, note=needs >= 3 seeds |
| P1:elevation_p50 | insufficient | age_myr=100, seeds=2, lines_mean=-5114, quad_mean=-4953, note=needs >= 3 seeds |
| P1:elevation_p95 | insufficient | age_myr=100, seeds=2, lines_mean=3400, quad_mean=2962, note=needs >= 3 seeds |
| P1:hc_volume_km3 | insufficient | age_myr=100, seeds=2, lines_mean=9.15e+09, quad_mean=9.31e+09, note=needs >= 3 seeds |
| P1:land_fraction | insufficient | age_myr=100, seeds=2, lines_mean=0.2147, quad_mean=0.2246, note=needs >= 3 seeds |
| P1:plates | insufficient | age_myr=100, seeds=2, lines_mean=21.5, quad_mean=21, note=needs >= 3 seeds |
| P1:sea_level_m | insufficient | age_myr=100, seeds=2, lines_mean=-1015, quad_mean=-817.7, note=needs >= 3 seeds |

## Checkpoints, seed 1

| seed 1 | 0 Myr lines | 0 Myr quad | 100 Myr lines | 100 Myr quad | 200 Myr lines | 200 Myr quad | 400 Myr lines | 400 Myr quad | 600 Myr lines | 600 Myr quad | 800 Myr lines | 800 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 14 | 14 | 16 | 17 | 27 | 26 | 43 | 36 | 53 | 35 | 55 | 43 |
| nodes | 16312 | 15575 | 16284 | 15948 | 17421 | 16194 | 18149 | 16166 | 20471 | 16312 | 21061 | 16084 |
| land fraction | 0.2784 | 0.2788 | 0.2973 | 0.2975 | 0.2994 | 0.251 | 0.2319 | 0.1309 | 0.2296 | 0.103 | 0.267 | 0.05756 |
| sea level (m) | 0 | 0 | -853.8 | -838 | -1313 | -1118 | -1942 | -1876 | -2606 | -2083 | -2657 | -2265 |
| Hc volume (km³) | 1.15e+10 | 1.15e+10 | 1.09e+10 | 1.07e+10 | 1.02e+10 | 9.82e+09 | 8.09e+09 | 7.3e+09 | 6.27e+09 | 6.31e+09 | 6.11e+09 | 5.43e+09 |
| continental Hc (km³) | 7.62e+09 | 7.66e+09 | 7.96e+09 | 7.62e+09 | 7.59e+09 | 6.24e+09 | 4.53e+09 | 2.6e+09 | 2.96e+09 | 1.39e+09 | 2.69e+09 | 8.84e+08 |
| elevation p05/p50/p95 (m) | -4770/-4356/2947 | -4769/-4355/2989 | -6120/-5030/4783 | -5810/-4825/4245 | -6191/-4629/5048 | -5898/-4738/4294 | -6307/-5130/4787 | -5955/-5188/2258 | -6210/-5138/3582 | -5975/-5265/650.3 | -6333/-5099/2988 | -5967/-5362/-2065 |
| uncovered | 0.03196 | 0.00836 | 0.0458 | 0.01115 | 0.07641 | 0.01431 | 0.09941 | 0.01809 | 0.1191 | 0.01809 | 0.1206 | 0.01888 |
| multiply covered | 0.01523 | 0.01004 | 0.00663 | 0.01215 | 0.01117 | 0.01696 | 0.00668 | 0.01916 | 0.00857 | 0.01681 | 0.00764 | 0.01843 |
| void | 0 | 0 | 0.00072 | 0 | 0.00201 | 0 | 0.00027 | 0 | 0.00034 | 0 | 0.00108 | 0 |
| nodes inside other plate | 0.0187 | 0.01637 | 0.0129 | 0.007148 | 0.02985 | 0.01136 | 0.03317 | 0.009712 | 0.03327 | 0.004414 | 0.03874 | 0.003544 |
| stacked | 0 | 0 | 0.01603 | 0 | 0.1498 | 0 | 0.1991 | 0 | 0.3911 | 0 | 0.4147 | 0 |
| anisotropic | 0.001042 | 0.000321 | 0.005895 | 0.002508 | 0.01458 | 0.005002 | 0.03229 | 0.008165 | 0.06179 | 0.008705 | 0.09463 | 0.01094 |
| row alignment | 0.6607 | 0.1263 | 0.09318 | -0.273 | 0.09374 | 0.06047 | 0.2161 | 0.07293 | 0.5047 | 0.04508 | 0.5155 | 0.0305 |
| thin | 0.000368 | 0.000385 | 0.005895 | 0.00627 | 0.009931 | 0.00636 | 0.02656 | 0.01052 | 0.03195 | 0.00987 | 0.05024 | 0.01418 |
| one-node lines | 0.01147 | — | 0.1409 | — | 0.1929 | — | 0.3258 | — | 0.3586 | — | 0.3993 | — |
| quad one-cell-thin | — | 0.00077 | — | 0.00997 | — | 0.0142 | — | 0.02134 | — | 0.01999 | — | 0.02443 |
| air temperature mean (°C) | 4.512 | 4.486 | 4.86 | 6.699 | 3.177 | 3.948 | 5.067 | -0.01078 | 9.993 | 9.038 | 9.713 | 17.57 |
| precipitation mean (mm) | 987.2 | 987.6 | 1003 | 1041 | 1065 | 1084 | 1146 | 1176 | 1209 | 1243 | 1232 | 1256 |
| river fraction of land | — | — | 0.175 | 0.158 | 0.1765 | 0.2036 | 0.2314 | 0.2377 | 0.1607 | 0.184 | 0.1803 | 0.1122 |
| HEALPix same node | 0.5007 | 0.5058 | 0.4981 | 0.5 | 0.4777 | 0.4961 | 0.4728 | 0.4965 | 0.4441 | 0.4974 | 0.4405 | 0.4921 |

## Checkpoints, seed 2

| seed 2 | 0 Myr lines | 0 Myr quad | 100 Myr lines | 100 Myr quad | 200 Myr lines | 200 Myr quad | 400 Myr lines | 400 Myr quad | 600 Myr lines | 600 Myr quad | 800 Myr lines | 800 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 18 | 18 | 27 | 25 | 47 | 35 | 59 | 44 | 60 | 49 | 67 | 46 |
| nodes | 16331 | 15491 | 17101 | 16054 | 18615 | 16406 | 20528 | 16066 | 20853 | 16256 | 22454 | 15763 |
| land fraction | 0.1401 | 0.1396 | 0.1321 | 0.1517 | 0.1132 | 0.1029 | 0.08161 | 0.02935 | 0.07339 | 5.98e-05 | 0.0902 | 0 |
| sea level (m) | 0 | 0 | -1177 | -797.5 | -1764 | -1095 | -2204 | -1598 | -2288 | -1805 | -2484 | -1691 |
| Hc volume (km³) | 8.89e+09 | 8.88e+09 | 7.43e+09 | 7.91e+09 | 5.88e+09 | 6.69e+09 | 4.48e+09 | 5.07e+09 | 3.7e+09 | 4.55e+09 | 3.71e+09 | 4.52e+09 |
| continental Hc (km³) | 4.27e+09 | 4.26e+09 | 3.58e+09 | 4.03e+09 | 2.51e+09 | 2.38e+09 | 1.26e+09 | 5.12e+08 | 8.88e+08 | 0 | 6.49e+08 | 0 |
| elevation p05/p50/p95 (m) | -5194/-4438/1652 | -5214/-4440/1648 | -6224/-5197/2018 | -5896/-5080/1679 | -6293/-5300/2066 | -5932/-5121/1220 | -6257/-5326/-266.6 | -5979/-5393/-4059 | -6283/-5286/-1099 | -5986/-5435/-4478 | -6253/-5299/-414.2 | -5982/-5431/-4450 |
| uncovered | 0.03441 | 0.00999 | 0.073 | 0.01328 | 0.1083 | 0.01846 | 0.1183 | 0.02064 | 0.1245 | 0.02131 | 0.1356 | 0.02189 |
| multiply covered | 0.1008 | 0.01114 | 0.01837 | 0.01625 | 0.00473 | 0.02013 | 0.00494 | 0.01958 | 0.00406 | 0.01998 | 0.00517 | 0.01993 |
| void | 0 | 0 | 0.00012 | 0 | 0.0017 | 0 | 0.00033 | 0 | 0.00048 | 0 | 0.00017 | 0 |
| nodes inside other plate | 0.08658 | 0.01607 | 0.03128 | 0.008721 | 0.02288 | 0.01036 | 0.02957 | 0.005851 | 0.03141 | 0.004368 | 0.04462 | 0.005202 |
| stacked | 0 | 0 | 0.1227 | 0 | 0.2688 | 0 | 0.3848 | 0 | 0.4122 | 0 | 0.4906 | 0 |
| anisotropic | 0.000857 | 0.000258 | 0.01275 | 0.004921 | 0.02235 | 0.009509 | 0.06606 | 0.01027 | 0.06613 | 0.01107 | 0.07767 | 0.01224 |
| row alignment | 0.3744 | -0.3999 | -0.07988 | 0.01681 | -0.07114 | -0.01802 | 0.4583 | 0.05223 | 0.3877 | -0.02759 | 0.3738 | -0.05817 |
| thin | 0.000306 | 0.000323 | 0.01257 | 0.006603 | 0.02251 | 0.01201 | 0.03619 | 0.01282 | 0.03515 | 0.01519 | 0.04926 | 0.0177 |
| one-node lines | 0.01746 | — | 0.1774 | — | 0.2291 | — | 0.3884 | — | 0.3877 | — | 0.4466 | — |
| quad one-cell-thin | — | 0.001162 | — | 0.01352 | — | 0.02072 | — | 0.02272 | — | 0.02565 | — | 0.02766 |
| air temperature mean (°C) | -3.448 | -3.523 | 2.7 | 3.592 | 0.2131 | 5.744 | 0.5058 | 13.76 | 1.34 | -6.521 | 8.794 | — |
| precipitation mean (mm) | 1133 | 1134 | 1172 | 1167 | 1202 | 1196 | 1254 | 1271 | 1260 | 1274 | 1252 | 1270 |
| river fraction of land | — | — | 0.2138 | 0.2132 | 0.2668 | 0.2106 | 0.1688 | 0.222 | 0.2953 | 1 | 0.1674 | 0 |
| HEALPix same node | 0.5007 | 0.5039 | 0.4824 | 0.4991 | 0.4609 | 0.4941 | 0.4469 | 0.4986 | 0.4399 | 0.4982 | 0.4255 | 0.4999 |

## Performance (mean seconds per step)

| s/step | seed 1 lines | seed 1 quad | seed 2 lines | seed 2 quad |
|---|---|---|---|---|
| step_total | 3.601 | 2.206 | 3.451 | 1.668 |
| deform_topology | 1.705 | 0.5196 | 1.759 | 0.4514 |
| climate_erosion_hydrology | 0.6935 | 0.6372 | 0.541 | 0.4056 |
| deform | 1.601 | 0.3398 | 1.657 | 0.2935 |
| faults | 0.7252 | 0.7955 | 0.6879 | 0.6226 |
| fluid_dynamics | 3.66e-06 | 2.88e-06 | 3.34e-06 | 1.51e-06 |
| gap_fill | 0.008574 | 0.1116 | 0.007001 | 0.1015 |
| magma_transport | 0.07951 | 0.05098 | 0.05937 | 0.01893 |
| overlap_tracking | 0.01731 | 0.03406 | 0.017 | 0.03012 |
| record_stats | 0.008041 | 0.003884 | 0.006602 | 0.002658 |
| resource_formation | 0.02393 | 0.007842 | 0.01966 | 0.003401 |
| sea_level | 0.03698 | 0.02829 | 0.03796 | 0.0247 |
| shift | 0.2908 | 0.1559 | 0.2968 | 0.1325 |
| stranded_basins | 1.73e-05 | 1.29e-05 | 1.1e-05 | 5.68e-06 |
| topology | 0.0774 | 0.03415 | 0.07764 | 0.02633 |
| volcanism | 0.03605 | 0.0062 | 0.04063 | 0.005156 |

Index builds per step (mean):

- seed 1 lines: line_row_lookup 102.0, plate_node_kdtree 102.1, plate_outline 109.8, plate_outline_kdtree 108.0, world_node_kdtree 1.0
- seed 1 quad: plate_node_kdtree 64.9, plate_outline 78.3, plate_outline_kdtree 66.2, quad_adjacency 27.3, quad_boundary_loops 47.3, world_node_kdtree 1.0
- seed 2 lines: line_row_lookup 142.0, plate_node_kdtree 142.2, plate_outline 153.7, plate_outline_kdtree 150.7, world_node_kdtree 1.0
- seed 2 quad: plate_node_kdtree 87.1, plate_outline 106.5, plate_outline_kdtree 88.9, quad_adjacency 36.2, quad_boundary_loops 65.5, world_node_kdtree 1.0

