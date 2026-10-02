# Line vs quad parity report

**Verdict: FAIL** — configs long; seeds 1, 2.

Reproduce:

    bin/debug/surface_parity.py paired --preset long --node-density 0.5 --checkpoints-myr 100,200,400,600,800 --audit-every 10 --samples 100000 --seeds 1,2 --jobs 2 --out ../analysis/issue249-final/stress

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
| H11:field_caps | pass | violations=0, surface=lines, seed=1 |

## Coverage and overlap

| gate | status | worst case |
|---|---|---|
| C1:uncovered | pass | age_myr=0, quad=0.00836, lines=0.03196, fail_above=0.03696, seed=1 |
| C2:void | pass | age_myr=0, quad=0, lines=0, fail_above=0.0005, seed=1 |
| C3:multiply_covered | fail | age_myr=200, quad=0.02156, warn_above=0.015, fail_above=0.02, surface=quad, seed=1 |
| C4:nodes_inside_other_plate | pass | age_myr=0, quad=0.01637, lines=0.0187, fail_above=0.0237, seed=1 |

## Conservation (actual node/cell areas)

| gate | status | worst case |
|---|---|---|
| K1:hc_volume_km3_drift | fail | age_myr=200, quad=-0.03193, lines=-0.1536, gap=0.1217, warn_above=0.05, fail_above=0.1, seed=1 |
| K2:continental_hc_volume_km3_drift | fail | age_myr=400, quad=-0.3804, lines=-0.5291, gap=0.1488, warn_above=0.05, fail_above=0.1, seed=1 |
| K3:quad_area_accounting | pass | age_myr=100, quad=1.003, expected=1.003, warn_above=0.01, seed=1 |
| K4:hm_volume_km3_drift | fail | age_myr=600, quad=-0.2353, lines=-0.3382, gap=0.1028, warn_above=0.05, fail_above=0.1, seed=1 |

## Mesh quality (quad)

| gate | status | worst case |
|---|---|---|
| M1:anisotropic_fraction | warn | age_myr=400, quad=0.01051, lines=0.03243, warn_above=0.01, seed=1 |
| M2:row_alignment | pass | age_myr=0, quad=0.1263, note=too few anisotropic nodes to have a direction, seed=1 |
| M3:thin_fraction | warn | age_myr=0, quad=0.000385, lines=0.000368, warn_above=0.000368, seed=1 |
| M4:aspect_gt_4 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M5:skew_gt_45 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M6:fragments | warn | age_myr=100, quad=3, lines=0, seed=1 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=1 |
| S2:stability_air_temperature_mean_c | pass | quad=0.2572, lines=0.1384, warn_above=0.3267, seed=1 |
| S2:stability_land_fraction | pass | quad=0.000692, lines=0.001053, warn_above=0.004107, seed=1 |
| S2:stability_ocean_temperature_mean_c | pass | quad=0.00726, lines=0.01015, warn_above=0.07029, seed=1 |
| S2:stability_precipitation_mean_mm | pass | quad=5.291, lines=5.448, warn_above=15.9, seed=1 |
| S2:stability_sea_level_m | pass | quad=6.619, lines=12.21, warn_above=25.42, seed=1 |

## Derived indexes: HEALPix vs KD-tree

| gate | status | worst case |
|---|---|---|
| X1:healpix_all | pass | age_myr=0, samples=100000, quad_same_node=0.5058, lines_same_node=0.5007, quad_distance_p95=1.112, lines_distance_p95=1.074, seed=1 |
| X1:healpix_antimeridian | pass | age_myr=0, samples=1668, quad_same_node=0.509, lines_same_node=0.473, quad_distance_p95=1.227, lines_distance_p95=1.141, seed=1 |
| X1:healpix_coast | pass | age_myr=0, samples=8727, quad_same_node=0.4994, lines_same_node=0.4906, quad_distance_p95=1.128, lines_distance_p95=1.098, seed=1 |
| X1:healpix_hole | insufficient | age_myr=0, samples=0, seed=1 |
| X1:healpix_plate_boundary | pass | age_myr=0, samples=12435, quad_same_node=0.4913, lines_same_node=0.485, quad_distance_p95=1.193, lines_distance_p95=1.163, seed=1 |
| X1:healpix_pole | pass | age_myr=0, samples=1520, quad_same_node=0.5441, lines_same_node=0.4783, quad_distance_p95=1.209, lines_distance_p95=1.128, seed=1 |

## Performance

| gate | status | worst case |
|---|---|---|
| R1:deform_topology_s_per_step | pass | quad=0.6803, lines=1.341, ratio=0.5075, warn_above=1.25, seed=1 |
| R2:step_total_s_per_step | pass | quad=2.286, lines=2.884, ratio=0.7925, warn_above=1.25, seed=1 |

## Ensemble parity

| gate | status | worst case |
|---|---|---|
| P1:continental_hc_volume_km3 | insufficient | age_myr=100, seeds=2, lines_mean=5.48e+09, quad_mean=5.97e+09, note=needs >= 3 seeds |
| P1:elevation_p05 | insufficient | age_myr=100, seeds=2, lines_mean=-6091, quad_mean=-5796, note=needs >= 3 seeds |
| P1:elevation_p50 | insufficient | age_myr=100, seeds=2, lines_mean=-5110, quad_mean=-4952, note=needs >= 3 seeds |
| P1:elevation_p95 | insufficient | age_myr=100, seeds=2, lines_mean=2760, quad_mean=3663, note=needs >= 3 seeds |
| P1:hc_volume_km3 | insufficient | age_myr=100, seeds=2, lines_mean=8.98e+09, quad_mean=9.62e+09, note=needs >= 3 seeds |
| P1:land_fraction | insufficient | age_myr=100, seeds=2, lines_mean=0.2156, quad_mean=0.2146, note=needs >= 3 seeds |
| P1:plates | insufficient | age_myr=100, seeds=2, lines_mean=23, quad_mean=23, note=needs >= 3 seeds |
| P1:sea_level_m | insufficient | age_myr=100, seeds=2, lines_mean=-1016, quad_mean=-693.4, note=needs >= 3 seeds |

## Checkpoints, seed 1

| seed 1 | 0 Myr lines | 0 Myr quad | 100 Myr lines | 100 Myr quad | 200 Myr lines | 200 Myr quad | 400 Myr lines | 400 Myr quad | 600 Myr lines | 600 Myr quad | 800 Myr lines | 800 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 14 | 14 | 17 | 20 | 31 | 28 | 49 | 37 | 49 | 43 | 57 | 50 |
| nodes | 16312 | 15575 | 16234 | 16052 | 17998 | 16212 | 18901 | 16371 | 19278 | 16292 | 19697 | 16045 |
| land fraction | 0.2784 | 0.2788 | 0.3063 | 0.2959 | 0.3112 | 0.244 | 0.2152 | 0.1773 | 0.2281 | 0.1242 | 0.2882 | 0.09622 |
| sea level (m) | 0 | 0 | -807.2 | -685 | -1569 | -896 | -2193 | -1379 | -2384 | -1796 | -2376 | -1914 |
| Hc volume (km³) | 1.15e+10 | 1.15e+10 | 1.07e+10 | 1.11e+10 | 9.7e+09 | 1.11e+10 | 7.4e+09 | 9.66e+09 | 6.62e+09 | 8.23e+09 | 6.38e+09 | 6.72e+09 |
| continental Hc (km³) | 7.62e+09 | 7.66e+09 | 7.7e+09 | 7.9e+09 | 6.8e+09 | 7.11e+09 | 3.59e+09 | 4.75e+09 | 2.87e+09 | 3.14e+09 | 3.18e+09 | 1.65e+09 |
| elevation p05/p50/p95 (m) | -4770/-4356/2947 | -4769/-4355/2989 | -6067/-5019/3888 | -5779/-4781/5412 | -6164/-4722/4849 | -5842/-4650/6338 | -6158/-5073/3787 | -5875/-5036/5835 | -6307/-5183/2901 | -5893/-5246/5209 | -6291/-5051/3361 | -5919/-5289/1234 |
| uncovered | 0.03196 | 0.00836 | 0.04848 | 0.01141 | 0.07714 | 0.01612 | 0.1066 | 0.0209 | 0.115 | 0.02289 | 0.1213 | 0.02193 |
| multiply covered | 0.01523 | 0.01004 | 0.0066 | 0.01422 | 0.00493 | 0.02156 | 0.00348 | 0.02023 | 0.00426 | 0.02064 | 0.00516 | 0.01948 |
| void | 0 | 0 | 1e-05 | 0 | 0.00033 | 0 | 0.00025 | 0 | 0.0033 | 0 | 0.00304 | 0 |
| nodes inside other plate | 0.0187 | 0.01637 | 0.01657 | 0.01134 | 0.01911 | 0.01813 | 0.0228 | 0.01063 | 0.025 | 0.008532 | 0.02848 | 0.007167 |
| stacked | 0 | 0 | 0.01404 | 0 | 0.1981 | 0 | 0.2813 | 0 | 0.3264 | 0 | 0.353 | 0 |
| anisotropic | 0.001042 | 0.000321 | 0.007022 | 0.005358 | 0.02789 | 0.006908 | 0.03243 | 0.01051 | 0.04103 | 0.01295 | 0.05574 | 0.01433 |
| row alignment | 0.6607 | 0.1263 | -0.02148 | 0.06032 | 0.517 | 0.0128 | 0.2924 | -0.04601 | 0.3641 | -0.06682 | 0.4006 | 0.05634 |
| thin | 0.000368 | 0.000385 | 0.006591 | 0.008472 | 0.01083 | 0.009746 | 0.02518 | 0.01472 | 0.02723 | 0.01848 | 0.0333 | 0.02107 |
| one-node lines | 0.01147 | — | 0.1564 | — | 0.2166 | — | 0.313 | — | 0.3251 | — | 0.3621 | — |
| quad one-cell-thin | — | 0.00077 | — | 0.01364 | — | 0.01906 | — | 0.0248 | — | 0.02848 | — | 0.03023 |
| air temperature mean (°C) | 4.512 | 4.486 | 6.507 | 4.322 | 4.607 | -3.33 | 7.901 | -3.11 | 10.47 | -2.201 | 11.8 | 6.433 |
| precipitation mean (mm) | 987.2 | 987.6 | 1034 | 1015 | 1107 | 1065 | 1185 | 1115 | 1187 | 1173 | 1212 | 1221 |
| river fraction of land | — | — | 0.1482 | 0.1627 | 0.15 | 0.245 | 0.1943 | 0.2779 | 0.1559 | 0.3304 | 0.167 | 0.3076 |
| HEALPix same node | 0.5007 | 0.5058 | 0.4986 | 0.4996 | 0.4742 | 0.4959 | 0.4647 | 0.4957 | 0.4575 | 0.4954 | 0.4521 | 0.4965 |

## Checkpoints, seed 2

| seed 2 | 0 Myr lines | 0 Myr quad | 100 Myr lines | 100 Myr quad | 200 Myr lines | 200 Myr quad | 400 Myr lines | 400 Myr quad | 600 Myr lines | 600 Myr quad | 800 Myr lines | 800 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 18 | 18 | 29 | 26 | 44 | 32 | 54 | 53 | 61 | 54 | 75 | 60 |
| nodes | 16331 | 15491 | 17255 | 15974 | 18109 | 16121 | 20022 | 16233 | 21572 | 16022 | 21307 | 16172 |
| land fraction | 0.1401 | 0.1396 | 0.125 | 0.1333 | 0.1321 | 0.1062 | 0.08913 | 0.06972 | 0.07346 | 0.05422 | 0.04879 | 0.01516 |
| sea level (m) | 0 | 0 | -1225 | -701.8 | -1584 | -895.9 | -2156 | -1207 | -2469 | -1241 | -2419 | -1523 |
| Hc volume (km³) | 8.89e+09 | 8.88e+09 | 7.23e+09 | 8.12e+09 | 6.55e+09 | 7.59e+09 | 4.41e+09 | 6.58e+09 | 3.2e+09 | 5.88e+09 | 3.66e+09 | 5.04e+09 |
| continental Hc (km³) | 4.27e+09 | 4.26e+09 | 3.26e+09 | 4.04e+09 | 2.99e+09 | 3.13e+09 | 1.5e+09 | 1.82e+09 | 6.74e+08 | 8.13e+08 | 5.43e+08 | 4.49e+08 |
| elevation p05/p50/p95 (m) | -5194/-4438/1652 | -5214/-4440/1648 | -6114/-5200/1632 | -5814/-5123/1915 | -6243/-5322/2144 | -5896/-5234/2364 | -6353/-5325/-333 | -5898/-5346/1432 | -6298/-5367/-242.5 | -5926/-5344/-960.9 | -6301/-5336/-2717 | -5906/-5405/-4345 |
| uncovered | 0.03441 | 0.00999 | 0.07343 | 0.01677 | 0.09135 | 0.02012 | 0.1137 | 0.02329 | 0.1338 | 0.02466 | 0.1361 | 0.02765 |
| multiply covered | 0.1008 | 0.01114 | 0.0155 | 0.01514 | 0.01209 | 0.01756 | 0.01119 | 0.01928 | 0.00386 | 0.01972 | 0.00462 | 0.0223 |
| void | 0 | 0 | 0 | 0 | 6e-05 | 0 | 0.00202 | 0 | 0.00254 | 0 | 0 | 0 |
| nodes inside other plate | 0.08658 | 0.01607 | 0.02521 | 0.005822 | 0.03617 | 0.005893 | 0.03486 | 0.006653 | 0.02874 | 0.004556 | 0.03628 | 0.002597 |
| stacked | 0 | 0 | 0.1376 | 0 | 0.1969 | 0 | 0.3685 | 0 | 0.482 | 0 | 0.4298 | 0 |
| anisotropic | 0.000857 | 0.000258 | 0.01008 | 0.005259 | 0.02949 | 0.009677 | 0.0433 | 0.01244 | 0.06091 | 0.01623 | 0.065 | 0.01818 |
| row alignment | 0.3744 | -0.3999 | 0.3351 | 0.07911 | 0.3449 | -0.04854 | 0.2916 | 0.2128 | 0.4419 | 0.114 | 0.3384 | -0.01548 |
| thin | 0.000306 | 0.000323 | 0.007882 | 0.00795 | 0.02093 | 0.01365 | 0.02797 | 0.02014 | 0.03676 | 0.02297 | 0.0429 | 0.02541 |
| one-node lines | 0.01746 | — | 0.1022 | — | 0.2555 | — | 0.3206 | — | 0.3831 | — | 0.3917 | — |
| quad one-cell-thin | — | 0.001162 | — | 0.01728 | — | 0.02463 | — | 0.02871 | — | 0.03189 | — | 0.03791 |
| air temperature mean (°C) | -3.448 | -3.523 | 2.966 | 0.1726 | -0.1016 | -2.641 | 1.631 | -10.62 | 8.026 | -6.357 | 11.61 | -11.04 |
| precipitation mean (mm) | 1133 | 1134 | 1182 | 1156 | 1196 | 1195 | 1245 | 1205 | 1272 | 1238 | 1263 | 1260 |
| river fraction of land | — | — | 0.2235 | 0.2474 | 0.203 | 0.3071 | 0.2268 | 0.3386 | 0.25 | 0.3175 | 0.2095 | 0.3668 |
| HEALPix same node | 0.5007 | 0.5039 | 0.4851 | 0.499 | 0.4742 | 0.4954 | 0.45 | 0.4956 | 0.4327 | 0.4983 | 0.4374 | 0.4959 |

## Performance (mean seconds per step)

| s/step | seed 1 lines | seed 1 quad | seed 2 lines | seed 2 quad |
|---|---|---|---|---|
| step_total | 2.884 | 2.286 | 2.691 | 1.747 |
| deform_topology | 1.341 | 0.6803 | 1.391 | 0.5793 |
| climate_erosion_hydrology | 0.5653 | 0.5742 | 0.4001 | 0.3589 |
| deform | 1.258 | 0.316 | 1.316 | 0.2646 |
| faults | 0.5957 | 0.769 | 0.541 | 0.6074 |
| fluid_dynamics | 2.33e-06 | 2.3e-06 | 1.39e-06 | 1.37e-06 |
| gap_fill | 0.006662 | 0.2977 | 0.005586 | 0.2615 |
| magma_transport | 0.05693 | 0.06132 | 0.04078 | 0.03476 |
| overlap_tracking | 0.01414 | 0.03361 | 0.01323 | 0.02912 |
| record_stats | 0.006761 | 0.003877 | 0.004961 | 0.002567 |
| resource_formation | 0.01661 | 0.007032 | 0.01247 | 0.002596 |
| sea_level | 0.03303 | 0.02571 | 0.03337 | 0.02357 |
| shift | 0.2335 | 0.1552 | 0.2313 | 0.13 |
| stranded_basins | 1.25e-05 | 1.62e-05 | 5.91e-06 | 5.34e-06 |
| topology | 0.06217 | 0.03302 | 0.05705 | 0.02408 |
| volcanism | 0.03227 | 0.007803 | 0.03223 | 0.006619 |

Index builds per step (mean):

- seed 1 lines: line_row_lookup 107.7, plate_node_kdtree 107.9, plate_outline 116.2, plate_outline_kdtree 114.5, world_node_kdtree 1.0
- seed 1 quad: plate_node_kdtree 95.9, plate_outline 103.6, plate_outline_kdtree 99.2, quad_adjacency 45.5, quad_boundary_loops 68.4, world_node_kdtree 1.0
- seed 2 lines: line_row_lookup 139.7, plate_node_kdtree 139.9, plate_outline 151.2, plate_outline_kdtree 148.3, world_node_kdtree 1.0
- seed 2 quad: plate_node_kdtree 121.9, plate_outline 130.3, plate_outline_kdtree 124.9, quad_adjacency 55.7, quad_boundary_loops 86.0, world_node_kdtree 1.0

