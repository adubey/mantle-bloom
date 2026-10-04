# Line vs quad parity report

**Verdict: FAIL** — configs long; seeds 1, 2.

Reproduce:

    bin/debug/surface_parity.py paired --preset long --node-density 0.5 --checkpoints-myr 100,200,400,600,800 --audit-every 10 --samples 100000 --seeds 1,2 --jobs 2 --out ../analysis/issue249-campaign3/stress

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
| C3:multiply_covered | warn | age_myr=600, quad=0.01788, warn_above=0.015, fail_above=0.02, surface=quad, seed=1 |
| C4:nodes_inside_other_plate | pass | age_myr=0, quad=0.01637, lines=0.0187, fail_above=0.0237, seed=1 |

## Conservation (actual node/cell areas)

| gate | status | worst case |
|---|---|---|
| K1:hc_volume_km3_drift | fail | age_myr=400, quad=-0.1913, lines=-0.3541, gap=0.1628, warn_above=0.05, fail_above=0.1, seed=1 |
| K2:continental_hc_volume_km3_drift | fail | age_myr=400, quad=-0.424, lines=-0.5291, gap=0.1052, warn_above=0.05, fail_above=0.1, seed=1 |
| K3:quad_area_accounting | pass | age_myr=100, quad=0.9992, expected=0.9988, warn_above=0.01, seed=1 |
| K4:hm_volume_km3_drift | fail | age_myr=200, quad=0.2073, lines=0.07341, gap=0.1339, warn_above=0.05, fail_above=0.1, seed=1 |

## Mesh quality (quad)

| gate | status | worst case |
|---|---|---|
| M1:anisotropic_fraction | warn | age_myr=600, quad=0.01333, lines=0.04103, warn_above=0.01, seed=1 |
| M2:row_alignment | pass | age_myr=0, quad=0.1263, note=too few anisotropic nodes to have a direction, seed=1 |
| M3:thin_fraction | warn | age_myr=0, quad=0.000385, lines=0.000368, warn_above=0.000368, seed=1 |
| M4:aspect_gt_4 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M5:skew_gt_45 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M6:fragments | warn | age_myr=600, quad=2, lines=0, seed=1 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=1 |
| S2:stability_air_temperature_mean_c | pass | quad=0.2497, lines=0.1384, warn_above=0.3267, seed=1 |
| S2:stability_land_fraction | pass | quad=0.000586, lines=0.001053, warn_above=0.004107, seed=1 |
| S2:stability_ocean_temperature_mean_c | pass | quad=0.005193, lines=0.01015, warn_above=0.07029, seed=1 |
| S2:stability_precipitation_mean_mm | pass | quad=5.797, lines=5.448, warn_above=15.9, seed=1 |
| S2:stability_sea_level_m | pass | quad=5.161, lines=12.21, warn_above=25.42, seed=1 |

## Derived indexes: HEALPix vs KD-tree

| gate | status | worst case |
|---|---|---|
| X1:healpix_all | pass | age_myr=0, samples=100000, quad_same_node=0.5058, lines_same_node=0.5007, quad_distance_p95=1.112, lines_distance_p95=1.074, seed=1 |
| X1:healpix_antimeridian | pass | age_myr=0, samples=1668, quad_same_node=0.509, lines_same_node=0.473, quad_distance_p95=1.227, lines_distance_p95=1.141, seed=1 |
| X1:healpix_coast | pass | age_myr=0, samples=8727, quad_same_node=0.4994, lines_same_node=0.4906, quad_distance_p95=1.128, lines_distance_p95=1.098, seed=1 |
| X1:healpix_hole | warn | age_myr=800, samples=333, quad_same_node=0.4955, lines_same_node=0.4204, quad_distance_p95=1.697, lines_distance_p95=1.296, seed=2 |
| X1:healpix_plate_boundary | pass | age_myr=0, samples=12435, quad_same_node=0.4913, lines_same_node=0.485, quad_distance_p95=1.193, lines_distance_p95=1.163, seed=1 |
| X1:healpix_pole | pass | age_myr=0, samples=1520, quad_same_node=0.5441, lines_same_node=0.4783, quad_distance_p95=1.209, lines_distance_p95=1.128, seed=1 |

## Performance

| gate | status | worst case |
|---|---|---|
| R1:deform_topology_s_per_step | pass | quad=0.5134, lines=1.045, ratio=0.4913, warn_above=1.25, seed=1 |
| R2:step_total_s_per_step | pass | quad=1.632, lines=2.243, ratio=0.7276, warn_above=1.25, seed=1 |

## Ensemble parity

| gate | status | worst case |
|---|---|---|
| P1:continental_hc_volume_km3 | insufficient | age_myr=100, seeds=2, lines_mean=5.48e+09, quad_mean=6.07e+09, note=needs >= 3 seeds |
| P1:elevation_p05 | insufficient | age_myr=100, seeds=2, lines_mean=-6091, quad_mean=-5838, note=needs >= 3 seeds |
| P1:elevation_p50 | insufficient | age_myr=100, seeds=2, lines_mean=-5110, quad_mean=-5021, note=needs >= 3 seeds |
| P1:elevation_p95 | insufficient | age_myr=100, seeds=2, lines_mean=2760, quad_mean=3804, note=needs >= 3 seeds |
| P1:hc_volume_km3 | insufficient | age_myr=100, seeds=2, lines_mean=8.98e+09, quad_mean=9.62e+09, note=needs >= 3 seeds |
| P1:land_fraction | insufficient | age_myr=100, seeds=2, lines_mean=0.2156, quad_mean=0.2161, note=needs >= 3 seeds |
| P1:plates | insufficient | age_myr=100, seeds=2, lines_mean=23, quad_mean=22.5, note=needs >= 3 seeds |
| P1:sea_level_m | insufficient | age_myr=100, seeds=2, lines_mean=-1016, quad_mean=-710.9, note=needs >= 3 seeds |

## Checkpoints, seed 1

| seed 1 | 0 Myr lines | 0 Myr quad | 100 Myr lines | 100 Myr quad | 200 Myr lines | 200 Myr quad | 400 Myr lines | 400 Myr quad | 600 Myr lines | 600 Myr quad | 800 Myr lines | 800 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 14 | 14 | 17 | 21 | 31 | 28 | 49 | 42 | 49 | 47 | 57 | 40 |
| nodes | 16312 | 15575 | 16234 | 16007 | 17998 | 16143 | 18901 | 16201 | 19278 | 16055 | 19697 | 16037 |
| land fraction | 0.2784 | 0.2788 | 0.3063 | 0.2933 | 0.3112 | 0.254 | 0.2152 | 0.1629 | 0.2281 | 0.1088 | 0.2882 | 0.07236 |
| sea level (m) | 0 | 0 | -807.2 | -750.8 | -1569 | -1030 | -2193 | -1584 | -2384 | -1886 | -2376 | -1950 |
| Hc volume (km³) | 1.15e+10 | 1.15e+10 | 1.07e+10 | 1.11e+10 | 9.7e+09 | 1.08e+10 | 7.4e+09 | 9.3e+09 | 6.62e+09 | 7.55e+09 | 6.38e+09 | 6.32e+09 |
| continental Hc (km³) | 7.62e+09 | 7.66e+09 | 7.7e+09 | 7.88e+09 | 6.8e+09 | 7.15e+09 | 3.59e+09 | 4.41e+09 | 2.87e+09 | 2.12e+09 | 3.18e+09 | 1e+09 |
| elevation p05/p50/p95 (m) | -4770/-4356/2947 | -4769/-4355/2989 | -6067/-5019/3888 | -5800/-4928/5469 | -6164/-4722/4849 | -5867/-4909/6338 | -6158/-5073/3787 | -5895/-5191/6241 | -6307/-5183/2901 | -5933/-5285/4205 | -6291/-5051/3361 | -5932/-5256/-280.1 |
| uncovered | 0.03196 | 0.00836 | 0.04848 | 0.01199 | 0.07714 | 0.01817 | 0.1066 | 0.02141 | 0.115 | 0.02419 | 0.1213 | 0.02146 |
| multiply covered | 0.01523 | 0.01004 | 0.0066 | 0.01078 | 0.00493 | 0.01483 | 0.00348 | 0.01496 | 0.00426 | 0.01788 | 0.00516 | 0.01353 |
| void | 0 | 0 | 1e-05 | 0 | 0.00033 | 0 | 0.00025 | 0 | 0.0033 | 0 | 0.00304 | 0 |
| nodes inside other plate | 0.0187 | 0.01637 | 0.01657 | 0.00681 | 0.01911 | 0.01059 | 0.0228 | 0.004506 | 0.025 | 0.00735 | 0.02848 | 0.002307 |
| stacked | 0 | 0 | 0.01404 | 0 | 0.1981 | 0 | 0.2813 | 0 | 0.3264 | 0 | 0.353 | 0 |
| anisotropic | 0.001042 | 0.000321 | 0.007022 | 0.003873 | 0.02789 | 0.005575 | 0.03243 | 0.00858 | 0.04103 | 0.01333 | 0.05574 | 0.009166 |
| row alignment | 0.6607 | 0.1263 | -0.02148 | 0.05008 | 0.517 | 0.04178 | 0.2924 | -0.05976 | 0.3641 | 0.1555 | 0.4006 | 0.04777 |
| thin | 0.000368 | 0.000385 | 0.006591 | 0.00656 | 0.01083 | 0.008672 | 0.02518 | 0.01222 | 0.02723 | 0.01626 | 0.0333 | 0.01228 |
| one-node lines | 0.01147 | — | 0.1564 | — | 0.2166 | — | 0.313 | — | 0.3251 | — | 0.3621 | — |
| quad one-cell-thin | — | 0.00077 | — | 0.007497 | — | 0.01456 | — | 0.01691 | — | 0.02199 | — | 0.02207 |
| air temperature mean (°C) | 4.512 | 4.486 | 6.507 | 3.908 | 4.607 | -0.3644 | 7.901 | -6.334 | 10.47 | -4.056 | 11.8 | 0.3103 |
| precipitation mean (mm) | 987.2 | 987.6 | 1034 | 1018 | 1107 | 1063 | 1185 | 1100 | 1187 | 1176 | 1212 | 1213 |
| river fraction of land | — | — | 0.1482 | 0.1715 | 0.15 | 0.2123 | 0.1943 | 0.3452 | 0.1559 | 0.3447 | 0.167 | 0.2581 |
| HEALPix same node | 0.5007 | 0.5058 | 0.4986 | 0.5007 | 0.4742 | 0.4914 | 0.4647 | 0.4979 | 0.4575 | 0.497 | 0.4521 | 0.497 |

## Checkpoints, seed 2

| seed 2 | 0 Myr lines | 0 Myr quad | 100 Myr lines | 100 Myr quad | 200 Myr lines | 200 Myr quad | 400 Myr lines | 400 Myr quad | 600 Myr lines | 600 Myr quad | 800 Myr lines | 800 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 18 | 18 | 29 | 24 | 44 | 35 | 54 | 49 | 61 | 52 | 75 | 59 |
| nodes | 16331 | 15491 | 17255 | 15809 | 18109 | 16086 | 20022 | 15910 | 21572 | 16165 | 21307 | 15851 |
| land fraction | 0.1401 | 0.1396 | 0.125 | 0.1388 | 0.1321 | 0.1102 | 0.08913 | 0.05369 | 0.07346 | 0.048 | 0.04879 | 0.02029 |
| sea level (m) | 0 | 0 | -1225 | -670.9 | -1584 | -810.1 | -2156 | -1337 | -2469 | -1255 | -2419 | -1507 |
| Hc volume (km³) | 8.89e+09 | 8.88e+09 | 7.23e+09 | 8.18e+09 | 6.55e+09 | 7.72e+09 | 4.41e+09 | 6.08e+09 | 3.2e+09 | 5.65e+09 | 3.66e+09 | 5.04e+09 |
| continental Hc (km³) | 4.27e+09 | 4.26e+09 | 3.26e+09 | 4.25e+09 | 2.99e+09 | 3.25e+09 | 1.5e+09 | 1.56e+09 | 6.74e+08 | 8.05e+08 | 5.43e+08 | 5.04e+08 |
| elevation p05/p50/p95 (m) | -5194/-4438/1652 | -5214/-4440/1648 | -6114/-5200/1632 | -5876/-5114/2139 | -6243/-5322/2144 | -5907/-5161/3692 | -6353/-5325/-333 | -5940/-5401/-741.2 | -6298/-5367/-242.5 | -5918/-5346/-1312 | -6301/-5336/-2717 | -5934/-5432/-4305 |
| uncovered | 0.03441 | 0.00999 | 0.07343 | 0.0151 | 0.09135 | 0.02025 | 0.1137 | 0.02561 | 0.1338 | 0.02433 | 0.1361 | 0.02725 |
| multiply covered | 0.1008 | 0.01114 | 0.0155 | 0.01246 | 0.01209 | 0.01398 | 0.01119 | 0.01709 | 0.00386 | 0.01587 | 0.00462 | 0.01607 |
| void | 0 | 0 | 0 | 0 | 6e-05 | 0 | 0.00202 | 0 | 0.00254 | 0 | 0 | 0 |
| nodes inside other plate | 0.08658 | 0.01607 | 0.02521 | 0.007401 | 0.03617 | 0.004289 | 0.03486 | 0.003394 | 0.02874 | 0.002413 | 0.03628 | 0.002019 |
| stacked | 0 | 0 | 0.1376 | 0 | 0.1969 | 0 | 0.3685 | 0 | 0.482 | 0 | 0.4298 | 0 |
| anisotropic | 0.000857 | 0.000258 | 0.01008 | 0.005566 | 0.02949 | 0.006217 | 0.0433 | 0.01125 | 0.06091 | 0.009774 | 0.065 | 0.01255 |
| row alignment | 0.3744 | -0.3999 | 0.3351 | 0.06938 | 0.3449 | 0.1262 | 0.2916 | 0.0283 | 0.4419 | -0.06545 | 0.3384 | -0.02365 |
| thin | 0.000306 | 0.000323 | 0.007882 | 0.006832 | 0.02093 | 0.009014 | 0.02797 | 0.01691 | 0.03676 | 0.01281 | 0.0429 | 0.01533 |
| one-node lines | 0.01746 | — | 0.1022 | — | 0.2555 | — | 0.3206 | — | 0.3831 | — | 0.3917 | — |
| quad one-cell-thin | — | 0.001162 | — | 0.01069 | — | 0.01517 | — | 0.02439 | — | 0.02066 | — | 0.02568 |
| air temperature mean (°C) | -3.448 | -3.523 | 2.966 | -0.434 | -0.1016 | -7.218 | 1.631 | -7.148 | 8.026 | 4.04 | 11.61 | -2.594 |
| precipitation mean (mm) | 1133 | 1134 | 1182 | 1158 | 1196 | 1173 | 1245 | 1247 | 1272 | 1244 | 1263 | 1264 |
| river fraction of land | — | — | 0.2235 | 0.255 | 0.203 | 0.2398 | 0.2268 | 0.3867 | 0.25 | 0.2918 | 0.2095 | 0.3625 |
| HEALPix same node | 0.5007 | 0.5039 | 0.4851 | 0.4994 | 0.4742 | 0.4971 | 0.45 | 0.4998 | 0.4327 | 0.4969 | 0.4374 | 0.496 |

## Performance (mean seconds per step)

| s/step | seed 1 lines | seed 1 quad | seed 2 lines | seed 2 quad |
|---|---|---|---|---|
| step_total | 2.243 | 1.632 | 2.735 | 1.735 |
| deform_topology | 1.045 | 0.5134 | 1.414 | 0.5981 |
| climate_erosion_hydrology | 0.4355 | 0.4087 | 0.4062 | 0.368 |
| deform | 0.981 | 0.2219 | 1.337 | 0.25 |
| faults | 0.4641 | 0.5226 | 0.5502 | 0.5687 |
| fluid_dynamics | 1.66e-06 | 1.36e-06 | 1.52e-06 | 1.71e-06 |
| gap_fill | 0.005647 | 0.2437 | 0.005744 | 0.2952 |
| magma_transport | 0.0444 | 0.03975 | 0.0415 | 0.03725 |
| overlap_tracking | 0.01098 | 0.02442 | 0.0135 | 0.02863 |
| record_stats | 0.00486 | 0.002868 | 0.00523 | 0.002569 |
| resource_formation | 0.01135 | 0.003769 | 0.01296 | 0.002847 |
| sea_level | 0.02762 | 0.02102 | 0.03368 | 0.02319 |
| shift | 0.183 | 0.114 | 0.2347 | 0.1269 |
| stranded_basins | 7.36e-06 | 8.57e-06 | 6.81e-06 | 5.54e-06 |
| topology | 0.04735 | 0.02339 | 0.05812 | 0.02423 |
| volcanism | 0.0248 | 0.005326 | 0.03286 | 0.006172 |

Index builds per step (mean):

- seed 1 lines: line_row_lookup 107.7, plate_node_kdtree 107.9, plate_outline 116.2, plate_outline_kdtree 114.5, world_node_kdtree 1.0
- seed 1 quad: plate_node_kdtree 93.2, plate_outline 99.6, plate_outline_kdtree 96.3, quad_adjacency 39.4, quad_boundary_loops 64.5, world_node_kdtree 1.0
- seed 2 lines: line_row_lookup 139.7, plate_node_kdtree 139.9, plate_outline 151.2, plate_outline_kdtree 148.3, world_node_kdtree 1.0
- seed 2 quad: plate_node_kdtree 116.6, plate_outline 123.2, plate_outline_kdtree 119.3, quad_adjacency 46.7, quad_boundary_loops 80.1, world_node_kdtree 1.0

