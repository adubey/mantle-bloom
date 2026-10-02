# Line vs quad parity report

**Verdict: FAIL** — configs long; seeds 1, 2, 3, 4, 5.

Reproduce:

    bin/debug/surface_parity.py paired --preset long --seeds 1,2,3,4,5 --surfaces quad --jobs 5 --out ../analysis/issue257/long

Tolerances and gate definitions: docs/surface-parity.md. `info` rows are line-baseline
findings; they never gate the quad surface.

## Hard invariants

| gate | status | worst case |
|---|---|---|
| H1:quad_topology | pass | violations=0, audits=41, surface=quad, seed=1 |
| H2:quad_neighbours | pass | violations=0, audits=41, surface=quad, seed=1 |
| H3:quad_folded | pass | violations=0, audits=41, surface=quad, seed=1 |
| H4:fields | pass | violations=0, audits=41, surface=lines, seed=1 |
| H5:frames | pass | violations=0, audits=41, surface=lines, seed=1 |
| H6:derived_caches | pass | violations=0, audits=41, surface=lines, seed=1 |
| H7:revisions | pass | violations=0, audits=41, surface=lines, seed=1 |
| H8:load_round_trip | pass | checks=6, surface=lines, seed=1 |
| H9:load_continuation | pass | surface=lines, seed=1 |
| H10:quad_no_stacked_nodes | pass | surface=quad, seed=1 |
| H11:field_caps | info (line baseline: 5 finding(s)) | violations=420, surface=lines, seed=1; first: {'detail': 'mantle_lithosphere_thickness_m: 17 values outside [2000, 240000] (range 147.928..240000)', 'kind': 'field_bounds', 'plate_id': 2, 'step': 10, 'surface': 'lines'} |

## Coverage and overlap

| gate | status | worst case |
|---|---|---|
| C1:uncovered | pass | age_myr=0, quad=0.005985, lines=0.02311, fail_above=0.02811, seed=1 |
| C2:void | pass | age_myr=0, quad=0, lines=0, fail_above=0.0005, seed=1 |
| C3:multiply_covered | pass | age_myr=0, quad=0.006475, warn_above=0.015, fail_above=0.02, surface=quad, seed=1 |
| C4:nodes_inside_other_plate | pass | age_myr=0, quad=0.009653, lines=0.01618, fail_above=0.02118, seed=1 |

## Conservation (actual node/cell areas)

| gate | status | worst case |
|---|---|---|
| K1:hc_volume_km3_drift | fail | age_myr=400, quad=-0.4656, lines=-0.3156, gap=0.1501, warn_above=0.05, fail_above=0.1, seed=1 |
| K2:continental_hc_volume_km3_drift | fail | age_myr=240, quad=-0.3883, lines=-0.09438, gap=0.2939, warn_above=0.05, fail_above=0.1, seed=1 |
| K3:quad_area_accounting | pass | age_myr=30, quad=0.9994, expected=0.9994, warn_above=0.01, seed=1 |
| K4:hm_volume_km3_drift | fail | age_myr=400, quad=0.006688, lines=-0.2088, gap=0.2155, warn_above=0.05, fail_above=0.1, seed=2 |

## Mesh quality (quad)

| gate | status | worst case |
|---|---|---|
| M1:anisotropic_fraction | pass | age_myr=0, quad=0.000159, lines=0.000521, warn_above=0.01, seed=1 |
| M2:row_alignment | pass | age_myr=0, quad=0.08457, note=too few anisotropic nodes to have a direction, seed=1 |
| M3:thin_fraction | warn | age_myr=0, quad=0.000191, lines=0.000184, warn_above=0.000184, seed=1 |
| M4:aspect_gt_4 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M5:skew_gt_45 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M6:fragments | warn | age_myr=400, quad=3, lines=0, seed=1 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=1 |
| S2:stability_air_temperature_mean_c | pass | quad=0.1087, lines=0.1332, warn_above=0.3163, seed=1 |
| S2:stability_land_fraction | pass | quad=0.001533, lines=0.002011, warn_above=0.006021, seed=1 |
| S2:stability_ocean_temperature_mean_c | pass | quad=0.01002, lines=0.01056, warn_above=0.07112, seed=1 |
| S2:stability_precipitation_mean_mm | pass | quad=5.693, lines=6.101, warn_above=17.2, seed=1 |
| S2:stability_sea_level_m | pass | quad=6.065, lines=10.04, warn_above=21.08, seed=1 |

## Derived indexes: HEALPix vs KD-tree

| gate | status | worst case |
|---|---|---|
| X1:healpix_all | pass | age_myr=0, samples=200000, quad_same_node=0.5515, lines_same_node=0.5572, quad_distance_p95=1.195, lines_distance_p95=1.18, seed=1 |
| X1:healpix_antimeridian | pass | age_myr=0, samples=3333, quad_same_node=0.5581, lines_same_node=0.5467, quad_distance_p95=1.283, lines_distance_p95=1.216, seed=1 |
| X1:healpix_coast | pass | age_myr=0, samples=12454, quad_same_node=0.5499, lines_same_node=0.5449, quad_distance_p95=1.212, lines_distance_p95=1.22, seed=1 |
| X1:healpix_hole | insufficient | age_myr=0, samples=0, seed=1 |
| X1:healpix_plate_boundary | pass | age_myr=0, samples=17562, quad_same_node=0.5399, lines_same_node=0.5446, quad_distance_p95=1.257, lines_distance_p95=1.233, seed=1 |
| X1:healpix_pole | pass | age_myr=0, samples=3038, quad_same_node=0.525, lines_same_node=0.5319, quad_distance_p95=1.353, lines_distance_p95=1.359, seed=1 |

## Ensemble parity

| gate | status | worst case |
|---|---|---|
| P1:continental_hc_volume_km3 | pass | age_myr=30, seeds=5, lines_mean=4.15e+09, quad_mean=4.2e+09, lines_sigma=2.14e+09, delta=4.97e+07, fail_above=4.28e+09 |
| P1:elevation_p05 | fail | age_myr=120, seeds=5, lines_mean=-6296, quad_mean=-5842, lines_sigma=114.7, delta=454, fail_above=229.5 |
| P1:elevation_p50 | pass | age_myr=30, seeds=5, lines_mean=-5097, quad_mean=-5086, lines_sigma=64.91, delta=10.84, fail_above=129.8 |
| P1:elevation_p95 | pass | age_myr=30, seeds=5, lines_mean=500.4, quad_mean=513.9, lines_sigma=2493, delta=13.47, fail_above=4986 |
| P1:hc_volume_km3 | pass | age_myr=30, seeds=5, lines_mean=8.11e+09, quad_mean=8.1e+09, lines_sigma=1.87e+09, delta=-6.01e+06, fail_above=3.74e+09 |
| P1:land_fraction | pass | age_myr=30, seeds=5, lines_mean=0.1453, quad_mean=0.1422, lines_sigma=0.11, delta=-0.003096, fail_above=0.22 |
| P1:plates | pass | age_myr=30, seeds=5, lines_mean=18.4, quad_mean=18.6, lines_sigma=2.408, delta=0.2, fail_above=4.817 |
| P1:sea_level_m | fail | age_myr=240, seeds=5, lines_mean=-1696, quad_mean=-994.5, lines_sigma=191.8, delta=701.4, fail_above=383.5 |

## Checkpoints, seed 1

| seed 1 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 14 | 14 | 15 | 15 | 16 | 16 | 19 | 20 | 26 | 30 | 41 | 43 |
| nodes | 32639 | 31493 | 32659 | 31600 | 32467 | 32011 | 33168 | 32464 | 36582 | 33055 | 39503 | 32792 |
| land fraction | 0.2795 | 0.2793 | 0.3138 | 0.3118 | 0.3138 | 0.3082 | 0.3115 | 0.284 | 0.3078 | 0.2028 | 0.35 | 0.09436 |
| sea level (m) | 0 | 0 | -688.6 | -631.7 | -768.6 | -671.9 | -1032 | -773.2 | -1829 | -1320 | -2214 | -2041 |
| Hc volume (km³) | 1.15e+10 | 1.15e+10 | 1.1e+10 | 1.11e+10 | 1.08e+10 | 1.08e+10 | 1.03e+10 | 1.05e+10 | 9.24e+09 | 8.65e+09 | 7.89e+09 | 6.16e+09 |
| continental Hc (km³) | 7.65e+09 | 7.65e+09 | 7.55e+09 | 7.67e+09 | 7.51e+09 | 7.6e+09 | 7.39e+09 | 7.19e+09 | 6.93e+09 | 4.68e+09 | 6.12e+09 | 1.92e+09 |
| elevation p05/p50/p95 (m) | -4768/-4354/2990 | -4768/-4353/2984 | -5479/-5009/2994 | -5476/-4986/3186 | -5771/-5055/3198 | -5654/-4983/3315 | -6115/-5054/3339 | -5790/-4752/4325 | -6322/-4894/4487 | -5864/-4959/3603 | -6314/-4387/5269 | -5898/-5358/423.9 |
| uncovered | 0.02311 | 0.005985 | 0.02747 | 0.007185 | 0.03253 | 0.007435 | 0.04503 | 0.00896 | 0.05169 | 0.01152 | 0.06657 | 0.01488 |
| multiply covered | 0.01406 | 0.006475 | 0.00236 | 0.006545 | 0.006555 | 0.007545 | 0.0046 | 0.01131 | 0.003345 | 0.01301 | 0.00532 | 0.01253 |
| void | 0 | 0 | 0 | 0 | 0.001105 | 0 | 0.00139 | 0 | 0.00052 | 0 | 0.00202 | 0 |
| nodes inside other plate | 0.01618 | 0.009653 | 0.007134 | 0.003165 | 0.01266 | 0.003405 | 0.01058 | 0.008409 | 0.01268 | 0.008259 | 0.01663 | 0.002714 |
| stacked | 0 | 0 | 0.000919 | 0 | 0.003727 | 0 | 0.06024 | 0 | 0.2166 | 0 | 0.3195 | 0 |
| anisotropic | 0.000521 | 0.000159 | 0.002051 | 0.001709 | 0.003758 | 0.001906 | 0.005548 | 0.002464 | 0.02135 | 0.003751 | 0.07394 | 0.006343 |
| row alignment | 0.5571 | 0.08457 | -0.07381 | 0.02076 | -0.3011 | -0.122 | -0.197 | -0.2557 | 0.4583 | 0.1135 | 0.807 | -0.05098 |
| thin | 0.000184 | 0.000191 | 0.00199 | 0.002658 | 0.004004 | 0.003811 | 0.006301 | 0.004035 | 0.008228 | 0.006262 | 0.01691 | 0.008325 |
| one-node lines | 0.009498 | — | 0.07475 | — | 0.1395 | — | 0.1534 | — | 0.2181 | — | 0.2385 | — |
| quad one-cell-thin | — | 0.000381 | — | 0.006171 | — | 0.006935 | — | 0.00881 | — | 0.01228 | — | 0.01534 |
| air temperature mean (°C) | 4.49 | 4.481 | 7.557 | 6.926 | 7.946 | 7.701 | 7.864 | 5.861 | 5.198 | 1.827 | 6.724 | 9.582 |
| precipitation mean (mm) | 986.1 | 987 | 1012 | 1010 | 1013 | 1000 | 1038 | 1048 | 1082 | 1134 | 1121 | 1233 |
| river fraction of land | — | — | 0.1056 | 0.114 | 0.1079 | 0.1176 | 0.1055 | 0.1416 | 0.1335 | 0.1769 | 0.1309 | 0.1471 |
| HEALPix same node | 0.5572 | 0.5515 | 0.5549 | 0.5502 | 0.5564 | 0.5497 | 0.5495 | 0.5493 | 0.5323 | 0.5455 | 0.5211 | 0.5449 |

## Checkpoints, seed 2

| seed 2 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 18 | 18 | 21 | 21 | 21 | 20 | 29 | 26 | 38 | 36 | 65 | 43 |
| nodes | 32634 | 31334 | 33331 | 31544 | 33498 | 31870 | 34629 | 32528 | 37402 | 32527 | 40632 | 32529 |
| land fraction | 0.1406 | 0.1403 | 0.1511 | 0.1469 | 0.1488 | 0.1438 | 0.1542 | 0.1336 | 0.1632 | 0.0588 | 0.08428 | 0.01557 |
| sea level (m) | 0 | 0 | -772.4 | -607.4 | -969.6 | -677.7 | -1289 | -731.2 | -1797 | -1327 | -2378 | -1603 |
| Hc volume (km³) | 8.93e+09 | 8.91e+09 | 8.38e+09 | 8.39e+09 | 7.95e+09 | 8.06e+09 | 7.01e+09 | 7.54e+09 | 5.87e+09 | 5.73e+09 | 4.09e+09 | 4.72e+09 |
| continental Hc (km³) | 4.29e+09 | 4.28e+09 | 4.28e+09 | 4.36e+09 | 4.13e+09 | 4.16e+09 | 3.79e+09 | 3.52e+09 | 3.33e+09 | 1.51e+09 | 1.54e+09 | 1.83e+08 |
| elevation p05/p50/p95 (m) | -5179/-4437/1657 | -5174/-4437/1652 | -5838/-5104/1754 | -5532/-5086/1813 | -6242/-5156/1929 | -5757/-5122/1695 | -6392/-5251/1906 | -5870/-5107/1566 | -6474/-5419/2469 | -5929/-5391/-841 | -6439/-5487/-233.5 | -5937/-5471/-4337 |
| uncovered | 0.02514 | 0.00755 | 0.03544 | 0.00868 | 0.04558 | 0.00884 | 0.06081 | 0.01103 | 0.07643 | 0.01507 | 0.09685 | 0.01637 |
| multiply covered | 0.09836 | 0.007265 | 0.01212 | 0.00897 | 0.01518 | 0.00921 | 0.0031 | 0.01044 | 0.0037 | 0.01249 | 0.00428 | 0.01351 |
| void | 0 | 0 | 0 | 0 | 0.0007 | 0 | 0.002305 | 0 | 0.00423 | 0 | 0.001205 | 0 |
| nodes inside other plate | 0.08234 | 0.01098 | 0.02223 | 0.004216 | 0.02155 | 0.005993 | 0.01135 | 0.004827 | 0.01684 | 0.003351 | 0.01966 | 0.002613 |
| stacked | 0 | 0 | 0.03585 | 0 | 0.07153 | 0 | 0.1464 | 0 | 0.2826 | 0 | 0.3952 | 0 |
| anisotropic | 0.000276 | 0.000287 | 0.00384 | 0.002853 | 0.004418 | 0.002322 | 0.008374 | 0.003443 | 0.0208 | 0.005257 | 0.03672 | 0.006886 |
| row alignment | 0.4286 | -0.4423 | -0.3053 | 0.2016 | 0.09466 | 0.05847 | 0.2076 | 0.1916 | 0.3881 | 0.06302 | 0.4304 | -0.07137 |
| thin | 0.000123 | 9.57e-05 | 0.00366 | 0.003138 | 0.005105 | 0.003263 | 0.007162 | 0.004611 | 0.0143 | 0.007747 | 0.02067 | 0.01036 |
| one-node lines | 0.01008 | — | 0.09871 | — | 0.09152 | — | 0.1366 | — | 0.2212 | — | 0.2851 | — |
| quad one-cell-thin | — | 0.000734 | — | 0.007704 | — | 0.006934 | — | 0.01024 | — | 0.01534 | — | 0.01977 |
| air temperature mean (°C) | -3.378 | -3.453 | 0.8642 | -0.1418 | 2.335 | 1.092 | 3 | 1.992 | 2.285 | 12.78 | -1.182 | 17.66 |
| precipitation mean (mm) | 1132 | 1132 | 1168 | 1171 | 1179 | 1184 | 1194 | 1189 | 1184 | 1236 | 1247 | 1274 |
| river fraction of land | — | — | 0.1487 | 0.1555 | 0.1725 | 0.1686 | 0.1826 | 0.1857 | 0.1859 | 0.1125 | 0.2442 | 0.09483 |
| HEALPix same node | 0.5526 | 0.5486 | 0.5474 | 0.5518 | 0.5443 | 0.5522 | 0.5356 | 0.5471 | 0.5251 | 0.5492 | 0.5092 | 0.5483 |

## Checkpoints, seed 3

| seed 3 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 18 | 18 | 20 | 21 | 27 | 24 | 32 | 27 | 55 | 30 | 52 | 37 |
| nodes | 32632 | 31213 | 32601 | 31537 | 33172 | 32160 | 36309 | 32736 | 38155 | 33062 | 40117 | 33031 |
| land fraction | 0.05561 | 0.05551 | 0.07847 | 0.07413 | 0.084 | 0.07705 | 0.06338 | 0.0627 | 0.06129 | 0.03523 | 0.04183 | 0.01801 |
| sea level (m) | 0 | 0 | -515.8 | -483 | -801.6 | -538.9 | -1492 | -712.1 | -1827 | -935.3 | -2124 | -1135 |
| Hc volume (km³) | 7.2e+09 | 7.2e+09 | 6.88e+09 | 6.9e+09 | 6.51e+09 | 6.77e+09 | 5.15e+09 | 5.97e+09 | 4.28e+09 | 5.22e+09 | 3.52e+09 | 4.95e+09 |
| continental Hc (km³) | 2.59e+09 | 2.59e+09 | 2.76e+09 | 2.79e+09 | 2.65e+09 | 2.67e+09 | 1.81e+09 | 1.79e+09 | 1.12e+09 | 9.38e+08 | 4.87e+08 | 5e+08 |
| elevation p05/p50/p95 (m) | -5222/-4494/40.48 | -5222/-4493/41.41 | -5807/-5127/130.8 | -5531/-5118/40.26 | -6222/-5210/387.8 | -5733/-5162/201.1 | -6396/-5383/-797.2 | -5880/-5222/-272.9 | -6427/-5503/-1305 | -5942/-5378/-1977 | -6493/-5562/-2636 | -5945/-5449/-4151 |
| uncovered | 0.02706 | 0.00765 | 0.04388 | 0.008495 | 0.05449 | 0.00947 | 0.0631 | 0.01107 | 0.08456 | 0.01247 | 0.0874 | 0.01401 |
| multiply covered | 0.01106 | 0.00755 | 0.005605 | 0.01032 | 0.004775 | 0.009135 | 0.0103 | 0.01002 | 0.006305 | 0.01052 | 0.00532 | 0.01229 |
| void | 0 | 0 | 0.001785 | 0 | 0.000425 | 0 | 0.001245 | 0 | 0.000125 | 0 | 0.003235 | 0 |
| nodes inside other plate | 0.01492 | 0.01157 | 0.01635 | 0.005676 | 0.01501 | 0.00482 | 0.02754 | 0.004154 | 0.02372 | 0.00251 | 0.02066 | 0.004238 |
| stacked | 0 | 0 | 0.009172 | 0 | 0.05493 | 0 | 0.2116 | 0 | 0.289 | 0 | 0.39 | 0 |
| anisotropic | 0.000552 | 0.000192 | 0.003374 | 0.002029 | 0.006512 | 0.002425 | 0.01063 | 0.004002 | 0.02941 | 0.004355 | 0.03507 | 0.005086 |
| row alignment | 0.7907 | 0.3694 | -0.132 | 0.04715 | -0.000998 | -0.07017 | 0.01312 | -0.001982 | 0.3071 | -0.03216 | 0.4714 | 0.1924 |
| thin | 6.13e-05 | 0.000192 | 0.003282 | 0.003044 | 0.007536 | 0.004322 | 0.01168 | 0.004552 | 0.01929 | 0.005596 | 0.018 | 0.006963 |
| one-node lines | 0.006969 | — | 0.08553 | — | 0.1146 | — | 0.1708 | — | 0.2617 | — | 0.2717 | — |
| quad one-cell-thin | — | 0.000609 | — | 0.007578 | — | 0.008364 | — | 0.009958 | — | 0.01155 | — | 0.01505 |
| air temperature mean (°C) | -9.087 | -9.163 | -3.114 | -3.634 | -8.534 | -9.615 | -1.442 | -7.176 | 4.876 | 14.42 | -0.3359 | -8.753 |
| precipitation mean (mm) | 1201 | 1201 | 1238 | 1235 | 1247 | 1241 | 1263 | 1257 | 1257 | 1260 | 1266 | 1256 |
| river fraction of land | — | — | 0.1492 | 0.1604 | 0.2252 | 0.2234 | 0.2209 | 0.1554 | 0.2163 | 0.1661 | 0.1636 | 0.3857 |
| HEALPix same node | 0.5533 | 0.556 | 0.554 | 0.5538 | 0.5485 | 0.5556 | 0.5323 | 0.5496 | 0.5278 | 0.5496 | 0.5128 | 0.5504 |

## Checkpoints, seed 4

| seed 4 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 17 | 17 | 19 | 19 | 23 | 20 | 29 | 27 | 43 | 38 | 47 | 45 |
| nodes | 32636 | 31134 | 32564 | 31357 | 33037 | 31721 | 34487 | 32578 | 36962 | 33014 | 38992 | 32515 |
| land fraction | 0.1417 | 0.1416 | 0.1614 | 0.1585 | 0.1686 | 0.159 | 0.1707 | 0.1596 | 0.1433 | 0.09205 | 0.1246 | 0.03317 |
| sea level (m) | 0 | 0 | -565.6 | -536.8 | -733.6 | -528.9 | -1087 | -479.6 | -1645 | -1052 | -1969 | -1417 |
| Hc volume (km³) | 8.49e+09 | 8.49e+09 | 8.16e+09 | 8.12e+09 | 7.97e+09 | 7.97e+09 | 7.47e+09 | 7.77e+09 | 6.42e+09 | 6.56e+09 | 5.45e+09 | 5.15e+09 |
| continental Hc (km³) | 4.13e+09 | 4.13e+09 | 4.22e+09 | 4.24e+09 | 4.3e+09 | 4.31e+09 | 4.31e+09 | 4.28e+09 | 3.48e+09 | 2.42e+09 | 2.59e+09 | 6.65e+08 |
| elevation p05/p50/p95 (m) | -5222/-4427/958.8 | -5222/-4427/962.3 | -5637/-5066/1180 | -5506/-5055/1058 | -6012/-5111/1359 | -5719/-5063/1063 | -6275/-5172/1543 | -5824/-5011/1085 | -6340/-5251/2621 | -5913/-5273/887.1 | -6383/-5373/174.1 | -5930/-5384/-3053 |
| uncovered | 0.02816 | 0.00778 | 0.03963 | 0.008025 | 0.05167 | 0.00896 | 0.05948 | 0.01094 | 0.07301 | 0.01418 | 0.07369 | 0.0155 |
| multiply covered | 0.01667 | 0.00729 | 0.002055 | 0.00819 | 0.00239 | 0.00877 | 0.00492 | 0.01012 | 0.00188 | 0.01341 | 0.002985 | 0.01263 |
| void | 0 | 0 | 9.5e-05 | 0 | 0.002415 | 0 | 0.000155 | 0 | 0.002095 | 0 | 0.000505 | 0 |
| nodes inside other plate | 0.01906 | 0.01092 | 0.007677 | 0.002169 | 0.01014 | 0.002806 | 0.01673 | 0.00396 | 0.01147 | 0.0063 | 0.01739 | 0.00326 |
| stacked | 0 | 0 | 0.007554 | 0 | 0.04852 | 0 | 0.118 | 0 | 0.252 | 0 | 0.3224 | 0 |
| anisotropic | 0.00046 | 0.000225 | 0.003654 | 0.001531 | 0.007083 | 0.002427 | 0.01061 | 0.003469 | 0.01943 | 0.005028 | 0.02626 | 0.006797 |
| row alignment | 0.3941 | 0.03748 | -0.1217 | -0.005632 | -0.1282 | -0.06129 | -0.0828 | -0.1176 | 0.344 | -0.0153 | 0.3528 | 0.1198 |
| thin | 0.000184 | 3.21e-05 | 0.005036 | 0.002424 | 0.007386 | 0.004193 | 0.00983 | 0.004205 | 0.01253 | 0.006088 | 0.01726 | 0.009473 |
| one-node lines | 0.009368 | — | 0.09133 | — | 0.1425 | — | 0.1804 | — | 0.1959 | — | 0.2572 | — |
| quad one-cell-thin | — | 0.000418 | — | 0.00641 | — | 0.008449 | — | 0.01038 | — | 0.01372 | — | 0.01805 |
| air temperature mean (°C) | -5.967 | -5.829 | -7.355 | -7.355 | -8.496 | -8.058 | -5.337 | -4.634 | -1.376 | 0.7354 | 4.615 | 6.405 |
| precipitation mean (mm) | 1146 | 1147 | 1237 | 1255 | 1168 | 1184 | 1185 | 1180 | 1172 | 1231 | 1196 | 1244 |
| river fraction of land | — | — | 0.1574 | 0.1556 | 0.1546 | 0.1546 | 0.155 | 0.1608 | 0.1764 | 0.2308 | 0.1725 | 0.243 |
| HEALPix same node | 0.5553 | 0.5487 | 0.5518 | 0.5488 | 0.5519 | 0.5489 | 0.5436 | 0.5494 | 0.5343 | 0.547 | 0.5221 | 0.5498 |

## Checkpoints, seed 5

| seed 5 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 16 | 16 | 17 | 17 | 18 | 18 | 23 | 24 | 44 | 36 | 49 | 35 |
| nodes | 32643 | 31550 | 32728 | 31544 | 32701 | 31555 | 33678 | 32091 | 38149 | 33052 | 40169 | 33011 |
| land fraction | 0.01461 | 0.01469 | 0.02181 | 0.01988 | 0.02562 | 0.02556 | 0.03194 | 0.03122 | 0.0558 | 0.02867 | 0.0821 | 0.01037 |
| sea level (m) | 0 | 0 | -451.6 | -432.5 | -487.1 | -405.5 | -700.3 | -354.8 | -1381 | -337.9 | -1606 | -660.4 |
| Hc volume (km³) | 6.19e+09 | 6.19e+09 | 6.11e+09 | 6.03e+09 | 6e+09 | 6.04e+09 | 5.62e+09 | 5.9e+09 | 4.48e+09 | 5.5e+09 | 4.09e+09 | 4.97e+09 |
| continental Hc (km³) | 1.86e+09 | 1.86e+09 | 1.96e+09 | 1.95e+09 | 2.02e+09 | 2.02e+09 | 1.96e+09 | 1.97e+09 | 1.64e+09 | 1.32e+09 | 1.39e+09 | 5.5e+08 |
| elevation p05/p50/p95 (m) | -5222/-4570/-4197 | -5222/-4570/-4198 | -5500/-5182/-3558 | -5511/-5187/-3527 | -5920/-5185/-3110 | -5706/-5207/-2836 | -6299/-5229/-2506 | -5844/-5221/-1872 | -6434/-5411/-1203 | -5897/-5223/-1824 | -6469/-5465/-494.2 | -5941/-5432/-3757 |
| uncovered | 0.02357 | 0.006985 | 0.03184 | 0.00752 | 0.03644 | 0.00841 | 0.05019 | 0.01027 | 0.07593 | 0.01406 | 0.07455 | 0.01383 |
| multiply covered | 0.03266 | 0.00718 | 0.00541 | 0.00821 | 0.00485 | 0.009195 | 0.005325 | 0.0104 | 0.005855 | 0.01155 | 0.00614 | 0.01201 |
| void | 0 | 0 | 0 | 0 | 0 | 0 | 0.0008 | 0 | 0.003855 | 0 | 0.001225 | 0 |
| nodes inside other plate | 0.0367 | 0.0109 | 0.0125 | 0.004216 | 0.01156 | 0.004849 | 0.01416 | 0.00402 | 0.02045 | 0.005113 | 0.02305 | 0.003241 |
| stacked | 0 | 0 | 0.001161 | 0 | 0.009939 | 0 | 0.08278 | 0 | 0.3035 | 0 | 0.3602 | 0 |
| anisotropic | 0.000429 | 9.51e-05 | 0.001681 | 0.002061 | 0.003119 | 0.001997 | 0.004602 | 0.004768 | 0.02197 | 0.005264 | 0.03164 | 0.005089 |
| row alignment | 0.3974 | -0.3597 | 0.08567 | 0.1791 | -0.3858 | -0.04844 | 0.01516 | 0.16 | 0.4266 | -0.07338 | 0.404 | 0.001211 |
| thin | 0.000123 | 6.34e-05 | 0.001558 | 0.003075 | 0.003089 | 0.00469 | 0.004691 | 0.006887 | 0.01371 | 0.005355 | 0.01387 | 0.007997 |
| one-node lines | 0.007371 | — | 0.05065 | — | 0.1036 | — | 0.106 | — | 0.2183 | — | 0.2729 | — |
| quad one-cell-thin | — | 0.000571 | — | 0.006689 | — | 0.007954 | — | 0.01134 | — | 0.01349 | — | 0.01387 |
| air temperature mean (°C) | -2.289 | -2.44 | -6.049 | -6.007 | -7.363 | -6.969 | -9.886 | -9.196 | -10.09 | 4.799 | -2.556 | -6.939 |
| precipitation mean (mm) | 1254 | 1254 | 1252 | 1256 | 1263 | 1262 | 1263 | 1260 | 1261 | 1262 | 1251 | 1265 |
| river fraction of land | — | — | 0.1357 | 0.1289 | 0.1629 | 0.1342 | 0.1667 | 0.2736 | 0.1962 | 0.1991 | 0.2063 | 0.2607 |
| HEALPix same node | 0.5552 | 0.5447 | 0.5529 | 0.5491 | 0.5526 | 0.5494 | 0.5464 | 0.5524 | 0.5257 | 0.5515 | 0.5186 | 0.5459 |

## Performance (mean seconds per step)

| s/step | seed 1 quad | seed 2 quad | seed 3 quad | seed 4 quad | seed 5 quad |
|---|---|---|---|---|---|
| step_total | 2.382 | 2.462 | 2.389 | 2.509 | 2.36 |
| deform_topology | 0.7525 | 0.8169 | 0.8107 | 0.8414 | 0.7829 |
| climate_erosion_hydrology | 0.5687 | 0.5551 | 0.5432 | 0.5448 | 0.5413 |
| deform | 0.3693 | 0.3868 | 0.3824 | 0.3921 | 0.3696 |
| faults | 0.6551 | 0.7217 | 0.6879 | 0.7441 | 0.6952 |
| fluid_dynamics | 1.72e-06 | 1.62e-06 | 1.83e-06 | 1.69e-06 | 1.74e-06 |
| gap_fill | 0.3172 | 0.3605 | 0.3595 | 0.3757 | 0.3443 |
| magma_transport | 0.1268 | 0.0779 | 0.05335 | 0.0845 | 0.04276 |
| overlap_tracking | 0.03356 | 0.03533 | 0.03545 | 0.03921 | 0.03595 |
| record_stats | 0.003701 | 0.003234 | 0.00326 | 0.003285 | 0.00321 |
| resource_formation | 0.009652 | 0.009371 | 0.009259 | 0.009029 | 0.008962 |
| sea_level | 0.04776 | 0.05489 | 0.05863 | 0.05435 | 0.05812 |
| shift | 0.2072 | 0.2114 | 0.2115 | 0.2163 | 0.2168 |
| stranded_basins | 1.42e-05 | 9.96e-06 | 8.85e-06 | 9.9e-06 | 8.04e-06 |
| topology | 0.03245 | 0.03431 | 0.03335 | 0.03437 | 0.03303 |
| volcanism | 0.009317 | 0.01018 | 0.01007 | 0.01004 | 0.009699 |

Index builds per step (mean):

- seed 1 quad: plate_node_kdtree 78.6, plate_outline 86.7, plate_outline_kdtree 80.8, quad_adjacency 43.2, quad_boundary_loops 58.9, world_node_kdtree 1.0
- seed 2 quad: plate_node_kdtree 90.7, plate_outline 101.7, plate_outline_kdtree 92.9, quad_adjacency 50.9, quad_boundary_loops 70.0, world_node_kdtree 1.0
- seed 3 quad: plate_node_kdtree 87.4, plate_outline 98.6, plate_outline_kdtree 90.4, quad_adjacency 49.2, quad_boundary_loops 67.7, world_node_kdtree 1.0
- seed 4 quad: plate_node_kdtree 93.1, plate_outline 105.7, plate_outline_kdtree 97.6, quad_adjacency 52.3, quad_boundary_loops 72.3, world_node_kdtree 1.0
- seed 5 quad: plate_node_kdtree 81.1, plate_outline 92.1, plate_outline_kdtree 85.3, quad_adjacency 44.7, quad_boundary_loops 62.5, world_node_kdtree 1.0

