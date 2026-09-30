# Line vs quad parity report

**Verdict: FAIL** — configs long; seeds 1, 2, 3, 4, 5.

Reproduce:

    bin/debug/surface_parity.py paired --preset long --seeds 1,2,3,4,5 --jobs 5 --out ../analysis/issue249-campaign/long

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
| H11:field_caps | warn (line baseline: 5 finding(s)) | violations=200, surface=quad, seed=1; first: {'detail': 'mantle_lithosphere_thickness_m: 3 values outside [2000, 240000] (range 1823.6..240000)', 'kind': 'field_bounds', 'plate_id': 2, 'step': 30, 'surface': 'quad'} |

## Coverage and overlap (quad vs paired line run)

| gate | status | worst case |
|---|---|---|
| C1:uncovered | pass | age_myr=0, quad=0.005985, lines=0.02311, fail_above=0.02811, seed=1 |
| C2:void | pass | age_myr=0, quad=0, lines=0, fail_above=0.0005, seed=1 |
| C3:multiply_covered | fail | age_myr=240, quad=0.01243, lines=0.003345, fail_above=0.008345, seed=1 |
| C4:nodes_inside_other_plate | pass | age_myr=0, quad=0.009653, lines=0.01618, fail_above=0.02118, seed=1 |

## Conservation (actual node/cell areas)

| gate | status | worst case |
|---|---|---|
| K1:hc_volume_km3_drift | fail | age_myr=400, quad=-0.4265, lines=-0.5417, gap=0.1152, warn_above=0.05, fail_above=0.1, seed=2 |
| K2:continental_hc_volume_km3_drift | fail | age_myr=240, quad=-0.2612, lines=-0.09438, gap=0.1668, warn_above=0.05, fail_above=0.1, seed=1 |
| K3:quad_area_accounting | pass | age_myr=30, quad=0.9898, expected=0.9897, warn_above=0.01, seed=1 |
| K4:hm_volume_km3_drift | fail | age_myr=400, quad=0.00123, lines=-0.2088, gap=0.21, warn_above=0.05, fail_above=0.1, seed=2 |

## Mesh quality (quad)

| gate | status | worst case |
|---|---|---|
| M1:anisotropic_fraction | pass | age_myr=0, quad=0.000159, lines=0.000521, warn_above=0.01, seed=1 |
| M2:row_alignment | pass | age_myr=0, quad=0.08457, note=too few anisotropic nodes to have a direction, seed=1 |
| M3:thin_fraction | warn | age_myr=0, quad=0.000191, lines=0.000184, warn_above=0.000184, seed=1 |
| M4:aspect_gt_4 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M5:skew_gt_45 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M6:fragments | warn | age_myr=240, quad=3, lines=0, seed=1 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=1 |
| S2:stability_air_temperature_mean_c | warn | quad=0.7868, lines=0.2848, warn_above=0.6196, seed=3 |
| S2:stability_land_fraction | warn | quad=0.004247, lines=0.000907, warn_above=0.003813, seed=5 |
| S2:stability_ocean_temperature_mean_c | pass | quad=0.01579, lines=0.01056, warn_above=0.07112, seed=1 |
| S2:stability_precipitation_mean_mm | pass | quad=5.78, lines=6.101, warn_above=17.2, seed=1 |
| S2:stability_sea_level_m | warn | quad=58.21, lines=10.04, warn_above=21.08, seed=1 |

## Derived indexes: HEALPix vs KD-tree

| gate | status | worst case |
|---|---|---|
| X1:healpix_all | pass | age_myr=0, samples=200000, quad_same_node=0.5515, lines_same_node=0.5572, quad_distance_p95=1.195, lines_distance_p95=1.18, seed=1 |
| X1:healpix_antimeridian | pass | age_myr=0, samples=3333, quad_same_node=0.5581, lines_same_node=0.5467, quad_distance_p95=1.283, lines_distance_p95=1.216, seed=1 |
| X1:healpix_coast | pass | age_myr=0, samples=12454, quad_same_node=0.5499, lines_same_node=0.5449, quad_distance_p95=1.212, lines_distance_p95=1.22, seed=1 |
| X1:healpix_hole | warn | age_myr=60, samples=131, quad_same_node=0.5038, lines_same_node=0.5613, quad_distance_p95=1.267, lines_distance_p95=1.261, seed=2 |
| X1:healpix_plate_boundary | pass | age_myr=0, samples=17562, quad_same_node=0.5399, lines_same_node=0.5446, quad_distance_p95=1.257, lines_distance_p95=1.233, seed=1 |
| X1:healpix_pole | pass | age_myr=0, samples=3038, quad_same_node=0.525, lines_same_node=0.5319, quad_distance_p95=1.353, lines_distance_p95=1.359, seed=1 |

## Performance

| gate | status | worst case |
|---|---|---|
| R1:deform_topology_s_per_step | pass | quad=0.7553, lines=2.074, ratio=0.3642, warn_above=1.25, seed=1 |
| R2:step_total_s_per_step | pass | quad=2.976, lines=4.308, ratio=0.6907, warn_above=1.25, seed=1 |

## Ensemble parity

| gate | status | worst case |
|---|---|---|
| P1:continental_hc_volume_km3 | pass | age_myr=30, seeds=5, lines_mean=4.15e+09, quad_mean=4.23e+09, lines_sigma=2.14e+09, delta=7.41e+07, fail_above=4.28e+09 |
| P1:elevation_p05 | fail | age_myr=120, seeds=5, lines_mean=-6296, quad_mean=-5870, lines_sigma=114.7, delta=425.5, fail_above=229.5 |
| P1:elevation_p50 | pass | age_myr=30, seeds=5, lines_mean=-5097, quad_mean=-5088, lines_sigma=64.91, delta=8.852, fail_above=129.8 |
| P1:elevation_p95 | pass | age_myr=30, seeds=5, lines_mean=500.4, quad_mean=560.6, lines_sigma=2493, delta=60.23, fail_above=4986 |
| P1:hc_volume_km3 | pass | age_myr=30, seeds=5, lines_mean=8.11e+09, quad_mean=8.09e+09, lines_sigma=1.87e+09, delta=-1.99e+07, fail_above=3.74e+09 |
| P1:land_fraction | pass | age_myr=30, seeds=5, lines_mean=0.1453, quad_mean=0.1441, lines_sigma=0.11, delta=-0.001261, fail_above=0.22 |
| P1:plates | pass | age_myr=30, seeds=5, lines_mean=18.4, quad_mean=19, lines_sigma=2.408, delta=0.6, fail_above=4.817 |
| P1:sea_level_m | fail | age_myr=240, seeds=5, lines_mean=-1696, quad_mean=-1094, lines_sigma=191.8, delta=601.9, fail_above=383.5 |

## Checkpoints, seed 1

| seed 1 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 14 | 14 | 15 | 15 | 16 | 16 | 19 | 19 | 26 | 27 | 41 | 32 |
| nodes | 32639 | 31493 | 32659 | 31306 | 32467 | 32056 | 33168 | 32668 | 36582 | 32734 | 39503 | 32907 |
| land fraction | 0.2795 | 0.2793 | 0.3138 | 0.3144 | 0.3138 | 0.3154 | 0.3115 | 0.3017 | 0.3078 | 0.2366 | 0.35 | 0.1635 |
| sea level (m) | 0 | 0 | -688.6 | -569.8 | -768.6 | -768.4 | -1032 | -1018 | -1829 | -1438 | -2214 | -1845 |
| Hc volume (km³) | 1.15e+10 | 1.15e+10 | 1.1e+10 | 1.11e+10 | 1.08e+10 | 1.09e+10 | 1.03e+10 | 1.05e+10 | 9.24e+09 | 9.31e+09 | 7.89e+09 | 7.57e+09 |
| continental Hc (km³) | 7.65e+09 | 7.65e+09 | 7.55e+09 | 7.69e+09 | 7.51e+09 | 7.69e+09 | 7.39e+09 | 7.35e+09 | 6.93e+09 | 5.65e+09 | 6.12e+09 | 3.03e+09 |
| elevation p05/p50/p95 (m) | -4768/-4354/2990 | -4768/-4353/2984 | -5479/-5009/2994 | -5473/-4976/3228 | -5771/-5055/3198 | -5667/-4971/3490 | -6115/-5054/3339 | -5795/-4851/4441 | -6322/-4894/4487 | -5935/-5059/4354 | -6314/-4387/5269 | -5939/-5207/2373 |
| uncovered | 0.02311 | 0.005985 | 0.02747 | 0.01451 | 0.03253 | 0.0067 | 0.04503 | 0.00751 | 0.05169 | 0.00987 | 0.06657 | 0.01029 |
| multiply covered | 0.01406 | 0.006475 | 0.00236 | 0.00421 | 0.006555 | 0.008245 | 0.0046 | 0.00954 | 0.003345 | 0.01243 | 0.00532 | 0.01379 |
| void | 0 | 0 | 0 | 0 | 0.001105 | 0 | 0.00139 | 0 | 0.00052 | 0 | 0.00202 | 0 |
| nodes inside other plate | 0.01618 | 0.009653 | 0.007134 | 0.003194 | 0.01266 | 0.004399 | 0.01058 | 0.006153 | 0.01268 | 0.007943 | 0.01663 | 0.008205 |
| stacked | 0 | 0 | 0.000919 | 0 | 0.003727 | 0 | 0.06024 | 0 | 0.2166 | 0 | 0.3195 | 0 |
| anisotropic | 0.000521 | 0.000159 | 0.002051 | 0.00099 | 0.003758 | 0.001809 | 0.005548 | 0.002296 | 0.02135 | 0.003452 | 0.07394 | 0.004163 |
| row alignment | 0.5571 | 0.08457 | -0.07381 | -0.3983 | -0.3011 | 0.009192 | -0.197 | 0.09361 | 0.4583 | -0.08133 | 0.807 | 0.04768 |
| thin | 0.000184 | 0.000191 | 0.00199 | 0.00214 | 0.004004 | 0.003369 | 0.006301 | 0.003122 | 0.008228 | 0.00388 | 0.01691 | 0.005257 |
| one-node lines | 0.009498 | — | 0.07475 | — | 0.1395 | — | 0.1534 | — | 0.2181 | — | 0.2385 | — |
| quad one-cell-thin | — | 0.000381 | — | 0.005718 | — | 0.006333 | — | 0.006428 | — | 0.008157 | — | 0.009299 |
| air temperature mean (°C) | 4.49 | 4.481 | 7.557 | 7.111 | 7.946 | 7.649 | 7.864 | 6.205 | 5.198 | 4.099 | 6.724 | 6.28 |
| precipitation mean (mm) | 986.1 | 987 | 1012 | 1009 | 1013 | 997.6 | 1038 | 1048 | 1082 | 1102 | 1121 | 1185 |
| river fraction of land | — | — | 0.1056 | 0.1183 | 0.1079 | 0.1182 | 0.1055 | 0.1424 | 0.1335 | 0.1458 | 0.1309 | 0.1621 |
| HEALPix same node | 0.5572 | 0.5515 | 0.5549 | 0.5504 | 0.5564 | 0.5503 | 0.5495 | 0.5475 | 0.5323 | 0.5466 | 0.5211 | 0.549 |

## Checkpoints, seed 2

| seed 2 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 18 | 18 | 21 | 20 | 21 | 22 | 29 | 29 | 38 | 39 | 65 | 49 |
| nodes | 32634 | 31334 | 33331 | 31095 | 33498 | 31900 | 34629 | 32490 | 37402 | 33048 | 40632 | 32962 |
| land fraction | 0.1406 | 0.1403 | 0.1511 | 0.1492 | 0.1488 | 0.1465 | 0.1542 | 0.1368 | 0.1632 | 0.09322 | 0.08428 | 0.03414 |
| sea level (m) | 0 | 0 | -772.4 | -561.5 | -969.6 | -788 | -1289 | -971.6 | -1797 | -1236 | -2378 | -1655 |
| Hc volume (km³) | 8.93e+09 | 8.91e+09 | 8.38e+09 | 8.38e+09 | 7.95e+09 | 8.15e+09 | 7.01e+09 | 7.57e+09 | 5.87e+09 | 6.36e+09 | 4.09e+09 | 5.11e+09 |
| continental Hc (km³) | 4.29e+09 | 4.28e+09 | 4.28e+09 | 4.4e+09 | 4.13e+09 | 4.22e+09 | 3.79e+09 | 3.63e+09 | 3.33e+09 | 1.93e+09 | 1.54e+09 | 6.25e+08 |
| elevation p05/p50/p95 (m) | -5179/-4437/1657 | -5174/-4437/1652 | -5838/-5104/1754 | -5533/-5089/1881 | -6242/-5156/1929 | -5784/-5133/1906 | -6392/-5251/1906 | -5905/-5172/1610 | -6474/-5419/2469 | -5946/-5270/955.8 | -6439/-5487/-233.5 | -5973/-5417/-3569 |
| uncovered | 0.02514 | 0.00755 | 0.03544 | 0.01838 | 0.04558 | 0.00809 | 0.06081 | 0.01039 | 0.07643 | 0.01446 | 0.09685 | 0.01423 |
| multiply covered | 0.09836 | 0.007265 | 0.01212 | 0.00445 | 0.01518 | 0.008895 | 0.0031 | 0.01171 | 0.0037 | 0.01389 | 0.00428 | 0.01363 |
| void | 0 | 0 | 0 | 0 | 0.0007 | 0 | 0.002305 | 0 | 0.00423 | 0 | 0.001205 | 0 |
| nodes inside other plate | 0.08234 | 0.01098 | 0.02223 | 0.00164 | 0.02155 | 0.00442 | 0.01135 | 0.005633 | 0.01684 | 0.005265 | 0.01966 | 0.003519 |
| stacked | 0 | 0 | 0.03585 | 0 | 0.07153 | 0 | 0.1464 | 0 | 0.2826 | 0 | 0.3952 | 0 |
| anisotropic | 0.000276 | 0.000287 | 0.00384 | 0.001383 | 0.004418 | 0.0021 | 0.008374 | 0.002709 | 0.0208 | 0.004539 | 0.03672 | 0.005977 |
| row alignment | 0.4286 | -0.4423 | -0.3053 | -0.05453 | 0.09466 | -0.2125 | 0.2076 | 0.04508 | 0.3881 | -0.04606 | 0.4304 | -0.02063 |
| thin | 0.000123 | 9.57e-05 | 0.00366 | 0.002541 | 0.005105 | 0.002226 | 0.007162 | 0.003447 | 0.0143 | 0.005719 | 0.02067 | 0.006674 |
| one-node lines | 0.01008 | — | 0.09871 | — | 0.09152 | — | 0.1366 | — | 0.2212 | — | 0.2851 | — |
| quad one-cell-thin | — | 0.000734 | — | 0.005821 | — | 0.006489 | — | 0.008249 | — | 0.01292 | — | 0.01474 |
| air temperature mean (°C) | -3.378 | -3.453 | 0.8642 | 0.1569 | 2.335 | 1.181 | 3 | 2.307 | 2.285 | 4.21 | -1.182 | 16.82 |
| precipitation mean (mm) | 1132 | 1132 | 1168 | 1164 | 1179 | 1175 | 1194 | 1190 | 1184 | 1212 | 1247 | 1265 |
| river fraction of land | — | — | 0.1487 | 0.1535 | 0.1725 | 0.1965 | 0.1826 | 0.2088 | 0.1859 | 0.1751 | 0.2442 | 0.1069 |
| HEALPix same node | 0.5526 | 0.5486 | 0.5474 | 0.55 | 0.5443 | 0.5525 | 0.5356 | 0.5467 | 0.5251 | 0.5455 | 0.5092 | 0.5496 |

## Checkpoints, seed 3

| seed 3 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 18 | 18 | 20 | 22 | 27 | 25 | 32 | 29 | 55 | 39 | 52 | 45 |
| nodes | 32632 | 31213 | 32601 | 31054 | 33172 | 32294 | 36309 | 32723 | 38155 | 33072 | 40117 | 32698 |
| land fraction | 0.05561 | 0.05551 | 0.07847 | 0.07717 | 0.084 | 0.08386 | 0.06338 | 0.05627 | 0.06129 | 0.04273 | 0.04183 | 0.0336 |
| sea level (m) | 0 | 0 | -515.8 | -464.6 | -801.6 | -748.1 | -1492 | -1008 | -1827 | -1101 | -2124 | -1191 |
| Hc volume (km³) | 7.2e+09 | 7.2e+09 | 6.88e+09 | 6.88e+09 | 6.51e+09 | 6.81e+09 | 5.15e+09 | 5.86e+09 | 4.28e+09 | 5.37e+09 | 3.52e+09 | 5.19e+09 |
| continental Hc (km³) | 2.59e+09 | 2.59e+09 | 2.76e+09 | 2.83e+09 | 2.65e+09 | 2.74e+09 | 1.81e+09 | 1.55e+09 | 1.12e+09 | 1.08e+09 | 4.87e+08 | 7.81e+08 |
| elevation p05/p50/p95 (m) | -5222/-4494/40.48 | -5222/-4493/41.41 | -5807/-5127/130.8 | -5585/-5129/84.96 | -6222/-5210/387.8 | -5801/-5185/287.6 | -6396/-5383/-797.2 | -5917/-5269/-840.1 | -6427/-5503/-1305 | -5949/-5325/-1418 | -6493/-5562/-2636 | -5967/-5398/-2616 |
| uncovered | 0.02706 | 0.00765 | 0.04388 | 0.02121 | 0.05449 | 0.009015 | 0.0631 | 0.01081 | 0.08456 | 0.01242 | 0.0874 | 0.01363 |
| multiply covered | 0.01106 | 0.00755 | 0.005605 | 0.00611 | 0.004775 | 0.01065 | 0.0103 | 0.01213 | 0.006305 | 0.0119 | 0.00532 | 0.01417 |
| void | 0 | 0 | 0.001785 | 0 | 0.000425 | 0 | 0.001245 | 0 | 0.000125 | 0 | 0.003235 | 0 |
| nodes inside other plate | 0.01492 | 0.01157 | 0.01635 | 0.004669 | 0.01501 | 0.00418 | 0.02754 | 0.005898 | 0.02372 | 0.003447 | 0.02066 | 0.003976 |
| stacked | 0 | 0 | 0.009172 | 0 | 0.05493 | 0 | 0.2116 | 0 | 0.289 | 0 | 0.39 | 0 |
| anisotropic | 0.000552 | 0.000192 | 0.003374 | 0.001739 | 0.006512 | 0.002106 | 0.01063 | 0.002781 | 0.02941 | 0.004142 | 0.03507 | 0.00526 |
| row alignment | 0.7907 | 0.3694 | -0.132 | -0.08311 | -0.000998 | -0.05464 | 0.01312 | -0.04059 | 0.3071 | 0.02868 | 0.4714 | 0.175 |
| thin | 6.13e-05 | 0.000192 | 0.003282 | 0.002351 | 0.007536 | 0.003344 | 0.01168 | 0.004462 | 0.01929 | 0.005745 | 0.018 | 0.007432 |
| one-node lines | 0.006969 | — | 0.08553 | — | 0.1146 | — | 0.1708 | — | 0.2617 | — | 0.2717 | — |
| quad one-cell-thin | — | 0.000609 | — | 0.005507 | — | 0.008082 | — | 0.009657 | — | 0.01149 | — | 0.01483 |
| air temperature mean (°C) | -9.087 | -9.163 | -3.114 | -4.145 | -8.534 | -10.08 | -1.442 | -4.943 | 4.876 | 2.165 | -0.3359 | 5.099 |
| precipitation mean (mm) | 1201 | 1201 | 1238 | 1234 | 1247 | 1236 | 1263 | 1263 | 1257 | 1260 | 1266 | 1252 |
| river fraction of land | — | — | 0.1492 | 0.1655 | 0.2252 | 0.2232 | 0.2209 | 0.1886 | 0.2163 | 0.1975 | 0.1636 | 0.1341 |
| HEALPix same node | 0.5533 | 0.556 | 0.554 | 0.5523 | 0.5485 | 0.5538 | 0.5323 | 0.5486 | 0.5278 | 0.5525 | 0.5128 | 0.5502 |

## Checkpoints, seed 4

| seed 4 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 17 | 17 | 19 | 21 | 23 | 23 | 29 | 32 | 43 | 36 | 47 | 45 |
| nodes | 32636 | 31134 | 32564 | 30883 | 33037 | 31546 | 34487 | 32209 | 36962 | 32625 | 38992 | 32655 |
| land fraction | 0.1417 | 0.1416 | 0.1614 | 0.1597 | 0.1686 | 0.1627 | 0.1707 | 0.1586 | 0.1433 | 0.116 | 0.1246 | 0.05123 |
| sea level (m) | 0 | 0 | -565.6 | -489.5 | -733.6 | -592.2 | -1087 | -732.1 | -1645 | -1175 | -1969 | -1534 |
| Hc volume (km³) | 8.49e+09 | 8.49e+09 | 8.16e+09 | 8.1e+09 | 7.97e+09 | 8.06e+09 | 7.47e+09 | 7.69e+09 | 6.42e+09 | 6.98e+09 | 5.45e+09 | 5.54e+09 |
| continental Hc (km³) | 4.13e+09 | 4.13e+09 | 4.22e+09 | 4.27e+09 | 4.3e+09 | 4.41e+09 | 4.31e+09 | 4.25e+09 | 3.48e+09 | 3.04e+09 | 2.59e+09 | 1.03e+09 |
| elevation p05/p50/p95 (m) | -5222/-4427/958.8 | -5222/-4427/962.3 | -5637/-5066/1180 | -5504/-5059/1109 | -6012/-5111/1359 | -5725/-5071/1145 | -6275/-5172/1543 | -5867/-5030/1249 | -6340/-5251/2621 | -5925/-5231/1238 | -6383/-5373/174.1 | -5967/-5359/-1477 |
| uncovered | 0.02816 | 0.00778 | 0.03963 | 0.02008 | 0.05167 | 0.009595 | 0.05948 | 0.01091 | 0.07301 | 0.01247 | 0.07369 | 0.01291 |
| multiply covered | 0.01667 | 0.00729 | 0.002055 | 0.00446 | 0.00239 | 0.008715 | 0.00492 | 0.01146 | 0.00188 | 0.0129 | 0.002985 | 0.01398 |
| void | 0 | 0 | 9.5e-05 | 0 | 0.002415 | 0 | 0.000155 | 0 | 0.002095 | 0 | 0.000505 | 0 |
| nodes inside other plate | 0.01906 | 0.01092 | 0.007677 | 0.001425 | 0.01014 | 0.003043 | 0.01673 | 0.00562 | 0.01147 | 0.007356 | 0.01739 | 0.005757 |
| stacked | 0 | 0 | 0.007554 | 0 | 0.04852 | 0 | 0.118 | 0 | 0.252 | 0 | 0.3224 | 0 |
| anisotropic | 0.00046 | 0.000225 | 0.003654 | 0.001878 | 0.007083 | 0.002219 | 0.01061 | 0.003167 | 0.01943 | 0.003556 | 0.02626 | 0.005941 |
| row alignment | 0.3941 | 0.03748 | -0.1217 | 0.02713 | -0.1282 | 0.02388 | -0.0828 | 0.0332 | 0.344 | -0.000946 | 0.3528 | 0.08282 |
| thin | 0.000184 | 3.21e-05 | 0.005036 | 0.002331 | 0.007386 | 0.002346 | 0.00983 | 0.003974 | 0.01253 | 0.004506 | 0.01726 | 0.008146 |
| one-node lines | 0.009368 | — | 0.09133 | — | 0.1425 | — | 0.1804 | — | 0.1959 | — | 0.2572 | — |
| quad one-cell-thin | — | 0.000418 | — | 0.00599 | — | 0.007798 | — | 0.009035 | — | 0.01192 | — | 0.0136 |
| air temperature mean (°C) | -5.967 | -5.829 | -7.355 | -7.412 | -8.496 | -7.864 | -5.337 | 0.5155 | -1.376 | 2.287 | 4.615 | 3.094 |
| precipitation mean (mm) | 1146 | 1147 | 1237 | 1252 | 1168 | 1175 | 1185 | 1184 | 1172 | 1198 | 1196 | 1234 |
| river fraction of land | — | — | 0.1574 | 0.1588 | 0.1546 | 0.1586 | 0.155 | 0.1756 | 0.1764 | 0.2212 | 0.1725 | 0.2458 |
| HEALPix same node | 0.5553 | 0.5487 | 0.5518 | 0.5491 | 0.5519 | 0.5481 | 0.5436 | 0.5492 | 0.5343 | 0.5507 | 0.5221 | 0.5506 |

## Checkpoints, seed 5

| seed 5 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 16 | 16 | 17 | 17 | 18 | 18 | 23 | 22 | 44 | 31 | 49 | 39 |
| nodes | 32643 | 31550 | 32728 | 31182 | 32701 | 31608 | 33678 | 32445 | 38149 | 32664 | 40169 | 32676 |
| land fraction | 0.01461 | 0.01469 | 0.02181 | 0.01979 | 0.02562 | 0.02655 | 0.03194 | 0.03552 | 0.0558 | 0.03988 | 0.0821 | 0.01725 |
| sea level (m) | 0 | 0 | -451.6 | -388.2 | -487.1 | -435.5 | -700.3 | -511.7 | -1381 | -518.1 | -1606 | -751.8 |
| Hc volume (km³) | 6.19e+09 | 6.19e+09 | 6.11e+09 | 6e+09 | 6e+09 | 6.05e+09 | 5.62e+09 | 5.93e+09 | 4.48e+09 | 5.68e+09 | 4.09e+09 | 5.25e+09 |
| continental Hc (km³) | 1.86e+09 | 1.86e+09 | 1.96e+09 | 1.95e+09 | 2.02e+09 | 2.04e+09 | 1.96e+09 | 2.06e+09 | 1.64e+09 | 1.58e+09 | 1.39e+09 | 8.87e+08 |
| elevation p05/p50/p95 (m) | -5222/-4570/-4197 | -5222/-4570/-4198 | -5500/-5182/-3558 | -5517/-5188/-3499 | -5920/-5185/-3110 | -5751/-5211/-2885 | -6299/-5229/-2506 | -5866/-5221/-1873 | -6434/-5411/-1203 | -5930/-5222/-1404 | -6469/-5465/-494.2 | -5973/-5369/-3155 |
| uncovered | 0.02357 | 0.006985 | 0.03184 | 0.01479 | 0.03644 | 0.00756 | 0.05019 | 0.0094 | 0.07593 | 0.01085 | 0.07455 | 0.01366 |
| multiply covered | 0.03266 | 0.00718 | 0.00541 | 0.00455 | 0.00485 | 0.00871 | 0.005325 | 0.0104 | 0.005855 | 0.01212 | 0.00614 | 0.01296 |
| void | 0 | 0 | 0 | 0 | 0 | 0 | 0.0008 | 0 | 0.003855 | 0 | 0.001225 | 0 |
| nodes inside other plate | 0.0367 | 0.0109 | 0.0125 | 0.003111 | 0.01156 | 0.005347 | 0.01416 | 0.005425 | 0.02045 | 0.005633 | 0.02305 | 0.003091 |
| stacked | 0 | 0 | 0.001161 | 0 | 0.009939 | 0 | 0.08278 | 0 | 0.3035 | 0 | 0.3602 | 0 |
| anisotropic | 0.000429 | 9.51e-05 | 0.001681 | 0.001058 | 0.003119 | 0.001329 | 0.004602 | 0.002681 | 0.02197 | 0.003827 | 0.03164 | 0.005111 |
| row alignment | 0.3974 | -0.3597 | 0.08567 | -0.1264 | -0.3858 | 0.04182 | 0.01516 | -0.01844 | 0.4266 | -0.07181 | 0.404 | 0.06811 |
| thin | 0.000123 | 6.34e-05 | 0.001558 | 0.001539 | 0.003089 | 0.002594 | 0.004691 | 0.003298 | 0.01371 | 0.005143 | 0.01387 | 0.007559 |
| one-node lines | 0.007371 | — | 0.05065 | — | 0.1036 | — | 0.106 | — | 0.2183 | — | 0.2729 | — |
| quad one-cell-thin | — | 0.000571 | — | 0.004105 | — | 0.006328 | — | 0.008507 | — | 0.01151 | — | 0.01503 |
| air temperature mean (°C) | -2.289 | -2.44 | -6.049 | -6.821 | -7.363 | -7.677 | -9.886 | -9.321 | -10.09 | 1.494 | -2.556 | 8.038 |
| precipitation mean (mm) | 1254 | 1254 | 1252 | 1255 | 1263 | 1266 | 1263 | 1256 | 1261 | 1256 | 1251 | 1267 |
| river fraction of land | — | — | 0.1357 | 0.1387 | 0.1629 | 0.1419 | 0.1667 | 0.2443 | 0.1962 | 0.2057 | 0.2063 | 0.1886 |
| HEALPix same node | 0.5552 | 0.5447 | 0.5529 | 0.5501 | 0.5526 | 0.5528 | 0.5464 | 0.553 | 0.5257 | 0.554 | 0.5186 | 0.5505 |

## Performance (mean seconds per step)

| s/step | seed 1 lines | seed 1 quad | seed 2 lines | seed 2 quad | seed 3 lines | seed 3 quad | seed 4 lines | seed 4 quad | seed 5 lines | seed 5 quad |
|---|---|---|---|---|---|---|---|---|---|---|
| step_total | 4.308 | 2.976 | 5.527 | 3.344 | 5.579 | 2.849 | 4.463 | 2.82 | 3.801 | 2.42 |
| deform_topology | 2.074 | 0.7553 | 2.934 | 0.9502 | 2.949 | 0.8185 | 2.183 | 0.7958 | 1.946 | 0.6411 |
| climate_erosion_hydrology | 0.8419 | 0.8005 | 0.8304 | 0.7603 | 0.8392 | 0.6831 | 0.7475 | 0.6897 | 0.6405 | 0.6126 |
| deform | 1.959 | 0.5263 | 2.785 | 0.6343 | 2.793 | 0.5497 | 2.058 | 0.5299 | 1.834 | 0.4337 |
| faults | 0.6835 | 0.8564 | 0.9167 | 1.1 | 0.9521 | 0.9005 | 0.7705 | 0.8862 | 0.6525 | 0.7783 |
| fluid_dynamics | 5.12e-06 | 3.03e-06 | 6.2e-06 | 1.01e-05 | 5.3e-06 | 2.5e-06 | 4.08e-06 | 2.83e-06 | 2.68e-06 | 2.05e-06 |
| gap_fill | 0.0159 | 0.1413 | 0.01713 | 0.2079 | 0.0172 | 0.1791 | 0.01454 | 0.1755 | 0.0131 | 0.1331 |
| magma_transport | 0.1534 | 0.1999 | 0.1728 | 0.1352 | 0.1304 | 0.08902 | 0.2002 | 0.1041 | 0.07344 | 0.07966 |
| overlap_tracking | 0.02255 | 0.04136 | 0.02882 | 0.05675 | 0.0284 | 0.04741 | 0.02431 | 0.04711 | 0.02064 | 0.03644 |
| record_stats | 0.008225 | 0.004862 | 0.009773 | 0.004291 | 0.01045 | 0.004032 | 0.00833 | 0.004101 | 0.006991 | 0.003494 |
| resource_formation | 0.02847 | 0.01316 | 0.0348 | 0.0135 | 0.03672 | 0.01192 | 0.02947 | 0.01214 | 0.02385 | 0.01057 |
| sea_level | 0.0687 | 0.05627 | 0.08511 | 0.06528 | 0.09755 | 0.06432 | 0.07664 | 0.05851 | 0.07679 | 0.06161 |
| shift | 0.407 | 0.28 | 0.4792 | 0.3028 | 0.4949 | 0.2658 | 0.3935 | 0.2581 | 0.3344 | 0.2236 |
| stranded_basins | 1.92e-05 | 1.82e-05 | 1.77e-05 | 1.31e-05 | 1.71e-05 | 1.09e-05 | 1.38e-05 | 1.17e-05 | 1.03e-05 | 9.64e-06 |
| topology | 0.0769 | 0.04638 | 0.1036 | 0.05121 | 0.1102 | 0.04236 | 0.08704 | 0.04337 | 0.07854 | 0.0378 |
| volcanism | 0.04125 | 0.008671 | 0.06125 | 0.01182 | 0.06516 | 0.01098 | 0.05131 | 0.009914 | 0.04462 | 0.008813 |

Index builds per step (mean):

- seed 1 lines: line_row_lookup 64.5, plate_node_kdtree 64.6, plate_outline 70.3, plate_outline_kdtree 69.3, world_node_kdtree 1.0
- seed 1 quad: plate_node_kdtree 49.7, plate_outline 63.8, plate_outline_kdtree 50.7, quad_adjacency 28.2, quad_boundary_loops 40.6, world_node_kdtree 1.0
- seed 2 lines: line_row_lookup 104.0, plate_node_kdtree 104.2, plate_outline 113.0, plate_outline_kdtree 110.4, world_node_kdtree 1.0
- seed 2 quad: plate_node_kdtree 76.3, plate_outline 99.5, plate_outline_kdtree 77.9, quad_adjacency 42.4, quad_boundary_loops 64.0, world_node_kdtree 1.0
- seed 3 lines: line_row_lookup 114.4, plate_node_kdtree 114.7, plate_outline 125.0, plate_outline_kdtree 121.7, world_node_kdtree 1.0
- seed 3 quad: plate_node_kdtree 71.5, plate_outline 93.9, plate_outline_kdtree 73.2, quad_adjacency 39.5, quad_boundary_loops 60.6, world_node_kdtree 1.0
- seed 4 lines: line_row_lookup 94.1, plate_node_kdtree 94.3, plate_outline 104.9, plate_outline_kdtree 100.3, world_node_kdtree 1.0
- seed 4 quad: plate_node_kdtree 71.2, plate_outline 93.2, plate_outline_kdtree 73.2, quad_adjacency 39.3, quad_boundary_loops 59.8, world_node_kdtree 1.0
- seed 5 lines: line_row_lookup 89.9, plate_node_kdtree 90.1, plate_outline 98.2, plate_outline_kdtree 95.5, world_node_kdtree 1.0
- seed 5 quad: plate_node_kdtree 57.3, plate_outline 74.0, plate_outline_kdtree 58.7, quad_adjacency 30.4, quad_boundary_loops 46.9, world_node_kdtree 1.0

