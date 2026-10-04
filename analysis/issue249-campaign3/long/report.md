# Line vs quad parity report

**Verdict: FAIL** — configs long; seeds 1, 2, 3, 4, 5.

Reproduce:

    bin/debug/surface_parity.py paired --preset long --seeds 1,2,3,4,5 --jobs 5 --out ../analysis/issue249-campaign3/long

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
| H11:field_caps | pass | violations=0, surface=lines, seed=1 |

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
| K1:hc_volume_km3_drift | fail | age_myr=240, quad=-0.06516, lines=-0.1837, gap=0.1185, warn_above=0.05, fail_above=0.1, seed=1 |
| K2:continental_hc_volume_km3_drift | fail | age_myr=240, quad=-0.1396, lines=-0.01223, gap=0.1274, warn_above=0.05, fail_above=0.1, seed=1 |
| K3:quad_area_accounting | pass | age_myr=30, quad=0.9982, expected=0.9982, warn_above=0.01, seed=1 |
| K4:hm_volume_km3_drift | fail | age_myr=240, quad=0.1254, lines=-0.01745, gap=0.1429, warn_above=0.05, fail_above=0.1, seed=1 |

## Mesh quality (quad)

| gate | status | worst case |
|---|---|---|
| M1:anisotropic_fraction | pass | age_myr=0, quad=0.000159, lines=0.000521, warn_above=0.01, seed=1 |
| M2:row_alignment | pass | age_myr=0, quad=0.08457, note=too few anisotropic nodes to have a direction, seed=1 |
| M3:thin_fraction | warn | age_myr=0, quad=0.000191, lines=0.000184, warn_above=0.000184, seed=1 |
| M4:aspect_gt_4 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M5:skew_gt_45 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M6:fragments | warn | age_myr=400, quad=6, lines=0, seed=1 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=1 |
| S2:stability_air_temperature_mean_c | pass | quad=0.2162, lines=0.1217, warn_above=0.2934, seed=1 |
| S2:stability_land_fraction | pass | quad=0.001105, lines=0.001618, warn_above=0.005236, seed=1 |
| S2:stability_ocean_temperature_mean_c | pass | quad=0.008222, lines=0.009419, warn_above=0.06884, seed=1 |
| S2:stability_precipitation_mean_mm | pass | quad=5.302, lines=5.341, warn_above=15.68, seed=1 |
| S2:stability_sea_level_m | pass | quad=8.549, lines=14.26, warn_above=29.51, seed=1 |

## Derived indexes: HEALPix vs KD-tree

| gate | status | worst case |
|---|---|---|
| X1:healpix_all | pass | age_myr=0, samples=200000, quad_same_node=0.5515, lines_same_node=0.5572, quad_distance_p95=1.195, lines_distance_p95=1.18, seed=1 |
| X1:healpix_antimeridian | pass | age_myr=0, samples=3333, quad_same_node=0.5581, lines_same_node=0.5467, quad_distance_p95=1.283, lines_distance_p95=1.216, seed=1 |
| X1:healpix_coast | pass | age_myr=0, samples=12454, quad_same_node=0.5499, lines_same_node=0.5449, quad_distance_p95=1.212, lines_distance_p95=1.22, seed=1 |
| X1:healpix_hole | warn | age_myr=120, samples=134, quad_same_node=0.5149, lines_same_node=0.561, quad_distance_p95=1.506, lines_distance_p95=1.177, seed=5 |
| X1:healpix_plate_boundary | pass | age_myr=0, samples=17562, quad_same_node=0.5399, lines_same_node=0.5446, quad_distance_p95=1.257, lines_distance_p95=1.233, seed=1 |
| X1:healpix_pole | pass | age_myr=0, samples=3038, quad_same_node=0.525, lines_same_node=0.5319, quad_distance_p95=1.353, lines_distance_p95=1.359, seed=1 |

## Performance

| gate | status | worst case |
|---|---|---|
| R1:deform_topology_s_per_step | pass | quad=0.7838, lines=1.431, ratio=0.5476, warn_above=1.25, seed=1 |
| R2:step_total_s_per_step | pass | quad=2.333, lines=3.052, ratio=0.7643, warn_above=1.25, seed=1 |

## Ensemble parity

| gate | status | worst case |
|---|---|---|
| P1:continental_hc_volume_km3 | pass | age_myr=30, seeds=5, lines_mean=4.15e+09, quad_mean=4.22e+09, lines_sigma=2.14e+09, delta=6.28e+07, fail_above=4.29e+09 |
| P1:elevation_p05 | fail | age_myr=60, seeds=5, lines_mean=-6094, quad_mean=-5739, lines_sigma=127.5, delta=355.1, fail_above=254.9 |
| P1:elevation_p50 | pass | age_myr=30, seeds=5, lines_mean=-5100, quad_mean=-5089, lines_sigma=63.5, delta=10.72, fail_above=127 |
| P1:elevation_p95 | pass | age_myr=30, seeds=5, lines_mean=492, quad_mean=546.8, lines_sigma=2516, delta=54.85, fail_above=5032 |
| P1:hc_volume_km3 | pass | age_myr=30, seeds=5, lines_mean=8.09e+09, quad_mean=8.12e+09, lines_sigma=1.86e+09, delta=2.71e+07, fail_above=3.72e+09 |
| P1:land_fraction | pass | age_myr=30, seeds=5, lines_mean=0.1454, quad_mean=0.1425, lines_sigma=0.1102, delta=-0.002869, fail_above=0.2204 |
| P1:plates | pass | age_myr=30, seeds=5, lines_mean=18.6, quad_mean=18.4, lines_sigma=2.608, delta=-0.2, fail_above=5.215 |
| P1:sea_level_m | fail | age_myr=240, seeds=5, lines_mean=-1496, quad_mean=-867.1, lines_sigma=286.6, delta=629, fail_above=573.2 |

## Checkpoints, seed 1

| seed 1 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 14 | 14 | 15 | 15 | 17 | 16 | 22 | 21 | 31 | 30 | 40 | 44 |
| nodes | 32639 | 31493 | 32657 | 31566 | 32506 | 31973 | 33096 | 32471 | 35801 | 32506 | 38640 | 33165 |
| land fraction | 0.2795 | 0.2793 | 0.3143 | 0.3125 | 0.314 | 0.3092 | 0.3104 | 0.2905 | 0.3488 | 0.2362 | 0.3557 | 0.1394 |
| sea level (m) | 0 | 0 | -712.5 | -609.7 | -810.9 | -679.4 | -1067 | -829.3 | -1466 | -1203 | -2037 | -1806 |
| Hc volume (km³) | 1.15e+10 | 1.15e+10 | 1.1e+10 | 1.11e+10 | 1.07e+10 | 1.1e+10 | 1.03e+10 | 1.1e+10 | 9.42e+09 | 1.08e+10 | 7.63e+09 | 8.78e+09 |
| continental Hc (km³) | 7.65e+09 | 7.65e+09 | 7.54e+09 | 7.7e+09 | 7.5e+09 | 7.72e+09 | 7.5e+09 | 7.78e+09 | 7.56e+09 | 6.58e+09 | 5.64e+09 | 3.5e+09 |
| elevation p05/p50/p95 (m) | -4768/-4354/2990 | -4768/-4353/2984 | -5537/-5016/3006 | -5473/-4981/3280 | -5948/-5061/3098 | -5677/-4994/3905 | -6240/-5054/3592 | -5820/-4925/6031 | -6229/-4156/5101 | -5906/-4939/6334 | -6204/-3912/4452 | -5905/-5279/5861 |
| uncovered | 0.02311 | 0.005985 | 0.02641 | 0.007625 | 0.03548 | 0.007175 | 0.04928 | 0.009645 | 0.04825 | 0.01218 | 0.06949 | 0.01623 |
| multiply covered | 0.01406 | 0.006475 | 0.002595 | 0.005805 | 0.008495 | 0.006185 | 0.004335 | 0.008235 | 0.2145 | 0.01171 | 0.1075 | 0.01042 |
| void | 0 | 0 | 0 | 0 | 0.000885 | 0 | 0.00418 | 0 | 0.00387 | 0 | 0.004115 | 0 |
| nodes inside other plate | 0.01618 | 0.009653 | 0.007135 | 0.002756 | 0.0163 | 0.003315 | 0.01218 | 0.004958 | 0.1357 | 0.009168 | 0.07081 | 0.002744 |
| stacked | 0 | 0 | 0.000429 | 0 | 0.004799 | 0 | 0.05892 | 0 | 0.2178 | 0 | 0.3226 | 0 |
| anisotropic | 0.000521 | 0.000159 | 0.002174 | 0.000982 | 0.004799 | 0.001376 | 0.005862 | 0.002371 | 0.01268 | 0.003292 | 0.02919 | 0.005397 |
| row alignment | 0.5571 | 0.08457 | 0.005663 | -0.1962 | -0.05471 | 0.1164 | -0.326 | -0.2864 | 0.4806 | -0.1725 | 0.4082 | 0.07087 |
| thin | 0.000184 | 0.000191 | 0.002205 | 0.001806 | 0.004276 | 0.003221 | 0.006375 | 0.002679 | 0.01193 | 0.004522 | 0.0183 | 0.005488 |
| one-node lines | 0.009498 | — | 0.06511 | — | 0.1246 | — | 0.1826 | — | 0.1545 | — | 0.3117 | — |
| quad one-cell-thin | — | 0.000381 | — | 0.004435 | — | 0.005473 | — | 0.005605 | — | 0.008768 | — | 0.01013 |
| air temperature mean (°C) | 4.49 | 4.481 | 7.632 | 6.992 | 7.903 | 6.919 | 7.503 | 3.03 | 4.105 | -3.885 | 7.546 | -11.91 |
| precipitation mean (mm) | 986.1 | 987 | 1013 | 1009 | 1010 | 1005 | 1038 | 1038 | 1065 | 1082 | 1139 | 1114 |
| river fraction of land | — | — | 0.1038 | 0.119 | 0.1032 | 0.117 | 0.1296 | 0.152 | 0.1523 | 0.2146 | 0.158 | 0.2712 |
| HEALPix same node | 0.5572 | 0.5515 | 0.5544 | 0.5504 | 0.5543 | 0.5507 | 0.5477 | 0.5479 | 0.5334 | 0.548 | 0.5214 | 0.5511 |

## Checkpoints, seed 2

| seed 2 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 18 | 18 | 21 | 20 | 24 | 19 | 29 | 26 | 45 | 40 | 66 | 52 |
| nodes | 32634 | 31334 | 33378 | 31465 | 33416 | 31764 | 35260 | 32272 | 36870 | 32890 | 40285 | 32873 |
| land fraction | 0.1406 | 0.1403 | 0.1513 | 0.1465 | 0.1336 | 0.1445 | 0.1412 | 0.1417 | 0.1254 | 0.11 | 0.09341 | 0.0498 |
| sea level (m) | 0 | 0 | -779.7 | -607.3 | -1033 | -701.9 | -1385 | -742.3 | -1690 | -954.6 | -2238 | -1374 |
| Hc volume (km³) | 8.93e+09 | 8.91e+09 | 8.4e+09 | 8.44e+09 | 7.53e+09 | 8.33e+09 | 6.7e+09 | 8.02e+09 | 5.74e+09 | 7.65e+09 | 3.94e+09 | 5.91e+09 |
| continental Hc (km³) | 4.29e+09 | 4.28e+09 | 4.3e+09 | 4.4e+09 | 3.67e+09 | 4.41e+09 | 3.56e+09 | 3.99e+09 | 3.06e+09 | 2.99e+09 | 1.36e+09 | 1.2e+09 |
| elevation p05/p50/p95 (m) | -5179/-4437/1657 | -5174/-4437/1652 | -5804/-5105/1770 | -5532/-5093/1863 | -6256/-5183/1597 | -5775/-5136/1908 | -6420/-5242/1620 | -5874/-5124/1769 | -6413/-5344/1974 | -5907/-5266/4110 | -6428/-5407/-442.4 | -5941/-5385/-1383 |
| uncovered | 0.02514 | 0.00755 | 0.03578 | 0.00922 | 0.04899 | 0.009215 | 0.05709 | 0.01098 | 0.07726 | 0.01505 | 0.09827 | 0.01843 |
| multiply covered | 0.09836 | 0.007265 | 0.01236 | 0.0066 | 0.01325 | 0.00676 | 0.00913 | 0.009745 | 0.01154 | 0.0106 | 0.004505 | 0.0135 |
| void | 0 | 0 | 0 | 0 | 0.000705 | 0 | 0.000245 | 0 | 0.001095 | 0 | 0.00141 | 0 |
| nodes inside other plate | 0.08234 | 0.01098 | 0.02187 | 0.00375 | 0.02008 | 0.003432 | 0.02155 | 0.006166 | 0.02829 | 0.004561 | 0.02681 | 0.00651 |
| stacked | 0 | 0 | 0.03994 | 0 | 0.07326 | 0 | 0.1616 | 0 | 0.2343 | 0 | 0.3775 | 0 |
| anisotropic | 0.000276 | 0.000287 | 0.003895 | 0.001367 | 0.004788 | 0.002267 | 0.007771 | 0.00282 | 0.02411 | 0.003983 | 0.03795 | 0.005993 |
| row alignment | 0.4286 | -0.4423 | -0.289 | -0.1641 | 0.1419 | 0.1745 | 0.1469 | 0.08646 | 0.3313 | 0.01397 | 0.3968 | -0.03743 |
| thin | 0.000123 | 9.57e-05 | 0.003715 | 0.00197 | 0.005207 | 0.002645 | 0.006722 | 0.003842 | 0.01679 | 0.004926 | 0.02445 | 0.007088 |
| one-node lines | 0.01008 | — | 0.09904 | — | 0.08882 | — | 0.1409 | — | 0.2394 | — | 0.3273 | — |
| quad one-cell-thin | — | 0.000734 | — | 0.004513 | — | 0.005258 | — | 0.007375 | — | 0.01012 | — | 0.01354 |
| air temperature mean (°C) | -3.378 | -3.453 | 0.5634 | -0.7211 | 2.76 | -1.371 | 4.149 | -0.7591 | 5.716 | -7.245 | 7.747 | -3.09 |
| precipitation mean (mm) | 1132 | 1132 | 1168 | 1169 | 1181 | 1179 | 1184 | 1191 | 1207 | 1178 | 1252 | 1239 |
| river fraction of land | — | — | 0.1435 | 0.1611 | 0.1711 | 0.1773 | 0.1733 | 0.18 | 0.1712 | 0.2126 | 0.178 | 0.2885 |
| HEALPix same node | 0.5526 | 0.5486 | 0.5474 | 0.5524 | 0.5442 | 0.551 | 0.5358 | 0.5469 | 0.5275 | 0.5531 | 0.5117 | 0.5466 |

## Checkpoints, seed 3

| seed 3 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 18 | 18 | 21 | 21 | 26 | 23 | 32 | 28 | 50 | 32 | 53 | 41 |
| nodes | 32632 | 31213 | 32625 | 31453 | 32906 | 32093 | 35745 | 32863 | 38164 | 32850 | 41852 | 32535 |
| land fraction | 0.05561 | 0.05551 | 0.07864 | 0.07493 | 0.08236 | 0.07897 | 0.06659 | 0.06943 | 0.05865 | 0.04962 | 0.04181 | 0.02356 |
| sea level (m) | 0 | 0 | -536 | -478.3 | -755.6 | -562.9 | -1380 | -708.2 | -1853 | -900.3 | -2258 | -1174 |
| Hc volume (km³) | 7.2e+09 | 7.2e+09 | 6.86e+09 | 6.91e+09 | 6.5e+09 | 7.05e+09 | 5.39e+09 | 6.33e+09 | 4.02e+09 | 5.92e+09 | 3.31e+09 | 5.2e+09 |
| continental Hc (km³) | 2.59e+09 | 2.59e+09 | 2.74e+09 | 2.79e+09 | 2.62e+09 | 2.95e+09 | 2.01e+09 | 2.1e+09 | 1.04e+09 | 1.63e+09 | 7.09e+08 | 5.38e+08 |
| elevation p05/p50/p95 (m) | -5222/-4494/40.48 | -5222/-4493/41.41 | -5855/-5132/151.2 | -5554/-5125/53.87 | -6192/-5204/368.1 | -5761/-5177/401.4 | -6339/-5353/-549.4 | -5873/-5257/109.2 | -6448/-5512/-1404 | -5939/-5392/-910.8 | -6491/-5590/-2867 | -5987/-5470/-4279 |
| uncovered | 0.02706 | 0.00765 | 0.04494 | 0.009745 | 0.06031 | 0.01052 | 0.05826 | 0.01149 | 0.08195 | 0.01404 | 0.08327 | 0.01563 |
| multiply covered | 0.01106 | 0.00755 | 0.005425 | 0.007525 | 0.005095 | 0.007455 | 0.00666 | 0.00778 | 0.00392 | 0.009965 | 0.005385 | 0.0106 |
| void | 0 | 0 | 0.0015 | 0 | 0.00508 | 0 | 0.00068 | 0 | 0.000825 | 0 | 0.000845 | 0 |
| nodes inside other plate | 0.01492 | 0.01157 | 0.01655 | 0.004229 | 0.01501 | 0.003147 | 0.02171 | 0.002952 | 0.01462 | 0.002648 | 0.0233 | 0.001076 |
| stacked | 0 | 0 | 0.01116 | 0 | 0.05333 | 0 | 0.1811 | 0 | 0.3157 | 0 | 0.4086 | 0 |
| anisotropic | 0.000552 | 0.000192 | 0.00328 | 0.001812 | 0.005379 | 0.002399 | 0.01111 | 0.003104 | 0.01845 | 0.003927 | 0.07438 | 0.005502 |
| row alignment | 0.7907 | 0.3694 | -0.2428 | 0.008403 | -0.05305 | -0.0341 | -0.1601 | 0.1542 | 0.3072 | -0.08415 | 0.6981 | 0.00345 |
| thin | 6.13e-05 | 0.000192 | 0.003954 | 0.002067 | 0.009755 | 0.002867 | 0.01234 | 0.003408 | 0.01242 | 0.004871 | 0.03238 | 0.006393 |
| one-node lines | 0.006969 | — | 0.1034 | — | 0.1099 | — | 0.2123 | — | 0.2042 | — | 0.3192 | — |
| quad one-cell-thin | — | 0.000609 | — | 0.004483 | — | 0.005578 | — | 0.006664 | — | 0.01023 | — | 0.01239 |
| air temperature mean (°C) | -9.087 | -9.163 | -2.935 | -3.652 | -8.748 | -13.08 | -3.941 | -9.61 | 4.527 | -3.539 | -5.604 | -14.6 |
| precipitation mean (mm) | 1201 | 1201 | 1239 | 1234 | 1246 | 1233 | 1256 | 1237 | 1253 | 1236 | 1259 | 1242 |
| river fraction of land | — | — | 0.1475 | 0.1573 | 0.2383 | 0.2043 | 0.1902 | 0.1635 | 0.1863 | 0.1716 | 0.2557 | 0.2818 |
| HEALPix same node | 0.5533 | 0.556 | 0.5556 | 0.5534 | 0.5491 | 0.5546 | 0.5359 | 0.5496 | 0.5233 | 0.5532 | 0.5106 | 0.5528 |

## Checkpoints, seed 4

| seed 4 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 17 | 17 | 19 | 19 | 24 | 21 | 29 | 30 | 39 | 34 | 48 | 43 |
| nodes | 32636 | 31134 | 32567 | 31243 | 33088 | 31613 | 34796 | 32576 | 35784 | 32473 | 38687 | 32567 |
| land fraction | 0.1417 | 0.1416 | 0.1611 | 0.159 | 0.1672 | 0.1601 | 0.1701 | 0.1531 | 0.1486 | 0.1289 | 0.1281 | 0.06749 |
| sea level (m) | 0 | 0 | -561 | -526.6 | -714.9 | -533.8 | -1055 | -555.6 | -1351 | -941 | -1988 | -1273 |
| Hc volume (km³) | 8.49e+09 | 8.49e+09 | 8.17e+09 | 8.12e+09 | 7.96e+09 | 7.98e+09 | 7.52e+09 | 8.05e+09 | 6.14e+09 | 8.41e+09 | 5.17e+09 | 6.2e+09 |
| continental Hc (km³) | 4.13e+09 | 4.13e+09 | 4.24e+09 | 4.24e+09 | 4.31e+09 | 4.34e+09 | 4.34e+09 | 4.51e+09 | 3.24e+09 | 4.01e+09 | 2.43e+09 | 1.33e+09 |
| elevation p05/p50/p95 (m) | -5222/-4427/958.8 | -5222/-4427/962.3 | -5657/-5065/1142 | -5508/-5062/1081 | -6065/-5101/1183 | -5734/-5075/1147 | -6261/-5139/1279 | -5827/-5053/1752 | -6323/-4991/1299 | -5939/-5283/5826 | -6449/-5386/727.7 | -5967/-5406/-274.1 |
| uncovered | 0.02816 | 0.00778 | 0.03876 | 0.00971 | 0.04596 | 0.01026 | 0.06164 | 0.0128 | 0.07119 | 0.01337 | 0.07251 | 0.016 |
| multiply covered | 0.01667 | 0.00729 | 0.002045 | 0.006415 | 0.005055 | 0.00669 | 0.01461 | 0.01009 | 0.002735 | 0.01218 | 0.01218 | 0.01152 |
| void | 0 | 0 | 0.00017 | 0 | 0.00027 | 0 | 0.00107 | 0 | 0.00099 | 0 | 0.00035 | 0 |
| nodes inside other plate | 0.01906 | 0.01092 | 0.007676 | 0.001184 | 0.01584 | 0.002214 | 0.02974 | 0.004942 | 0.01417 | 0.007391 | 0.03332 | 0.005404 |
| stacked | 0 | 0 | 0.006233 | 0 | 0.03962 | 0 | 0.1348 | 0 | 0.201 | 0 | 0.2911 | 0 |
| anisotropic | 0.00046 | 0.000225 | 0.003163 | 0.001824 | 0.008281 | 0.002341 | 0.01058 | 0.002793 | 0.0218 | 0.003357 | 0.03236 | 0.004974 |
| row alignment | 0.3941 | 0.03748 | -0.3379 | 0.07459 | -0.1497 | 0.07542 | 0.02997 | 0.02543 | 0.2699 | 0.09947 | 0.5594 | -0.0857 |
| thin | 0.000184 | 3.21e-05 | 0.004882 | 0.002433 | 0.008372 | 0.003669 | 0.009829 | 0.003868 | 0.01336 | 0.004527 | 0.01515 | 0.004913 |
| one-node lines | 0.009368 | — | 0.09032 | — | 0.1679 | — | 0.1766 | — | 0.2224 | — | 0.2353 | — |
| quad one-cell-thin | — | 0.000418 | — | 0.005473 | — | 0.006295 | — | 0.008964 | — | 0.008992 | — | 0.01093 |
| air temperature mean (°C) | -5.967 | -5.829 | -7.379 | -7.288 | -7.705 | -8.102 | -2.663 | -7.381 | 3.256 | -12.7 | 9.181 | 0.405 |
| precipitation mean (mm) | 1146 | 1147 | 1238 | 1251 | 1170 | 1176 | 1184 | 1176 | 1226 | 1150 | 1196 | 1203 |
| river fraction of land | — | — | 0.1637 | 0.1583 | 0.1628 | 0.1571 | 0.1538 | 0.1816 | 0.1629 | 0.2347 | 0.1734 | 0.2612 |
| HEALPix same node | 0.5553 | 0.5487 | 0.5532 | 0.5505 | 0.5507 | 0.5494 | 0.5374 | 0.5476 | 0.5359 | 0.5475 | 0.5211 | 0.5523 |

## Checkpoints, seed 5

| seed 5 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 16 | 16 | 17 | 17 | 19 | 17 | 24 | 23 | 44 | 30 | 59 | 46 |
| nodes | 32643 | 31550 | 32711 | 31505 | 32615 | 31545 | 33372 | 32093 | 36755 | 32863 | 39871 | 32920 |
| land fraction | 0.01461 | 0.01469 | 0.02159 | 0.01965 | 0.02522 | 0.02612 | 0.03185 | 0.03645 | 0.0503 | 0.06263 | 0.06085 | 0.05628 |
| sea level (m) | 0 | 0 | -457.7 | -429.4 | -495.4 | -426.2 | -665.6 | -385.8 | -1120 | -336.8 | -1509 | -387.3 |
| Hc volume (km³) | 6.19e+09 | 6.19e+09 | 6.08e+09 | 6.04e+09 | 5.94e+09 | 6.04e+09 | 5.58e+09 | 6.07e+09 | 4.94e+09 | 6.65e+09 | 4.53e+09 | 6.35e+09 |
| continental Hc (km³) | 1.86e+09 | 1.86e+09 | 1.94e+09 | 1.95e+09 | 2.02e+09 | 2.03e+09 | 2.02e+09 | 2.13e+09 | 1.92e+09 | 2.14e+09 | 1.64e+09 | 1.68e+09 |
| elevation p05/p50/p95 (m) | -5222/-4570/-4197 | -5222/-4570/-4198 | -5532/-5183/-3610 | -5511/-5187/-3544 | -6011/-5187/-3178 | -5750/-5216/-2868 | -6309/-5240/-2554 | -5883/-5222/-1649 | -6382/-5345/-1107 | -5948/-5326/259.3 | -6401/-5389/-1256 | -5948/-5373/262.5 |
| uncovered | 0.02357 | 0.006985 | 0.03218 | 0.007625 | 0.04061 | 0.00874 | 0.05264 | 0.01043 | 0.07284 | 0.01306 | 0.08543 | 0.01638 |
| multiply covered | 0.03266 | 0.00718 | 0.00512 | 0.006345 | 0.002615 | 0.007125 | 0.005555 | 0.00887 | 0.00297 | 0.00999 | 0.003715 | 0.01214 |
| void | 0 | 0 | 0 | 0 | 5e-06 | 0 | 0.00017 | 0 | 0.000155 | 0 | 0.00088 | 0 |
| nodes inside other plate | 0.0367 | 0.0109 | 0.01241 | 0.002857 | 0.008616 | 0.003836 | 0.01244 | 0.004892 | 0.01681 | 0.003682 | 0.01547 | 0.003858 |
| stacked | 0 | 0 | 0.001345 | 0 | 0.01263 | 0 | 0.06844 | 0 | 0.2377 | 0 | 0.3602 | 0 |
| anisotropic | 0.000429 | 9.51e-05 | 0.001529 | 0.001524 | 0.004293 | 0.001268 | 0.004195 | 0.002493 | 0.01861 | 0.00353 | 0.03812 | 0.005711 |
| row alignment | 0.3974 | -0.3597 | -0.05293 | 0.1199 | -0.1406 | -0.02741 | 0.2143 | 0.06737 | 0.1823 | 0.03979 | 0.5579 | -0.0267 |
| thin | 0.000123 | 6.34e-05 | 0.001498 | 0.002 | 0.004906 | 0.002568 | 0.005064 | 0.003864 | 0.01521 | 0.004138 | 0.01706 | 0.007625 |
| one-node lines | 0.007371 | — | 0.04772 | — | 0.1111 | — | 0.09831 | — | 0.248 | — | 0.2341 | — |
| quad one-cell-thin | — | 0.000571 | — | 0.003555 | — | 0.005516 | — | 0.006762 | — | 0.00849 | — | 0.01367 |
| air temperature mean (°C) | -2.289 | -2.44 | -6.466 | -6.276 | -7.534 | -6.809 | -10.37 | -9.617 | -7.978 | -24.33 | -2.054 | -23.09 |
| precipitation mean (mm) | 1254 | 1254 | 1252 | 1255 | 1263 | 1263 | 1265 | 1246 | 1266 | 1219 | 1255 | 1205 |
| river fraction of land | — | — | 0.1389 | 0.1436 | 0.1583 | 0.1359 | 0.1711 | 0.2535 | 0.2274 | 0.3226 | 0.2056 | 0.3004 |
| HEALPix same node | 0.5552 | 0.5447 | 0.5523 | 0.5513 | 0.552 | 0.5514 | 0.5476 | 0.5524 | 0.5316 | 0.5443 | 0.5161 | 0.5454 |

## Performance (mean seconds per step)

| s/step | seed 1 lines | seed 1 quad | seed 2 lines | seed 2 quad | seed 3 lines | seed 3 quad | seed 4 lines | seed 4 quad | seed 5 lines | seed 5 quad |
|---|---|---|---|---|---|---|---|---|---|---|
| step_total | 3.052 | 2.333 | 3.811 | 2.504 | 3.974 | 2.962 | 3.884 | 2.945 | 3.668 | 2.708 |
| deform_topology | 1.431 | 0.7838 | 1.986 | 0.8836 | 2.115 | 1.061 | 1.901 | 1.012 | 1.826 | 0.9256 |
| climate_erosion_hydrology | 0.6093 | 0.5572 | 0.5991 | 0.5489 | 0.6037 | 0.6342 | 0.6805 | 0.6634 | 0.6133 | 0.6245 |
| deform | 1.345 | 0.3768 | 1.872 | 0.3883 | 2.007 | 0.4912 | 1.788 | 0.4658 | 1.728 | 0.4561 |
| faults | 0.4921 | 0.606 | 0.6462 | 0.6939 | 0.6538 | 0.8104 | 0.6515 | 0.7988 | 0.6439 | 0.7256 |
| fluid_dynamics | 2.22e-06 | 1.64e-06 | 2.67e-06 | 1.87e-06 | 2.45e-06 | 2.73e-06 | 3.16e-06 | 3.99e-06 | 3.12e-06 | 2.48e-06 |
| gap_fill | 0.01169 | 0.3431 | 0.01273 | 0.4243 | 0.01248 | 0.4844 | 0.013 | 0.4624 | 0.01202 | 0.3967 |
| magma_transport | 0.1207 | 0.1186 | 0.09658 | 0.09524 | 0.08181 | 0.1008 | 0.1462 | 0.1346 | 0.08015 | 0.07972 |
| overlap_tracking | 0.01685 | 0.03297 | 0.02044 | 0.03699 | 0.0213 | 0.04286 | 0.02016 | 0.04164 | 0.01888 | 0.03649 |
| record_stats | 0.006086 | 0.003422 | 0.007211 | 0.003248 | 0.007159 | 0.003984 | 0.007809 | 0.004157 | 0.006415 | 0.00363 |
| resource_formation | 0.01978 | 0.009741 | 0.02404 | 0.009578 | 0.02456 | 0.01263 | 0.02549 | 0.01205 | 0.02158 | 0.01067 |
| sea_level | 0.05441 | 0.04601 | 0.06978 | 0.05271 | 0.07676 | 0.06217 | 0.06928 | 0.05869 | 0.07369 | 0.06139 |
| shift | 0.2852 | 0.1983 | 0.3312 | 0.2054 | 0.3556 | 0.2637 | 0.3518 | 0.2496 | 0.3558 | 0.2652 |
| stranded_basins | 1.38e-05 | 2.62e-05 | 9.98e-06 | 2e-05 | 8.79e-06 | 2.38e-05 | 1.41e-05 | 2.46e-05 | 1.45e-05 | 1.47e-05 |
| topology | 0.05766 | 0.03084 | 0.08032 | 0.03405 | 0.07397 | 0.04232 | 0.07973 | 0.04171 | 0.06758 | 0.03624 |
| volcanism | 0.03049 | 0.008841 | 0.04668 | 0.009622 | 0.05054 | 0.01211 | 0.04667 | 0.01094 | 0.04351 | 0.01006 |

Index builds per step (mean):

- seed 1 lines: line_row_lookup 73.5, plate_node_kdtree 73.6, plate_outline 79.2, plate_outline_kdtree 78.0, world_node_kdtree 1.0
- seed 1 quad: plate_node_kdtree 75.8, plate_outline 84.5, plate_outline_kdtree 77.9, quad_adjacency 40.8, quad_boundary_loops 57.4, world_node_kdtree 1.0
- seed 2 lines: line_row_lookup 105.4, plate_node_kdtree 105.7, plate_outline 114.5, plate_outline_kdtree 111.9, world_node_kdtree 1.0
- seed 2 quad: plate_node_kdtree 95.3, plate_outline 106.4, plate_outline_kdtree 98.0, quad_adjacency 49.6, quad_boundary_loops 72.6, world_node_kdtree 1.0
- seed 3 lines: line_row_lookup 113.0, plate_node_kdtree 113.3, plate_outline 122.8, plate_outline_kdtree 120.4, world_node_kdtree 1.0
- seed 3 quad: plate_node_kdtree 91.5, plate_outline 103.5, plate_outline_kdtree 94.5, quad_adjacency 47.5, quad_boundary_loops 70.9, world_node_kdtree 1.0
- seed 4 lines: line_row_lookup 91.4, plate_node_kdtree 91.6, plate_outline 100.4, plate_outline_kdtree 97.1, world_node_kdtree 1.0
- seed 4 quad: plate_node_kdtree 87.5, plate_outline 99.4, plate_outline_kdtree 91.3, quad_adjacency 45.9, quad_boundary_loops 67.7, world_node_kdtree 1.0
- seed 5 lines: line_row_lookup 100.6, plate_node_kdtree 100.8, plate_outline 110.5, plate_outline_kdtree 108.0, world_node_kdtree 1.0
- seed 5 quad: plate_node_kdtree 79.6, plate_outline 91.0, plate_outline_kdtree 83.3, quad_adjacency 41.4, quad_boundary_loops 62.1, world_node_kdtree 1.0

