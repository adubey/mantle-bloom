# Line vs quad parity report

**Verdict: FAIL** — configs long; seeds 1, 2, 3, 4, 5.

Reproduce:

    bin/debug/surface_parity.py paired --preset long --seeds 1,2,3,4,5 --jobs 5 --out ../analysis/issue249-final/long

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
| C3:multiply_covered | warn | age_myr=240, quad=0.01579, warn_above=0.015, fail_above=0.02, surface=quad, seed=2 |
| C4:nodes_inside_other_plate | pass | age_myr=0, quad=0.009653, lines=0.01618, fail_above=0.02118, seed=1 |

## Conservation (actual node/cell areas)

| gate | status | worst case |
|---|---|---|
| K1:hc_volume_km3_drift | fail | age_myr=400, quad=-0.1628, lines=-0.3381, gap=0.1753, warn_above=0.05, fail_above=0.1, seed=1 |
| K2:continental_hc_volume_km3_drift | fail | age_myr=240, quad=-0.1872, lines=-0.01223, gap=0.175, warn_above=0.05, fail_above=0.1, seed=1 |
| K3:quad_area_accounting | pass | age_myr=30, quad=1, expected=1, warn_above=0.01, seed=1 |
| K4:hm_volume_km3_drift | fail | age_myr=400, quad=-0.03762, lines=-0.3125, gap=0.2749, warn_above=0.05, fail_above=0.1, seed=1 |

## Mesh quality (quad)

| gate | status | worst case |
|---|---|---|
| M1:anisotropic_fraction | pass | age_myr=0, quad=0.000159, lines=0.000521, warn_above=0.01, seed=1 |
| M2:row_alignment | pass | age_myr=0, quad=0.08457, note=too few anisotropic nodes to have a direction, seed=1 |
| M3:thin_fraction | warn | age_myr=0, quad=0.000191, lines=0.000184, warn_above=0.000184, seed=1 |
| M4:aspect_gt_4 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M5:skew_gt_45 | pass | age_myr=0, quad=0, warn_above=0.001, seed=1 |
| M6:fragments | warn | age_myr=30, quad=1, lines=8, seed=1 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=lines, seed=1 |
| S2:stability_air_temperature_mean_c | pass | quad=0.1201, lines=0.1217, warn_above=0.2934, seed=1 |
| S2:stability_land_fraction | pass | quad=0.00081, lines=0.001618, warn_above=0.005236, seed=1 |
| S2:stability_ocean_temperature_mean_c | pass | quad=0.006929, lines=0.009419, warn_above=0.06884, seed=1 |
| S2:stability_precipitation_mean_mm | pass | quad=5.066, lines=5.341, warn_above=15.68, seed=1 |
| S2:stability_sea_level_m | pass | quad=5.624, lines=14.26, warn_above=29.51, seed=1 |

## Derived indexes: HEALPix vs KD-tree

| gate | status | worst case |
|---|---|---|
| X1:healpix_all | pass | age_myr=0, samples=200000, quad_same_node=0.5515, lines_same_node=0.5572, quad_distance_p95=1.195, lines_distance_p95=1.18, seed=1 |
| X1:healpix_antimeridian | pass | age_myr=0, samples=3333, quad_same_node=0.5581, lines_same_node=0.5467, quad_distance_p95=1.283, lines_distance_p95=1.216, seed=1 |
| X1:healpix_coast | pass | age_myr=0, samples=12454, quad_same_node=0.5499, lines_same_node=0.5449, quad_distance_p95=1.212, lines_distance_p95=1.22, seed=1 |
| X1:healpix_hole | warn | age_myr=120, samples=165, quad_same_node=0.5273, lines_same_node=0.5777, quad_distance_p95=1.144, lines_distance_p95=1.29, seed=2 |
| X1:healpix_plate_boundary | pass | age_myr=0, samples=17562, quad_same_node=0.5399, lines_same_node=0.5446, quad_distance_p95=1.257, lines_distance_p95=1.233, seed=1 |
| X1:healpix_pole | pass | age_myr=0, samples=3038, quad_same_node=0.525, lines_same_node=0.5319, quad_distance_p95=1.353, lines_distance_p95=1.359, seed=1 |

## Performance

| gate | status | worst case |
|---|---|---|
| R1:deform_topology_s_per_step | pass | quad=0.7428, lines=1.46, ratio=0.5087, warn_above=1.25, seed=1 |
| R2:step_total_s_per_step | pass | quad=2.334, lines=3.089, ratio=0.7554, warn_above=1.25, seed=1 |

## Ensemble parity

| gate | status | worst case |
|---|---|---|
| P1:continental_hc_volume_km3 | pass | age_myr=30, seeds=5, lines_mean=4.15e+09, quad_mean=4.22e+09, lines_sigma=2.14e+09, delta=6.85e+07, fail_above=4.29e+09 |
| P1:elevation_p05 | fail | age_myr=60, seeds=5, lines_mean=-6094, quad_mean=-5723, lines_sigma=127.5, delta=371.5, fail_above=254.9 |
| P1:elevation_p50 | pass | age_myr=30, seeds=5, lines_mean=-5100, quad_mean=-5087, lines_sigma=63.5, delta=12.57, fail_above=127 |
| P1:elevation_p95 | pass | age_myr=30, seeds=5, lines_mean=492, quad_mean=553.5, lines_sigma=2516, delta=61.49, fail_above=5032 |
| P1:hc_volume_km3 | pass | age_myr=30, seeds=5, lines_mean=8.09e+09, quad_mean=8.12e+09, lines_sigma=1.86e+09, delta=2.6e+07, fail_above=3.72e+09 |
| P1:land_fraction | pass | age_myr=30, seeds=5, lines_mean=0.1454, quad_mean=0.1422, lines_sigma=0.1102, delta=-0.003163, fail_above=0.2204 |
| P1:plates | pass | age_myr=30, seeds=5, lines_mean=18.6, quad_mean=18.4, lines_sigma=2.608, delta=-0.2, fail_above=5.215 |
| P1:sea_level_m | fail | age_myr=240, seeds=5, lines_mean=-1496, quad_mean=-792.9, lines_sigma=286.6, delta=703.2, fail_above=573.2 |

## Checkpoints, seed 1

| seed 1 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 14 | 14 | 15 | 15 | 17 | 17 | 22 | 20 | 31 | 27 | 40 | 36 |
| nodes | 32639 | 31493 | 32657 | 31632 | 32506 | 31976 | 33096 | 32749 | 35801 | 33126 | 38640 | 32721 |
| land fraction | 0.2795 | 0.2793 | 0.3143 | 0.3119 | 0.314 | 0.3039 | 0.3104 | 0.2817 | 0.3488 | 0.231 | 0.3557 | 0.1985 |
| sea level (m) | 0 | 0 | -712.5 | -639.9 | -810.9 | -700.6 | -1067 | -811.3 | -1466 | -1085 | -2037 | -1490 |
| Hc volume (km³) | 1.15e+10 | 1.15e+10 | 1.1e+10 | 1.11e+10 | 1.07e+10 | 1.1e+10 | 1.03e+10 | 1.12e+10 | 9.42e+09 | 1.05e+10 | 7.63e+09 | 9.65e+09 |
| continental Hc (km³) | 7.65e+09 | 7.65e+09 | 7.54e+09 | 7.71e+09 | 7.5e+09 | 7.77e+09 | 7.5e+09 | 7.79e+09 | 7.56e+09 | 6.22e+09 | 5.64e+09 | 4.73e+09 |
| elevation p05/p50/p95 (m) | -4768/-4354/2990 | -4768/-4353/2984 | -5537/-5016/3006 | -5473/-4986/3330 | -5948/-5061/3098 | -5662/-4981/4403 | -6240/-5054/3592 | -5783/-4777/6336 | -6229/-4156/5101 | -5884/-4766/6337 | -6204/-3912/4452 | -5905/-5140/5968 |
| uncovered | 0.02311 | 0.005985 | 0.02641 | 0.00697 | 0.03548 | 0.008125 | 0.04928 | 0.00904 | 0.04825 | 0.01076 | 0.06949 | 0.01249 |
| multiply covered | 0.01406 | 0.006475 | 0.002595 | 0.007395 | 0.008495 | 0.0081 | 0.004335 | 0.01106 | 0.2145 | 0.01353 | 0.1075 | 0.01418 |
| void | 0 | 0 | 0 | 0 | 0.000885 | 0 | 0.00418 | 0 | 0.00387 | 0 | 0.004115 | 0 |
| nodes inside other plate | 0.01618 | 0.009653 | 0.007135 | 0.004205 | 0.0163 | 0.004597 | 0.01218 | 0.00858 | 0.1357 | 0.01081 | 0.07081 | 0.008099 |
| stacked | 0 | 0 | 0.000429 | 0 | 0.004799 | 0 | 0.05892 | 0 | 0.2178 | 0 | 0.3226 | 0 |
| anisotropic | 0.000521 | 0.000159 | 0.002174 | 0.001897 | 0.004799 | 0.00197 | 0.005862 | 0.002443 | 0.01268 | 0.004981 | 0.02919 | 0.00382 |
| row alignment | 0.5571 | 0.08457 | 0.005663 | 0.002959 | -0.05471 | 0.2231 | -0.326 | -0.05709 | 0.4806 | -0.09235 | 0.4082 | 0.02773 |
| thin | 0.000184 | 0.000191 | 0.002205 | 0.003667 | 0.004276 | 0.003972 | 0.006375 | 0.003328 | 0.01193 | 0.00643 | 0.0183 | 0.005532 |
| one-node lines | 0.009498 | — | 0.06511 | — | 0.1246 | — | 0.1826 | — | 0.1545 | — | 0.3117 | — |
| quad one-cell-thin | — | 0.000381 | — | 0.007208 | — | 0.007412 | — | 0.008183 | — | 0.01072 | — | 0.01232 |
| air temperature mean (°C) | 4.49 | 4.481 | 7.632 | 6.52 | 7.903 | 5.89 | 7.503 | 1.32 | 4.105 | -0.8958 | 7.546 | -2.807 |
| precipitation mean (mm) | 986.1 | 987 | 1013 | 1008 | 1010 | 993.9 | 1038 | 1046 | 1065 | 1094 | 1139 | 1138 |
| river fraction of land | — | — | 0.1038 | 0.1143 | 0.1032 | 0.1318 | 0.1296 | 0.1608 | 0.1523 | 0.183 | 0.158 | 0.2179 |
| HEALPix same node | 0.5572 | 0.5515 | 0.5544 | 0.5506 | 0.5543 | 0.5492 | 0.5477 | 0.5504 | 0.5334 | 0.5512 | 0.5214 | 0.5457 |

## Checkpoints, seed 2

| seed 2 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 18 | 18 | 21 | 20 | 24 | 21 | 29 | 26 | 45 | 36 | 66 | 50 |
| nodes | 32634 | 31334 | 33378 | 31561 | 33416 | 31858 | 35260 | 32492 | 36870 | 32796 | 40285 | 32768 |
| land fraction | 0.1406 | 0.1403 | 0.1513 | 0.1464 | 0.1336 | 0.147 | 0.1412 | 0.1446 | 0.1254 | 0.1192 | 0.09341 | 0.04271 |
| sea level (m) | 0 | 0 | -779.7 | -616.4 | -1033 | -686.4 | -1385 | -672.5 | -1690 | -1006 | -2238 | -1427 |
| Hc volume (km³) | 8.93e+09 | 8.91e+09 | 8.4e+09 | 8.43e+09 | 7.53e+09 | 8.38e+09 | 6.7e+09 | 8.22e+09 | 5.74e+09 | 8.09e+09 | 3.94e+09 | 5.81e+09 |
| continental Hc (km³) | 4.29e+09 | 4.28e+09 | 4.3e+09 | 4.41e+09 | 3.67e+09 | 4.44e+09 | 3.56e+09 | 4.2e+09 | 3.06e+09 | 3.36e+09 | 1.36e+09 | 1.08e+09 |
| elevation p05/p50/p95 (m) | -5179/-4437/1657 | -5174/-4437/1652 | -5804/-5105/1770 | -5537/-5089/1876 | -6256/-5183/1597 | -5776/-5123/1864 | -6420/-5242/1620 | -5876/-5090/2338 | -6413/-5344/1974 | -5917/-5319/5251 | -6428/-5407/-442.4 | -5938/-5411/-2622 |
| uncovered | 0.02514 | 0.00755 | 0.03578 | 0.008205 | 0.04899 | 0.008495 | 0.05709 | 0.01095 | 0.07726 | 0.01328 | 0.09827 | 0.01542 |
| multiply covered | 0.09836 | 0.007265 | 0.01236 | 0.009035 | 0.01325 | 0.0101 | 0.00913 | 0.01152 | 0.01154 | 0.01579 | 0.004505 | 0.01371 |
| void | 0 | 0 | 0 | 0 | 0.000705 | 0 | 0.000245 | 0 | 0.001095 | 0 | 0.00141 | 0 |
| nodes inside other plate | 0.08234 | 0.01098 | 0.02187 | 0.004087 | 0.02008 | 0.006749 | 0.02155 | 0.006463 | 0.02829 | 0.009696 | 0.02681 | 0.003265 |
| stacked | 0 | 0 | 0.03994 | 0 | 0.07326 | 0 | 0.1616 | 0 | 0.2343 | 0 | 0.3775 | 0 |
| anisotropic | 0.000276 | 0.000287 | 0.003895 | 0.002155 | 0.004788 | 0.002919 | 0.007771 | 0.004001 | 0.02411 | 0.005641 | 0.03795 | 0.005493 |
| row alignment | 0.4286 | -0.4423 | -0.289 | 0.2856 | 0.1419 | 0.2076 | 0.1469 | 0.01081 | 0.3313 | 0.1307 | 0.3968 | 0.09789 |
| thin | 0.000123 | 9.57e-05 | 0.003715 | 0.003295 | 0.005207 | 0.003547 | 0.006722 | 0.006217 | 0.01679 | 0.006373 | 0.02445 | 0.008759 |
| one-node lines | 0.01008 | — | 0.09904 | — | 0.08882 | — | 0.1409 | — | 0.2394 | — | 0.3273 | — |
| quad one-cell-thin | — | 0.000734 | — | 0.006685 | — | 0.007785 | — | 0.01074 | — | 0.01448 | — | 0.01645 |
| air temperature mean (°C) | -3.378 | -3.453 | 0.5634 | -0.4713 | 2.76 | -0.8927 | 4.149 | -0.7687 | 5.716 | -11.03 | 7.747 | -17.17 |
| precipitation mean (mm) | 1132 | 1132 | 1168 | 1171 | 1181 | 1177 | 1184 | 1174 | 1207 | 1171 | 1252 | 1255 |
| river fraction of land | — | — | 0.1435 | 0.1585 | 0.1711 | 0.1809 | 0.1733 | 0.2084 | 0.1712 | 0.2705 | 0.178 | 0.2905 |
| HEALPix same node | 0.5526 | 0.5486 | 0.5474 | 0.551 | 0.5442 | 0.5522 | 0.5358 | 0.546 | 0.5275 | 0.5454 | 0.5117 | 0.5503 |

## Checkpoints, seed 3

| seed 3 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 18 | 18 | 21 | 21 | 26 | 23 | 32 | 28 | 50 | 35 | 53 | 47 |
| nodes | 32632 | 31213 | 32625 | 31545 | 32906 | 32194 | 35745 | 32820 | 38164 | 32645 | 41852 | 32915 |
| land fraction | 0.05561 | 0.05551 | 0.07864 | 0.07441 | 0.08236 | 0.07875 | 0.06659 | 0.06569 | 0.05865 | 0.05011 | 0.04181 | 0.02345 |
| sea level (m) | 0 | 0 | -536 | -488.1 | -755.6 | -552 | -1380 | -745.8 | -1853 | -925.5 | -2258 | -1123 |
| Hc volume (km³) | 7.2e+09 | 7.2e+09 | 6.86e+09 | 6.9e+09 | 6.5e+09 | 7.08e+09 | 5.39e+09 | 6.67e+09 | 4.02e+09 | 5.92e+09 | 3.31e+09 | 5.24e+09 |
| continental Hc (km³) | 2.59e+09 | 2.59e+09 | 2.74e+09 | 2.79e+09 | 2.62e+09 | 2.98e+09 | 2.01e+09 | 2.31e+09 | 1.04e+09 | 1.29e+09 | 7.09e+08 | 5.85e+08 |
| elevation p05/p50/p95 (m) | -5222/-4494/40.48 | -5222/-4493/41.41 | -5855/-5132/151.2 | -5532/-5119/37.66 | -6192/-5204/368.1 | -5739/-5164/404.8 | -6339/-5353/-549.4 | -5883/-5267/3.116 | -6448/-5512/-1404 | -5940/-5384/-925 | -6491/-5590/-2867 | -5933/-5426/-4157 |
| uncovered | 0.02706 | 0.00765 | 0.04494 | 0.008285 | 0.06031 | 0.009335 | 0.05826 | 0.01126 | 0.08195 | 0.01424 | 0.08327 | 0.01653 |
| multiply covered | 0.01106 | 0.00755 | 0.005425 | 0.01041 | 0.005095 | 0.00935 | 0.00666 | 0.01153 | 0.00392 | 0.0121 | 0.005385 | 0.01392 |
| void | 0 | 0 | 0.0015 | 0 | 0.00508 | 0 | 0.00068 | 0 | 0.000825 | 0 | 0.000845 | 0 |
| nodes inside other plate | 0.01492 | 0.01157 | 0.01655 | 0.005643 | 0.01501 | 0.004162 | 0.02171 | 0.006155 | 0.01462 | 0.002665 | 0.0233 | 0.00319 |
| stacked | 0 | 0 | 0.01116 | 0 | 0.05333 | 0 | 0.1811 | 0 | 0.3157 | 0 | 0.4086 | 0 |
| anisotropic | 0.000552 | 0.000192 | 0.00328 | 0.002092 | 0.005379 | 0.002609 | 0.01111 | 0.004022 | 0.01845 | 0.005269 | 0.07438 | 0.007686 |
| row alignment | 0.7907 | 0.3694 | -0.2428 | 0.08543 | -0.05305 | 0.08246 | -0.1601 | -0.03678 | 0.3072 | -0.04961 | 0.6981 | -0.06585 |
| thin | 6.13e-05 | 0.000192 | 0.003954 | 0.003012 | 0.009755 | 0.003199 | 0.01234 | 0.006277 | 0.01242 | 0.006984 | 0.03238 | 0.01082 |
| one-node lines | 0.006969 | — | 0.1034 | — | 0.1099 | — | 0.2123 | — | 0.2042 | — | 0.3192 | — |
| quad one-cell-thin | — | 0.000609 | — | 0.007703 | — | 0.008014 | — | 0.01164 | — | 0.01467 | — | 0.01975 |
| air temperature mean (°C) | -9.087 | -9.163 | -2.935 | -3.606 | -8.748 | -13.45 | -3.941 | -18.31 | 4.527 | -1.647 | -5.604 | -10.03 |
| precipitation mean (mm) | 1201 | 1201 | 1239 | 1235 | 1246 | 1233 | 1256 | 1232 | 1253 | 1238 | 1259 | 1249 |
| river fraction of land | — | — | 0.1475 | 0.1608 | 0.2383 | 0.2106 | 0.1902 | 0.222 | 0.1863 | 0.2406 | 0.2557 | 0.3509 |
| HEALPix same node | 0.5533 | 0.556 | 0.5556 | 0.5538 | 0.5491 | 0.5551 | 0.5359 | 0.5463 | 0.5233 | 0.5521 | 0.5106 | 0.5484 |

## Checkpoints, seed 4

| seed 4 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 17 | 17 | 19 | 19 | 24 | 20 | 29 | 28 | 39 | 37 | 48 | 48 |
| nodes | 32636 | 31134 | 32567 | 31356 | 33088 | 31746 | 34796 | 32677 | 35784 | 32935 | 38687 | 32632 |
| land fraction | 0.1417 | 0.1416 | 0.1611 | 0.1586 | 0.1672 | 0.1582 | 0.1701 | 0.1573 | 0.1486 | 0.1361 | 0.1281 | 0.0971 |
| sea level (m) | 0 | 0 | -561 | -534.1 | -714.9 | -526.8 | -1055 | -482.2 | -1351 | -762.7 | -1988 | -1063 |
| Hc volume (km³) | 8.49e+09 | 8.49e+09 | 8.17e+09 | 8.12e+09 | 7.96e+09 | 7.98e+09 | 7.52e+09 | 7.95e+09 | 6.14e+09 | 8.73e+09 | 5.17e+09 | 7.54e+09 |
| continental Hc (km³) | 4.13e+09 | 4.13e+09 | 4.24e+09 | 4.24e+09 | 4.31e+09 | 4.33e+09 | 4.34e+09 | 4.42e+09 | 3.24e+09 | 4.22e+09 | 2.43e+09 | 2.63e+09 |
| elevation p05/p50/p95 (m) | -5222/-4427/958.8 | -5222/-4427/962.3 | -5657/-5065/1142 | -5506/-5055/1055 | -6065/-5101/1183 | -5714/-5060/1095 | -6261/-5139/1279 | -5795/-5001/1198 | -6323/-4991/1299 | -5874/-5151/5647 | -6449/-5386/727.7 | -5906/-5279/3949 |
| uncovered | 0.02816 | 0.00778 | 0.03876 | 0.008305 | 0.04596 | 0.009085 | 0.06164 | 0.01104 | 0.07119 | 0.01297 | 0.07251 | 0.01656 |
| multiply covered | 0.01667 | 0.00729 | 0.002045 | 0.008235 | 0.005055 | 0.008715 | 0.01461 | 0.01192 | 0.002735 | 0.0155 | 0.01218 | 0.0146 |
| void | 0 | 0 | 0.00017 | 0 | 0.00027 | 0 | 0.00107 | 0 | 0.00099 | 0 | 0.00035 | 0 |
| nodes inside other plate | 0.01906 | 0.01092 | 0.007676 | 0.002232 | 0.01584 | 0.002992 | 0.02974 | 0.006916 | 0.01417 | 0.0109 | 0.03332 | 0.006098 |
| stacked | 0 | 0 | 0.006233 | 0 | 0.03962 | 0 | 0.1348 | 0 | 0.201 | 0 | 0.2911 | 0 |
| anisotropic | 0.00046 | 0.000225 | 0.003163 | 0.001818 | 0.008281 | 0.002551 | 0.01058 | 0.003703 | 0.0218 | 0.005374 | 0.03236 | 0.007048 |
| row alignment | 0.3941 | 0.03748 | -0.3379 | 0.01907 | -0.1497 | 0.1034 | 0.02997 | -0.05317 | 0.2699 | 0.04081 | 0.5594 | 0.06797 |
| thin | 0.000184 | 3.21e-05 | 0.004882 | 0.002775 | 0.008372 | 0.003938 | 0.009829 | 0.003795 | 0.01336 | 0.006953 | 0.01515 | 0.01039 |
| one-node lines | 0.009368 | — | 0.09032 | — | 0.1679 | — | 0.1766 | — | 0.2224 | — | 0.2353 | — |
| quad one-cell-thin | — | 0.000418 | — | 0.006889 | — | 0.008946 | — | 0.009885 | — | 0.0143 | — | 0.0175 |
| air temperature mean (°C) | -5.967 | -5.829 | -7.379 | -7.338 | -7.705 | -8.213 | -2.663 | -6.745 | 3.256 | -11.01 | 9.181 | -11.18 |
| precipitation mean (mm) | 1146 | 1147 | 1238 | 1254 | 1170 | 1183 | 1184 | 1172 | 1226 | 1156 | 1196 | 1153 |
| river fraction of land | — | — | 0.1637 | 0.1556 | 0.1628 | 0.1607 | 0.1538 | 0.1729 | 0.1629 | 0.192 | 0.1734 | 0.2502 |
| HEALPix same node | 0.5553 | 0.5487 | 0.5532 | 0.5488 | 0.5507 | 0.5496 | 0.5374 | 0.5483 | 0.5359 | 0.5515 | 0.5211 | 0.5524 |

## Checkpoints, seed 5

| seed 5 | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 16 | 16 | 17 | 17 | 19 | 18 | 24 | 25 | 44 | 32 | 59 | 42 |
| nodes | 32643 | 31550 | 32711 | 31544 | 32615 | 31535 | 33372 | 32220 | 36755 | 32969 | 39871 | 32629 |
| land fraction | 0.01461 | 0.01469 | 0.02159 | 0.01978 | 0.02522 | 0.0253 | 0.03185 | 0.03494 | 0.0503 | 0.036 | 0.06085 | 0.02246 |
| sea level (m) | 0 | 0 | -457.7 | -434.5 | -495.4 | -419.3 | -665.6 | -363.2 | -1120 | -185.4 | -1509 | -494.5 |
| Hc volume (km³) | 6.19e+09 | 6.19e+09 | 6.08e+09 | 6.03e+09 | 5.94e+09 | 6.04e+09 | 5.58e+09 | 6.09e+09 | 4.94e+09 | 6.25e+09 | 4.53e+09 | 5.56e+09 |
| continental Hc (km³) | 1.86e+09 | 1.86e+09 | 1.94e+09 | 1.95e+09 | 2.02e+09 | 2.03e+09 | 2.02e+09 | 2.13e+09 | 1.92e+09 | 1.96e+09 | 1.64e+09 | 1.16e+09 |
| elevation p05/p50/p95 (m) | -5222/-4570/-4197 | -5222/-4570/-4198 | -5532/-5183/-3610 | -5512/-5188/-3532 | -6011/-5187/-3178 | -5723/-5211/-2910 | -6309/-5240/-2554 | -5881/-5221/-1609 | -6382/-5345/-1107 | -5922/-5201/-411.6 | -6401/-5389/-1256 | -5928/-5329/-2538 |
| uncovered | 0.02357 | 0.006985 | 0.03218 | 0.00752 | 0.04061 | 0.00856 | 0.05264 | 0.01056 | 0.07284 | 0.01214 | 0.08543 | 0.01498 |
| multiply covered | 0.03266 | 0.00718 | 0.00512 | 0.008215 | 0.002615 | 0.00882 | 0.005555 | 0.01027 | 0.00297 | 0.01272 | 0.003715 | 0.01247 |
| void | 0 | 0 | 0 | 0 | 5e-06 | 0 | 0.00017 | 0 | 0.000155 | 0 | 0.00088 | 0 |
| nodes inside other plate | 0.0367 | 0.0109 | 0.01241 | 0.004216 | 0.008616 | 0.00463 | 0.01244 | 0.00478 | 0.01681 | 0.006794 | 0.01547 | 0.002912 |
| stacked | 0 | 0 | 0.001345 | 0 | 0.01263 | 0 | 0.06844 | 0 | 0.2377 | 0 | 0.3602 | 0 |
| anisotropic | 0.000429 | 9.51e-05 | 0.001529 | 0.002061 | 0.004293 | 0.002759 | 0.004195 | 0.003445 | 0.01861 | 0.003882 | 0.03812 | 0.006926 |
| row alignment | 0.3974 | -0.3597 | -0.05293 | 0.1791 | -0.1406 | -0.07852 | 0.2143 | -0.0458 | 0.1823 | 0.0551 | 0.5579 | 0.06734 |
| thin | 0.000123 | 6.34e-05 | 0.001498 | 0.003075 | 0.004906 | 0.00482 | 0.005064 | 0.005587 | 0.01521 | 0.00637 | 0.01706 | 0.008918 |
| one-node lines | 0.007371 | — | 0.04772 | — | 0.1111 | — | 0.09831 | — | 0.248 | — | 0.2341 | — |
| quad one-cell-thin | — | 0.000571 | — | 0.006689 | — | 0.008784 | — | 0.009994 | — | 0.01147 | — | 0.01652 |
| air temperature mean (°C) | -2.289 | -2.44 | -6.466 | -6.234 | -7.534 | -7.4 | -10.37 | -9.688 | -7.978 | -7.478 | -2.054 | -0.2478 |
| precipitation mean (mm) | 1254 | 1254 | 1252 | 1256 | 1263 | 1261 | 1265 | 1251 | 1266 | 1250 | 1255 | 1248 |
| river fraction of land | — | — | 0.1389 | 0.1282 | 0.1583 | 0.1478 | 0.1711 | 0.2645 | 0.2274 | 0.2757 | 0.2056 | 0.2265 |
| HEALPix same node | 0.5552 | 0.5447 | 0.5523 | 0.5491 | 0.552 | 0.5494 | 0.5476 | 0.5508 | 0.5316 | 0.5546 | 0.5161 | 0.55 |

## Performance (mean seconds per step)

| s/step | seed 1 lines | seed 1 quad | seed 2 lines | seed 2 quad | seed 3 lines | seed 3 quad | seed 4 lines | seed 4 quad | seed 5 lines | seed 5 quad |
|---|---|---|---|---|---|---|---|---|---|---|
| step_total | 3.089 | 2.334 | 3.992 | 2.533 | 4.184 | 2.703 | 4.636 | 3.149 | 4.79 | 3.308 |
| deform_topology | 1.46 | 0.7428 | 2.092 | 0.858 | 2.248 | 0.9376 | 2.347 | 1.089 | 2.46 | 1.113 |
| climate_erosion_hydrology | 0.6047 | 0.5535 | 0.6148 | 0.5451 | 0.6133 | 0.5597 | 0.7531 | 0.6611 | 0.7298 | 0.7001 |
| deform | 1.371 | 0.3636 | 1.975 | 0.41 | 2.135 | 0.4398 | 2.205 | 0.5068 | 2.321 | 0.5198 |
| faults | 0.505 | 0.668 | 0.679 | 0.7302 | 0.6964 | 0.8083 | 0.7843 | 0.9126 | 0.858 | 1.014 |
| fluid_dynamics | 2.19e-06 | 2.02e-06 | 2.49e-06 | 1.66e-06 | 3.02e-06 | 2.16e-06 | 4.63e-06 | 3.53e-06 | 5.59e-06 | 3.36e-06 |
| gap_fill | 0.01217 | 0.312 | 0.0128 | 0.3753 | 0.01281 | 0.4174 | 0.01562 | 0.4879 | 0.01552 | 0.4946 |
| magma_transport | 0.1229 | 0.1097 | 0.1012 | 0.1103 | 0.0826 | 0.08517 | 0.1721 | 0.1305 | 0.1124 | 0.09322 |
| overlap_tracking | 0.01724 | 0.03259 | 0.02136 | 0.03768 | 0.0223 | 0.04079 | 0.02513 | 0.05027 | 0.02609 | 0.05165 |
| record_stats | 0.006362 | 0.0034 | 0.007602 | 0.00331 | 0.007779 | 0.003399 | 0.00914 | 0.004055 | 0.01 | 0.00421 |
| resource_formation | 0.01913 | 0.008705 | 0.02412 | 0.00892 | 0.02529 | 0.008587 | 0.02749 | 0.01073 | 0.02902 | 0.01163 |
| sea_level | 0.05552 | 0.04596 | 0.07181 | 0.05377 | 0.08042 | 0.06094 | 0.08229 | 0.0619 | 0.09464 | 0.07577 |
| shift | 0.2816 | 0.1923 | 0.3462 | 0.2121 | 0.3717 | 0.227 | 0.398 | 0.2658 | 0.432 | 0.2823 |
| stranded_basins | 1.25e-05 | 2.14e-05 | 1.46e-05 | 1.47e-05 | 9.65e-06 | 1.28e-05 | 3.71e-05 | 2.41e-05 | 3.56e-05 | 1.48e-05 |
| topology | 0.05957 | 0.03462 | 0.08366 | 0.03495 | 0.07836 | 0.03962 | 0.1007 | 0.04371 | 0.09745 | 0.04675 |
| volcanism | 0.0313 | 0.008425 | 0.04993 | 0.00985 | 0.0539 | 0.01119 | 0.05706 | 0.01239 | 0.05865 | 0.01207 |

Index builds per step (mean):

- seed 1 lines: line_row_lookup 73.5, plate_node_kdtree 73.6, plate_outline 79.2, plate_outline_kdtree 78.0, world_node_kdtree 1.0
- seed 1 quad: plate_node_kdtree 72.7, plate_outline 82.1, plate_outline_kdtree 75.6, quad_adjacency 43.8, quad_boundary_loops 56.1, world_node_kdtree 1.0
- seed 2 lines: line_row_lookup 105.4, plate_node_kdtree 105.7, plate_outline 114.5, plate_outline_kdtree 111.9, world_node_kdtree 1.0
- seed 2 quad: plate_node_kdtree 92.4, plate_outline 104.9, plate_outline_kdtree 95.5, quad_adjacency 54.2, quad_boundary_loops 71.8, world_node_kdtree 1.0
- seed 3 lines: line_row_lookup 113.0, plate_node_kdtree 113.3, plate_outline 122.8, plate_outline_kdtree 120.4, world_node_kdtree 1.0
- seed 3 quad: plate_node_kdtree 92.9, plate_outline 105.5, plate_outline_kdtree 95.7, quad_adjacency 54.6, quad_boundary_loops 72.9, world_node_kdtree 1.0
- seed 4 lines: line_row_lookup 91.4, plate_node_kdtree 91.6, plate_outline 100.4, plate_outline_kdtree 97.1, world_node_kdtree 1.0
- seed 4 quad: plate_node_kdtree 94.5, plate_outline 107.7, plate_outline_kdtree 98.8, quad_adjacency 55.9, quad_boundary_loops 73.9, world_node_kdtree 1.0
- seed 5 lines: line_row_lookup 100.6, plate_node_kdtree 100.8, plate_outline 110.5, plate_outline_kdtree 108.0, world_node_kdtree 1.0
- seed 5 quad: plate_node_kdtree 82.4, plate_outline 93.5, plate_outline_kdtree 86.3, quad_adjacency 46.3, quad_boundary_loops 63.7, world_node_kdtree 1.0

