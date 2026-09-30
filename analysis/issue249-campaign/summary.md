# Parity campaign summary

| set | runs | verdict |
|---|---:|---|
| long | 10 | fail |
| repro-equivalent | 2 | fail |
| repro-lines | 1 | pass |
| stress | 4 | fail |

## Gate status by set

| gate | long | repro-equivalent | repro-lines | stress |
|---|---|---|---|---|
| H1:quad_topology | pass | pass | — | pass |
| H2:quad_neighbours | pass | pass | — | pass |
| H3:quad_folded | pass | pass | — | pass |
| H4:fields | pass | pass | pass | pass |
| H5:frames | pass | pass | pass | pass |
| H6:derived_caches | pass | pass | pass | pass |
| H7:revisions | pass | pass | pass | pass |
| H8:load_round_trip | pass | pass | pass | pass |
| H9:load_continuation | pass | pass | pass | pass |
| H10:quad_no_stacked_nodes | pass | pass | — | pass |
| H11:field_caps | warn (+5 line) | warn (+1 line) | info (+1 line) | warn (+2 line) |
| C1:uncovered | pass | pass | — | pass |
| C2:void | pass | pass | — | pass |
| C3:multiply_covered | fail | fail | — | fail |
| C4:nodes_inside_other_plate | pass | pass | — | pass |
| K1:hc_volume_km3_drift | fail | warn | — | warn |
| K2:continental_hc_volume_km3_drift | fail | fail | — | fail |
| K3:quad_area_accounting | pass | pass | — | pass |
| K4:hm_volume_km3_drift | fail | fail | — | fail |
| M1:anisotropic_fraction | pass | pass | — | warn |
| M2:row_alignment | pass | pass | — | pass |
| M3:thin_fraction | warn | pass | — | warn |
| M4:aspect_gt_4 | pass | pass | — | pass |
| M5:skew_gt_45 | pass | pass | — | pass |
| M6:fragments | warn | warn | — | warn |
| S1:climate_hydrology_finite | pass | pass | pass | pass |
| S2:stability_air_temperature_mean_c | warn | warn | — | warn |
| S2:stability_land_fraction | warn | pass | — | pass |
| S2:stability_ocean_temperature_mean_c | pass | pass | — | pass |
| S2:stability_precipitation_mean_mm | pass | pass | — | pass |
| S2:stability_sea_level_m | warn | warn | — | warn |
| X1:healpix_all | pass | pass | — | pass |
| X1:healpix_antimeridian | pass | pass | — | pass |
| X1:healpix_coast | pass | pass | — | warn |
| X1:healpix_hole | warn | insufficient | — | insufficient |
| X1:healpix_plate_boundary | pass | pass | — | pass |
| X1:healpix_pole | pass | pass | — | pass |
| R1:deform_topology_s_per_step | pass | pass | — | pass |
| R2:step_total_s_per_step | pass | pass | — | pass |
| P1:continental_hc_volume_km3 | pass | insufficient | — | insufficient |
| P1:elevation_p05 | fail | insufficient | — | insufficient |
| P1:elevation_p50 | pass | insufficient | — | insufficient |
| P1:elevation_p95 | pass | insufficient | — | insufficient |
| P1:hc_volume_km3 | pass | insufficient | — | insufficient |
| P1:land_fraction | pass | insufficient | — | insufficient |
| P1:plates | pass | insufficient | — | insufficient |
| P1:sea_level_m | fail | insufficient | — | insufficient |

## long: seed means (seeds 1, 2, 3, 4, 5)

| metric | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 16.6 | 16.6 | 18.4 | 19 | 21 | 20.8 | 26.4 | 26.2 | 41.2 | 34.4 | 50.8 | 42 |
| land fraction | 0.1264 | 0.1263 | 0.1453 | 0.1441 | 0.1482 | 0.147 | 0.1463 | 0.1378 | 0.1463 | 0.1057 | 0.1366 | 0.05994 |
| sea level (m) | 0 | 0 | -598.8 | -494.7 | -752.1 | -666.4 | -1120 | -848.2 | -1696 | -1094 | -2058 | -1395 |
| Hc drift since 0 Myr | 0 | 0 | -0.04081 | -0.0443 | -0.07301 | -0.05343 | -0.1629 | -0.1124 | -0.2935 | -0.1991 | -0.4129 | -0.3097 |
| uncovered | 0.02541 | 0.00719 | 0.03565 | 0.01779 | 0.04414 | 0.008192 | 0.05572 | 0.009805 | 0.07232 | 0.01201 | 0.07981 | 0.01294 |
| multiply covered | 0.03456 | 0.007152 | 0.00551 | 0.004756 | 0.00675 | 0.009043 | 0.005648 | 0.01105 | 0.004217 | 0.01265 | 0.004809 | 0.01371 |
| nodes inside other plate | 0.03384 | 0.0108 | 0.01318 | 0.002808 | 0.01418 | 0.004278 | 0.01607 | 0.005745 | 0.01703 | 0.005929 | 0.01948 | 0.00491 |
| stacked | 0 | 0 | 0.01093 | 0 | 0.03773 | 0 | 0.1238 | 0 | 0.2687 | 0 | 0.3574 | 0 |
| anisotropic | 0.000447 | 0.000192 | 0.00292 | 0.00141 | 0.004978 | 0.001913 | 0.007954 | 0.002727 | 0.02259 | 0.003903 | 0.04073 | 0.00529 |
| thin | 0.000135 | 0.000115 | 0.003105 | 0.00218 | 0.005424 | 0.002776 | 0.007932 | 0.003661 | 0.01361 | 0.004999 | 0.01734 | 0.007014 |
| one-node lines | 0.008657 | — | 0.08019 | — | 0.1183 | — | 0.1495 | — | 0.223 | — | 0.2651 | — |
| quad one-cell-thin | — | 0.000542 | — | 0.005428 | — | 0.007006 | — | 0.008375 | — | 0.0112 | — | 0.0135 |
| quad hole loops | — | 0 | — | 0 | — | 0.4 | — | 0.8 | — | 0.8 | — | 2.2 |
| air temperature (°C) | -3.246 | -3.281 | -1.619 | -2.222 | -2.822 | -3.358 | -1.16 | -1.047 | 0.179 | 2.851 | 1.453 | 7.867 |
| precipitation (mm) | 1144 | 1144 | 1181 | 1183 | 1174 | 1170 | 1188 | 1188 | 1191 | 1206 | 1216 | 1240 |
| river fraction of land | — | — | 0.1393 | 0.147 | 0.1646 | 0.1677 | 0.1661 | 0.1919 | 0.1817 | 0.1891 | 0.1835 | 0.1675 |
| HEALPix same node | 0.5547 | 0.5499 | 0.5522 | 0.5504 | 0.5507 | 0.5515 | 0.5415 | 0.549 | 0.529 | 0.5499 | 0.5168 | 0.55 |
| HEALPix p95 distance (s) | 1.18 | 1.201 | 1.203 | 1.212 | 1.213 | 1.196 | 1.224 | 1.197 | 1.226 | 1.196 | 1.234 | 1.198 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 4.736 | 2.882 | 0.6085 |
| deform_topology | 2.417 | 0.7922 | 0.3277 |
| deform | 2.285 | 0.5348 | 0.234 |
| topology | 0.09125 | 0.04422 | 0.4846 |
| gap_fill | 0.01557 | 0.1674 | 10.75 |
| overlap_tracking | 0.02494 | 0.04581 | 1.837 |
| faults | 0.7951 | 0.9042 | 1.137 |
| shift | 0.4218 | 0.2661 | 0.6308 |
| climate_erosion_hydrology | 0.7799 | 0.7092 | 0.9094 |
| magma_transport | 0.146 | 0.1216 | 0.8325 |

## repro-equivalent: seed means (seeds 804913535)

| metric | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 352 Myr lines | 352 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|
| plates | 9 | 9 | 11 | 11 | 21 | 21 | 32 | 28 | 37 | 31 |
| land fraction | 0.2737 | 0.2741 | 0.3132 | 0.3095 | 0.2676 | 0.2173 | 0.1965 | 0.129 | 0.1946 | 0.0998 |
| sea level (m) | 0 | 0 | -616.5 | -526.1 | -1289 | -1352 | -2156 | -1884 | -2632 | -2135 |
| Hc drift since 0 Myr | 0 | 0 | -0.02313 | -0.02533 | -0.1334 | -0.1864 | -0.3858 | -0.4007 | -0.5003 | -0.4567 |
| uncovered | 0.0202 | 0.00533 | 0.0233 | 0.01206 | 0.04532 | 0.008345 | 0.07395 | 0.0108 | 0.06213 | 0.0126 |
| multiply covered | 0.01241 | 0.00524 | 0.004405 | 0.00508 | 0.007875 | 0.01162 | 0.003755 | 0.01209 | 0.00247 | 0.01197 |
| nodes inside other plate | 0.01452 | 0.007697 | 0.0127 | 0.005894 | 0.02094 | 0.01038 | 0.01751 | 0.005662 | 0.01213 | 0.005218 |
| stacked | 0 | 0 | 0.004203 | 0 | 0.06922 | 0 | 0.2102 | 0 | 0.3153 | 0 |
| anisotropic | 0.000306 | 9.39e-05 | 0.001767 | 0.000784 | 0.004421 | 0.002617 | 0.02327 | 0.003734 | 0.04686 | 0.003125 |
| thin | 6.13e-05 | 0 | 0.001828 | 0.001035 | 0.0046 | 0.002769 | 0.01497 | 0.005111 | 0.02425 | 0.004247 |
| one-node lines | 0 | — | 0.05845 | — | 0.1028 | — | 0.2794 | — | 0.2224 | — |
| quad one-cell-thin | — | 0.000219 | — | 0.003386 | — | 0.007852 | — | 0.01108 | — | 0.01077 |
| quad hole loops | — | 0 | — | 0 | — | 0 | — | 0 | — | 1 |
| air temperature (°C) | 8.022 | 8.051 | 5.29 | 4.991 | 0.376 | 0.5143 | 1.538 | 2.254 | -10.65 | -7.998 |
| precipitation (mm) | 1004 | 1004 | 1035 | 1042 | 1068 | 1128 | 1184 | 1213 | 1151 | 1224 |
| river fraction of land | — | — | 0.1103 | 0.1068 | 0.1735 | 0.1837 | 0.1582 | 0.2295 | 0.1807 | 0.1868 |
| HEALPix same node | 0.5548 | 0.5569 | 0.5571 | 0.5524 | 0.5491 | 0.5531 | 0.5327 | 0.5494 | 0.5272 | 0.5521 |
| HEALPix p95 distance (s) | 1.176 | 1.177 | 1.179 | 1.184 | 1.225 | 1.167 | 1.279 | 1.184 | 1.199 | 1.19 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 4.361 | 2.871 | 0.6583 |
| deform_topology | 2.133 | 0.7623 | 0.3573 |
| deform | 1.996 | 0.5243 | 0.2626 |
| topology | 0.09874 | 0.04925 | 0.4988 |
| gap_fill | 0.01577 | 0.1473 | 9.34 |
| overlap_tracking | 0.02248 | 0.04146 | 1.844 |
| faults | 0.7016 | 0.8799 | 1.254 |
| shift | 0.3933 | 0.2633 | 0.6696 |
| climate_erosion_hydrology | 0.838 | 0.79 | 0.9427 |
| magma_transport | 0.1383 | 0.08812 | 0.6371 |

## repro-lines: seed means (seeds 804913535)

| metric | 352.4 Myr lines | 352.9 Myr lines | 353.4 Myr lines |
|---|---|---|---|
| plates | 35 | 35 | 35 |
| land fraction | 0.3892 | 0.3924 | 0.3983 |
| sea level (m) | -826.2 | -831.5 | -837.7 |
| Hc drift since 0 Myr | 0 | -0.002385 | -0.008671 |
| uncovered | 0.0232 | 0.02341 | 0.02261 |
| multiply covered | 0.004735 | 0.00392 | 0.004505 |
| nodes inside other plate | 0.01379 | 0.01257 | 0.01383 |
| stacked | 0.275 | 0.2761 | 0.2767 |
| anisotropic | 0.04662 | 0.04726 | 0.04814 |
| thin | 0.05789 | 0.05775 | 0.05702 |
| one-node lines | 0.293 | 0.3036 | 0.3092 |
| quad one-cell-thin | — | — | — |
| quad hole loops | — | — | — |
| air temperature (°C) | -7.643 | -7.687 | -7.675 |
| precipitation (mm) | 1021 | 1027 | 1027 |
| river fraction of land | 0.1134 | 0.1131 | 0.1132 |
| HEALPix same node | 0.5346 | 0.5358 | 0.534 |
| HEALPix p95 distance (s) | 1.174 | 1.17 | 1.173 |

## stress: seed means (seeds 1, 2)

| metric | 0 Myr lines | 0 Myr quad | 100 Myr lines | 100 Myr quad | 200 Myr lines | 200 Myr quad | 400 Myr lines | 400 Myr quad | 600 Myr lines | 600 Myr quad | 800 Myr lines | 800 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 16 | 16 | 21.5 | 21 | 37 | 30.5 | 51 | 40 | 56.5 | 42 | 61 | 44.5 |
| land fraction | 0.2093 | 0.2092 | 0.2147 | 0.2246 | 0.2063 | 0.1769 | 0.1568 | 0.08015 | 0.1515 | 0.05152 | 0.1786 | 0.02878 |
| sea level (m) | 0 | 0 | -1015 | -817.7 | -1539 | -1107 | -2073 | -1737 | -2447 | -1944 | -2571 | -1978 |
| Hc drift since 0 Myr | 0 | 0 | -0.1081 | -0.08868 | -0.2231 | -0.1963 | -0.3948 | -0.3971 | -0.5186 | -0.4696 | -0.5248 | -0.5098 |
| uncovered | 0.03319 | 0.009175 | 0.0594 | 0.01222 | 0.09238 | 0.01639 | 0.1089 | 0.01937 | 0.1218 | 0.0197 | 0.1281 | 0.02039 |
| multiply covered | 0.05803 | 0.01059 | 0.0125 | 0.0142 | 0.00795 | 0.01854 | 0.00581 | 0.01937 | 0.006315 | 0.0184 | 0.006405 | 0.01918 |
| nodes inside other plate | 0.05264 | 0.01622 | 0.02209 | 0.007934 | 0.02637 | 0.01086 | 0.03137 | 0.007781 | 0.03234 | 0.004391 | 0.04168 | 0.004373 |
| stacked | 0 | 0 | 0.06936 | 0 | 0.2093 | 0 | 0.292 | 0 | 0.4016 | 0 | 0.4527 | 0 |
| anisotropic | 0.00095 | 0.00029 | 0.009322 | 0.003715 | 0.01846 | 0.007255 | 0.04917 | 0.009218 | 0.06396 | 0.009889 | 0.08615 | 0.01159 |
| thin | 0.000337 | 0.000354 | 0.009234 | 0.006437 | 0.01622 | 0.009184 | 0.03138 | 0.01167 | 0.03355 | 0.01253 | 0.04975 | 0.01594 |
| one-node lines | 0.01447 | — | 0.1592 | — | 0.211 | — | 0.3571 | — | 0.3732 | — | 0.4229 | — |
| quad one-cell-thin | — | 0.000966 | — | 0.01174 | — | 0.01746 | — | 0.02203 | — | 0.02282 | — | 0.02605 |
| quad hole loops | — | 0 | — | 0 | — | 2.5 | — | 3.5 | — | 1 | — | 5 |
| air temperature (°C) | 0.5321 | 0.4815 | 3.78 | 5.146 | 1.695 | 4.846 | 2.786 | 6.877 | 5.667 | 1.258 | 9.254 | 17.57 |
| precipitation (mm) | 1060 | 1061 | 1088 | 1104 | 1133 | 1140 | 1200 | 1223 | 1234 | 1259 | 1242 | 1263 |
| river fraction of land | — | — | 0.1944 | 0.1856 | 0.2217 | 0.2071 | 0.2001 | 0.2299 | 0.228 | 0.592 | 0.1738 | 0.05609 |
| HEALPix same node | 0.5007 | 0.5048 | 0.4902 | 0.4996 | 0.4693 | 0.4951 | 0.4599 | 0.4975 | 0.442 | 0.4978 | 0.433 | 0.496 |
| HEALPix p95 distance (s) | 1.077 | 1.119 | 1.12 | 1.11 | 1.161 | 1.106 | 1.166 | 1.106 | 1.172 | 1.09 | 1.188 | 1.119 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 3.526 | 1.937 | 0.5493 |
| deform_topology | 1.732 | 0.4855 | 0.2804 |
| deform | 1.629 | 0.3167 | 0.1944 |
| topology | 0.07752 | 0.03024 | 0.3901 |
| gap_fill | 0.007788 | 0.1065 | 13.68 |
| overlap_tracking | 0.01716 | 0.03209 | 1.87 |
| faults | 0.7066 | 0.7091 | 1.003 |
| shift | 0.2938 | 0.1442 | 0.4909 |
| climate_erosion_hydrology | 0.6172 | 0.5214 | 0.8447 |
| magma_transport | 0.06944 | 0.03496 | 0.5034 |
