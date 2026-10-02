# Parity campaign summary

| set | runs | verdict |
|---|---:|---|
| issue147 | 2 | fail |
| long | 10 | fail |
| repro-equivalent | 2 | fail |
| repro-lines | 1 | pass |
| stress | 4 | fail |

## Gate status by set

| gate | issue147 | long | repro-equivalent | repro-lines | stress |
|---|---|---|---|---|---|
| H1:quad_topology | pass | pass | pass | — | pass |
| H2:quad_neighbours | pass | pass | pass | — | pass |
| H3:quad_folded | pass | pass | pass | — | pass |
| H4:fields | pass | pass | pass | pass | pass |
| H5:frames | pass | pass | pass | pass | pass |
| H6:derived_caches | pass | pass | pass | pass | pass |
| H7:revisions | pass | pass | pass | pass | pass |
| H8:load_round_trip | pass | pass | pass | pass | pass |
| H9:load_continuation | pass | pass | pass | pass | pass |
| H10:quad_no_stacked_nodes | pass | pass | pass | — | pass |
| H11:field_caps | pass | pass | pass | info (+1 line) | pass |
| C1:uncovered | pass | pass | pass | — | pass |
| C2:void | pass | pass | pass | — | pass |
| C3:multiply_covered | pass | warn | pass | — | fail |
| C4:nodes_inside_other_plate | pass | pass | pass | — | pass |
| K1:hc_volume_km3_drift | pass | fail | fail | — | fail |
| K2:continental_hc_volume_km3_drift | pass | fail | warn | — | fail |
| K3:quad_area_accounting | pass | pass | pass | — | pass |
| K4:hm_volume_km3_drift | fail | fail | warn | — | fail |
| M1:anisotropic_fraction | pass | pass | pass | — | warn |
| M2:row_alignment | pass | pass | pass | — | pass |
| M3:thin_fraction | warn | warn | pass | — | warn |
| M4:aspect_gt_4 | pass | pass | pass | — | pass |
| M5:skew_gt_45 | pass | pass | pass | — | pass |
| M6:fragments | warn | warn | warn | — | warn |
| S1:climate_hydrology_finite | pass | pass | pass | pass | pass |
| S2:stability_air_temperature_mean_c | pass | pass | pass | — | pass |
| S2:stability_land_fraction | pass | pass | pass | — | pass |
| S2:stability_ocean_temperature_mean_c | pass | pass | pass | — | pass |
| S2:stability_precipitation_mean_mm | pass | pass | pass | — | pass |
| S2:stability_sea_level_m | pass | pass | pass | — | pass |
| X1:healpix_all | pass | pass | pass | — | pass |
| X1:healpix_antimeridian | pass | pass | pass | — | pass |
| X1:healpix_coast | pass | pass | pass | — | pass |
| X1:healpix_hole | insufficient | warn | insufficient | — | insufficient |
| X1:healpix_plate_boundary | pass | pass | pass | — | pass |
| X1:healpix_pole | pass | pass | pass | — | pass |
| R1:deform_topology_s_per_step | pass | pass | pass | — | pass |
| R2:step_total_s_per_step | pass | pass | pass | — | pass |
| P1:continental_hc_volume_km3 | insufficient | pass | insufficient | — | insufficient |
| P1:elevation_p05 | insufficient | fail | insufficient | — | insufficient |
| P1:elevation_p50 | insufficient | pass | insufficient | — | insufficient |
| P1:elevation_p95 | insufficient | pass | insufficient | — | insufficient |
| P1:hc_volume_km3 | insufficient | pass | insufficient | — | insufficient |
| P1:land_fraction | insufficient | pass | insufficient | — | insufficient |
| P1:plates | insufficient | pass | insufficient | — | insufficient |
| P1:sea_level_m | insufficient | fail | insufficient | — | insufficient |

## issue147: seed means (seeds 0)

| metric | 0 Myr lines | 0 Myr quad | 10 Myr lines | 10 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad |
|---|---|---|---|---|---|---|---|---|
| plates | 19 | 19 | 28 | 27 | 31 | 29 | 38 | 30 |
| land fraction | 0.153 | 0.153 | 0.1737 | 0.1627 | 0.1677 | 0.1543 | 0.156 | 0.1456 |
| sea level (m) | 0 | 0 | -373.9 | -384.2 | -695.5 | -542 | -865.9 | -547.3 |
| Hc drift since 0 Myr | 0 | 0 | -0.0102 | -0.0272 | 0.01493 | -0.01741 | 0.02745 | 0.001031 |
| uncovered | 0.01328 | 0.003935 | 0.02469 | 0.0058 | 0.02671 | 0.0064 | 0.03074 | 0.00723 |
| multiply covered | 0.01601 | 0.0038 | 0.00614 | 0.005195 | 0.001435 | 0.00658 | 0.001725 | 0.00611 |
| nodes inside other plate | 0.01771 | 0.005561 | 0.007473 | 0.002548 | 0.004076 | 0.003599 | 0.005954 | 0.002176 |
| stacked | 0 | 0 | 0.008937 | 0 | 0.07282 | 0 | 0.1344 | 0 |
| anisotropic | 0.000145 | 3.19e-05 | 0.001019 | 0.001006 | 0.004386 | 0.001867 | 0.0137 | 0.002263 |
| thin | 3.06e-05 | 2.39e-05 | 0.001226 | 0.002476 | 0.002352 | 0.003503 | 0.003608 | 0.004058 |
| one-node lines | 0.002784 | — | 0.0404 | — | 0.08951 | — | 0.152 | — |
| quad one-cell-thin | — | 0.00016 | — | 0.005144 | — | 0.006535 | — | 0.007679 |
| quad hole loops | — | 0 | — | 7 | — | 13 | — | 6 |
| air temperature (°C) | 5.78 | 5.798 | 3.999 | 4.211 | 0.3401 | -0.2379 | -2.967 | -5.288 |
| precipitation (mm) | 1096 | 1096 | 1129 | 1137 | 1126 | 1140 | 1153 | 1160 |
| river fraction of land | — | — | 0.08387 | 0.09501 | 0.1018 | 0.1178 | 0.1135 | 0.1325 |
| HEALPix same node | 0.5559 | 0.5458 | 0.5557 | 0.5471 | 0.5477 | 0.5443 | 0.5414 | 0.5463 |
| HEALPix p95 distance (s) | 1.183 | 1.211 | 1.196 | 1.214 | 1.207 | 1.22 | 1.205 | 1.227 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 8.451 | 3.767 | 0.4457 |
| deform_topology | 1.562 | 1.364 | 0.8735 |
| deform | 1.339 | 0.5454 | 0.4073 |
| topology | 0.1416 | 0.08523 | 0.602 |
| gap_fill | 0.03938 | 0.6647 | 16.88 |
| overlap_tracking | 0.04148 | 0.06877 | 1.658 |
| faults | 0.628 | 0.6846 | 1.09 |
| shift | 0.4021 | 0.3149 | 0.7832 |
| climate_erosion_hydrology | 0.8317 | 0.742 | 0.8921 |
| magma_transport | 4.721 | 0.4461 | 0.0945 |

## long: seed means (seeds 1, 2, 3, 4, 5)

| metric | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 16.6 | 16.6 | 18.6 | 18.4 | 22 | 19.8 | 27.2 | 25.4 | 41.8 | 33.4 | 53.2 | 44.6 |
| land fraction | 0.1264 | 0.1263 | 0.1454 | 0.1422 | 0.1445 | 0.1426 | 0.144 | 0.1368 | 0.1463 | 0.1145 | 0.136 | 0.07684 |
| sea level (m) | 0 | 0 | -609.4 | -542.6 | -762 | -577 | -1111 | -615 | -1496 | -792.9 | -2006 | -1120 |
| Hc drift since 0 Myr | 0 | 0 | -0.04228 | -0.04037 | -0.08569 | -0.04088 | -0.1635 | -0.05289 | -0.2926 | -0.06536 | -0.4194 | -0.1998 |
| uncovered | 0.02541 | 0.00719 | 0.03561 | 0.007857 | 0.04627 | 0.00872 | 0.05578 | 0.01057 | 0.0703 | 0.01268 | 0.08179 | 0.0152 |
| multiply covered | 0.03456 | 0.007152 | 0.00551 | 0.008658 | 0.006903 | 0.009018 | 0.008058 | 0.01126 | 0.04713 | 0.01393 | 0.02665 | 0.01378 |
| nodes inside other plate | 0.03384 | 0.0108 | 0.01313 | 0.004077 | 0.01517 | 0.004626 | 0.01952 | 0.006579 | 0.04192 | 0.008173 | 0.03394 | 0.004713 |
| stacked | 0 | 0 | 0.01182 | 0 | 0.03673 | 0 | 0.121 | 0 | 0.2413 | 0 | 0.352 | 0 |
| anisotropic | 0.000447 | 0.000192 | 0.002808 | 0.002004 | 0.005508 | 0.002562 | 0.007902 | 0.003523 | 0.01913 | 0.005029 | 0.0424 | 0.006195 |
| thin | 0.000135 | 0.000115 | 0.003251 | 0.003165 | 0.006503 | 0.003895 | 0.008065 | 0.005041 | 0.01394 | 0.006622 | 0.02147 | 0.008883 |
| one-node lines | 0.008657 | — | 0.08113 | — | 0.1205 | — | 0.1621 | — | 0.2137 | — | 0.2855 | — |
| quad one-cell-thin | — | 0.000542 | — | 0.007035 | — | 0.008188 | — | 0.01009 | — | 0.01313 | — | 0.01651 |
| quad hole loops | — | 0 | — | 0.4 | — | 1.4 | — | 1 | — | 3.6 | — | 5.8 |
| air temperature (°C) | -3.246 | -3.281 | -1.717 | -2.226 | -2.665 | -4.814 | -1.063 | -6.839 | 1.925 | -6.412 | 3.363 | -8.287 |
| precipitation (mm) | 1144 | 1144 | 1182 | 1185 | 1174 | 1170 | 1185 | 1175 | 1203 | 1182 | 1220 | 1209 |
| river fraction of land | — | — | 0.1395 | 0.1435 | 0.1668 | 0.1664 | 0.1636 | 0.2057 | 0.18 | 0.2324 | 0.1942 | 0.2672 |
| HEALPix same node | 0.5547 | 0.5499 | 0.5526 | 0.5507 | 0.55 | 0.5511 | 0.5409 | 0.5483 | 0.5303 | 0.551 | 0.5162 | 0.5494 |
| HEALPix p95 distance (s) | 1.18 | 1.201 | 1.2 | 1.198 | 1.219 | 1.195 | 1.224 | 1.199 | 1.233 | 1.187 | 1.238 | 1.195 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 4.138 | 2.805 | 0.6779 |
| deform_topology | 2.121 | 0.948 | 0.4469 |
| deform | 2.001 | 0.448 | 0.2239 |
| topology | 0.08395 | 0.03993 | 0.4757 |
| gap_fill | 0.01378 | 0.4175 | 30.29 |
| overlap_tracking | 0.02242 | 0.0426 | 1.9 |
| faults | 0.7046 | 0.8266 | 1.173 |
| shift | 0.3659 | 0.2359 | 0.6448 |
| climate_erosion_hydrology | 0.6631 | 0.6039 | 0.9107 |
| magma_transport | 0.1182 | 0.1058 | 0.8946 |

## repro-equivalent: seed means (seeds 804913535)

| metric | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 352 Myr lines | 352 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|
| plates | 9 | 9 | 11 | 11 | 21 | 25 | 32 | 34 | 45 | 42 |
| land fraction | 0.2737 | 0.2741 | 0.3143 | 0.3078 | 0.2703 | 0.2299 | 0.1804 | 0.1481 | 0.1211 | 0.09691 |
| sea level (m) | 0 | 0 | -590.6 | -545 | -1379 | -1163 | -2429 | -1658 | -2760 | -2013 |
| Hc drift since 0 Myr | 0 | 0 | -0.02386 | -0.01836 | -0.1319 | -0.02366 | -0.4189 | -0.2535 | -0.5519 | -0.3503 |
| uncovered | 0.0202 | 0.00533 | 0.02459 | 0.005855 | 0.04796 | 0.01003 | 0.06039 | 0.01336 | 0.06531 | 0.01549 |
| multiply covered | 0.01241 | 0.00524 | 0.004545 | 0.00793 | 0.01961 | 0.01206 | 0.00424 | 0.01349 | 0.005185 | 0.01342 |
| nodes inside other plate | 0.01452 | 0.007697 | 0.01404 | 0.00726 | 0.0427 | 0.009476 | 0.01467 | 0.006405 | 0.02196 | 0.005069 |
| stacked | 0 | 0 | 0.003967 | 0 | 0.08112 | 0 | 0.2461 | 0 | 0.2436 | 0 |
| anisotropic | 0.000306 | 9.39e-05 | 0.002044 | 0.001582 | 0.004858 | 0.002369 | 0.01594 | 0.00613 | 0.03052 | 0.006307 |
| thin | 6.13e-05 | 0 | 0.004516 | 0.001831 | 0.007053 | 0.003614 | 0.01021 | 0.008265 | 0.01318 | 0.007996 |
| one-node lines | 0 | — | 0.06045 | — | 0.08988 | — | 0.1823 | — | 0.2408 | — |
| quad one-cell-thin | — | 0.000219 | — | 0.004995 | — | 0.008959 | — | 0.0133 | — | 0.01629 |
| quad hole loops | — | 0 | — | 0 | — | 1 | — | 5 | — | 3 |
| air temperature (°C) | 8.022 | 8.051 | 5.354 | 5.126 | 0.7672 | -8.147 | -1.416 | -6.035 | -12.73 | -14.43 |
| precipitation (mm) | 1004 | 1004 | 1031 | 1046 | 1093 | 1091 | 1199 | 1170 | 1215 | 1174 |
| river fraction of land | — | — | 0.1078 | 0.1119 | 0.1358 | 0.1992 | 0.1926 | 0.2466 | 0.2232 | 0.3145 |
| HEALPix same node | 0.5548 | 0.5569 | 0.5588 | 0.5545 | 0.5447 | 0.5529 | 0.53 | 0.5537 | 0.5327 | 0.5514 |
| HEALPix p95 distance (s) | 1.176 | 1.177 | 1.178 | 1.177 | 1.209 | 1.168 | 1.218 | 1.185 | 1.196 | 1.193 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 4.148 | 3.028 | 0.7301 |
| deform_topology | 2.023 | 1.016 | 0.5025 |
| deform | 1.906 | 0.4902 | 0.2572 |
| topology | 0.07971 | 0.04238 | 0.5317 |
| gap_fill | 0.01505 | 0.438 | 29.11 |
| overlap_tracking | 0.02255 | 0.04584 | 2.033 |
| faults | 0.6769 | 0.8688 | 1.283 |
| shift | 0.3782 | 0.2585 | 0.6835 |
| climate_erosion_hydrology | 0.7756 | 0.6961 | 0.8976 |
| magma_transport | 0.1299 | 0.1009 | 0.7765 |

## repro-lines: seed means (seeds 804913535)

| metric | 352.4 Myr lines | 352.9 Myr lines | 353.4 Myr lines |
|---|---|---|---|
| plates | 35 | 35 | 35 |
| land fraction | 0.3892 | 0.3924 | 0.3982 |
| sea level (m) | -826.2 | -831.5 | -838 |
| Hc drift since 0 Myr | 0 | -0.002392 | -0.008526 |
| uncovered | 0.0232 | 0.02343 | 0.0226 |
| multiply covered | 0.004735 | 0.00392 | 0.00448 |
| nodes inside other plate | 0.01379 | 0.01257 | 0.01383 |
| stacked | 0.275 | 0.2761 | 0.2767 |
| anisotropic | 0.04662 | 0.04726 | 0.04814 |
| thin | 0.05789 | 0.05774 | 0.05703 |
| one-node lines | 0.293 | 0.3036 | 0.3096 |
| quad one-cell-thin | — | — | — |
| quad hole loops | — | — | — |
| air temperature (°C) | -7.643 | -7.687 | -7.674 |
| precipitation (mm) | 1021 | 1027 | 1027 |
| river fraction of land | 0.1134 | 0.1131 | 0.1132 |
| HEALPix same node | 0.5346 | 0.5358 | 0.534 |
| HEALPix p95 distance (s) | 1.174 | 1.17 | 1.172 |

## stress: seed means (seeds 1, 2)

| metric | 0 Myr lines | 0 Myr quad | 100 Myr lines | 100 Myr quad | 200 Myr lines | 200 Myr quad | 400 Myr lines | 400 Myr quad | 600 Myr lines | 600 Myr quad | 800 Myr lines | 800 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 16 | 16 | 23 | 23 | 37.5 | 30 | 51.5 | 45 | 55 | 48.5 | 66 | 55 |
| land fraction | 0.2093 | 0.2092 | 0.2156 | 0.2146 | 0.2216 | 0.1751 | 0.1522 | 0.1235 | 0.1508 | 0.08921 | 0.1685 | 0.05569 |
| sea level (m) | 0 | 0 | -1016 | -693.4 | -1576 | -896 | -2174 | -1293 | -2426 | -1519 | -2398 | -1719 |
| Hc drift since 0 Myr | 0 | 0 | -0.1252 | -0.05925 | -0.2088 | -0.0883 | -0.4292 | -0.2093 | -0.5314 | -0.3112 | -0.5158 | -0.4243 |
| uncovered | 0.03319 | 0.009175 | 0.06095 | 0.01409 | 0.08425 | 0.01812 | 0.1101 | 0.0221 | 0.1244 | 0.02378 | 0.1287 | 0.02479 |
| multiply covered | 0.05803 | 0.01059 | 0.01105 | 0.01468 | 0.00851 | 0.01956 | 0.007335 | 0.01976 | 0.00406 | 0.02018 | 0.00489 | 0.02089 |
| nodes inside other plate | 0.05264 | 0.01622 | 0.02089 | 0.00858 | 0.02764 | 0.01201 | 0.02883 | 0.008641 | 0.02687 | 0.006544 | 0.03238 | 0.004882 |
| stacked | 0 | 0 | 0.07581 | 0 | 0.1975 | 0 | 0.3249 | 0 | 0.4042 | 0 | 0.3914 | 0 |
| anisotropic | 0.00095 | 0.00029 | 0.008553 | 0.005308 | 0.02869 | 0.008293 | 0.03787 | 0.01148 | 0.05097 | 0.01459 | 0.06037 | 0.01626 |
| thin | 0.000337 | 0.000354 | 0.007236 | 0.008211 | 0.01588 | 0.0117 | 0.02658 | 0.01743 | 0.032 | 0.02072 | 0.0381 | 0.02324 |
| one-node lines | 0.01447 | — | 0.1293 | — | 0.2361 | — | 0.3168 | — | 0.3541 | — | 0.3769 | — |
| quad one-cell-thin | — | 0.000966 | — | 0.01546 | — | 0.02184 | — | 0.02675 | — | 0.03019 | — | 0.03407 |
| quad hole loops | — | 0 | — | 2 | — | 3 | — | 2 | — | 5.5 | — | 4 |
| air temperature (°C) | 0.5321 | 0.4815 | 4.736 | 2.247 | 2.253 | -2.986 | 4.766 | -6.866 | 9.247 | -4.279 | 11.7 | -2.301 |
| precipitation (mm) | 1060 | 1061 | 1108 | 1086 | 1152 | 1130 | 1215 | 1160 | 1230 | 1206 | 1238 | 1240 |
| river fraction of land | — | — | 0.1859 | 0.205 | 0.1765 | 0.276 | 0.2106 | 0.3082 | 0.2029 | 0.3239 | 0.1882 | 0.3372 |
| HEALPix same node | 0.5007 | 0.5048 | 0.4918 | 0.4993 | 0.4742 | 0.4956 | 0.4573 | 0.4957 | 0.4451 | 0.4968 | 0.4447 | 0.4962 |
| HEALPix p95 distance (s) | 1.077 | 1.119 | 1.115 | 1.114 | 1.138 | 1.107 | 1.156 | 1.101 | 1.211 | 1.106 | 1.179 | 1.11 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 2.788 | 2.016 | 0.7233 |
| deform_topology | 1.366 | 0.6298 | 0.461 |
| deform | 1.287 | 0.2903 | 0.2256 |
| topology | 0.05961 | 0.02855 | 0.4789 |
| gap_fill | 0.006124 | 0.2796 | 45.66 |
| overlap_tracking | 0.01368 | 0.03136 | 2.292 |
| faults | 0.5684 | 0.6882 | 1.211 |
| shift | 0.2324 | 0.1426 | 0.6135 |
| climate_erosion_hydrology | 0.4827 | 0.4665 | 0.9665 |
| magma_transport | 0.04886 | 0.04804 | 0.9832 |
