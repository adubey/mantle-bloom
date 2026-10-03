# Parity campaign summary

| set | runs | verdict |
|---|---:|---|
| issue147 | 2 | fail |
| issue147-clean | 2 | fail |
| long | 10 | fail |
| repro-equivalent | 2 | fail |
| repro-saves | 2 | warn |
| stress | 4 | fail |

## Gate status by set

| gate | issue147 | issue147-clean | long | repro-equivalent | repro-saves | stress |
|---|---|---|---|---|---|---|
| H1:quad_topology | pass | pass | pass | pass | pass | pass |
| H2:quad_neighbours | pass | pass | pass | pass | pass | pass |
| H3:quad_folded | pass | pass | pass | pass | pass | pass |
| H4:fields | pass | pass | pass | pass | pass | pass |
| H5:frames | pass | pass | pass | pass | pass | pass |
| H6:derived_caches | pass | pass | pass | pass | pass | pass |
| H7:revisions | pass | pass | pass | pass | pass | pass |
| H8:load_round_trip | pass | pass | pass | pass | pass | pass |
| H9:load_continuation | pass | pass | pass | pass | pass | pass |
| H10:quad_no_stacked_nodes | pass | pass | pass | pass | pass | pass |
| H11:field_caps | pass | pass | pass | pass | info (+1 line) | pass |
| C1:uncovered | pass | pass | pass | pass | pass | pass |
| C2:void | pass | pass | pass | pass | pass | pass |
| C3:multiply_covered | pass | pass | pass | pass | pass | warn |
| C4:nodes_inside_other_plate | pass | pass | pass | pass | pass | pass |
| K1:hc_volume_km3_drift | pass | pass | fail | fail | pass | fail |
| K2:continental_hc_volume_km3_drift | pass | pass | fail | warn | pass | fail |
| K3:quad_area_accounting | pass | pass | pass | pass | pass | pass |
| K4:hm_volume_km3_drift | fail | fail | fail | fail | pass | fail |
| M1:anisotropic_fraction | pass | pass | pass | pass | pass | warn |
| M2:row_alignment | pass | pass | pass | pass | pass | pass |
| M3:thin_fraction | pass | pass | warn | pass | pass | warn |
| M4:aspect_gt_4 | pass | pass | pass | pass | pass | pass |
| M5:skew_gt_45 | pass | pass | pass | pass | pass | pass |
| M6:fragments | pass | pass | warn | warn | warn | warn |
| S1:climate_hydrology_finite | pass | pass | pass | pass | pass | pass |
| S2:stability_air_temperature_mean_c | pass | pass | pass | warn | pass | pass |
| S2:stability_land_fraction | pass | pass | pass | pass | pass | pass |
| S2:stability_ocean_temperature_mean_c | pass | pass | pass | pass | pass | pass |
| S2:stability_precipitation_mean_mm | pass | pass | pass | pass | pass | pass |
| S2:stability_sea_level_m | pass | pass | pass | pass | pass | pass |
| X1:healpix_all | pass | pass | pass | pass | pass | pass |
| X1:healpix_antimeridian | pass | pass | pass | pass | pass | pass |
| X1:healpix_coast | pass | pass | pass | pass | pass | pass |
| X1:healpix_hole | insufficient | insufficient | warn | insufficient | pass | warn |
| X1:healpix_plate_boundary | pass | pass | pass | pass | pass | pass |
| X1:healpix_pole | pass | pass | pass | pass | pass | pass |
| R1:deform_topology_s_per_step | warn | pass | pass | pass | pass | pass |
| R2:step_total_s_per_step | pass | pass | pass | pass | pass | pass |
| P1:continental_hc_volume_km3 | insufficient | insufficient | pass | insufficient | insufficient | insufficient |
| P1:elevation_p05 | insufficient | insufficient | fail | insufficient | insufficient | insufficient |
| P1:elevation_p50 | insufficient | insufficient | pass | insufficient | insufficient | insufficient |
| P1:elevation_p95 | insufficient | insufficient | pass | insufficient | insufficient | insufficient |
| P1:hc_volume_km3 | insufficient | insufficient | pass | insufficient | insufficient | insufficient |
| P1:land_fraction | insufficient | insufficient | pass | insufficient | insufficient | insufficient |
| P1:plates | insufficient | insufficient | pass | insufficient | insufficient | insufficient |
| P1:sea_level_m | insufficient | insufficient | fail | insufficient | insufficient | insufficient |

## issue147: seed means (seeds 0)

| metric | 0 Myr lines | 0 Myr quad | 10 Myr lines | 10 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad |
|---|---|---|---|---|---|---|---|---|
| plates | 19 | 19 | 28 | 27 | 31 | 28 | 38 | 31 |
| land fraction | 0.153 | 0.153 | 0.1737 | 0.1637 | 0.1677 | 0.1556 | 0.156 | 0.1499 |
| sea level (m) | 0 | 0 | -373.9 | -383.6 | -695.5 | -570.5 | -865.9 | -613.6 |
| Hc drift since 0 Myr | 0 | 0 | -0.0102 | -0.02604 | 0.01493 | -0.0179 | 0.02745 | -0.002964 |
| uncovered | 0.01328 | 0.003935 | 0.02469 | 0.005625 | 0.02671 | 0.00624 | 0.03074 | 0.00705 |
| multiply covered | 0.01601 | 0.0038 | 0.00614 | 0.00425 | 0.001435 | 0.00476 | 0.001725 | 0.00518 |
| nodes inside other plate | 0.01771 | 0.005561 | 0.007473 | 0.002135 | 0.004076 | 0.002616 | 0.005954 | 0.002047 |
| stacked | 0 | 0 | 0.008937 | 0 | 0.07282 | 0 | 0.1344 | 0 |
| anisotropic | 0.000145 | 3.19e-05 | 0.001019 | 0.00044 | 0.004386 | 0.00112 | 0.0137 | 0.001418 |
| thin | 3.06e-05 | 2.39e-05 | 0.001226 | 0.000896 | 0.002352 | 0.001616 | 0.003608 | 0.002843 |
| one-node lines | 0.002784 | — | 0.0404 | — | 0.08951 | — | 0.152 | — |
| quad one-cell-thin | — | 0.00016 | — | 0.002623 | — | 0.003184 | — | 0.004364 |
| quad hole loops | — | 0 | — | 0 | — | 1 | — | 1 |
| air temperature (°C) | 5.78 | 5.798 | 3.999 | 4.32 | 0.3401 | -0.3207 | -2.967 | -5.172 |
| precipitation (mm) | 1096 | 1096 | 1129 | 1135 | 1126 | 1139 | 1153 | 1156 |
| river fraction of land | — | — | 0.08387 | 0.09397 | 0.1018 | 0.1212 | 0.1135 | 0.1364 |
| HEALPix same node | 0.5559 | 0.5458 | 0.5557 | 0.5472 | 0.5477 | 0.5457 | 0.5414 | 0.5465 |
| HEALPix p95 distance (s) | 1.183 | 1.211 | 1.196 | 1.213 | 1.207 | 1.221 | 1.205 | 1.231 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 10.19 | 8.38 | 0.8226 |
| deform_topology | 2.016 | 2.93 | 1.453 |
| deform | 1.755 | 1.487 | 0.8471 |
| topology | 0.1669 | 0.1572 | 0.9419 |
| gap_fill | 0.04361 | 1.165 | 26.72 |
| overlap_tracking | 0.05 | 0.1201 | 2.402 |
| faults | 0.7641 | 1.203 | 1.575 |
| shift | 0.5636 | 0.975 | 1.73 |
| climate_erosion_hydrology | 0.9973 | 1.454 | 1.458 |
| magma_transport | 5.488 | 1.471 | 0.268 |

## issue147-clean: seed means (seeds 0)

| metric | 0 Myr lines | 0 Myr quad | 10 Myr lines | 10 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad |
|---|---|---|---|---|---|---|---|---|
| plates | 19 | 19 | 28 | 27 | 31 | 28 | 38 | 31 |
| land fraction | 0.153 | 0.153 | 0.1737 | 0.1637 | 0.1677 | 0.1556 | 0.156 | 0.1499 |
| sea level (m) | 0 | 0 | -373.9 | -383.6 | -695.5 | -570.5 | -865.9 | -613.6 |
| Hc drift since 0 Myr | 0 | 0 | -0.0102 | -0.02604 | 0.01493 | -0.0179 | 0.02745 | -0.002964 |
| uncovered | 0.01328 | 0.003935 | 0.02469 | 0.005625 | 0.02671 | 0.00624 | 0.03074 | 0.00705 |
| multiply covered | 0.01601 | 0.0038 | 0.00614 | 0.00425 | 0.001435 | 0.00476 | 0.001725 | 0.00518 |
| nodes inside other plate | 0.01771 | 0.005561 | 0.007473 | 0.002135 | 0.004076 | 0.002616 | 0.005954 | 0.002047 |
| stacked | 0 | 0 | 0.008937 | 0 | 0.07282 | 0 | 0.1344 | 0 |
| anisotropic | 0.000145 | 3.19e-05 | 0.001019 | 0.00044 | 0.004386 | 0.00112 | 0.0137 | 0.001418 |
| thin | 3.06e-05 | 2.39e-05 | 0.001226 | 0.000896 | 0.002352 | 0.001616 | 0.003608 | 0.002843 |
| one-node lines | 0.002784 | — | 0.0404 | — | 0.08951 | — | 0.152 | — |
| quad one-cell-thin | — | 0.00016 | — | 0.002623 | — | 0.003184 | — | 0.004364 |
| quad hole loops | — | 0 | — | 0 | — | 1 | — | 1 |
| air temperature (°C) | 5.78 | 5.798 | 3.999 | 4.32 | 0.3401 | -0.3207 | -2.967 | -5.172 |
| precipitation (mm) | 1096 | 1096 | 1129 | 1135 | 1126 | 1139 | 1153 | 1156 |
| river fraction of land | — | — | 0.08387 | 0.09397 | 0.1018 | 0.1212 | 0.1135 | 0.1364 |
| HEALPix same node | 0.5559 | 0.5458 | 0.5557 | 0.5472 | 0.5477 | 0.5457 | 0.5414 | 0.5465 |
| HEALPix p95 distance (s) | 1.183 | 1.211 | 1.196 | 1.213 | 1.207 | 1.221 | 1.205 | 1.231 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 8.518 | 4.068 | 0.4775 |
| deform_topology | 1.583 | 1.36 | 0.859 |
| deform | 1.358 | 0.5351 | 0.3941 |
| topology | 0.1436 | 0.08609 | 0.5997 |
| gap_fill | 0.03967 | 0.671 | 16.91 |
| overlap_tracking | 0.04231 | 0.0677 | 1.6 |
| faults | 0.6358 | 0.6206 | 0.976 |
| shift | 0.4111 | 0.3133 | 0.762 |
| climate_erosion_hydrology | 0.8442 | 0.7661 | 0.9076 |
| magma_transport | 4.733 | 0.7909 | 0.1671 |

## long: seed means (seeds 1, 2, 3, 4, 5)

| metric | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 60 Myr lines | 60 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 400 Myr lines | 400 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 16.6 | 16.6 | 18.6 | 18.4 | 22 | 19.2 | 27.2 | 25.6 | 41.8 | 33.2 | 53.2 | 45.2 |
| land fraction | 0.1264 | 0.1263 | 0.1454 | 0.1425 | 0.1445 | 0.1438 | 0.144 | 0.1382 | 0.1463 | 0.1175 | 0.136 | 0.0673 |
| sea level (m) | 0 | 0 | -609.4 | -530.3 | -762 | -580.8 | -1111 | -644.2 | -1496 | -867.1 | -2006 | -1203 |
| Hc drift since 0 Myr | 0 | 0 | -0.04228 | -0.04007 | -0.08569 | -0.04401 | -0.1635 | -0.06736 | -0.2926 | -0.06403 | -0.4194 | -0.2196 |
| uncovered | 0.02541 | 0.00719 | 0.03561 | 0.008785 | 0.04627 | 0.009183 | 0.05578 | 0.01107 | 0.0703 | 0.01354 | 0.08179 | 0.01653 |
| multiply covered | 0.03456 | 0.007152 | 0.00551 | 0.006538 | 0.006903 | 0.006843 | 0.008058 | 0.008943 | 0.04713 | 0.01089 | 0.02665 | 0.01164 |
| nodes inside other plate | 0.03384 | 0.0108 | 0.01313 | 0.002955 | 0.01517 | 0.003189 | 0.01952 | 0.004782 | 0.04192 | 0.00549 | 0.03394 | 0.003918 |
| stacked | 0 | 0 | 0.01182 | 0 | 0.03673 | 0 | 0.121 | 0 | 0.2413 | 0 | 0.352 | 0 |
| anisotropic | 0.000447 | 0.000192 | 0.002808 | 0.001502 | 0.005508 | 0.00193 | 0.007902 | 0.002716 | 0.01913 | 0.003618 | 0.0424 | 0.005515 |
| thin | 0.000135 | 0.000115 | 0.003251 | 0.002055 | 0.006503 | 0.002994 | 0.008065 | 0.003532 | 0.01394 | 0.004597 | 0.02147 | 0.006301 |
| one-node lines | 0.008657 | — | 0.08113 | — | 0.1205 | — | 0.1621 | — | 0.2137 | — | 0.2855 | — |
| quad one-cell-thin | — | 0.000542 | — | 0.004492 | — | 0.005624 | — | 0.007074 | — | 0.00932 | — | 0.01213 |
| quad hole loops | — | 0 | — | 0 | — | 1 | — | 0.8 | — | 0.8 | — | 3 |
| air temperature (°C) | -3.246 | -3.281 | -1.717 | -2.189 | -2.665 | -4.488 | -1.063 | -4.867 | 1.925 | -10.34 | 3.363 | -10.46 |
| precipitation (mm) | 1144 | 1144 | 1182 | 1184 | 1174 | 1171 | 1185 | 1177 | 1203 | 1173 | 1220 | 1201 |
| river fraction of land | — | — | 0.1395 | 0.1479 | 0.1668 | 0.1583 | 0.1636 | 0.1861 | 0.18 | 0.2312 | 0.1942 | 0.2806 |
| HEALPix same node | 0.5547 | 0.5499 | 0.5526 | 0.5516 | 0.55 | 0.5514 | 0.5409 | 0.5489 | 0.5303 | 0.5492 | 0.5162 | 0.5496 |
| HEALPix p95 distance (s) | 1.18 | 1.201 | 1.2 | 1.199 | 1.219 | 1.198 | 1.224 | 1.195 | 1.233 | 1.2 | 1.238 | 1.202 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 3.678 | 2.69 | 0.7315 |
| deform_topology | 1.852 | 0.9331 | 0.5039 |
| deform | 1.748 | 0.4356 | 0.2492 |
| topology | 0.07185 | 0.03703 | 0.5154 |
| gap_fill | 0.01238 | 0.4222 | 34.1 |
| overlap_tracking | 0.01953 | 0.03819 | 1.956 |
| faults | 0.6175 | 0.7269 | 1.177 |
| shift | 0.3359 | 0.2364 | 0.7039 |
| climate_erosion_hydrology | 0.6212 | 0.6057 | 0.975 |
| magma_transport | 0.1051 | 0.1058 | 1.007 |

## repro-equivalent: seed means (seeds 804913535)

| metric | 0 Myr lines | 0 Myr quad | 30 Myr lines | 30 Myr quad | 120 Myr lines | 120 Myr quad | 240 Myr lines | 240 Myr quad | 352 Myr lines | 352 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|
| plates | 9 | 9 | 11 | 11 | 21 | 25 | 32 | 32 | 45 | 38 |
| land fraction | 0.2737 | 0.2741 | 0.3143 | 0.3077 | 0.2703 | 0.2277 | 0.1804 | 0.1462 | 0.1211 | 0.1116 |
| sea level (m) | 0 | 0 | -590.6 | -554.8 | -1379 | -1236 | -2429 | -1733 | -2760 | -1970 |
| Hc drift since 0 Myr | 0 | 0 | -0.02386 | -0.01998 | -0.1319 | -0.03946 | -0.4189 | -0.2857 | -0.5519 | -0.3815 |
| uncovered | 0.0202 | 0.00533 | 0.02459 | 0.006155 | 0.04796 | 0.01046 | 0.06039 | 0.01264 | 0.06531 | 0.01424 |
| multiply covered | 0.01241 | 0.00524 | 0.004545 | 0.00648 | 0.01961 | 0.01013 | 0.00424 | 0.0111 | 0.005185 | 0.0101 |
| nodes inside other plate | 0.01452 | 0.007697 | 0.01404 | 0.006311 | 0.0427 | 0.007051 | 0.01467 | 0.006929 | 0.02196 | 0.003202 |
| stacked | 0 | 0 | 0.003967 | 0 | 0.08112 | 0 | 0.2461 | 0 | 0.2436 | 0 |
| anisotropic | 0.000306 | 9.39e-05 | 0.002044 | 0.00087 | 0.004858 | 0.001923 | 0.01594 | 0.003873 | 0.03052 | 0.004453 |
| thin | 6.13e-05 | 0 | 0.004516 | 0.000964 | 0.007053 | 0.003052 | 0.01021 | 0.004569 | 0.01318 | 0.005703 |
| one-node lines | 0 | — | 0.06045 | — | 0.08988 | — | 0.1823 | — | 0.2408 | — |
| quad one-cell-thin | — | 0.000219 | — | 0.002425 | — | 0.00586 | — | 0.008774 | — | 0.01089 |
| quad hole loops | — | 0 | — | 0 | — | 0 | — | 1 | — | 2 |
| air temperature (°C) | 8.022 | 8.051 | 5.354 | 4.894 | 0.7672 | -8.555 | -1.416 | -3.375 | -12.73 | -9.614 |
| precipitation (mm) | 1004 | 1004 | 1031 | 1043 | 1093 | 1071 | 1199 | 1190 | 1215 | 1195 |
| river fraction of land | — | — | 0.1078 | 0.1121 | 0.1358 | 0.1755 | 0.1926 | 0.2661 | 0.2232 | 0.2288 |
| HEALPix same node | 0.5548 | 0.5569 | 0.5588 | 0.5551 | 0.5447 | 0.5511 | 0.53 | 0.5474 | 0.5327 | 0.5494 |
| HEALPix p95 distance (s) | 1.176 | 1.177 | 1.178 | 1.176 | 1.209 | 1.165 | 1.218 | 1.198 | 1.196 | 1.207 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 2.642 | 2.023 | 0.7655 |
| deform_topology | 1.263 | 0.6813 | 0.5396 |
| deform | 1.184 | 0.3115 | 0.2632 |
| topology | 0.05362 | 0.03067 | 0.5721 |
| gap_fill | 0.01061 | 0.3096 | 29.19 |
| overlap_tracking | 0.01496 | 0.02953 | 1.975 |
| faults | 0.4471 | 0.5508 | 1.232 |
| shift | 0.2158 | 0.1632 | 0.7561 |
| climate_erosion_hydrology | 0.5146 | 0.4994 | 0.9705 |
| magma_transport | 0.09105 | 0.06371 | 0.6997 |

## repro-saves: seed means (seeds 804913535)

| metric | 352.4 Myr lines | 352.4 Myr quad | 352.9 Myr lines | 352.9 Myr quad | 353.4 Myr lines | 353.4 Myr quad |
|---|---|---|---|---|---|---|
| plates | 35 | 35 | 35 | 35 | 35 | 36 |
| land fraction | 0.3892 | 0.3313 | 0.3924 | 0.3288 | 0.3982 | 0.3291 |
| sea level (m) | -826.2 | -826.2 | -831.5 | -810.5 | -838 | -815.9 |
| Hc drift since 0 Myr | 0 | 0 | -0.002392 | -0.000749 | -0.008526 | -0.001449 |
| uncovered | 0.0232 | 0.00829 | 0.02343 | 0.00591 | 0.0226 | 0.00612 |
| multiply covered | 0.004735 | 0.00713 | 0.00392 | 0.004805 | 0.00448 | 0.00465 |
| nodes inside other plate | 0.01379 | 0.01165 | 0.01257 | 0.002717 | 0.01383 | 0.002417 |
| stacked | 0.275 | 0 | 0.2761 | 0 | 0.2767 | 0 |
| anisotropic | 0.04662 | 0.001235 | 0.04726 | 0.000641 | 0.04814 | 0.000672 |
| thin | 0.05789 | 0.001389 | 0.05774 | 0.000841 | 0.05703 | 0.00088 |
| one-node lines | 0.293 | — | 0.3036 | — | 0.3096 | — |
| quad one-cell-thin | — | 0.002971 | — | 0.002987 | — | 0.003228 |
| quad hole loops | — | 38 | — | 3 | — | 2 |
| air temperature (°C) | -7.643 | -7.643 | -7.687 | -8.07 | -7.674 | -8.052 |
| precipitation (mm) | 1021 | 1021 | 1027 | 1024 | 1027 | 1025 |
| river fraction of land | 0.1134 | — | 0.1131 | 0.1243 | 0.1132 | 0.1232 |
| HEALPix same node | 0.5346 | 0.5544 | 0.5358 | 0.5557 | 0.534 | 0.5567 |
| HEALPix p95 distance (s) | 1.174 | 1.179 | 1.17 | 1.172 | 1.172 | 1.174 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 7.958 | 5.824 | 0.7318 |
| deform_topology | 2.953 | 1.737 | 0.5884 |
| deform | 2.675 | 0.7376 | 0.2757 |
| topology | 0.1853 | 0.09648 | 0.5208 |
| gap_fill | 0.03521 | 0.81 | 23.01 |
| overlap_tracking | 0.0569 | 0.09329 | 1.64 |
| faults | 0.7586 | 0.7029 | 0.9265 |
| shift | 0.7078 | 0.4167 | 0.5888 |
| climate_erosion_hydrology | 1.123 | 0.949 | 0.8452 |
| magma_transport | 2.003 | 1.824 | 0.9106 |

## stress: seed means (seeds 1, 2)

| metric | 0 Myr lines | 0 Myr quad | 100 Myr lines | 100 Myr quad | 200 Myr lines | 200 Myr quad | 400 Myr lines | 400 Myr quad | 600 Myr lines | 600 Myr quad | 800 Myr lines | 800 Myr quad |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| plates | 16 | 16 | 23 | 22.5 | 37.5 | 31.5 | 51.5 | 45.5 | 55 | 49.5 | 66 | 49.5 |
| land fraction | 0.2093 | 0.2092 | 0.2156 | 0.2161 | 0.2216 | 0.1821 | 0.1522 | 0.1083 | 0.1508 | 0.07838 | 0.1685 | 0.04633 |
| sea level (m) | 0 | 0 | -1016 | -710.9 | -1576 | -919.8 | -2174 | -1461 | -2426 | -1571 | -2398 | -1728 |
| Hc drift since 0 Myr | 0 | 0 | -0.1252 | -0.05873 | -0.2088 | -0.09644 | -0.4292 | -0.2535 | -0.5314 | -0.3534 | -0.5158 | -0.4413 |
| uncovered | 0.03319 | 0.009175 | 0.06095 | 0.01355 | 0.08425 | 0.01921 | 0.1101 | 0.02351 | 0.1244 | 0.02426 | 0.1287 | 0.02436 |
| multiply covered | 0.05803 | 0.01059 | 0.01105 | 0.01162 | 0.00851 | 0.0144 | 0.007335 | 0.01603 | 0.00406 | 0.01688 | 0.00489 | 0.0148 |
| nodes inside other plate | 0.05264 | 0.01622 | 0.02089 | 0.007105 | 0.02764 | 0.007441 | 0.02883 | 0.00395 | 0.02687 | 0.004881 | 0.03238 | 0.002163 |
| stacked | 0 | 0 | 0.07581 | 0 | 0.1975 | 0 | 0.3249 | 0 | 0.4042 | 0 | 0.3914 | 0 |
| anisotropic | 0.00095 | 0.00029 | 0.008553 | 0.00472 | 0.02869 | 0.005896 | 0.03787 | 0.009915 | 0.05097 | 0.01155 | 0.06037 | 0.01086 |
| thin | 0.000337 | 0.000354 | 0.007236 | 0.006696 | 0.01588 | 0.008843 | 0.02658 | 0.01456 | 0.032 | 0.01453 | 0.0381 | 0.01381 |
| one-node lines | 0.01447 | — | 0.1293 | — | 0.2361 | — | 0.3168 | — | 0.3541 | — | 0.3769 | — |
| quad one-cell-thin | — | 0.000966 | — | 0.009093 | — | 0.01486 | — | 0.02065 | — | 0.02132 | — | 0.02388 |
| quad hole loops | — | 0 | — | 0 | — | 0.5 | — | 2.5 | — | 1.5 | — | 2 |
| air temperature (°C) | 0.5321 | 0.4815 | 4.736 | 1.737 | 2.253 | -3.791 | 4.766 | -6.741 | 9.247 | -0.0082 | 11.7 | -1.142 |
| precipitation (mm) | 1060 | 1061 | 1108 | 1088 | 1152 | 1118 | 1215 | 1173 | 1230 | 1210 | 1238 | 1239 |
| river fraction of land | — | — | 0.1859 | 0.2133 | 0.1765 | 0.2261 | 0.2106 | 0.366 | 0.2029 | 0.3182 | 0.1882 | 0.3103 |
| HEALPix same node | 0.5007 | 0.5048 | 0.4918 | 0.5 | 0.4742 | 0.4943 | 0.4573 | 0.4989 | 0.4451 | 0.497 | 0.4447 | 0.4965 |
| HEALPix p95 distance (s) | 1.077 | 1.119 | 1.115 | 1.118 | 1.138 | 1.113 | 1.156 | 1.11 | 1.211 | 1.1 | 1.179 | 1.119 |

| s/step (seed mean) | lines | quad | quad / lines |
|---|---:|---:|---:|
| step_total | 2.489 | 1.683 | 0.6764 |
| deform_topology | 1.229 | 0.5557 | 0.452 |
| deform | 1.159 | 0.236 | 0.2036 |
| topology | 0.05273 | 0.02381 | 0.4515 |
| gap_fill | 0.005695 | 0.2694 | 47.3 |
| overlap_tracking | 0.01224 | 0.02652 | 2.167 |
| faults | 0.5072 | 0.5457 | 1.076 |
| shift | 0.2089 | 0.1205 | 0.5766 |
| climate_erosion_hydrology | 0.4209 | 0.3883 | 0.9227 |
| magma_transport | 0.04295 | 0.0385 | 0.8964 |
