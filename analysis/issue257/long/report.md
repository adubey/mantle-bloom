# Line vs quad parity report

**Verdict: PASS** — configs long; seeds 1, 2, 3, 4, 5.

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
| H4:fields | pass | violations=0, audits=41, surface=quad, seed=1 |
| H5:frames | pass | violations=0, audits=41, surface=quad, seed=1 |
| H6:derived_caches | pass | violations=0, audits=41, surface=quad, seed=1 |
| H7:revisions | pass | violations=0, audits=41, surface=quad, seed=1 |
| H8:load_round_trip | pass | checks=6, surface=quad, seed=1 |
| H9:load_continuation | pass | surface=quad, seed=1 |
| H10:quad_no_stacked_nodes | pass | surface=quad, seed=1 |
| H11:field_caps | pass | violations=0, surface=quad, seed=1 |

## Coverage and overlap

| gate | status | worst case |
|---|---|---|
| C3:multiply_covered | pass | age_myr=0, quad=0.006475, warn_above=0.015, fail_above=0.02, surface=quad, seed=1 |

## Climate and hydrology

| gate | status | worst case |
|---|---|---|
| S1:climate_hydrology_finite | pass | surface=quad, seed=1 |

## Checkpoints, seed 1

| seed 1 | 0 Myr quad | 30 Myr quad | 60 Myr quad | 120 Myr quad | 240 Myr quad | 400 Myr quad |
|---|---|---|---|---|---|---|
| plates | 14 | 15 | 16 | 20 | 30 | 43 |
| nodes | 31493 | 31600 | 32011 | 32464 | 33055 | 32792 |
| land fraction | 0.2793 | 0.3118 | 0.3082 | 0.284 | 0.2028 | 0.09436 |
| sea level (m) | 0 | -631.7 | -671.9 | -773.2 | -1320 | -2041 |
| Hc volume (km³) | 1.15e+10 | 1.11e+10 | 1.08e+10 | 1.05e+10 | 8.65e+09 | 6.16e+09 |
| continental Hc (km³) | 7.65e+09 | 7.67e+09 | 7.6e+09 | 7.19e+09 | 4.68e+09 | 1.92e+09 |
| elevation p05/p50/p95 (m) | -4768/-4353/2984 | -5476/-4986/3186 | -5654/-4983/3315 | -5790/-4752/4325 | -5864/-4959/3603 | -5898/-5358/423.9 |
| uncovered | 0.005985 | 0.007185 | 0.007435 | 0.00896 | 0.01152 | 0.01488 |
| multiply covered | 0.006475 | 0.006545 | 0.007545 | 0.01131 | 0.01301 | 0.01253 |
| void | 0 | 0 | 0 | 0 | 0 | 0 |
| nodes inside other plate | 0.009653 | 0.003165 | 0.003405 | 0.008409 | 0.008259 | 0.002714 |
| stacked | 0 | 0 | 0 | 0 | 0 | 0 |
| anisotropic | 0.000159 | 0.001709 | 0.001906 | 0.002464 | 0.003751 | 0.006343 |
| row alignment | 0.08457 | 0.02076 | -0.122 | -0.2557 | 0.1135 | -0.05098 |
| thin | 0.000191 | 0.002658 | 0.003811 | 0.004035 | 0.006262 | 0.008325 |
| one-node lines | — | — | — | — | — | — |
| quad one-cell-thin | 0.000381 | 0.006171 | 0.006935 | 0.00881 | 0.01228 | 0.01534 |
| air temperature mean (°C) | 4.481 | 6.926 | 7.701 | 5.861 | 1.827 | 9.582 |
| precipitation mean (mm) | 987 | 1010 | 1000 | 1048 | 1134 | 1233 |
| river fraction of land | — | 0.114 | 0.1176 | 0.1416 | 0.1769 | 0.1471 |
| HEALPix same node | 0.5515 | 0.5502 | 0.5497 | 0.5493 | 0.5455 | 0.5449 |

## Checkpoints, seed 2

| seed 2 | 0 Myr quad | 30 Myr quad | 60 Myr quad | 120 Myr quad | 240 Myr quad | 400 Myr quad |
|---|---|---|---|---|---|---|
| plates | 18 | 21 | 20 | 26 | 36 | 43 |
| nodes | 31334 | 31544 | 31870 | 32528 | 32527 | 32529 |
| land fraction | 0.1403 | 0.1469 | 0.1438 | 0.1336 | 0.0588 | 0.01557 |
| sea level (m) | 0 | -607.4 | -677.7 | -731.2 | -1327 | -1603 |
| Hc volume (km³) | 8.91e+09 | 8.39e+09 | 8.06e+09 | 7.54e+09 | 5.73e+09 | 4.72e+09 |
| continental Hc (km³) | 4.28e+09 | 4.36e+09 | 4.16e+09 | 3.52e+09 | 1.51e+09 | 1.83e+08 |
| elevation p05/p50/p95 (m) | -5174/-4437/1652 | -5532/-5086/1813 | -5757/-5122/1695 | -5870/-5107/1566 | -5929/-5391/-841 | -5937/-5471/-4337 |
| uncovered | 0.00755 | 0.00868 | 0.00884 | 0.01103 | 0.01507 | 0.01637 |
| multiply covered | 0.007265 | 0.00897 | 0.00921 | 0.01044 | 0.01249 | 0.01351 |
| void | 0 | 0 | 0 | 0 | 0 | 0 |
| nodes inside other plate | 0.01098 | 0.004216 | 0.005993 | 0.004827 | 0.003351 | 0.002613 |
| stacked | 0 | 0 | 0 | 0 | 0 | 0 |
| anisotropic | 0.000287 | 0.002853 | 0.002322 | 0.003443 | 0.005257 | 0.006886 |
| row alignment | -0.4423 | 0.2016 | 0.05847 | 0.1916 | 0.06302 | -0.07137 |
| thin | 9.57e-05 | 0.003138 | 0.003263 | 0.004611 | 0.007747 | 0.01036 |
| one-node lines | — | — | — | — | — | — |
| quad one-cell-thin | 0.000734 | 0.007704 | 0.006934 | 0.01024 | 0.01534 | 0.01977 |
| air temperature mean (°C) | -3.453 | -0.1418 | 1.092 | 1.992 | 12.78 | 17.66 |
| precipitation mean (mm) | 1132 | 1171 | 1184 | 1189 | 1236 | 1274 |
| river fraction of land | — | 0.1555 | 0.1686 | 0.1857 | 0.1125 | 0.09483 |
| HEALPix same node | 0.5486 | 0.5518 | 0.5522 | 0.5471 | 0.5492 | 0.5483 |

## Checkpoints, seed 3

| seed 3 | 0 Myr quad | 30 Myr quad | 60 Myr quad | 120 Myr quad | 240 Myr quad | 400 Myr quad |
|---|---|---|---|---|---|---|
| plates | 18 | 21 | 24 | 27 | 30 | 37 |
| nodes | 31213 | 31537 | 32160 | 32736 | 33062 | 33031 |
| land fraction | 0.05551 | 0.07413 | 0.07705 | 0.0627 | 0.03523 | 0.01801 |
| sea level (m) | 0 | -483 | -538.9 | -712.1 | -935.3 | -1135 |
| Hc volume (km³) | 7.2e+09 | 6.9e+09 | 6.77e+09 | 5.97e+09 | 5.22e+09 | 4.95e+09 |
| continental Hc (km³) | 2.59e+09 | 2.79e+09 | 2.67e+09 | 1.79e+09 | 9.38e+08 | 5e+08 |
| elevation p05/p50/p95 (m) | -5222/-4493/41.41 | -5531/-5118/40.26 | -5733/-5162/201.1 | -5880/-5222/-272.9 | -5942/-5378/-1977 | -5945/-5449/-4151 |
| uncovered | 0.00765 | 0.008495 | 0.00947 | 0.01107 | 0.01247 | 0.01401 |
| multiply covered | 0.00755 | 0.01032 | 0.009135 | 0.01002 | 0.01052 | 0.01229 |
| void | 0 | 0 | 0 | 0 | 0 | 0 |
| nodes inside other plate | 0.01157 | 0.005676 | 0.00482 | 0.004154 | 0.00251 | 0.004238 |
| stacked | 0 | 0 | 0 | 0 | 0 | 0 |
| anisotropic | 0.000192 | 0.002029 | 0.002425 | 0.004002 | 0.004355 | 0.005086 |
| row alignment | 0.3694 | 0.04715 | -0.07017 | -0.001982 | -0.03216 | 0.1924 |
| thin | 0.000192 | 0.003044 | 0.004322 | 0.004552 | 0.005596 | 0.006963 |
| one-node lines | — | — | — | — | — | — |
| quad one-cell-thin | 0.000609 | 0.007578 | 0.008364 | 0.009958 | 0.01155 | 0.01505 |
| air temperature mean (°C) | -9.163 | -3.634 | -9.615 | -7.176 | 14.42 | -8.753 |
| precipitation mean (mm) | 1201 | 1235 | 1241 | 1257 | 1260 | 1256 |
| river fraction of land | — | 0.1604 | 0.2234 | 0.1554 | 0.1661 | 0.3857 |
| HEALPix same node | 0.556 | 0.5538 | 0.5556 | 0.5496 | 0.5496 | 0.5504 |

## Checkpoints, seed 4

| seed 4 | 0 Myr quad | 30 Myr quad | 60 Myr quad | 120 Myr quad | 240 Myr quad | 400 Myr quad |
|---|---|---|---|---|---|---|
| plates | 17 | 19 | 20 | 27 | 38 | 45 |
| nodes | 31134 | 31357 | 31721 | 32578 | 33014 | 32515 |
| land fraction | 0.1416 | 0.1585 | 0.159 | 0.1596 | 0.09205 | 0.03317 |
| sea level (m) | 0 | -536.8 | -528.9 | -479.6 | -1052 | -1417 |
| Hc volume (km³) | 8.49e+09 | 8.12e+09 | 7.97e+09 | 7.77e+09 | 6.56e+09 | 5.15e+09 |
| continental Hc (km³) | 4.13e+09 | 4.24e+09 | 4.31e+09 | 4.28e+09 | 2.42e+09 | 6.65e+08 |
| elevation p05/p50/p95 (m) | -5222/-4427/962.3 | -5506/-5055/1058 | -5719/-5063/1063 | -5824/-5011/1085 | -5913/-5273/887.1 | -5930/-5384/-3053 |
| uncovered | 0.00778 | 0.008025 | 0.00896 | 0.01094 | 0.01418 | 0.0155 |
| multiply covered | 0.00729 | 0.00819 | 0.00877 | 0.01012 | 0.01341 | 0.01263 |
| void | 0 | 0 | 0 | 0 | 0 | 0 |
| nodes inside other plate | 0.01092 | 0.002169 | 0.002806 | 0.00396 | 0.0063 | 0.00326 |
| stacked | 0 | 0 | 0 | 0 | 0 | 0 |
| anisotropic | 0.000225 | 0.001531 | 0.002427 | 0.003469 | 0.005028 | 0.006797 |
| row alignment | 0.03748 | -0.005632 | -0.06129 | -0.1176 | -0.0153 | 0.1198 |
| thin | 3.21e-05 | 0.002424 | 0.004193 | 0.004205 | 0.006088 | 0.009473 |
| one-node lines | — | — | — | — | — | — |
| quad one-cell-thin | 0.000418 | 0.00641 | 0.008449 | 0.01038 | 0.01372 | 0.01805 |
| air temperature mean (°C) | -5.829 | -7.355 | -8.058 | -4.634 | 0.7354 | 6.405 |
| precipitation mean (mm) | 1147 | 1255 | 1184 | 1180 | 1231 | 1244 |
| river fraction of land | — | 0.1556 | 0.1546 | 0.1608 | 0.2308 | 0.243 |
| HEALPix same node | 0.5487 | 0.5488 | 0.5489 | 0.5494 | 0.547 | 0.5498 |

## Checkpoints, seed 5

| seed 5 | 0 Myr quad | 30 Myr quad | 60 Myr quad | 120 Myr quad | 240 Myr quad | 400 Myr quad |
|---|---|---|---|---|---|---|
| plates | 16 | 17 | 18 | 24 | 36 | 35 |
| nodes | 31550 | 31544 | 31555 | 32091 | 33052 | 33011 |
| land fraction | 0.01469 | 0.01988 | 0.02556 | 0.03122 | 0.02867 | 0.01037 |
| sea level (m) | 0 | -432.5 | -405.5 | -354.8 | -337.9 | -660.4 |
| Hc volume (km³) | 6.19e+09 | 6.03e+09 | 6.04e+09 | 5.9e+09 | 5.5e+09 | 4.97e+09 |
| continental Hc (km³) | 1.86e+09 | 1.95e+09 | 2.02e+09 | 1.97e+09 | 1.32e+09 | 5.5e+08 |
| elevation p05/p50/p95 (m) | -5222/-4570/-4198 | -5511/-5187/-3527 | -5706/-5207/-2836 | -5844/-5221/-1872 | -5897/-5223/-1824 | -5941/-5432/-3757 |
| uncovered | 0.006985 | 0.00752 | 0.00841 | 0.01027 | 0.01406 | 0.01383 |
| multiply covered | 0.00718 | 0.00821 | 0.009195 | 0.0104 | 0.01155 | 0.01201 |
| void | 0 | 0 | 0 | 0 | 0 | 0 |
| nodes inside other plate | 0.0109 | 0.004216 | 0.004849 | 0.00402 | 0.005113 | 0.003241 |
| stacked | 0 | 0 | 0 | 0 | 0 | 0 |
| anisotropic | 9.51e-05 | 0.002061 | 0.001997 | 0.004768 | 0.005264 | 0.005089 |
| row alignment | -0.3597 | 0.1791 | -0.04844 | 0.16 | -0.07338 | 0.001211 |
| thin | 6.34e-05 | 0.003075 | 0.00469 | 0.006887 | 0.005355 | 0.007997 |
| one-node lines | — | — | — | — | — | — |
| quad one-cell-thin | 0.000571 | 0.006689 | 0.007954 | 0.01134 | 0.01349 | 0.01387 |
| air temperature mean (°C) | -2.44 | -6.007 | -6.969 | -9.196 | 4.799 | -6.939 |
| precipitation mean (mm) | 1254 | 1256 | 1262 | 1260 | 1262 | 1265 |
| river fraction of land | — | 0.1289 | 0.1342 | 0.2736 | 0.1991 | 0.2607 |
| HEALPix same node | 0.5447 | 0.5491 | 0.5494 | 0.5524 | 0.5515 | 0.5459 |

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

