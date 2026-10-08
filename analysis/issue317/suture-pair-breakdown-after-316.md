# Issue #317 suture Hm by plate pair

Replay commit `e11a04be19c5a19f7fa19a7c5e39e453622f1f19`, final checkpoint at 130.5 Myr.
All volumes are exact-area Hm volumes, in million km³ (M km³).

The final ledger contains 4,282 connected donor fronts across 18 single-neighbor bilateral pairs and 9 multi-neighbor candidate groups. Each direction is donor plate → neighbor plate.

For a single-neighbor pair, both donor directions are shown separately. A front whose adjacent cells identify multiple candidate plates remains grouped under its full candidate set; its volume is not divided among candidates. The saved ledger aggregates fronts by this key, so it reports front counts and grouped volumes rather than individual spatial front IDs.

## Single-neighbor bilateral pairs

| Donor → neighbor | Fronts | Donor Hm | Placed on survivors | Unplaced / delaminated |
|---|---:|---:|---:|---:|
| 3->21 | 729 | 1,513.074 | 921.095 | 591.980 |
| 22->1 | 475 | 1,055.733 | 954.342 | 101.391 |
| 21->3 | 431 | 818.661 | 817.327 | 1.334 |
| 1->22 | 327 | 786.793 | 782.978 | 3.815 |
| 2->23 | 425 | 468.155 | 415.451 | 52.705 |
| 4->1 | 217 | 411.841 | 403.309 | 8.533 |
| 1->4 | 253 | 360.567 | 345.774 | 14.793 |
| 1->23 | 123 | 347.841 | 315.069 | 32.772 |
| 23->1 | 115 | 280.728 | 266.078 | 14.650 |
| 9->1 | 227 | 251.942 | 247.153 | 4.789 |
| 20->21 | 139 | 210.965 | 207.731 | 3.234 |
| 20->3 | 100 | 150.927 | 149.512 | 1.415 |
| 1->9 | 200 | 146.655 | 146.655 | 0.000 |
| 22->2 | 87 | 144.477 | 144.477 | 0.000 |
| 21->20 | 106 | 125.532 | 125.532 | 0.000 |
| 3->20 | 113 | 106.667 | 106.363 | 0.304 |
| 2->22 | 61 | 89.756 | 89.756 | 0.000 |
| 23->2 | 23 | 38.578 | 38.578 | 0.000 |
| 2->1 | 22 | 36.455 | 36.455 | 0.000 |
| 1->2 | 26 | 32.234 | 32.234 | 0.000 |
| 9->22 | 16 | 12.543 | 12.543 | 0.000 |
| 22->9 | 14 | 8.704 | 8.704 | 0.000 |
| 23->9 | 3 | 3.504 | 3.504 | 0.000 |
| 9->23 | 4 | 3.334 | 3.334 | 0.000 |
| 21->22 | 6 | 0.349 | 0.349 | 0.000 |
| 22->21 | 7 | 0.322 | 0.322 | 0.000 |
| 21->23 | 4 | 0.187 | 0.187 | 0.000 |
| 3->4 | 1 | 0.168 | 0.168 | 0.000 |
| 23->21 | 4 | 0.132 | 0.132 | 0.000 |
| 20->22 | 1 | 0.110 | 0.110 | 0.000 |
| 21->1 | 2 | 0.094 | 0.094 | 0.000 |
| 4->3 | 3 | 0.093 | 0.093 | 0.000 |
| 9->20 | 1 | 0.049 | 0.049 | 0.000 |

## Multi-neighbor candidate groups

| Donor → candidate neighbors | Fronts | Donor Hm | Placed on survivors | Unplaced / delaminated |
|---|---:|---:|---:|---:|
| 21->3|20 | 5 | 9.215 | 8.494 | 0.721 |
| 3->20|21 | 5 | 5.336 | 4.526 | 0.810 |
| 23->1|9 | 1 | 2.937 | 1.271 | 1.666 |
| 1->9|23 | 1 | 2.931 | 2.931 | 0.000 |
| 1->9|22 | 1 | 2.321 | 2.321 | 0.000 |
| 20->3|21 | 1 | 0.848 | 0.848 | 0.000 |
| 9->1|23 | 1 | 0.434 | 0.434 | 0.000 |
| 9->1|22 | 1 | 0.054 | 0.054 | 0.000 |
| 3->15|20 | 1 | 0.030 | 0.030 | 0.000 |

## Reconciliation

| Measure | By pair/candidate groups | Recorded total |
|---|---:|---:|
| Fronts | 4,282 | 4,282 |
| Donor Hm (M km³) | 7,431.277 | 7,431.277 |
| Placed Hm (M km³) | 6,596.365 | 6,596.365 |
| Unplaced / delaminated Hm (M km³) | 834.912 | 834.912 |
