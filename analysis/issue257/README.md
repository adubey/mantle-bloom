# Issue #257: whole-world sums weighted by exact quad cell areas

`long/` holds the quad runs from
`paired --preset long --seeds 1,2,3,4,5 --surfaces quad --jobs 5 --out ../analysis/issue257/long`,
on this branch rebased onto main after #259 (gap filling every step) and #256 (Hm caps). The
change does not touch line worlds (line plates keep the nominal area), so the line baseline is
`../issue249-campaign/long`. `long-vs-249-lines/` is `compare` run over the quad JSONs here
plus that campaign's `seed*-lines.json`; the inputs are not duplicated. To reproduce, copy both
sets into one directory and run `surface_parity.py compare` on it.

## S2 `stability_sea_level_m`

Passes on every seed. Worst seed 1: quad 6.1 m/step against lines 10.0 m (limit 21.1 m). The
#249 campaign had quad at 58–76 m/step. That drop is #259's fix; exact areas keep it.

## Quad sea level by run (m)

Second-half std of per-step sea-level change in the last column.

| seed | run | 60 Myr | 120 | 240 | 400 | std Δ |
|---|---|---:|---:|---:|---:|---:|
| 1 | #249, nominal, older main | −768 | −1018 | −1438 | −1845 | 58.2 |
| 1 | #259, nominal, main | −766 | −968 | −1652 | −2112 | 5.5 |
| 1 | this, exact | −672 | −773 | −1320 | −2041 | 6.1 |
| 2 | #249 | −788 | −972 | −1236 | −1655 | 76.3 |
| 2 | #259 | −732 | −778 | −1266 | −1748 | 5.0 |
| 2 | this | −678 | −731 | −1327 | −1603 | 4.1 |
| 3 | #249 | −748 | −1008 | −1101 | −1191 | 72.7 |
| 3 | this | −539 | −712 | −935 | −1135 | 4.3 |
| 4 | #249 | −592 | −732 | −1175 | −1534 | 72.0 |
| 4 | #259 | −609 | −658 | −1223 | −1603 | 5.5 |
| 4 | this | −529 | −480 | −1052 | −1417 | 5.2 |
| 5 | #249 | −436 | −512 | −518 | −752 | 73.9 |
| 5 | #259 | −420 | −412 | −343 | −736 | 6.4 |
| 5 | this | −406 | −355 | −338 | −660 | 6.0 |

(#259 kept no seed 3. Its runs predate #256, so the #259 → this comparison is not fully
isolated.)

Exact areas mostly raise quad sea level, as the bias measurement in #257 predicts: the
real / nominal ocean-volume ratio falls over a run, so a budget held in nominal units loses
real water as the world ages. That moves quad away from lines on P1 `sea_level_m`: at
240 Myr, quad mean is −994 m against lines −1696 m (σ 192 m), up from −1094 m in the #249
campaign. The remaining gap is the shallower quad ocean floor (#254), which this change
uncovers rather than causes.

## Other gates

Every gate that fails here (K1, K2, K4, P1 `elevation_p05`, P1 `sea_level_m`) also failed in
the #249 campaign. C3 `multiply_covered`, which failed there, now passes (#255).
