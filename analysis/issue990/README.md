# Issue 990 · Lake–glacier boundary

[#990](https://github.com/adubey/mantle-bloom/issues/990) implements
[RFC 0015 D4](../../rfc/0015-frozen-water-in-physical-units.md#d4-lakes-that-freeze-and-lakes-under-ice-with-rfc-0014).

Principles 5, 6, and 8 govern the boundary. Wet snowfall was counted as glacier
snowmelt even though it entered the lake supply directly. The lake solve used bare
terrain under glacier ice. D4 assigns snowfall once, balances water-equivalent mass
across the cap handoff, and excludes visibly glaciated cells from lake storage.

## Closing-ice decision

The measurements in this section are the original D4 diagnosis, before the PR review fixes.
The controlled reruns and updated decision measure are recorded in the review follow-up below.

[G1](../issue766/README.md) counted liquid under all closing glacier ice. For seed
997271774 its median was 20%, an upper bound on the one-step lag. To isolate the lag,
the probe runs D4 with opening glacier geometry and measures the ratio of liquid
lake volume on cells newly reaching 10 m of glacier depth to total liquid lake
volume. Volumes use accounting cell area. The decision rule is a maximum over 1%
on any diagnosis step, so the probe stops on the first crossing.

The probe started at 87.4 Myr on seed 997271774 with ice ages off and the production
glacier march enabled. The matching D3 diagnosis save is
`mantle-bloom-seed997271774-87400000y.bloom` with SHA-256
`07971ee25cdbc364cd887ea8f4fffc69ceed49b1933d09ea38e3864ea95328b8`.
The opening-geometry measurement ran at `57652f9574624ce10fa90219bca4715759ffbb01`.

| Geometry | Step | Liquid lake volume (km³) | Newly covered liquid (km³) | Share |
| --- | ---: | ---: | ---: | ---: |
| Opening glacier | 87.4–87.65 Myr | 770,477.39 | 13,558.60 | 1.7598% |
| Closing glacier | 87.4–87.65 Myr | 1,589,473.64 | 0 | 0% |

This triggered closing glacier geometry after the ice solve. The precipitation
partition still uses opening wet area. The lake hierarchy is built once, on the
closing surface, and the opening hypsometry is reused when glacier depth is unchanged.

Run the lag arm from the current branch with:

```bash
backend/.venv/bin/python bin/debug/measure_d4_lake_lag.py <save> \
    --opening-ice --myr 10 --step-years 250000 --out <json>
```

At the measured commit, opening geometry was the default and the flag was absent.
The confirmation at `e1a1f4ba51f681db7b91bc258120d3dae8fbdb45` used the same save and
one 250 kyr step with the production closing geometry. It found zero liquid lake
water under any visible glacier. Total liquid storage rose by 818,996.26 km³
between the two geometry arms. This is the whole-system response to changing
hydraulic geometry, not the volume intercepted under newly covered cells. The
original report did not attribute this change to the altered spill hierarchy,
ice-dammed basins, or glacier-inflated margin relief, so it is not evidence that
any one mechanism caused the increase. This is one step on one diagnostic world.
The two replay wall times are from a shared machine with other active simulations.

On the closing step the lake and land residuals were +204 and −2,690 m³ of water
and the glacier residual was −718 m³ of ice. No cap transfer occurred on that
step. The focused tests exercise both handoffs with nonzero transfers and compare
their lake, land, glacier ice, and glacier water ledgers.

## Review follow-up preregistration

Before rerunning the diagnostic after PR review, compare three one-step arms from
the same 87.4 Myr save, seed 997271774 and 250 kyr schedule: opening glacier
geometry with rock-bed relief, closing geometry with the old glacier-surface
relief, and closing geometry with rock-bed relief. Report liquid volume under
closing visible ice, on one-cell glacier margins, and elsewhere, plus total
storage and ledger residuals. The primary acceptance measure stays the
preregistered maximum share of liquid under newly covered ice, with a 1%
threshold. The relief counterfactual tests whether bed-based sub-cell relief
removes margin inflation; the three-arm totals will show whether that mechanism
explains the earlier +818,996 km³ whole-world change. These are one-step
diagnostic comparisons, not validation across seeds or step sizes.

The reruns used the same save and 250 kyr step on the reviewed branch. The
tectonic mesh changed from 132,407 to 132,248 nodes during the step, so the
opening glacier mask was mapped to closing nodes by nearest opening-state cell
centre before measuring newly covered cells. The opening arm crossed the 1%
threshold on its first step. Values below are km³ of liquid water; residuals
are in m³ in their respective ledger units.

| Geometry | Total lake water | Under closing ice | On glacier margin | Elsewhere | Newly covered share | Lake residual | Land residual | Ice residual |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Opening ice, rock-bed relief | 768,728.96 | 10,312.99 | 26,742.40 | 731,673.57 | 1.3416% | −445 | −74 | +1,516 |
| Closing ice, glacier-surface relief | 1,526,613.27 | 0 | 125,708.85 | 1,400,904.42 | 0% | −2 | −1,415 | +194 |
| Closing ice, rock-bed relief | 1,525,681.34 | 0 | 124,819.54 | 1,400,861.80 | 0% | −953 | −1,470 | +1,516 |

Replacing glacier-surface sub-cell relief with rock-bed relief lowered total
storage by 931.93 km³ (0.06% of the rock-bed total), far below the 756,952 km³
difference between the opening- and closing-geometry arms. Thus margin relief
inflation does not explain the large whole-system response. This comparison
does not separate the remaining effects of spill hierarchy, ice-dammed basins,
and closing-surface routing. The opening-geometry result also confirms that
liquid under newly covered cells is above the preregistered 1% threshold on this
diagnostic step; the production closing-geometry result removes it.

[The #1005 follow-up](../issue1005/README.md) repeated both rock-bed geometry
arms after restoring lapse-driven snowfall. The opening overlap reached
1.7628%, and closing geometry kept liquid out from under visible ice.
