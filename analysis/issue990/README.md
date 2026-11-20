# Issue 990 · Lake–glacier boundary

[#990](https://github.com/adubey/mantle-bloom/issues/990) implements
[RFC 0015 D4](../../rfc/0015-frozen-water-in-physical-units.md#d4-lakes-that-freeze-and-lakes-under-ice-with-rfc-0014).

Principles 5, 6, and 8 govern the boundary. Wet snowfall was counted as glacier
snowmelt even though it entered the lake supply directly. The lake solve used bare
terrain under glacier ice. D4 assigns snowfall once, balances water-equivalent mass
across the cap handoff, and excludes visibly glaciated cells from lake storage.

## Closing-ice decision

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
between the two geometry arms. This is one step on one diagnostic world. The two
replay wall times are from a shared machine with other active simulations.

On the closing step the lake and land residuals were +204 and −2,690 m³ of water
and the glacier residual was −718 m³ of ice. No cap transfer occurred on that
step. The focused tests exercise both handoffs with nonzero transfers and compare
their lake, land, glacier ice, and glacier water ledgers.
