# Issue 1005 · D4 review follow-up

[Issue #1005](https://github.com/adubey/mantle-bloom/issues/1005) follows the
[final D4 reviews](https://github.com/adubey/mantle-bloom/pull/992#issuecomment-7065743558).
Principles 5, 6 and 8 require consistent units, an outlet for all displaced lake
water, and snowfall that responds to the trial glacier surface. The regressions
check dry-share lapse feedback, preservation of already-partitioned wet snow,
cap eligibility with nonzero precipitation, overflow delivery from a merged
ice-only closed root, and lake phase-budget closure. The D4 probe also checks
each interval against its own opening mesh and ice mask (Principle 9).

## Replay plan

Replay one 250 kyr step from seed 997271774 at 87.4 Myr with ice ages off and
the glacier march on, matching the [D4 diagnosis](../issue990/README.md).
Use `mantle-bloom-seed997271774-87400000y.bloom` with SHA-256
`07971ee25cdbc364cd887ea8f4fffc69ceed49b1933d09ea38e3864ea95328b8`.
Use the corrected diagnostic in both backends on `origin/main` and the fix
branch. For each backend compare opening and closing glacier geometry, with
rock-bed subcell relief in both arms. The 250 kyr interval ends at 87.65 Myr.

The primary metric is the share of global liquid lake volume under newly
visible ice. D4 calls for closing geometry if the opening arm exceeds 1%.
The prior reviewed opening result was 1.3416%. Restoring lapse feedback may
increase advance and the exposed overlap, at a scale of tenths of a percentage
point or more. The closing arm should hold zero liquid under visible ice in
either backend. Compare total and regional volumes and the lake, land and ice
budget residuals. These four controlled one-step replays do not estimate
variability across seeds or intervals.

## Regression evidence

- A seasonal bed-frozen lake with a fully wet opening cell, a 0 °C annual mean,
  and 1 °C amplitude receives 650 mm/yr of precipitation. Only 325 mm/yr
  falls as snow on the lake. The previous handoff transferred the cap even
  though that snowfall could not sustain glacier ice. The corrected handoff
  keeps the cap in the lake budget at 650 mm/yr and transfers it at
  1,000 mm/yr, where the frozen share can sustain it. Lake precipitation plus
  glacier input accounts for all precipitation in both cases.
- A merged, ice-only closed root previously booked overflow but exported
  none of it. The two-leaf regression now routes all 26 × 10⁹ m³ of opening
  storage and external supply through one leaf sink and closes the lake budget.
- For `q = 6 m²/yr`, `loss = 10⁻¹³ m/yr` and a tenth of a year, adding
  `10⁻¹⁴ m/yr` of snow previously made the ice about 4.2 mm thinner. The
  corrected result stays within `10⁻¹³ m` of the Stefan trajectory and the
  freezing, melting and retained snowfall terms sum to the ice storage change.
- In the two-step lag fixture the correct newly covered volumes are 3 and
  5 m³. The old probe reported 3 and 8 m³ from the cached state and 5 and
  10 m³ from the compact state. Both save paths now give 3 and 5 m³.

## Baseline result

The opening-geometry arm on `353a37470fb7f821687687b6107c9cf3ab3f18f8`
reproduced the prior diagnosis with the corrected probe. The comparison uses
the probe at `32701054` and the baseline worktree's Python environment.

```bash
<baseline-python> bin/debug/measure_d4_lake_lag.py <save> \
    --backend <baseline-backend> --opening-ice --myr 0.25 --step-years 250000 \
    --out <capture>
```

| Backend and geometry | Total liquid (km³) | Under closing ice (km³) | Glacier margin (km³) | Elsewhere (km³) | Newly covered share |
| --- | ---: | ---: | ---: | ---: | ---: |
| Main, opening geometry | 768,728.96 | 10,312.99 | 26,742.40 | 731,673.57 | 1.3416% |

All liquid under closing ice in this arm was on newly covered cells. The lake and
land residuals were −511 and −738 m³ of water, and the glacier residual was
+1,516 m³ of ice. Cap transfer was zero. The step took 50.73 s and peak resident
memory was 0.84 GiB with other simulations running. This arm verifies the
baseline above the 1% threshold. The remaining arms will measure the change
from restoring snowfall feedback.
