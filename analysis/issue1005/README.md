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
