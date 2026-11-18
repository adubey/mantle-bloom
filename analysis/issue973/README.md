# Issue 973: RFC 0015 D3 implementation

Implements the accepted D3 direction in [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md#d3-glacial-erosion): integrate Cook et al. (2020)'s physical basal-sliding abrasion law through the D2 ice march, update a working bed between intervals, and transfer accumulated rock once through the existing terrain and material ledgers.

**Principles:** #5 (consistent units) and #8 (each process follows its own physics). The prior closing-ice-times-world-step estimate was then clipped with the shared neighbor-drop limit, creating a 66–73% gap between 0.25 Myr and 1 Myr steps on G3's diagnosis saves. The invariant is that each march interval contributes the physical law times interval duration, the next interval receives the bed lowered by the previous interval's abrasion, no held interval omits continuing abrasion, and the final material transfer remains within crustal headroom and closes the existing sediment ledger.

## Pre-registered replay

Before replaying, expect the integrated abrasion total at 0.25 Myr and 1 Myr to differ materially less than G3's prior 66–73% gap; a gap under 20% is the working target, while a remaining gap above 50% means D3 has not removed the dominant timestep sensitivity. Expect applied/raw to rise above G3's 0.0045–0.210 range where the removed material has available crustal headroom; any remaining loss should be explainable by headroom, mobile cover and craton resistance. Runtime is a measured cost rather than a pass threshold; compare each D3 arm with the D2 baseline under the same save, horizon and step schedule, and report slowdown alongside the interval count.

Controls: Cook et al. (2020), coefficient 1e-4 m/yr, exponent 0.65, glacier_erosion_multiplier = 1, ice ages disabled, same 10 Myr horizon and deterministic save for each paired run. Use each save at 0.25 Myr × 40 steps and 1 Myr × 10 steps. D2 baseline runs use origin/main; D3 runs use the issue branch. Do not tune coefficients after seeing these results.

| Save | Age | SHA-256 |
|---|---:|---|
| mantle-bloom-seed649658380-6400000y.bloom | 6.4 Myr | 06167cb4e45fac7141e98df2cf091999ffc40b36be4c0da3bd94703b6b31131f |
| mantle-bloom-seed997271774-87400000y.bloom | 87.4 Myr | 07971ee25cdbc364cd887ea8f4fffc69ceed49b1933d09ea38e3864ea95328b8 |

### Runs

Pending preregistered D2 and D3 replays.
