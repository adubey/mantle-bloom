# Issue #368: glacier mass balance in physical units

Glacier mass balance runs in model units (decision B in #353). Snowfall builds ice at
`GLACIER_ACCUMULATION_RATE` = 0.02 m per Myr per mm/yr, while the physical exchange is
1000/917 mm of ice per mm of water, or 1,090.5 m per Myr per mm/yr: 54,526× more.
Geothermal basal melt and `WATER_TO_ICE` (frozen river and lake water) are both defined
through the same exchange. Surface melt, `GLACIER_MELT_RATE_M_PER_MYR` = 400 m/Myr per 20 °C
of year-averaged positive degrees, is model-scale too. The abrasion coefficient (430) was
calibrated against the model-unit sliding speed (#353).

This analysis asks whether converting to physical units would leave glacier extent roughly
where it is (the balance points and depth guard don't change) or change ice volume,
meltwater, erosion and depth-guard shedding enough to need an RFC.

Frozen river and lake ice are out of scope here: RFC 0014 (PR #523) covers them, and PR #504
fixed lake freezing.

## Conversion assumptions

| Quantity | Model now | Physical | Factor |
|---|---:|---:|---:|
| Accumulation (m ice / Myr per mm/yr) | 0.02 | 1,090.5 (ρ_w/ρ_i = 1000/917) | 54,526 |
| Basal melt, G = 0.06 W/m² (m/Myr) | 0.113 | 6,182 (G / ρ_i L_f) | 54,526 |
| `WATER_TO_ICE` (m ice per m water) | 2e-5 | 1.0905 | 54,526 |
| Surface melt (m ice / Myr per °C of year-averaged positive degrees) | 20 | 2.79 × 10⁶ | 139,408 |
| Abrasion (m / Myr per m/yr of sliding) | 430 (calibrated) | 100 (ė = 10⁻⁴ u_b) | — |
| Meltwater routed per m³ of ice | 1 m³ | 0.917 m³ (by mass) | 0.917 |

Physical surface melt uses a degree-day factor of 7 mm water equivalent per °C per day, which
falls inside Hock's (2003) range for ice and above the snow range. One °C of year-averaged
positive degrees is 365.25 positive degree-days per year. Physical melt is therefore 2.56×
stronger relative to accumulation than the model's current ratio, so a full conversion is
not the same as scaling both by one factor. The abrasion value is the linear law of Humphrey
& Raymond (1994), the upper end of the coefficients summarised by Cook et al. (2020).

## Design (preregistered)

All arms run the same backend (`origin/main` at `30766af`, branch
`claude/issue-368-glacier-mass-balance`, with no model changes). `bin/debug/measure_ice_sheet.py`
patches the constants in-process, so arm `k1` is the unmodified model and serves as the
baseline. No baseline worktree is needed.

Saves (SHA-256):
- seed 997271774 at 87.4 Myr, `08ad9461…ed762b` (#353's abrasion calibration replay);
- seed 649658380 at 6.4 Myr, `ebbfb5da…73f7b3` (#503's young, glaciating world).

Each arm replays 10 Myr in 0.25 Myr steps with 1 Myr checkpoints.

| Arm | Accumulation, basal, `WATER_TO_ICE` × | Surface melt × | Abrasion | Meltwater by mass |
|---|---:|---:|---:|---|
| `k1` | 1 | 1 | 430 | no |
| `k10` | 10 | 10 | 430 | no |
| `k100` | 100 | 100 | 430 | no |
| `k1000` | 1,000 | 1,000 | 430 | no |
| `kphys` | 54,526 | 54,526 | 430 | no |
| `phys-neutral` | 54,526 | 139,408 | 430 / 54,526 | yes |
| `phys` | 54,526 | 139,408 | 100 | yes |

The `k*` arms test proportional scaling: the same balance points, reached faster. The
`phys-neutral` arm shows a full mass-balance conversion with erosion rescaled to keep the
current calibration. The `phys` arm adds physical abrasion. Step check: `kphys` and `phys`
also run on seed 997271774 at 1 Myr steps.

Runs are serial and wait until no other replay, probe, unit or stress run is active.

### Predictions

1. **Time to balance.** At k = 1, typical accumulation (~6 m/Myr) needs hundreds of Myr to
   build a kilometres-thick sheet, so 10 Myr of `k1` is transient. From k ≈ 100 up, the ice
   should reach equilibrium within the window. Prediction: ice-covered land area at 10 Myr
   differs by < 10% between `k1000` and `kphys`, and by > 10% between `k1` and `kphys`.
2. **Meltwater** in balance equals accumulation, so it should scale roughly with k: `kphys`
   routed meltwater at 10 Myr is within a factor of 3 of 54,526 × `k1`'s accumulation.
3. **Erosion.** Raw abrasion scales with sliding, so it should scale roughly with k. If the
   Hc-headroom and detachment limits bind, applied abrasion grows less than raw. Prediction:
   applied/raw falls below 0.5 at `kphys`.
4. **Depth guard.** Cap steps and shed volume rise with k.
5. **Physical melt.** `phys-neutral` has less ice area than `kphys`, since its melt is 2.56×
   stronger relative to accumulation.
6. **Step size.** `kphys` and `phys` at 1 Myr steps match their 0.25 Myr runs within 5% in ice
   area and volume.

Primary metrics: ice-covered land area and ice volume at 10 Myr. Secondary: sea level,
routed meltwater, ocean inflow, raw and applied glacial erosion, cap steps, cap nodes and
shed volume, ledger residual, step time and peak RSS.

### Decision rule

An RFC is warranted if, on either save, `phys-neutral` or `phys` moves ice area, ice
volume, sea level or applied glacial erosion by more than 25% from `k1` at 10 Myr, **or** if
the depth guard sheds more than 10% of routed meltwater, **or** if applied/raw abrasion falls
below 0.5 (the caps are doing the model's work). Otherwise the conversion is an ordinary
issue: change the constants, keep the balance points, and recheck the guard.

## Results

Pending.
