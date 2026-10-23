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

Physical surface melt uses a degree-day factor of 7 mm water equivalent per °C per day. That
is above the 3.5 ± 1.4 and 4.6 ± 1.4 mm d⁻¹ K⁻¹ Braithwaite (2008, doi:10.3189/002214308785836968)
finds for melting snow at 66 glaciers' equilibrium lines, as an ice factor should be, but the
value itself is a design choice: it is not checked against a published ice range here (RFC
0015's gate G2 sources it). One °C of year-averaged positive degrees is 365.25 positive
degree-days per year. Physical melt is therefore 2.56× stronger relative to accumulation than
the model's current ratio, so a full conversion is not the same as scaling both by one
factor. The abrasion value is the dimensionless erosion rate of order 1.0 × 10⁻⁴ that Humphrey
& Raymond (1994, doi:10.3189/S0022143000012429) measured on Variegated Glacier, and the mean
coefficient of Cook et al.'s (2020, doi:10.1038/s41467-020-14583-8) 38-glacier compilation.

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

All 16 runs completed (exit 0) on `30766af`, one at a time. The numbers are in
[`compare.json`](compare.json): per-checkpoint ice volume, ice cover and sea level for every
run, plus ledger flows, meltwater, erosion and depth-guard counts. Neither save has ice-age
forcing (`ice_age_period_years` = 0), so every change in the ice comes from inside the model.
The ice ledger closes in every run (worst residual 1.75 × 10⁵ m³ per step against ice volumes
of 10¹⁵–10¹⁸ m³).

Command, per arm (flags from the design table):

```bash
cd backend
.venv/bin/python ../bin/debug/measure_ice_sheet.py <local-save> --myr 10 --step-years 250000 \
    --checkpoint-myr 1 [--accumulation-scale K] [--melt-scale M] \
    [--abrasion-coefficient A] [--meltwater-density] --out <scratch-output>
```

### The ice does not settle at physical rates

Ice volume (million km³) at each 1 Myr checkpoint:

| Run | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 997 `k1` | 5.7 | 5.6 | 5.6 | 5.7 | 5.8 | 5.8 | 5.9 | 5.9 | 6.0 | 6.0 | 6.0 |
| 997 `k100` | 5.7 | 11.1 | 16.1 | 17.5 | 18.1 | 18.6 | 18.4 | 19.2 | 18.9 | 18.9 | 18.3 |
| 997 `k1000` | 5.7 | 75 | 148 | 208 | 285 | 354 | 438 | 506 | 525 | 719 | 488 |
| 997 `kphys` | 5.7 | 748 | 115 | 177 | 1,016 | 349 | 524 | 73 | 156 | 105 | 125 |
| 997 `phys` | 5.7 | 158 | 447 | 113 | 491 | 665 | 111 | 129 | 78 | 101 | 657 |
| 649 `k1` | 3.6 | 3.6 | 3.6 | 3.6 | 3.6 | 3.6 | 3.6 | 3.7 | 3.7 | 3.7 | 3.8 |
| 649 `k1000` | 3.6 | 4.5 | 4.2 | 4.6 | 4.2 | 4.4 | 5.4 | 4.6 | 5.1 | 4.7 | 4.6 |
| 649 `kphys` | 3.6 | 94 | 15.5 | 5.1 | 89 | 23 | 9.9 | 85 | 73 | 37 | 86 |
| 649 `phys` | 3.6 | 91 | 18.0 | 5.1 | 86 | 23 | 7.5 | 84 | 26 | 10.2 | 83 |

- Up to k = 100 the ice behaves as predicted: on seed 997271774 it grows to a new level
  (about 18–19 M km³, 13% of land) by 4 Myr and stays there. `k1`'s ice is still the save's
  inherited ice, far from its own balance.
- From k = 1,000 it doesn't settle on seed 997271774, where `k1000` grows without levelling
  off (seed 649658380's `k1000` holds near 4.6 M km³). At physical rates,
  volume swings by 10–15× between checkpoints, and sea level by up to 3 km. On seed 649658380
  the swing repeats about every 3 Myr.
- **The result depends on step size.** On seed 997271774, `kphys` ends at 125 M km³ with
  0.25 Myr steps and 927 M km³ with 1 Myr steps, and `phys` at 657 and 154 M km³. A
  physical-rate sheet takes about 10 kyr to respond (a few km of ice at ~300 m/kyr of
  accumulation), far shorter than the step. Sea level, the climate's ice surface and frozen
  river ice all reach the ice only once per step, a step late. The swings look like the
  overshoot of those lagged couplings, not a physical cycle; this analysis doesn't isolate
  which coupling drives them.

### 10 Myr means (checkpoints 1–10)

Seed 997271774:

| Run | Ice (M km³) | Ice cover (% of land) | Sea level (m) | Meltwater (km³/yr) | Ocean inflow (m³/s) | Applied glacial erosion (km³) | Applied / raw | Guard-clipped share of meltwater |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `k1` | 5.9 | 5.3 | −508 | 0.38 | 920,000 | 3.16 M | 1.42 | 0 |
| `k10` | 8.6 | 7.6 | −514 | 3.9 | 922,000 | 4.79 M | 1.50 | 4% |
| `k100` | 17.5 | 12.4 | −536 | 44 | 909,000 | 16.9 M | 0.85 | 6% |
| `k1000` | 375 | 25.1 | −1,559 | 2,600 | 779,000 | 60.8 M | 0.28 | 2% |
| `kphys` | 339 | 16.4 | −1,455 | 134,000 | 2,090,000 | 30.0 M | 0.009 | 17% |
| `phys-neutral` | 118 | 14.3 | −818 | 52,000 | 1,218,000 | 0.028 M | 2.68 | 35% |
| `phys` | 295 | 16.4 | −1,307 | 60,000 | 1,277,000 | 5.46 M | 0.040 | 31% |
| `kphys`, 1 Myr steps | 498 | 21.3 | −1,938 | 144,000 | 2,423,000 | 10.4 M | 0.004 | 15% |
| `phys`, 1 Myr steps | 219 | 17.5 | −1,108 | 65,000 | 1,326,000 | 2.01 M | 0.011 | 30% |

Seed 649658380:

| Run | Ice (M km³) | Ice cover (% of land) | Sea level (m) | Meltwater (km³/yr) | Ocean inflow (m³/s) | Applied glacial erosion (km³) | Applied / raw | Guard-clipped share of meltwater |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `k1` | 3.7 | 2.8 | −146 | 0.33 | 732,000 | 0.87 M | 1.57 | 0 |
| `k10` | 3.8 | 3.0 | −146 | 3.4 | 736,000 | 1.12 M | 1.56 | 5% |
| `k100` | 3.4 | 2.9 | −145 | 33 | 738,000 | 2.33 M | 1.54 | 9% |
| `k1000` | 4.6 | 7.5 | −149 | 325 | 729,000 | 5.12 M | 0.36 | 11% |
| `kphys` | 51.8 | 5.7 | −299 | 20,000 | 866,000 | 6.49 M | 0.007 | 45% |
| `phys-neutral` | 43.0 | 5.6 | −276 | 16,000 | 819,000 | 0.005 M | 3.73 | 57% |
| `phys` | 43.3 | 5.6 | −277 | 16,000 | 816,000 | 1.02 M | 0.067 | 57% |

Applied / raw exceeds 1 at low k because loose sediment cover erodes up to
`MOBILE_COVER_ERODIBILITY_FACTOR` times faster than rock (erosion's `applied_scale`). The
script's raw is the abrasion law alone. The guard-clipped share is ice shed at the depth
guard (`capped_m3`) plus frozen river ice clipped there (`frozen_river_capped_m3`), over
routed meltwater. Ice shed alone stays below 0.01% of meltwater in every run, and the clipped
frozen river ice is almost all of it.

### What it shows

- **The balance points don't survive the conversion.** On every measure the decision rule
  names, `phys-neutral` and `phys` land far outside 25% of `k1`: mean ice volume 20–51× `k1`'s on
  seed 997271774 and 12× on seed 649658380, mean sea level 130–800 m lower, ice cover 2–3×.
  On seed 997271774 the converted ice peaks at 665 M km³ (`phys`), about 8× the Last Glacial
  Maximum's ~80 M km³. Since the result also changes with step size, no
  single number here is "the" physical ice volume.
- **The caps do the work.** Raw abrasion with the current coefficient rises about 1,500× at
  `kphys` while applied rises only 9×: the neighbor-drop and Hc-headroom limits clip 99% of
  it (applied / raw 0.007–0.009). Physical abrasion (`phys`) is clipped 93–96%. With the
  current calibration preserved (`phys-neutral`), glacial erosion collapses by 100–200×.
  Sliding grew about 120–260× under `phys-neutral`'s mass balance and 1,500× under `kphys`'s,
  not 54,526×, so dividing the coefficient by the accumulation factor over-corrects.
- **Meltwater grows faster than k.** At physical rates it is 15–17× more than k × `k1`'s
  accumulation, because the ice covers more ground and snow melting on bare ground (60–96% of
  meltwater) grows with it. Ocean inflow rises 11–127%.
- **The depth guard clips frozen river ice, not sheet ice.** The physical `WATER_TO_ICE`
  makes rivers that freeze build kilometres of ice in a step, and the guard clips back that ice
  at 17–57% of meltwater's volume. This crosses into RFC 0014's ground (frozen
  river ice).
- **A water leak today.** At frozen nodes `compute_hydrology` drops all precipitation from
  physical runoff (`liquid_precip`), but the ice gains only 1/54,526 of it. On `k1` that is
  roughly 7,100 km³/yr (about 225,000 m³/s, a quarter of ocean inflow) on seed 997271774 and
  about 1,200 km³/yr on seed 649658380, missing from the river budget. Sea level is
  unaffected, since eustasy counts the stored ice. This bears on #355.

### Predictions

| # | Prediction | Outcome |
|---|---|---|
| 1 | `k1000` vs `kphys` ice cover within 10%; `k1` vs `kphys` more than 10% apart | First not met (25.4% vs 12.6% and 7.1% vs 10.5% at 10 Myr); second met |
| 2 | `kphys` meltwater within 3× of 54,526 × `k1` accumulation | Not met: 15–17× |
| 3 | Applied/raw below 0.5 at `kphys` | Met: 0.007–0.009 |
| 4 | Cap steps and shed volume rise with k | Met: 0–2 of 40 steps at `k1`, 40 of 40 at `kphys` |
| 5 | `phys-neutral` has less ice cover than `kphys` | Met on seed 997271774 (14.3% vs 16.4% mean); within the swings on seed 649658380 (5.6% vs 5.7%) |
| 6 | 1 Myr steps within 5% of 0.25 Myr | Not met: 4–7× apart in volume at 10 Myr |

### Recommendation

**Write an RFC.** The rule's ice-volume, sea-level and applied/raw criteria are all met. The
guard criterion is met only through frozen river ice; ice shed from sheets is negligible.
Changing the constants alone gives a mean ice volume of 43–498 M km³ that changes with step size. The RFC would
need:

1. a coupling that stays stable when the ice responds in ~10 kyr: sub-stepping the ice
   against sea level and the ice-surface climate, or solving for the equilibrium sheet
   directly;
2. a physical abrasion law and the role of the neighbor-drop and Hc-headroom limits, which
   currently clip 93–99% of it;
3. frozen river ice under physical `WATER_TO_ICE`, settled with RFC 0014;
4. the snowfall water leak at frozen nodes, which a physical exchange fixes and which #355
   depends on.

### Runtime

Mean step time was 4.9–8.9 s, and peak RSS 1.6–1.9 GB per run. Seed 997271774's `k1`,
`kphys`, `phys-neutral`, `phys` and `k10` ran while another replay and a 1.6 GB probe were
running (the user cleared the machine). Their step times are not from an idle machine. The
other 11 runs started only after the queue's busy check found no other replay, probe or test
running.
