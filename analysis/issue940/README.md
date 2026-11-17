# Issue 940: RFC 0015 gate G3, glacial erosion and its limits

Gate G3 of [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md#diagnosis-gates) for
[#940](https://github.com/adubey/mantle-bloom/issues/940). Tracking:
[#544](https://github.com/adubey/mantle-bloom/issues/544).

**Principles:** #5 (consistent units) and #8 (each process follows its own physics). Replacing
`origin/main`'s model-unit abrasion coefficient (`430 * MODEL_UNIT_SCALE`) with physical
literature rates (`K_G = 1e-4`, `l = 0.65` or `1.0`) increases raw abrasion by three orders of
magnitude, causing `erosion.py`'s shared per-step cap and scaling (`applied_scale`) to remove
79–99.5% of raw abrasion and creating a 66–73% step-size gap in applied glacial erosion between
0.25 Myr and 1.0 Myr steps.

**Verdict: D3's within-step bound must be designed explicitly before implementation.** Under the
G3 decision rule (`applied / raw < 0.5`), both `(K_G, l)` pairs fail the clip-ratio threshold on
every deciding save and step size (`applied / raw` = `0.048–0.210` for Cook et al. 2020;
`0.0045–0.049` for Humphrey & Raymond 1994). The `l = 2` option (Herman et al., 2015) was
dropped because its coefficient was not verified directly from the paper.

## Method and runs

Replays were run with [gate_g3_glacial_erosion.py](../../bin/debug/gate_g3_glacial_erosion.py)
and [run_940_g3.sh](../../bin/debug/run_940_g3.sh) on `origin/main` (`bf624090`, default
`trend-1k` march, ice ages off) for 10 Myr total at 0.25 Myr (40 steps) and 1.0 Myr (10 steps):

- **Deciding saves:** `mantle-bloom-seed649658380-6400000y.bloom` (6.4 Myr, peak RSS 1.11–1.22 GB)
  and `mantle-bloom-seed997271774-87400000y.bloom` (87.4 Myr, peak RSS 1.38–1.68 GB).
- **Arms (`glacier_erosion_multiplier = 1.0`):**
  - `main`: `GLACIER_ABRASION_COEFFICIENT = 430 * MODEL_UNIT_SCALE` (`≈ 0.007886` m/Myr per `(m/yr)`), `l = 1.0`
  - `cook2020`: `K_G = 1e-4` m/yr per `(m/yr)^0.65` (`100.0` m/Myr), `l = 0.65` (Cook et al., 2020)
  - `humphrey1994`: `K_G = 1e-4` m/yr per `(m/yr)^1.0` (`100.0` m/Myr), `l = 1.0` (Humphrey & Raymond, 1994)
- **Knob sensitivity (`cook2020`, 1.0 Myr × 10 steps):** `World.glacier_erosion_multiplier` at
  `0.25×` and `4×`, both full (abrasion + flattening) and abrasion-only (`flatten_multiplier = 1.0`).

Full per-step JSONs and logs are in `~/mantle-bloom-captures/issue-940/` with `SHA256SUMS`;
compact summary metrics are in [`g3.json`](g3.json).

## Results

### 1. Abrasion volume, clip ratio, step gap, and flattening (10 Myr)

Note: `applied / raw` is the effective abrasion ratio from `erosion.py`'s shared
`applied_scale = erosion_amount / raw_erosion_total`, which reflects the combined per-step
terrain/headroom cap (`min(max(drop_to_lowest_neighbor_m, march_cut_m), removable_m)`),
mobile-cover erodibility (`_strip_mobile_cover`), craton resistance (`craton_keep`), and shared
competition with the other subaerial erosion terms (`rain + river + weathering + seismic`), not
the neighbor-drop clip in isolation. On `main`, `applied / raw > 1` (`1.97–4.10`) because raw
abrasion is small enough that the cap rarely binds while mobile-cover erodibility amplifies
cover entrainment up to `MOBILE_COVER_ERODIBILITY_FACTOR`.

| Save | Arm | Step (Myr) | Raw abrasion (km³) | Applied abrasion (km³) | Applied / Raw | Step gap (applied) | Flattening (km³) | Flatten / Applied |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| 649658380 | `main` | 0.25 | 4,217 | 17,274 | 4.0964 | 0.2% | 1,982,806 | 114.78 |
| 649658380 | `main` | 1.0 | 4,384 | 17,304 | 3.9469 | — | 1,650,974 | 95.41 |
| 649658380 | `cook2020` | 0.25 | 12,977,604 | 2,720,198 | **0.2096** | **66.5%** | 1,715,729 | 0.63 |
| 649658380 | `cook2020` | 1.0 | 12,833,599 | 912,146 | **0.0711** | — | 1,567,935 | 1.72 |
| 649658380 | `humphrey1994` | 0.25 | 57,074,332 | 2,774,295 | **0.0486** | **65.7%** | 1,749,997 | 0.63 |
| 649658380 | `humphrey1994` | 1.0 | 66,930,332 | 951,464 | **0.0142** | — | 1,567,811 | 1.65 |
| 997271774 | `main` | 0.25 | 169,374 | 457,849 | 2.7032 | 23.6% | 16,240,047 | 35.47 |
| 997271774 | `main` | 1.0 | 177,328 | 349,649 | 1.9718 | — | 12,203,415 | 34.90 |
| 997271774 | `cook2020` | 0.25 | 216,012,038 | 33,883,494 | **0.1569** | **70.4%** | 16,298,441 | 0.48 |
| 997271774 | `cook2020` | 1.0 | 211,032,797 | 10,041,063 | **0.0476** | — | 12,228,636 | 1.22 |
| 997271774 | `humphrey1994` | 0.25 | 2,258,460,285 | 39,875,153 | **0.0177** | **72.6%** | 16,266,853 | 0.41 |
| 997271774 | `humphrey1994` | 1.0 | 2,420,986,990 | 10,922,151 | **0.0045** | — | 12,297,338 | 1.13 |

Key observations:
1. **Clip saturation:** Under `cook2020`, the neighbor-drop clip removes **79.0–92.9%** of raw
   abrasion on `649658380` and **84.3–95.2%** on `997271774`. Under `humphrey1994`, it removes
   **95.1–99.5%**. Because both laws saturate the same per-step neighbor-drop cap, their applied
   erosion volumes are nearly identical (`2.72M` vs `2.77M km³` on `649658380`; `33.88M` vs
   `39.88M km³` on `997271774` at 0.25 Myr) despite `humphrey1994` having 4.4–10.5× higher raw
   demand.
2. **Step-size dependence:** Because the neighbor-drop cap resets each step, taking 4× as many
   steps (0.25 Myr vs 1.0 Myr) increases applied glacial erosion by **~3×** (`66–73%` step gap)
   under both physical laws. By contrast, closing `land_ice_area_km2` and
   `land_volume_above_sea_km3` agree within 0.1–8.5% and 1.2–1.7%, respectively.
3. **Comparison with `origin/main`:** At 0.25 Myr, `cook2020` applied abrasion is **157×** `main`
   on `649658380` and **74×** `main` on `997271774`.
4. **Flattening beside abrasion:** Under `main`, flattening is **35–115×** applied abrasion. Under
   `cook2020`, flattening is comparable to applied abrasion (**0.48–0.63×** at 0.25 Myr;
   **1.22–1.72×** at 1.0 Myr).

### 2. `World.glacier_erosion_multiplier` sensitivity (`cook2020`, 1.0 Myr × 10 steps)

| Save | Multiplier | Flattening mult | Land ice area (M km²) | Δ Ice area | Land vol (M km³) | Δ Land vol | Applied abrasion (km³) | Flattening (km³) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 649658380 | 0.25× | 0.25× | 2.98 | +2.18% | 89.73 | +0.18% | 872,805 (-4.3%) | 431,225 (-72.5%) |
| 649658380 | 0.25× | 1.00× | 2.91 | -0.37% | 89.58 | +0.01% | 830,691 (-8.9%) | 1,582,900 (+1.0%) |
| 649658380 | 1.00× | 1.00× | 2.92 | 0.00% | 89.57 | 0.00% | 912,146 (0.0%) | 1,567,935 (0.0%) |
| 649658380 | 4.00× | 1.00× | 2.97 | +1.80% | 89.59 | +0.01% | 933,856 (+2.4%) | 1,559,583 (-0.5%) |
| 649658380 | 4.00× | 4.00× | 2.72 | -6.89% | 89.24 | -0.37% | 860,587 (-5.7%) | 4,893,319 (+212.1%) |
| 997271774 | 0.25× | 0.25× | 13.72 | -3.39% | 260.59 | +0.02% | 7,936,259 (-21.0%) | 3,188,320 (-73.9%) |
| 997271774 | 0.25× | 1.00× | 13.66 | -3.83% | 260.40 | -0.06% | 8,075,538 (-19.6%) | 12,249,564 (+0.2%) |
| 997271774 | 1.00× | 1.00× | 14.20 | 0.00% | 260.55 | 0.00% | 10,041,063 (0.0%) | 12,228,636 (0.0%) |
| 997271774 | 4.00× | 1.00× | 13.87 | -2.32% | 260.50 | -0.02% | 11,379,878 (+13.3%) | 12,323,524 (+0.8%) |
| 997271774 | 4.00× | 4.00× | 14.00 | -1.42% | 260.13 | -0.16% | 12,346,522 (+23.0%) | 44,488,158 (+263.8%) |

- **Flattening's share:** With abrasion saturated by the neighbor-drop clip, `glacier_erosion_multiplier`
  barely moves applied abrasion on `649658380` (`-4.3%` to `-5.7%` across a 16× span) and moves it
  modestly on `997271774` (`-21%` to `+23%`). Almost all of the knob's effect on **land volume**
  (`+0.18%` to `-0.37%` on `649658380`, `+0.02%` to `-0.16%` on `997271774`) comes from
  **flattening** (with flattening held at `1×`, land volume moves by `≤ 0.06%`). On `649658380`,
  flattening also drives the `-6.89%` ice-cap area reduction at `4×`; on `997271774`, ice-cap area
  shifts (`-1.4%` to `-3.8%`) are non-monotonic and within margin root-find variance.
