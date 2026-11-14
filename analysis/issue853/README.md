# Issue #853: lake solver result against slice length

Status: first pass. One save, one step from the start state, one anchor. Not a full answer to #853.

Driver: [`bin/debug/measure_lake_slice_length.py`](../../bin/debug/measure_lake_slice_length.py) on
`origin/main` `b72e4461`. Save: `seed910211954-225.1myr.bloom` (RFC 0016 window D start, 225.1 Myr;
bucket `legacy_sha256` `6a8d6f90...`, the registered hash), ice-age cycle off, anchor
[`target-anchor-window-d.json`](../rfc0016/target-anchor-window-d.json). Raw data:
`fixed-5k.json`, `fixed-25k.json`, `refreshed-5k.json` here (full captures and `SHA256SUMS` in
`~/mantle-bloom-captures/issue-853/`).

## Method

`hydrology` calls `lakes.solve_lakes_prognostic` with the default `substeps=1`, so the lake slice is
the whole step. The driver reloads the save per row and takes one `step_world`.

- **Fixed climate:** the driver forces `substeps = ceil(step / slice)`. Climate, rivers and terrain are
  computed before the solve, so they are identical in every row; only the lake integration is sliced.
- **Refreshed climate:** the same 5 kyr span run as `5000 / slice` ordinary steps, so everything is
  recomputed every slice. This changes more than climate (plates, erosion, ice, rivers also advance
  per slice), so it bounds the refresh effect rather than isolating climate.
- **Anchor lake:** the lake component at the anchor cell, or the nearest lake when the anchor is dry
  (same rule as `replay_rfc0016_river_water.py`). The anchor is dry in every row. "Merged" means that
  lake covers over 3M km2 (the issue's merged state is 6 to 9M km2; the small basin is under 1M).

## Anchor lake volume after one step (km3)

Fixed climate, 5 kyr step:

| Slice (yr) | 5000 | 2500 | 1000 | 500 | 200 | 100 |
|---|---|---|---|---|---|---|
| Anchor lake km3 | 348,607 | 356,508 | 366,435 | 381,263 | 389,495 | 392,063 |
| Change from previous | | +2.3% | +2.8% | +4.0% | +2.2% | +0.7% |

Fixed climate, 25 kyr step:

| Slice (yr) | 25000 | 5000 | 1000 | 200 |
|---|---|---|---|---|
| Anchor lake km3 | 392,527 | 408,890 | 410,447 | 410,649 |

Climate refreshed, 5 kyr span:

| Slice (yr) | 5000 | 2500 | 1000 | 500 | 200 |
|---|---|---|---|---|---|
| Anchor lake km3 | 348,607 | 357,277 | 410,823 | 490,308 | 541,291 |
| Lake area M km2 | 0.26 | 0.26 | 0.28 | 0.33 | 0.34 |
| World lake km3 | 3.61M | 4.37M | 4.79M | 4.85M | 4.99M |

Lake ledger residuals are 2e-12 or smaller (normalized) in every row.

## Findings

1. **At fixed climate the lake solver converges.** The anchor lake rises 12.5% from one slice to 100 yr
   slices at a 5 kyr step, with the change shrinking to 0.7% from 200 to 100 yr; the 25 kyr step
   moves 0.05% from 1,000 to 200 yr slices. Slices of about 500 yr or shorter are within about 3%
   of the limit at 5 kyr, and 1,000 yr is within 0.5% at 25 kyr. So the slice length is a real but
   modest effect, and it is not what makes results swing.
2. **Refreshing climate every slice matters more and is not converged at 200 yr.** The anchor lake
   grows 55% from 5000 to 200 yr, still rising 10% from 500 to 200 yr, and world lake volume grows
   38%. This arm also changes terrain, ice and rivers per slice, so this run does not say which
   input drives it.
3. **The issue's large spread and the merged state do not reproduce on current main.** The issue's
   0.29M to 2.29M km3 values and the merged inland sea come from the pinned RFC 0016 source; #811
   already found they do not occur on current main (#757, #758). Here the anchor lake is
   0.35M to 0.54M km3 over 0.26M to 0.34M km2 in every row, never merged. So whether D/E merged
   survives at converged slices cannot be answered from this save: there is no merged state to test.
4. **Not examined:** more than one step, other anchors or saves (windows E and A), whether the
   slice-length effect grows over a multi-step run, and why 10 kyr did not merge in the old data.

## Recommendation

No change to the lake solver's slicing for convergence: it is converged at fixed climate at slices
of 500 yr and below, at a cost of linear runtime in slices (65 s against 51 s at 50 slices). The
open dependence on step length comes from the refresh interval, so the next issue is to separate
which per-slice input (climate, rivers and river storage, ice, terrain) drives the growth in the
refreshed arm, starting from this driver, and to replace the D/E windows as #811 already recommends.
Whether #769's drain is a step artefact is not answered here; the owner decides on closing it.
