# Issue 917: RFC 0015 D6, sea level: debit ice above flotation

Implementation of [RFC 0015 D6](../../rfc/0015-frozen-water-in-physical-units.md#d6-the-ocean)
(sea level half) for [#917](https://github.com/adubey/mantle-bloom/issues/917). Tracking:
[#544](https://github.com/adubey/mantle-bloom/issues/544).

**Principle:** #1, conservation. Sea level comes from a conserved water budget, so every cubic
metre booked as ice must be water that has left the ocean.

## Pre-registered comparison

- **Saves:** the [#548](../issue548/README.md) diagnosis saves (converted to `.bloom` in B2):
  - `649658380` at 6.4 Myr (`mantle-bloom-seed649658380-6400000y.bloom`, deciding, ice ages off)
  - `997271774` at 87.4 Myr (`mantle-bloom-seed997271774-87400000y.bloom`, deciding, ice ages off)
  - `374464662` at 132.6 Myr (`mantle-bloom-seed374464662-132600000y.bloom`, report only, 8 Myr ice-age period)
- **Protocol:** 10 steps at 0.25 Myr (primary) and 1 Myr (secondary) with the default `trend-1k`
  glacier march, comparing `origin/main` (`011f288d`) against the branch.
- **Primary metric & threshold:**
  - On the deciding saves, the change in sea level from `origin/main` is expected to be positive
    and small (within G4's measured over-debit of ≤ 0.08 m, plus coupled-run drift).
  - After the change, the replay's own `sea_level_model_m` matches `measure_flotation_debit.py`'s
    generalised arm `h_generalised_m` (to within eustasy's solver tolerance, ~1 m).
  - Land and glacier water ledger residuals remain closed / unchanged relative to `origin/main`.
  - Ice volume is reported alongside sea level per step and at the end.
