# RFC 0014 D3: water-surface temperature

## Before/after replay preregistration

**Question.** Does applying the previous liquid lake level's hypsometric water-surface elevation change lake ice and water-budget outcomes in a production step, while preserving the existing lake water ledger bounds?

**Principle and invariant.** The previous behavior used bed-level temperature for seasonal lake losses and ice phase changes, contrary to design principles 1 (state the model boundary) and 3 (aggregate seasonal processes over long steps). The code-level invariant is checked by `test_lake_temperature_uses_lagged_water_surface_and_connected_ocean_mask`: enclosed below-sea-level lakes use lapse-adjusted air temperature at the prior water surface; only connected ocean uses ocean temperature; a newly wet cell has no prior-level adjustment.

**Comparison.** Replay the preregistered held-out save for seed 997271774 (start age 130.5 Myr) for ten 100 kyr steps using `bin/debug/gate_lake_ice_heldout.py`, once at the pre-change `origin/main` commit and once on the D3 branch. The save SHA-256 is `a43962a149a0ec96f1dbee3e6c675be98d50422c47e0455e062aa72a7d1da4d7`. Set `ice_age_period_years = 0` as the driver requires. Run each arm alone, serially, with no tests or other replay running.

**Expected direction.** Raising the temperature reference from the bed to a lake's water surface applies cooling with elevation, which should reduce ice-free fraction and open-season evaporation in previously wet cells. The resulting lake ice volume may increase and liquid storage/overflow may change. The primary diagnostic is the endpoint and per-step lake ice / liquid state; ledger residual must remain within the script's `0.1%` positive-supply acceptance bound, physical ice must remain ≤2 m, and ice must remain within local liquid-plus-ice inventory. This is a directional diagnostic, not a calibrated target; no pass threshold is set for the magnitude of the difference.

**Noise and limits.** This is one held-out save and a short 1 Myr interval, so it does not estimate seed variability or long-term behavior. The saved world represents the same baseline for both arms; no parameters are tuned from this comparison. Full replay output stays in scratch storage. The checked-in compact summary will record both commits, checksum, command, endpoint values, runtime and peak RSS.
