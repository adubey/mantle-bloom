# RFC 0014 D3: water-surface temperature

## Before/after replay preregistration

**Question.** Does applying the previous liquid lake level's hypsometric water-surface elevation change lake ice and water-budget outcomes in a production step, while preserving the existing lake water ledger bounds?

**Principle and invariant.** The previous behavior used bed-level temperature for seasonal lake losses and ice phase changes, contrary to design principles 1 (state the model boundary) and 3 (aggregate seasonal processes over long steps). The code-level invariant is checked by `test_lake_temperature_uses_lagged_water_surface_and_connected_ocean_mask`: enclosed below-sea-level lakes use lapse-adjusted air temperature at the prior water surface; only connected ocean uses ocean temperature; a newly wet cell has no prior-level adjustment.

**Comparison.** Replay the preregistered held-out save for seed 997271774 (start age 130.5 Myr) for ten 100 kyr steps using `bin/debug/gate_lake_ice_heldout.py`, once at the pre-change `origin/main` commit and once on the D3 branch. The save SHA-256 is `a43962a149a0ec96f1dbee3e6c675be98d50422c47e0455e062aa72a7d1da4d7`. Set `ice_age_period_years = 0` as the driver requires. Run each arm alone, serially, with no tests or other replay running.

**Expected direction.** Raising the temperature reference from the bed to a lake's water surface applies cooling with elevation, which should reduce ice-free fraction and open-season evaporation in previously wet cells. The resulting lake ice volume may increase and liquid storage/overflow may change. The primary diagnostic is the endpoint and per-step lake ice / liquid state; ledger residual must remain within the script's `0.1%` positive-supply acceptance bound, physical ice must remain ≤2 m, and ice must remain within local liquid-plus-ice inventory. This is a directional diagnostic, not a calibrated target; no pass threshold is set for the magnitude of the difference.

**Noise and limits.** This is one held-out save and a short 1 Myr interval, so it does not estimate seed variability or long-term behavior. The saved world represents the same baseline for both arms; no parameters are tuned from this comparison. Full replay output stays in scratch storage. The checked-in compact summary will record both commits, checksum, command, endpoint values, runtime and peak RSS.

## Held-out replay results

The same held-out save was replayed serially on `origin/main` (`8a69bab69c32114c9a42a321d0b80c2db3d9c393`) and the D3 implementation (`a1ab8eb1df558a3ef495308604781993983101e2`). Both arms ran the exact ten-step 100 kyr held-out driver; the baseline completed in 64.22 s (peak RSS 2,145 MB) and the D3 arm in 65.10 s (peak RSS 1,938 MB).

| Endpoint at 131.5 Myr | `origin/main` | D3 | Change |
|---|---:|---:|---:|
| Lake ice water equivalent | 108.278 km³ | 238.329 km³ | +130.050 km³ (+120.1%) |
| Liquid lake volume | 250,172.190 km³ | 265,116.883 km³ | +14,944.692 km³ (+6.0%) |
| Wet lake area | 1,242,839.933 km² | 1,288,309.427 km² | +45,469.494 km² (+3.7%) |
| Maximum physical ice | 1.999 m | 1.905 m | -0.093 m (-4.7%) |
| Shallow bedfast area | 0.000 % | 0.756 % | +0.756 % |
| Deep bedfast area | 0.000 % | 0.000 % | +0.000 % |
| Max per-step ledger residual / positive supply | 1.1e-10% | 1.9e-10% | both pass 0.1% gate |

Both arms passed all held-out checks: physical ice ≤2 m (within floating-point tolerance), nonnegative reservoirs, ice within local liquid-plus-ice inventory, deep bedfast area ≤1%, and lake ledger residual ≤0.1% of positive supply. The observed ice-water-equivalent volume rose by about 111% on this one save. This is a clear model response on the tested trajectory, not a calibrated estimate of typical change; one seed and one 1 Myr interval do not quantify seed-to-seed or long-term variability. The per-step series, both commits, checksums, runtime, RSS and acceptance results are in [`d3-replay.json`](d3-replay.json). Full driver outputs remain in scratch storage.

Reproduction command (replace `<local-save>` and `<scratch-output>` with local paths):

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/gate_lake_ice_heldout.py \
  <local-save> --myr 1 --step-years 100000 --checkpoint-myr 0.1 \
  --out <scratch-output>
```
