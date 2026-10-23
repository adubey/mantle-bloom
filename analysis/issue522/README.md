# Issue #522: insolation in W/m² and humidity in g/kg

## Preregistration

Written and committed before any run.

**Change.** Insolation becomes annual-mean top-of-atmosphere flux in W/m². Temperature is
interpolated through the existing temperature-by-latitude calibration. Humidity becomes
near-surface specific humidity in g/kg: 80% of Clausius–Clapeyron saturation over open water,
capped at 16.2 g/kg. Every constant in the old unitless humidity index was converted at its
reference point, where 1.0 is air over 20 °C water, or 11.56 g/kg.

**Design principle.** Explicit physical units (principle 5). The change is meant to make the
fields readable, not to move the model. The metrics below test that it doesn't move it much.

**Runs.** `bin/debug/measure_climate_fields.py`, 10 Myr in 250 kyr steps, checkpoints every
2 Myr. Two saves:

| Save | SHA-256 |
|---|---|
| seed 997271774 at 87.4 Myr (reference) | `08ad9461d9585e92be93a58d81e326c7e4b685ae7d462135d01b0e9931ed762b` |
| seed 447276134 at 116.3 Myr (held out) | `1ccae7159f4a6326831ba032078a408f1a692ca4b59ad315e798d7103c71c513` |

Baseline: `origin/main` at `59d2ef3` in a detached worktree. Branch: `claude/issue-522-climate-units`.
The runs go one at a time, with nothing else running.

**Primary measurement: the t = 0 checkpoint.** Both arms compute a fresh climate snapshot of
the same loaded world, so their difference is the climate change alone, with no replay
divergence and no noise. The 2–10 Myr checkpoints are secondary. Small changes there grow
through erosion, ice and topography, so a gap is only suggestive unless it grows with time.

**Expected at t = 0** (from the curves; old index × 11.56 vs new g/kg):

| Metric | Expected | Threshold to act on |
|---|---|---|
| Air temperature by 15° band, land and ocean | Within 0.7 °C (interpolation error of the calibration table) | Any band > 1 °C |
| Ocean humidity, 0–15° (SST ~26 °C) | ~5% higher (old 1.32 → 15.3 vs cap 16.2) | — |
| Ocean humidity, 45–75° (SST 5–12 °C) | 10–25% higher (old floor 0.3 removed, Clausius–Clapeyron steeper below 20 °C) | — |
| Ocean humidity, 75–90° (SST ≲ 2 °C) | ~5–10% lower | — |
| Mean land precipitation | Within ±5% | > 10% |
| Land precipitation by band | Within ±15%; the 45–75° bands wetter | Any band > 25% |
| Land share of precipitation | Within ±1 percentage point | > 2 pp |
| Köppen group shares of land | Within ±3 pp | > 5 pp |

**Secondary at 10 Myr.** Ice cap (% of land) and land fraction: flag a change of more than
2 pp of land. The direction follows the 45–90° precipitation change and isn't predicted.

If a primary metric crosses its threshold, retune the precipitation coefficient or the marine
relative humidity before opening the PR, or report why not.
