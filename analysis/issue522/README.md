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

## Results

Four runs, one at a time with nothing else running: `origin/main` at `59d2ef3` and the branch at
`454a1f3`, on each save. The command was the same for both arms, with
`--backend <baseline-worktree>/backend` for `origin/main`:

```bash
backend/.venv/bin/python bin/debug/measure_climate_fields.py <local-save> \
  --myr 10 --step-years 250000 --checkpoint-myr 2 --out <scratch-output>.json
```

The compact numbers are in [`comparison.json`](comparison.json).

### t = 0 (the primary measurement: same world, climate change only)

| Metric | 997271774 main | branch | Δ | 447276134 main | branch | Δ | Threshold |
|---|---|---|---|---|---|---|---|
| Land air temp (°C) | 7.22 | 7.19 | −0.03 | 7.28 | 7.25 | −0.03 | band > 1 °C |
| Worst band, land air temp (°C) | | | −0.28 (75–90°) | | | −0.26 (75–90°) | > 1 °C |
| Ocean humidity (g/kg) | 11.93 | 12.44 | +4.3% | 11.85 | 12.35 | +4.2% | — |
| Land precipitation (mm/yr) | 648 | 659 | +1.7% | 823 | 838 | +1.8% | > 10% |
| Ocean precipitation (mm/yr) | 1471 | 1528 | +3.9% | 1442 | 1502 | +4.1% | — |
| Land share of precipitation (%) | 15.23 | 14.96 | −0.27 pp | 22.39 | 22.00 | −0.39 pp | > 2 pp |
| Worst Köppen group shift (pp of land) | | | −0.52 (B) | | | −0.70 (B) | > 5 pp |

Ocean humidity and land precipitation by latitude band (both hemispheres together; Δ is branch vs main):

| Band | Δ ocean humidity (997 / 447) | Δ land precip (997 / 447) |
|---|---|---|
| 0–15° | +5% / +5% | +3% / +2% |
| 15–30° | +7% / +7% | +2% / +3% |
| 30–45° | +2% / +2% | +3% / +3% |
| 45–60° | −2% / −3% | 0% / 0% |
| 60–75° | +6% / +6% | +2% / +2% |
| 75–90° | +17% / +18% | −5% / −1% |

### 10 Myr later (secondary)

| Metric | 997271774 main | branch | 447276134 main | branch |
|---|---|---|---|---|
| Land (% of surface) | 30.14 | 30.18 | 33.86 | 33.84 |
| Ice cap (% of land) | 5.56 | 5.48 | 5.23 | 5.19 |
| Land precipitation (mm/yr) | 650 | 677 (+4.1%) | 841 | 866 (+3.0%) |
| Land share of precipitation (%) | 15.81 | 15.88 | 23.10 | 22.92 |
| Worst Köppen group shift (pp) | | −1.38 (B) | | −1.01 (B) |
| Land precip, 75–90° (mm/yr) | 96.5 | 83.3 (−14%) | 114.9 | 110.8 (−4%) |

### What it means

- **Temperature didn't move.** The calibration table reproduces the old profile at Earth's tilt:
  every band is within 0.3 °C, slightly cooler poleward of 60°, which is the table's
  interpolation error.
- **Precipitation moved less than expected.** Land precipitation rose about 2% at t = 0,
  inside the ±5% expected. No land band moved more than 5%, well inside the 25% threshold.
  The 45–75° bands didn't get the 10–25% rise predicted: their SSTs (~10–16 °C) sit where the
  old linear index and Clausius–Clapeyron nearly agree. The prediction assumed the old 0.3
  floor bit there, but it binds only below 6 °C.
- **Ocean humidity rose ~4%, so ocean precipitation rose ~4%.** That comes from 15–30° SSTs
  (~24 °C), where Clausius–Clapeyron runs above the old line, and from polar water, where the
  old 0.3 floor was below what Clausius–Clapeyron gives at 4 °C. That rain falls on open
  ocean, outside the land water budget, so it reaches lakes, rivers and erosion only as the
  ~2% more land precipitation above.
- **Downstream barely moved.** After 10 Myr, ice cap and land area differ by under 0.2 pp
  at every checkpoint, with no trend. Köppen shares shift by up to 1.4 pp, mostly fewer arid
  cells. Land precipitation stays 1.6–4.1% higher at every checkpoint on both saves.
- **One persistent regional change: polar land is a little drier.** On save 997271774 the
  75–90° land band is 5% drier at t = 0, and 12–17% drier from the first step on (about
  13 mm/yr on a ~90 mm/yr band); save 447276134 shows −1 to −4%. It doesn't grow, and it's
  well inside the 25% band threshold, but it isn't replay noise. A likely cause, not checked:
  lakes in near-freezing air evaporate less than under the old 0.3 floor (Clausius–Clapeyron
  gives 3.0 g/kg at 0 °C against the floor's 3.47). Ice cap didn't respond.
- No primary metric crossed its threshold, so no retune was made.

**Runtime and memory.** Mean step time was 6.5–7.3 s on both arms, except two windows (main
447276134 at 8–10 Myr, 23 s; branch 447276134 at 0–2 Myr, 15.5 s). Those look like machine
load rather than the change, since the climate code path costs the same in both arms. Peak RSS
was 1.2–1.6 GB.
