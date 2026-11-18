# Issue 958: RFC 0021 G0, targets and harness

G0 of [RFC 0021](../../rfc/0021-plate-force-calibration.md#diagnosis-gates) fixes the calibration targets and builds
the scorer before any force correction is probed. It changes no model code. Tracking issue #957.

## Files

| File | What |
|---|---|
| `development-seeds.json` | RFC 0018's 20 development worlds (#888), as #942's census ran them, and the five stage 1 worlds |
| `heldout-seeds.json` | RFC 0021's 20 held-out seeds, frozen by `bin/debug/freeze_rfc0021_heldout.py` |
| `settings.json` | the G0 settings: the rifting baseline with the brake on and off, and the great-circle splitter |
| `earth-plate-speed.json` | T1's Earth figure from PB2002, `bin/debug/earth_plate_speed.py` |
| `t2-census.json` | per plate of each generated development world, B.4b's slab-attached length and its band |
| `t2-threshold.json` | T2's slab-attached share, fixed by the owner before G1 |
| `g0-baseline-brake-off.json` | the scored 20-world brake-off baseline |
| `g0-step0.5-brake-off.json` | the scored 0.5 Myr step check on three stage 1 worlds |

## Rifting baseline

#955's `deep-marks` arm (`analysis/issue888/rift-lines-baseline/arms.json`): stage 1 setting 0 of #888 with
ocean-floor strength from age, `k_cap` 0.14 and `mark_power` 2.8. It ran on PR #950 at `6bda6906`, with
`run_rift_lines_arms.py` and `score_rift_lines.py` at `bd7b2c3b`, under run labels
`issue888/rift-lines-baseline/deep-marks/<seed>` for the five stage 1 worlds. Brake off, it met 6 of RFC 0018's 10
targets and split 7.8 times per world per 100 Myr. Here it is the `rifting-baseline-brake-off` setting and, with the
brake on, `rifting-baseline`.

## T1: Earth's mean plate speed

Bird's PB2002 model (52 plates, Pacific-fixed Euler poles; Bird 2003, doi:10.1029/2001GC000252) on a 100,000-point
Fibonacci grid, each point assigned to the one plate whose outline winds around it (every point is assigned exactly
once):

| Frame | Area-weighted mean speed (cm/yr) |
|---|---:|
| No net rotation | **3.76** (RMS 4.39) |
| Pacific fixed | 6.06 |

The net rotation removed is 0.649°/Myr, against 0.641°/Myr for NUVEL-1A's Pacific plate in its no-net-rotation
frame, so the frame transformation checks out. The Pacific moves at 6.5 cm/yr and covers 21% of the surface.
Africa, Antarctica, North America and Eurasia (together 43%) move at 1.3–2.8 cm/yr.

3.76 cm/yr is below the RFC's working 4–5 cm/yr. **The owner moved T1 to 3–5 cm/yr in the no-net-rotation frame
(2026-11-18).** The breakdown now removes each sample's net rotation from the model's plates as well
(`speed_nnr_cm_per_yr`, `net_rotation_cm_per_yr`), and the scorer reports T1 both ways.

## T2: what counts as slab-attached

The RFC proposed B.4b's geometric test on the generation snapshots: down-going crust in the boundary band, at a
contact that isn't diverging, with only the older side counting where ocean meets ocean. On the 20 development
worlds' 282 generated plates it doesn't separate anything:

| Share of a plate's band slab-attached under B.4b, at generation | p10 | p25 | p50 | p75 | p90 |
|---|---:|---:|---:|---:|---:|
| All plates | 0.40 | 0.51 | 0.63 | 0.74 | 0.84 |

| Threshold | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 |
|---|---:|---:|---:|---:|---:|
| Share of plate area counted slab-attached | 99% | 98% | 96% | 92% | 82% |

Continental-majority plates (137 of 282) have a median of 69%, higher than the oceanic ones. Two things cause it:

- **Generation leaves every plate at zero rotation.** The world's motion solve at generation
  (`world.py`, for `assign_initial_rock_ages`) dates the ocean floor and is then discarded, so "not diverging"
  passes at every contact. Any-rate closing and the 1 cm/yr gate are both zero everywhere.
- **Ocean floor ages don't step at contacts.** Across ocean–ocean contacts that aren't diverging, the two sides'
  ages differ by a median 3.3 Myr (p90 26 Myr), at a median 47 Myr old, on world 106533960. So the older-side rule
  is close to a coin toss.

Stepping doesn't rescue it. On world 106533960 after 5 and 10 Myr, the test marks a median 49% of each band, and
any-rate closing 40–43%. Their rank correlations with plate speed are −0.10 and −0.25 (B.4b) and 0.22 and 0.28
(closing). About half of every plate's boundary faces converging ocean floor.

**The owner chose (2026-11-18): a plate is slab-attached when it has a slab of its own.** That means oceanic floor
it lost under its neighbours in the last 10 Myr (`collision_polarity`'s consumption evidence, `ROLE_LOWER`) is at
least a set share of its area. The breakdown records it per plate as `g0.consumed_10myr_km2`.

On the 20-world baseline below (7,634 plate samples from 10 Myr on), 99% of plates lost some floor in the window, the
median share is 15% of the plate's area, and the share's rank correlation with no-net-rotation speed is 0.51:

| Share | 5% | 10% | 15% | 20% | 30% | 50% |
|---|---:|---:|---:|---:|---:|---:|
| Plate area counted slab-attached | 68% | 32% | 18% | 11% | 5% | 2% |
| Speed ratio, slab-attached to the rest | 1.52 | 1.58 | 1.59 | 1.61 | 1.66 | 1.70 |

**The owner fixed the share at 10% (2026-11-18)**, which puts about Earth's share of the surface on slab-attached
plates. The ratio is under 2.5 at every share, so the choice does not decide the baseline's verdict. (The scorer's
1.54 differs from this table's 1.58 because it also counts the 0 Myr sample.) The owner also asked to think about seeding this at
generation with a pseudo-history, so that worlds start with slabs (see tracking issue #957).

## Held-out set

`heldout-seeds.json`: SHA-256 of `mantle-bloom/rfc0021/v1/heldout/{index:02d}`, first four bytes big-endian masked
with `0x7fffffff`. Indices 0–19 gave 20 distinct seeds, none zero, none in RFC 0018's preregistration and none a
development world. `freeze_rfc0021_heldout.py --check` re-derives the set and compares it with the file.
`calibrate_plate_forces.py` refuses these seeds, and RFC 0018's held-out and reserve seeds, unless given `--heldout`
(G3).

## Harness

`bin/debug/calibrate_plate_forces.py` runs one setting on one world (`run`), a seed list in series (`queue`), and
pools a setting's worlds into one JSON (`score`):

- **forces:** `plate_torque_breakdown.sample` every 10 Myr, which now also records the G0 measures:
  - T2: consumed floor, and B.4b's length for comparison;
  - T4: area-weighted basal traction `μ/d_s |v_mantle − v_plate|`;
  - T5: ridge push over its shared front length at the pushing nodes, and the plate's oldest oceanic floor;
  - no-net-rotation speeds;
- **RFC 0008:** `probe_craton_contacts.py`'s per-step contact rows, from which the score takes collision-time drops
  on jammed pairs, the craton overlap at fusion, the mean largest-plate share and the plate count. The probe's tap on
  `collision_drag_coefficients` predated #760's brake gains and failed on them; it now passes them through, and its
  brake figure includes them. Every step also logs each plate's `|ω| R` and area, for T6's pre- and post-contact
  speeds;
- **RFC 0018:** `run_rift_lines_d1.py`'s split record and recorder, scored by `rift_line_recorder.pool` against the
  production splitter's rate on the same worlds;
- **closure:** the solver's residual over its driving torque, its maximum over every sample.

A 10 Myr smoke run of world 106533960, brake off, gave:
- 2.9 cm/yr, with mantle flow at 76% of the driving torque;
- a solver closure of 1.8×10⁻¹⁵;
- 18.6 plates on average.

## Runs

In series, plates-only, 1 Myr steps, 200 Myr, sampled every 10 Myr, under `~/mantle-bloom-captures/issue958/`:

1. `repro942/great-circle-brake-off/`: #942's two brake-off worlds, to check the harness against #942's numbers.
2. `g0/rifting-baseline-brake-off/`: the 20 development worlds.
3. `g0-step0.5/rifting-baseline-brake-off/`: three stage 1 worlds at 0.5 Myr.

All three finished without failures; every solver sample converged, with closure at most 6.5×10⁻¹⁵.

**Reproduction.** On #942's two brake-off worlds the scorer matches `analysis/issue942/torque-breakdown-brake-off.json`
to four decimals: continental 2.0335 and oceanic 2.5569 cm/yr, mantle flow 56.1% and slab pull 41.7% of the
driving torque, sunk-slab drag 49.9% and basal drag 44.0% of the resisting.

**20-world baseline** (`g0-baseline-brake-off.json`):

| Target | Baseline | Band | |
|---|---:|---|---|
| T1 mean speed, no net rotation | 1.88 cm/yr (1.90 with it) | 3–5 | fails |
| T2 speed ratio, slab-attached (32% of area) to the rest | 1.54 | ≥ 2.5 | fails |
| T3 slab pull share of driving torque | 55% | ≥ 50% | passes |
| T4 basal traction under plates without slabs | 0.22 MPa | ≤ 1 | passes |
| T5 ridge push, oceanic floor ≥ 50 Myr | 4.5×10¹⁰ N/m | 1–4×10¹² | fails |
| T7 fastest to slowest world | 2.0 (1.32–2.64 cm/yr) | ≤ 3 | passes |
| T8 contacts | 122 fusions, 1 collision-time drop on a jammed pair, largest plate 23%, 18.7 plates | not worse than `main` | for G2 |
| T9 RFC 0018 targets met | 5; 9.0 splits per world per 100 Myr | | split rate needs the great-circle rate |

Continental plates move at 1.62 cm/yr and oceanic ones at 2.17. Sunk-slab drag (51%) and basal drag (45%) are the
resisting torque, and ridge push is 0.6% of the driving torque. 30% of the subducting length carries no pull. T6's
pre- and post-contact speeds are in each run's `breakdown.json`, not yet scored.

**0.5 Myr step check** (`g0-step0.5-brake-off.json`, against the same three worlds at 1 Myr):

| | 1 Myr | 0.5 Myr | Change |
|---|---:|---:|---:|
| T1 (cm/yr, no net rotation) | 1.89 | 1.98 | +5.1% |
| T2 ratio | 1.88 | 1.79 | −4.7% |
| T3 | 48% | 49% | +2.0% |
| T4 (MPa) | 0.91 | 0.47 | −49% |
| T5 (10¹⁰ N/m) | 3.76 | 3.89 | +3.5% |
| T7 | 1.44 | 1.49 | +3.9% |
| Splits (106533960, 1097624328, 1397817354) | 23, 15, 10 | 38, 16, 19 | +52% |
| Mean plate count | 17.5 | 19.2 | +9.8% |

T1 and T2 are inside the step check's 15% and 20%. The rift-line splitter is not: it splits half as often again at
half the step, so its rate is a property of the 1 Myr step, and every RFC 0021 and RFC 0018 split-rate comparison has
to use the same step.

## Reproduce

```bash
python3 bin/debug/freeze_rfc0021_heldout.py --check
python3 -I bin/debug/earth_plate_speed.py --plates <dir>/PB2002_plates.dig.txt --poles <dir>/PB2002_poles.dat.txt \
    --points 100000 --out analysis/issue958/earth-plate-speed.json
backend/.venv/bin/python bin/debug/calibrate_plate_forces.py t2-census --seeds-file analysis/issue958/development-seeds.json \
    --out analysis/issue958/t2-census.json
backend/.venv/bin/python bin/debug/calibrate_plate_forces.py queue --settings analysis/issue958/settings.json \
    --setting rifting-baseline-brake-off --seeds-file analysis/issue958/development-seeds.json --out ~/mantle-bloom-captures/issue958/g0
backend/.venv/bin/python bin/debug/calibrate_plate_forces.py score --runs ~/mantle-bloom-captures/issue958/g0/rifting-baseline-brake-off \
    --production-rate <great-circle rate on the same worlds> --t2-share <fixed share>
```
