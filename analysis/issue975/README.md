# Issue 975: RFC 0021 G1, each force correction alone

G1 of [RFC 0021](../../rfc/0021-plate-force-calibration.md#diagnosis-gates) (tracking #957). It uses the switches from #971
(PR #974) and the T2 scorer from PR #970. Every run is plates-only, with rift lines at the rifting baseline (#955's
`deep-marks`), 200 Myr, 1 Myr steps, sampled every 10 Myr, on the five stage 1 worlds (`stage1-seeds.json`).
`settings.json` holds each correction alone, brake on and off.

## Baseline (brake off, the five stage 1 worlds)

From G0's runs of the same worlds (`analysis/issue958/`). With every switch off the solve is bit-identical with and
without #971 (`test_default_forms_leave_the_motion_solve_bit_identical`), so they stand for the baseline here.

| Measure | Baseline |
|---|---:|
| T1 mean speed, no net rotation | 1.80 cm/yr |
| T2 speed ratio, slab-attached to the rest | 1.76 |
| T3 slab pull share of driving torque | 54.5% |
| T4 basal traction under plates without slabs | 0.80 MPa |
| T5 ridge push, ocean floor ≥ 50 Myr | 4.3×10¹⁰ N/m |
| Net slab work, share of pull's power | +0.5% (drag trench-normal −57%, trench-parallel −38%) |
| Slab pull per metre of subducting length | 3.0×10¹² N/m, on a 49 km column |
| Closing down-going length that pulls | 47% (69% counting nodes that pull nothing, G0's measure) |
| Ridge push share of driving torque | 0.6% |
| Sunk-slab drag share of resisting torque | 52% |
| Splits per world per 100 Myr | 7.8 |

## Predictions (posted on #975 before any run, 2026-11-18)

Brake off, against the baseline above. The adoption tests are the owner's sharpened ones (#975). With the brake on,
the predictions are the same directions, with smaller speed changes.

| Run | Its own test | Predicted | Other predictions |
|---|---|---|---|
| B.1 | T5 within 20% of the column calculation (about 5×10¹¹) | T5 4–5.3×10¹¹: **passes** | ridge push 5–15% of drive; T1 2.0–2.1; T2 1.6–1.8, since it pushes every spreading plate; splits within ±20% |
| B.2 | net slab work ≥ +20% | +20% to +40%: **passes** | sunk-slab drag 35–45% of resistance; T1 2.1–2.6; T2 1.9–2.3; splits up 0–30% |
| B.3 | pull per metre up ≥ 30% | +30% to +70%: **passes, near the bar** | T3 60–65%; T1 2.0–2.3; T2 1.8–2.0; T5 unchanged |
| B.4a | closing length pulling ≥ 90% | 82–88%: **fails** | about 31% of closing length moves under 1 cm/yr and gets about half weight at `v_gate` 1 cm/yr; T1 within +10% |
| B.4b | closing length pulling ≥ 90% | 55–70% of the total: **fails**; 100% of the down-going side by construction | it drops the younger ocean side; T1 within ±10% |

If B.4a fails as predicted, the reason is `v_gate`, which G2 searches over 0.25–2 cm/yr. B.4a would then be carried
into G2 only if it raises T1 or T2.

## Review fixes (PR #979, 2026-11-19)

Codex's review found that B.2's drag acted on every subducting node, where RFC 0021 sums it over nodes with positive
pull, and that the "closing length pulled" measure counted nodes that pull nothing. Both are fixed in #974
(`slab_drag_nodes`, `pulling_nodes`). Muse and Codex also found that the force-constant overrides ran before world
generation (so a corner changed the starting ocean-floor ages; now applied after it), that a resumed run was matched on
its directory name only (now on step length, cadence and the resolved setting), and that knob values were unchecked.
Every brake-off setting and the four brake-on settings that contain B.2 were re-run, 200 Myr, 1 Myr steps, five stage 1
worlds, into a new capture directory; the old runs are superseded. The brake-on runs without B.2 (`rifting-baseline`,
`b1`, `b3`, `b4a`, `b4b`) were not re-run: their dynamics don't depend on the fix, but their "closing length pulled"
still uses the old count (the first table of the old README), so read that row for them as an upper bound.

The baseline is now the stage 1 worlds' own brake-off run (`rifting-baseline-brake-off.json`), not G0's 20-world file.

## Results, brake off (five stage 1 worlds)

| Measure | Baseline | B.1 | B.2 | B.3 | B.4a | B.4b |
|---|---:|---:|---:|---:|---:|---:|
| T1 NNR speed, cm/yr | 1.80 | 1.95 | 2.53 | 2.68 | 1.66 | 1.43 |
| T2 ratio (+10 Myr) | 1.76 | 1.52 | 1.79 | 1.62 | 1.61 | 1.80 |
| T3 pull share of drive | 54.5% | 51.8% | 56.5% | 77.7% | 58.6% | 49.0% |
| T4 MPa (see note) | 0.80 | 1.28 | 0.81 | 0.53 | 0.33 | - |
| T5 N/m | 4.33e+10 | 4.03e+11 | 5.56e+10 | 5.86e+10 | 4.26e+10 | 3.83e+10 |
| T7 fastest/slowest world | 1.44 | 1.52 | 1.62 | 1.29 | 1.62 | 1.68 |
| Net slab work | +0.5% | -5.8% | +24.0% | +17.7% | -5.2% | -15.3% |
| Pull per metre, N/m | 2.98e+12 | 2.93e+12 | 2.80e+12 | 4.89e+12 | 2.44e+12 | 2.43e+12 |
| Closing length pulled | 47% | 48% | 53% | 54% | 56% | 42% |
| Ridge push share of drive | 0.6% | 6.3% | 0.8% | 0.6% | 0.6% | 0.5% |
| Splits per world per 100 Myr | 7.8 | 8.1 | 12.6 | 29.4 | 11.3 | 14.8 |
| Jammed collision-time drops (T8) | 0 | 0 | 1 | 1 | 0 | 0 |
| Plates (mean) | 19.0 | 18.6 | 21.5 | 35.2 | 22.2 | 25.3 |

Against each adoption test:

| Correction | Test | Result |
|---|---|---|
| B.1 | T5 within 20% of the column calculation, taken as the rounded 5×10¹¹ N/m (owner, 2026-11-19; the exact column value is 5.31×10¹¹, which 4.03×10¹¹ misses by 24%) | 4.03×10¹¹, 19% under: **pass** (T2 falls to 1.52 and T4 rises to 1.28) |
| B.2 | net slab work at least +20% | +24.0%: **pass** (T1 2.53; splits 12.6 against 7.8) |
| B.3 | pull per metre of subducting length up at least 30% | 4.89×10¹² against 2.98×10¹², +64%: **pass**, with 29.4 splits and 35 plates |
| B.4a | closing length pulling at least 90% | 56%: **fail** |
| B.4b | the same | 42%: **fail** |

Before the fix B.2 read +16.7% and was adopted by the owner as an exception; it now passes outright.

**T4 note.** T4 is the mean basal traction under plates with no subducting node, and B.4 changes which nodes count as
subducting, so each B.4 setting averages over a different set of plates. `b4b-brake-off`, `b4a` and `b4b` (brake on),
and `b123-b4a` have no such plate and score nothing ("-"); the B.4 settings' T4 values are not comparable with each
other or with the baseline's.

## Results, combined B.1 + B.2 + B.3

Brake off:

| Measure | b123 | + B.4a | + B.4b |
|---|---:|---:|---:|
| T1 NNR speed, cm/yr | 4.25 | 4.51 | 4.13 |
| T2 ratio (+10 Myr) | 1.59 | 1.46 | 1.61 |
| T3 pull share of drive | 60.0% | 61.2% | 54.4% |
| T4 MPa (see note) | 1.52 | 1.78 | 1.53 |
| T5 N/m | 1.90e+12 | 1.93e+12 | 1.96e+12 |
| T7 fastest/slowest world | 1.54 | 1.51 | 1.53 |
| Net slab work | +29.0% | +26.0% | +28.8% |
| Pull per metre, N/m | 4.38e+12 | 4.04e+12 | 4.00e+12 |
| Closing length pulled | 56% | 60% | 43% |
| Ridge push share of drive | 19.9% | 20.4% | 23.1% |
| Splits per world per 100 Myr | 11.8 | 14.8 | 13.9 |
| Jammed collision-time drops (T8) | 3 | 6 | 3 |
| Plates (mean) | 17.5 | 20.9 | 20.2 |

Brake on:

| Measure | B.2 | b123 | + B.4a | + B.4b |
|---|---:|---:|---:|---:|
| T1 NNR speed, cm/yr | 0.66 | 1.02 | 1.34 | 1.28 |
| T2 ratio (+10 Myr) | 2.45 | 2.62 | 2.14 | 1.65 |
| T3 pull share of drive | 14.7% | 51.4% | 65.5% | 63.7% |
| T4 MPa (see note) | 0.41 | 0.42 | - | 0.37 |
| T5 N/m | 1.35e+10 | 1.44e+12 | 1.47e+12 | 1.50e+12 |
| T7 fastest/slowest world | 3.24 | 2.76 | 1.94 | 1.35 |
| Net slab work | +65.0% | +71.6% | +68.3% | +51.7% |
| Pull per metre, N/m | 3.85e+12 | 7.36e+12 | 4.21e+12 | 4.77e+12 |
| Closing length pulled | 8% | 24% | 50% | 42% |
| Ridge push share of drive | 0.1% | 9.5% | 11.2% | 10.5% |
| Splits per world per 100 Myr | 3.6 | 8.9 | 17.0 | 25.3 |
| Jammed collision-time drops (T8) | 0 | 0 | 1 | 1 |
| Plates (mean) | 17.5 | 21.7 | 29.4 | 40.0 |

Brake off, the combined corrections reach T1 (4.25 cm/yr), T3 (60%) and T5 (1.90×10¹²) at the default constants, but T2
is 1.59 (the baseline's is 1.76; target at least 2.5) and T4 is 1.52 MPa (limit 1), and three collision-time drops occur
while jammed (the baseline has none). Brake on, T1 is 1.02 and T2 2.62 (the brake slows the plates without slabs, so T2
there says little about the slab forces). Solver closure is at most 2.0×10⁻¹⁴ of the driving torque in every run, read from the read-only solves at the 10 Myr samples.

## Second review round (2026-11-20)

Fixed in the harness and #974 after Codex's and Muse's second reviews; none of them moves a G1 result, so nothing was
re-run:

- The resume check also compares the run length and whether contacts were captured (a finished 10 Myr run, or one
  made with `--no-contacts`, is no longer reused for a 200 Myr run with contacts).
- B.2's drag set is read from the actual pull (`pulling_nodes`: buoyant, nonzero `f_pull`, nonzero B.4a ramp weight), so
  G2's `f_pull` and `v_gate` ranges keep drag and pull on the same nodes. At G1's defaults it is the same set as before.
- The breakdown measures B.2's drag work on the same trench normal the solve uses (the earlier figure used the chord
  and understated it slightly), and `down_going_closing_past_gate` is now taken inside the closing set. The checked-in
  scores predate both: B.2's net slab work (+24.0%) is a touch higher than it would read, and brake-on `b4b`'s
  `closing_down_going_length_past_gate` reads 139.7% for the stated reason. The B.4 verdicts use
  `closing_down_going_length_pulled`, which was already inside the closing set.
- RFC B.4b's wording now matches the code (older side where ocean meets ocean, not collision polarity).

## Not done here

- RFC 0021's G1 asks for what alternative D would change in bathymetry and land fraction, and for the ledgers to close.
  These plates-only probes score neither (`run.json`'s `ledger_relative_gap_max` is a placeholder). Both are deferred
  to the close-out issue (#980), which should close the ledgers before any correction becomes a default; G2's
  preregistration should say so.

## Still to run

The best-case check (four `corner-*` settings in `settings.json`: combined corrections, brake off,
`ASTHENOSPHERE_VISCOSITY_PA_S` 2e19, `f_pull` 1.0 (the default), `SLAB_MANTLE_VISCOSITY_PA_S` 1e20 or 1e21,
`SUBDUCTION_LSINK_M` 150 or 660 km) and the 0.5 Myr repeat of the combined setting on three worlds, in #980.
