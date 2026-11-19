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
| Closing down-going length that pulls | 69% |
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

## Results, brake off (five stage 1 worlds, scored 2026-11-19)

| Measure | Baseline | B.1 | B.2 | B.3 | B.4a | B.4b |
|---|---:|---:|---:|---:|---:|---:|
| T1 NNR speed, cm/yr | 1.80 | 1.95 | 2.25 | 2.68 | 1.66 | 1.43 |
| T2 ratio (+10 Myr) | 1.76 | 1.52 | 1.97 | 1.62 | 1.61 | 1.80 |
| T3 pull share of drive | 54.5% | 51.8% | 52.6% | 77.7% | 58.6% | 49.0% |
| T4 MPa | 0.80 | 1.28 | 0.74 | 0.53 | 0.33 | not scored |
| T5 N/m | 4.3e10 | 4.0e11 | 4.7e10 | 5.9e10 | 4.3e10 | 3.8e10 |
| Net slab work | +0.5% | -5.8% | +16.7% | +17.7% | -5.2% | -15.3% |
| Pull per metre, N/m | 2.98e12 | 2.93e12 | 2.79e12 | 4.89e12 | 2.44e12 | 2.43e12 |
| Closing length pulled | 69% | 71% | 73% | 80% | 83% | 73% of total |
| Splits per world per 100 Myr | 7.8 | 8.1 | 8.1 | 29.4 | 11.3 | 14.8 |
| Plates (mean) | 19.0 | 18.6 | 19.6 | 35.2 | 22.2 | 25.3 |

Against each adoption test: B.1 T5 4.0e11 is 19% under 5e11, so it passes (T2 and T4 move the wrong way);
B.2 net slab work +16.7% misses +20% (predicted +20 to +40%); B.3 pull per metre +64% passes, but the split rate
rises 3.8x and the plate count to 35; B.4a pulls 83% of the closing length, fails as predicted; B.4b pulls 73% of
the total (predicted 55 to 70%), fails as predicted.

Brake-on runs are still running.

## Decision (owner, 2026-11-19)

B.2 is adopted although its net slab work (+16.7%) missed the +20% bar: it raises T1 (2.25 cm/yr) and T2 (1.97)
without raising the split rate, and the bar was a prediction, not a derived limit. Adopted so far: B.1 and B.3
(pass), B.2 (owner). B.4a and B.4b failed and go to G2 only through the gate-speed search.
