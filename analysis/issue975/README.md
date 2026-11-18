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

## Results

Running.
