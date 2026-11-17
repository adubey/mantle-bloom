# Issue #942: what sets plate speed

Status: in progress. First pass at deliverable 2 (the torque breakdown and standalone speeds) on two worlds. The speed census on the 20
development worlds (deliverable 1) and the `contact_drag_multiplier` sweep (deliverable 3) are not done.

## Torque breakdown, two worlds

**Finding.** Plates are carried by the prescribed mantle flow and held back by the contact brake. Slab pull supplies
about 15% of the driving torque and ridge push about 0.1%. The 1 cm/yr motion gate lets slab pull act on about 11%
of the down-going boundary length that is closing.

**Run.** `bin/debug/plate_torque_breakdown.py` at `bf624090`, worlds 106533960 and 1097624328, production
`generate_world(seed)` defaults and the production splitter, plates-only (`simulate_climate_biomes = False`),
200 Myr at 1 Myr steps, brake on. Every 10 Myr the script repeats the motion solve read-only and takes each torque's
component along the plate's solved rotation axis: positive drives the motion, negative resists it. Shares below are
of the summed driving or resisting torque over all plate-samples (1,047), pooled by
`bin/debug/summarize_torque_breakdown.py` into [`torque-breakdown.json`](torque-breakdown.json). The predictions
were posted on the issue before the run.

| | 106533960 | 1097624328 | Pooled |
|---|---:|---:|---:|
| Continental crust speed (cm/yr) | 1.37 | 0.59 | 0.98 |
| Oceanic crust speed (cm/yr) | 1.47 | 0.77 | 1.12 |
| Plate area moving under 1 cm/yr (by `\|ω\| R`) | 16% | 79% | 47% |
| **Driving:** mantle flow | 83% | 73% | 80% |
| slab pull | 12% | 21% | 15% |
| ridge push | 0.1% | 0.1% | 0.1% |
| **Resisting:** contact brake | 59% | 61% | 59% |
| basal drag | 27% | 25% | 26% |
| sunk-slab drag | 6% | 8% | 7% |
| blocked edges (residual) | 6% | 5% | 6% |
| Plate area where the brake is the largest resistance | 66% | 57% | 61% |
| Closing down-going boundary length past the 1 cm/yr gate | 12% | 11% | 11% |
| Opening boundary length past the 1 cm/yr gate | 23% | 7% | 15% |
| Slab pull of all closing length / applied slab pull | 3.1 | 3.5 | 3.3 |

Small shares not listed (a term acting on the other side on a few plates) are in the JSON. "Mantle flow" and
"basal drag" are the two halves of one coupling, `b − K ω`: the asthenosphere pulls a plate toward the local
prescribed flow (`mantle.flow_at`, centres of 2–6 cm/yr). Its net effect is to drive.

**Against the predictions.**

| # | Prediction | Result |
|---|---|---|
| 1 | Ridge push under 1% of driving torque | Holds: 0.1% |
| 2 | Basal drag is the largest resistance on over half the plate area, the brake second | Fails: the brake is largest on 61%, basal drag on 28% |
| 3 | Mantle flow supplies more driving torque than slab pull on plates under 1 cm/yr | Holds: 77% against 9% |
| 4 | On plates under 1 cm/yr, under half the closing down-going length passes the gate and ungated pull is at least twice the applied | Holds: 7% passes; ratio 4.2 |
| 5 | 1097624328 has a smaller gated share than 106533960 | Fails: 11% against 12%, within the 10-point margin |

**What this does and does not show.**

- Slab pull is a minor load in the current model, for three compounding reasons: its magnitude (4.9×10¹² N/m at the
  40 km oceanic reference), the gate (about a ninth of the closing length), and the down-going share of the boundary.
  Without the gate it would be about 3.3 times larger, still under the mantle-flow term.
- The gate applies to fast plates too (14% of closing length passes on plates over 1 cm/yr), so it is not only a
  trap for slow plates. Prediction 5 failing says the gate does not explain why one world is slower than the other.
- Prediction 2 failing changes the reading of the brake: it is the largest resistance, not an addition to basal
  drag. That agrees with the #888 arms, where turning it off raised speed 2–3 times.
- A breakdown gives the balance at the solved speed. It does not give the speed under a changed term; the next
  section estimates that from each snapshot.

**Validity.** The solver converged at every sample. Speeds are near the #888 brake-on arm (1.14/1.35 and 0.41/0.56
cm/yr, rift-line prototype). The closure check (terms summing to zero) could only be read where there are no
blocked edges, since the residual is the blocked-edge term; at 2 Myr with none it was 10⁻¹⁵ of the driving torque.
Both worlds had blocked edges (up to 31 and 22).

**Limits.** Two worlds, plates-only stepping, one sample per 10 Myr.

## Speeds a plate would have alone

**Finding.** Mantle flow alone would carry plates at about 3.7 cm/yr (4.8 and 2.5 in the two worlds), close to
Earth's mean, and they move at about 1.0. Contacts between plates remove about 70% of that speed. The prescribed
flow is not the ceiling that keeps plates slow.

**Run.** The same two worlds and schedule, rerun with three more measures per plate from each snapshot's own terms;
the driving and resisting shares are identical to the first run. Predictions 6–8 were posted before the run.

| Area-weighted crust speed (cm/yr) | 106533960 | 1097624328 | Pooled |
|---|---:|---:|---:|
| As solved (brake on) | 1.42 | 0.68 | 1.05 |
| Mantle flow only (`K_basal ω = b`) | 4.82 | 2.51 | 3.67 |
| Own loads, no contact (adds slab pull, ridge push, sunk-slab drag) | 4.36 | 2.46 | 3.41 |
| Prescribed flow under the plates (mean `\|flow\|`) | 5.24 | 3.08 | 4.16 |
| #888 arm with `contact_drag_multiplier = 0`, continental / oceanic | 2.06 / 2.91 | 1.29 / 1.58 | |

| # | Prediction | Result |
|---|---|---|
| 6 | Mantle-flow-only speed under 2.5 cm/yr (against: 4 or more) | Not met: 3.67 pooled. One world is at the 2.5 line and the other, at 4.8, counts against |
| 7 | Own-loads-no-contact within 30% of mantle-flow-only | Holds: 0.93 times. Slab pull and sunk-slab drag together slow plates slightly |
| 8 | Flow under plates at least 1.5 times the mantle-flow-only plate speed (against: under 1.2) | Counts against: 1.13. A rigid plate keeps nearly all of the flow speed under it |

**What this shows.**

- Plate speed in the current model is the prescribed flow speed (convection centres of 2–6 cm/yr,
  `MANTLE_FLOW_REFERENCE_RATE` = 4 cm/yr) less what contacts take away. The slab terms change it by under 10%.
- Continental crust would be as fast as oceanic crust or faster under its own loads (3.57 against 3.25 cm/yr
  pooled). In the solved motion and in every #888 arm it is slower, so the continental deficit comes from contacts.
- The check against the #888 brake-off arm does not close. Own loads with no contact give 2.5–4.7 cm/yr, and the
  runs with the brake multiplier at zero gave 1.3–2.9. Setting `contact_drag_multiplier = 0` does not remove all
  contact: the blocked-edge contact force (`_blocked_contact_response`) is not scaled by it, and neighbours still
  couple through it. The brake-off runs also have their own history (more splits). Which of these holds the
  remaining factor of about 1.5–2 has not been measured; a breakdown of a brake-off run would show it.

**Limits.** Snapshot estimates on the brake-on history, each plate solved alone. They leave out what faster plates
would do to the gate, the boundaries and plate shapes.

**Reproduce.**

```bash
backend/.venv/bin/python bin/debug/plate_torque_breakdown.py --out <dir> --seeds 106533960 1097624328
backend/.venv/bin/python bin/debug/summarize_torque_breakdown.py <dir> --out analysis/issue942/torque-breakdown.json
```

Raw outputs of the rerun (SHA-256): `106533960.json` `cc41c208251ae77e3be2fcce809dc4ed7a7632a38d580200864f1942f99612ab`,
`1097624328.json` `c576fa804c45a7864a5a172c00ef17b47e4d7a2671ac9f863a70659cce1b331d`.
