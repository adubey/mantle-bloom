# Issue #942: what sets plate speed

Status: in progress. First pass at deliverable 2 (the torque breakdown, standalone speeds and a brake-off run) on two worlds. The speed census on the 20
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
  couple through it. The brake-off runs also have their own history (more splits). The next section measures
  which of these holds the remaining factor of about 1.5–2: neither. It is the slab terms.

**Limits.** Snapshot estimates on the brake-on history, each plate solved alone. They leave out what faster plates
would do to the gate, the boundaries and plate shapes.

## Brake off

**Finding.** With the contact brake off, plates reach 2.3 cm/yr, and what holds them below the 3.6 cm/yr that mantle
flow would give is sunk-slab drag, not the blocked-edge contact force. Once plates move fast enough for boundaries
to pass the 1 cm/yr gate, subduction switches on, and its drag outweighs its pull: slabs are a net brake. A retune of
the contact brake alone therefore tops out near 2.3 cm/yr on these two worlds, about half of Earth's mean.

**Run.** The same two worlds and schedule with `contact_drag_multiplier = 0` from the start (the #888 brake-off
knob, on the production splitter), pooled into [`torque-breakdown-brake-off.json`](torque-breakdown-brake-off.json).
Predictions 9 and 10 were posted before the run.

| | 106533960 | 1097624328 | Pooled | Pooled, brake on |
|---|---:|---:|---:|---:|
| Solved crust speed (cm/yr) | 2.75 | 1.87 | 2.31 | 1.05 |
| continental / oceanic | 2.50 / 2.98 | 1.56 / 2.14 | 2.03 / 2.56 | 0.98 / 1.12 |
| Own loads, no contact (cm/yr) | 2.91 | 2.01 | 2.46 | 3.41 |
| Mantle flow only (cm/yr) | 4.85 | 2.33 | 3.59 | 3.67 |
| **Driving:** mantle flow | 64% | 44% | 56% | 80% |
| slab pull | 35% | 53% | 42% | 15% |
| ridge push | 0.5% | 0.7% | 0.6% | 0.1% |
| **Resisting:** sunk-slab drag | 51% | 48% | 50% | 7% |
| basal drag | 43% | 45% | 44% | 26% |
| blocked edges (residual) | 3% | 5% | 4% | 6% |
| contact brake | 0 | 0 | 0 | 59% |
| Closing down-going boundary length past the 1 cm/yr gate | 66% | 66% | 66% | 11% |

| # | Prediction | Result |
|---|---|---|
| 9 | Blocked edges are over 40% of the resisting torque and solved speed is under 70% of own-loads speed (against: under 15%, or over 90%) | Counts against: 4%, and solved speed is 94% of own-loads speed |
| 10 | Own-loads-no-contact speed stays within 25% of its brake-on-history value of 3.41 cm/yr (against: under 2.4) | Not met, by a little: 2.46, which is 28% lower |

**What this shows.**

- The three steps from the prescribed flow to the solved speed, on these two worlds: mantle flow alone gives about
  3.6 cm/yr; the slab terms take it to about 2.5; the contact brake takes it to about 1.05.
- With the brake on, the slab terms are nearly absent, because few boundaries close faster than 1 cm/yr. So the
  brake-on estimate of own-loads speed (3.41 cm/yr) was too high: it was taken with subduction mostly switched off.
  With the brake off, two thirds of the closing length passes the gate.
- Slab pull and sunk-slab drag are both larger with the brake off, and the drag wins: own-loads speed (2.46) is
  below mantle-flow-only speed (3.59). Sunk-slab drag resists a subducting node's motion in every direction, at
  `SLAB_MANTLE_VISCOSITY_PA_S` = 10²¹ Pa·s, while slab pull acts only toward the trench and pulls from trenches on
  different sides of a plate partly cancel. Whether this is the cause of the net drag has not been tested.
- This reverses the earlier reading of the gate. Opening it does not free plates; in the current force model it
  adds a net resistance.
- On Earth, plates with attached slabs are the fastest. Here the slab slows its plate, and speed follows the
  prescribed flow.

**Validity.** The solver converged at every sample, and the contact-brake term is exactly zero on every plate.
Solved speeds are near the #888 brake-off arm (2.06 / 2.91 and 1.29 / 1.58 cm/yr), somewhat higher.

**Limits.** Two worlds, plates-only stepping. No collision measure was taken, so this says nothing about what
loosening the brake costs (deliverable 3).

**Reproduce.**

```bash
backend/.venv/bin/python bin/debug/plate_torque_breakdown.py --out <dir> --seeds 106533960 1097624328
backend/.venv/bin/python bin/debug/summarize_torque_breakdown.py <dir> --out analysis/issue942/torque-breakdown.json
# brake off: add --world-knobs '{"contact_drag_multiplier": 0.0}' to the first command
```

Raw outputs of the rerun (SHA-256): `106533960.json` `cc41c208251ae77e3be2fcce809dc4ed7a7632a38d580200864f1942f99612ab`,
`1097624328.json` `c576fa804c45a7864a5a172c00ef17b47e4d7a2671ac9f863a70659cce1b331d`. Brake off: `106533960.json`
`463002e3551357b02e8d8a5591b1aa77e0e3c8e07487f034ee1816f2b1f91d98`, `1097624328.json`
`fd2f3a4477c321e7bf676c9010f37491a9e9bb47d0e25f60dde0a79d24850a58`.
