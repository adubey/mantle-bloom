# Issue #942: what sets plate speed

Status: in progress. First pass at deliverable 2 (the torque breakdown, standalone speeds, a brake-off run and the slab terms' work budget) on two worlds. The speed census on the 20
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

## Why slabs slow their plates

**Finding.** With the brake off, sunk-slab drag takes more work out of the plates than slab pull puts in, for three
reasons of similar size: 40% of the drag is on trench-parallel motion, which pull cannot offset; 34% of the
subducting length is dragged but not pulled, because its column is buoyant; and the drag on trench-normal motion
alone is already two thirds of the pull's work. Without the trench-parallel part the slab terms would net positive.

**Run.** The brake-off run repeated with a work budget over each plate's subducting boundary (its driving and
resisting shares are identical to the earlier run). Work is force times the node's own velocity at the solved
rotation, given as a share of the work slab pull does where plates move toward their trench. Predictions 11–14 were
posted before the run.

| Pooled over plate-samples | 106533960 | 1097624328 | Pooled |
|---|---:|---:|---:|
| Pull work where the plate moves toward its trench | +100% | +100% | +100% |
| Pull work where the plate moves away from its trench | −11% | −5% | −8% |
| Drag on trench-normal motion | −77% | −55% | −68% |
| Drag on trench-parallel motion | −53% | −37% | −46% |
| **Net slab work** | **−41%** | **+4%** | **−22%** |
| Subducting length dragged but not pulled (buoyant column) | 35% | 34% | 34% |
| Subducting length where the plate moves away from its trench | 26% | 17% | 23% |
| Subducting length past the pull-to-drag terminal speed | 25% | 25% | 25% |
| Pull torque surviving cancellation across the plate | 62% | 64% | 63% |
| Mean slab pull (10¹² N/m; reference 4.9) | 2.84 | 2.95 | 2.89 |
| Mantle lithosphere thickness at the trench (km; reference 40) | 50.8 | 52.8 | 51.7 |

| # | Prediction | Result |
|---|---|---|
| 11 | Trench-parallel motion takes at least a third of the drag dissipation | Holds: 40% |
| 12 | Less than 60% of the pull torque survives cancellation | Not met, by a little: 63% |
| 13 | Mean pull at the trench is below the 4.9×10¹² N/m reference, because trench lithosphere is thin or buoyant | Holds for the number (59% of the reference) and for buoyancy; thinness is wrong, the lithosphere there is 52 km |
| 14 | Slab work nets negative | Holds pooled (−22%); one world is at +4% |

**What this shows.**

- Sunk-slab drag is isotropic: it resists a subducting node's whole velocity. Slab pull acts only toward the
  trench. With drag on trench-normal motion only, the pooled net would be about +24% in place of −22%, at these
  speeds.
- A third of the subducting length gets drag and no pull. `subducting_boundary_mask` includes continental nodes
  behind a slab (`slab_attached_continent`), whose buoyant columns pull nothing (`slab_pull_force_n` clips at zero),
  and the drag has no such clip. How much of the drag dissipation falls on that length was not measured.
- The pull that is applied is partly lost: 37% cancels across the plate, and 8% of its work is negative where the
  neighbour overrides a plate moving the other way.
- A quarter of the subducting length is past its terminal speed, where drag exceeds pull even on normal motion.
  At `SLAB_MANTLE_VISCOSITY_PA_S` = 10²¹ Pa·s over a 100 km shear zone and a 250 km slab, that speed is 6.2 cm/yr for
  the 40 km reference column and lower for buoyant ones.

**Limits.** Two worlds, plates-only stepping, brake off. A work budget at the solved speed; it does not give the
speed a changed slab law would produce.

## The mantle flow reference rate

A check of `MANTLE_FLOW_REFERENCE_RATE` (4 cm/yr) against the literature, with no run.

- **Origin.** Added on 2026-08-15 (`ddb78dc2`) with the first motion model, which fitted each plate's rotation to
  the flow, so it was in effect the target plate speed. It predates the torque engine, sunk-slab drag and the
  contact brake, and has not been recalibrated.
- **Magnitude: plausible.** Centres of 2–6 cm/yr give a mean flow of 4.2 cm/yr under plates. With plates about
  3 cm/yr slower, the basal coefficient of 5×10¹⁴ Pa·s/m gives a traction near 0.5 MPa (an estimate from the speed
  difference), inside published basal tractions of about 0.1–1 MPa (Bird, Liu & Rucker 2008 report no more than
  1 MPa under slabless plates). The viscosity of 5×10¹⁹ Pa·s is at the low end of what global models need.
- **Role: not Earth-like.** The flow supplies 80% of the driving torque with the brake on and 56% with it off. On
  Earth slab pull and slab suction are the main drivers (about 60% and 40% in Conrad & Lithgow-Bertelloni 2004),
  and basal traction mostly resists. The model's centres are random and fixed, unrelated to where the slabs are;
  on Earth the flow under plates is largely set up by the slabs.
- **Consequences.** A plate dragged by the flow cannot outrun it, so the constant is a ceiling on plate speed. The
  random draw also sets each world's speed: mantle-flow-only speed is 4.8 cm/yr in one world and 2.5 in the other.
- **Retune.** Doubling the rate would probably reach Earth-like speeds, with driving tractions over 1 MPa and
  plates tied more closely to a fixed pattern. Slab-attached plates would still not be the fast ones.

## Diagnosis so far

On two worlds: plate speed is the prescribed mantle flow (about 3.6 cm/yr), less a net drag from the slab terms
(to about 2.5), less the contact brake (to about 1.05). Ridge push is under 1% of the drive in every run. The
brake alone cannot bring speeds to Earth's: with it off they are 2.3 cm/yr here and 1.55 / 2.12 cm/yr continental /
oceanic over five worlds in #888. Earth-like speeds, and slab-attached plates being the fast ones, need a change to
the driving forces, which is RFC-sized under the issue's own rule.

**Reproduce.**

```bash
backend/.venv/bin/python bin/debug/plate_torque_breakdown.py --out <dir> --seeds 106533960 1097624328
backend/.venv/bin/python bin/debug/summarize_torque_breakdown.py <dir> --out analysis/issue942/torque-breakdown.json
# brake off: add --world-knobs '{"contact_drag_multiplier": 0.0}' to the first command
```

Raw outputs of the rerun (SHA-256): `106533960.json` `cc41c208251ae77e3be2fcce809dc4ed7a7632a38d580200864f1942f99612ab`,
`1097624328.json` `c576fa804c45a7864a5a172c00ef17b47e4d7a2671ac9f863a70659cce1b331d`. Brake off: `106533960.json`
`463002e3551357b02e8d8a5591b1aa77e0e3c8e07487f034ee1816f2b1f91d98`, `1097624328.json`
`fd2f3a4477c321e7bf676c9010f37491a9e9bb47d0e25f60dde0a79d24850a58`. Brake off with the slab budget:
`106533960.json` `0d53b9291641061a14f2beaebfbd9ccd27e86c74f0f040511a9786950a2c46b3`, `1097624328.json`
`aeb9af59573b5ddbbbddfcd59e5cd3c02eb00f1dc061e07cc450ddaef8e0ae7f`.
