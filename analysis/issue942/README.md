# Issue #942: what sets plate speed

Status: in progress. First pass at deliverable 2 (the torque breakdown) on two worlds. The speed census on the 20
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
- A breakdown gives the balance at the solved speed. It does not give the speed under a changed term. In
  particular the speed mantle flow alone would give each plate (`K⁻¹ b`), the ceiling when the flow is the main
  driver, has not been computed.

**Validity.** The solver converged at every sample. Speeds are near the #888 brake-on arm (1.14/1.35 and 0.41/0.56
cm/yr, rift-line prototype). The closure check (terms summing to zero) could only be read where there are no
blocked edges, since the residual is the blocked-edge term; at 2 Myr with none it was 10⁻¹⁵ of the driving torque.
Both worlds had blocked edges (up to 31 and 22).

**Limits.** Two worlds, plates-only stepping, one sample per 10 Myr.

**Reproduce.**

```bash
backend/.venv/bin/python bin/debug/plate_torque_breakdown.py --out <dir> --seeds 106533960 1097624328
backend/.venv/bin/python bin/debug/summarize_torque_breakdown.py <dir> --out analysis/issue942/torque-breakdown.json
```

Raw outputs (SHA-256): `106533960.json` `e12cee4cbf265bc6448418622c4900b87b9a6cae1159880697de503061c7162d`, `1097624328.json` `0b32789ab2368b8b3b3a165b67b6f0172f8652fb59aebdeb2d99780ce0d1fc5e`.
