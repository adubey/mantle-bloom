# RFC 0008 evidence: craton sinks on current main

Background numbers for [RFC 0008](../../rfc/0008-craton-collision-jamming.md) and the 2026-10-14
revision of [RFC 0002](../../rfc/0002-craton-trench-attrition.md). The gate itself is
below (Phase 1 gate).

`bin/debug/compare_cratons.py --seeds 7,11,2024 --variants on --years 1e9 --every 100e6 --workers 1`
on origin/main at 1b4f686, default density, 5 Myr steps. Numbers are in
[craton_sinks_1gyr.json](craton_sinks_1gyr.json).

This baseline includes #358, which counts the split cooldown in years (2 Myr) instead of steps.
At 5 Myr steps the old cooldown was 20 steps, 100 Myr, so these runs split about 3.5 times as
often as the runs behind #378, #385 and #388 (104–115 splits and 27–55 merges per Gyr, against
32–37 splits and 9–15 merges). Numbers from before #358 at 5 Myr steps aren't comparable.

# Phase 1 gate (issue #395)

The gate in RFC 0008's Design: how often cratons jam plate contacts, whether jams lose
collision time, how much slab pull comes from buoyant columns, and whether craton contacts
close more slowly, before any solution is built.

## Method

`bin/debug/probe_craton_contacts.py --seed <S> --years 1e9`, seeds 7, 11 and 2024, default
density, 5 Myr steps, one at a time, on the gate branch (origin/main 3b9d9c1 plus the probe
and the diagnostic `run_eligible` mask). The probe taps the step in its own process and
feeds nothing back: each seed ends with exactly the craton volume of the baseline above
(679,376,267 km³ on seed 7). Its per-cell loss totals agree with the ledger sinks to within
2%. The per-seed summary is in [gate.json](gate.json).

Definitions are the RFC's (see the probe's docstring):
- a contact node is **jammed** when it is a retreat candidate holding craton (on an active
  polarized front, the frozen `run_eligible` lower-plate cells), or when it and the
  neighbour's nearest node are both cratonic. A pair is jammed with at least 4 jammed nodes
  and closing ≥ −1 cm/yr;
- craton loss is booked to the nearest contact within 500 km (`CRATON_BRAKE_REACH_KM`).
  A fusion's loss goes to the fused pair;
- **J** is the share of all craton destruction near a jammed contact whose pair didn't fuse
  within 50 Myr.

## Results (0–1 Gyr)

| | Seed 7 | Seed 11 | Seed 2024 |
|---|---|---|---|
| Craton destroyed, % of start | 86% | 66% | 68% |
| Near any contact / near a jammed contact, % of destruction | 98% / 93% | 97% / 92% | 97% / 90% |
| **J** | **22%** | **57%** | **40%** |
| Jammed share of continental collision length | 44% | 56% | 37% |
| Collision time dropped: all / jammed pairs | 106 / 28 | 54 / 22 | 111 / 53 |
| Collision time lost on jammed pairs | 455 Myr | 270 Myr | 955 Myr |
| **Buoyant share of slab pull** (cratonic) | **6.4%** (4.4%) | **4.5%** (2.8%) | **4.0%** (3.2%) |
| Closing rate, converging continental contacts: craton in reach / clear (cm/yr) | 1.70 / 1.59 | 1.84 / 1.12 | 1.48 / 1.61 |
| Brake today at jammed contacts, median (N/m) | 7.8 × 10⁴ | 4.9 × 10⁴ | 2.5 × 10⁴ |
| RFC's contact drag at jammed contacts, p10 / p50 / p90 (10¹² N/m) | 0.26 / 0.94 / 2.2 | 0.24 / 0.84 / 1.9 | 0.22 / 0.77 / 1.8 |
| Fusions: roll / forced | 38 / 17 | 15 / 12 | 29 / 7 |
| Fusions of a pair jammed at the time | 53 of 55 | 25 of 27 | 31 of 36 |
| Craton lost in fusions, % of start | 54% | 18% | 31% |
| Craton-on-craton overlap at fusion, summed, % of start | 90% | 75% | 53% |

## Findings

- **Jams are common.** 37–56% of converging continental contact length is craton-jammed,
  far above the falsification threshold (5%). J is 22–57%, so the loss acceptance
  criterion stands.
- **The merge cap is the delamination sink.** Craton lost in fusions (54%, 18% and 31% of
  the starting craton) equals each seed's `delaminated_m3` share in the baseline (54%, 18%,
  31%). Root delamination contributes almost nothing over 1 Gyr. Nearly every fusion is
  of a pair jammed at the time, and the cratons they stack sum to 53–90% of the starting
  craton. RFC 0001's placement fix has to land before or with solution 6, as the RFC says.
- **Jams reset collision time.** 26–48% of collision-time drops are of jammed pairs, which
  together lose 270–955 Myr of accumulated time per seed. Solution 5 targets exactly this.
- **Buoyant columns pull.** Continental-coded or cratonic subducting nodes supply 4.0–6.4%
  of slab-pull torque. That's above the 1% line, so solution 1's resistance is tuned, not
  only its keel fix.
- **Craton contacts don't close more slowly today.** With craton in reach, continental
  contacts close 7% faster, 65% faster and 8% slower than without. That's expected with
  a brake about 10⁸ times too weak (median 2.5–7.8 × 10⁴ N/m against Copley's 5 × 10¹²).
  The RFC's rule for solution 2 compares rates *after* solution 1, so it's decided then.
- **Sizing the contact drag.** At the contacts' own closing rates (about 1.5 cm/yr, against
  v_ref = 4 cm/yr) and jam factors, the RFC's drag gives a median of 0.8–0.9 × 10¹² N/m.
  That's just under the acceptance band of 10¹²–10¹³. The p90 (1.8–2.2 × 10¹²) is inside it.
  Slowing will lower the rates further, so the force check is read after solutions 0 and 1.

## Gate outcome

| Gate row | Result | Build |
|---|---|---|
| Any | — | 0 (required), 1, 5 and 6 |
| Closing rate at craton contacts within 10% of craton-free, after 1 | Not yet measurable: 1 not built. Baseline shows no braking | Re-measure after 0 and 1 |
| Buoyant-pull share under 1% on every seed | No: 4.0–6.4% | 1, with its resistance tuned |
| J under 5% on every seed | No: 22–57% | Keep the loss criterion |
