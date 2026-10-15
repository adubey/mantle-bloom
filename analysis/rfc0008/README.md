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

# Solution 0: contact drag (issue #399)

Collision drag in Pa·s/m, sized to 5 × 10¹² N/m of trench at 4 cm/yr, with each contested
node braked by max(relief factor, craton strength). Same probe and seeds, origin/main 22a2f3c
against the branch, both with the probe's plate statistics. Per-seed summaries are in
[solution0.json](solution0.json).

| Seed 7 / 11 / 2024 | origin/main | Solution 0 |
|---|---|---|
| Brake at jammed contacts, median (N/m) | 7.8 × 10⁴ / 4.9 × 10⁴ / 2.5 × 10⁴ | 3.0 × 10¹¹ / 1.4 × 10¹¹ / 1.2 × 10¹¹ |
| Closing rate, continental contacts with craton in reach (cm/yr) | 1.70 / 1.84 / 1.48 | 1.51 / 1.81 / 1.52 |
| Mean plate speed, area-weighted (cm/yr) | 3.65 / 3.35 / 3.33 | 3.33 / 3.13 / 3.21 |
| Mean plate count | 46.7 / 44.1 / 44.6 | 44.8 / 46.6 / 43.4 |
| Largest plate, mean / max share of surface | 6.6 / 18% · 7.1 / 31% · 7.7 / 22% | 6.6 / 18% · 6.8 / 31% · 7.6 / 22% |
| Collision-time drops: all (jammed pairs) | 106 (28) / 54 (22) / 111 (53) | 188 (128) / 67 (57) / 87 (59) |
| Fusions, roll + forced | 55 / 27 / 36 | 45 / 15 / 24 |
| Craton-on-craton overlap at fusion, % of start | 90% / 75% / 53% | 52% / 32% / 18% |
| Craton lost in fusions, % of start | 54% / 18% / 31% | 39% / 4% / 12% |
| Craton left at 1 Gyr, % of start | 17% / 40% / 43% | 37% / 52% / 75% |

- **The brake is real now, but light.** Its median at jammed contacts rose about 10⁶-fold,
  to 1.2–3.0 × 10¹¹ N/m, still under the RFC's 10¹²–10¹³ band. g comes only from relief and
  craton-on-craton contacts. The jam factor j, which covers buoyant crust at the trench,
  arrives with solution 1.
- **Plates slow a little.** Mean plate speed falls 3–9%, far from the 30% falsification
  limit. Plate count moves under 6%, and the largest plate is unchanged. Closing at contacts
  with craton in reach changes by −12%, −2% and +3%, short of the RFC's 25% (a target for the
  RFC as a whole).
- **Craton survives mostly through fewer fusions.** Craton left at 1 Gyr rises from 17–43% to
  37–75% of the start, and nearly all of the gain is the fusion (`delaminated_m3`) sink.
  Destruction outside fusions moves −7%, +11% and −30%. Fusions are fewer because braked pairs
  stall: collision time is dropped on jammed pairs 2.6–4.6 times as often on seeds 7 and 11.
  Solution 5 keeps that time, so part of this gain will come back as fusions. Craton loss
  should be read again after 5 and RFC 0001's placement fix.

# Solution 1: buoyant columns (issue #400)

Slab pull per subducting node scales with its column's net negative buoyancy (a full-strength
cratonic keel is neutral), and continental crust forced down at a contact that isn't
diverging jams the contact drag with j (15 km → 0, 30 km → 1). Same probe and seeds, at
5 Myr steps to 1 Gyr. The baseline is solution 0's replays ([solution0.json](solution0.json)):
no model code changed on origin/main between them and this branch. Per-seed summaries are in
[solution1.json](solution1.json).

| Seed 7 / 11 / 2024 | Solution 0 (origin/main) | Solution 1 |
|---|---|---|
| Buoyant share of slab pull | 5.4% / 5.3% / 4.7% | 0.01% / 0.01% / 0.01% |
| Brake at jammed contacts, median (N/m) | 3.0 × 10¹¹ / 1.4 × 10¹¹ / 1.2 × 10¹¹ | 3.4 × 10¹¹ / 2.4 × 10¹¹ / 1.7 × 10¹¹ |
| Brake at jammed contacts, p90 (N/m) | 7.7 × 10¹¹ / 2.8 × 10¹¹ / 3.0 × 10¹¹ | 9.2 × 10¹¹ / 9.3 × 10¹¹ / 4.6 × 10¹¹ |
| Closing rate, continental contacts, craton in reach / clear (cm/yr) | 1.51 / 1.61 · 1.81 / 1.60 · 1.52 / 1.37 | 1.50 / 1.27 · 1.41 / 1.28 · 1.20 / 1.23 |
| Closing rate, all contacts, craton in reach / clear (cm/yr) | 1.42 / 1.50 · 1.46 / 1.47 · 1.46 / 1.48 | 1.36 / 1.45 · 1.41 / 1.45 · 1.34 / 1.47 |
| Mean plate speed, area-weighted (cm/yr) | 3.33 / 3.13 / 3.21 | 3.00 / 2.92 / 2.98 |
| Mean plate count | 44.8 / 46.6 / 43.4 | 46.3 / 43.8 / 43.0 |
| Largest plate, mean / max share of surface | 6.6 / 18% · 6.8 / 31% · 7.6 / 22% | 6.7 / 19% · 6.7 / 31% · 6.5 / 22% |
| Collision-time drops: all (jammed pairs) | 188 (128) / 67 (57) / 87 (59) | 231 (168) / 43 (26) / 54 (22) |
| Fusions, roll + forced | 45 / 15 / 24 | 36 / 11 / 11 |
| Craton-on-craton overlap at fusion, % of start | 52% / 32% / 18% | 43% / 31% / 20% |
| Craton lost in fusions, % of start | 39% / 4% / 12% | 42% / 2% / 2% |
| Craton destroyed outside fusions, % of start | 29% / 52% / 26% | 24% / 43% / 20% |
| Craton left at 1 Gyr, % of start | 37% / 52% / 75% | 38% / 67% / 95% |
| J | 34% / 82% / 58% | 31% / 82% / 79% |

- **The keel fix works.** Buoyant columns supplied 4.7–5.4% of slab pull; they now supply
  0.01% (thin continental crust on a thick root, which still sinks).
- **Plates slow further, within bounds.** Mean plate speed falls another 7–10% (11–18% against
  the pre-solution-0 baseline of 3.65 / 3.35 / 3.33 cm/yr), under the RFC's 30%
  falsification limit. Plate count moves under 6% and the largest plate is unchanged.
- **The brake at the probe's craton-jammed contacts rises 1.1–1.7× at the median** and 1.2–3.3×
  at p90, still under the 10¹²–10¹³ N/m band. The jam acts where thick continental crust is
  forced down; the probe's "jammed" contacts are where craton blocks retreat, and many of
  those have a thin-crusted or upper-plate side with j = 0.
- **Solution 2's trigger is met.** Craton contacts don't close more slowly than craton-free
  ones: continental pairs with craton in reach close +18%, +10% and −2% against clear pairs,
  and all pairs −6%, −3% and −9%, within 10% on every seed by the all-pairs measure. Closing
  falls at contacts of both kinds, not specifically near cratons. By the RFC's table this
  selects solution 2 (the proximity brake).
- **Craton survives better on seeds 11 and 2024** (67% and 95% left, from 52% and 75%), with
  fewer fusions and less loss outside them. Seed 7 is flat (38%), and its collision-time drops
  on jammed pairs rise from 128 to 168, which solution 5 addresses.
- **Cost.** `contact_jam_factor` (including `forced_down_mask`) took 0.28% of a 100 Myr seed 7
  run at density 1 under cProfile; the buoyancy and drag changes 0.04%. Not measured at
  density 4.
