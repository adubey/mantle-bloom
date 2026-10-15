# RFC 0009 gate G1: does a craton mosaic give rifts a way around? (#426)

[RFC 0009](../../rfc/0009-craton-mosaic.md) proposes seeding each continent as a mosaic of
cratonic cores separated by belts that never become craton. Gate G1 asks whether the mosaic
lets RFC 0010's minimum cut route a rift between the cores. Its outcomes were fixed in the
RFC before this run.

**Outcome: neither building branch.** At the starting knobs (1,500 km core spacing, 300 km
belts), pooled over seeds 7, 11 and 2024, L = 24% and B = 0.91. The cheaper branch needs
L ≥ 20% *and* B ≤ 0.75. No seed passes it on its own. Band exposure is worse than on today's
layout with formation off (B = 0.77). Under the RFC's rule this is the falsification case:
don't build RFC 0010's cut. Belts also don't last. By 500 Myr they keep 5%, 19% and 34% of
their t = 0 area, under a third on two seeds, so under the RFC's other rule 0009 waits for
the sink RFCs (0001, 0002, 0005). Seed 2024's 34% means the stronger falsification, "every
seed under a third", is narrowly missed.

## Method

`bin/debug/probe_craton_mosaic.py` replays fresh worlds on origin/main `02d5145` (the probe
commit, `8ffdade`, adds only the script), with every default: 5 Myr steps to 1 Gyr and
100 Myr checkpoints. It runs one seed at a time and changes no saved state.

- **Mosaic seeding** (RFC 0009 solution 1). Per plate, the eligible columns are today's
  seeding host with the 500 km formation margin. Poisson-disc nuclei are placed
  `CORE_SPACING_KM` apart, from an RNG keyed on (seed, plate). Each column goes to its
  nearest nucleus by graph hops, within its connected piece of the interior. Columns within
  `BELT_WIDTH_KM / 2` of a partition boundary are belts. Everything else becomes craton, as
  seeding does today.
- **Belts** are carried by a probe-only HISTORY surface field, registered before the app
  loads. A merged or coarsened cell stays a belt if any contributor was one, which is the
  remap rule the RFC proposes.
- **Formation off** (`CRATON_FORMATION_MYR` = ∞), so the seeded layout is the only craton.
- **Every cut** that goes ahead (split or failed rift) is evaluated by
  `prototype_mincut_split._evaluate` at a 1,000 km corridor and floor 0.05, with all of RFC
  0010's admissibility gates. Strength is RFC 0010's: craton or keel, but
  `RIFT_BELT_STRENGTH` = 0 on a belt.
- **L**: the share of craton-crossing flow cuts that have an admissible minimum cut whose
  line touches no craton (no craton node at either end of a crossed link).
- **B**: craton area in the ±2.5-spacing band, using the minimum cut where it is
  admissible and the flow cut otherwise, divided by the same sum for the flow cuts.
- **Control** (`--layout blob`): today's 875 km seeding, also with formation off, so the
  two layouts differ only in where the seeded craton sits. RFC 0010's own control
  ([analysis/rfc0010](../rfc0010/README.md)) had formation on: L = 17% and B = 0.73.
- **Sensitivity** (`--belt-width-km 600`): the top of the RFC's belt-width sweep. This is
  not the gate.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py --out $S/g1_mosaic.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py --layout blob --out $S/g1_blob.json
PYTHONPATH=backend backend/.venv/bin/python bin/debug/probe_craton_mosaic.py --belt-width-km 600 --out $S/g1_mosaic600.json
```

Runs took 730–870 s each. Generation took 0.8–1.8 s with either layout. Compact numbers are
in [g1.json](g1.json).

## Results

### The gate (0–1 Gyr, 1,000 km corridor)

| Layout | Seed 7 L / B | Seed 11 L / B | Seed 2024 L / B | Pooled L / B | Craton-crossing cuts |
|---|---|---|---|---|---|
| **Mosaic, 300 km belts (gate)** | 9% / 0.94 | 37% / 0.83 | 33% / 0.98 | **24% / 0.91** | 53, 43, 27 |
| Today's layout, formation off (control) | 14% / 0.80 | 21% / 0.69 | 9% / 0.77 | 15% / 0.77 | 63, 33, 22 |
| Mosaic, 600 km belts (sensitivity) | 29% / 0.90 | 42% / 0.75 | 41% / 0.75 | 38% / 0.81 | 48, 59, 34 |

- **The mosaic raises L but not B.** Minimum cuts do run along belts. In one inspected seed 7
  cut, 37 of the 83 nodes on the cut were belt and 2 were craton. A clear line is also
  more common than on today's layout (24% against 15%). But the band a failed rift thins is
  ±312 km, wider than a 300 km belt. A cut down the middle of a belt still has cores on both
  sides of its band, so B stays near 1. Belts would have to be wider than the band, about
  625 km, before following one empties the band. Even at 600 km, B is 0.81 pooled, and only
  seeds 11 and 2024 reach 0.75.
- **The mosaic exposes more craton than today's layout.** It seeds with the 500 km
  formation margin, not 875 km, so it has more craton at t = 0 (125, 92 and 59 M km² against
  104, 82 and 46), and more of each plate is craton or belt.
- **Early cuts do better than late ones.** Over the first 300 Myr, the 300 km mosaic on
  seed 11 reaches L = 50% and B = 0.78. After 300 Myr the belts are mostly gone and both
  measures drift back toward the control (from the per-cut rows, which aren't checked in).

### Seeding (t = 0)

| | Seed 7 | Seed 11 | Seed 2024 |
|---|---|---|---|
| Seeded interior (M km²) | 150 | 111 | 71 |
| Belt share of it, 300 km / 600 km | 17% / 32% | 17% / 33% | 17% / 32% |
| Interior pieces with ≥ 2 nuclei, of which with ≥ 2 cores | 10 / 10 | 4 / 4 | 6 / 6 |
| Craton components (≥ 1 M km²), mosaic / today | 81 (61) / 23 (8) | 49 (37) / 10 (4) | 42 (28) / 12 (4) |

The seeding meets the RFC's t = 0 criteria: every piece large enough for two nuclei has at
least two cores, and belts are 10–30% of the seeded interior at the starting knobs.

### Belt survival (300 km belts)

Belt area counts belt-marked genuine continental columns without craton. "Marked on craton"
counts belt-marked cells that also carry craton, after a remap blended the two.

| Myr | Seed 7 belt (of t = 0) | Seed 11 | Seed 2024 | Seed 7 marked on craton (M km²) |
|---|---|---|---|---|
| 0 | 25.0 (100%) | 19.1 (100%) | 12.0 (100%) | 0.0 |
| 100 | 16.2 (65%) | 12.6 (66%) | 10.6 (89%) | 5.0 |
| 200 | 7.7 (31%) | 8.4 (44%) | 8.6 (72%) | 11.4 |
| 300 | 4.9 (20%) | 6.4 (34%) | 7.5 (63%) | 8.7 |
| 500 | 1.4 (5%) | 3.6 (19%) | 4.0 (34%) | 9.0 |
| 1,000 | 0.2 (1%) | 0.2 (1%) | 0.9 (8%) | 6.7 |

Genuine continental area falls too, by 69%, 43% and 33% at 500 Myr on the same runs. Belts
fall faster than that, by two routes:
- **Sinks.** Belts are non-cratonic host, so they are thinned out of the host window or
  destroyed like the rest of it (#409).
- **Merges onto craton.** When a quad merge or coarsening combines a belt cell with a
  craton cell, the result keeps both the mark (HISTORY) and part of the craton (extensive).
  On seed 7, 9 M km² of marked area sits on craton by 500 Myr. In the probe that cell still
  scores belt strength, but it counts as craton for L and B. The RFC's proposed remap rule
  ("keep a mark when any contributor has one") would make the same blend, and so would
  formation's exclusion of marked cells. The implementation would need to decide which
  state wins.

At 600 km, belts keep 7%, 15% and 26% of their area at 500 Myr: wider belts don't last
proportionally longer.

## What the outcome means

- **By the pre-registered rule**, G1 lands in the "Otherwise" branch: don't build RFC 0010's
  minimum cut. The belt-survival rule independently says the mosaic waits for the sink RFCs.
- **The falsification's stated reading** is "the flow puts regime boundaries through
  cratons". It fits L less well than B. L does rise with a mosaic (15% → 24%, or 38% with
  600 km belts), so the cut's geometry was part of the obstacle. B fails mainly because the
  operational band is wider than a belt. That is a property of how the model thins a failed
  rift (`FAILED_RIFT_BAND_MULT`), not of where the flow puts its boundaries.
- **Belts vanish within a few hundred Myr whatever their width.** No layout fix survives
  that. This matches #409's finding that non-cratonic host is what the model's sinks remove.

## What this does not show

- Whether a mosaic would pass G1 once the sinks are repaired (RFCs 0001, 0002, 0005), with
  belts added at collisions and rifts (solution 2's writers, not in this probe), or with
  formation on. The probe seeds belts only at generation.
- Whether B would fall with a failed-rift band scaled to the weak zone it follows. G1 fixed
  the band at today's `FAILED_RIFT_BAND_MULT`.
- How much rifting loss a mosaic saves after the split (#388). As in RFC 0010's probe, stretch
  donation and divergent thinning aren't followed.
