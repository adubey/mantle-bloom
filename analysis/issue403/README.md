# Issue #403: RFCs 0001–0007 on the year-based split pace

#358 counts plate age, and so the split cooldown, in years (`SPLIT_MIN_AGE_YEARS` = 2 Myr).
Before it, the cooldown was 20 steps: 100 Myr at the 5 Myr steps of the tuning runs, and
200 Myr for most of the 19.2 Gyr run behind RFCs 0001–0007 (0.1 Myr steps to 21.8 My, then
10 Myr). Every number those RFCs quote was measured at that slower split pace. This is the
baseline they are re-stated against.

## Method

On origin/main at **2562f63**, default density, 5 Myr steps, one run at a time:

- `bin/debug/compare_cratons.py --seeds 7,11,2024 --variants on --years 1e9 --every 100e6 --workers 1 --save-dir <dir>`
- `bin/debug/compare_cratons.py --seeds 456862899 --variants on --years 3e9 --every 100e6 --workers 1 --save-dir <dir>`
- `bin/debug/rebaseline_summary.py`, which adds `ledger_history.py`'s tables from each final
  save's `stats_history` (one sample per step) to the checkpoint rows. Its output is
  [baseline.json](baseline.json).

`--save-dir` is new in this change: it writes each run's final world so its whole history
can be read. The saves aren't checked in.

The three tuning seeds reproduce RFC 0008's baseline at 1b4f686
(`analysis/rfc0008/craton_sinks_1gyr.json`) exactly: the same craton volume at 1 Gyr, sink
shares, events and plate counts. Nothing merged between 1b4f686 and 2562f63 moves them.

**The regenerated seed 456862899 isn't the 19.2 Gyr save's world.** World generation has
changed since that run. It starts with 3.9 × 10⁹ km³ of continental material, against the
save's 7.9 × 10⁹, 82 M km² of genuine continental area and 15% land. Its numbers below are
the current model's behaviour on that seed, not a re-run of the save.

## Topology

| | Seed 7 | Seed 11 | Seed 2024 | 456862899, Gyr 0 / 1 / 2 |
|---|---|---|---|---|
| Splits per Gyr | 115 | 109 | 104 | 104 / 77 / 80 |
| Merges per Gyr | 55 | 27 | 36 | 12 / 1 / 2 |
| Failed rifts per Gyr | 84 | 80 | 85 | 94 / 79 / 58 |
| Plates: start / mean / end | 20 / 46 / 51 | 9 / 43 / 57 | 11 / 44 / 60 | 9 / 53 / 58 |
| Split cuts with craton in the band | 39% | 24% | 30% | 8% (3 Gyr) |

Before #358, on the same tuning seeds, #385's baseline had 28–37 splits and 15.5–18.6 mean
plates per Gyr (`analysis/issue385/`).

## Craton sinks

Share of the craton volume at the window start. For the tuning seeds the window is the whole
run (0–1 Gyr).

| Sink | Seed 7 | Seed 11 | Seed 2024 | 456862899, Gyr 0 / 1 / 2 |
|---|---|---|---|---|
| Craton at start (M km³) | 4,098 | 3,170 | 1,724 | 1,449 / 355 / 239 |
| `delaminated_m3` | **54%** | 18% | **31%** | **28%** / 0% / 0% |
| `eroded_m3` | 9.7% | **21%** | 10% | 21% / 3.8% / 9.9% |
| `rifted_m3` | 13% | 9.0% | 16% | 9.0% / 7.1% / 12% |
| `subducted_m3` | 5.9% | 14% | 9.3% | 15% / **22%** / **33%** |
| `collision_reworked_m3` | 1.3% | 2.1% | 1.3% | 1.8% / 0.9% / 2.1% |
| Total | 84% | 65% | 67% | 75% / 34% / 57% |
| `formed_m3`, M km³ (0–500 / 500–1,000 Myr) | 34 / 0 | 162 / 0 | 167 / 0 | 0.7 / 0, and 0 after |

Rifting by path, as a share of the start craton: divergent 3.2%, 2.8%, 3.6%; stretch 9.0%,
5.6%, 11%; failed 0.9%, 0.6%, 1.3% (seeds 7, 11, 2024). On seed 456862899: 2.9% / 2.5% /
5.6% divergent, 5.5% / 4.7% / 6.2% stretch, 0.7% / 0% / 0% failed.

- **`delaminated_m3` is the largest sink while continents are large** (seeds 7 and 2024,
  and the first Gyr of 456862899). RFC 0008's gate found it is the merge cap.
- **Merge losses come in large steps.** Steps that book more than 1% of the live craton to
  `delaminated_m3`: 20, 7 and 5 on the tuning seeds, and 6 on 456862899 (all before 420 My).
  The largest single step takes 17% (seed 7, 365 My), 13% (seed 11, 565 My), 14% (seed 2024,
  715 My) and 14% (456862899, 125 My) of the live craton. A step can hold more than one
  merge, so these are upper bounds per merge.
- **`subducted_m3` takes over once merges stop.** On 456862899, merges fall to 1–2 per Gyr
  after the first Gyr, `delaminated_m3` books nothing, and subduction is the largest sink
  (22% and 33% per Gyr), as the cratons shrink to fragments.
- **`eroded_m3` is as large as rifting**, 10–21% in the first Gyr, and no RFC targets it (#413).
- **Formation stops within 500 Myr** on every run. Seed 456862899 forms 0.7 M km³ once, at
  255 My.
- **Craton share of genuine continental area rises from loss, not formation**: 40%, 51%,
  38% and 45% at the start; 87%, 92%, 85% and 89% at 500 Myr; 94–97% at 1 Gyr on the tuning
  seeds. Non-cratonic continental crust disappears faster than craton.
- Ledgers balance: craton balance error at most 1.2 km³.

## Continental material

Fluxes in km³/yr (`continental_material_m`).

| Run, window (Myr) | Inventory at ends (10⁹ km³) | Juvenile | All sinks | `deeply_subducted` | `collision_subducted` | `rift_thinned` |
|---|---|---|---|---|---|---|
| 7, 0–100 | 12.0 → 10.9 | 10.7 | 22.6 | 11.7 | 5.8 | 2.7 |
| 7, 100–500 | 10.9 → 3.7 | 8.3 | 27.9 | 12.3 | 13.1 | 2.0 |
| 7, 500–1,000 | 3.7 → 1.7 | 5.0 | 9.8 | 6.4 | 2.3 | 1.0 |
| 11, 0–100 | 6.9 → 6.4 | 10.4 | 16.8 | 15.6 | 0.0 | 1.0 |
| 11, 100–500 | 6.4 → 3.8 | 7.3 | 15.0 | 12.6 | 1.4 | 1.0 |
| 11, 500–1,000 | 3.8 → 1.9 | 4.2 | 9.0 | 7.4 | 1.0 | 0.6 |
| 2024, 0–100 | 6.5 → 6.3 | 6.8 | 9.8 | 8.3 | 0.0 | 1.3 |
| 2024, 100–500 | 6.3 → 4.6 | 4.9 | 10.2 | 6.4 | 2.1 | 1.4 |
| 2024, 500–1,000 | 4.6 → 1.6 | 3.1 | 9.7 | 6.4 | 2.5 | 0.8 |
| 456862899, 0–100 | 3.9 → 3.4 | 7.0 | 12.3 | 10.3 | 0.9 | 1.0 |
| 456862899, 100–500 | 3.4 → 1.0 | 3.4 | 10.1 | 7.6 | 1.8 | 0.7 |
| 456862899, 500–1,000 | 1.0 → 0.57 | 1.7 | 3.0 | 2.9 | 0.0 | 0.2 |
| 456862899, 1,000–2,000 | 0.57 → 0.44 | 1.3 | 1.9 | 1.7 | 0.0 | 0.1 |
| 456862899, 2,000–3,000 | 0.44 → 0.19 | 1.1 | 1.6 | 1.4 | 0.0 | 0.1 |

- The inventory falls **73–86% in the first Gyr** on every run (the 19.2 Gyr save: 47%).
  Genuine continental area falls from 82–260 M km² to 10–27 M km², and land from 15–45% to
  2–7% of the planet.
- Sinks run at 9.8–28 km³/yr while the inventory is above 1.5 × 10⁹ km³ (the save: 8–11).
- `deeply_subducted_m3` is still most of the sinks (44–95% per window), but
  `collision_subducted_m3`, which books the merge's and the suture's unplaced remainder, is
  47% of seed 7's sinks over 100–500 Myr and 21–26% on seed 2024 after 100 Myr.
- Juvenile additions are 3–11 km³/yr while continents are large, and 1.1–1.7 km³/yr on
  456862899 after 500 Myr (the save: 0.4 over its last 7 Gyr).
- Continental material on oceanic hosts at 1 Gyr: 16%, 12%, 21% and 6%.

## What changes for each RFC

Recorded in each RFC, with the pre-#358 figures kept and labelled.

| RFC | Outcome |
|---|---|
| 0001 | Premise stronger: merges are the largest craton sink on two of three seeds, 5–20 merge-sized losses per Gyr, the largest 13–17% of live craton. The gate no longer needs the save's 3,052 My event: seed 7 at 365 My is a regenerable substitute. Recommendation unchanged. |
| 0002 | `subducted_m3` is the 3rd or 4th sink in the first Gyr, but the largest once merges stop (22–33% per Gyr on 456862899). The recommendation stands; the "today" figure is now 15–33%, not 6–18%. |
| 0003 | Rifting 9.0–16% on the tuning seeds and 7.1–12% per Gyr on 456862899 (was "4–6% today"). Phase 1b (#388, PR #392) already re-selected the branch on this pace. |
| 0004 | **Gate outcome changes.** Formation still stops, but the no-ratchet criterion (craton share of genuine continental area 30–65%) already fails on main, from losing non-cratonic crust, not from formation. Starting share is 38–51%, not 35%. The criterion, and solution 5's trigger, need revisiting: #411. |
| 0005 | **Gate scope changes.** Sinks are 2–3× the save's; P is −73% to −86% in the first Gyr. `collision_subducted_m3` is a large channel on two seeds and needs a place in the split. The acceptance criterion "while the inventory is above 5 × 10⁹ km³" never applies on the regenerated 456862899 (3.9 × 10⁹ at start). #412. |
| 0006 | Baselines re-stated. The world now reaches a near continent-free state (2% land) within 1 Gyr, not over 12 Gyr. Recommendation unchanged. |
| 0007 | Baseline re-stated (juvenile 1.1–1.7 km³/yr once continents are gone). Recommendation unchanged. |
| — | `eroded_m3`, 10–21% of craton in the first Gyr, has no RFC: #413. |

## Not measured

The step-counted cadences (`DEFRAG_INTERVAL_STEPS`, `GAP_FILL_INTERVAL_STEPS`,
`MAGMA_TRANSPORT_INTERVAL_STEPS`, the once-per-step merge roll) weren't varied. All runs here
use 5 Myr steps, so they are internally consistent, but they would move at other step
lengths.
