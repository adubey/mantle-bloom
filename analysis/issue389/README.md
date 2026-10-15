# Issue #389: steering the split cut off cratons (RFC 0003 solution 1)

[RFC 0003](../../rfc/0003-rifting-around-cratons.md) solution 1 tilts a plate split's
flow-derived cut by up to 20° onto the candidate whose band crosses the least lithospheric
strength (craton, or a keel thicker than 100 km), as long as the tilted cut still separates
the same flow regimes. This A/B measures it against origin/main on the RFC's tuning seeds.

**Outcome: the RFC's "no admissible alternative" falsification case holds.** In 158
flow cuts that crossed craton on the three seeds, a craton-free cut passing the
admissibility tests existed in 19 (12%) at the default 20° tilt and 10% misfit, and in at
most 30 (19%) anywhere in the RFC's sweep ranges (tilt 10–30°, misfit 5–20%). Steering cut
the share of splits crossing craton by 13–21% relative, against the 75% criterion. Mean
rifting loss fell 22% (8.55% to 6.70% of starting craton per Gyr), but seed-to-seed scatter
between variants is as large as that drop (see Noise).

## Method

`bin/debug/steer_split_ab.py` wraps `rebaseline_rfcs._measure` (the #403 harness):
`world.generate_world(seed)` with all defaults, cratons on, 5 Myr steps to 1 Gyr, 100 Myr
checkpoints, seeds 7, 11 and 2024, one replay at a time. Its read-only taps add:

- each successful split's daughters' RMS flow-fit residual against their set poles;
- time in `maybe_split_plate`, and in the steering itself;
- a **gate census** for every cut that went ahead (split or failed rift): did the
  flow-derived cut's band hold craton, did the chosen cut's, and, for tilts 10, 20, 30 and
  45°, did a craton-free candidate exist and pass each admissibility test, with misfit caps
  of 5, 10 and 20%. The census recomputes the steering's arrays on the run's own trajectory;
  it does not change it.

Runs:

| Variant | Model commit | Overrides |
|---|---|---|
| main | `02d35b9` (origin/main) | none |
| steered | `168b060` | none |
| steered, band ×2 | `168b060` | `SPLIT_CUT_BAND_MULT` = 5.0 (scoring band only) |

Commands (scratch outputs; `S` is a scratch directory):

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/steer_split_ab.py \
    --out $S/branch.json --save-dir $S/saves-branch
PYTHONPATH=backend backend/.venv/bin/python bin/debug/steer_split_ab.py \
    --set SPLIT_CUT_BAND_MULT=5.0 --out $S/band2x.json --save-dir $S/saves-band2x
cp bin/debug/steer_split_ab.py ../mantle-bloom-origin-main/bin/debug/
(cd ../mantle-bloom-origin-main && PYTHONPATH=backend backend/.venv/bin/python \
    bin/debug/steer_split_ab.py --out $S/main.json --save-dir $S/saves-main)
```

The baseline needs the copy inside the origin/main worktree: `compare_cratons.py` puts its
own tree's `backend` first on `sys.path`, overriding `PYTHONPATH`. A first baseline run
without the copy imported the branch's model; it is discarded, but it reproduced the steered
run exactly (every ledger, cut count, plate count and residual), as did two later steered runs
with the census taps, so steered replays are deterministic. (Save-file hashes differ between
repeats; the saves aren't byte-stable, so they aren't the check.)

Compact numbers: [steering_ab.json](steering_ab.json).

## Results (0–1 Gyr)

| Seed | Variant | Rifted, % of start craton | Div / stretch / failed | Splits through craton | Craton share of split-band area | Splits | Mean plates | Post-split residual, cm/yr |
|---|---|---|---|---|---|---|---|---|
| 7 | main | 9.61% | 1.34 / 7.46 / 0.81 | 39 of 107 (36.4%) | 25.6% | 107 | 43.5 | 1.96 |
| 7 | steered | 8.17% | 1.14 / 6.41 / 0.63 | 34 of 110 (30.9%) | 18.8% | 110 | 47.1 | 1.83 |
| 7 | band ×2 | 10.11% | 1.33 / 8.32 / 0.46 | 35 of 107 (32.7%) | 20.1% | 107 | 48.7 | 1.95 |
| 11 | main | 8.06% | 1.95 / 5.70 / 0.41 | 24 of 95 (25.3%) | 11.7% | 95 | 45.5 | 1.51 |
| 11 | steered | 6.67% | 1.63 / 4.72 / 0.32 | 20 of 100 (20.0%) | 9.3% | 100 | 43.1 | 1.65 |
| 11 | band ×2 | 5.69% | 1.16 / 3.85 / 0.68 | 20 of 98 (20.4%) | 9.0% | 98 | 43.8 | 1.51 |
| 2024 | main | 7.97% | 1.28 / 5.87 / 0.82 | 29 of 103 (28.2%) | 9.3% | 103 | 42.8 | 1.63 |
| 2024 | steered | 5.26% | 1.30 / 3.44 / 0.52 | 23 of 94 (24.5%) | 6.8% | 94 | 42.0 | 1.73 |
| 2024 | band ×2 | 7.22% | 1.61 / 5.06 / 0.55 | 21 of 95 (22.1%) | 5.7% | 95 | 44.4 | 1.73 |

Means: rifting 8.55% (main), 6.70% (steered, −22%), 7.67% (band ×2, −10%). Splits through
craton 30.0%, 25.1% (−16% relative) and 25.1%. Post-split residual 1.70, 1.74 (+2%) and 1.73
cm/yr. Craton balance error is at most 0.02 km³ in every run.

Origin/main is already well below #403's `2562f63` rates (13.17%, 8.97%, 15.62%; mean
12.59%). #410 (continental contact drag) merged in between and moved the trajectories, so
this A/B compares against `02d35b9`, not the RFC's pinned numbers.

### Gate census (steered, all three seeds)

Of 585 cuts that went ahead, 158 flow-derived cuts had craton in their band; 144 chosen cuts
still did. For those 158:

| Max tilt | Craton-free candidate exists | …passes size and fit-rank | …passes pole separation | …passes residual | Admissible, misfit ≤ 5% / 10% / 20% |
|---|---|---|---|---|---|
| 10° | 32 | 21 | 32 | 32 | 11 / 17 / 21 |
| **20°** | 64 | 38 | 64 | 52 | 9 / **19** / 27 |
| 30° | 101 | 49 | 100 | 67 | 11 / 17 / 30 |
| 45° | 139 | 46 | 138 | 72 | 8 / 17 / 30 |

Within 20°, three cuts in five have no craton-free candidate at all. A wider cone finds one
in 88% of cases, but such a cut moves more than 20% of the plate to the other side, or leaves
a half under `SPLIT_MIN_NODES`; the residual test removes some more. Pole separation almost
never binds. The cratons sit where the flow's
regime boundary runs, so the gates that keep the cut faithful to the flow are the ones that
stop it missing the craton. This is the RFC's "no admissible alternative" case: "the flow
field puts regime boundaries through cratons. The fix belongs in the mantle model, not the
cut."

Steering still lowers the craton area inside split bands by 21–27% (it often slides the band
toward the craton's edge), which is a plausible route for the drop in stretch loss.

### Noise

Each variant is one trajectory per seed. Seed 2024 rifts 5.26% with the default scoring band
and 7.22% with a wider one, though the two variants move a similar number of cuts. Seed 7
rises above main with the wider band. A spread of ±2 percentage points between variants that
differ only in where a handful of cuts fall is as large as the mean drop, so −22% is not
established as the effect of steering. More seeds or repeat runs with perturbed seeds would be
needed to resolve it.

### Acceptance criteria (RFC 0003)

| Criterion | Result |
|---|---|
| Rifting ≤ 2% of live craton per Gyr | Not met: 6.70% mean (steered). Solution 3 was always needed |
| Split cuts through craton −75% | **Not met**: −16% mean (−15%, −21%, −13%) |
| Splits and mean plates within 10% | Met: splits +3%, +5%, −9%; plates +8%, −5%, −2% |
| Post-split residual within +10% | Met: +2% mean (−7%, +9%, +7% per seed) |
| < 100 ms per attempt at 10⁵ nodes | Met: 25.6 ms, synthetic 10⁵-node cap, full path (strength excluded, 0.5 ms) |
| Step time within 2% | Met at default density: steering took 0.25–0.34 s of 750–800 s per run (< 0.05%); `maybe_split_plate` 3.6–4.0 s in both variants. Density 4 not measured |
| Ledgers balance; deterministic | Met (see Method) |
| Seed 456862899 to 3 Gyr | Not run: held pending the decision below |

## What this means

Solution 1, as designed, can't reach its own cut-share criterion on the current model. Phase
1b's F = 42% was an upper bound that assumed a steered cut could miss the craton; the census
shows it usually can't without abandoning the flow's regime boundary. Options:

1. Don't build solution 1. Record the falsification in RFC 0003 and size solution 3 to the
   whole loss (8.55%/Gyr on origin/main).
2. Keep steering for its partial effect (less craton in the band), and size solution 3 to the
   remainder. That needs more seeds first, to show the effect is above noise.
3. Pursue the RFC's own redirect: why the flow's regime boundaries sit under cratons
   (mantle model).
