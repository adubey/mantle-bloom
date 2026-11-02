# Enclosed ocean-floor basins as ocean: survey and replay comparison

Findings for [#758](https://github.com/adubey/mantle-bloom/issues/758). It classifies an enclosed
below-sea-level body as ocean by what it is, rather than by its share of the largest ocean, and
measures what that does to #754's windows D and E and to one other mature world.

**Principles:** 1, because the model boundary between ocean and lake was a share of the largest
component (5% of its nodes), not a property of the body, and 10, because a world-scale feature
near that share changes class with small changes in elevation or sea level.

**The change.** `hydrology.connected_ocean_mask` keeps the largest below-sea-level component as
ocean. Any other component is ocean when at least half its floor area is oceanic crust and it
covers at least 5 million km² (`ENCLOSED_OCEAN_MIN_AREA_KM2`). A component holding any of last
step's ocean stays ocean down to 2.5 million km² (`ENCLOSED_OCEAN_KEEP_AREA_KM2`). Lake water on a
cell that becomes ocean is exported to the ocean by the lake solve, as before; lake ice on such a
cell now melts into it first, so it leaves with the water. Eustasy then fills the basin to sea
level from the conserved budget. See
[Which water below sea level is ocean](../../docs/simulation-model/hydrology.md#which-water-below-sea-level-is-ocean).

## Survey: where the threshold falls

`bin/debug/survey_enclosed_basins.py` lists every connected below-sea-level component of at least
20 cells in a save, with its area, oceanic-crust share by area, median rock age
(`node_created_years`), depth and sill. The sill is the lowest level at which the body would join
the largest component, from a minimax flood over the hydrology graph at cell centres, so it
ignores breach notches. Nine saves, on `origin/main` `38c2460d`:

```bash
backend/.venv/bin/python bin/debug/survey_enclosed_basins.py <save> [<save> ...] \
  --min-nodes 20 --out <scratch>/survey.json
```

The largest enclosed bodies, excluding each save's main ocean
([`survey-enclosed-basins.json`](survey-enclosed-basins.json) has every body):

| Save | Age (Myr) | Area (M km²) | Share of main ocean | Oceanic share | Median rock age (Myr) | Floor below sea (m) | Sill above sea (m) |
|---|---:|---:|---:|---:|---:|---:|---:|
| seed910211954 (D) | 225.1 | **15.14** | 4.73% | 0.99 | 25.8 | 3,448 | 1,380 |
| seed910211954 (E) | 224.0 | **15.11** | 4.72% | 0.99 | 24.5 | 3,444 | 1,162 |
| seed628977350 | 226.7 | 5.01 | 1.11% | 0.08 | 226.7 | 522 | 631 |
| seed317739385 (S) | 201.6 | 2.37 | 0.53% | 0.00 | 201.6 | 367 | 1,945 |
| seed635056485 | 452.8 | 1.92 | 0.42% | 1.00 | 23.3 | 2,385 | 7,562 |
| seed910211954 (D) | 225.1 | 0.89 | 0.26% | 0.00 | 225.1 | 988 | 658 |
| seed373936961 | 244.4 | 0.68 | 0.17% | 1.00 | 7.5 | 1,341 | 1,909 |

The other saves' largest enclosed bodies are under 0.4 million km². Save SHA-256s:

| Save | SHA-256 |
|---|---|
| `seed910211954-225.1myr` (D) | `6a8d6f90959b3cc4235a48ebbd6ddb0373041f79db131d9ac9f370dd6806bff7` |
| `seed910211954-224myr` (E) | `4045b915a2d3d9e534a84089bba7277bbd488690e51d70f4e928f58c8cf674b1` |
| `seed317739385-201.6myr` (S) | `ad8d34d60ac89db97efdfd722bc3db0287db22d8917288871502fe9a7747f94e` |
| `mantle-bloom-seed373936961-244400000y` | `9980103ce9bd7733e4671bc3bec0d48359b45a8e457bc76a73dd769f8840f342` |
| `mantle-bloom-seed635056485-452800000y` | `ec27d5077bdf43b3a6f0339b5726fa4d110da9eb1f143f96a3004d79f309b045` |
| `mantle-bloom-seed628977350-226700000y` | `707d35d312ce9b5bdf1cc78c64ad880e2e3d072c80458b54a289b04e89e4acb7` |
| `mantle-bloom-seed447276134-244400000y` | `7b1490776f0d9994da68f79e2c27f686c90d0997b8e2ae4c1a5f550e83dfb7fb` |
| `mantle-bloom-seed295569319-180000000y` | `1dd920a60a5a1ff6db4f0ad8430f1467f23de5d8f55d11402fc0ed9a72089840` |
| `mantle-bloom-seed997271774-130500000y` | `a43962a149a0ec96f1dbee3e6c675be98d50422c47e0455e062aa72a7d1da4d7` |

What this shows:

- **D's basin stands alone.** No other enclosed body floored by oceanic crust exceeds
  1.9 million km², and D's is 15.1 million, about the Arctic Ocean's size. 5 million km² sits in
  that gap. It is twice the Mediterranean, larger than any of Earth's marginal seas.
- **Crust separates the next-largest bodies.** The 5.0 million km² body on seed 628977350 is 92%
  continental crust, and S's 2.4 million km² body is all continental: drowned inland basins, not
  ocean basins. That save's sea level is +4,856 m, which is unusual, and wasn't investigated.
- **D's basin isn't a hidden strait.** Its lowest pass is 1.2–1.4 km above sea level at cell
  centres. #754 found the hierarchy's notched rim between about −1,070 and +330 m in the steps it
  checked. So the rule doesn't claim a real connection: it treats an ocean-scale basin with the
  ocean model, whose surface sits at sea level, because the lake model gives it no circulation or
  ocean heat. Filling it from the world ocean is the cost.
- Filling D's basin to today's sea level takes 45.4 million km³. Over the 332 million km² main
  ocean that is about 137 m of sea level, before hypsometry and the water main already holds in
  the basin.

## Replay comparison: preregistration

Committed before any replay.

**Arms.**

| Arm | Tree | Commit |
|---|---|---|
| main | Detached `origin/main` | `38c2460d` |
| branch | `claude/issue-758-enclosed-ocean-basins` | the commit carrying this preregistration |

**Windows.** 100 kyr steps with the ice-age cycle off, as in #754. Each runs 5 Myr (50 steps),
which is also the scan for class-change frequency. D and E are #754's registered checkpoints of
seed 910211954, and S is its registered seed-317739385 checkpoint, the "other mature save": its
enclosed bodies are all continental or small, so S is a control the rule shouldn't touch. SHA-256s
are in the survey table above.

**Command.** One arm at a time, each after a `compute-budget` check (`--run sim`):

```bash
cd backend
.venv/bin/python ../bin/debug/replay_enclosed_ocean.py <checkpoint.mbworld> --myr 5 \
  --step-years 100000 --backend <arm tree>/backend --out <scratch>/<window>-<arm>.json
```

**Metrics.** From `bin/debug/replay_enclosed_ocean.py`, per step:

- **Class changes:** a connected group of cells below sea level both before and after a step
  whose `is_ocean` flipped in it, covering at least 200,000 km² (the sea tier's size), counted
  by direction. Shoreline moves aren't class changes.
- **Ocean:** node count and area. **Sea level**, and the eustasy closure: ocean volume at the
  solved level plus trapped water, minus the conserved total, against the bisection tolerance
  (1 m over one mean cell).
- **Lakes:** whole-world liquid and ice volume. **Ledgers:** lake residual over the larger of
  external supply and overflow out; land residual over rain plus meltwater.
- **Enclosed bodies** of at least 200,000 km², with area, oceanic share and class.
- **Precipitation over nearby land:** the cos-latitude-weighted mean over land climate cells
  (not ocean, and no visible lake) within 1,000 km of the largest enclosed body's cells at the
  start but more than 100 km from them, so the basin floor itself (land on main, ocean on the
  branch) isn't in the sample. On S, the largest enclosed body is the 2.4 million km² continental
  basin.
- **Runtime:** median step seconds and peak RSS.

**Predictions and thresholds.**

| # | Prediction | Supported if |
|---|---|---|
| 1 | D's and E's basin becomes ocean at step 1 and stays ocean | Branch: the basin (the 15 million km² body) is ocean at every step from 1 to 50 in D and E. Main: it is never ocean |
| 2 | Sea level falls by the basin's fill | Branch minus main sea level at step 1 between −100 and −180 m in D and E |
| 3 | Ocean grows by the basin | Branch ocean node count above main's by 3,000–4,500 at step 1 in D and E |
| 4 | Ledgers close at reclassification and after | Every branch step: normalized lake and land residuals at most 1e-9, and eustasy closure within its tolerance |
| 5 | No flicker from the new rule | Branch, D and E: exactly one class change into ocean of at least 10 million km² (step 1) and none of the basin back to lake. Branch total class changes no more than main's in each window |
| 6 | Whole-world lake volume falls | Branch whole-world lake water (liquid plus ice) below main's at step 50 in D and E |
| 7 | Nearby land gets wetter | Branch mean nearby-land precipitation over steps 1–50 above main's in D and E |
| 8 | S is untouched | S: branch sea level within 1 m of main at every step, and the same class-change count |
| 9 | Runtime is unchanged | Median step time within 10% of main's in each window |

Sizes: prediction 2 is 45.4 million km³ over about 332 million km², 137 m, less the water main's
basin holds, which swings by millions of km³ (about 10 m of sea level) step to step. Prediction 3
is the basin's 4,046 cells, less cells left dry by the lower sea. No size is predicted for 7; its
sign is the test. Main's D and E windows swing every step, so 6 compares step 50 only, and its
whole-world lake volume isn't otherwise interpreted.

**Variability and limits.** One seed carries D and E, and E's window contains D's start, so D and
E aren't independent. Main's lake volume swings by about 2.6 million km³ a step in E (#754), which
is about 8 m of sea level, so sea-level gaps under about 10 m aren't meaningful. #757's open PR
(#761) changes the same lakes and lake-ice handling; whichever of the two lands second reruns
these numbers.
