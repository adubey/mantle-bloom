# Issue #521: rebalance both halves when a merged lake splits

`lakes._resolve` used to write both children of a splitting lake full to the saddle. This
branch resolves each child at its own balance instead. The replays below compare the branch
with `origin/main` on two saves.

## Preregistration

Written before any replay was run.

- **Saves.**
  - `seed308161711-121.5Myr` (SHA-256 `c5572514…214d6e0`, 128,340 nodes), the #519 evidence
    save, tracking sea 26 and lake 13.
  - `seed374464662-132.6Myr` (SHA-256 `0f8c50c9…a4279c`), a held-out mature save, tracking
    its three largest bodies by saved volume.
- **Runs.** 10 × 100 kyr steps from each save, on `origin/main` `59d2ef3` and on this branch,
  one at a time, using [`bin/debug/measure_lake_split_reset.py`](../../bin/debug/measure_lake_split_reset.py).
- **Primary metric:** global visible lake volume per step, and its max/min over the 10
  steps. #519's H7 check found the saddle reset adding 4.7–5.5 M km³ per step, about the
  whole visible lake volume. #519 knockout K2b (this change *plus* K1, the frozen rule
  judged from wet cells only) cut the max/min swing of global lake area from 2.4× to 1.4×.
- **Expected direction.** On the evidence save, branch lake volume lower at every step and
  the swing in lake area and volume smaller. This run lacks K1: a frozen merged parent still
  dries before it splits, but each child now gets its own frozen check, so children without a
  cold cell keep their own water. The swing could therefore stay above K2b's 1.4×.
  - Sea 26's largest wet extent well under main's 1,713 nodes (K2b: ≤ 5).
  - Lake 13 no longer resting on its saddle, and declining or steady instead of alternating.
  - Frozen-sink inflow lower (K2b: about 75% lower); ocean discharge within 2% of main.
  - Split events higher (K2b: 42–78 per step against 12–17), merges similar.
  - Step runtime within 10% of main.
- **Decision.** The fix holds if branch lake volume stays below main's at every step on the
  evidence save, sea 26 never floods to more than 100 nodes, lake 13 isn't on its saddle at
  the end, and ocean discharge moves by less than 2%. A merge/split count that grows step after
  step without levelling off would count as a returned oscillation.
- **Variability.** Each arm is deterministic from the save, so there is no seed noise. The
  baseline's own step-to-step swing (#519: 2.4× in area) is the scale against which branch
  changes are judged. The held-out save guards against a result specific to one world.

## Results

Both saves, 10 × 100 kyr, `origin/main` `59d2ef3` against branch `33f4079`, each run alone
on the machine. The later review commit only changes a comment in `lakes.py`, so these runs
still apply. Per-step numbers are in [`compare.json`](compare.json).

```bash
cd backend
.venv/bin/python -I ../bin/debug/measure_lake_split_reset.py <local-save> --steps 10 \
    [--track 26 13] --backend <tree>/backend --out <scratch-output>.json
```

### Evidence save `seed308161711-121.5Myr` (sea 26, lake 13)

| Measure | main | branch |
|---|---:|---:|
| Global lake volume, min–max (max/min) | 1.36–5.98 M km³ (4.4×) | 0.40–0.76 M km³ (1.9×) |
| Global lake area, min–max (max/min) | 6.8–28.2 M km² (4.1×) | 1.5–3.6 M km² (2.3×) |
| Sea 26, largest wet extent | 1,712 nodes (step 3) | 43 nodes (step 2) |
| Lake 13, steps on its saddle; volume at step 10 | 9 of 10; 0.54 M km³, after 4.4–5.2 M km³ | 0; 118 km³, falling every step from 0.22 M km³ |
| Frozen-sink inflow, median steps 2–10 | 106,700 m³/s | 9,700 m³/s |
| Ocean discharge, median steps 2–10 | 448,100 m³/s | 357,100 m³/s (−20%) |
| Merge / split events per step | 0–12 / 14–19 | 2–17 / 46–104 |
| Step runtime, mean; peak RSS | 30.0 s; 769 MB | 14.1 s; 1,248 MB |

### Held-out save `seed374464662-132.6Myr` (three largest bodies: 5, 10, 21)

| Measure | main | branch |
|---|---:|---:|
| Global lake volume, min–max (max/min) | 5.32–5.37 M km³ (1.01×) | 0.97–1.16 M km³ (1.19×) |
| Global lake area, min–max (max/min) | 6.6–8.1 M km² (1.22×) | 1.7–2.2 M km² (1.30×) |
| Lake 5, steps on saddle; volume at step 10 | 10; 3.24 M km³ (evaporation 2.1–3.2× inflow) | 0; 0.07 M km³ |
| Lake 10, steps on saddle; volume at step 10 | 1; 1.43 M km³ | 8; 0.90 M km³, rising |
| Lake 21, steps on saddle; volume at step 10 | 10; 0.31 M km³ | 1; 39 km³ |
| Frozen-sink inflow, median steps 2–10 | 19,200 m³/s | 60,000 m³/s |
| Ocean discharge, median steps 2–10 | 963,700 m³/s | 893,100 m³/s (−7%) |
| Merge / split events per step | 0–2 / 9–12 | 0–2 / 18–74, falling |
| Step runtime, mean; peak RSS | 6.9 s; 1,245 MB | 6.2 s; 1,294 MB |

### What the numbers mean

- **Most of main's visible lake water was the reset's.** On both saves, branch lake volume is
  80–90% lower from the first step on, and lower than main's at every step. This matches the
  #519 H7 estimate of 4.7–5.5 M km³ added per step on the evidence save, about the whole
  visible volume. The bodies main held on their saddles while they evaporated 2–3× their
  inflow (lake 13, lake 5, lake 21) drain once nothing refills them.
- **Flicker.** Sea 26 no longer floods (43 nodes at most, against 1,712), and the evidence
  save's swing in lake area falls from 4.1× to 2.3×. That's above K2b's 1.4×, as preregistered
  could happen without K1. On the held-out save the swing rises slightly (1.22× → 1.30×), mostly
  from the first step's drop; over steps 2–10 it's 1.12×.
- **Lake 10 rests on a saddle in 8 of 10 branch steps.** This is the churn the review found: a
  fed child fills to its own rim, the saddle, reads as merged next step, and splits again. Its
  volume rises by its own inflow (0.69 → 0.90 M km³); nothing is reset into it.
  `test_a_child_fed_to_the_saddle_re_splits_each_step_without_creating_water` pins this case,
  and #520 part 1 addresses the false merge behind it.
- **Ocean discharge fails the preregistered 2% bound** (−20% and −7%), because the branch
  has less rain on land, not because water goes missing in routing. My first explanation was
  that the reset's water had left basins near their rims, so runoff spilled to the ocean. A
  rerun of the evidence pair with inflow to unfrozen sinks added ruled that out: that inflow
  falls too (median 448,200 → 201,900 m³/s), so less water reaches the ocean *and* the land
  sinks. A second rerun added the step's land precipitation and Budyko runoff, read from the
   call in :

  | Evidence save, median steps 2–10 | main | branch |
  |---|---:|---:|
  | Land liquid precipitation | 2,902,000 m³/s | 1,701,000 m³/s (−41%) |
  | Land runoff | 1,229,000 m³/s | 635,000 m³/s (−48%) |
  | Ocean discharge | 448,100 m³/s | 357,100 m³/s (−20%) |
  | Unfrozen-sink inflow | 448,200 m³/s | 201,900 m³/s (−55%) |

  The two arms are equal at step 1, when both still start from the saved climate. They part
  from step 2, once the branch's lake area has fallen from about 22 to about 2 M km². The
  extra lakes main kept were an evaporation source for the climate, and the rain from it fed the
  runoff. With them gone, land rain and runoff roughly halve. A smaller share of the reduced
  runoff now ends in land sinks, so ocean discharge drops by less. The reruns repeat the first
  run's ocean discharge to the m³/s, so the arms are deterministic. Land precipitation wasn't
  measured on the held-out save, where lake area falls from about 7 to about 2 M km² and ocean
  discharge falls 7%.
- **Frozen-sink inflow** falls 91% on the evidence save, as K2b found, but rises about 3× on
  the held-out save. That rise is not explained yet.
- **Events and cost.** Split events rise 3–5× without climbing across steps, and merges are
  unchanged. Step runtime is lower on both saves. Peak RSS is 62% higher on the evidence save
  and 4% higher on the held-out one; the snapshots are member-sized copies, so they don't
  explain a 480 MB rise, and that difference is unattributed.

### Against the decision rule

| Criterion | Result |
|---|---|
| Branch lake volume below main at every step (evidence) | Met |
| Sea 26 never floods past 100 nodes | Met (43) |
| Lake 13 not on its saddle at the end | Met |
| Ocean discharge within 2% | **Not met** (−20%, −7%): less land rain once the reset's lakes stop evaporating, above |
| Merge/split counts don't grow step after step | Met |
| Runtime within 10% | Met (faster) |
