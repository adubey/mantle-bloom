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

Pending.
