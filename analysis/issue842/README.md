# Issue #842: trend-adaptive march as the D2 default (round 3)

Round 3 of the #822 march-cost work (rounds 1–2: `analysis/issue822/`, PR #839). It asks
whether the trend-adaptive march at a 1 kyr longest substep (preset `trend-adapt-1k`)
can replace the fixed-250-yr D2 march as the production default. Fixed-250-yr D2 stays
selectable as the reference. The tests are RFC 0015's own criteria plus #765's
invariant. The owner decides the revision; this README only reports.

## Preregistration

Written and committed before any round-3 run. Data already seen: rounds 1–2 on gen5 (seed
5, density 0.25), gen30 (seed 30, 8 plates) and the diagnosis saves 997271774 @ 87.4 Myr
and 649658380 @ 6.4 Myr, up to 4 steps each. `trend_window = 30` was chosen from a gen5
trace, and the 1 kyr ceiling was chosen by the owner from #741. Nothing has been run on
seed 374464662 or on the held-out set. Two tool smoke tests ran before this commit: one
`trend-adapt-1k` step of 649658380 with ice ages at 8 Myr (E1 runs that save with ice ages
off), and the ensemble mode on gen5, which is not in round 3's set.

**Arms:** `d2` (the production march) and `trend-adapt-1k`. E3 also runs `trend-adapt`
(4 kyr ceiling) for context; it has no part in the decision rule.

**Tree:** branch `claude/issue-842-adaptive-march-default`, stacked on PR #839. The arms
differ only by `--preset`. `bin/debug/march_cli.py` adds `--ice-age-period-myr`; saves
store `ice_age_period_years = 0`, and nothing earlier set it in replays.

### Runs

| Run | What | Reads |
|---|---|---|
| E1 diagnosis replays | `replay_789_march.py`, 10 × 0.25 Myr: 997271774 @ 87.4 and 649658380 @ 6.4 (ice ages off), 374464662 @ 132.6 (`--ice-age-period-myr 8`) | Per-step glacier and land ledger closure; march and step time; ice volume and area after 10 steps vs `d2` |
| E2 step robustness | 997271774: 10 Myr as 40 × 0.25 Myr and as 10 × 1 Myr, both arms | RFC criterion: ice volume and area at 10 Myr agree within 10% between step lengths |
| E3 stability ensemble | `probe_822_march_cost.py --ensemble`: an unperturbed member plus 12 opening-ice members (±0.01–0.5%). Run on 997 after 3 warm-up steps (step 3, the long one), 649 after 1 (step 1) and gen30 (step 0). Arms: `d2`, `trend-adapt-1k`, `trend-adapt` | #765 invariant: largest single-substep volume move per member; members that **flip**, i.e. whose closing ice differs from the same arm's unperturbed member by > 2% (4× the largest perturbation) |
| E4 convergence reference | `d2 --substep-years 125` and `62.5` on gen30 (1 step) and 997 (4 steps, step 3 read) | How far 250 yr `d2` is from itself at shorter substeps, against `trend-adapt-1k`'s distance from 250 yr `d2` |
| E5 held-out acceptance | As E1, on 447276134 @ 116.3 Myr (`--ice-age-period-myr 8`), 628977350 @ 226.7 and 295569319 @ 180.0 (from B2; the RFC lists 116.4 and 180.1, these are the nearest saves) | Same as E1. Never used for any knob or threshold choice |
| E6 cost targets | #822's monster tests (`test_sample_is_ocean_and_is_sea_healpix_mode_agree_closely_with_kdtree`, `test_spawning_is_deterministic_per_seed`) with the march on, both arms | #822's acceptance: each under ~60 s |

### Metrics and expected direction

- **Cost (primary):** median march time per step, and median whole-step time, as a ratio
  to `d2`. Expected 0.2–0.7× on E1/E5 steps that don't settle quickly, and ~1.0× on steps
  that settle in under ~30 substeps (both arms share the pair test there).
- **Ledgers:** glacier closure |residual| / accumulation on every step (RFC: ≤ 0.1%); land
  closure ≤ 0.1%. Expected: land at round-off on both arms. Glacier ≤ `d2`'s own, since a
  trend hold averages a whole window. `d2` itself was 0.4–2% on some #789 steps.
- **Ice after 10 steps (E1, E5):** relative difference in closing ice volume and area vs
  `d2`. Expected within ~2%. Trajectories diverge step by step, so this is a drift bound,
  not a pointwise one.
- **Step gap (E2):** |V(0.25 Myr) − V(1 Myr)| / V at 10 Myr, and the same for area.
  Expected: `trend-adapt-1k`'s gap is close to `d2`'s (#765 measured 9.71% for `d2` on
  this save with drifting noise; this run uses main's climate).
- **Ensemble (E3):** expected no single-substep move > 20% in any arm. Expected 0 flips for
  `d2` and `trend-adapt-1k`; `trend-adapt` (4 kyr) may flip on gen30.
- **Convergence (E4):** expected `d2`'s 250 → 125 yr difference on gen30 to be similar to
  round 2's −7.3% ice, and much smaller on 997 step 3.

### Decision rule ("propose the flip")

All of these must hold for `trend-adapt-1k`:

1. **E1 and E5, every step:** glacier closure ≤ max(1e-3, `d2`'s on the same step + 1e-4),
   and land closure ≤ 1e-3.
2. **E1 and E5:** median march time ≤ 0.6× `d2` over the steps where `d2` marches more
   than 30 substeps. If no such step exists on a save, the save only has to stay ≤ 1.1×.
3. **E1 and E5:** closing ice volume after 10 steps within 2% of `d2`. Where E4 shows a
   larger self-difference for `d2` (250 vs 125 yr), within that difference instead.
4. **E2:** step gap ≤ 10%, or ≤ `d2`'s own gap.
5. **E3:** no member with a single-substep move > 20%, and no more flips than `d2`'s
   count plus 2.

If any part fails, the round reports it and the owner decides. There is no rerun with
changed thresholds, and no flip. E6 is reported against #822's acceptance and doesn't gate
the flip.

### Controls and variability

- One run at a time (`bin/compute_budget.py check --run sim` first).
- The accuracy metrics are deterministic for a given tree, save and preset.
- Wall time varies with load: a ratio within ~10% of 1 counts as no change.
- A `d2` parity check (round 1's gen5 closing depths, bit for bit) runs first, to show the
  round-3 tree's `d2` is unchanged.
- Outputs go to `~/mantle-bloom-captures/issue-842/` with `SHA256SUMS`. The compact
  numbers go in `analysis/issue842/round3.json`.
