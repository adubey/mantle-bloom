# Issue #1001: stop forming new cratons

## Question

Cratonic crust grows on long runs because `cratons.update` cratonises any continental interior that has been quiet for `CRATON_FORMATION_MYR` (200 Myr). On Earth essentially no craton has formed since about 1 Ga, and cratons are destroyed at roughly 7–25% of their area per Gyr.

The question is whether turning formation off by default (`World.form_new_cratons = False`) is enough to make the cratonic budget net-negative at an Earth-like rate. A second question is whether that reduces the plate slowdown that the growing cratons cause.

## Earth rates used

The issue has the sources.

- **Formation:** about 0 since about 1 Ga.
- **Destruction, long-term:** about 14% of the original craton area per 2 Gyr, about 7% per Gyr (Liu et al. 2018, JGR, doi:10.1029/2017JB015372). This is as reported in search summaries; the full text could not be fetched.
- **Destruction, Phanerozoic:** roughly 25% per Gyr. This is an order-of-magnitude sum of the roots destroyed or thinned in the last about 200 Myr: eastern North China, western Wyoming, southern Kaapvaal, the Angolan Shield, western Man-Léo and Tanzania (Celli et al. 2020, Nat Commun 11:92).

## Model rates before the change

These come from the RFC 0012 G2 baseline at commit `2279d7e`: four seeds, 500 Myr, 1 Myr steps, taken from [`../rfc0012/g2_craton_ledger_divergence.json`](../rfc0012/g2_craton_ledger_divergence.json).

| | Per Gyr |
|---|---:|
| Formed | 26–50% of the starting craton volume |
| Destroyed | 12–20% |
| Net | +8% to +36% |

## Pre-registration (written before any run)

**Arms.** Each seed is generated at node density 1 and run with 1 Myr steps to 500 Myr, with full stepping and the glacier march off (the UI default).

- **on:** the branch with `form_new_cratons = True`. It runs the same code path as `origin/main`, because the flag only gates the formation write.
- **off:** the branch default.

The off arm resumes from the on arm's 190 Myr save. That is valid only while `formed_m3` is still 0 there, and the script refuses the save otherwise. If formation starts earlier, the off arm runs from generation instead.

**Seeds.** 326494039 and 111065103: two G2 seeds, one with mid formation and one with the most.

**Script.** [`bin/debug/craton_budget_trajectory.py`](../../bin/debug/craton_budget_trajectory.py), snapshot every 10 Myr.

**Primary metrics and expectations:**

1. **Cumulative `formed_m3` at 500 Myr.**
   - on: hundreds of M km³, as in G2.
   - off: exactly 0.
2. **Net live-craton change per Gyr, 0→500 Myr.**
   - on: positive, about +5% to +40% per Gyr.
   - off: negative, between −5% and −30% per Gyr, which is Earth's range. If off falls outside that range, the existing sinks are off-Earth and need their own issue. This change does not retune them.
3. **Craton share of genuine continental area.**
   - on: rises after about 200 Myr.
   - off: falls on net over 200→500 Myr.
4. **Plate speed, area-weighted |ω×r|, all nodes and by crust kind, averaged over the 10 Myr samples from 210–500 Myr.**
   - Expected: off ≥ on, of uncertain size; a few to about 20% is plausible.
   - Noise is the spread of the 10 Myr samples within an arm. The comparison counts only if both seeds move the same way by more than about two standard errors of the per-arm mean, treating samples as roughly independent. Otherwise it is inconclusive.

**Controls.** Balance error stays near 0 in both arms. The arms match exactly up to the fork. Peak RSS and runtime are recorded.

## Results

Pending.
