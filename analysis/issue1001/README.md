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

Runs were on 2026-11-21 at branch commit `d2985549`, two seeds, 500 Myr at 1 Myr steps, 10 Myr snapshots. The numbers are in [`comparison.json`](comparison.json).

**Controls.** Both 190 Myr saves had `formed_m3 = 0`, and the off arm's first snapshot reproduced the on arm's 190 Myr volume and speed exactly. The worst |balance error| was 1e-12 M km³. Peak RSS was 0.62 GB, and each arm took 19–31 min.

### Pre-registered metrics

| Metric | Seed 326494039, on | Seed 326494039, off | Seed 111065103, on | Seed 111065103, off | Expected | Met? |
|---|---:|---:|---:|---:|---|---|
| Formed by 500 Myr (M km³) | 1,114 | **0** | 782 | **0** | on: hundreds; off: 0 | yes |
| Destroyed by 500 Myr (M km³) | 347 | 305 | 244 | 232 | | |
| Net craton volume, % per Gyr | +33.8 | **−13.5** | +30.6 | **−13.2** | on: +5 to +40; off: −5 to −30 | yes |
| Net craton **area**, % per Gyr | +53.2 | **+0.8** | +40.3 | **−0.03** | (not pre-registered) | |
| Share of genuine continental area, 200 → 500 Myr | 0.542 → 0.601 | 0.446 → **0.477** | 0.585 → 0.622 | 0.505 → **0.517** | on rises; off falls | on yes; **off no** |
| Speed, all nodes, 210–500 Myr mean ± SE (cm/yr) | 0.923 ± 0.030 | 0.945 ± 0.029 | 0.676 ± 0.017 | 0.674 ± 0.016 | off ≥ on, by > 2 SE on both seeds | **inconclusive** |
| Off − on speed: all / continental / oceanic | | +2.3% / +1.0% / +3.3% (z 0.5 / 0.2 / 0.8) | | −0.3% / +1.6% / −1.1% (z −0.1 / 0.2 / −0.4) | | |

### Reading

1. **Formation off works as intended.** No craton forms, and the cratonic volume falls by about 13% per Gyr on both seeds, inside the pre-registered Earth range.

2. **That loss is thinning, not area loss.** The off arms' craton area stays flat to within 1% over 500 Myr. The mean cratonic thickness falls from 45.0 to 41.8 km and from 34.5 to 32.3 km. Of the volume lost, erosion books 54% and 39% on the two seeds, collision rework 16% and 41%, and rifting 29% and 20%. A column leaves the craton only when its `craton_crust_m` reaches exactly 0, and none of the sinks take it there.

   Earth's 7–25% per Gyr is area: keels removed, leaving ordinary continental crust. Comparing the model's volume sinks against it, as the issue and this pre-registration did, compared unlike quantities. By area, the model has no Earth-like destruction.

   The off-arm share still rises, from 0.446 to 0.477 and from 0.505 to 0.517. That is because genuine continental area shrinks under a constant craton area: 227 → 212 and 201 → 197 M km².

3. **New cratons don't explain the slowdown over 500 Myr.** By 500 Myr the on arms carry 27% and 20% more craton area than the off arms. Even so, the 210–500 Myr speeds agree within about 1 SE on both seeds, continental and oceanic nodes alike. This was with the craton brake at its default, 1×.

   Both arms slow down together after the fork. Seed 326494039 goes from 1.35 cm/yr at 200 Myr to 0.75–0.83 at 500 Myr; seed 111065103 from 0.94 to 0.64–0.66. Something other than newly formed craton drives that.

### Exploratory: do the existing cratons slow the plates? (not pre-registered)

Reading 3 left open whether cratons matter for speed at all, so one extra arm was run after the results above were in. It used seed 326494039, the seed that slows most. The arm resumes the same 190 Myr save with `form_new_cratons = False` and sets `craton_crust_m = 0` on every node. That removes every craton, and with them the craton brake and the jam, then the arm steps to 500 Myr. The ledger is not rebalanced, so its balance error means nothing in this arm.

The setup is the same as the off arm, plus:

```python
for plate in world.plates:
    if plate.node_count():
        plate.set_fields_on_plate(craton_crust_m=np.zeros(plate.node_count()))
```

| Seed 326494039, 210–500 Myr | off (cratons kept) | no craton | Change |
|---|---:|---:|---:|
| Speed, all nodes (cm/yr) | 0.945 ± 0.029 | 1.118 ± 0.022 | +18% (z 4.8) |
| Continental | 0.793 ± 0.032 | 0.955 ± 0.024 | +20% (z 4.1) |
| Oceanic | 1.115 ± 0.028 | 1.286 ± 0.024 | +15% (z 4.6) |

On this seed the existing cratons cost about 18% of plate speed, a clear effect. It comes from the cratons the world started with. Craton formed after 200 Myr added nothing measurable (reading 3), probably because a new craton forms at least 500 km inside its plate, away from the fronts that jam.

The no-craton arm still slows from 1.29 cm/yr at 190 Myr to about 0.95 at 490 Myr, so most of the secular slowdown has another cause. Over the same span the plate count rises from 31 to 41. That cause belongs in its own issue. It is one seed and one exploratory arm.

Per the pre-registration, a destruction rate outside Earth's range means the sinks need their own issue, and this change does not retune them. The follow-up is a decratonization process that removes craton area, not crust.

