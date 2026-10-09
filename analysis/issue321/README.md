# Issue #321: terrane docking

Replays of two saves with continental terranes on oceanic plates, 200 × 100 ky steps each,
`origin/main` at `10b657b` against this branch. Seed `997271774` (the #310/#320 replay) has
no terranes on oceanic plates at 87.4 Myr, so it can't exercise docking.

```bash
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_terrane_docking.py \
  ~/Downloads/mantle-bloom-seed936513024-90000000y.mbworld --myr 20 --out dock-936.json
```

The compact comparison, with checkpoints every 2 Myr, is [`compare.json`](compare.json).
The branch runs are from before the review fixes. Those fixes only touch partial docks,
docked craton and a few edge cases, and the replays recorded no partial docks and only
0.5 M km³ of docked craton.

## Seed 936513024, 90.0 → 110.0 Myr

9.3 M km² of terrane on oceanic plates at the start; no cratons in this save.

| Measure | main | branch |
|---|---:|---:|
| Terrane area on oceanic plates at 110 Myr | 1.24 M km² | 1.11 M km² |
| Continental area | 440.3 M km² | 442.3 M km² |
| Continental Hc | 14,762 M km³ | 14,817 M km³ |
| Land area | 289.4 M km² | 281.2 M km² |
| Relocations | 0 | 0 |
| Mean step | 10.1 s | 9.1 s |

- Most terrane crust already docked on `main`, through #320's polarized fronts: 2,173 terrane
  cells docked through a frozen front on the branch, against 36 by contact.
- 206 fronts (28 M km³ of Hc) had no continental contact and stayed on their carrier; none
  relocated. No partial docks.
- The `carrier` fallback tier decided 27 fronts that the motion tier decided on `main`.
- Land: the branch ends 8 M km² (2.8%) lower, about #315's run-to-run noise. The gap
  swings: −21 M km² at 96 Myr, +2 M km² at 100 Myr. The 96 Myr dip is not explained.
- The two runs ran side by side, so their timings are comparable only roughly.

## Seed 910211954, 220.2 → 240.4 Myr

44.6 M km² of continental crust on nominally oceanic plates, much of it likely a continent
on a plate labelled oceanic rather than small terranes. 485 M km³ of it is cratonic.

| Measure | main | branch |
|---|---:|---:|
| Terrane area on oceanic plates at 240.4 Myr | 19.93 M km² | 17.80 M km² |
| Continental area | 188.46 M km² | 188.66 M km² |
| Land area | 157.61 M km² | 159.24 M km² |
| Relocations | 159 | 196 |
| Mean step | 18.4 s | 20.8 s |

- 1,872 fronts docked (495 M km³): 1,769 cells through frozen fronts, 397 by contact, of
  which 345 onto another oceanic plate's terrane.
- 4,308 fronts (1,079 M km³) found no continental crust at their contact and took the old
  own-plate path. Docking is not the main fate of this seed's "terranes".
- Docked craton was 0.5 M km³. Terrane craton fell from 485 to 133 M km³, mostly through the
  own-plate path, which reworks it (#332).
- The `carrier` tier decided no fronts: these carriers are mostly continent, which the
  area test excludes.
- Runtime is 13% higher on the branch. Both ran alone. At #320's 1.6 ms per transfer, the
  docking itself costs seconds over the run; the likely cause is the trajectories diverging,
  not tested.

## Ledgers and accounts

- Continental material, craton and mobile cover close to ≤ 10⁻⁴ km³ in all four runs.
- Hm: a terrane's Hm is now booked to the carrier's `oceanic_and_deep_subduction` sink,
  not `suture_hm_subducted_m3`, and its fronts no longer appear in `hm_suture_budget`. On
  seed 936 that moved ~1,100 M km³ between the two accounts over 20 Myr; the total removed
  is unchanged. Suture-Hm totals compared with #319/#320 replays drop by that much for
  accounting reasons alone.

## Follow-ups

- #331: dock terranes by handing their cells to the overrider, which removes partial docking.
- #332: carry craton through continental collision transfer, as terranes now do.
