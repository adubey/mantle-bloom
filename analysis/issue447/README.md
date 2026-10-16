# Issue #447: boundary physics from each node's crust (RFC 0009)

## Summary

On 18 fresh-world runs (seeds 7, 11 and 2024; control, mixed A and mixed B; 500 Myr), deciding
the 22 boundary sites from each node's crust instead of its plate's majority type moves the
model a lot.

- **Continents last longer.** Genuine continental area at 500 Myr in the control runs is
  30.5 / 20.9 / 22.9% of the sphere, against 17.1 / 13.8 / 16.5% on origin/main. Less
  continent goes down trenches (`deeply_subducted_m3` falls 14–33% in the controls), and
  juvenile additions stay within −24% to +2%.
- **Fewer plates, merges and splits.** Active plates at 500 Myr in the controls fall from
  47 / 45 / 43 to 36 / 37 / 35, and splits from 74 / 60 / 56 to 49 / 47 / 46.
- **Ledgers balance as before.** The worst relative balance error over every checkpoint is
  1.1e-3 on the branch and 1.3e-3 on origin/main (the continental-material drift G2 traced to
  #299). The craton ledger closes.
- **G2 rerun.** The 100→500 Myr rule passes on every seed. The juvenile-additions gap between
  the control and the mixed runs narrows on all three seeds (from 782 / 1,209 / 959 to
  462 / 631 / 92 M km³), and so does the deep-subduction gap. The 500 Myr continental-area gap
  narrows on seed 11 (4.2 to 3.4 points) and stays inside the noise on seeds 7 and 2024. G2 as
  written still fails on seed 11's area. It now also fails on splits and on overlap at 200 Myr
  (the control has fewer splits and more overlap than the mixed runs), and still on plate
  count at 500 Myr.
- **No runtime cost.** Rerun alone, the branch matches origin/main step for step. The full
  branch set took longer only because something else was using the machine (see Runtime).

## What was run

- **Probe:** `bin/debug/probe_mixed_plates.py` from PR #446 (G2), unchanged, copied into each
  tree's `bin/debug/` for the run so it loads that tree's backend. Defaults: seeds 7, 11 and
  2024; control, mixed A and mixed B; 5 Myr steps to 500 Myr; checkpoints every 100 Myr. See
  [analysis/rfc0009](../rfc0009/README.md) for what the mixed runs are.
- **Baseline:** origin/main `788e654`, a detached worktree.
- **Branch:** `9ad2aca` (`claude/issue-447-node-crust-rules`; the analysis commits after it
  change no model code).
- The two sets ran one after the other, nothing else of mine running alongside.

```bash
# in each tree, with probe_mixed_plates.py copied into bin/debug/
PYTHONPATH=backend /usr/bin/time -l backend/.venv/bin/python bin/debug/probe_mixed_plates.py --out $S/g2_<tree>.json
python3 bin/debug/compare_g2.py $S/g2_main.json $S/g2_branch.json --out analysis/issue447/g2_ab.json
```

Compact results are in [g2_ab.json](g2_ab.json). A 100 Myr rerun of seed 2024 mixed B on the
branch reproduced the set's numbers exactly.

## Results at 500 Myr

Per seed: control / mixed A / mixed B. Areas are percent of the sphere; ledger accounts are
cumulative, M km³.

| Measure | Seed 7, main | Seed 7, branch | Seed 11, main | Seed 11, branch | Seed 2024, main | Seed 2024, branch |
|---|---|---|---|---|---|---|
| Genuine continental area | 17.1 / 21.2 / 17.6 | 30.5 / 25.7 / 29.5 | 13.8 / 19.0 / 17.0 | 20.9 / 24.3 / 24.2 | 16.5 / 18.6 / 16.9 | 22.9 / 20.8 / 22.7 |
| Change 100→500 Myr | −20.8 / −20.5 / −22.8 | −11.5 / −16.8 / −13.4 | −10.1 / −9.5 / −11.9 | −5.5 / −3.6 / −4.3 | −7.4 / −6.4 / −9.1 | −0.7 / −2.1 / −1.3 |
| Juvenile additions | 2,512 / 1,638 / 1,823 | 2,550 / 2,036 / 2,140 | 3,437 / 2,085 / 2,371 | 2,609 / 2,150 / 1,807 | 2,219 / 1,134 / 1,387 | 1,824 / 1,541 / 1,922 |
| Deeply subducted | 4,291 / 2,927 / 2,919 | 3,683 / 3,180 / 3,148 | 5,449 / 3,687 / 3,697 | 3,643 / 2,304 / 2,105 | 2,815 / 1,410 / 1,988 | 2,106 / 2,064 / 1,999 |
| Active plates | 47 / 50 / 49 | 36 / 40 / 40 | 45 / 48 / 50 | 37 / 45 / 45 | 43 / 48 / 51 | 35 / 48 / 41 |
| Merges | 40 / 33 / 38 | 27 / 44 / 43 | 13 / 13 / 12 | 9 / 16 / 25 | 13 / 13 / 6 | 16 / 21 / 9 |
| Splits | 74 / 62 / 63 | 49 / 66 / 62 | 60 / 54 / 58 | 47 / 55 / 55 | 56 / 53 / 44 | 46 / 57 / 44 |
| Plates mixed (5–95%) | 0.34 / 0.40 / 0.41 | 0.47 / 0.47 / 0.65 | 0.22 / 0.44 / 0.40 | 0.38 / 0.49 / 0.47 | 0.35 / 0.29 / 0.35 | 0.49 / 0.40 / 0.56 |

Across 0→500 Myr the branch's control loses continent more slowly from the first checkpoint
on (seed 7: 51.0 → 42.0% at 100 Myr against 37.9% on main), so the change isn't one event.
These runs don't tell which of the 22 sites carries it. Two that act on the continent's
edge directly are the likeliest: ocean floor on a continental plate now subducts and retreats
on its own, so the continent stops being consumed with it, and a continent behind its own
slab is conserved as a terrane rather than subducted. That is a hypothesis, not a measurement.

### G2 rerun

The gap is mean(mixed) − control. The noise floor is max(|mixed A − mixed B|, floor), with a
1-point floor for area and none for the ledger accounts, and a gap fails beyond 2 × noise.

| Gap (noise) | Seed 7, main → branch | Seed 11 | Seed 2024 |
|---|---|---|---|
| Genuine area at 500 Myr, points | +2.3 (3.6) → −2.9 (3.7) | **+4.2 (2.0) → +3.4 (1.0)** | +1.3 (1.7) → −1.2 (1.9) |
| Area change 100→500 Myr, points | −0.9 (2.3) → −3.6 (3.5) | −0.5 (2.4) → +1.6 (1.0) | −0.4 (2.7) → −1.0 (1.0) |
| Juvenile additions, M km³ | −782 (185) → −462 (104) | −1,209 (286) → −631 (343) | −959 (253) → −92 (381) |
| Deeply subducted, M km³ | −1,368 (8) → −519 (32) | −1,757 (10) → −1,438 (200) | −1,116 (578) → −74 (65) |

Bold fails G2's rule. The metrics G2 judges on two or more seeds:

| | origin/main `788e654` | Branch |
|---|---|---|
| Failing G2 metrics | active plates at 500 Myr; overlap at 100 Myr | active plates at 500 Myr; overlap at 200 Myr; splits |

On origin/main the area gap now fails on seed 11 only (G2 measured it at `652c2c0`, before
#435, failing on 11 and 2024). The branch's new failures run the same direction as its
area gap: the controls split less and overlap more than the mixed runs. The owner's G2
decision already reads the controls and the mixed runs as different starting layouts, so
these don't by themselves point at a defect in the 22 sites.

## Runtime

| | origin/main | Branch |
|---|---|---|
| Wall time, 9 runs | 3,206 s | 7,130 s |
| Instructions retired | 57.0 × 10¹² | 54.3 × 10¹² |
| Cycles | 13.42 × 10¹² | 13.48 × 10¹² |
| Peak memory footprint | 1.07 GB | 0.76 GB |

The branch did the same CPU work, but its wall time doubled because something else was using
the machine during that set: seed 2024's mixed runs ran 4× slower than origin/main from the
first 100 Myr. Rerun alone under cProfile, the first 100 Myr of seed 2024 mixed B took 52.8 s
on the branch and 51.0 s on origin/main, with the same profile.
