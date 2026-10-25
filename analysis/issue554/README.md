# Issue 554: RFC 0015 gate G0, what makes the implicit ice solve swing

Gate G0 of [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md) (tracking #544): at
D1's physical mass-balance rates, the end-of-step ice solve swings from step to step and
depends on the step length (#368). G0 tests which of the solve's lagged inputs, or its root
choice, causes that.

**Result: the swing is in the solve itself, as the RFC's outcome rule puts it.**
- No single input fix or root rule passes: 0 of 27 runs of the RFC's arms (9 arms × 3 saves),
  and 0 of 12 diagnostic runs.
- Two of the in-step fixed points don't converge: membership and yield stress (f), and
  deflection (c′).
- What does work is marching each body to equilibrium in short substeps within the step, with
  every fast input recomputed between substeps. That is the RFC's alternative D, run with
  adaptive substeps.
  - It passes on seed 649658380.
  - It passes on seed 997271774 once the climate's per-step noise is frozen; with the noise,
    the swing is 27.4%, against a 25% limit.
  - On the ice-age save, its step gap passes but the ice runs away to about 15× the Last
    Glacial Maximum's excess ice. That points at the ice-age forcing and the sea-level
    feedback, not at the solver.

## Method

The prototype is on gate branch `claude/issue-554-g0-ice-swing` (`ea9a187`). It is not for
merging.
- `ice_sheet.G0` holds the switches, and every switch off is `origin/main`.
- D1 is prototyped as specified, with frictional heat lagged one step and carried on each node:
  - accumulation by mass;
  - degree-day melt of the year's snow, then bare ice (`DDF_snow` 4, `DDF_ice` 8 from #548);
  - the basal energy balance `G + Q_f − k·max(0, T_pmp − T_surf)/H`.
- Physical arms also:
  - use D3's starting abrasion law (`K_G` 1e-4, `l` 0.65);
  - route meltwater by mass;
  - drop the ice-age sea-ice cap, which D6 removes.

`bin/debug/replay_g0.py` (on the gate branch) replays one arm. The gate's metrics come from
[`score_g0.py`](../../bin/debug/score_g0.py):
- the *swing*: the largest step-to-step |ΔV|/V after 2 Myr;
- the *step gap*: |V(0.25 Myr steps) − V(1 Myr steps)|/V at the end;
- on the ice-age save, the second forcing cycle's range of V over the first's, and V's lag
  behind the forcing.

Saves (#548's diagnosis set):
- seed 997271774 at 87.4 Myr and seed 649658380 at 6.4 Myr, each 10 Myr with ice ages off;
- seed 374464662 at 132.6 Myr, 18 Myr with ice ages at an 8 Myr period.

Every arm ran at 0.25 and 1 Myr steps, one replay at a time, with peak RSS around 1–1.3 GB. Since #550
(RFC 0014 D2(b)), `origin/main` has no frozen river or lake ice, so arm (d) is `origin/main`
and every arm includes it.

Scores for every run: [`g0-scores.json`](g0-scores.json).

## The RFC's arms

| Arm | 997271774: swing 0.25 / 1 Myr, gap | 649658380: swing 0.25 / 1 Myr, gap | 374464662 (ice ages): gap, peak V |
|---|---|---|---|
| `origin/main` (model units) | 0.3% / 1.2%, 0.6% | 0.2% / 0.6%, 1.3% | 4%, 1.9 M km³ (doesn't follow the forcing) |
| (a) D1, today's inputs | 86% / 40%, 5.0× | 295% / 19%, 4% | 93%, 820 M km³ |
| (b) sea level held | 419% / 99%, 3.8× | >2,000× / 274%, 100% | 1.7×, 332 M km³ |
| (b′) sea level within the step | 121% / 99%, 1.7× | 753% / 506%, 3.4× | 168×, 997 M km³ |
| (c) deflection held | 63% / 91%, 0.9% | 1,309% / 22%, 59% | 34%, 998 M km³ |
| (c′) deflection within the step | 259% / 107%, 25% | 119% / 367%, 99% | 46%, 970 M km³ |
| (e) seasonal snow fraction | 274% / 73%, 3.5× | 2,599% / 4,235%, 3% | 12.6×, 766 M km³ |
| (f) membership, yield stress at the solved thickness | 150% / 93%, 34% | 41% / >7,000×, 23% | 97%, 712 M km³ |
| (h1) largest stable root | 844% / 742%, 95% | 665% / 977×, 22% | 52%, 1,039 M km³ |
| (h2) smallest stable root | 279% / 54%, 16% | 529× / 2,989%, 115% | 39.5×, 681 M km³ |

A swing of "2,000×" means the ice fell almost to nothing in one step and then regrew.

Arms (g), (i1), (i2), the leave-one-out arms and (j) weren't run. (g) nests three in-step
loops, two of which never converged, so it would take minutes per step and days for the set.
The single-input arms and the in-step iterations already show what the outcome rule asks.

### What the arms show

- **Fast ice on a knife edge.** On seed 997271774, accumulation (~8,000 km³/yr) and surface
  melt (~7,700) nearly cancel. One 0.25 Myr step of the difference is several times the ice
  volume, and the ice's response time (~10 kyr) is far shorter than the step. So each step's
  closing volume is the equilibrium of that step's inputs.
- **One step-lagged inputs that push back.** Several of those inputs come from last step's ice
  and push against it. A thawed bed softens the yield stress, giving a lower, warmer profile.
  Thicker ice deflects the bed further. A bigger sliding body melts more at its base. A
  one-step-lagged negative feedback through a knife-edge equilibrium overshoots: the 0.25 Myr
  series alternate high and low from step to step.
- **Removing one feedback isn't enough.** Diagnostic arms each remove one feedback, and all
  fail on all three saves:
  - (y), yield stress fixed at the hard-bed value, comes closest (seed 997271774: 29% / 707%,
    gap 15%);
  - (q), no frictional heat;
  - (x), no glacial erosion;
  - (n), climate noise frozen.
- **The in-step fixed points inherit the overshoot.** Membership and yield stress (f) and
  deflection (c′) hit their 5-iteration limit on every step. Sea level (b′), solved by
  bracketing, converges in 2 evaluations to ≤ 0.2 m.
- **Multiple stable roots** occur on 10–93% of steps in every arm. Neither fixed root rule
  helps.

## March to equilibrium (alternative D)

The owner's design direction on #554: within each step, march each body forward in short
physical-time substeps. Between substeps, every fast input is recomputed: membership and
yield stress, frictional heat, the snow switch, the mass-balance surface on the current ice's
load-equilibrium deflection, and sea level. Once the ice settles, hold it for the rest of the
step. Only the slow forcing is held for the whole step.

| Arm | What it does |
|---|---|
| (m) | 1 kyr substeps; settled when the means of successive substep pairs differ by < 0.1%; at most 100 kyr marched, then held |
| (m2) | (m) with the lagged inputs damped (mean of the last two substeps), adaptive substeps, no march cap, settled at < 0.02% |
| (m3) | (m2) settled at 0.1%, halving the substep only for alternations > 1%, holding the mean rates of the last six substeps |

| Arm | 997271774: swing, gap | 649658380: swing, gap | 374464662: gap, cycle growth, lags, V |
|---|---|---|---|
| (m) | 8.5% / 8.2%, 5.6% | 6.8% / 4.3%, 1.0% | 99% (0.25 Myr runs away to 772 M km³; 1 Myr held to 4–13 by the cap) |
| (m3) | 27.4% / 10%, 3.2% | 9.4% / 6.5%, 0.3% | 9.2%, 3.07, −0.25 and 0.5 Myr, 597–796 M km³ |
| (m3), climate noise frozen | 14.9% / 2.4%, 2.9% | | |

(m2) was too slow to run in full, at 155–191 s per step, and gave the same volumes as (m).

- **The overshoot is gone.** The marched series are smooth. Seed 997271774 settles near
  13 M km³ (2.3× its model-unit ice) and seed 649658380 near 2.3 M km³, at both step sizes.
- **The ledger closes where the march runs the whole step.** On seeds 997271774 and 649658380,
  (m3) rarely settles at 10 kyr substeps, so it marches the whole step. Its residual there is at
  most 4 × 10⁻⁵ of a step's accumulation, inside the RFC's 0.1%. Where a step settles and the
  alternating state is held, the hold's booking leaves more: up to 0.6% for (m3) on the
  ice-age save, and up to 1.1% for (m). The hold needs a booking that nets the alternation to
  zero, or the march has to run to the end of the step.
- **A small alternation remains.** Between substeps it is 0.1–1% at 1–10 kyr substeps (±0.4% at
  500 yr), probably thin margins responding within centuries. Damping didn't remove it.
- **The rest of seed 997271774's swing is the climate.** Its one 27% drop at ~95.6 Myr happens
  on the first substep of a step: the march then settles at the new level. The model redraws
  its climate noise (humidity, moisture-flux and current textures, from `(seed,
  elapsed_years)`) every step, and fast ice now follows it. With the noise frozen at the
  opening draw, the largest swing is 14.9%. The redraw is a climate design question for
  RFC 0015 to raise, not a solver one.
- **Cost.** (m3) took 40 s per step at 0.25 Myr steps and 160–200 s at 1 Myr steps on seed
  997271774, against ~10 s for (a). Substeps halve whenever the alternation exceeds 1%.

## The ice-age save runs away

Under the march, both step sizes take seed 374464662's ice from 1.9 to 600–800 M km³ in the
first step (the save opens near the glacial maximum, intensity 0.95). Sea level falls to
−2,800 m, and the ice stays above ~600 M km³ through the interglacials. At 0.25 Myr steps it
collapses twice, to 51 and 31 M km³ within half a million years, and regrows. That pattern
looks like the marine-grounding feedback running in both directions.

The march only lets the ice reach the equilibrium D1's rates imply, so this belongs to the
forcing and to sea level (G4), not to G0:
- the model's uniform 8 °C ice-age cooling;
- sea level free to fall by kilometres;
- no ice shelves or marine retreat.

For scale, the Last Glacial Maximum's excess was ~52 M km³ and ~134 m.

## Other findings

- **Sliding speeds aren't physical.** At 0.25 Myr steps on seed 997271774:
  - 18–22% of the ice area slides faster than 100 m/yr, 3.7% faster than 1 km/yr, and 0.1–0.2%
    faster than 10 km/yr;
  - the fastest 1% of nodes are ~200 m thick, slide ~7 km/yr and each carries ~150 km³/yr of
    balance flux.

  These speeds come from the balance-flux velocity at thin margins, the 10 m thickness floor in
  `_diagnose_flow`. They feed frictional heat and abrasion, but (q) and (x) show they don't
  drive the swing.
- **`origin/main` drifts.** Its sea level falls 117 m over 10 Myr on seed 997271774 (−408 to
  −525 m) while its ice barely changes. The march adds ~20 m for its extra ~7 M km³. That drift
  isn't from the ice.
- **Erosion during the hold is deferred.** At its equilibrium flow the ice erodes the bed for
  the whole step. If that changes the bed enough to move the ice's equilibrium, the hold
  misses a slow feedback (#554).

## Outcome and next steps

By the RFC's outcome rule, "if (g), (i1) and (i2) all fail, the instability is in the solve
itself: D2 substeps all bodies (alternative D), and the RFC is revised." Every single-input arm
fails, and the in-step loops (g) would combine don't converge. The march does what the rule
calls for, and it passes where the climate isn't redrawn under it.

For the RFC revision (`rfc-authoring`):
1. **D2 becomes a march.** Each body evolves in physical-time substeps, with all fast inputs
   recomputed, until it settles or the step ends. The owner's slice-by-slice flow version
   would take velocity from shallow-ice flow; that replaces the plastic equilibrium shape and
   fixes the sliding speeds.
2. **Climate noise per step.** A climate that redraws its noise every step moves fast ice with
   it. Either the noise persists over longer times, or the RFC accepts it as forcing.
3. **G4 and the ice-age forcing.** Physical rates under the uniform 8 °C cooling glaciate to
   ~15× the LGM.
4. **G1** (response times, substep costs) would measure the march directly.

## Reproduce

```bash
git switch claude/issue-554-g0-ice-swing   # in its own worktree
cd backend
.venv/bin/python ../bin/debug/replay_g0.py <save> --arm <arm> --myr 10 --step-years 250000 \
    --checkpoint-myr 0.25 --out <scratch>/<arm>-<seed>-250000.json
# The ice-age save: --myr 18 --ice-age-period-years 8e6. The noise freeze: --freeze-climate-noise.
python ../bin/debug/score_g0.py --fine <...-250000.json> --coarse <...-1000000.json>
```
