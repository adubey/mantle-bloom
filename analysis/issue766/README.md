# Issue 766: RFC 0015 gate G1 on the march

Gate G1 of [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md) for
[#766](https://github.com/adubey/mantle-bloom/issues/766). Tracking:
[#544](https://github.com/adubey/mantle-bloom/issues/544). Feeds the D2 revision,
[#765](https://github.com/adubey/mantle-bloom/issues/765).

## Findings

1. **G1's pass rule, as written:** seed 649658380 passes. Seed 997271774 meets the falsification
   criterion by 0.08 points (20.08% of volume over 0.1 Δt), all from one body on one step: the
   equatorial plateau on the 97.4 Myr drop step, where losing ice grows instead of decaying. The
   ice-age save meets it (57%) because its runaway sheets don't decay.
2. **Response times are short:** volume-weighted median τ of 438 yr and 1.0 kyr on the two saves
   without ice ages. τ/4 is under the prototype's 250 yr floor on every save.
3. **No cap converges within 1% for every body,** even 500 yr → 250 yr. The failures are small
   bodies, 0.7–5.9% of volume at the smallest caps.
4. **Short substeps are cheaper, not dearer.** At a 250 yr longest substep the march settles on
   119 of 120 steps and holds; at 10 kyr it settles on 35 and oscillates through the rest. Median
   step time falls to 0.32–0.51×. The ice tracks the 10 kyr run without ice ages; the ice-age save
   still runs away.
5. **Lake water under closing ice** is likely past the RFC's 1% trigger on two saves (a median of
   20% and 5.8%, an upper bound), so D4 would rebuild the lake hierarchy after the solve.

For #765: a substep limit at the 250 yr floor fits the +15% runtime target with room to spare,
and removes the march's oscillation. It doesn't fix the ice-age runaway, which stays with G0 and
G4.

**Principle:** #3, fast processes in aggregate over long steps. A fast body's state should be set
by its forcing, not by how the step is integrated.

**Owner's decision (2026-11-02):** run G1 on the march, arm (m3), not on D2's end-of-step solve,
which G0 replaced.

## Preregistration

Written before any run in this report. Two smoke runs on seed 649658380 (one step at 6.9 Myr)
shaped the method and are reported below as such. They are the only data seen.

### Tree and saves

G0 prototype `ea9a187` with #575's drifting noise (`a94047d`), plus the probe on gate branch
`claude/issue-766-g1-march`. Arm (m3), 0.25 Myr steps, draw 0 (each world's own noise), 10 Myr from
each diagnosis save:

| Save | Ice ages | Probed steps (closing time) |
|---|---|---|
| Seed 997271774 at 87.4 Myr | off | 89.9, 92.4, 94.9, 97.4 Myr |
| Seed 649658380 at 6.4 Myr | off | 8.9, 11.4, 13.9, 16.4 Myr |
| Seed 374464662 at 132.6 Myr | 8 Myr period | 135.1, 137.6, 140.1, 142.6 Myr |

### What the probe measures, per step

The real march runs as usual. From the same inputs, the probe then re-marches the step:

- **τ, the response time.** The opening ice is scaled by D2's perturbations, ±5% and ±10%. Each
  scaled opening is marched for 10 kyr at fixed 125 yr substeps beside an unperturbed one, and
  each body's τ is when its gap first falls to 1/e. The full march can't resolve τ: in the smoke
  runs the gap fell below 1/e within the first 1–2 kyr substep, and scatter within the step swamped
  its closing ratio. Resolution is 125 yr; a body whose gap never falls to 1/e in 10 kyr has τ
  over 10 kyr.
- **Turnover `V/melt`**, the opening volume over the step's melt rate, as the RFC's `H/|b|`.
- **Convergence,** the RFC's direct check: the same step with the longest substep at 10 kyr (as
  run), 5 kyr, 1 kyr, 500 yr and 250 yr. Each cap's closing volume per body is compared with the
  next smaller cap's.
- **Lake water under closing ice** as a share of all lake water, and **frictional heat's melt
  equivalent** on warm-based ice beside the basal melt.

Bodies are connected components of the opening ice over the hydrology neighbour graph. Those over
10⁴ km² are measured one by one; the rest are pooled and not measured.

### Rules

1. **Body τ** is the median of its four perturbations' e-folding times. Bodies are weighted by
   opening volume, pooled over a save's probed steps.
2. **G1's pass rule, as written:** at least 95% of ice volume has τ ≤ 0.1 Δt (25 kyr). A body with
   τ over 10 kyr counts as failing. Expected: pass on every save, since the smoke runs gave
   125–1,625 yr.
3. **Falsification, as written:** more than 20% of ice volume over 25 kyr.
4. **The substep limit** for #765, from D2's rule of substeps at most τ/4. It is the save's
   volume-weighted 5th-percentile τ over 4, rounded down to 1,000, 500, 250 or 125 yr. Expected:
   250 yr or less on every save. That's below the prototype's 250 yr floor, which would then bind.
5. **Convergence** passes at a cap when every measured body's closing volume is within 1% of the
   next smaller cap's. Expected: fail at 10 kyr (the smoke runs gave gaps of 2.7–16%). Unknown
   below 1 kyr: the smoke runs gave 0.25–3.6%, with the march's settling and hold included.
6. **Cost.** Each save is replayed again for 10 Myr, without probes, at the smaller of rule 4's
   limit and the largest cap that passes rule 5, floored at 250 yr. The median step time is
   compared with the probe run's unprobed steps. The RFC's target is +15% over `origin/main`;
   this ratio is reported beside it and doesn't decide it.

**Reported, not ruled on:** the ice-age save also stands as G0's missing rerun with drifting
noise. Its runaway (600–800 M km³ in #554) is reported, but G1 doesn't close G0.

**The RFC's other G1 reports are not applicable to the march,** because they belong to the endpoint
solve's loops: the sea-level bracket's `|r|`, its transition steps and bodies changing marine
status within it, root-enumeration cost, and membership iterations.

### Smoke runs (method only)

Seed 649658380, one step closing at 6.9 Myr, longest substep 10 kyr:

| Body | Opening (k km³) | Turnover | τ, ±5% and ±10% | Cap 10→5 kyr | 1 kyr→500 yr | 500→250 yr |
|---|---:|---:|---|---:|---:|---:|
| 1 (79° N) | 2,769 | 14.6 kyr | 0.75–1.6 kyr | −2.7% | −0.4% | −0.7% |
| 2 (57° N) | 144 | 332 yr | 125–500 yr | −6.1% | −3.6% | −0.3% |
| 3 (61° N) | 48 | 172 yr | 125–250 yr | −16.4% | +1.6% | −1.6% |

## Results

### What was run

| Run | Save (SHA-256) | Longest substep | Steps |
|---|---|---|---|
| Probes | seed 997271774 at 87.4 Myr (`08ad9461…`) | 10 kyr, as (m3) | 40 × 0.25 Myr, 4 probed |
| Probes | seed 649658380 at 6.4 Myr (`ebbfb5da…`) | 10 kyr | 40, 4 probed |
| Probes | seed 374464662 at 132.6 Myr (`0f8c50c9…`), ice ages at 8 Myr | 10 kyr | 40, 4 probed |
| Cost | the same three | 250 yr | 40 each, no probes |

Tree: gate branch `claude/issue-766-g1-march` at `884fe85` (prototype `ea9a187`, drifting noise
`a94047d`, #741's probes, and `g1_march.py`).

```bash
# From backend/ in the gate-branch worktree:
.venv/bin/python ../bin/debug/replay_g0.py <local-save> --arm m3 --step-years 250000 --checkpoint-myr 0.25 \
    --myr 10 [--ice-age-period-years 8e6] --g1-myr <four closing times> --capture-dir <scratch-output>/<save> \
    --out <scratch-output>/<save>.json
.venv/bin/python ../bin/debug/replay_g0.py <local-save> --arm m3 --step-years 250000 --checkpoint-myr 0.25 \
    --myr 10 [--ice-age-period-years 8e6] --set march_max_substep_years=250 --out <scratch-output>/cost-<save>.json
# From this branch's root:
backend/.venv/bin/python bin/debug/summarize_g1_march.py <scratch-output> --out analysis/issue766/g1.json
```

Each run ran alone and peaked at 1.0–1.6 GB. Each probe took 2–6 minutes. `g1.json` (per body)
and `cost.json` are the compact results. Raw outputs are in B2 as
`issue766-runs-2026-11-02.tar.gz` (SHA-256
`a09d48d52247fb0fb10ddc4a070bc36b69d4e5c1b76d3510e48b6da0a55b8461`).

### τ and the pass rule

| Save | Bodies × steps | τ, volume-weighted p5 / median / p95 | Volume over 0.1 Δt | Rule | Turnover, median |
|---|---:|---|---:|---|---:|
| Seed 997271774 | 25 | 250 yr / 438 yr / over 10 kyr | 20.08% | Falsified | 2.0 kyr |
| Seed 649658380 | 8 | 688 yr / 1.0 kyr / 1.06 kyr | 0% | Pass | 12.2 kyr |
| Seed 374464662 | 25 | 500 yr / over 10 kyr / over 10 kyr | 57.4% | Falsified | 3.7 kyr |

- **Seed 997271774's failure is one body-step.** On the 97.4 Myr step, the plateau body (10.5 M
  km³ at 7.5° N) absorbs added ice (τ = 750 yr for +5% and +10%), but a 5% or 10% loss grows:
  the gap ends the 10 kyr window at 1.17 and 2.60 times its opening value. That is an unstable
  branch, the tipping point #741 found. On the other three steps the same body has τ =
  375–438 yr. That body-step alone is 20.08% of the save's probed volume.
- **On the ice-age save,** the two main sheets (260–510 M km³) don't decay within 10 kyr on several
  steps. They're the runaway (rule note above), not a resolution limit.
- **τ is much shorter than turnover:** 0.4–1 kyr against 2–12 kyr, so the mass-balance feedback,
  not throughput, sets how fast a body answers.

### Convergence

Closing volume per body against the next smaller cap. "Over 1%" counts bodies and their share of
probed volume. A body the larger cap melts away (closing under 1% of its opening) is counted apart.

| Save | 10 kyr → 5 kyr | 5 kyr → 1 kyr | 1 kyr → 500 yr | 500 → 250 yr |
|---|---|---|---|---|
| Seed 997271774 | 11 of 25, 14.3% | 10, 13.9% | 9, 9.2% | 8, 5.9% |
| Seed 649658380 | 1 of 8, 0.5% | 3, 2.4% | 3, 2.4% | 2, 1.9% |
| Seed 374464662 | 6 + 3 vanished, 9.8% | 12 + 1, 10.5% | 8 + 1, 0.5% | 8 + 1, 0.7% |

No cap passes on any save. The largest bodies mostly stay within 1% at every cap; the plateau body
does except on the drop step. The small ones keep moving by 2–5% even at 250 yr. Which state each
settles in is what moves, not the integration.

### Cost

Rule 6 gives 250 yr: rule 4's limit is 125 yr, no cap converges, and the floor is 250 yr.

| Save | Step time, 250 yr / 10 kyr | Substeps per step | Settled steps | End volume (M km³) | Largest step change after 2 Myr |
|---|---|---|---|---|---|
| Seed 997271774 | 15.9 / 30.9 s (0.51×) | 15.5 / 45.5 | 40 / 2 of 40 | 13.36 / 12.90 | 2.9% / 3.1% |
| Seed 649658380 | 13.2 / 41.8 s (0.32×) | 14 / 72.5 | 40 / 0 | 2.11 / 2.11 | 3.2% / 5.9% |
| Seed 374464662 | 22.5 / 65.7 s (0.34×) | 13 / 48.5 | 39 / 33 | 629 / 658 | 667% / 5.1% |

The 10 kyr step times exclude probed steps. Other work on the machine differed between runs, so
substeps per step is the steadier measure, and it agrees. The RFC's +15% is against `origin/main`,
which these runs don't measure; the march at 250 yr costs less than the march as G0 left it.

**Why it's cheaper:** with substeps no longer than τ, the march reaches its settling tolerance in
about 4 kyr and holds the rest of the step. At 10 kyr each substep overshoots, the two-substep
alternation never damps, and it marches the whole step. That oscillation is the within-step
scatter #741 found behind the +8.3% and −10.0% steps.

**The ice-age save** goes from 1.9 to 802 M km³ on its first step at 10 kyr. At 250 yr it reaches
50 M km³ at 2 Myr, then jumps 667% in one step, to the same 600–700 M km³. The limit delays the
runaway; it doesn't stop it.

### Other G1 reports

- **Lake water under visible closing ice,** as a share of lake volume (median / max over the run):
  seed 997271774 20% / 83% (the 83% is the opening, before any march step), seed 649658380 0.1%
  / 1.2%, seed 374464662 5.8% / 22%. It counts water under any ice, so it bounds D4's lag from
  above.
- **Frictional heat on warm-based ice,** as melt, over basal melt: 0.93–0.98 on seed 997271774,
  0.52–0.95 on seed 649658380, 0.72–0.83 on the ice-age save.
- **Not applicable to the march:** the sea-level bracket's `|r|`, transition steps and
  marine-status changes, root-enumeration cost, membership iterations.

### Deviations

- A body whose τ is past the 10 kyr window counts as over 0.1 Δt (25 kyr), as rule 2 said, though
  it may lie between 10 and 25 kyr. All such cases are the drop step's plateau body and the
  ice-age runaway, where the gap grows rather than decays slowly.
- The summary script first divided by the closing volume of bodies the larger cap melted away,
  giving gaps of 10¹⁷%. It now counts those bodies as failing, apart.
