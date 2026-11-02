# Issue 766: RFC 0015 gate G1 on the march

Gate G1 of [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md) for
[#766](https://github.com/adubey/mantle-bloom/issues/766). Tracking:
[#544](https://github.com/adubey/mantle-bloom/issues/544). Feeds the D2 revision,
[#765](https://github.com/adubey/mantle-bloom/issues/765).

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
