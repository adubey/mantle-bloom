# Issue 912: RFC 0015 gate G4, sea level

Gate G4 of [RFC 0015](../../rfc/0015-frozen-water-in-physical-units.md) for
[#912](https://github.com/adubey/mantle-bloom/issues/912). Tracking:
[#544](https://github.com/adubey/mantle-bloom/issues/544).

**Principle:** #1, conservation. Sea level comes from a conserved water budget, so every cubic
metre booked as ice must be water that has left the ocean. Today's debit books grounded marine
ice in full, including the part that sits in ocean capacity the solve also counts.

**Result: adopt D6's simple form.** On both deciding saves, the simple form and the generalised
correction never differ by more than 0.01 m of sea level at 0.25 Myr steps (0.016 m at 1 Myr),
against a threshold of 1 m. Today's debit and D6's differ by at most 0.08 m on those saves. On
the ice-age save, which is reported and doesn't decide, they differ by up to 2.4 m, and the
simple and generalised forms by up to 0.26 m.

## Method

The method was fixed in [#912's first comment](https://github.com/adubey/mantle-bloom/issues/912#issuecomment-6902313139)
before any run. In short, three sea levels are solved after each step on the same final state
(elevations `z`, ice `H`, accounting areas `A`, budget `V_total`, and other trapped water `L`
for lakes, lake ice and river water), with `w = 0.917 H`:

| Arm | Equation solved for `h` |
|---|---|
| **T, today** | `C(h) = V_total − L − Σ A·w` |
| **S, D6 simple** | `C(h) = V_total − L − Σ A·max(0, w − D(h_ref))`, with `h_ref` the arm's own level from the step before |
| **G, generalised** | `Σ_mask(h) A·max(0, (h − z) − w) = V_total − L − Σ A·w` |

`C(h)` is eustasy's capacity over the seeded connected-ocean mask, and `D(h) = max(0, h − z)` on
that mask (0 elsewhere). Arm G is the liquid-water balance: Goelzer et al. (2020) Eqs. 5, 13 and
14 (with the 2022 corrigendum's outer `max`), with the reference level set to the level being
solved and the bed to the current one. Bed and sea-level changes from any cause therefore enter
without a lag or a reference state. Arm S lags `D` by one step because that is the level
`update_sea_level` knows when it runs.

The replay runs today's rule. S and G are diagnostics, not fed back. Ice is classified at the
replay's level as land (outside the connected ocean), grounded marine (`w > D`) or floating
(`w ≤ D`).

**Gate statistic:** the largest `|h_S − h_G|` over the 10 steps on each deciding save. The
threshold is 1 m, on both deciding saves at 0.25 Myr.

## Runs

```bash
OUT=~/mantle-bloom-captures/issue-912 bin/debug/run_912_g4.sh
# each run is:
PYTHONPATH=backend backend/.venv/bin/python bin/debug/measure_flotation_debit.py \
    <save> --steps 10 --myr {0.25|1} [--ice-age-period-years 8e6] --out <captures>/<seed>-<myr>.json
```

- **Commit:** `0d4505ee` (`origin/main`, no model changes). The march is the default,
  `trend-1k` (#906).
- **Step size:** 0.25 Myr is primary and 1 Myr second, 10 steps each.
- **Saves:** the [#548](../issue548/README.md) diagnosis saves, with SHA-256 checked before each
  run:

| Seed | Age | SHA-256 | Ice ages | Role |
|---|---:|---|---|---|
| 997271774 | 87.4 Myr | `08ad9461…ed762b` | off | deciding |
| 649658380 | 6.4 Myr | `ebbfb5da…73f7b3` | off | deciding |
| 374464662 | 132.6 Myr | `0f8c50c9…a4279c` | 8 Myr period | report only (#907) |

- **Captures:** full outputs and logs are in `~/mantle-bloom-captures/issue-912/`, with
  `SHA256SUMS`. [`g4.json`](g4.json) holds the compact per-step results and the capture
  checksums.
- **Resources:** peak RSS was 0.95–1.53 GB per run. The first run (649658380 at 0.25 Myr)
  overlapped another agent's run, with the owner's override. The rest ran alone.
- **Check:** arm T re-solves the replay's own level to within 0.0001 m on every step, so the
  differences below are not solver error.

## Results

`T − S` is today's level minus D6's. It is negative because today's debit over-counts
grounded marine ice and so puts sea level too low. `S − G` is the generalised correction.

### Gate statistic

| Save | Step | max \|S − G\| (m) | max \|T − S\| (m) | End: T − S (m) | End: S − G (m) | Sea level, open → end (m) |
|---|---:|---:|---:|---:|---:|---:|
| 997271774 | 0.25 | **0.0095** | 0.063 | −0.0434 | 0.0014 | −408.2 → −473.0 |
| 997271774 | 1 | 0.0164 | 0.080 | −0.0115 | 0.0006 | −408.2 → −536.1 |
| 649658380 | 0.25 | **0.0001** | 0.0043 | −0.0001 | 0.0000 | −100.5 → −111.4 |
| 649658380 | 1 | 0.0003 | 0.0033 | −0.0024 | −0.0001 | −100.5 → −104.5 |
| 374464662 (ice ages) | 0.25 | 0.259 | 2.385 | −0.0724 | −0.0019 | −355.3 → −338.8 |
| 374464662 (ice ages) | 1 | 0.054 | 1.661 | −0.0952 | −0.0070 | −355.3 → −320.5 |

### Ice by class

Ice volume in M km³ of ice, and area in M km², at the opening state, at the step with the
most grounded marine ice, and at the end.

| Save | Step | State | Land | Grounded marine | Floating |
|---|---:|---|---:|---:|---:|
| 997271774 | — | open | 5.65 / 40.49 | 0 / 0 | 0.0000 / 0.005 |
| 997271774 | 0.25 | step 1 (most marine) | 15.62 / 14.92 | 0.095 / 0.167 | 0 / 0 |
| 997271774 | 0.25 | end | 14.29 / 14.67 | 0.045 / 0.102 | 0.0026 / 0.007 |
| 997271774 | 1 | step 1 (most marine) | 14.22 / 14.78 | 0.088 / 0.156 | 0.0001 / 0.004 |
| 997271774 | 1 | end | 12.26 / 14.62 | 0.017 / 0.040 | 0.0001 / 0.003 |
| 649658380 | — | open | 3.56 / 21.40 | 0.0069 / 0.028 | 0.0001 / 0.016 |
| 649658380 | 0.25 | step 3 (most marine) | 2.53 / 3.45 | 0.0050 / 0.012 | 0.0001 / 0.008 |
| 649658380 | 0.25 | end | 2.45 / 3.38 | 0 / 0 | 0.0001 / 0.008 |
| 649658380 | 1 | step 2 (most marine) | 2.32 / 3.46 | 0.0082 / 0.021 | 0 / 0 |
| 649658380 | 1 | end | 2.18 / 3.17 | 0.0022 / 0.008 | 0.0003 / 0.013 |
| 374464662 | — | open | 1.53 / 26.02 | 0.33 / 0.38 | 0 / 0 |
| 374464662 | 0.25 | step 0 (most marine) | 12.06 / 17.77 | 2.96 / 3.26 | 0.0002 / 0.012 |
| 374464662 | 0.25 | end | 4.73 / 8.79 | 0.11 / 0.17 | 0.0001 / 0.004 |
| 374464662 | 1 | step 6 (most marine) | 12.11 / 18.45 | 2.29 / 2.90 | 0.0096 / 0.055 |
| 374464662 | 1 | end | 5.82 / 10.63 | 0.14 / 0.23 | 0.0002 / 0.017 |

Every step's classes, debits and levels are in [`g4.json`](g4.json).

## Findings

1. **Gate outcome: adopt D6's simple form.** The largest `|S − G|` on a deciding save is
   0.0095 m at 0.25 Myr and 0.016 m at 1 Myr, two orders of magnitude under 1 m. The step-size
   dependence the method set aside for an owner decision doesn't arise.
2. **D6's change is small without ice ages.** Grounded marine ice is under 1% of ice volume on
   both deciding saves (at most 0.095 M km³ against 14–16 M km³ on land), so today's
   over-debit moves sea level by at most 0.08 m. Floating ice is negligible: at most
   0.0026 M km³.
3. **The simple form's residual comes from lagging `D` by a step.** `|S − G|` is largest on
   the step where sea level moves most: on 997271774 that is step 0, which falls 30 m. On the
   ice-age save it peaks at 0.26 m when grounded marine ice reaches 3 M km³ over 3.3 M km². If
   D6 computes `D` at each candidate level inside the bisection, the debit is arm G exactly and
   this residual goes. That costs one extra `max` per node per candidate, and the RFC's design
   doesn't depend on it.
4. **Under ice ages the over-debit reaches metres (report only).** At the cycle's peak, today's
   rule puts sea level 2.4 m too low (0.25 Myr) or 1.7 m too low (1 Myr). Ice ages are out of
   RFC 0015's scope (#907).
5. **Density is the larger unmodelled term (report only, not in the rule).** The model has one
   water density, so Goelzer's Eq. 10 (ocean water at 1028 kg m⁻³ against fresh meltwater) is
   zero in both arms. Evaluated separately over the ocean area, it is 0.16–0.22 m on 649658380
   and 0.86–1.10 m on 997271774. That is larger than the flotation correction on these saves.
   It concerns the model's single water density, not D6's form, so it is recorded here and not
   carried into D6.
6. **The opening states aren't comparable with the stepped states.** The saves were captured
   before D1–D2, and the first step re-grows the ice: 997271774 goes from 5.65 M km³ over
   40.5 M km² of thin cover to 15.8 M km³ over 14.8 M km². The per-step results, not the
   opening state, are the gate's evidence.
