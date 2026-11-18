# RFC 0019 G3: three-round bounded tuning plan

**Status: three tuning rounds completed 2026-11-18; G3 failed and the owner authorized an exception to continue RFC 0019.** The
previously opened held-out worlds are spent and remain acceptance-only: they will not be
used to choose a parameter or evaluate a successor. This plan allows at most three further
tuning rounds. A round is one preregistered parameter question and its paired arms; repeated
reports of the same arm do not create another round. Stop early if a round rules out the
candidate family. No candidate passed after the three rounds. The owner explicitly authorized
RFC continuation with G3 recorded as failed and the inland-rain mismatch treated as a known
modeling simplification; this does not count as passing G3 or selecting a candidate for
production.

## Decision basis

The desired output is a visually plausible rainfall pattern, especially rain belts, rain
shadows and dry interiors. Cold-land precipitation remains a formal positive-supply check;
cold-inland amount is a reported tradeoff rather than a tuning target. Water closure,
nonnegative finite fields, density convergence and the existing runtime limit remain hard
constraints. The governing principle is **conservation**: rainfall may be redistributed by
the sink field, but it must not create water. Pattern decisions also respect **scale
separation**: changes to coarse regional delivery must be assessed before fine placement.

The current frozen `p=0.25` candidate is closed, density-convergent, and passed the approved
rain-belt, ridge-direction and cold-coast checks on the two original held-outs. It failed the
dry-interior check in 5/16 held-out combinations (all CFD), with ratios 0.613–0.835 against
the 0.6 limit. The earlier tuning screen shows that the finite-volume (FV) solver has lower
worst-case CFD interior/coast ratios than wind arrows on A/D/E; both pass ridge-direction
checks there. Use FV as the single candidate solver in this plan, with diagnostic and CFD
winds retained as separate modes. This is a tuning scope choice, not a claim that solver
selection alone fixes the known failure.

The proposed change follows the **conservation** principle by moving a fixed source through
a different positive rainout-rate field in the stationary solve. Its testable invariant is
that total placed precipitation equals total source within `1e-8`, with nonnegative finite
column water and rain. It follows **scale separation** by changing only the coarse sink
distribution, with transport, source, cap, ridge interception and native placement fixed.

## New acceptance-only registry

The original two held-outs have been opened. Two replacement acceptance worlds are registered
here before any of their climate fields are inspected. Selection used only the private bucket
metadata (save name, age and SHA-1), excluding all RFC 0019 tuning saves and the two spent
held-outs. The listed exact saves are not reused as tuning cases.

| Save | Simulated age | Bucket SHA-1 |
|---|---:|---|
| `mantle-bloom-seed66121723-26800000y.bloom` | 26.8 Myr | `c434815348a65596b3234987b4d279d294b063d5` |
| `mantle-bloom-seed544349906-422900000y.bloom` | 422.9 Myr | `584ca71983ac95555eaf7544a774af07aa52fc3f` |

The initial candidate `seed447276134-244400000y` was replaced after peer review found its
geography in a prior basin survey. `seed66121723-26800000y` and `seed544349906-422900000y`
have no prior record in this repository; they were chosen by age stratum and bucket hash only.
Reverify each digest against the loaded bytes before the final acceptance evaluation. If a
registered save is unavailable or fails an input-integrity check, report that and stop; do
not substitute another save after seeing candidate results. Use the existing predeclared G3
geographic eligibility masks and all approved thresholds unchanged. New acceptance results
are a final holdout evaluation, not a fourth tuning round. If a candidate fails, do not revise
it using these results.

## Tuning cases and common controls

Use only existing registered tuning captures A, D, E and the cold source, at density 4 first,
with diagnostic and CFD winds, using FV. A, D and E are three checkpoints from one seed, so
report each checkpoint but treat them as one independent world family; the cold save is the
second family. If a round selects a plausible setting, confirm its density-1 behavior and
cold-coast supply on these same tuning captures before freeze.

Keep the `p=0.25` ridge interception, source arrays, transport faces, column cap, rain
placement, masks and `9.1 day` baseline rainout constant fixed in every sink-profile arm.
For each sink profile, pair its `p=0.25` arm with a matching `p=0` arm solely to measure the
ridge-direction and contrast-gain checks. The alpha-zero `p=0.25` arm must reproduce the
existing frozen `ridge_025` result to `rtol=1e-10`. Report absolute near-coast and inland
rain as well as ratios. To prevent ratio-only wins, a tuning arm is ineligible if it reduces
absolute `>1,500 km` inland rain below `80%` of the frozen `p=0.25, alpha=0` candidate or
changes its land rain share by more than `2 percentage points`. These are **tuning safeguards**, not
new formal G3 acceptance limits. Save common-scale maps for a visual check that a numeric
shadow improvement has not erased recognizable rain-shadow structure.

For choosing among eligible arms, apply conservation, finite/nonnegative, cap and placement
checks first. An arm that fails a hard check or either tuning safeguard is excluded from the
ranking. A pattern-failing but otherwise eligible arm may advance to the next declared tuning
round; it is not freeze-eligible and cannot be reported as passing G3. Then minimize the
worst normalized shortfall from the approved rain-belt (`1.5`),
ridge-direction (`60%` plus positive median contrast gain), and inland/coast (`0.6`) limits
across each world family, checkpoint and wind mode. Cold-coast precipitation must be positive
on eligible cold land. Break ties by lower absolute change in land rain share, then stronger
family-balanced median ridge contrast gain (first median within the seed family, then across
families). Legacy RMS, folded bands, Köppen shares, cold-inland rain and
land-rain share versus the current climate remain diagnostics, not score terms. A case is not
eligible for a metric only under the geography rule already registered in
`g3-preregistration.md`; eligibility cannot be chosen from rainfall output.

## The three rounds

### Round 1 outcome

Round 1 completed on all four registered tuning captures, both wind modes, at density 4.
All arms passed the hard budget/placement checks. For eligible arms the belt, ridge-shadow,
and cold-coast pattern checks passed; the inland/coast ratio was the only pattern shortfall.

The coast-focused `alpha=1, L=750 km` arm is eligible to advance: its worst normalized
shortfall was `1.5841`, versus `2.0451` for the fresh same-run `ridge_025` baseline and
`1.7023` for its equal-integrated uniform-land control. Its worst inland/coast ratio was
`1.550` on the cold capture with CFD winds, against the `0.6` limit. It retains at least
`91.8%` of baseline inland rain and changes land rain share by at most `0.73` percentage
points. Thus it improves the error but does not pass G3.

The `alpha=3` coast arm had the lowest raw pattern shortfall (`1.0049`, worst inland/coast
ratio `1.203`) but is excluded: several A/D/E diagnostic cases retained less than `80%` of
baseline inland rain (minimum `76.9%`), and the cold CFD case shifted land rain share by
`2.01` percentage points. Its ratio-only gain does not qualify. The equal-integral uniform
arms are eligible but score worse than their coast-focused counterparts. No Round 1 arm is
freeze-eligible. Advance `alpha=1` to Round 2 and test only `L=500` and `1,500 km`.

### Reproducibility correction before Round 1

The first Round 1 attempt stopped before writing a tuning result. Its check compared the
fresh `p=0.25` arm with `ridge-strengths-with-distribution.json`, which has no embedded save,
code-commit or capture hash. That file reported finite-volume RMS values of `1.9610`
(diagnostic wind) and `2.0587` (CFD). A fresh standard-path reproduction from the registered
density-4 A capture (save SHA-256 `1f167a8d9e0edf06c3c52b23598dc2f776adc2927f843beb992b5193ad30fd67`,
capture-results SHA-256 `439336b38bd8af9e0eeabe2083b3d88f5919f4f627ded03ebe33ed9f5ac69519`,
code commit `8ea19e3af229d48a9d9c7339780951ce23fb5932`) reproduced `0.664751` and
`0.672209`, respectively. The cause of the old report's mismatch is unknown; without
provenance it is not a valid numerical baseline.

The coast tuning compares each treatment with its freshly computed `p=0.25, alpha=0` arm
and its matched `p=0` sink control. Selection uses only these within-run comparisons plus
hard conservation, cap, placement, and fixed-geography checks. The separate standard-path
run documents the old report discrepancy; it is not used as an external numeric oracle,
because its full regional metrics lack provenance. Each tuning result fingerprints the
dirty-worktree state, all diagnostic and solver source files, dependency constraints, and
Python/NumPy/SciPy versions; it also records this protocol's version identifier. Cross-round
comparisons are valid only when these fingerprints match. No parameter, mask, threshold, or
solver choice changed.
This is a pre-experiment measurement-provenance correction, not a tuning round. The initial
attempt produced no scored arm and does not consume one of the three rounds.

For fixed geography eligibility, each regional sample must cover at least `0.1%` of the
planetary surface. A ridge-shadow sample must include at least 100 valid crest/lee crossings
across at least 25 distinct crests. These thresholds were set from geography alone before
Round 1 candidate results; an absent metric may be skipped only when this support test fails.

### Round 1 — coast focus versus added land sink

Use `L=750 km`. For each native land cell compute `x=exp(-d/L)` from the archived coast
distance in km, then area-average `x` within each land reservoir to get `p_i` in `[0,1]`.
Apply `lambda_i = lambda0 * (1 + alpha * p_i)` on land and `lambda_i=lambda0` on ocean.
Do **not** exponentiate the already aggregated `p_i` again. Compare:

| Arm | Ridge interception | Land sink |
|---|---:|---|
| Control | `p=0` | `alpha=0` |
| Frozen candidate | `p=0.25` | `alpha=0` |
| Coast 1 pair | `p=0` and `0.25` | `alpha=1`, coast profile |
| Uniform matched 1 pair | `p=0` and `0.25` | uniform land multiplier equal to Coast 1's area-integrated land sink |
| Coast 3 pair | `p=0` and `0.25` | `alpha=3`, coast profile |
| Uniform matched 3 pair | `p=0` and `0.25` | uniform land multiplier equal to Coast 3's area-integrated land sink |

The matched uniform arm distinguishes coast focusing from simply increasing total land
rainout. Measure ridge gains against the matched `p=0` arm at the same sink; compare the
resulting `p=0.25` maps to the frozen `alpha=0, p=0.25` candidate for inland and land-share
guardrails. Select Coast 1 or Coast 3 only if it beats both the frozen candidate and its
matched-uniform alternative on the minimax pattern score without violating the tuning
safeguards. If neither
coast arm improves the minimax score over control, stop early; do not use the midpoint on a
mechanism already ruled out. Round 3 is reserved for a coast arm that improved the score but
still misses a formal tuning pattern threshold after the length check.

### Round 2 — one coast-distance scale

If Round 1 identifies a coast-focused winner, hold its `alpha` fixed and compare `L=500 km`
and `L=1,500 km` against its `L=750 km` result. All other inputs and controls remain fixed.
Select the eligible length with the best minimax score. If neither length improves the
`L=750 km` score, retain `L=750 km` only if it remains eligible; otherwise reject this sink
family.

Round 2 completed at both lengths for `alpha=1`. `L=1,500 km` scored `1.5257`, better than
`L=750 km` (`1.5841`) and `L=500 km` (`1.6400`); its matched uniform control scored `1.6170`.
All hard gates and tuning safeguards passed. Its worst inland/coast ratio remains `1.515`
on the cold capture with CFD winds, while retaining at least `88.0%` of baseline inland rain
and changing land rain share by at most `0.82` percentage points. Rain belts, ridge direction,
and positive cold-coast rain continue to pass. No candidate is freeze-eligible. Use
`alpha=1, L=1,500 km` for the final strength-midpoint round.

### Round 3 — one bounded strength midpoint

If the best Round 1/2 coast-focused candidate improves the minimax score but still misses an
approved tuning pattern threshold, keep its selected `L` and test only `alpha=2` against the
same-L `alpha=1` candidate and its matched uniform sink. Round 1 ruled out `alpha=3` on the
inland-retention safeguard, so do not spend this round retesting it at another length. Compare
the `alpha=2` score to the Round 2 `alpha=1, L=1,500 km` score only if their source
fingerprints match. If no coast-focused candidate improved the
Round 1 minimax score, do not try another mechanism under this plan; record G3 failed and
prepare the owner's exception. No additional parameter or fallback wind multiplier may be
introduced in this round.

If any candidate satisfies all tuning safeguards and pattern thresholds before Round 3, stop
tuning and freeze it. After freeze, verify new-held-out byte hashes and perform the registered
acceptance evaluation once. If no candidate passes the new held-outs after at most three
rounds, mark G3 failed, preserve the failed measurements, and add an explicit owner-approved
exception to unblock RFC continuation; do not imply that the gate passed.

## Three-round disposition: G3 failed by owner exception

All three authorized tuning rounds are complete. Round 1 found `alpha=1, L=750 km` eligible
with minimax score `1.5841` versus `2.0451` for its baseline; `alpha=3` was excluded by the
inland-retention safeguard. Round 2 selected the better `alpha=1, L=1,500 km` result at
`1.5257` (versus `1.6400` at 500 km and `1.5841` at 750 km). Its worst inland/coast rain
ratio was `1.515` on the cold case with CFD winds, against the formal maximum `0.6`.
Round 3 tested `alpha=2, L=1,500 km`: minimax improved to `1.1674`, and the worst
inland/coast ratio improved to `1.300`, but minimum inland-rain retention fell to `77.5%`,
below the `80%` safeguard. It is ineligible. Thus `alpha=1, L=1,500 km` is the best eligible
tuning result, but it still misses the dry-interior criterion by a factor of about 2.5 at its
worst case. None was frozen, and no new held-out acceptance run was warranted because no
candidate met the preregistered tuning eligibility and pattern requirements.

The owner explicitly authorized the exception on 2026-11-18: RFC 0019 may continue to G4
and implementation planning with the inland-rain mismatch documented as a modeling
simplification. The G3 gate remains **failed**, and the exception waives only the sequencing
block. Do not report any three-round candidate as a G3 pass or as production-selected. The
owner's stated target remains visually useful rain belts, rain shadows and dry interiors on
synthetic worlds, accepting approximate intensities and positions while the single-column
model's known errors are tracked for later work.

The round score files are private captures:

- Round 1: `~/mantle-bloom-captures/issue913/g3-round1-score.json`
- Round 2: `~/mantle-bloom-captures/issue913/g3-round2-score.json`
- Round 3: `~/mantle-bloom-captures/issue913/g3-round3-score.json`

All three scores used source fingerprint
`593332fee198422d0e50e76bd86063a7f16a60d40521794d874cf500d0091d4f`. Reproduction,
case outputs, caveats and per-round score details are retained in those files and the
respective per-case JSON files in the same private capture directory.

## Required implementation and audit checks

- Aggregate `exp(-distance/L)` from field orientation into the reservoir orientation using
  `grid.native_ids` and `grid.native_area`; test a synthetic coast/inland pair so a missing
  row flip fails.
- Require finite nonnegative distances and `p_i` in `[0,1]`; require ocean sink to equal
  `RAIN_RATE`, land sink to lie in `[RAIN_RATE, (1+alpha)*RAIN_RATE]`, and alpha zero to be
  bit-identical to the uniform control.
- Audit each arm's raw global and summed cell residuals at `1e-8`, finite/nonnegative mass
  and rain, mass cap, overflow, exact native placement totals and source equality. The solve
  may not repair a residual by adding rain.
- Refuse any non-control sink profile on a held-out capture before loading its `.npz` arrays;
  retain the existing exact tuning-save SHA checks. Every output records save digest, density,
  wind mode, solver, alpha, length, commit, dirty-worktree state, runtime versions, and hashes
  of every diagnostic and solver source file used. The scorer requires identical fingerprints
  within a round. Compare score files across rounds only when their fingerprints and registered
  save digests match.
- The scorer excludes any arm failing finite/nonnegative mass, cap, overflow, rain, global or
  summed-cell closure, native placement, or either tuning safeguard. Pattern checks remain
  separate: an otherwise eligible candidate can advance to
  another tuning round while missing them, but cannot freeze unless each required pattern has
  geographic support and passes. If geography supports a metric but its value is undefined,
  assign that pattern its full normalized shortfall; skip it only when the fixed geography test
  fails. Rank only arms passing hard gates and safeguards. The family-balanced
  tie-break prevents the three A/D/E checkpoints from outweighing the independent cold-source
  family.

No production climate code is changed in these tuning rounds. Any later wind diagnostic or
CFD work belongs to the separate wind-transport handoff.
