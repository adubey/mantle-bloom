# RFC 0018: Plates split along weak paths that are pulling apart

- **Status:** Accepted (revised 2026-11-08 after G2 sample-size review; revised 2026-11-06 after G2 non-worsening criterion review; revised 2026-11-06 after G2 metric review; revised 2026-11-04 after G1 gate review; owner confirmed 2026-11-03)
- **Date:** 2026-11-03
- **Related:** #770, #797, #802, #846; RFCs [0009](0009-craton-mosaic.md) (cratons-first worlds), [0010](0010-rifts-follow-weak-lithosphere.md) (accepted predecessor; withdrawal follows gates)

## Revision note

**2026-11-08 — G0–G2 diagnosis results recorded.** Revised G0 and the coupled-probability
G1 diagnostics pass on the 20 preregistered development worlds. G2 passes under the
owner-directed 10-pair amendment merged in [PR #845](https://github.com/adubey/mantle-bloom/pull/845):
all first 10 frozen pairs completed, zero had `D_s < 0`, and they produced 255 combined actual
splits. The G2 amendment was made after those ten outcomes were observed, so it is post-result
evidence with lower precision and supports non-worsening on the evaluated pairs, not positive
efficacy. Production implementation is tracked in [#846](https://github.com/adubey/mantle-bloom/issues/846). Gate details and compact evidence are in
[`analysis/issue797/README.md`](../analysis/issue797/README.md). Production implementation,
implementation-level ledger checks, G3 and the remaining acceptance criteria are still open.

**2026-11-08 — G2 sample-size amendment.** The owner requested reducing the G2 held-out
sample from 20 pairs to 10 because the graph runs are taking a long time. The first ten
seeds in the frozen registry are the evaluation set; all ten paired runs are complete. Keep
the existing non-worsening tolerance at 10% of the sample, rounded to one seed: at most one
of ten may have `D_s < 0`. Keep the combined 50-actual-split requirement and all other
simulation settings and outcome definitions unchanged. The ten pairs produced 255 actual
splits (157 baseline, 98 graph), zero tiny splits in either arm, and `D_s = 0` for every pair.
Thus the amended criteria pass on these ten pairs.

This is a post-result, owner-directed amendment: the first ten outcomes were observed before
the sample size was reduced. It is not preregistered or confirmatory, and the smaller sample
has less precision. The result supports non-worsening on these ten pairs; it does not show
positive efficacy. The ten-pair results and amendment are recorded in
[`analysis/issue797/g2-sample-size-amendment.json`](../analysis/issue797/g2-sample-size-amendment.json).

| Gate | Previous criterion | Revised criterion | Reason |
|---|---|---|---|
| G2 sample and non-worsening threshold | Complete all 20 frozen held-out pairs; at most two may have `D_s < 0`; at least 50 combined actual split events | Complete the first 10 seeds in the frozen registry; at most one may have `D_s < 0`; at least 50 combined actual split events | The owner requested a shorter evaluation because graph runs take a long time; preserving the 10% tolerance gives a limit of one worsened pair out of ten. The change was made after the first ten outcomes were observed. |

**2026-11-06 — G2 non-worsening criterion amendment.** This revision is proposed in
[PR #824](https://github.com/adubey/mantle-bloom/pull/824). The owner approved replacing the
positive-effect requirements after the first four held-out pairs had `D_s = 0` in every pair:
neither method produced a tiny breakup (66 baseline splits and 38 graph splits, zero tiny
events in either arm). The previous amendment required a positive mean
`D_s`, a paired bootstrap interval wholly above zero, and improvement on at least two seeds.
At that stage, the criterion passed when no more than two of the 20 held-out seeds had `D_s < 0`;
`D_s = 0` counts as non-worsening. This permits a result with no improvement if the graph
method also avoids adding tiny breakups. This is an owner-selected tolerance, not a geological
rate threshold. The seed list and simulation settings remain fixed; all four observed pairs
remain in the analysis. Because the criterion changed after four outcomes were observed,
label the result as a second amended G2 analysis, not as preregistered. Continue to report the
mean `D_s`, its 95% paired seed-bootstrap interval, and the number of seeds with `D_s > 0` as
descriptive evidence, but they are no longer pass conditions.

| Gate | Previous criterion | Revised criterion | Reason |
|---|---|---|---|
| G2 outcome threshold | Mean `D_s > 0`, two-sided 95% paired seed-bootstrap interval above zero, and `D_s > 0` on at least two seeds | No more than two of the 20 held-out seeds have `D_s < 0`; `D_s = 0` is non-worsening. Complete all 20 pairs and meet the existing 50-event requirement. | The first four pairs showed zero tiny events in both methods. Requiring positive improvement and a positive interval would reject ties even though the owner's goal is to avoid making tiny breakups worse. |

**2026-11-06 — G2 tiny-breakup endpoint amendment.** This revision merged in
[PR #823](https://github.com/adubey/mantle-bloom/pull/823). The owner approved this change after
the first two held-out seed pairs exposed a zero-baseline failure: both methods had zero tiny
breakups, so the relative reduction divided by a zero baseline share and was undefined. The
relative share also conditions on the number of splits, which can differ between methods
under the same 200 Myr exposure. The two pairs had 42 baseline splits and 19 graph splits,
with zero tiny breakups in either arm (`S_baseline = 0`, `S_graph = 0`, so the share
difference is 0 percentage points); their paired tiny-count differences are both zero.
These are preliminary outcomes from the G2 diagnostic run associated with [#797](https://github.com/adubey/mantle-bloom/issues/797). The seed
list and simulation protocol remain fixed. Because the endpoint changed after two pairs were
observed, report the revised result as an amended G2 analysis, retain the originally
preregistered share-based result as a secondary analysis, and do not describe the revised
endpoint as preregistered. Keep both completed pairs in all summaries and uncertainty
calculations.

| Gate | Earlier criterion | Revised criterion | Reason |
|---|---|---|---|
| G2 primary outcome | At least one-third relative reduction in the share of actual splits with a smaller daughter under 7% of parent area, and improvement on at least two seeds | For each seed, `D_s = T_baseline,s - T_graph,s`, the difference in tiny-breakup counts over 200 Myr. Pass when the mean paired difference is positive, its two-sided 95% paired seed-bootstrap interval excludes zero in the beneficial direction, and `D_s > 0` for at least two seeds. Use 10,000 resamples of whole seed pairs. | Relative reduction is undefined at a zero baseline share. Paired count differences are defined even when an arm has no splits or no tiny splits, and measure event burden at the fixed world exposure. |
| G2 supporting measures | Tiny-breakup share and daughter mass distribution | Report the original relative share reduction, percentage-point difference `S_baseline - S_graph`, and smaller-daughter area and mass distributions as secondary measures. | Preserve comparison with #770 while showing conditional split quality and daughter sizes. |
| G3 timestep comparison | Tiny-breakup share within 20% of the normal-step result | Keep the 20% relative comparison when the normal-step share is nonzero; if it is zero, require zero tiny breakups at half step. Report paired counts and shares. | Avoid dividing by a zero normal-step share. |

The current outcome threshold in the later 2026-11-06 amendment, further amended on
2026-11-08 for sample size, supersedes the intermediate positive-effect criterion in this
table. The intermediate endpoint and criterion were both
changed transparently after held-out outcomes had been observed.

The G2 endpoint is an amended analysis because it was selected after two seed-pair outcomes
were observed. The original relative-share result remains reported for transparency. The
seed list and the minimum of 50 actual split events across the seed set, with its frozen
reserve rule, remain unchanged.

**2026-11-04 — G1 gate review.** This revision proposed [PR #810](https://github.com/adubey/mantle-bloom/pull/810). The previous G1 required at least 10 observed opportunities
and a usable candidate in at least two thirds of trigger-positive opportunities. The 2/3 rate
was a project-selected usability cutoff, not a geological estimate. Geological rift systems
are not independent draws from the simulation's per-trigger opportunities, and observations
show both successful and failed systems whose outcomes depend on evolving force balance,
weakening and competition. The literature does not establish a universal success fraction for
this denominator; the mechanisms separating successful and failed rifts remain debated
[[Brune et al., 2023]](#ref-brune2023), and mapped failed systems illustrate the variety of
incomplete rift histories [[Elling et al., 2022]](#ref-elling2022).

The old G1 rate criterion is removed. G1 now checks structural invariants and the stated
probability model's behavior on controlled cases, while development-world outcome counts are
reported descriptively with uncertainty. No aggregate per-opportunity fraction is called a
geological pass threshold. G1 does not claim geological rate calibration; that would need a
separate mapping between simulated opportunity exposure and geological rift episodes.

## Problem

When a plate splits today, the model divides it with a great circle. The line is chosen from
the difference between two mantle-flow patterns, but it does not follow the plate's weak
places or check whether the two new plates would move apart. It can also leave a very small
piece.

Two recent efforts did not solve this. RFC 0010's earlier minimum-cost cut did not route
around cratons reliably. Issue #770 added inherited pseudo-sutures to otherwise unchanged
cratons-first worlds, but the agreed routing measure showed no improvement in tiny
breakups: the 0.2 and 0.4 route-cost trials produced more tiny outcomes overall, and no seed
improved ([#770](https://github.com/adubey/mantle-bloom/issues/770),
[0.4 trial](https://github.com/adubey/mantle-bloom/issues/770#issuecomment-6617164261)).
Lowering the cost of a weak line alone was not enough.

A new graph-cut probe made actual daughter plates in three saved worlds. One split passed
the existing size and motion checks, leaving daughters with 40.2% and 59.8% of the parent
mass. Its net opening was only +0.013 cm/yr. The other two worlds produced no split: their
candidate daughters repeatedly failed the existing minimum-size or pole-separation checks.
This is encouraging geometry evidence, but it is too small a sample and too weak an opening
to show that the method works. The saved worlds are not available for download; their names
and compact results are in the [graph-cut replay record](../analysis/exploration-graph-rift-algorithm/applied-split-replay.json).

**Design principles:** the great-circle cut is a geometry shortcut that does not account
for local weakness or outward motion, which conflicts with principle 8 (let each process
follow its own physics). This proposal uses cell area and mass (principles 4 and 5) and
requires multiple seeds, measured limits and visible results (principle 9).

The current accepted RFC 0010 proposal detects four kinds of rift setting and routes a growing
rift around cratons. This proposal replaces the great-circle split choice and the independent
split-success roll with a graph-cut method whose candidate availability, opening and
resistance contribute to one stochastic decision. It retains minimum plate age and cooldown.
The project owner accepted this RFC on 2026-11-03, authorizing the diagnosis gates. If its
gates pass, RFC 0010 will be withdrawn after its open gates, issues and dependent references
are settled. Until that follow-up is complete, RFC 0010 remains the active accepted proposal.

## Proposed solutions

**What nature tells us.** Rifts can reuse old weak zones, but their paths and growth depend on
the direction and amount of extension; some stop before the continent breaks apart
[[Vauchez et al., 1997]](#ref-vauchez1997) [[Brune et al., 2023]](#ref-brune2023). Geological
inventories document failed systems alongside rifts that reached seafloor spreading
[[Elling et al., 2022]](#ref-elling2022). The outcomes depend on changing force balance,
weakening, strengthening and competition between branches; no universal success fraction
follows from these observations, and the transition mechanisms remain debated. This supports
representing breakup as conditional and probabilistic, with weakness and extension coupled. It
does not support a fixed 2/3 success target per simulation trigger, a particular graph-cut
method or a mass exponent.

1. **Keep the great-circle cut.** It is simple and makes two clean sides, but ignores weak
   paths and the motion that would follow a split.
2. **Keep RFC 0010's detected, routed rifts.** That plan has detailed setting and route
   gates, but the #770 routing experiments found that weaker pseudo-sutures did not reduce
   tiny breakups. Continuing to tune their route cost is not the proposal tested here.
3. **Choose a graph cut using weakness and opening (recommended).** Treat each surface cell
   as a node and each shared edge as a possible break. Generate several low-cost cuts through
   sampled weak edges. Reject cuts that do not form a connected path across the plate or
   would make the sides move together. Among the rest, prefer cuts whose daughter motions
   differ more from the parent's, weighted by mass, and prefer more balanced daughter masses.
   This makes weakness, opening and split size part of the
   same choice while leaving the starting world unchanged.

## Design

### One split decision

The minimum plate age and cooldown remain hard gates. Motion separation and candidate
availability contribute to a stochastic breakup decision rather than acting as independent
deterministic pass/fail checks. The initial cratons-first world, including land, water, plate
ownership, crust and mantle thickness, is not changed. The algorithm reads the existing Hc
and Hm fields to estimate column mass; it does not change Hm. The successful split,
failed-rift thinning and craton-loss diagnostics must all follow the same graph path, so no
part of the event keeps using a different great-circle line.

This model divides the whole plate at the split event. It does not grow a crack over time or
join separate crack segments. Those questions were central to the earlier RFC 0010 design
and are outside this proposal.

For an eligible plate:

1. Cluster its cell velocities into the two motion groups already used by the split logic.
   Fit a likely rotation for each group. No shortest-path search is needed: only the split
   path is found, and that path is undirected.
2. Give each cell a structural resistance to breaking. More massive columns and cratons
   resist more; thinning, faults, volcanoes, plume heat, belts and sutures lower resistance.
   Motion fit is evaluated separately for each candidate cut, because it depends on which
   daughter group the cell would stay with. These are bounded relative scores, not measured
   rock strengths; no cell has zero or negative resistance.
3. Propose candidate cuts from sampled weak edges that connect different motion groups.
   Candidate proposals have a documented probability distribution; the probability of
   proposing no structurally valid cut is part of breakup availability. For each proposed
   cut, keep connected daughters with one dominant edge-to-edge trace as hard topology
   requirements. At least 90% of the cut length must belong to that trace; report the
   remaining fragments. Key the candidate stream to world seed, plate and event index so
   replay and timestep behavior can be checked without coupling it to unrelated random
   processes.
4. Measure signed, shared-edge-length-weighted opening `O` and candidate resistance `R` for
   each valid candidate. Use them to define a prior over that candidate's breakup
   probability. A weaker candidate must require less outward opening to reach the same
   breakup probability. This is a physical decision model, not a 95% confidence interval on
   a mean opening: do not use the old `1.96 * RMS residual` rule as its probability mapping.
   The exact prior and its numerical scale are not selected by current evidence; freeze them
   before G1 and test the declared behavior on controlled cases.
5. Use the proposal distribution and candidate probabilities to determine an event-level
   breakup probability. Draw once for breakup; if breakup is selected, draw the cut from
   the candidates conditional on their proposal and breakup probabilities. Do not apply
   the legacy independent 0.55 success chance. The potential-momentum score may rank or
   weight candidate proposals, but it is not physical momentum and must not override the
   coupled breakup probability.

For each cell `n`, estimate its mass as cell area times the mass per area of its crust and
mantle lithosphere: `m_n = A_n * (rho_c * Hc_n + rho_m * Hm_n)` (kg). The fitted rotations
are angular rates in radians/year-equivalent, as in `mantle.flow_at`. For unit position
vector `p_n` and planet radius `R`, convert the fitted daughter and parent rotations
`omega_new` and `omega_old` to surface velocities with `v_new(n) = R * (omega_new x p_n)`
and `v_old(n) = R * (omega_old x p_n)` (m/year). The potential momentum score before the
balance factor is:

```text
P = sum over cells n of m_n * |v_new(n) - v_old(n)|   [kg m/year]
```

This is a mass-weighted velocity-change ranking score, not physical momentum; it omits the
time and force history needed to represent momentum transfer. The radius factor is constant
across candidates and cancels in a within-plate ranking. Per-year units are retained, so no
seconds conversion is applied.

Let `f` be the smaller daughter's share of the parent mass. A soft balance factor is
`B = 4 f (1 - f)`, which ranges from zero for an extremely small daughter to one for an
even split. Rank candidates by:

```text
score = P * B^mu
```

`mu` only changes how strongly the score prefers a more balanced split. The first trial value
is `mu = 1`; compare `mu = 0, 0.5, 1, 2` during tuning.
No rock-fracture experiment or geology result justifies `mu > 0` or any particular value.
It is an explicit model choice, and only helps if it improves held-out simulation outcomes.

The offline prototype instead ranked cuts by reduced-mass-weighted net opening. The new score
follows the owner's stated idea of weighting each side's change from its old motion by its
mass, while outward opening and resistance feed the coupled breakup probability. The
development gate compares both scores on tuning seeds, so the additional balance factor must
earn its place before the winning version is frozen for held-out runs.

### Starting values and limits

The first tuning pass uses the prototype's 12 sampled cuts per eligible plate and its
dimensionless resistance rules. The entries below define one resistance calculation; they
are design choices, not measured strengths. For cell `n`, let `q_n` be its column areal mass
divided by the plate median, and let `r_n` be its velocity residual from the fitted rotation
of the candidate daughter that contains it. Let `r_med` be the median residual for that
candidate. Define:

```text
M_n = 0.6 + 0.4 * clamp(q_n, 0, 4)
C_n = 2.5 if cell n is craton, otherwise 1
F_n = 1 / (1 + clamp(r_n / r_med, 0, 4))
W_n = clamp(product of active weakening multipliers, 0.12, 1)
R_n = max(0.02, M_n * C_n * F_n * W_n)
edge_cost(j, k) = (edge_length / median_edge_length) * (R_j + R_k) / 2
```

An active weakening signal multiplies resistance by one minus its listed reduction. Apply
the weakening multipliers first, clamp their product to `0.12–1.0`, multiply by `M`, `C`
and `F`, then apply the `0.02` floor. For example, with median column mass, craton present,
residual equal to `r_med`, and all weakening signals active, their product is
`0.65 * 0.55 * 0.65 * 0.65 * 0.60 * 0.75 = 0.068`, clamped to `0.12`; thus
`R = max(0.02, 1 * 2.5 * 0.5 * 0.12) = 0.15`. If both cells on a median-length edge have
that resistance, its edge cost is `0.15`.

If `r_med` is zero, use `F_n = 1` for every cell in that candidate. The resistance is
candidate-specific because each cell's residual is measured against the daughter rotation
of the side it would remain with. If G1 fails, the predeclared follow-up sensitivity check
varies each factor one at a time on the development seeds. For `M`, `C` and `F`, use
`X_lambda = 1 + lambda * (X - 1)` for `lambda = 0.5` and `1.5`. For each weakening
reduction `delta`, test `lambda * delta` at those same two values of `lambda`, capped at
`0.95`. Keep the `0.12` weakening clamp and `0.02` resistance floor fixed. Record all
results; select any revised values before rerunning G1, and never use G2 seeds to tune them.

| Signal | First trial | Meaning |
|---|---:|---|
| Column mass | `0.6 + 0.4 * clamp(mass / plate-median, 0, 4)` | A thicker or denser column costs more to cut |
| Craton | `C = 2.5` when present, otherwise `1` | A craton is harder to cut, not an absolute no-go zone |
| Motion fit | `F = 1 / (1 + clamp(residual / median-residual, 0, 4))` | Residual is to the fitted rotation of the candidate side containing the cell; a poor fit lowers resistance to breaking from that side |
| Belt / suture | Multiply by `0.65` / `0.55` when present | These are reductions of `0.35` / `0.45` |
| Volcano / plume heat | Multiply by `0.65` each when present | Each is a reduction of `0.35` |
| Fault / thinning | Multiply by `0.60` / `0.75` when present | These are reductions of `0.40` / `0.25` |
| Weakening limit | Clamp product of active multipliers to `0.12–1.0` | Applied before the other resistance factors |
| Resistance floor | `max(0.02, ...)` | Applied after all other factors |
| Candidate count | `12` | Number of sampled weak edges; tune only on development seeds |
| Balance exponent | `mu = 1`; sweep `0, 0.5, 1, 2` | Design choice; not a measured physical value |

The motion-residual reduction is candidate- and assignment-specific: the fit is measured
against the daughter group a cell would remain with, never against whichever of the two
groups happens to fit it better. All score factors stay bounded. Area is used when converting
column mass to a plate or daughter total; raw node count is not a substitute for mass.

A failed breakup attempt thins the existing failed-rift band along its proposed graph path.
The removed 0.55 production chance is not multiplied into the new event probability.

### What earlier work contributes

| Earlier work | Carries forward | Does not establish |
|---|---|---|
| RFC 0010's #480 cut trials | Avoid relying on one static cheapest path; test connected daughters and craton cost | That a different graph score works |
| RFC 0010 Gates D1–D3 and #699 | The 200 Myr history, plate-count range, craton-loss budget, timestep and runtime checks | That the old detected-setting and route gates are needed for this split-only change |
| #770 geometry and maps | Preserve the generated world and inspect land, water, cratons, belts and sutures | That cheaper pseudo-sutures improve actual breakups; the 0.2 and 0.4 route trials failed |
| Three graph-cut replays | The method can make a plausible 40/60 split | A reliable rate or strong opening; only one split occurred and its opening was tiny |

Results from every row above are development evidence. None of its seeds may be counted in
G2's held-out sample.

### Diagnosis gates

Use the three saved worlds above only to develop the method. They have already shaped the
proposal and cannot count as held-out evidence. The #770 seeds, the 20-world pseudo-suture
census and all other seeds used to choose graph weights or `mu` are also tuning evidence.

1. **G0 — Make the geometry and score measurable.** On 20 predeclared development seeds,
   record every eligible split attempt, sampled edge, candidate rejection, daughter mass
   share, cut length, trace continuity, signed opening `O`, candidate resistance `R`,
   proposal outcomes, score and runtime. Use area-derived column mass. Keep every failed
   world-generation seed in the report; add only predeclared replacements to reach 20
   constructible baseline worlds, following the #770 census practice. Record fitted-motion
   residual structure as diagnostic information, not as a confidence interval on mean
   opening. Do not impose an opening floor derived from `1.96 * RMS residual`; instead,
   define and freeze the `O,R`-to-breakup-probability mapping before G1. The pilot's tiny
   `+0.013 cm/yr` opening remains a descriptive calibration case, not a binary pass.
   Compare the potential-momentum score with the prototype's reduced-mass-times-opening
   score and compare `mu = 0, 0.5, 1, 2`; resistance values remain fixed at the starting
   values in the table. Freeze the score, `mu`, candidate proposal rule and the `O,R`-to-prior
   mapping before G1. If G1 reveals a failed monotonic response or structurally invalid
   applied split, diagnose the mechanism before changing weights; use the predeclared
   one-at-a-time resistance sweep only on development seeds, then rerun G0/G1. Report
   uncertainty by resampling worlds, not individual opportunities. Do not tune against G2
   seeds.
2. **G1 — Check the coupled breakup model.** Before the gate, freeze the event-level
   `P_break` mapping and candidate proposal rule. A single uniform draw decides breakup from
   that probability; do not apply the legacy independent 55% success roll. The model must
   include candidate availability and couple outward motion to resistance so that, with
   other inputs held fixed, more opening cannot lower breakup probability and weaker
   resistance cannot raise the opening required for the same probability. Validate the
   probability-to-draw implementation on controlled cases with known probabilities. Every
   applied split must produce two connected daughters, a dominant edge-to-edge trace, and
   pass the existing hard daughter-size and ownership checks. Treat pole separation and the
   motion difference as probabilistic inputs to `P_break`, not independent hard pass/fail
   filters. On the 20 development seeds, report every eligible opportunity, candidate
   proposal, rejection reason, `P_break`, random draw, outcome, and per-world outcome rates
   with uncertainty. There is no fixed usable-opportunity percentage for passing G1: the
   geological record does not define one for this denominator. This gate tests internal
   probability behavior and structural safety, not Earth's breakup frequency.
3. **G2 — Compare actual breakups on held-out worlds.** Freeze a 20-seed registry before
   running; the amended gate evaluates its first 10 seeds. Compare the production
   great-circle method and graph-cut method from identical generated worlds, settings, step
   schedule and random streams for 200 Myr. Use the same
   A tiny breakup is an actual split whose smaller daughter has less than 7% of the parent
   plate's area, preserving the event definition used in #770. For each seed `s`, let
   `T_baseline,s` and `T_graph,s` be the tiny-breakup counts over the fixed 200 Myr horizon,
   and define `D_s = T_baseline,s - T_graph,s`. Under the 2026-11-08 owner-directed amendment,
   the graph method passes only if the first 10
   pairs in the frozen registry are complete, no more than one has `D_s < 0`, and the existing
   combined 50-actual-split event threshold is met. This sample-size reduction was made after
   the first ten outcomes were observed; describe the result as post-result amended evidence,
   not as a preregistered confirmatory result. A tie (`D_s = 0`) counts as
   non-worsening. Report the mean `D_s`, its two-sided 95% paired seed-bootstrap interval, and
   the number of seeds with `D_s > 0` descriptively; resample whole seed pairs 10,000 times
   and never resample individual split events as independent observations. For method `m`,
   calculate its overall tiny-breakup share as
   `S_m = sum_s(T_m,s) / sum_s(B_m,s)`, where `B_m,s` is its actual split count. Report the
   original relative share reduction only when `S_baseline > 0` and both shares are defined;
   report the percentage-point difference `S_baseline - S_graph` when both shares are
   defined. If a method makes no actual splits, report its share as undefined; the primary
   count difference remains defined. A zero baseline tiny-breakup count provides no positive
   reduction to estimate; evaluate non-worsening with the held-out seed sign criterion above.
   Also report per-seed counts and
   smaller-daughter area and mass distributions. The comparison must include at least 50
   actual split events across the seed set; if it does not, continue only with seeds from a
   reserve list frozen before the run. Do not infer success from a favorable image or
   routing-only proxy, and do not
   impose a hard 30% daughter-mass floor.
4. **G3 — Check cratons, stability and cost.** The existing production split-success gate
   uses `round(elapsed_years)` in its random key. That is a known step-dependent carryover;
   keep it for the first comparison and judge its effect with this timestep check rather
   than changing it here. Preserve the RFC 0010 limits unless the
   owner changes them: rifting-related craton loss at or below 2% of live craton per Gyr,
   and 10–40 plates at every output. Repeat representative held-out runs at half the step
   size on three predeclared held-out seeds; plate counts and craton loss must stay within 20%
   of normal-step results. Compare tiny-breakup shares within 20% when the normal-step share
   is nonzero; when it is zero, the half-step run must also have zero tiny breakups. Report
   the paired tiny-breakup counts and shares. Median full-simulation runtime must remain within 2% of the paired current
   production baseline. Record any seed where the model cannot finish or maintain its
   ledgers.

### Acceptance criteria

- [x] G0 passes on the 20 development seeds, with all seed failures and rejection reasons
      kept in the report; see `analysis/issue797/README.md` and `g0-prior-summary.json`.
- [x] G1 validates candidate availability and the coupled breakup probability on controlled
      cases, including the required monotonic responses; development-world event rates and
      uncertainty are reported without a fixed per-opportunity success fraction in
      `analysis/issue797/README.md` and `g1-results.json`.
- [x] G2 passes on the first 10 seeds in the frozen held-out registry: no more than one seed has
      `D_s < 0`, ties count as non-worsening, and the event-count requirement is met. The
      sample-size amendment was made after these outcomes were observed. Report
      the mean paired difference, its 95% seed-bootstrap interval, improved-seed count,
      `S_baseline - S_graph` in percentage points and the original relative share reduction
      when defined. Report all three G2 amendments transparently; each was adopted after
      outcomes were observed. Production implementation and remaining checks are tracked in
      [#846](https://github.com/adubey/mantle-bloom/issues/846).
- [ ] Every applied split makes exactly two connected plates, assigns each parent cell to
      one daughter, and preserves crust and mantle ledgers. Report daughter mass shares and
      cut trace fragments; no fixed 30% mass minimum is added.
- [ ] G3 passes: craton loss at most 2% per Gyr, plate count 10–40 at every output,
      timestep metrics within 20%, and median runtime within 2% of baseline.
- [ ] Old saves load and can continue stepping. The method adds no saved fields unless a
      later revision explains why they are needed and how old saves initialize them.
- [ ] Unit tests cover connected daughters, boundary-to-boundary traces, probability bounds
      and monotonic responses to opening and weakness, balanced-score ordering, area-weighted
      mass, failed-rift thinning along the graph path, exact ownership and old-save loading.
- [ ] Show five predeclared held-out world snapshots at the first split, or show the
      no-split result for any seed with no split in 200 Myr. Use a shared projection and
      legend. Show land, water, cratons, belts, sutures, the graph-cut trace and both daughter
      plate outlines so a reader can inspect the full crack and resulting plates.
- [ ] Document the new split choice in `docs/simulation-model/tectonics.md` and update the
      architecture module map if the implementation adds or moves a module.

### Falsification criteria

- **Weakness does not guide real splits.** If the chosen cuts do not follow the weak marks
  more often than matched control cuts, the edge-cost hypothesis is unsupported. Revisit the
  cost signals; do not claim pseudo-sutures helped solely because they appear under a cut.
- **Mass weighting worsens tiny-breakup burden.** If more than one of the first ten evaluated
  held-out seeds has `D_s < 0` after the weights and `mu` have been frozen on development
  seeds, the method
  exceeds the owner's non-worsening tolerance; reject or revise the score. Do not tune on the
  held-out seeds.
- **The split score rewards closing or broken cuts.** If high-scoring cuts often close, break
  into disconnected traces, fail daughter checks or create extra fragments, the score is
  incomplete. Revise it before a production proposal.
- **The breakup probability ignores the coupled controls.** If controlled cases show that
  increasing outward opening lowers `P_break`, or weakening resistance raises the required
  opening at a fixed probability, the proposed probability model is internally inconsistent;
  revise it before G2.
- **The result depends strongly on step size.** If G3's 20% limit fails, make the candidate
  sampling and event timing stable across timestep changes or reject the design.
- **Craton safety or bookkeeping fails.** If the 2% craton-loss bound, connected ownership,
  or material ledgers fail, stop implementation and diagnose the path before changing the
  bound.

## Detailed design (provisional touchpoints)

- `backend/app/merge_split.py`: replace great-circle cut selection inside
  `maybe_split_plate`; retain minimum age and cooldown, replace deterministic motion and the
  independent success roll with one coupled event probability and keyed stochastic draw.
  Pass the selected daughter mask and crossed-edge trace through split, failed-rift and
  diagnostic handling.
- `backend/app/plates.py` and `backend/app/sparse_quad_patch.py`: add a mask-based split
  helper that validates daughter sizes, calls `_plates_from_node_masks()` and resets both
  daughter ages. `maybe_split_plate` uses this helper instead of `Plate.split(cut_normal, ...)`;
  leave the existing great-circle `Plate.split` signature in place for any other callers.
- Add one trace-distance helper that measures each cell's angular distance to the polyline
  through crossed-edge midpoints and returns its distance and in-band mask. Keep the current
  band width (`FAILED_RIFT_BAND_MULT * line_spacing_rad`). Change
  `Plate.apply_failed_rift` to receive the distance array and band width so its thinning
  taper follows the trace. Change `cratons.record_rift_cut` to receive the in-band mask.
  The success-path belt call already accepts a mask, so pass it the same trace-derived mask;
  no `cratons.mark_belt` signature change is needed. These operations share one trace-band
  definition, avoiding any remaining great-circle band in this split path.
- `backend/app/mantle.py`: reuse the existing flow sampling and Euler-fit functions; do not
  add shortest-path work.
- `backend/unit_tests/test_merge_split.py`, `backend/unit_tests/test_sparse_quad_patch.py`,
  `backend/unit_tests/test_cratons.py` and `backend/unit_tests/test_persistence.py`:
  cover candidate validity, conservation, timestep-stable selection, and old saves.
- `docs/simulation-model/tectonics.md`: explain when the graph cut is evaluated, how it is
  chosen, and what the relative cost and balance score do not claim.

No persistent state is expected. Each eligible plate needs one sparse partition solve and
12 candidate checks; the full work grows with the plate's cells and shared edges. Measure
the actual end-to-end cost at production resolution; do not assume sparse graph work is free.

## Assumptions and uncertainty

- The model's per-cell Hc and Hm are adequate proxies for how much material a cut separates.
  They are read-only inputs to the score; neither field is changed by this RFC.
- Existing belt, pseudo-suture, collision-suture, fault, volcano, plume and thinning marks
  can guide a possible cut. Their weights are not laboratory-calibrated strengths.
- Generated pseudo-sutures are not in the `origin/main` base used for this draft; their
  marker is still part of #770 work. Use that distinct signal only after it lands or is
  added through an accepted implementation plan. This RFC does not change world generation.
- Fitted daughter rotations approximate the motion after separation well enough to estimate
  outward opening as an input to the breakup prior. This needs held-out replay evidence.
- A mass-weighted velocity change can rank split options. It is a heuristic score, not
  angular momentum, a force law or a conservation ledger.
- A balance exponent may reduce very small daughters, but no evidence establishes its value.
  Its value is chosen on development seeds and judged only on separate held-out seeds.
- The #770 routing results establish that lowering suture route cost alone did not reduce
  tiny breakups in that router. They do not test this graph-cut selection method.
- Earth's catalog of successful and failed rifts does not yet provide a calibrated rate per
  simulation trigger. Geological systems, model opportunities and their exposure times are
  not equivalent units; development-world split fractions are descriptive until that mapping
  is made.
- The three applied graph-cut splits are an exploratory demonstration, not a success-rate
  estimate. One had only +0.013 cm/yr net opening and two made no valid split.

## References

- <a id="ref-vauchez1997"></a>Vauchez, A., Barruol, G., & Tommasi, A. (1997). Why do
  continents break-up parallel to ancient orogenic belts? *Terra Nova*, 9, 62–66.
  [doi:10.1111/j.1365-3121.1997.tb00003.x](https://doi.org/10.1111/j.1365-3121.1997.tb00003.x)
- <a id="ref-brune2023"></a>Brune, S., Kolawole, F., Olive, J.-A., et al. (2023).
  Geodynamics of continental rift initiation and evolution. *Nature Reviews Earth & Environment*,
  4, 235–253. [doi:10.1038/s43017-023-00391-3](https://doi.org/10.1038/s43017-023-00391-3)
- <a id="ref-elling2022"></a>Elling, R., Stein, S., Stein, C. A., & Gefeke, K. (2022).
  Three major failed rifts in central North America: Similarities and differences. *GSA Today*,
  32(6), 4–11. [doi:10.1130/GSATG518A.1](https://doi.org/10.1130/GSATG518A.1)
