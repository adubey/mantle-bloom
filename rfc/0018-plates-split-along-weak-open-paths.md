# RFC 0018: Plates split along weak paths that are pulling apart

- **Status:** Accepted (revision proposed 2026-11-14 after #848 runtime evidence and #876 rift-candidate exploration; revision proposed 2026-11-09 after spectral-weighting, runtime and failed-rift review; revised 2026-11-08 after G2 sample-size review; revised 2026-11-06 after G2 non-worsening criterion review; revised 2026-11-06 after G2 metric review; revised 2026-11-04 after G1 gate review; owner confirmed 2026-11-03)
- **Date:** 2026-11-03
- **Related:** #487, #699, #770, #797, #802, #836, #846, #848, #850, #851, #866, #876; RFC [0007](0007-mantle-plumes-and-plateaus.md) (plumes); RFCs [0009](0009-craton-mosaic.md) (cratons-first worlds), [0010](0010-rifts-follow-weak-lithosphere.md) (accepted predecessor; withdrawal follows gates)

## Revision note

**2026-11-14 — damage bands replace the eigensolve candidate search (revision proposed after #848 and #876).**
The previous version chose a split by a weighted-Laplacian graph cut: up to twelve candidate cuts per
plate, each refined by repeated sparse eigensolves, filtered by structural checks, scored by an
`O/R/S` breakup prior and drawn as one keyed stochastic event. Failed attempts thinned a band along
the drawn path.

The evidence against keeping it:

- **Runtime ([#848](https://github.com/adubey/mantle-bloom/issues/848), [`analysis/issue848/`](../analysis/issue848/README.md)).**
  Eigensolves are about 85% of search time. A perfect exact skip of opportunities with no valid
  candidate still leaves about 1.84 times the baseline step time. Lower refinement caps change
  partitions and `P_break`, and the one exact shortcut saves 16% of solves.
- **Validity is geometric ([#876](https://github.com/adubey/mantle-bloom/issues/876), [`analysis/issue876/`](../analysis/issue876/README.md)).**
  Structural validity is unrelated to extensional load relative to strength (AUC 0.47–0.56), so
  valid cuts are bookkeeping candidates, not physical ones. Valid cuts carry about 7% of column
  capacity under one step's load.
- **Damage needs softening ([#699](https://github.com/adubey/mantle-bloom/issues/699), [`analysis/issue487/`](../analysis/issue487/README.md)).**
  Damage that does not weaken rock fragments plates. The absolute strength threshold depends on an
  unconstrained background stress, so it must be calibrated.

The owner then set the goal and rules for this revision in [#876](https://github.com/adubey/mantle-bloom/issues/876#issuecomment-6874579897):
efficient, realistic-looking rifting rather than new science; parameters with literature defaults
and a tuning search; loads from mantle drag, slab pull and plumes; bands that start about
100–150 km wide and grow; physical candidates that can start anywhere; and gates on realism and
runtime. This revision replaces the candidate search, prior and draws with rift strain stored in
each cell. Strain grows where load exceeds a softened strength and thins the band in place, and a
split follows the cheapest path through the remaining strength once the intact ligament fails. It
keeps the minimum age, cooldown, mask-based split helper, trace-band craton diagnostics, tiny-breakup
definition, held-out sample and the G3 limits.

| Gate | Previous criterion | Proposed criterion | Reason |
|---|---|---|---|
| G0 | Candidate geometry, opening, availability and solve counts on 20 development worlds for direct-cost weighting | **D0:** calibrate the cheap directional load tensor against the thin-shell diagnostic: manufactured sign cases, then sign agreement, principal direction, rank correlation and magnitude ratio on development worlds; otherwise switch to periodic shell solves | The candidate search is removed; the load proxy is the new uncertain input |
| G1 | `O/R/S` prior monotonicity and path-distribution checks | **D1:** predeclared parameter search on development worlds against the realism targets, including the #699 fragmentation check | The stochastic prior is removed; tuning replaces probability calibration, as the owner directed |
| G2 | Ten revised held-out pairs (registry 12–20 plus reserve `1479220500`); at most one `D_s < 0`; at least 50 combined actual splits, with the frozen reserve-list fallback | **D2:** same ten seeds, same tiny-breakup non-worsening limit, same combined 50-split count and reserve fallback, plus every realism target | Owner gates on realism; the held-out sample is still unrun for any revised method |
| G3 | Runtime within 2%; craton loss at most 2% per Gyr; 10–40 plates; half-step plate count, craton loss and tiny-breakup share within 20% (zero-share rule) | **D3:** all accepted limits unchanged, with breakup count and rift duration added to the half-step comparison | Owner gates on runtime; limits unchanged |

| Item | Verdict | Reason |
|---|---|---|
| G0–G2 results for either weighting | Historical | They evaluate the removed candidate method |
| G3 | Re-scoped as D3 | Limits unchanged; run after implementation |
| [#846](https://github.com/adubey/mantle-bloom/issues/846) production implementation | Blocked, re-scope | Implement the damage-band method after D0 and D1 |
| [#797](https://github.com/adubey/mantle-bloom/issues/797) tracking | Re-scope | Replace the G0–G3 checklist with D0–D3 |
| [#848](https://github.com/adubey/mantle-bloom/issues/848) | Close after a return comment | Every measured implementation and tested exact shortcut stays far above G3's limit; untested equivalent optimizations are not ruled out, but the method they would speed up is replaced. Before closing, post a consolidated return on #848: the runtime findings, the failed-rift finding (the experimental callback never applied failed-rift thinning; failed rifts are now bands left in the rock) and the G0–G3 to D0–D3 gate change |
| [#836](https://github.com/adubey/mantle-bloom/issues/836), [#866](https://github.com/adubey/mantle-bloom/issues/866) | Addressed in part | Interior localization is the proposed fix for wide subsidence; margin `_open_rift` is unchanged |
| [#699](https://github.com/adubey/mantle-bloom/issues/699) | Evidence | Its fragmentation result motivates softening; no further gate depends on it |
| RFC 0010 | Unchanged | Remains accepted until this RFC's gates and lifecycle decision settle |

**2026-11-09 — weak-path weighting, failed-rift outcomes and runtime evidence (revision proposed in #848).**
The accepted method used `conductance = 1 / edge_cost` in its weighted graph Laplacian even though
`edge_cost` is the resistance to cutting. Since the Laplacian quadratic form is
`xᵀLx = Σ wᵢⱼ(xᵢ - xⱼ)²`, this mapping penalizes separating the low-resistance edges that the
proposal intends to cut. In the controlled comparison in
[#850](https://github.com/adubey/mantle-bloom/issues/850), the inverse mapping crossed both
planted weak links in 0/100 cycles; direct cost crossed both in 100/100. Across 68 frozen G0
opportunities, direct cost shared no candidate partition fingerprints with the inverse mapping,
changed the selected partition in all 68, and selected the same anchor in 9/68. Direct cost
therefore matches the weak-link direction better, but it is not a validated final method: its
selected development cuts had median net opening `−0.059 cm/yr`, so geometry, availability and
the accepted `O/R/S` prior must be rechecked before held-out runs. The direct-cost snapshot
comparison made 6,072 sparse eigensolves in 259.6 s versus 6,481 in 311.6 s for inverse cost;
these are diagnostic timings, not full-simulation runtime evidence. [PR #851](https://github.com/adubey/mantle-bloom/pull/851)
records the reproducible comparison.

The first ten G2 pairs had median wall times of 8,363.5 s for the graph callback and 505.2 s
for baseline (16.25×), but the run record lacks sufficient isolation metadata to call this a
controlled G3 failure. The graph callback also returned on a failed breakup without applying
the RFC-required thinning and craton diagnostics along a selected path. This means G2 remains
evidence for the experimental callback and the old candidate method only; it does not establish
runtime or full failed-rift behavior for the revised method. The accepted G3 runtime limit
remains unchanged until a controlled paired measurement exists. No output-preserving runtime
optimization has been demonstrated, so the candidate method's ability to meet the 2% limit
remains unknown.

The inverse weighting conflicts with principle 8 because it steers the partition away from
the weak paths the breakup process is intended to use. The solve cost raises principle 2's
runtime constraint, and principle 9 requires controlled evidence before calling G3 failed or
passed.

This revision proposes `edge_cost` as the positive, dimensionless Laplacian adjacency weight,
preserving the cost meaning: lower-cost edges carry a smaller separation penalty. It also makes
the event outcome explicit. For candidate `i`, let `qᵢ` be its proposal probability and `pᵢ` its
frozen breakup prior; let `q₀` be the probability of proposing no structurally valid path. Draw
breakup with `P_break = Σᵢ qᵢpᵢ`. On breakup, retain the conditional path distribution
`qᵢpᵢ / P_break`. On no breakup, choose a failed path with conditional mass
`qᵢ(1−pᵢ) / (1−P_break)`, or choose no path with mass `q₀ / (1−P_break)`. The same path then
drives the split or failed-rift thinning and craton diagnostics. The no-path outcome makes no
geological mark. This is the candidate-first proposal-and-outcome model written in event-first
form; it preserves the successful path distribution and gives a failed attempt a path without
applying a second independent success roll. Event and path draws use separate keyed streams.

| Item | Verdict | Reason |
|---|---|---|
| G0 | Must rerun | The Laplacian weighting changes candidate partitions and selected traces; remeasure geometry, candidate availability, opening, score and runtime on development worlds. |
| G1 | Must rerun | Revalidate the frozen `O/R/S` prior and structural checks against the revised candidate distribution, and test failed-path selection and no-path outcomes. |
| G2 | Must rerun | The existing pass is for the inverse-weight experimental callback and omits failed-rift effects. Use only still-unrun held-out seeds for the revised method after the revised sample is frozen. |
| G3 | Blocked | Run after the method and failed-rift semantics are implemented; use controlled paired full simulations. Keep the 2% runtime limit pending valid measurements. |
| #846 | Blocked | Keep production implementation paused until the revision and required G0/G1 gates are accepted and pass. |
| #847 | Unchanged as historical evidence | It records the old method's G0–G2 results; it is not validation for this revision. |
| #850 / #851 | Unchanged as completed diagnostic evidence | The controlled weighting result stands; it motivates a revised candidate but does not validate it. |
| #797 | Re-scope after revision decision | Update the gate checklist and implementation dependency after the owner decides on this revision. |
| RFC 0010 | Blocked | It remains accepted until the revised RFC 0018 gates and lifecycle decision are settled. |

The first ten G2 seeds remain excluded from tuning and revised-method validation. Of registry
seeds 11–20, nine remain wholly unrun. Seed `205240614` has a completed baseline-only run and
should not be described as unseen; this revision excludes it from a clean revised-method
sample. To retain ten pairs, the proposed new sample is the nine untouched registry seeds
12–20 plus the first frozen reserve seed, `1479220500`, subject to the owner's approval of this
revised gate. The remaining 29 reserve seeds stay available under the preregistered reserve
rule. No seed in the proposed sample is used for development or tuning.

| Gate | Previously accepted evidence or criterion | Proposed revision | Reason |
|---|---|---|---|
| G0 | Passed on 20 development worlds for inverse weighting | Repeat on 20 development worlds with direct-cost weighting; compare with uniform control and report candidate geometry, opening, availability and solve cost | Direct cost changed the selected partition in all 68 measured opportunities, so the old geometry and runtime evidence does not apply. |
| G1 | Passed controlled probability and structural checks for the old candidate flow | Repeat the `O/R/S` checks and validate success-path, failed-path and no-path distributions | The candidate paths change, and the old callback had no failed-path outcome. |
| G2 | The amended gate evaluated the first ten registry seeds; at most one of ten could have `D_s < 0` and the pair set needed at least 50 actual splits | For the revised method, evaluate seeds 12–20 plus frozen reserve seed `1479220500`; retain the one-of-ten and 50-split limits | The first ten have already evaluated the old method. Seed 11 (`205240614`) has a baseline-only run, leaving nine wholly unseen registry seeds. This sample change needs owner approval before use. |
| G3 | Median runtime within 2% of paired production baseline; existing craton, plate-count and timestep limits | Keep all limits. Use paired, serialized, same-host and same-settings full runs with profiling separate | Existing wall times are not controlled enough to justify changing the 2% limit. |

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
piece. A plate's interior never stretches before it splits. Extension happens only at plate
margins, where `quad_tectonics._open_rift` thins the cells within two rows of newly inserted
cells ([#876 step 3](https://github.com/adubey/mantle-bloom/issues/876#issuecomment-6874155833)).
Rifts therefore do not form, grow, fail or leave scars inside plates. They appear as
instantaneous cuts.

Two earlier fixes did not solve this. RFC 0010's minimum-cost cut did not route around cratons
reliably. Issue #770's inherited pseudo-sutures did not reduce tiny breakups: the 0.2 and 0.4
route-cost trials produced more tiny outcomes and no seed improved
([#770](https://github.com/adubey/mantle-bloom/issues/770)).

The graph-cut method that this RFC first proposed carries a serious runtime risk and does not
measure physics:

- **Runtime.** Eigensolves are about 85% of its search time, with up to 97 sparse solves per
  decision. 94.6% of mature-world opportunities have no structurally valid candidate, yet a
  perfect exact skip of those still leaves about 1.84 times the baseline step time against G3's
  2% limit. Lower refinement caps change partitions and `P_break`. The only exact shortcut found,
  stopping at the first repeated partition, saves 16% of solves. Untested output-preserving
  optimizations (warm-started solves, candidate reuse across steps) are not ruled out, but none
  tested comes near the limit
  ([`analysis/issue848/README.md`](../analysis/issue848/README.md)).
- **Validity is geometry, not load.** On 125 opportunities (1,491 candidates), neither a
  drag-based load map nor a thin-shell stress map separates valid from invalid cuts (rank AUC
  0.47–0.56). Plate size does (0.86). Breakup probability among valid cuts does follow load
  (AUC 0.89 on the thin-shell map) ([`analysis/issue876/README.md`](../analysis/issue876/README.md)).
- **One step's load is far below rock strength.** Along valid cuts the median thin-shell load is
  about 7% of column capacity. RFC 0010 Gate 2 found about 12% on its cuts
  ([`analysis/issue487/README.md`](../analysis/issue487/README.md)). Intact rock does not fail
  under a single step's load.

The offline damage experiments for RFC 0010 Gate 3b
([#699](https://github.com/adubey/mantle-bloom/issues/699)) let damage accumulate under slab pull
and mantle currents. They produced broad damage that fragmented plates into hundreds or
thousands of pieces, never two daughters. The diagnosis recorded there was that damaged rock
did not weaken, so a crack could not concentrate load; failures that were not kept also healed.
Issue [#866](https://github.com/adubey/mantle-bloom/issues/866)'s wide below-sea-level basins,
which #836 traces to margin thinning, show the same missing feedback from the other side:
extension spread too wide.

**Design principles:** the great-circle cut and the instantaneous graph cut are geometry
shortcuts that do not let rifting follow its own physics (principle 8). The graph cut's runtime
conflicts with principle 2. This proposal keeps extension, weakening and thinning in physical
units (principles 4 and 5) and judges the result on multiple seeds and visible maps
(principle 9).

The project owner's goal for this revision ([#876 decisions](https://github.com/adubey/mantle-bloom/issues/876#issuecomment-6874579897))
is a rifting process that runs efficiently and makes realistic-looking planets. It is not a
proof of new rift mechanics. Rifting on Earth is well documented, so the model should reproduce
its look and timing cheaply. Weakly constrained quantities are parameters, with literature
defaults where they exist and a tuning search against Earth targets.

## Proposed solutions

**What nature tells us.**

- Rifts reuse old weak zones. Of eighteen Atlantic and Indian Ocean conjugate margin pairs, all
  but the Madagascar–Seychelles–India system opened along former sutures
  [[Buiter & Torsvik, 2014]](#ref-buiter2014). Continents tend to break up parallel to ancient
  orogenic belts [[Vauchez et al., 1997]](#ref-vauchez1997).
- Rifts have a long slow phase, under 10 mm/yr full rate, then an abrupt increase in divergence
  before continental rupture. Brune et al. attribute the increase to the nonlinear decay of rift
  strength [[Brune et al., 2016]](#ref-brune2016). In model terms, weakening feeds back on
  extension.
- Rifts localize. Extension may concentrate in a narrow zone or spread over a region wider than
  the lithosphere is thick, depending on crustal thickness and heat flow [[Buck, 1991]](#ref-buck1991).
  Numerical rift models produce strong localization by letting frictional strength decrease
  with accumulated strain [[Huismans & Beaumont, 2003]](#ref-huismans2003).
- Many rifts fail. In failed rifts the force balance becomes dominated by lithospheric
  strengthening, and activity stops before rupture [[Brune et al., 2023]](#ref-brune2023).
  Geological inventories document failed systems beside rifts that reached seafloor spreading
  [[Elling et al., 2022]](#ref-elling2022).
- Thinned lithosphere cools back toward its original state on a thermal time constant of about
  62.8 Myr for 125 km lithosphere [[McKenzie, 1978]](#ref-mckenzie1978).
- Plumes help. Twelve of the eighteen margins are linked to large igneous provinces, with
  breakup at the same time as the province [[Buiter & Torsvik, 2014]](#ref-buiter2014).
  Courtillot et al. argue that a plume is a prerequisite for opening a major ocean basin
  [[Courtillot et al., 1999]](#ref-courtillot1999). Some margins broke up without one.

**Alternatives.**

1. **Keep the great-circle cut.** Simple and clean, but it ignores weakness, opening and rift
   history, and makes interiors rigid until the cut.
2. **Keep RFC 0010's detected, routed rifts.** Detailed setting and route gates, but #770 found
   that weaker pseudo-sutures did not reduce tiny breakups, and it shares the instantaneous-cut
   limitation.
3. **The graph-cut candidate search this RFC first proposed.** It chooses a well-formed cut, but
   its validity checks measure geometry, its eigensolves cost several times the step budget, and
   no exact shortcut brings it within G3 (#848).
4. **Damage bands that grow, localize, fail or break through (recommended).** Every cell carries
   accumulated rift strain. Loads from mantle drag, slab pull and plumes strain cells whose load
   exceeds their strength. Strain weakens the rock, so extension concentrates into narrow bands
   that thin, subside and widen over time. A band whose remaining intact ligament fails under the
   concentrated load splits the plate along the cheapest path through the remaining strength.
   A band whose load stops leaves a scar that heals partly and can reopen. All operations are
   local or a single path search, so no eigensolve is needed.

Alternative 4 matches the documented behavior: slow then fast, along inherited weakness, often
failing, aided by plumes. It turns the #699 lesson into the central design choice: damage must
soften rock so that bands localize. It also addresses principle 2's runtime constraint.

## Design

### Breakup bookkeeping kept from the accepted method

- Minimum plate age (`SPLIT_MIN_AGE_YEARS`) and split cooldown remain hard gates.
- Every applied split goes through a mask-based split helper that calls
  `_plates_from_node_masks()`, validates daughter sizes and ownership, and resets daughter ages.
- The success path's craton-loss diagnostic (`cratons.record_rift_cut`) and belt marking use one
  trace-derived band around the split path. No great-circle band remains in the split path.
- The generated starting world is not changed.

The candidate search, its `O/R/S` breakup prior, proposal distribution, momentum score, keyed
breakup draws and great-circle failed-rift thinning are removed. Failure is now a band that stops
growing, not a sampled event.

### State in the rock

Two new per-cell surface fields move with the rock:

- `rift_strain` (dimensionless, ≥ 0): accumulated extensional strain. It never decays and
  records the permanent, fault-like part of the damage.
- `rift_warm_strain` (dimensionless, 0 ≤ `rift_warm_strain` ≤ `rift_strain`): the share of that
  strain whose thermal weakening has not yet healed. It decays with time constant `τ_heal`.

Both are intensive (area-weighted on remap). Old saves load with both at their zero defaults
(`surface_fields.SURFACE_FIELDS`), so a pre-revision world starts with no rift damage. Its belts,
faults and plume heat still mark inherited weakness.

### Each step, for each eligible plate

**1. Load, `N` (N/m, extensional, depth-integrated).** Each cell carries a symmetric 2×2
membrane-resultant tensor `𝐍` in its tangent plane, summed from cheap directional terms that
use forces the model already computes. Tension is positive. The extensional load is the largest
principal value, `N = max(0, λ_max(𝐍))`, and its principal direction gives the band's strike
(perpendicular to the largest tension). Opposing contributions cancel, so compression cannot
register as extension.

- **Basal drag:** traction `τ = c_drag (v_mantle - v_plate)` (Pa, from
  `torque.basal_drag_force_vectors_n` divided by cell area). Along each shared edge with unit
  tangent `t̂`, the signed traction difference `Δτ = (τ_j - τ_i)·t̂` (positive pulls the ends
  apart, negative pushes them together) times a coupling length `L_drag` adds
  `Δτ L_drag (t̂ ⊗ t̂)` to both end cells, weighted by edge length.
- **Slab pull:** each subducting boundary cell's slab-pull force
  (`torque.slab_pull_force_vectors_n`) divided by its boundary length is a line load `F_sp`
  (N/m) along direction `ê`, the tangent of the pull. A slab pulls its plate toward the trench,
  so the plate carries it as tension along `ê`. It reaches into the plate as
  `F_sp exp(-d / L_pull) (ê ⊗ ê)`, where `d` is the distance to that boundary cell. Contributions
  from separate trenches add as tensors, so two slabs pulling a plate apart add tension and
  their components along the same line combine by direction.
- **Collision and ridge push:** contact drag (`torque.contact_drag_force_vectors_n`) and ridge
  push (`torque.ridge_push_force_vectors_n`) enter the same way, decaying with `L_pull`, with
  negative (compressive) sign along their force direction. This keeps colliding margins from
  rifting.
- **Plumes:** a cell with plume heat adds an isotropic tension `N_plume f_plume 𝐈`, standing in
  for spreading above a plume head. `f_plume = min(1, plume_heat_myr / step_myr)` is the share
  of the step the heat is active, as in `cratons.heated_strength`, so a countdown that ends
  mid-step contributes only its active share.

`N` is a calibrated estimate of the thin-shell extensional resultant `N_ext`. `L_drag` and `L_pull` are calibrated in D0 against the thin-shell membrane
diagnostic (`thin_shell.solve_thin_shell`), which solves the same per-node force field exactly
but costs about 1.6 s per plate. The cheap proxy runs every step; the shell solve runs only
offline for calibration (alternative B under D0).

**2. Strength, `N_eff` (N/m).** It is `N_eff = k_cap N_cap(kind, Hc, Hm, ΔT_moho) × M × W`, where:

- `N_cap` is the column's integrated extensional capacity to normal-fault yield, from the
  RFC 0010 Gate 2 yield envelope (`analysis/issue487/`, hydrostatic reference, `K0 = 0.8`). It
  is precomputed once on a grid of crust kind, Hc, Hm and Moho thermal lag, and interpolated.
  The 100 km continental anchor is 6.28 TN/m.
- `k_cap` is a single capacity scale. Gate 2 showed that the unconstrained background stress
  ratio `K0` moves the over-capacity share from 3% to 58% of cells, so the absolute scale cannot
  be derived and is calibrated instead.
- `M` is the product of the accepted inherited-weakness multipliers, clamped to `[0.12, 1]`: belt
  0.65, suture 0.55, volcano 0.65, plume heat 0.65 and fault 0.60. Plume heat applies its
  reduction over the active share of the step: `1 - 0.35 f_plume`. The `thinning` multiplier is
  dropped, because thinning now enters through Hc, Hm and `W`.
- `W = 1 - min(s_max, s_fault φ(ε̄) + s_warm φ(ε̄_warm))` is strain softening, with
  `φ(ε) = min(1, ε / ε_soft)`. `W` lies in `[1 - s_max, 1]`. `ε̄` and `ε̄_warm` are the two strain
  fields averaged over a disk of diameter `w_0` around the cell (area-weighted). This non-local
  average gives the band a physical width scale: softening cannot concentrate into a zone
  narrower than about `w_0`, whatever the node spacing. Local softening alone would make band
  width follow the grid.

**3. Strain.** A cell strains where its load exceeds its strength:

```text
dε/dt = ε̇_0 · max(0, N / N_eff - 1)^n        [1/yr]
```

`ε̇_0` sets the slow-phase rate. With a band about 125 km wide and Brune et al.'s slow phase
under 10 mm/yr full rate, `ε̇ ≈ 0.01 m/yr / 1.25e5 m = 8e-8 /yr`. As strain softens a cell,
`N_eff` falls and `dε/dt` rises, which is the slow-to-fast transition. Strain increments
`rift_strain` and `rift_warm_strain` equally; `rift_warm_strain` then decays as
`d(warm)/dt = -warm / τ_heal`.

Integration is timestep-robust. The rate depends nonlinearly on `ε` through `W`, so each step
sub-steps the cell's update with `Δε ≤ Δε_max` per sub-step, and the decay uses its exact
exponential. A cell's strain over a step therefore does not depend on how the step is divided,
within the sub-step tolerance, which D3 tests.

**4. Thinning, subsidence and relief.** Pure-shear thinning over the step scales the cell's Hc
and Hm by `exp(-Δε)` [[McKenzie, 1978]](#ref-mckenzie1978), so the band subsides isostatically
into a rift valley. The model's cells cannot gain area inside a plate, so the removed volume
moves to the band's flanks. The flank-transfer helper adapts `_allocate_stretch`'s volume and
tracer rule in reverse. At margins, `_open_rift` thins pre-existing donors and gives their volume
to newly inserted cells. Here the straining cells are the donors, and the recipients are cells
within `RIFT_FLANK_ROWS` hops that did not strain this step, in proportion to their area.
Straining cells never receive, so volume is not shuffled within the band. Hc and Hm move together
with the same tracers `_allocate_stretch` moves: continental material, craton crust, mobile cover
and its continental share, and the plume-head and arc tracers. Donor and recipient elevations are
then recomputed isostatically. Volumes are conserved by construction, and the flanks rise
slightly, like rift shoulders.

Interior band cells do not erupt before breakup. Thinning stops at the floor
`RIFT_CRITICAL_THICKNESS_M` for crust and `MIN_MANTLE_LITHOSPHERE_THICKNESS_M` for mantle
lithosphere. A cell that reaches the crust floor is latched broken by setting
`rift_strain = max(rift_strain, ε_break)`. `lithosphere_plate._erupt_melted_nodes` resets a
melting land cell to the 35 km continental reference column, which would re-thicken the band, so
it is not called here. New oceanic crust appears only after breakup, when the daughters separate
and the existing margin path (`_open_rift`) fills the gap. Normal-fault
relief comes from the existing fault machinery: the spawn weight in
`faults._maybe_spawn_faults` gains a term proportional to the local strain rate, so active
bands grow sub-parallel normal faults, striking along the band's principal direction.

**5. Physical candidates (diagnostic).** A physical candidate is a connected component of cells
with `rift_strain ≥ ε_cand`. The two strain fields cannot recover a candidate's history: onset
time, onset width, time since last strain or rate history. A diagnostic recorder in the gate
scripts therefore tracks candidates, outside saved state:

- **Sampling:** at every step, record each component's cells, area, length (its longest geodesic
  extent along the principal strike), width (area divided by length), mean strain rate and the
  share of its length on inherited weakness.
- **Identity:** a component keeps the identity of the previous step's component it overlaps most
  (by area). When components merge, the merged one keeps the oldest identity and the others are
  recorded as merged. When one splits, the largest piece keeps the identity and the rest start
  new ones with the parent's onset time.
- **Onset** is the first sample of an identity. **Onset width** is its width there. A candidate
  is *active* while any of its cells strain, *failed* once none has strained for `t_fail`, and
  *successful* when its path splits the plate.
- **Censoring:** candidates still active at the end of a run are censored. They are counted, but
  excluded from duration and failure-share estimates, which use survival estimates over censored
  and completed candidates. Candidates already present when a run starts (from a loaded save)
  have unknown onset and are excluded from onset statistics.

Production code needs only the two fields; the recorder runs in the D1–D3 harness.

**6. Breakthrough and split.** A cell is *broken* once `rift_strain ≥ ε_break`, including cells
latched there on reaching the crust floor (step 4). For a physical candidate with broken cells,
find the cheapest path that runs from the plate boundary through the band and back to the
boundary on the far side: two Dijkstra searches (`scipy.sparse.csgraph`) from the band's broken
cells, one toward each end of the band's long axis. Broken cells cost zero and intact cells cost
`N_eff × ds`, so the path runs through the band and bridges gaps through the weakest intact rock. Another band's
broken cells also cost zero, so rifts can link. The intact part of the path, the ligament, must
carry the band's whole load. In integrated force form, with `ds_i` the path length through cell
`i` (half the centre-to-centre distances to its two path neighbours, so adaptive resolution does
not bias the sum), the plate splits when

```text
Σ_path N_i ds_i  ≥  Σ_intact N_eff,i ds_i
```

A fully broken path has an empty right-hand side, so it splits whenever the band still carries
load (`Σ_path N_i ds_i > 0`). A fully broken band whose load has stopped stays an unsplit, failed
rift.

The path is an ordered chain of cells, each sharing an edge with the next, whose two ends are
plate-boundary cells. Remove the path cells from the plate's adjacency graph; the remaining
cells must form exactly two connected components. Each path cell then joins the component with
which it shares the most edge length, with ties going to the side whose fitted rotation it
matches best. The split is applied only if both daughters are connected after that assignment
and pass the existing hard daughter-size and ownership checks. Otherwise the band keeps
growing. The rule is deterministic. Randomness enters through the model's sampled plume events and stochastic faults,
consistent with the owner's rule to threshold endogenous causes and use probability only for
exogenous ones.

### Starting values

| Knob | Start | Source | Tuning |
|---|---:|---|---|
| `ε̇_0` | `8e-8 /yr` | Brune et al. 2016's slow phase (< 10 mm/yr) over a 125 km band | D1 search |
| `n` | `1` | Design choice | D1 search over 1–3 |
| `ε_soft` | `1.0` | Design choice; linear softening form of Huismans & Beaumont 2003 | D1 search over 0.5–2 |
| `s_fault`, `s_warm`, `s_max` | `0.40`, `0.40`, `0.88` | `s_fault` matches the accepted fault reduction; `s_max` matches the accepted 0.12 clamp | D1 search on `s_warm` |
| `τ_heal` | `62.8 Myr` | McKenzie 1978 thermal time constant | Fixed; sensitivity at 30 and 120 Myr |
| `k_cap` | `1.0` | #487 reference envelope | D1 search over 0.3–3 (the `K0` bracket) |
| `L_drag`, `L_pull` | fitted | D0 fit to thin-shell `N_ext` | D0 only |
| `N_plume` | `0.5 TN/m` | Design choice, under a tenth of the 100 km continental capacity (6.28 TN/m) | D1 search over 0–2 TN/m |
| `w_0` | `125 km` | Owner's onset width (100–150 km) | D1 resolution check |
| `ε_cand` | `0.1` | Design choice | Diagnostic only |
| `ε_break` | `2.0` | Design choice; about `β = 7.4`, which takes 35 km crust below 5 km | Fixed |
| `t_fail` | `5 Myr` | Design choice | Diagnostic only |
| `RIFT_FLANK_ROWS` | `2` | Matches `K_NEIGHBOUR_ROWS_FOR_MASS_CONSERVATION` | Fixed |
| `Δε_max` | `0.05` | Design choice | D3 timestep check |
| Inherited multipliers | accepted RFC 0018 values | Accepted design choices | Fixed |

At the default node density the line spacing is 62.5 km, so a band at the `w_0 = 125 km` scale
is about two cells wide. That is the owner's chosen onset width
([#876 decisions](https://github.com/adubey/mantle-bloom/issues/876#issuecomment-6874579897)),
set by the non-local softening average rather than by the grid. Bands widen as softened flank
cells reach their own threshold. D1 checks the width at two node densities.

Scale check: with load at twice strength and no softening, `dε/dt = ε̇_0` and reaching `ε_break`
takes 25 Myr. Softening shortens this, so
first-guess durations fall within the tens of Myr that Earth's rifts take
[[Brune et al., 2016]](#ref-brune2016). D1 tunes the knobs against the targets below; it does not
treat these first guesses as validated.

### Runtime and saves

Per plate per step, the load uses existing per-node forces and a nearest-boundary distance query.
The strength is a table interpolation, the strain update is per cell, and the components search
is linear in cells. The path search runs only for candidates with broken cells. The expected cost
is milliseconds per plate, small against a production step, but under plates-only stepping it
could be a noticeable share. D3 measures it on production settings. The thin-shell solve is not
on the step path. Saves gain two float surface fields, and old saves load with zeros.

### Diagnosis gates

Development worlds are the 20 RFC 0018 development worlds
([`analysis/issue797/`](../analysis/issue797/README.md)) and the #848 development seeds. The
held-out sample is the ten-seed revised-method sample already set aside: registry seeds 12–20
plus the first frozen reserve seed `1479220500`, none of which has been run for a revised method.
No held-out seed is used for calibration or tuning.

1. **D0: calibrate the cheap load.** On development worlds, fit `L_drag` and `L_pull` so the
   cheap tensor matches the thin-shell resultant under the same forces. First check
   manufactured cases with known answers: one-sided slab pull, two slabs pulling apart, a slab
   opposed by collision contact, and pure collision (compression only). The proxy must give the
   shell's sign (tension or compression) in each. Then, on development worlds over cells where
   the shell is loaded, require:
   - **sign agreement** (extension versus none) on at least 80% of cells;
   - a **principal-direction difference** under 30° on at least 70% of extensional cells;
   - **Spearman correlation** of `N` with `N_ext` of at least 0.5;
   - a **median ratio** `N / N_ext` between 0.5 and 2 in the top decile, by crust kind.

   If any of these fails, use alternative B: run the thin-shell solve every `K` steps
   (start `K = 5 Myr`), hold `𝐍` between solves, and measure its cost in D3.
2. **D1: tune the mechanism on development worlds.** Run 200 Myr with plates-only stepping
   (`simulate_climate_biomes = False`) across a predeclared parameter search (Latin hypercube
   over the D1 knobs, as in #313). Score each setting against the realism targets below and
   freeze the best before D2. Also record the #699 failure mode for every run: the largest
   damaged share of a plate and the number of damaged components at the time of each split.
   Repeat the frozen setting on five development seeds at half the default node density. Onset
   width and median onset-to-breakup duration must stay within 25% of the default-density
   values, or band width is following the grid.
3. **D2: realism on held-out worlds.** With the frozen knobs, compare the damage-band method with
   the production great-circle method on the ten held-out seeds for 200 Myr, from identical
   generated worlds, settings and step schedules. The method passes when it meets every realism
   target, and no more than one of ten pairs has more tiny breakups than its baseline
   (`D_s = T_baseline,s - T_rift,s < 0`; the accepted tiny-breakup definition, smaller daughter
   under 7% of parent area). The pairs must contain at least 50 actual splits across both arms
   combined, the accepted rule. If they do not, add seeds in order from the frozen reserve list,
   excluding `1479220500`, until they do. Report all target values with seed-bootstrap intervals.
4. **D3: runtime and timestep.** Keep the accepted G3 limits: median runtime within 2% of the
   paired production baseline on three held-out seeds, from controlled serial runs with profiling
   kept separate; craton loss at most 2% per Gyr; and 10–40 plates at every output. Repeat at half
   the step size on three predeclared held-out seeds. Plate count and craton loss must stay
   within 20% of the normal-step results, as accepted. So must the tiny-breakup share; when the
   normal-step share is zero, the half-step run must also have zero tiny breakups. Breakup count
   and median rift duration are added, also within 20%.

**Realism targets** (median over physical candidates unless stated):

| Target | Pass range | Basis |
|---|---|---|
| Onset width (first time `ε ≥ ε_cand`) | 100–150 km | Owner design choice |
| Width growth | Successful bands are wider at breakup than at onset | Owner design choice |
| Onset-to-breakup duration | 10–50 Myr | Tens of Myr [[Brune et al., 2016]](#ref-brune2016); range is a design choice |
| Slow-to-fast transition | Strain rate in the last quarter before breakup at least 3 times the first quarter | Abrupt acceleration [[Brune et al., 2016]](#ref-brune2016); factor is a design choice |
| Breakups along inherited weakness | At least two thirds of successful paths have most of their length on belts, sutures, faults or earlier rift strain, where such marks exist on the plate | 17 of 18 margins [[Buiter & Torsvik, 2014]](#ref-buiter2014); threshold lowered for the model's incomplete suture signal |
| Failed rifts | Between 20% and 80% of physical candidates fail | Failure is common [[Brune et al., 2023]](#ref-brune2023) [[Elling et al., 2022]](#ref-elling2022); range is a design choice |
| No super seas | No below-sea-level basin made by interior rift thinning is wider than three onset band widths (about 375 km) across its narrow axis, outside crust thicker and hotter than the plate median | Narrow rifts dominate; wide rifts need thick, hot crust [[Buck, 1991]](#ref-buck1991); #866 |
| Localization | At each split, damaged cells form at most three components above 1% of plate area, and the damaged share stays under 30% of plate area | #699's fragmentation failure mode; thresholds are design choices |
| Plume link | Reported: share of breakups with plume heat on the path | 12 of 18 margins [[Buiter & Torsvik, 2014]](#ref-buiter2014); descriptive only |
| Splits happen | At least 50 actual splits across both arms combined, with the frozen reserve-list fallback | Accepted event-count rule |

**Pictures.** For five held-out seeds, show map sequences of one successful and one failed rift,
from onset to breakup or to failure: land, water, cratons, belts, rift strain and plate outlines,
in a shared projection and legend.

### Acceptance criteria

- [ ] D0 records the fitted `L_drag` and `L_pull` with their agreement statistics, or the switch
      to alternative B with its cadence.
- [ ] D1 records the parameter search, the scores and the frozen knob values, using development
      seeds only.
- [ ] D2 passes on the ten held-out seeds: every realism target is met, no more than one pair
      has `D_s < 0`, and both arms together have at least 50 actual splits, with the reserve-list
      fallback.
- [ ] D3 passes: median runtime within 2% of baseline, craton loss at most 2% per Gyr, 10–40
      plates at every output, and half-step plate count, craton loss, tiny-breakup share (with
      the zero-share rule), breakup count and rift duration within 20%.
- [ ] Every applied split makes exactly two connected plates, assigns each parent cell to one
      daughter, and preserves crust and mantle ledgers. Interior thinning conserves crust,
      continental material and tracer volumes.
- [ ] Old saves load with zero rift strain and continue stepping (`test_persistence.py`).
- [ ] Unit tests cover the bounds of the strength factors, softening monotonicity, timestep-robust
      strain integration, flank volume and tracer conservation, healing, the broken latch (no
      interior eruption), plume step fractions, directional load signs (compression never
      registers as extension), the integrated ligament rule, a fully broken edge-to-edge band that
      splits, path validity and two-daughter splits.
- [ ] The pictures above are in the gate report.
- [ ] `docs/simulation-model/tectonics.md` documents rift damage, bands and breakup, and
      `docs/architecture.md` lists any new module.

### Falsification criteria

- **Softening does not localize.** If, across the D1 search, damage still fragments plates (the
  localization target fails in most runs at every softening setting), the premise that softening
  alone localizes bands is wrong. Add a tip-growth rule (strain concentrated at band ends, as in
  RFC 0010 Gate 3a) before D2, or reject the method.
- **The cheap load cannot stand in for stress, and the shell solve is too costly.** If D0 falls
  back to alternative B and D3 shows that even a coarse cadence exceeds the runtime limit, revise
  the load model.
- **No knob setting looks right.** If no D1 setting meets most realism targets on development
  worlds, the mechanism is missing a process. Diagnose which target fails before adding knobs.
- **Rifts ignore weakness.** If successful paths follow inherited weakness no more often than
  matched random edge-to-edge paths, the strength multipliers are not doing their job. Revisit
  them.
- **Timestep dependence.** If D3's half-step comparison fails, make the strain integration
  stable or reject the design.
- **Craton safety or bookkeeping fails.** If craton loss exceeds 2% per Gyr or a ledger does not
  close, stop and diagnose before changing a bound.

## Detailed design (provisional touchpoints)

- `backend/app/rift_damage.py` (new): load proxy, strength table, strain update, healing,
  flank transfer, physical-candidate diagnostics, breakthrough path and ligament rule.
- `backend/app/surface_fields.py`: add `rift_strain` and `rift_warm_strain` to
  `SURFACE_FIELDS` as intensive float fields with zero default.
- `backend/app/merge_split.py`: `maybe_split_plate` keeps its age and cooldown gates and calls
  `rift_damage` for the breakthrough decision and path, replacing the great-circle cut, its flow
  gates and the 0.55 success roll. Remove its `apply_failed_rift` call: failed rifts are bands
  that stop growing.
- `backend/app/plates.py` and `backend/app/sparse_quad_patch.py`: the mask-based split helper
  around `_plates_from_node_masks()`, as in the accepted design.
- `backend/app/quad_tectonics.py`: reuse the `_allocate_stretch` volume rule for flank transfer.
  Interior band cells do not call `lithosphere_plate._erupt_melted_nodes` (see step 4).
- `backend/app/torque.py`: reuse `basal_drag_force_vectors_n`, `slab_pull_force_vectors_n`,
  `contact_drag_force_vectors_n` and `ridge_push_force_vectors_n` for the load tensor; no change
  to plate-motion forces.
- `backend/app/faults.py`: add a strain-rate term to the spawn weight in `_maybe_spawn_faults`.
- `backend/app/cratons.py`: change `record_rift_cut` to take the path's band mask instead of a
  cut normal; belt marking already accepts a mask.
- Strength table: precomputed from the Gate 2 envelope (`bin/debug/diagnose_gate2_column_envelope.py`)
  as a small checked-in data file, or computed at import; decide in implementation.
- `bin/debug/` (new scripts): D0 calibration against `thin_shell`, D1 search driver, realism-target
  scorer and map sequences. `bin/debug/compare_rift_load_maps.py` is the starting point for D0.
- Tests: `test_rift_damage.py` (new), `test_merge_split.py`, `test_sparse_quad_patch.py`,
  `test_cratons.py`, `test_faults.py`, `test_persistence.py`.
- Docs: `docs/simulation-model/tectonics.md`, `docs/architecture.md`.

## Assumptions and uncertainty

- **The load proxy is a calibrated stand-in.** Drag and pull are scalar estimates of the largest
  extensional resultant, not a stress tensor. Their orientation is implicit. D0 tests them
  against the thin-shell diagnostic, which itself assumes linear elasticity about a zero-load
  reference (#487).
- **Capacity scale.** `k_cap` absorbs the unconstrained background stress, pore pressure and
  rheology (#487 Gate 2), so the absolute threshold is a calibration, not a prediction.
- **Softening.** The linear softening form follows numerical rift models; its amount and strain
  scale are design choices tuned for appearance and timing.
- **Fixed-grid extension.** Cells inside a plate cannot gain area, so moving band volume to its
  flanks stands in for the band widening. It keeps volumes exact but places thickening on the
  shoulders. An alternative is to hold the volume until breakup and release it at the new margins.
- **Healing.** One thermal time constant stands in for lithospheric cooling and strengthening.
  The permanent part stands in for fault damage.
- **Plume load.** `N_plume` is a design choice; plume weakening already exists (RFC 0007).
- **Sutures.** Generated pseudo-sutures are still #770 work. Until a per-cell suture field lands,
  belts carry inherited weakness, so the inherited-weakness target is lowered.
- **Oceanic columns.** About 7 km of oceanic crust reaches the `RIFT_CRITICAL_THICKNESS_M` floor
  and latches broken after a strain of only about 0.34, against about 1.95 for 35 km continental
  crust. Oceanic plates stay
  hard to break only through their higher mature-oceanic capacity in the #487 envelope. D1 reports
  oceanic and continental breakups separately; if oceanic interiors break far more often, use the
  strain threshold `ε_break` alone for oceanic cells.
- **Directional proxy.** The load tensor sums uniaxial contributions with decay lengths. It is
  not an equilibrium stress field. D0's sign, direction and magnitude checks against the
  thin-shell solve decide whether it is good enough or periodic shell solves are needed.
- **Width scale.** The non-local softening disk `w_0` sets band width directly, which is a
  regularization choice, not a prediction. D1 checks that width does not follow the node spacing.
- **Not geological rates.** Realism targets judge appearance and timing on model worlds. They are
  not calibrated rift-frequency estimates.

## References

- <a id="ref-brune2016"></a>Brune, S., Williams, S. E., Butterworth, N. P., & Müller, R. D.
  (2016). Abrupt plate accelerations shape rifted continental margins. *Nature*, 536, 201–204.
  [doi:10.1038/nature18319](https://doi.org/10.1038/nature18319)
- <a id="ref-brune2023"></a>Brune, S., Kolawole, F., Olive, J.-A., et al. (2023).
  Geodynamics of continental rift initiation and evolution. *Nature Reviews Earth & Environment*,
  4, 235–253. [doi:10.1038/s43017-023-00391-3](https://doi.org/10.1038/s43017-023-00391-3)
- <a id="ref-buck1991"></a>Buck, W. R. (1991). Modes of continental lithospheric extension.
  *Journal of Geophysical Research*, 96(B12), 20161–20178.
  [doi:10.1029/91JB01485](https://doi.org/10.1029/91JB01485)
- <a id="ref-buiter2014"></a>Buiter, S. J. H., & Torsvik, T. H. (2014). A review of Wilson Cycle
  plate margins: A role for mantle plumes in continental break-up along sutures? *Gondwana
  Research*, 26, 627–653. [doi:10.1016/j.gr.2014.02.007](https://doi.org/10.1016/j.gr.2014.02.007)
- <a id="ref-courtillot1999"></a>Courtillot, V., Jaupart, C., Manighetti, I., Tapponnier, P., &
  Besse, J. (1999). On causal links between flood basalts and continental breakup. *Earth and
  Planetary Science Letters*, 166, 177–195.
  [doi:10.1016/S0012-821X(98)00282-9](https://doi.org/10.1016/S0012-821X(98)00282-9)
- <a id="ref-elling2022"></a>Elling, R., Stein, S., Stein, C. A., & Gefeke, K. (2022).
  Three major failed rifts in central North America: Similarities and differences. *GSA Today*,
  32(6), 4–11. [doi:10.1130/GSATG518A.1](https://doi.org/10.1130/GSATG518A.1)
- <a id="ref-huismans2003"></a>Huismans, R. S., & Beaumont, C. (2003). Symmetric and asymmetric
  lithospheric extension: Relative effects of frictional-plastic and viscous strain softening.
  *Journal of Geophysical Research*, 108(B10), 2496.
  [doi:10.1029/2002JB002026](https://doi.org/10.1029/2002JB002026)
- <a id="ref-mckenzie1978"></a>McKenzie, D. (1978). Some remarks on the development of
  sedimentary basins. *Earth and Planetary Science Letters*, 40, 25–32.
  [doi:10.1016/0012-821X(78)90071-7](https://doi.org/10.1016/0012-821X(78)90071-7)
- <a id="ref-vauchez1997"></a>Vauchez, A., Barruol, G., & Tommasi, A. (1997). Why do
  continents break-up parallel to ancient orogenic belts? *Terra Nova*, 9, 62–66.
  [doi:10.1111/j.1365-3121.1997.tb00003.x](https://doi.org/10.1111/j.1365-3121.1997.tb00003.x)
