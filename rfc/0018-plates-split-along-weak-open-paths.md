# RFC 0018: Plates split along weak paths that are pulling apart

- **Status:** Draft
- **Date:** 2026-11-03
- **Related:** #770; RFCs [0009](0009-craton-mosaic.md) (cratons-first worlds), [0010](0010-rifts-follow-weak-lithosphere.md) (current accepted rift proposal; proposed replacement)

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
rift around cratons. The new proposal is narrower: keep the current split trigger and replace
only how it chooses the dividing line. This RFC is intended to supersede RFC 0010 if its
gates pass and the owner accepts it. Until then, RFC 0010 remains the accepted proposal.

## Proposed solutions

**What nature tells us.** Rifts can reuse old weak zones, but their paths and growth depend on
the direction and amount of extension; some stop before the continent breaks apart
[[Vauchez et al., 1997]](#ref-vauchez1997) [[Brune et al., 2023]](#ref-brune2023). This
supports considering both weak rock and how the two sides would move. It does not support a
particular graph-cut method or a mass exponent.

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

The production split trigger, minimum plate age, cooldown and existing daughter motion checks
stay in place for the first implementation. The initial cratons-first world, including land,
water, plate ownership, crust and mantle thickness, is not changed. The algorithm reads the
existing Hc and Hm fields to estimate column mass; it does not change Hm. The successful
split, failed-rift thinning and craton-loss diagnostics must all follow the same graph path,
so no part of the event keeps using a different great-circle line.

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
3. Sample weak edges that connect different motion groups. For each sampled edge, find a
   low-cost cut that crosses it. The cut cost combines structural resistance with a
   candidate-specific motion-fit term: compare each cell's velocity with the fitted rotation
   of the daughter side that contains it. A good fit means the cell resists breaking away
   from that group; a poor fit means it offers less resistance. Since the term depends on
   the proposed cut, evaluate it for each candidate partition rather than storing one
   motion-resistance value per cell for all cuts. Keep the candidate only if each side is
   connected and one dominant cut trace reaches the plate edge at two places. At least 90%
   of the cut length must belong to that trace; report the remaining fragments. Use a
   separate, deterministic sample keyed by world seed, plate and edge so this choice does
   not consume random numbers used by unrelated simulation processes.
4. Measure signed opening along each candidate cut, using shared-edge length as the weight.
   Reject candidates whose net opening is closing, or whose opening is no larger than the
   fit's estimated uncertainty. Also report the share of cut length that opens.
5. Rank the remaining candidates by their potential momentum score, then split on the best
   candidate if it passes the existing production checks. The score ranks options; it is
   not a claim that the simulation conserves or releases physical momentum at the cut.

For each cell `n`, estimate its mass as cell area times the mass per area of its crust and
mantle lithosphere: `m_n = A_n * (rho_c * Hc_n + rho_m * Hm_n)` (kg). Let `v_new(n)` and
`v_old(n)` be its velocity under the fitted daughter and parent rotations (m/s). The
potential momentum score before the balance factor is:

```text
P = sum over cells n of m_n * |v_new(n) - v_old(n)|   [kg m/s]
```

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
mass, while keeping outward opening as a separate required check. The held-out gate compares
both scores on tuning seeds, so the additional balance factor must earn its place before the
winning version is frozen for held-out runs.

### Starting values and limits

The first tuning pass uses the prototype's 12 sampled cuts per eligible plate and its
dimensionless resistance rules. They are starting points only:

| Signal | First trial | Meaning |
|---|---:|---|
| Column mass | `0.6 + 0.4 * clamp(mass / plate-median, 0, 4)` | A thicker or denser column costs more to cut |
| Craton | Multiply resistance by `2.5` | A craton is harder to cut, not an absolute no-go zone |
| Motion fit | Divide resistance by `1 + clamp(residual / median-residual, 0, 4)` | For each candidate, residual is measured against the fitted rotation of the daughter side containing the cell; a poor fit lowers resistance to breaking from that side |
| Belt / suture | Reduce resistance by `0.35` / `0.45` | Inherited weak marks lower the cost |
| Volcano / plume heat | Reduce resistance by `0.35` each | Existing volcanic or heated cells lower the cost |
| Fault / thinning | Reduce resistance by `0.40` / `0.25` | Older faults count; thinner columns lower the cost |
| Weakening limit | Combined reduction clamped to `0.12–1.0` | Avoid zero or negative resistance |
| Resistance floor | `0.02` | Avoid zero edge cost |
| Candidate count | `12` | Number of sampled weak edges; tune only on development seeds |
| Balance exponent | `mu = 1`; sweep `0, 0.5, 1, 2` | Design choice; not a measured physical value |

The motion-residual reduction is candidate- and assignment-specific: the fit is measured
against the daughter group a cell would remain with, never against whichever of the two
groups happens to fit it better. The edge-cost floor from the prototype is retained for the
first tuning pass. All score factors stay bounded. Area is used when converting column
mass to a plate or daughter total; raw node count is not a substitute for mass.

Keep the existing 0.55 split-success chance for the first paired comparison so this change
isolates the cut choice. A failed split thins the existing failed-rift band along the chosen
graph path. A separate revision would be needed to remove that chance or change the failed
rift response.

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
   share, cut length, trace continuity, net opening, opening uncertainty, score and runtime.
   Use area-derived column mass. Keep every failed world-generation seed in the report; add
   only predeclared replacements to reach 20 constructible baseline worlds, following the
   #770 census practice. Use fit residuals to define and freeze an opening floor so that a
   tiny positive value like the pilot's +0.013 cm/yr does not pass merely because it is
   above zero. Compare the potential-momentum score with the prototype's reduced-mass-times-
   opening score, and compare `mu = 1` with `mu = 0`; tune only on this set, then freeze the
   winning score and weights.
2. **G1 — Check that the candidate can make a usable split.** Require that, in the frozen
   development run, at least 10 plate opportunities are observed and at least two thirds
   of opportunities that pass the current split
   trigger produce at least one candidate that is connected on both sides, has a dominant
   edge-to-edge trace, passes existing daughter checks and has outward opening above its
   uncertainty. If not, stop and revise the graph construction or opening measure before a
   long replay. Report separately the failures from minimum daughter size, pole separation,
   trace breaks, closing motion and the stochastic success chance. The three-world pilot
   passed this only once in three worlds and therefore does not pass G1.
3. **G2 — Compare actual breakups on held-out worlds.** Freeze a new 20-seed list before
   running. Compare the production great-circle method and graph-cut method from identical
   generated worlds, settings, step schedule and random streams for 200 Myr. Use the same
   tiny-breakup measure agreed in #770: the share of actual breakups whose smaller daughter
   has less than 7% of the parent plate's area. The graph method must reduce that share by at
   least one third overall and improve at least two seeds. The comparison must include at
   least 50 actual split events across the seed set; if it does not, continue only with seeds
   from a reserve list frozen before the run. Report counts and per-seed results;
   do not infer success from a favorable image or from a routing-only proxy. Report the
   smaller daughter mass distribution too, but do not impose a hard 30% floor.
4. **G3 — Check cratons, stability and cost.** Preserve the RFC 0010 limits unless the
   owner changes them: rifting-related craton loss at or below 2% of live craton per Gyr,
   and 10–40 plates at every output. Repeat representative held-out runs at half the step
   size on three predeclared held-out seeds; plate counts, tiny-breakup share and craton loss
   must stay within 20% of the normal step results. Median full-simulation runtime must remain within 2% of the paired current
   production baseline. Record any seed where the model cannot finish or maintain its
   ledgers.

### Acceptance criteria

- [ ] G0 and G1 pass on the 20 development seeds, with all seed failures and rejection
      reasons kept in the report.
- [ ] G2 passes on 20 separate, predeclared held-out seeds: at least one-third fewer
      breakups with a smaller daughter under 7% of parent area overall, with improvement on
      at least two seeds.
- [ ] Every applied split makes exactly two connected plates, assigns each parent cell to
      one daughter, and preserves crust and mantle ledgers. Report daughter mass shares and
      cut trace fragments; no fixed 30% mass minimum is added.
- [ ] G3 passes: craton loss at most 2% per Gyr, plate count 10–40 at every output,
      timestep metrics within 20%, and median runtime within 2% of baseline.
- [ ] Old saves load and can continue stepping. The method adds no saved fields unless a
      later revision explains why they are needed and how old saves initialize them.
- [ ] Unit tests cover connected daughters, boundary-to-boundary traces, positive opening,
      balanced-score ordering, area-weighted mass, failed-rift thinning along the graph path,
      exact ownership and old-save loading.
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
- **Mass weighting does not reduce tiny breakups.** If G2 fails after the weights and `mu`
  have been frozen on development seeds, reject this score. Do not tune on the held-out seeds.
- **The split score rewards closing or broken cuts.** If high-scoring cuts often close, break
  into disconnected traces, fail daughter checks or create extra fragments, the score is
  incomplete. Revise it before a production proposal.
- **The result depends strongly on step size.** If G3's 20% limit fails, make the candidate
  sampling and event timing stable across timestep changes or reject the design.
- **Craton safety or bookkeeping fails.** If the 2% craton-loss bound, connected ownership,
  or material ledgers fail, stop implementation and diagnose the path before changing the
  bound.

## Detailed design (provisional touchpoints)

- `backend/app/merge_split.py`: replace great-circle cut selection inside
  `maybe_split_plate`; keep the existing trigger, cooldown and success chance initially.
  Route failed-rift thinning through the selected cut rather than the old great circle.
- `backend/app/plates.py`: use `adjacency()`, `surface_nodes().area_m2` and
  `_plates_from_node_masks()` to score cell edges and build the two daughters.
- `backend/app/sparse_quad_patch.py`: let failed-rift thinning follow the graph trace while
  keeping the existing material accounting.
- `backend/app/cratons.py`: keep the existing rift-loss diagnostic, but measure the graph
  trace rather than a great-circle band.
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
- Fitted daughter rotations approximate the motion after separation well enough to screen
  outward opening. This needs held-out replay evidence.
- A mass-weighted velocity change can rank split options. It is a heuristic score, not
  angular momentum, a force law or a conservation ledger.
- A balance exponent may reduce very small daughters, but no evidence establishes its value.
  Its value is chosen on development seeds and judged only on separate held-out seeds.
- The #770 routing results establish that lowering suture route cost alone did not reduce
  tiny breakups in that router. They do not test this graph-cut selection method.
- The three applied graph-cut splits are an exploratory demonstration, not a success-rate
  estimate. One had only +0.013 cm/yr net opening and two made no valid split.

## References

- <a id="ref-vauchez1997"></a>Vauchez, A., Barruol, G., & Tommasi, A. (1997). Why do
  continents break-up parallel to ancient orogenic belts? *Terra Nova*, 9, 62–66.
  [doi:10.1111/j.1365-3121.1997.tb00003.x](https://doi.org/10.1111/j.1365-3121.1997.tb00003.x)
- <a id="ref-brune2023"></a>Brune, S., Kolawole, F., Olive, J.-A., et al. (2023).
  Geodynamics of continental rift initiation and evolution. *Nature Reviews Earth & Environment*,
  4, 235–253. [doi:10.1038/s43017-023-00391-3](https://doi.org/10.1038/s43017-023-00391-3)
