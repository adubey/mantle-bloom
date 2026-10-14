# RFC 0003: Rifts go around cratons

- **Status:** Draft
- **Date:** 2026-10-13
- **Related:** #274 (craton rift and stretch resistance), #340, #378 (phase 1 gate); RFCs
  [0002](0002-craton-trench-attrition.md), [0004](0004-proto-cratons.md)

## Problem

A plate splits when one rigid rotation explains its mantle flow poorly
(`merge_split.maybe_split_plate`). The velocity field is clustered into two regimes, and
the plate is cut along the great circle equidistant from the two clusters' centroids
(`_refine_split_cut`). Nothing in that choice looks at the lithosphere: the cut goes
wherever the flow puts it, straight through a craton if one is there. Then:

- **If the rift succeeds** (`RIFT_SUCCESS_PROBABILITY` = 0.55), the cut becomes a divergent
  boundary. Divergent thinning in `lithosphere_plate.deform_columns` is cut by
  `CRATON_RIFT_RESISTANCE` (0.75 at full strength), so a craton on the new boundary still
  thins at a quarter of the normal rate. Rifted footprint taken from cratonic donors is cut
  the same way (`CRATON_STRETCH_RESISTANCE`, 0.75), in `quad_tectonics._allocate_stretch`.
  Over many steps the boundary still eats through the craton.
- **If the rift fails** (45% of attempts), `apply_failed_rift` (`sparse_quad_patch.py`)
  thins a tapered band reaching `FAILED_RIFT_BAND_MULT` = 2.5 cells either side of the cut
  by up to `FAILED_RIFT_THINNING_FRACTION` = 10%. It
  applies **no craton resistance at all**: the craton thins in proportion
  (`cratons.thin_with_column`) and books to `rifted_m3`.

In the run behind RFCs 0001–0007, rifting was the second-largest craton sink. The run is
seed 456862899, 19.2 Gyr, saved locally as `mantle-bloom-seed456862899-19191800000y.mbworld`;
the file is on the author's machine and **not available for download**. Numbers come from
`bin/debug/ledger_history.py` on that save.

- `rifted_m3` totalled 555 M km³, **21%** of all craton destruction.
- It took 6%, 4%, 6%, 5% and 10% of the live craton per Gyr over the first five Gyr.
- It was steady: no single 10 Myr sample took more than 1% of the live craton.

The run is a stress test of the model's long-term behaviour, not a calibration interval. It
shows that rifting wears cratons down steadily. It doesn't say how often cratons should
rift, and `rifted_m3` lumps three code paths (divergent thinning, stretch donation and
failed rifts), so it doesn't say which path matters most either.

## Proposed solutions

**What happens in nature.** Rifts localise in weak lithosphere and go around strong blocks.
- The East African Rift splits into two branches around the Tanzania craton. Numerical
  models reproduce this from a single plume beneath the craton, whose keel deflects the
  melt and strain to its edges [[Koptev et al., 2015]](#ref-koptev2015).
- Continental breakup follows old orogenic belts again and again. This is the Wilson cycle:
  the Atlantic reopened along the Caledonian–Appalachian suture
  [[Wilson, 1966]](#ref-wilson1966). Two explanations are offered, the mechanical
  anisotropy those belts leave in the mantle lithosphere
  [[Tommasi & Vauchez, 2001]](#ref-tommasi2001) and their thinner, weaker lithosphere
  [[Audet & Bürgmann, 2011]](#ref-audet2011).
- Cratons survive in part because weaker mobile belts around them take the strain
  [[Lenardic et al., 2003]](#ref-lenardic2003).
- Rifting through a craton needs the craton to be weakened first, by heating or melt
  metasomatism [[Foley, 2008]](#ref-foley2008).

Should proto-cratons resist rifting too? Mostly not: young collision belts are the
opposite of rift-resistant, since they are the weak zones breakup follows. Rift resistance
should come with maturity, once a column becomes craton (RFC 0004). A newly formed belt
shouldn't get it.

1. **Steer the split cut onto weak lithosphere.** Among great circles close to the
   flow-derived cut, choose the one that crosses the least lithospheric strength, as long
   as it still separates the two velocity regimes. Strength combines craton strength with a
   lithosphere-thickness term, so thin lithosphere is preferred, as breakup prefers it at
   old orogens [[Audet & Bürgmann, 2011]](#ref-audet2011). This removes most rifts through
   cratons at their source, for tens of ms per split attempt (Design).
2. **Give failed rifts the same craton resistance as active ones.** Scale
   `apply_failed_rift`'s thinning by `1 − CRATON_RIFT_RESISTANCE × strength`, matching
   `deform_columns`. This is a consistency fix.
3. **Raise rift resistance.** Raise `CRATON_RIFT_RESISTANCE` and
   `CRATON_STRETCH_RESISTANCE` from 0.75 toward 0.9. Cheap, but it treats the symptom: a
   boundary through a craton would still exist, only slower.
4. **Plume-style weakening.** Make craton rifting depend on a weakening state (heat or
   metasomatism), after [[Foley, 2008]](#ref-foley2008). The model has no plumes yet;
   [RFC 0007](0007-mantle-plumes-and-plateaus.md) adds them, with a `plume_heat_myr`
   field that cancels craton rift resistance and makes heated cells weak in this RFC's
   cut score.

**Recommendation:** first a gate that splits `rifted_m3` by path, so the fix targets the
path that matters. That gate has run (Design, "Phase 1 outcome"): stretch donation and
divergent thinning take 83–90% of craton rifting, failed rifts 10–17%, and only 12–18% of
cuts cross craton. So:
- do 2 now. It is a consistency fix, but it acts on at most 10–17% of the sink;
- run a second gate (phase 1b) to measure how much of the divergent and stretch loss
  happens on boundaries a craton-crossing cut created. Only that share is what 1 can remove;
- build 1 only if phase 1b says that share is large. Otherwise, 3 or a fix for boundary
  migration is the main lever.

4 comes with RFC 0007.

## Design

**Phase 1: split the rifting sink (gate).** Book `rifted_m3` in three sub-accounts:
- `rifted_divergent_m3`, from `lithosphere_plate.deform_columns`;
- `rifted_stretch_m3`, from the stretch donors in `quad_tectonics._open_rift` (their shares
  come from `_allocate_stretch`);
- `rifted_failed_m3`, from the failed-rift branch in `merge_split.maybe_split_plate`.

Record, per split and per failed rift, the craton area inside the cut's band. Run the
tuning seeds. The results decide the order of work: solution 1 helps the first two paths,
and solution 2 only the third.

**Phase 1 outcome (#378, PR #380).** The sub-accounts above, plus per-outcome counts of
cuts, cuts with craton in the failed-rift band, and the band's total and cratonic area.
They are diagnostics outside the balance and don't change the model. Seeds 7, 11 and 2024
ran to 1 Gyr at default density with `bin/debug/compare_cratons.py`; the numbers are in
`analysis/issue378/`.

| Seed | Rifted, % of start craton per Gyr | Divergent | Stretch | Failed | Splits / failed rifts with craton in band |
|---|---|---|---|---|---|
| 7 | 4.0% | 32% | 51% | 17% | 5 of 33 / 4 of 29 |
| 11 | 4.3% | 28% | 58% | 14% | 6 of 37 / 5 of 31 |
| 2024 | 10.7% | 41% | 50% | 10% | 4 of 32 / 5 of 28 |

Path columns are shares of `rifted_m3`.
- Stretch donation is the largest path, then divergent thinning. Failed rifts cost 0.4–1.0%
  of live craton per Gyr, so removing them entirely would leave 3.4–9.7%/Gyr.
- Reaching 2%/Gyr means cutting rifting loss by about half on seeds 7 and 11, by 81% on seed
  2024, and by 68% on the three-seed mean (6.3%/Gyr).
- Only 12–18% of cuts have craton in their band, and craton is 4–18% of the band area.
- What phase 1 can't say is how much divergent and stretch loss comes from boundaries those
  craton-crossing cuts created. A craton can also reach a divergent boundary without a cut:
  cratons are seeded 875 km and form 500 km inside their plate's margin, so one reaches a
  boundary only if the margin is consumed or a convergent boundary later turns divergent.
  That is the falsification case "rift cuts aren't the path", and phase 1 neither confirms
  nor rules it out.

**Phase 1b: attribute boundary loss to its origin (gate).** Before solution 1:
- At each split, mark the nodes in the cut's band of both daughters with whether the band
  held craton. Over the next 1 Gyr, book their divergent and stretch craton loss by that
  mark: cut through craton, cut clear of craton, or never in a cut band.
- Do this in a replay-only probe in `bin/debug/`, so no saved state changes. It can carry
  the mark as a surface field that exists only in the probe's process.
- Run it on seeds 7, 11 and 2024 to 1 Gyr.
- Let F be the share of `rifted_m3` that solution 1 could remove: all failed-rift loss
  (a failed rift loses craton only if its cut crosses one), plus divergent and stretch loss
  on nodes marked "cut through craton".
  - F ≥ 50%: build solution 1 with 2, as designed below.
  - 25% ≤ F < 50%: steering can't reach the target alone. Build 1 only together with 3,
    and size 3 to the remainder.
  - F < 25%: the "rift cuts aren't the path" falsification case holds. Don't build 1;
    revise this RFC toward boundary migration or solution 3.

**Strength proxy (a cut score only).** Per node, dimensionless and bounded to [0, 1]:

`S = clip(max(c, k), 0, 1)`, where

- `c = cratons.strength(craton_crust_m)`, already in [0, 1] (full at 25 km of craton);
- `k = clip((Hm / REFERENCE_HM_CONTINENTAL_M − 1) / (KEEL_FULL_RATIO − 1), 0, 1)`, rising
  from 0 at the 100 km reference to 1 at `KEEL_FULL_RATIO` × reference;
- RFC 0007's plume heat, once it exists, multiplies S by `1 − PLUME_CRATON_WEAKENING` on
  heated cells (`plume_heat_myr > 0`).

Taking the maximum rather than the sum means a cratonic, thick-keeled column scores 1, not
2. S is used **only** to rank candidate cuts; it never enters a thinning factor.

**Candidate cuts.** Start from the refined normal n₀ of `_refine_split_cut`.
- Draw `SPLIT_CUT_CANDIDATES` normals in a cone of half-angle `SPLIT_CUT_MAX_TILT_DEG`
  around n₀, on a deterministic Fibonacci spiral (no RNG, so replays are unchanged per
  seed).
- Score each by the area-weighted strength in its band, `Σ S·A` over nodes with
  |p·n| < sin(`FAILED_RIFT_BAND_MULT` × spacing). This is the same band a failed rift
  would thin.
- **Keep only cuts that still separate the flow regimes.** A candidate is admissible only
  if all of these hold:
  - it puts at most `SPLIT_CUT_MAX_MISFIT` of the plate's area on the other side from the
    velocity clustering;
  - its two halves' best-fit Euler poles are still at least `pole_separation_threshold`
    apart. That is the *relaxed* threshold `maybe_split_plate` already computes,
    `SPLIT_MIN_POLE_SEPARATION × (1 − relax)`, passed in so steered and unsteered cuts
    face the same gate;
  - the RMS residual of the two halves' fits, per node and unweighted as
    `merge_split._fit_residual_rms` computes it, is at most
    `1 + SPLIT_CUT_MAX_RESIDUAL_INCREASE` times that of n₀'s halves.

  Minimising strength alone could steer the cut into a different tectonic configuration;
  these tests stop that.
- Among admissible candidates, pick the lowest score. Ties, or no admissible candidate,
  keep n₀.

The cut then goes through `Plate.split` exactly as now, and the min-size gate is re-checked
with the new halves.

**Thinning factors, bounded.** Wherever craton resistance meets a thinning factor, the
factor is clamped:
- `resisted = clip(CRATON_RIFT_RESISTANCE × cratons.strength(craton), 0, 1)`, and with RFC
  0007, `× (1 − PLUME_CRATON_WEAKENING × heated)`, with heated = 1 while
  `plume_heat_myr > 0` and 0 otherwise, still clipped to [0, 1];
- failed rift (2): `factor = clip(1 − FAILED_RIFT_THINNING_FRACTION × taper × (1 −
  resisted), 0, 1)`;
- the same clamp in `deform_columns`, where `1 − resisted` scales the thinning delta.

Units are unchanged: a fraction of Hc and Hm.

**Batched fits.** Refitting each candidate's halves with `mantle.fit_euler_pole` builds an
N × 3 least-squares system per half. A synthetic check (in review) of 64 candidates at 10⁵
nodes put that at ~0.4 s per attempt, against ~10 ms for scoring. So the fits are batched
through moments instead:
- For unit positions p and flow velocities v, the least-squares pole for any half H of
  the nodes solves `M_H ω = b_H`, with `M_H = Σ_H (I − p pᵀ)` and `b_H = Σ_H p × v`.
- Its residual sum of squares is `Σ_H |v|² − ω · b_H`. Both identities were checked
  numerically against `fit_euler_pole` and `_fit_residual_rms` (agreement to ~1e-14).
- These are sums of per-node moments (6 of `p pᵀ`, 3 of `p × v`, 1 of `|v|²`, 1 count).
- One product of the N × 11 moment matrix with an N × K side-indicator matrix gives every
  candidate's half-A sums. Half B is the total minus half A.
- A 3 × 3 solve per half then gives the poles, separations and RMS residuals for all K
  candidates.
- *Rank handling:* a half whose `M_H` has smallest eigenvalue below
  `SPLIT_FIT_MIN_EIGENVALUE × count` (a tiny or near-degenerate half) makes the candidate
  inadmissible rather than producing an ill-conditioned pole.
- A unit test checks the batched poles and residuals against `fit_euler_pole` and
  `_fit_residual_rms` to 1e-9 relative.

**Runtime.** Splits are attempted at most once per plate per step, and only after
`SPLIT_MIN_AGE_STEPS`.
- At 64 candidates × ~10⁵ nodes, scoring and the batched moments are each one N × K
  product: ~7 × 10⁷ multiply-adds, tens of ms. Side indicators are built in chunks of 16
  candidates to keep memory near 13 MB.
- The budget is under 100 ms per attempt at 10⁵ nodes, nothing on steps without one, and
  under 2% of step time overall.
- The whole path (scoring, indicators, moments, solves) is profiled before merging, not
  just the score.

**Save compatibility.** No new state; only behaviour changes.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `SPLIT_CUT_MAX_TILT_DEG` | 20° | Geometric: tilting the cut plane by 20° moves the cut by up to ~2,200 km, more than the width of most cratons (1,000–2,000 km), so a craton straddling the flow cut can usually be missed | Sweep 10–30° |
| `SPLIT_CUT_CANDIDATES` | 64 | Resolution of the cone search (about 4.5° spacing at 20°) | Fixed |
| `SPLIT_CUT_MAX_MISFIT` | 10% of plate area | Design choice: keeps the cut faithful to the flow regimes | Sweep 5–20% |
| `SPLIT_CUT_MAX_RESIDUAL_INCREASE` | 0.10 | Design choice: a steered cut may fit the flow at most 10% worse | Sweep 0.05–0.25 |
| `SPLIT_FIT_MIN_EIGENVALUE` | 1e-3 | Numerical: rejects halves too small or too narrow to constrain a pole | Fixed |
| `KEEL_FULL_RATIO` | 2.0 (a 200 km keel) | Cratonic keels reach ~200–250 km [[Pearson et al., 2021]](#ref-pearson2021) | Fixed |
| `CRATON_RIFT_RESISTANCE`, `CRATON_STRETCH_RESISTANCE` | 0.75 (unchanged) | #274 | Raise only if solution 3 is needed |

Tune with `bin/debug/compare_cratons.py` on seeds 7, 11 and 2024 to 1 Gyr, plus a
regenerated seed 456862899 to 3 Gyr; run the replays one at a time. Watch:
- `rifted_m3` per Gyr;
- splits and failed rifts per Gyr, and plate count;
- the fraction of split cuts whose band holds craton.

The tilt and misfit sweeps pick the smallest values that meet the acceptance criteria, so
the flow keeps as much control of breakup as possible.

## Detailed design (provisional touchpoints)

- `backend/app/merge_split.py`
  - New `_steer_cut_off_strong_lithosphere(points, areas, strength, velocities, mask_a,
    cut_normal, band_sin, pole_separation_threshold)`, returning `(mask_a, cut_normal)`,
    with the admissibility tests (misfit, pole separation, residual).
  - New `_batched_half_fits(points, velocities, normals)`: the moment-based poles and RMS
    residuals for every candidate's two halves.
  - `maybe_split_plate` calls the steering after `_refine_split_cut`. It passes strength
    from the plate's `craton_crust_m` and `mantle_lithosphere_thickness_m`, and its relaxed
    `pole_separation_threshold`.
  - The failed-rift branch passes the steered normal to `apply_failed_rift`, and books
    `rifted_failed_m3` (phase 1).
  - New constants `SPLIT_CUT_MAX_TILT_DEG`, `SPLIT_CUT_CANDIDATES`,
    `SPLIT_CUT_MAX_MISFIT`, `SPLIT_CUT_MAX_RESIDUAL_INCREASE`, `SPLIT_FIT_MIN_EIGENVALUE`
    and `KEEL_FULL_RATIO`, with comments citing this RFC.
- `backend/app/sparse_quad_patch.py`, `apply_failed_rift`: scale `factor` by craton
  strength (read `craton_crust_m`), for both Hc and Hm, clamped to [0, 1].
- `backend/app/lithosphere_plate.py`, `deform_columns`: clamp `resisted` and book
  `rifted_divergent_m3`.
- `backend/app/quad_tectonics.py`, `_open_rift`: book `rifted_stretch_m3` (done, #378).
- `backend/app/cratons.py`: the three `rifted_*` sub-accounts, as diagnostics outside the
  balance (their sum is `rifted_m3`); update the Resistance paragraph of the module
  docstring.
- Tests:
  - `backend/unit_tests/test_merge_split.py`:
    - a synthetic plate whose flow-derived cut crosses a craton patch gets a steered cut
      that misses it, when a candidate within the tilt does;
    - with no craton, the cut is unchanged;
    - when every candidate fails an admissibility test (misfit, pole separation or
      residual), n₀ is kept;
    - a candidate that would merge the two flow regimes into one half is rejected even
      when it scores lowest;
    - `_batched_half_fits` matches `mantle.fit_euler_pole` and `_fit_residual_rms` per
      candidate to 1e-9 relative, and rejects a degenerate half;
    - with `relax` near 1 (a large, stressed plate), the steered cut is checked against
      the relaxed pole threshold, not the base constant;
    - S stays in [0, 1] for craton plus keel, and every thinning factor stays in [0, 1]
      for full-strength craton, zero strength and full plume weakening (a property test
      over a grid of inputs).
  - `backend/unit_tests/test_sparse_quad_patch.py`, next to the existing
    `apply_failed_rift` tests: a full-strength craton in the band thins by 25% of the
    non-cratonic amount, and `rifted_m3` books exactly the craton lost.
- Docs: `docs/simulation-model.md`, the split description in "Merge and split" and the
  Cratons section ("Resistance").

## Assumptions and uncertainty

- **Which path matters is measured, but not where its loss comes from.** Phase 1 found
  stretch and divergent loss dominant. Whether they run on boundaries born from
  craton-crossing cuts is unmeasured until phase 1b.
- **Three seeds to 1 Gyr.** Seed 2024 rifts 2.5× more than the others with a similar share
  of craton-crossing cuts, so the spread between seeds is wide.
- **Keel thickness stands in for lithospheric strength.** Real strength depends on
  temperature, composition and water content as well as thickness. The proxy is used only
  to *rank* cuts, which limits the damage if it's wrong.
  [[Audet & Bürgmann, 2011]](#ref-audet2011) support thin, weak lithosphere as the
  preferred breakup site.
- **The flow decides the rift.** Steering assumes the flow's regime boundary is right
  within ~20° and only the exact line is free. If the flow field itself is the problem
  (cuts through cratons because the convection centres sit under them), steering can't
  fix it.
- **Weakening comes later.** Without RFC 0007, nothing weakens a craton, so the model
  will under-produce the rare craton rifts that do happen.

## Acceptance criteria

- [x] **Gate (phase 1):** `rifted_m3` is split into its three paths, and on the tuning
      seeds each path's share and the craton area in cut bands are reported in
      `analysis/`, before solutions 1 and 2 are implemented (#378).
- [ ] **Gate (phase 1b):** F, the share of `rifted_m3` from craton-crossing cuts, is
      reported per seed in `analysis/`, and this RFC is updated with the branch it selects,
      before solution 1 is implemented.
- [ ] Craton loss to rifting is at most **2% of live craton per Gyr** over the first 3 Gyr
      on seed 456862899 (4–6% today), and on average across seeds 7, 11 and 2024. This
      is a model design target, not a measured Earth rate.
- [ ] The share of split cuts whose band crosses craton drops by at least 75% against
      `origin/main` on the same seeds.
- [ ] Splits per Gyr and mean plate count change by less than 10%, and the mean fit
      residual of post-split plates rises by less than 10%. Breakup still follows the
      flow.
- [ ] Failed rifts thin cratonic columns by `1 − resisted` of the normal amount, with
      every factor in [0, 1] (unit tests).
- [ ] A split attempt at 10⁵ nodes, full path, takes under 100 ms; step time changes by
      less than 2% at density 4.
- [ ] Ledgers balance; replays with the same seed are deterministic.
- [ ] Unit tests above pass; affected tests, ruff and `bin/quality_check.sh` are green.

## Falsification criteria

- **Rift cuts aren't the path.** If phase 1b shows most craton rifting is stretch donation
  or divergent thinning along boundaries that formed elsewhere and later migrated onto
  cratons (F < 25%), steering the cut won't help. Revise toward boundary migration, or
  toward solution 3. Phase 1 already found these two paths at 83–90% of the sink.
- **No admissible alternative.** If, in most splits that cross craton, no candidate within
  the tilt cone passes the admissibility tests, the flow field puts regime boundaries
  through cratons. The fix belongs in the mantle model, not the cut.
- **Steering changes the tectonics.** If plate counts, split frequency or fit residuals
  move by more than the limits above even at the smallest tilt that helps, the mechanism
  is too invasive. Fall back to solution 3.

## References

- <a id="ref-audet2011"></a>Audet, P., & Bürgmann, R. (2011). Dominant role of tectonic
  inheritance in supercontinent cycles. *Nature Geoscience*, 4, 184–187.
  [doi:10.1038/ngeo1080](https://doi.org/10.1038/ngeo1080)
- <a id="ref-foley2008"></a>Foley, S. F. (2008). Rejuvenation and erosion of the cratonic
  lithosphere. *Nature Geoscience*, 1, 503–510.
  [doi:10.1038/ngeo261](https://doi.org/10.1038/ngeo261)
- <a id="ref-koptev2015"></a>Koptev, A., Calais, E., Burov, E., Leroy, S., & Gerya, T.
  (2015). Dual continental rift systems generated by plume–lithosphere interaction.
  *Nature Geoscience*, 8, 388–392. [doi:10.1038/ngeo2401](https://doi.org/10.1038/ngeo2401)
- <a id="ref-lenardic2003"></a>Lenardic, A., Moresi, L.-N., & Mühlhaus, H. (2003).
  Longevity and stability of cratonic lithosphere: Insights from numerical simulations of
  coupled mantle convection and continental tectonics. *Journal of Geophysical Research:
  Solid Earth*, 108(B6), 2303.
  [doi:10.1029/2002JB001859](https://doi.org/10.1029/2002JB001859)
- <a id="ref-pearson2021"></a>Pearson, D. G., Scott, J. M., Liu, J., Schaeffer, A., Wang,
  L. H., van Hunen, J., Szilas, K., Chacko, T., & Kelemen, P. B. (2021). Deep continental
  roots and cratons. *Nature*, 596, 199–210.
  [doi:10.1038/s41586-021-03600-5](https://doi.org/10.1038/s41586-021-03600-5)
- <a id="ref-tommasi2001"></a>Tommasi, A., & Vauchez, A. (2001). Continental rifting
  parallel to ancient collisional belts: An effect of the mechanical anisotropy of the
  lithospheric mantle. *Earth and Planetary Science Letters*, 185, 199–210.
  [doi:10.1016/S0012-821X(00)00350-2](https://doi.org/10.1016/S0012-821X(00)00350-2)
- <a id="ref-wilson1966"></a>Wilson, J. T. (1966). Did the Atlantic close and then
  re-open? *Nature*, 211, 676–681. [doi:10.1038/211676a0](https://doi.org/10.1038/211676a0)
