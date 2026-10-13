# RFC 0003: Rifts go around cratons

- **Status:** Draft
- **Date:** 2026-10-13
- **Related:** #274 (craton rift and stretch resistance), #340; RFCs
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
  thins a tapered band 2.5 cells wide by up to `FAILED_RIFT_THINNING_FRACTION` = 10%. It
  applies **no craton resistance at all**: the craton thins in proportion
  (`cratons.thin_with_column`) and books to `rifted_m3`.

In the run behind RFCs 0001–0006, rifting was the second-largest craton sink. The run is
seed 456862899, 19.2 Gyr, saved locally as `mantle-bloom-seed456862899-19191800000y.mbworld`;
the file is on the author's machine and **not available for download**. Numbers come from
`bin/debug/ledger_history.py` on that save.

- `rifted_m3` totalled 555 M km³, **21%** of all craton destruction.
- It took 6%, 4%, 6%, 5% and 10% of the live craton per Gyr over the first five Gyr.
- It was steady: no single 10 Myr sample took more than 1% of the live craton.

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

On the user's question of proto-cratons: young collision belts are the opposite of
rift-resistant. They are the weak zones breakup follows. Rift resistance should come with
maturity, which RFC 0004 handles. A newly formed belt shouldn't get it.

1. **Steer the split cut onto weak lithosphere.** Among great circles close to the
   flow-derived cut, choose the one that crosses the least lithospheric strength, as long
   as it still separates the two velocity regimes. Strength combines craton strength with a
   lithosphere-thickness term, so thin belts and old sutures are preferred, after
   [[Audet & Bürgmann, 2011]](#ref-audet2011). This removes most rifts through cratons at
   their source, at a cost of a few dot products per split.
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

**Recommendation:** 1 and 2 together. Re-measure, and use 3 only if steady loss remains
above target. 4 comes with RFC 0007.

## Design

**Strength proxy.** Per node, in units of the reference continental lithosphere:

`S = cratons.strength(craton_crust_m) + max(0, Hm / REFERENCE_HM_CONTINENTAL_M − 1)`

- The first term is 0–1 (full at 25 km of craton).
- The second adds strength for keels thicker than the 100 km continental reference.
- A thinned belt (Hm below reference) or oceanic lithosphere scores 0.
- S is dimensionless, built from thicknesses in metres.

**Candidate cuts.** Start from the refined normal n₀ of `_refine_split_cut`.
- Draw `SPLIT_CUT_CANDIDATES` normals in a cone of half-angle `SPLIT_CUT_MAX_TILT_DEG`
  around n₀, on a deterministic Fibonacci spiral (no RNG, so replays are unchanged per
  seed).
- Score each by the area-weighted strength in its band, `Σ S·A` over nodes with
  |p·n| < sin(`FAILED_RIFT_BAND_MULT` × spacing). This is the same band a failed rift
  would thin.
- Reject a candidate that puts more than `SPLIT_CUT_MAX_MISFIT` of the plate's area on
  the other side from the velocity clustering. A cut that no longer separates the flow
  regimes is not the rift the flow asked for.
- Pick the lowest score. Ties keep n₀.

The cut then goes through `Plate.split` exactly as now, and the min-size and
pole-separation gates are re-checked with the new halves.

**Failed rift (2).** `factor = 1 − FAILED_RIFT_THINNING_FRACTION × taper × (1 −
CRATON_RIFT_RESISTANCE × strength)`. Units unchanged: a fraction of Hc and Hm.

**Runtime.** Splits are attempted at most once per plate per step, and only after
`SPLIT_MIN_AGE_STEPS`. With 64 candidates × ~10⁵ nodes, scoring is one 64 × N matrix
product, about 10 ms per attempt and nothing on steps without one.

**Save compatibility.** No new state; only behaviour changes.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `SPLIT_CUT_MAX_TILT_DEG` | 20° | Geometric: tilting the cut plane by 20° moves the cut by up to ~2,200 km, more than the width of most cratons (1,000–2,000 km), so a craton straddling the flow cut can usually be missed | Sweep 10–30° |
| `SPLIT_CUT_CANDIDATES` | 64 | Resolution of the cone search (about 4.5° spacing at 20°) | Fixed |
| `SPLIT_CUT_MAX_MISFIT` | 10% of plate area | Keeps the cut faithful to the flow regimes | Sweep 5–20% |
| `CRATON_RIFT_RESISTANCE`, `CRATON_STRETCH_RESISTANCE` | 0.75 (unchanged) | #274 | Raise only if solution 3 is needed |

Tune with `bin/debug/compare_cratons.py` on seeds 7, 11 and 2024 to 1 Gyr, plus a
regenerated seed 456862899 to 3 Gyr; run the replays one at a time. Watch:
- `rifted_m3` per Gyr;
- splits and failed rifts per Gyr, and plate count;
- the fraction of split cuts whose band holds craton.

The tilt and misfit sweeps pick the smallest values that meet the acceptance criteria, so
the flow keeps as much control of breakup as possible.

## Detailed design

- `backend/app/merge_split.py`
  - New `_steer_cut_off_strong_lithosphere(points, areas, strength, mask_a, cut_normal,
    band_sin)`, returning `(mask_a, cut_normal)`.
  - `maybe_split_plate` calls it after `_refine_split_cut`, with strength from the plate's
    `craton_crust_m` and `mantle_lithosphere_thickness_m`.
  - New constants `SPLIT_CUT_MAX_TILT_DEG`, `SPLIT_CUT_CANDIDATES` and
    `SPLIT_CUT_MAX_MISFIT`, with comments citing this RFC.
  - The failed-rift branch passes the same steered normal to `apply_failed_rift`.
- `backend/app/sparse_quad_patch.py`, `apply_failed_rift`: scale `factor` by craton
  strength (read `craton_crust_m`), for both Hc and Hm.
- `backend/app/cratons.py`: update the Resistance paragraph of the module docstring.
- Tests:
  - `backend/unit_tests/test_merge_split.py`: a synthetic plate whose flow-derived cut
    crosses a craton patch gets a steered cut that misses it when a candidate within the
    tilt does. With no craton, the cut is unchanged. When every candidate breaks the misfit
    limit, n₀ is kept.
  - `backend/unit_tests/test_sparse_quad_patch.py`, next to the existing
    `apply_failed_rift` tests: a full-strength craton in the band thins by 25% of the non-cratonic amount,
    and `rifted_m3` books exactly the craton lost.
- Docs: `docs/simulation-model.md`, the split description in "Merge and split" and the
  Cratons section ("Resistance").

## Acceptance criteria

- [ ] Craton loss to rifting is at most **2% of live craton per Gyr** over the first 3 Gyr
      on seed 456862899 (4–6% today), and on average across seeds 7, 11 and 2024.
- [ ] The share of split cuts whose band crosses craton drops by at least 75% against
      `origin/main` on the same seeds.
- [ ] Splits per Gyr and mean plate count change by less than 10%.
- [ ] Failed rifts thin cratonic columns by `1 − CRATON_RIFT_RESISTANCE × strength` of the
      normal amount (unit test).
- [ ] Ledgers balance; replays with the same seed are deterministic.
- [ ] Unit tests above pass; affected tests, ruff and `bin/quality_check.sh` are green.

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
- <a id="ref-tommasi2001"></a>Tommasi, A., & Vauchez, A. (2001). Continental rifting
  parallel to ancient collisional belts: An effect of the mechanical anisotropy of the
  lithospheric mantle. *Earth and Planetary Science Letters*, 185, 199–210.
  [doi:10.1016/S0012-821X(00)00350-2](https://doi.org/10.1016/S0012-821X(00)00350-2)
- <a id="ref-wilson1966"></a>Wilson, J. T. (1966). Did the Atlantic close and then
  re-open? *Nature*, 211, 676–681. [doi:10.1038/211676a0](https://doi.org/10.1038/211676a0)
