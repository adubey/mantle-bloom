# RFC 0002: Cratons are rarely subducted: margin-dependent attrition

- **Status:** Draft
- **Date:** 2026-10-13
- **Related:** #340 (attrition at trenches), #177 (margin retreat budget), #320 (crust
  transfer); RFCs [0001](0001-merge-cap-craton-loss.md),
  [0005](0005-continental-recycling-balance.md)

## Problem

Since #340, a craton cell is never consumed whole at a plate boundary. Its exposed edge
wears back at a trench instead (`cratons.attrite`): each exposed craton cell loses
`CRATON_ATTRITION_KM2_PER_MYR` = 40 km²/Myr of crust per unit of edge length, which is the
cross-section a non-accreting margin loses to subduction erosion
[[Scholl & von Huene, 2007]](#ref-scholl2007).

A column of Hc = 35 km under a 40 km²/Myr loss retreats at 40/35 ≈ **1.1 km/Myr**, at
any node density. That is a realistic rate for an *erosive* margin: the Mariana forearc
retreats about 1 km/Myr and northern Chile about 3 km/Myr
[[Clift & Vannucchi, 2004]](#ref-clift2004). The problem is that the model applies it at
every trench a craton edge touches, for as long as it touches it. A craton edge parked at
a trench retreats **about 1,100 km per Gyr**, and a craton 1,000–2,000 km across is gone
in 1–2 Gyr.

The run behind RFCs 0001–0006 shows it. The run is seed 456862899, 19.2 Gyr, saved locally
as `mantle-bloom-seed456862899-19191800000y.mbworld`; the file is on the author's machine
and **not available for download**. Numbers come from `bin/debug/ledger_history.py` on that
save.

- `subducted_m3` was the largest craton sink: 1,128 M km³, **42%** of all craton
  destruction.
- It was steady rather than lumpy: no single 10 Myr sample took more than 0.5% of the live
  craton.
- Per Gyr, it took 6%, 18%, 7%, 12% and 10% of the live craton in the first five Gyr. Later
  windows ran at 20–77%, as the cratons shrank and had more edge per unit volume.

Geologically this is backwards. Archean cratons have kept most of their extent for 2.5–3.5
Gyr [[Pearson et al., 2021]](#ref-pearson2021). Craton destruction is rare: the eastern
North China Craton is the textbook case, and it lost its *mantle root* under more than 100
Myr of Pacific subduction while its crust stayed at the surface
[[Zhu et al., 2012]](#ref-zhu2012).

The ledger also lumps three paths into `subducted_m3`:
- `cratons.attrite`, edge wear at trenches;
- the craton on cells `quad_tectonics._retreat` removes down a trench, which should be ~0
  because `cratons.retreat_allowed` keeps craton cells from retreating;
- the lost share at polarized collision fronts (`crust_transfer`, #320).

Until they are split, we can't say for certain how much is attrition.

## Proposed solutions

**What controls subduction erosion.** Margins are either *accretionary*, growing an
accretionary prism from offscraped trench sediment, or *erosive*, losing forearc crust
down the trench. [[von Huene & Scholl, 1991]](#ref-vonhuene1991) and
[[Clift & Vannucchi, 2004]](#ref-clift2004) found the split is set mainly by trench
sediment and convergence rate. Clift & Vannucchi find tectonic erosion favoured where
convergence exceeds 6 cm/yr *and* the trench's sedimentary cover is under 1 km. Accretion
favours slow convergence (below 7.6 cm/yr) and/or trench sediment over 1 km. Erosion happens at the forearc, in the upper plate's crust
above the slab, and erosive episodes alternate with accretion, collisions and ridge
subduction. Cratons survive because they sit behind younger, weaker mobile belts that take
the deformation [[Lenardic et al., 2003]](#ref-lenardic2003). Where cratonic lithosphere is
destroyed, the root is removed by refertilisation, heating and delamination
[[Foley, 2008]](#ref-foley2008) [[Lee et al., 2011]](#ref-lee2011)
[[Zhu et al., 2012]](#ref-zhu2012). Wholesale removal of the crust is not how it goes.

1. **Split the ledger first.** Add `attrited_m3` (from `cratons.attrite`) and
   `collision_lost_m3` (crust_transfer's lost and unplaced shares) to the craton ledger.
   `subducted_m3` keeps only what `_retreat` takes down a trench, and an assertion or
   diagnostic flags it as non-zero. Record attrition separately for continental and
   oceanic host plates, since only the continental case is capped by the #177 arc budget.
2. **Attrition only at erosive margins.** Classify each trench contact as accretionary or
   erosive from the subducting plate's trench sediment and the convergence rate, using the
   thresholds of [[Clift & Vannucchi, 2004]](#ref-clift2004). An accretionary contact wears
   no craton. With RFC 0005, it also keeps the offscraped sediment on the upper plate. An
   erosive contact keeps today's 40 km²/Myr.
3. **Erode the root, not the crust.** At an erosive contact on a cratonic margin, attrition
   first thins the mantle lithosphere (Hm) of the exposed cells toward the continental
   reference. Only once the keel is gone does the cell lose cratonic status
   (`decratonised_m3`); its crust stays continental. This matches the North China history
   [[Zhu et al., 2012]](#ref-zhu2012) and keeps continental crust that today vanishes.
4. **A plain resistance factor.** Scale attrition by
   `1 − CRATON_ATTRITION_RESISTANCE × strength`, as erosion, rifting and stretching already
   do. It is simple, but has no physical basis for any particular value, so it's listed
   only as the fallback.

**Recommendation:** 1 then 2. Solution 2 has thresholds taken straight from the literature,
uses fields the model already carries, and also benefits RFC 0005, which needs the same
classification for sediment. Re-measure after 2. If the attrition share is still above the
target, add 3. It is the more physical answer, but changes how cratons die and needs its
own calibration.

## Design

**Units.**
- Attrition rate: km²/Myr (= m²/yr) of cross-section per unit edge length. A cell of area A
  loses `rate × years / sqrt(A)` metres of craton per step.
- Trench sediment: metres of sediment on the subducting plate at the contact.
  `mobile_cover_m` is the unconsolidated cover, and continental material on an oceanic host
  (`continental_material_m`) is consolidated sediment; trench fill is their sum.
- Convergence: the node's normal closing rate, `rheology.normal_closing_rate_m_per_s`,
  converted to cm/yr for the threshold.

**Classification (2).** Per node in the trench set of `_retreat` on the *overriding*
(craton) plate:
- read the sediment at its nearest neighbour node, through
  `inputs.neighbor_node_index` / `neighbor_plate_id`, exactly as
  `torque.contact_craton_m` reads the neighbour's craton;
- the node is erosive when sediment < `EROSIVE_MAX_TRENCH_SEDIMENT_M` *and* closing rate
  > `EROSIVE_MIN_CONVERGENCE_CM_PER_YR`, the joint condition of
  [[Clift & Vannucchi, 2004]](#ref-clift2004), and accretionary otherwise;
- attrition runs only on erosive nodes.

The rule is binary, as in the source. A smooth ramp across ±25% of each threshold would
avoid flicker between steps; it is listed as a tuning choice, not a default.

**Root-first variant (3), if needed.**
- At an erosive cratonic node, Hm thins at `CRATON_KEEL_EROSION_M_PER_MYR` until it reaches
  `REFERENCE_HM_CONTINENTAL_M` (100 km). The mantle removed books to the Hm ledger
  (`hm_ledger`, account `oceanic_and_deep_subduction`).
- Below that, the cell's craton converts to ordinary continental crust, booked to the new
  craton account `decratonised_m3`. Hc and `continental_material_m` are untouched.
- The thinning rate is set so a cratonic margin loses its keel over ~100 Myr of continuous
  erosive subduction, the North China timescale [[Zhu et al., 2012]](#ref-zhu2012):
  about 1 km/Myr for a keel 100 km above reference.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `EROSIVE_MAX_TRENCH_SEDIMENT_M` | 1,000 m | [[Clift & Vannucchi, 2004]](#ref-clift2004) | Fixed unless sediment thickness in the model is off by a large factor; check its distribution at trench contacts first |
| `EROSIVE_MIN_CONVERGENCE_CM_PER_YR` | 6 cm/yr | [[Clift & Vannucchi, 2004]](#ref-clift2004) | Fixed; 7.6 cm/yr (their upper bound for accretion) is the sensitivity case |
| `CRATON_ATTRITION_KM2_PER_MYR` | 40 (unchanged) | [[Scholl & von Huene, 2007]](#ref-scholl2007) | Unchanged |
| `CRATON_KEEL_EROSION_M_PER_MYR` (variant 3 only) | 1,000 m/Myr | [[Zhu et al., 2012]](#ref-zhu2012) timescale | Sweep 300–3,000 |

How tuning works:
- Use `bin/debug/compare_cratons.py` on seeds 7, 11 and 2024 to 1 Gyr, plus a regenerated
  seed 456862899 to 3 Gyr. Run the replays one at a time.
- Report each craton sink per Gyr, the fraction of trench contact length classed
  accretionary, and genuine continental area.
- Both regimes are common on Earth. [[von Huene & Scholl, 1991]](#ref-vonhuene1991)
  counted about 19,000 km of non-accretionary margin, and
  [[Clift & Vannucchi, 2004]](#ref-clift2004) compiled both kinds around the Pacific and
  Indian oceans. If either regime covers less than 20% of the model's trench length, the
  model's trench sediment or convergence rates are off, and RFC 0005 should look at that
  before these thresholds move.

**Runtime.** Classification reads one field per trench node from a neighbour, with indices
already computed for the step: O(trench nodes), well under 1 ms per plate. Variant 3 adds an
O(trench nodes) Hm update. There's no new search.

**Save compatibility.** No new per-node field. New craton ledger accounts (`attrited_m3`,
`collision_lost_m3`, and `decratonised_m3` for variant 3) are backfilled at 0.0 by
`cratons.ensure_ledger`. A pre-change save keeps its history in `subducted_m3`.

## Detailed design

- `backend/app/cratons.py`
  - Add `attrited_m3`, `collision_lost_m3` (and `decratonised_m3`) to
    `CRATON_SINK_ACCOUNTS`.
  - `attrite` books to `attrited_m3` and takes an `erosive` mask (or the classification
    inputs).
  - For variant 3, a `wear_keel` helper next to `attrite`.
  - Update the module docstring's Destruction paragraph.
- `backend/app/quad_tectonics.py`, `_retreat`: build the erosive mask for the `trench` set
  before calling `cratons.attrite`; split the attrition booking by host plate type. The
  #177 budget logic is unchanged.
- `backend/app/torque.py`: add a `contact_trench_sediment_m(plate, neighbours, inputs)`
  next to `contact_craton_m`, reading `mobile_cover_m + continental_material_m` at the
  nearest neighbour node, counting the latter only on oceanic hosts.
- `backend/app/lithosphere_plate.py`: expose the closing rate in the `BoundaryContext` if
  `_retreat` doesn't already have it.
- `backend/app/crust_transfer.py`: book the lost and unplaced craton shares to
  `collision_lost_m3` instead of `subducted_m3`.
- `backend/app/stats.py` and `cratons.diagnostics`: report the new accounts and the
  accretionary fraction of trench length.
- `bin/debug/compare_cratons.py`: report the new accounts per Gyr.
- Tests:
  - `backend/unit_tests/test_cratons.py`: attrition skips accretionary contacts; erosive
    contacts lose 40 km²/Myr × edge; the ledger balances with the new accounts; an old
    ledger without them loads.
  - `backend/unit_tests/test_torque.py`: `contact_trench_sediment_m` reads the neighbour's
    nearest node.
- Docs: `docs/simulation-model.md`, Cratons ("Resistance" and the destruction ledger), and
  the boundary-retreat description in "Quad-surface deformation".

## Acceptance criteria

- [ ] Ledger split: on a 1 Gyr replay, `subducted_m3` is below 0.1% of the initial craton
      volume, and `attrited_m3` and `collision_lost_m3` account for the rest of what was
      booked there before.
- [ ] Craton loss to attrition is at most **3% of live craton per Gyr** over the first
      3 Gyr on seed 456862899 (6–18% today), and on average across seeds 7, 11 and 2024.
- [ ] Each regime covers at least 20% of trench contact length, averaged over a run.
- [ ] Genuine continental area at 1 Gyr is no lower than on `origin/main` for each seed.
- [ ] Ledgers balance (craton within 1e-9 of live volume; continental material closes).
- [ ] Step time on a density-4 world changes by less than 2% (`docs/profiling.md` method).
- [ ] Unit tests above pass; affected tests, ruff and `bin/quality_check.sh` are green.
- [ ] Saves from before the change load.

## References

- <a id="ref-clift2004"></a>Clift, P., & Vannucchi, P. (2004). Controls on tectonic
  accretion versus erosion in subduction zones: Implications for the origin and recycling
  of the continental crust. *Reviews of Geophysics*, 42, RG2001.
  [doi:10.1029/2003RG000127](https://doi.org/10.1029/2003RG000127)
- <a id="ref-foley2008"></a>Foley, S. F. (2008). Rejuvenation and erosion of the cratonic
  lithosphere. *Nature Geoscience*, 1, 503–510.
  [doi:10.1038/ngeo261](https://doi.org/10.1038/ngeo261)
- <a id="ref-lee2011"></a>Lee, C.-T. A., Luffi, P., & Chin, E. J. (2011). Building and
  destroying continental mantle. *Annual Review of Earth and Planetary Sciences*, 39,
  59–90. [doi:10.1146/annurev-earth-040610-133505](https://doi.org/10.1146/annurev-earth-040610-133505)
- <a id="ref-lenardic2003"></a>Lenardic, A., Moresi, L.-N., & Mühlhaus, H. (2003).
  Longevity and stability of cratonic lithosphere: Insights from numerical simulations of
  coupled mantle convection and continental tectonics. *Journal of Geophysical Research:
  Solid Earth*, 108(B6), 2303.
  [doi:10.1029/2002JB001859](https://doi.org/10.1029/2002JB001859)
- <a id="ref-pearson2021"></a>Pearson, D. G., Scott, J. M., Liu, J., Schaeffer, A., Wang,
  L. H., van Hunen, J., Szilas, K., Chacko, T., & Kelemen, P. B. (2021). Deep continental
  roots and cratons. *Nature*, 596, 199–210.
  [doi:10.1038/s41586-021-03600-5](https://doi.org/10.1038/s41586-021-03600-5)
- <a id="ref-scholl2007"></a>Scholl, D. W., & von Huene, R. (2007). Crustal recycling at
  modern subduction zones applied to the past—Issues of growth and preservation of
  continental basement crust, mantle geochemistry, and supercontinent reconstruction.
  *Geological Society of America Memoir* 200, 9–32.
  [doi:10.1130/2007.1200(02)](https://doi.org/10.1130/2007.1200(02))
- <a id="ref-vonhuene1991"></a>von Huene, R., & Scholl, D. W. (1991). Observations at
  convergent margins concerning sediment subduction, subduction erosion, and the growth of
  continental crust. *Reviews of Geophysics*, 29, 279–316.
  [doi:10.1029/91RG00969](https://doi.org/10.1029/91RG00969)
- <a id="ref-zhu2012"></a>Zhu, R., Xu, Y., Zhu, G., Zhang, H., Xia, Q., & Zheng, T.
  (2012). Destruction of the North China Craton. *Science China Earth Sciences*, 55,
  1565–1587. [doi:10.1007/s11430-012-4516-y](https://doi.org/10.1007/s11430-012-4516-y)
