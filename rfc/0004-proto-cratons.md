# RFC 0004: Proto-cratons stabilised by collision

- **Status:** Draft
- **Date:** 2026-10-13
- **Related:** #274 (craton formation); RFCs [0001](0001-merge-cap-craton-loss.md),
  [0003](0003-rifting-around-cratons.md), [0006](0006-intra-oceanic-arcs.md)

## Problem

Today a craton forms in one way only (`cratons.update`). A column must:
- be genuine continental crust, with Hc in [`CRATON_HOST_MIN_HC_M`,
  `CRATON_HOST_MAX_HC_M`] = [28 km, 50 km];
- sit at least `CRATON_FORMATION_MARGIN_KM` = 500 km inside its plate's edge and any
  non-genuine cell;
- stay that way for `CRATON_FORMATION_MYR` = 200 Myr.

Its whole Hc then becomes cratonic. In practice only the deep interiors of large
continents ever qualify, and those are mostly craton already.

In the run behind RFCs 0001–0006, **craton formation stopped for good at about 500 Myr**.
The run is seed 456862899, 19.2 Gyr, saved locally as
`mantle-bloom-seed456862899-19191800000y.mbworld`; the file is on the author's machine and
**not available for download**. Numbers come from `bin/debug/ledger_history.py` on that
save.

- `formed_m3` reached 357 M km³ by 500 Myr (201 M km³ of it in one sample at 202 My, the
  first time the 200 Myr clock completed). It then stayed at exactly that value for the
  remaining 18.7 Gyr.
- By 500 Myr cratons made up 73% of genuine continental area, and 80–88% for most of the
  rest of the run. The non-cratonic crust that remained was all near margins.

With no source, every craton sink in RFCs 0001–0003 is a one-way drain.

The user's proposal is a *proto-craton* stage: crust that can mature into craton, with
continental collisions as a main route.

## Proposed solutions

**What happens in nature.**
- Cratons are continental cores underlain by thick (about 200–250 km), cold and chemically
  depleted mantle keels [[Pearson et al., 2021]](#ref-pearson2021). Depletion makes the
  keel buoyant enough to offset its cold, dense temperature, which is why it lasts
  [[Jordan, 1978]](#ref-jordan1978).
- Those keels formed mostly in the Archean, by high-degree melting of a hotter mantle
  [[Herzberg & Rudnick, 2012]](#ref-herzberg2012). Lithospheric mantle formed later is
  progressively less depleted [[Griffin et al., 2009]](#ref-griffin2009).
- The depleted lithosphere was then *thickened and stabilised by lateral compression*,
  i.e. collisions [[Capitanio et al., 2020]](#ref-capitanio2020)
  [[Lee et al., 2011]](#ref-lee2011).
- Continents also grow outward: accretionary and collisional belts weld onto older cores
  and become part of the stable shield once they stop deforming
  [[Cawood et al., 2009]](#ref-cawood2009). Laurentia's Proterozoic belts around its
  Archean nuclei are the classic case.

So the proposal is sound, with two caveats:
- A freshly thickened orogen is hot and weak. It collapses, delaminates or rifts before it
  cools, so stabilisation follows the *end* of collision, not the collision itself.
- After the Archean, new cratons should be rarer and less depleted
  [[Griffin et al., 2009]](#ref-griffin2009).

1. **Welded growth.** A quiet genuine-continental column within `CRATON_WELD_REACH_KM` of
   existing craton keeps a formation clock running even when it is closer than 500 km to
   the plate edge. It needs only `CRATON_WELD_MARGIN_KM` of interior. Cratons then grow at
   their edges by absorbing the belts welded onto them.
2. **Collisional keel.** A quiet genuine-continental column whose mantle lithosphere was
   thickened past `PROTO_CRATON_MIN_HM_M` also qualifies, wherever it is. In the model, Hm
   thickens under convergent shortening (`rheology.apply_convergent_deformation`), and
   continental Hm does not relax back the way oceanic Hm does. This is the collision path
   of [[Capitanio et al., 2020]](#ref-capitanio2020).
3. **Name the stage.** A *proto-craton* is a non-cratonic column whose formation clock is
   running on any of the three paths (interior, welded, keel). It gets no extra mechanical
   resistance: young belts are weak (see RFC 0003). It is reported in stats and shown in
   the existing `craton` debug view, so formation can be watched.
4. **Secular decline (optional).** Scale the cratonic share given at formation by
   `CRATON_FORMATION_DEPLETION`, to mimic less-depleted post-Archean keels. The model
   seeds its starting cratons 2.5 Gyr old, so the simulated era is post-Archean.

**Recommendation:** 1, 2 and 3, with 4 left at 1.0 (off) until a run shows cratons
ratcheting up to most of the continental area. Path 1 alone restarts formation on the
margins of existing cratons. Path 2 lets collisions far from any craton seed new ones,
which is what the user asked for.

## Design

**Quantities and units.**
- Hc, Hm and craton thicknesses in metres.
- Distances in km, converted to graph hops with `cratons._hops` (62.5 km per hop at
  density 4, 125 km at density 1).
- The formation clock `stable_continental_myr` in Myr; it already exists and is a CLOCK
  field.

**Formation rule.** `quiet` currently means `host & (Hc ≤ CRATON_HOST_MAX_HC_M) &
(interior ≥ hops(500 km))`. It becomes `host & (Hc ≤ CRATON_HOST_MAX_HC_M) & (A | B | C)`:

- **A (interior):** `interior ≥ hops(CRATON_FORMATION_MARGIN_KM)`, as today.
- **B (welded):** `interior ≥ hops(CRATON_WELD_MARGIN_KM)` and
  `hops_to_craton ≤ hops(CRATON_WELD_REACH_KM)`. Here `hops_to_craton` is
  `quad_tectonics.hop_distance(plate, craton > 0, max_hops)`.
- **C (keel):** `interior ≥ hops(CRATON_WELD_MARGIN_KM)` and `Hm ≥ PROTO_CRATON_MIN_HM_M`.

The interior test for all three still excludes plate edges and non-genuine cells, so an
active margin never qualifies. The 200 Myr clock and the Hc window are unchanged: a belt
must relax back to ≤ 50 km of crust and stay quiet for 200 Myr, so stabilisation follows
the end of the orogeny. On formation the craton is
`CRATON_FORMATION_DEPLETION × Hc`, 1.0 by default.

**Runtime.** Path B adds one bounded BFS per plate per step, out to
`hops(CRATON_WELD_REACH_KM)` = 8 hops at density 4. That is the same kind and size of work
as the interior BFS `cratons.update` already does. Path C is one vector comparison.
Expected cost: under 1% of step time. Acceptance requires measuring it.

**Save compatibility.** No new per-node field: the clock, Hm and craton fields exist.
Proto-craton status is derived each step. The new stats keys are additive.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `CRATON_WELD_REACH_KM` | 500 km | One belt width: orogens welded onto shields are a few hundred km wide [[Cawood et al., 2009]](#ref-cawood2009) | Sweep 250–1,000 km |
| `CRATON_WELD_MARGIN_KM` | 250 km | Half the interior margin: keeps active margins and their arcs out | Sweep 125–500 km |
| `PROTO_CRATON_MIN_HM_M` | 150,000 m (1.5 × the 100 km reference) | Toward the 200–250 km depth of cratonic keels [[Pearson et al., 2021]](#ref-pearson2021) | First measure the Hm distribution in quiet former collision belts on replays, then set to the value that ~10–30% of them reach |
| `CRATON_FORMATION_DEPLETION` | 1.0 | Off; see solution 4 | Lower only if the ratchet criterion fails |
| `CRATON_FORMATION_MYR` | 200 Myr (unchanged) | #274 | Unchanged |

Tune with `bin/debug/compare_cratons.py` on seeds 7, 11 and 2024 to 2 Gyr (formation needs
time to show), plus a regenerated seed 456862899 to 3 Gyr; run the replays one at a time.
Report `formed_m3` per Gyr by path (A, B, C), craton share of genuine continental area,
and proto-craton area.

## Detailed design

- `backend/app/cratons.py`
  - New constants as above, with comments citing this RFC.
  - `update`: compute paths A, B and C, combine them into `quiet`, and apply
    `CRATON_FORMATION_DEPLETION` at formation.
  - Book `formed_m3` as now, plus per-path diagnostic counters kept on the ledger as
    `formed_interior_m3`, `formed_welded_m3` and `formed_keel_m3`. These are
    sub-accounts of `formed_m3`, outside the balance.
  - New `proto_craton_mask(plate)`, shared by `update`, `diagnostics` and rendering.
  - `diagnostics`: add `proto_craton_area_m2` and the per-path formed volumes.
  - Update the module docstring's Formation paragraph.
- `backend/app/stats.py`: `proto_craton_area_km2` in the stats row.
- `backend/app/render_image.py`: the `craton` debug view (`_render_craton_view`) already
  colours quiet interiors by their clock's progress. It takes its mask from
  `proto_craton_mask`, so welded and keel proto-cratons show too.
  `frontend/src/legendData.ts`'s `craton` legend relabels "Quiet interior near
  cratonising" as "Proto-craton near cratonising".
- `bin/debug/compare_cratons.py` and `bin/debug/ledger_history.py`: report per-path
  formation.
- Tests (`backend/unit_tests/test_cratons.py`):
  - a quiet belt 300 km from the edge, next to a craton, forms after 200 Myr (path B), and
    the same belt with no craton nearby does not;
  - a quiet column with Hm ≥ the threshold forms (path C);
  - an active margin never forms;
  - ledger balance is unchanged by the sub-accounts.
- Docs: `docs/simulation-model.md`, the Cratons section ("Formation"), and the render-view
  description if the craton view changes.

## Acceptance criteria

- [ ] On seeds 7, 11 and 2024 to 2 Gyr, and seed 456862899 to 3 Gyr, `formed_m3` increases
      in every 500 Myr window in which genuine continental area is at least 25% of its
      initial value. Today it stops after 500 Myr.
- [ ] Paths B and C each contribute in at least two of the four runs.
- [ ] No ratchet: craton share of genuine continental area stays between 30% and 65% at
      every 500 Myr checkpoint (it starts at 35%).
- [ ] Step time changes by less than 2% at density 4.
- [ ] Unit tests above pass; ledgers balance; affected tests, ruff and
      `bin/quality_check.sh` are green.

## References

- <a id="ref-capitanio2020"></a>Capitanio, F. A., Nebel, O., & Cawood, P. A. (2020).
  Thermochemical lithosphere differentiation and the origin of cratonic mantle. *Nature*,
  588, 89–94. [doi:10.1038/s41586-020-2976-3](https://doi.org/10.1038/s41586-020-2976-3)
- <a id="ref-cawood2009"></a>Cawood, P. A., Kröner, A., Collins, W. J., Kusky, T. M.,
  Mooney, W. D., & Windley, B. F. (2009). Accretionary orogens through Earth history.
  *Geological Society, London, Special Publications*, 318, 1–36.
  [doi:10.1144/SP318.1](https://doi.org/10.1144/SP318.1)
- <a id="ref-griffin2009"></a>Griffin, W. L., O'Reilly, S. Y., Afonso, J. C., & Begg,
  G. C. (2009). The composition and evolution of lithospheric mantle: A re-evaluation and
  its tectonic implications. *Journal of Petrology*, 50, 1185–1204.
  [doi:10.1093/petrology/egn033](https://doi.org/10.1093/petrology/egn033)
- <a id="ref-herzberg2012"></a>Herzberg, C., & Rudnick, R. (2012). Formation of cratonic
  lithosphere: An integrated thermal and petrological model. *Lithos*, 149, 4–15.
  [doi:10.1016/j.lithos.2012.01.010](https://doi.org/10.1016/j.lithos.2012.01.010)
- <a id="ref-jordan1978"></a>Jordan, T. H. (1978). Composition and development of the
  continental tectosphere. *Nature*, 274, 544–548.
  [doi:10.1038/274544a0](https://doi.org/10.1038/274544a0)
- <a id="ref-lee2011"></a>Lee, C.-T. A., Luffi, P., & Chin, E. J. (2011). Building and
  destroying continental mantle. *Annual Review of Earth and Planetary Sciences*, 39,
  59–90. [doi:10.1146/annurev-earth-040610-133505](https://doi.org/10.1146/annurev-earth-040610-133505)
- <a id="ref-pearson2021"></a>Pearson, D. G., Scott, J. M., Liu, J., Schaeffer, A., Wang,
  L. H., van Hunen, J., Szilas, K., Chacko, T., & Kelemen, P. B. (2021). Deep continental
  roots and cratons. *Nature*, 596, 199–210.
  [doi:10.1038/s41586-021-03600-5](https://doi.org/10.1038/s41586-021-03600-5)
