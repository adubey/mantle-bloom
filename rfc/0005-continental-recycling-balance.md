# RFC 0005: Balance continental recycling at trenches

- **Status:** Draft
- **Date:** 2026-10-13 (re-baselined 2026-10-14 on the year-based split pace, #403)
- **Related:** #120 (land fraction decline), #177 (margin retreat budget), #276 and #320
  (collision sinks), #297 (mobile cover), #358 (split cooldown in years), #403 (re-baseline),
  #412 (gate scope); RFCs [0002](0002-craton-trench-attrition.md),
  [0006](0006-intra-oceanic-arcs.md)

## Problem

Continental crust in the model leaves the surface much faster than it is made, and the
continents drain away. The continental-material ledger (`continental_ledger.py`) books
every source and sink.

In the run behind RFCs 0001–0007, the ledger shows the imbalance from the first step. The
run is seed 456862899, 19.2 Gyr, saved locally as
`mantle-bloom-seed456862899-19191800000y.mbworld`; the file is on the author's machine and
**not available for download**. Numbers come from `bin/debug/ledger_history.py` on that
save. The save predates #358: it ran post-#340 code (exact commit unknown) with the split
cooldown counted in steps, 200 Myr for most of the run.

| Window (Myr) | Juvenile (km³/yr) | Remelted returns (km³/yr) | All sinks (km³/yr) | of which `deeply_subducted` | Mean inventory (M km³) | Share on oceanic hosts at window end |
|---|---|---|---|---|---|---|
| 0–100 | 7.3 | 0.7 | 11.3 | 10.1 | 7,851 | 6% |
| 100–500 | 2.6 | 1.1 | 8.5 | 7.7 | 6,460 | 18% |
| 500–1,000 | 4.3 | 1.2 | 8.4 | 7.4 | 4,638 | 30% |
| 1,000–2,000 | 2.8 | 0.7 | 5.2 | 4.9 | 3,228 | 25% |
| 2,000–4,000 | 2.1 | 1.1 | 3.5 | 3.1 | 1,378 | 35% |
| 12,000–19,200 | 0.4 | 0.3 | 0.8 | 0.7 | 128 | 72% |

The world started with 7.9 × 10⁹ km³ of continental material, almost exactly Earth's
8.0 × 10⁹ km³ [[Clift & Vannucchi, 2004]](#ref-clift2004), so its fluxes are on the same
scale as Earth's.

**On the current model** (origin/main 2562f63, 5 Myr steps; `analysis/issue403/`), the drain
is two to three times faster. Fluxes in km³/yr:

| Run | Start inventory (10⁹ km³) | Window (Myr) | Juvenile | All sinks | `deeply_subducted` | `collision_subducted` | Inventory change |
|---|---|---|---|---|---|---|---|
| Seed 7 | 12.0 | 0–100 / 100–500 / 500–1,000 | 10.7 / 8.3 / 5.0 | 22.6 / 27.9 / 9.8 | 11.7 / 12.3 / 6.4 | 5.8 / 13.1 / 2.3 | −86% in 1 Gyr |
| Seed 11 | 6.9 | same | 10.4 / 7.3 / 4.2 | 16.8 / 15.0 / 9.0 | 15.6 / 12.6 / 7.4 | 0.0 / 1.4 / 1.0 | −73% |
| Seed 2024 | 6.5 | same | 6.8 / 4.9 / 3.1 | 9.8 / 10.2 / 9.7 | 8.3 / 6.4 / 6.4 | 0.0 / 2.1 / 2.5 | −75% |
| Seed 456862899 | 3.9 | same | 7.0 / 3.4 / 1.7 | 12.3 / 10.1 / 3.0 | 10.3 / 7.6 / 2.9 | 0.9 / 1.8 / 0.0 | −85% |

- The regenerated seed 456862899 isn't the save's world: generation has changed, and it
  starts with half the save's continental material. The tuning seeds start at 6.5–12.0 ×
  10⁹ km³.
- Sinks are 9.8–28 km³/yr while continents are large, 2–6 times Earth's recycling estimate
  below (the save: about twice).
- `deeply_subducted_m3` is still most of it (44–95% per window), but
  `collision_subducted_m3`, where merges and sutures book their unplaced remainder, is 47%
  of seed 7's sinks over 100–500 Myr and 21–26% of seed 2024's after 100 Myr. That account
  is outside this RFC's channel split; RFC 0001 is about the same merges.

The run is a stress test of the model's long-term behaviour, not a calibration interval.
The first Gyr's fluxes are compared with Earth's below only because that is when the
model's continents are Earth-sized. Later windows describe a degenerate world.

**Earth's budget.** For the Cenozoic, [[Clift et al., 2009]](#ref-clift2009) estimate these
sinks:

| Sink | km³/yr |
|---|---|
| Sediment subduction | 1.65 |
| Subducted eroded forearc | 1.3 |
| Lower-crust delamination | ~1.1 |
| Collision losses | 0.4 |
| Material fixed in oceanic crust | ~0.4 |
| **Total** | **≈ 4.9** |

Arcs (about 2.0–3.8 km³/yr) and accreted oceanic plateaus (about 1.1) replace it, and the
crust has been in rough equilibrium through the Phanerozoic
[[Clift et al., 2009]](#ref-clift2009) [[Stern & Scholl, 2010]](#ref-stern2010).

**These are sanity checks, not targets.** Clift et al.'s figures are global estimates for
specific processes, mostly over the Cenozoic, with large uncertainties. Three kinds of
number must also be kept apart:
- **gross production**: what a process makes (arc magma, eroded sediment);
- **net preservation**: what survives as continental crust;
- **model-ledger fluxes**: what an account books in a step.

The comparison here is ledger sinks against Earth's recycling estimates: both are
material leaving the continental inventory, so they are comparable, with that caveat. In the
pre-#358 save, the model's new crust, 2.6–4.3 km³/yr after the first 100 Myr, is within
Earth's range for arc production. Its sinks, 8–11 km³/yr, are about twice Earth's recycling
estimate, and nearly all of the excess is `deeply_subducted_m3`. On origin/main 2562f63,
juvenile additions are 3–11 km³/yr and the sinks 2–6 times Earth's, with
`collision_subducted_m3` a second large channel on two of the four runs.

That account mixes several channels:
- continental material, mostly sediment, on oceanic plates consumed at trenches
  (`quad_tectonics._retreat`);
- continental margins consumed under an oceanic neighbour, which on continental plates is
  capped 1:1 by arc production (#177);
- craton attrition (`cratons.attrite`; RFC 0002).

Because the #177 cap holds margin consumption to arc output, **the prime suspect is
sediment**: continental material eroded onto ocean floor and carried into trenches. The
share of continental material sitting on oceanic hosts rose from 0 to 30% in the first Gyr,
which supports that. But the ledger can't confirm it until the channels are split.

## Proposed solutions

**How sediment moves on Earth.**
- Rivers delivered about 14 Gt/yr of sediment to the coast before human influence
  [[Syvitski et al., 2005]](#ref-syvitski2005). At 2.65 t/m³ that is about 5.3 km³/yr of
  rock.
- About 1.65 km³/yr of sediment is subducted [[Clift et al., 2009]](#ref-clift2009). The
  two figures have different scopes (modern rivers versus a Cenozoic global budget), so
  their ratio is not a measured share of sediment kept on margins, and this RFC doesn't
  use it as one.
- At the trench, accretion is inefficient. A median of only 17% (range 7–37%) of the
  sediment arriving is scraped off into the accretionary prism; the rest is subducted
  [[Clift & Vannucchi, 2004]](#ref-clift2004). Erosive and accretionary margins are told
  apart by trench sediment and convergence rate
  [[von Huene & Scholl, 1991]](#ref-vonhuene1991).
- Some subducted sediment returns: it melts into arc magmas, and buoyant subducted crust
  can relaminate to the base of the overriding crust
  [[Hacker et al., 2011]](#ref-hacker2011).

1. **Measure the channels (required first).** Split `deeply_subducted_m3` into:
   - `sediment_subducted_m3`: material on oceanic hosts, consumed at a trench;
   - `margin_subducted_m3`: continental hosts consumed under an oceanic neighbour, the
     #177 channel;
   - `craton_attrited_m3`: from `cratons.attrite`.

   Add the flux diagnostics defined in Design (G, E, D, T, S, A, P). These tell us whether
   the excess is *erosion* (too much crust eroded), *export* (too much sediment reaching
   the deep ocean) or *consumption* (too much of what reaches the trench subducted). No
   fix is chosen until this measurement is in.
2. **Correct the dominant flux.** Which fix depends on the measurement:
   - Erosion high (G well above Earth's sediment delivery to the coast): erosion
     calibration, which is out of scope here, so a separate RFC.
   - Export high (most of E ends on oceanic hosts): a sediment-routing change. One
     candidate is an accommodation-limited priority for continental-host cells (Design).
     It needs independent justification before it's chosen, from the passive-margin
     sedimentation literature and from the model's own shelf geometry; the Earth budget
     above doesn't establish where sediment should settle.
   - Consumption high (S/T near 1 at margins the classification calls accretionary): 3.
3. **Accretionary prisms.** At contacts RFC 0002 classes as accretionary, transfer
   `OFFSCRAPE_FRACTION` of the arriving sediment onto the overriding plate's nearest cells
   as continental material instead of subducting it.
4. **Sediment return through arcs.** Book a fraction of subducted sediment as
   `remelted_relaminated_returns_m3` on the overriding arc band, after
   [[Hacker et al., 2011]](#ref-hacker2011). It has a large effect and is hard to
   constrain, so it is kept as a last resort.

**After the re-baseline (#403).** `collision_subducted_m3` is up to 47% of the sinks on
one tuning seed, so the gate's measurement should report it beside the three channels of
1, and the gate should say whether it or `deeply_subducted_m3` is the dominant excess
before a branch of 2 is chosen. If it is, the fix is RFC 0001's merge placement, not this
RFC's (#412).

**Recommendation:** 1 as a diagnosis gate, then the branch of 2 that the measurement
points to. Do 3 alongside RFC 0002, which needs the same contact classification; it only
moves a median 17% of trench sediment, but it is cheap and grounded. Do 4 only if the
budget still fails once the rest is in.

## Design

**Units.** Ledger accounts are m³ and reports are km³/yr. Sediment thicknesses are metres
of rock-equivalent crust, since every sediment field is a share of Hc; Earth figures given
in mass are converted at 2.65 t/m³.

**Channel split (1).** In `quad_tectonics._retreat`, the `subducted` set is booked in one
call today. Book three ways instead:
- cells with `~continental` (oceanic hosts) → `sediment_subducted_m3`;
- continental cells on a continental plate → `margin_subducted_m3`;
- attrition as RFC 0002 books it.

`deeply_subducted_m3` stays as their sum for old saves and dashboards: the new keys are
sub-accounts, and the balance uses the parent only.

**Flux definitions (1).** All in km³/yr of continental material (the
`continental_material_m` tracer), so each is like-for-like with the others:

| Symbol | Flux | Where measured | Earth sanity value (scope) |
|---|---|---|---|
| G | Gross denudation: continental material eroded from continental hosts | `erosion.py`, removal step | ~5.3 km³/yr of sediment reaching the coast, pre-human [[Syvitski et al., 2005]](#ref-syvitski2005) (modern rivers; excludes dissolved load and glacial input) |
| E | Erosion product deposited in the sea, on any host | `erosion.py` marine settling paths | — |
| D | The part of E deposited on oceanic hosts | the same paths, split by host type | — |
| T | Trench delivery: continental material on oceanic cells consumed at trenches | `quad_tectonics._retreat` | — |
| S | Sediment subducted: T minus offscraped | `sediment_subducted_m3` | 1.65 km³/yr [[Clift et al., 2009]](#ref-clift2009) (Cenozoic, global) |
| A | Offscraped and accreted | solution 3 | median 17% of T [[Clift & Vannucchi, 2004]](#ref-clift2004) (modern margins) |
| P | Net preservation: change in the live inventory | `continental_ledger.inventories` | ~0 over the Phanerozoic [[Clift et al., 2009]](#ref-clift2009) |

E, D and the stored marine inventory are recorded per step in a world-level diagnostic
(not pickled). D/E and T/D are the model's own ratios: there is no established Earth value
for them, and none is assumed.

**Accommodation-limited routing (a candidate for 2, if export dominates).** Sediment
reaches the sea floor by several paths in `erosion.py`:
- `_spread_marine_sediment`, for submarine and sea-cliff erosion: inverse-distance over
  the nearest lower ocean nodes within `MARINE_SPREAD_RANGE_RAD`;
- `_spread_beach_sediment`;
- settling of the river-routed pool where it meets the sea.

None of them looks at the host's crust type. The flux diagnostic in (1) says which path
carries most continental material onto oceanic hosts, and the change goes there first.
That path's receiver weights give continental-host cells first claim, up to an
accommodation limit: the room between the cell's depth and `SHELF_ACCOMMODATION_DEPTH_M`,
a physical depth for the shelf, slope and upper rise. Only the remainder goes to oceanic
hosts. Each path conserves mass through `np.add.at`, and that is kept.

Before this candidate is chosen, the RFC must be revised with:
- an independent justification for the depth limit and the host priority;
- a comparison against a physically simpler alternative, transport-limited spreading by
  slope, so sediment runs out until the gradient is too low to carry it, regardless of
  host.

**Offscraping (3).** The model's sediment fields overlap: `mobile_cover_m` is a share of
Hc, and `mobile_cover_continental_m` is a share of both it and `continental_material_m`.
So the transfer is written over disjoint parts, and each physical volume moves once.

For a consumed trench cell on an oceanic plate (oceanic-coded) in the accretionary class,
with f = `OFFSCRAPE_FRACTION`, the parts are:
- L = `mobile_cover_m`: all loose cover, of which Lc = `mobile_cover_continental_m` is
  continental;
- C = `continental_material_m − Lc`: consolidated continental sediment.

This is the same disjoint column RFC 0002 uses for classification
(`torque.contact_trench_sediment_m`).

**Decision: offscraped sediment becomes consolidated prism crust on the upper plate.** An
accretionary prism is deformed, dewatered sediment, not loose cover, and this matches what
`mobile_cover.accreted_m3` already means: cover metamorphosed into accreted crust, a sink
of the *cover* inventory only. The alternative, keeping it as loose cover on the receiver,
would be a conservative cover transfer with no sink; it's rejected because prisms aren't
loose. Concretely, before the cell is removed:

| Quantity | Donor (trench cell) | Receivers (upper plate) |
|---|---|---|
| Hc volume | − f·(L + C) | + f·(L + C), spread by the placement |
| `continental_material_m` volume | − f·(Lc + C) | + f·(Lc + C), in the same proportions |
| `mobile_cover_m`, `mobile_cover_continental_m` | − f·L, − f·Lc | unchanged |
| Mobile-cover ledger | f·L booked once to `accreted_m3` | — |
| Continental ledger | nothing: the material stays live | — |

The cell's remainder is then consumed as now: its continental material books to
`sediment_subducted_m3`, and its cover to the mobile-cover `subducted_m3`.

The receivers come from the upper-plate placement in `crust_transfer.py` (#320), exposed
as a function rather than duplicated. Each receiver gets Hc and continental material in
the same proportion, so the material never exceeds Hc anywhere. On an oceanic upper plate
(an ocean–ocean front, RFC 0006), the receivers are oceanic-coded cells, and the material
stays continental material on oceanic hosts until it docks or subducts.

**Runtime.**
- The split and the flux diagnostics are O(cells) sums over masks already built,
  negligible.
- Accommodation-limited routing, if chosen, reweights receivers the spreading paths
  already query: O(sources × neighbours), with no new tree query.
- Offscraping is O(trench cells).
- Budget: under 3% of step time in total.

**Save compatibility.** New ledger keys are backfilled at 0.0 by
`continental_ledger.ensure_initialized`. No new per-node field: the flux diagnostic isn't
pickled. Old saves keep their history in `deeply_subducted_m3`.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `SHELF_ACCOMMODATION_DEPTH_M` (candidate only) | 2,000 m | Calibration assumption: the foot of the continental slope and upper rise; needs independent justification before use | Sweep 1,000–4,000 m against the acceptance targets on S and P, not against a fixed export share |
| `OFFSCRAPE_FRACTION` | 0.17 | Global median [[Clift & Vannucchi, 2004]](#ref-clift2004) | Fixed; 0.07 and 0.37 are the sensitivity bounds from the same source |
| `ocean_deposition_multiplier` (existing live knob) | 1.0 | — | Unchanged; not a substitute for 2, since below 1.0 it books the sediment it declines to a sink |

Tune on seeds 7, 11 and 2024 to 1 Gyr, plus a regenerated seed 456862899 to 3 Gyr; run the
replays one at a time. Report G, E, D, T, S, A and P, juvenile additions and total sinks,
and genuine continental area and land fraction. The `replay-ab` skill does the
before/after against a detached `origin/main` worktree.

## Detailed design (provisional touchpoints)

- `backend/app/continental_ledger.py`: add `sediment_subducted_m3`, `margin_subducted_m3`
  and `craton_attrited_m3` as documented sub-accounts of `deeply_subducted_m3`, outside
  `balance_error_m3`. `inventories` reports them.
- `backend/app/quad_tectonics.py`, `_retreat`: book the split. For (3), call the offscrape
  transfer before `plate.remove_cells`.
- `backend/app/cratons.py`, `attrite`: book `craton_attrited_m3` alongside
  `deeply_subducted_m3`.
- `backend/app/hydroclimate/erosion.py`: the G, E and D diagnostics across the removal
  step and all marine settling paths. Only if the export branch of 2 is chosen:
  host-aware routing in whichever of `_spread_marine_sediment`, `_spread_beach_sediment`
  and river-mouth settling dominates.
- `backend/app/crust_transfer.py`: expose the upper-plate placement as a function callable
  for offscraped trench sediment.
- `backend/app/hydroclimate/mobile_cover.py`: book offscraped cover once as
  `accreted_m3`, the existing sink.
- `backend/app/stats.py`: G, E, D, T, S, A and P in the stats row.
- `bin/debug/ledger_history.py`: print the new channels when present.
- Tests:
  - `test_continental_ledger.py`: the sub-accounts sum to the parent, and an old ledger
    loads.
  - `test_quad_tectonics.py`:
    - the three-way split on a synthetic trench;
    - a cross-plate offscrape case: a trench cell with known L, Lc and C, consumed under a
      continental upper plate. Check that:
      - the donor loses exactly f·(L + C) of Hc and f·(Lc + C) of material;
      - the receivers gain exactly that, cell by cell, with material ≤ Hc everywhere;
      - the mobile-cover ledger books f·L once to `accreted_m3` and the remainder to
        `subducted_m3`;
      - the continental ledger books only the remainder, to `sediment_subducted_m3`;
      - both ledgers close.
  - `test_erosion.py`: the E/D split sums to E and matches the host masks. If routing
    changes, the total is conserved.
- Docs: `docs/simulation-model.md` ("Erosion", "Quad-surface deformation"),
  `docs/api-reference.md` if stats keys are listed there.

## Assumptions and uncertainty

- **Sediment is only the prime suspect.** The channel split may show margin consumption
  or attrition dominating instead.
- **Earth values are scoped.** Every Earth figure here is a global, process-specific,
  mostly Cenozoic estimate with large uncertainty. None says where sediment should settle
  in this model. They bound magnitudes; they are not fitting targets.
- **The tracer may mislead.** `continental_material_m` on oceanic hosts mixes sediment
  with continental fragments. The channel split separates them only by host type.
- **Offscraping efficiency doesn't travel well.** The 17% median comes from modern margins
  with modern sediment lithologies. It may not hold for the model's sediment, which has
  no lithology.

## Acceptance criteria

- [ ] **Gate:** channel split and flux diagnostics are in place, and a written measurement
      (in `analysis/`) of G, E, D, T, S, A and P on the tuning seeds names the dominant
      excess, before any fix under 2–4 is chosen.
- [ ] With the chosen fixes, over the first 2 Gyr on seed 456862899 and the 1 Gyr runs of
      seeds 7, 11 and 2024 (sanity ranges from scoped Earth estimates, not fits):
  - [ ] total continental sinks are 3–7 km³/yr while the inventory is above 5 × 10⁹ km³
        (Earth's recycling estimate ≈ 4.9). On origin/main 2562f63 they are 9.8–28
        km³/yr on the tuning seeds (8–11 in the pre-#358 save). The regenerated seed
        456862899 starts at 3.9 × 10⁹ km³, so this criterion never applies to it, and it
        applies to seeds 11 and 2024 only for their first few hundred Myr (#412);
  - [ ] S, `sediment_subducted_m3`, is at most 2.5 km³/yr (Earth 1.65);
  - [ ] P, the continental inventory, changes by less than 15% per Gyr. On origin/main
        2562f63 it is −73% to −86% in the first Gyr on all four runs (−47% in the pre-#358
        save).
- [ ] Juvenile additions are unchanged within 10%: this RFC is about sinks, not sources.
- [ ] Both ledgers close; step time changes by less than 3%.
- [ ] Unit tests above pass; affected tests, ruff and `bin/quality_check.sh` are green.
- [ ] Saves from before the change load.

## Falsification criteria

- **Sediment isn't the sink.** If S is under a third of `deeply_subducted_m3` on the tuning
  seeds, the prime suspect is wrong. Re-scope to the dominant channel: margin
  consumption (the #177 cap) or attrition (RFC 0002).
- **The budget isn't the cause of the decline.** If S and the other sinks are brought
  within the sanity ranges and P still falls by more than 15% per Gyr, the loss is not at
  trenches. Look for the unledgered change (`balance_error_m3`) or a source that's too
  small (RFCs 0006 and 0007).
- **Offscraping can't matter.** If A is under 0.1 km³/yr on every seed at the 0.37
  sensitivity bound, solution 3 isn't worth its cross-plate write. Drop it.

## References

- <a id="ref-clift2004"></a>Clift, P., & Vannucchi, P. (2004). Controls on tectonic
  accretion versus erosion in subduction zones: Implications for the origin and recycling
  of the continental crust. *Reviews of Geophysics*, 42, RG2001.
  [doi:10.1029/2003RG000127](https://doi.org/10.1029/2003RG000127)
- <a id="ref-clift2009"></a>Clift, P. D., Vannucchi, P., & Morgan, J. P. (2009). Crustal
  redistribution, crust–mantle recycling and Phanerozoic evolution of the continental
  crust. *Earth-Science Reviews*, 97, 80–104.
  [doi:10.1016/j.earscirev.2009.10.003](https://doi.org/10.1016/j.earscirev.2009.10.003)
- <a id="ref-hacker2011"></a>Hacker, B. R., Kelemen, P. B., & Behn, M. D. (2011).
  Differentiation of the continental crust by relamination. *Earth and Planetary Science
  Letters*, 307, 501–516.
  [doi:10.1016/j.epsl.2011.05.024](https://doi.org/10.1016/j.epsl.2011.05.024)
- <a id="ref-stern2010"></a>Stern, R. J., & Scholl, D. W. (2010). Yin and yang of
  continental crust creation and destruction by plate tectonic processes. *International
  Geology Review*, 52, 1–31.
  [doi:10.1080/00206810903332322](https://doi.org/10.1080/00206810903332322)
- <a id="ref-syvitski2005"></a>Syvitski, J. P. M., Vörösmarty, C. J., Kettner, A. J., &
  Green, P. (2005). Impact of humans on the flux of terrestrial sediment to the global
  coastal ocean. *Science*, 308, 376–380.
  [doi:10.1126/science.1109454](https://doi.org/10.1126/science.1109454)
- <a id="ref-vonhuene1991"></a>von Huene, R., & Scholl, D. W. (1991). Observations at
  convergent margins concerning sediment subduction, subduction erosion, and the growth of
  continental crust. *Reviews of Geophysics*, 29, 279–316.
  [doi:10.1029/91RG00969](https://doi.org/10.1029/91RG00969)
