# RFC 0010: Rifts follow weak lithosphere

- **Status:** Accepted 2026-10-18
- **Date:** 2026-10-15
- **Related:** #481 (tracking), #480 (implementation), #472 / PR #479 (RFC 0009 G3), #389 (RFC 0003 tracking), #378, #385, #388 (RFC 0003 gates), PR #418 (RFC
  0003 solution 1); RFCs [0003](0003-rifting-around-cratons.md) (superseded by this RFC),
  [0009](0009-craton-mosaic.md) (craton mosaic; its gate G1 gates this RFC),
  [0007](0007-mantle-plumes-and-plateaus.md) (plume weakening)

## Problem

A splitting plate breaks along a line the mantle flow alone decides, so continents break
wherever the flow puts that line, through cratons and strong lithosphere as readily as
through weak belts.

**How a plate splits today** (`merge_split.maybe_split_plate`, origin/main `1227c36`):
1. A plate at least 2 Myr old (`SPLIT_MIN_AGE_YEARS`) whose flow one rigid rotation fits
   poorly has its node velocities clustered into two regimes (`kmeans2`).
2. `_refine_split_cut` places a great circle midway between the two clusters' centroids
   and iterates until the halves stop changing. The great circle is there for a
   computational reason: the raw cluster labels aren't spatially contiguous, and a plane
   through the centre always gives two clean halves.
3. If the halves' poles are far enough apart, the plate either splits along that circle
   (`RIFT_SUCCESS_PROBABILITY` = 0.55) or fails to rift, and `apply_failed_rift` thins a
   band `FAILED_RIFT_BAND_MULT` = 2.5 spacings either side of it.

The cut never looks at the lithosphere. On origin/main `02d35b9`, 25–36% of successful
splits cross craton on seeds 7, 11 and 2024 (PR #418); rifting takes 8.0–9.6% of the
starting craton per Gyr, and stretch donation and divergent thinning along the new
boundaries are most of that loss (#378, #388).

**RFC 0003 tried to steer the cut and failed.**
- **Tilting the great circle** (PR #418): of 158 flow cuts that crossed craton, a
  craton-free candidate passing the admissibility gates existed for 19 (12%) within 20°,
  and at most 30 (19%) anywhere in the sweep ranges. A great circle can only move as a
  whole, so missing a craton sweeps a large strip of plate to the other side.
- **Letting the cut bend** (this RFC's probe, below): a minimum cut through a strength
  graph follows the flow as well as the great circle does, but still crosses craton in
  most cases. Of 135 craton-crossing cuts on the same seeds, 13 (10%) get a craton-free
  admissible minimum cut within a 1,000 km corridor, and at most 16 (12%) at any corridor
  up to 2,500 km. Wider corridors find craton-free routes only by slicing a coastal sliver
  off the plate, which the size gate rejects ([analysis/rfc0010](../analysis/rfc0010/README.md)).

The bent cut fails for a different reason than the tilted one: craton fills each
continent's interior, leaving no seam to follow. RFC 0009 addresses that. This RFC
designs the cut that can use the seams once they exist.

What this evidence doesn't show: whether minimum cuts would avoid cratons on a mosaic
(RFC 0009's gate G1 measures it), or how much rifting loss that would remove.

## Proposed solutions

**What happens in nature.**
- Rifts localise in weak lithosphere and go around strong blocks. The East African Rift
  splits into two branches around the Tanzania craton; numerical models reproduce this from
  a single plume beneath the craton, whose keel deflects melt and strain to its edges
  [[Koptev et al., 2015]](#ref-koptev2015).
- Breakup follows old orogenic belts. The Atlantic reopened along the
  Caledonian–Appalachian suture [[Wilson, 1966]](#ref-wilson1966). Rifts run parallel to
  pre-existing orogenic fabric, explained by mechanical anisotropy frozen into the
  lithospheric mantle [[Vauchez et al., 1997]](#ref-vauchez1997)
  [[Tommasi & Vauchez, 2001]](#ref-tommasi2001), or by those belts' thinner, weaker
  lithosphere [[Audet & Bürgmann, 2011]](#ref-audet2011). Structures frozen into the
  mantle lithosphere by earlier suturing behave as latent plate boundaries
  [[Heron et al., 2016]](#ref-heron2016).
- Cratons survive in part because weaker belts around them take the strain
  [[Lenardic et al., 2003]](#ref-lenardic2003). Rifting through a craton needs it to be
  weakened first, by heating or melt metasomatism [[Foley, 2008]](#ref-foley2008).

So the flow decides *that* a plate breaks and which two parts separate; the lithosphere's
strength decides *where* the break runs.

1. **Minimum cut through a strength graph (recommended).** Keep the trigger and the
   two-regime clustering. Replace the great circle with the cheapest boundary between the
   two regimes across the plate's node graph, where cutting a link costs its length times
   the lithosphere's strength. The cut bends around strong blocks and follows weak ones,
   and still separates the regimes by construction.
2. **Nucleate and propagate.** A rift starts where tension over strength peaks and grows
   through the weakest neighbours until it reaches the plate's edge (a split) or stalls (a
   failed rift, replacing the success coin flip). Closer to how rifts form, and pairs with
   RFC 0007's plumes (three-armed rifts over plume heads). A larger change to the split
   timing and frequency; a follow-up once 1 works.
3. **Stress-field trigger.** Compute intraplate stress from the mismatch between basal
   drag and plate motion, and split where it exceeds strength. The most physical, and the
   largest change: it moves *when* plates split, not only where. Out of scope.
4. **Raise craton resistance (RFC 0003's solution 3).** Raise `CRATON_RIFT_RESISTANCE` and
   `CRATON_STRETCH_RESISTANCE`. Cheap, but treats the symptom. Kept as the fallback, sized
   to whatever loss remains after 1.
5. **Plume weakening (RFC 0007).** Heated lithosphere counts as weak in the cut cost and
   in craton resistance. Lands with RFC 0007, through this RFC's strength function.

**Recommendation:** 1, built only if RFC 0009's cratons-first G3 selects a building branch
(line clearance L ≥ 20% and band exposure B ≤ 0.75). G3 selected the fallback branch; size
4 to the remaining loss in the acceptance runs. Solution 2 is the follow-up; 5 comes with
RFC 0007.

## Design

**Gate (resolved by RFC 0009 G3).** G1 on the mosaic inside today's plates selected no
building branch. The later G3 replayed the measures on the cratons-first layout and found
that layout alone did not halve cut-band craton share on every seed. Its pooled L = 27.6%
and B = 0.732 select the fallback branch: build solution 1 and tune solution 4 against the
matched RFC 0009 baseline. G3's recorded result and method are in [#472 and PR
#479](https://github.com/adubey/mantle-bloom/issues/472). Recheck the cut-routing selection
if RFC 0012 changes split initiation before this implementation lands.

**Probe evidence (no model change).** `bin/debug/prototype_mincut_split.py` taps
`cratons.record_rift_cut`, so it sees every cut that goes ahead without changing the
trajectory, and computes the minimum cut for it. It ran seeds 7, 11 and 2024 to 1 Gyr on
a model identical to origin/main `1227c36` (5 Myr steps, default density and controls, one
at a time); [analysis/rfc0010](../analysis/rfc0010/README.md) has the method and
[mincut_probe.json](../analysis/rfc0010/mincut_probe.json) the numbers.
Of 542 cuts, 135 crossed craton. Floor 0.05; L counts admissible cuts whose line touches
no craton, and B is the craton left in the operational band (the minimum cut where
admissible, else the flow cut) over the flow cuts', pooled:

| Corridor | Band craton-free, admissible | L | B | Median residual ratio | Median length ratio |
|---|---|---|---|---|---|
| 500 km | 7 of 135 | 26 (19%) | 0.85 | 1.001–1.004 | 0.87–0.89 |
| 1,000 km | 13 | 23 (17%) | 0.73 | 1.007–1.016 | 0.80–0.83 |
| 1,500 km | 16 | 29 (21%) | 0.69 | 1.015–1.028 | 0.75–0.78 |
| 2,500 km | 13 | 21 (16%) | 0.78 | 1.065–1.093 | 0.62–0.64 |

- The minimum cut follows the flow at least as well as the great circle: within 1,000 km
  the residual rises about 1%, and the boundary (summed chord length of cut links) is
  17–20% shorter.
- It trims the craton in the operational band by about a quarter at 1,000 km, but rarely
  clears it: the size gate rejects most craton-free routes, which cut off coastal slivers.
- On today's layout it falls outside both of G1's building branches, as it should: this is
  G1's control.
- It costs at most 22 ms per cut, 33 ms with the fits.

**Strength.** Per node, dimensionless, in [0, 1]:
`S = clip(max(c, k, b), 0, 1)`, where
- `c = cratons.strength(craton_crust_m)`;
- `k = clip((Hm / REFERENCE_HM_CONTINENTAL_M − 1) / (KEEL_FULL_RATIO − 1), 0, 1)`, 0 at the
  100 km reference and 1 at a 200 km keel;
- with RFC 0009, a belt-marked column scores `RIFT_BELT_STRENGTH` (default 0) whatever its
  keel, so old sutures and rifts are the cheapest path; and with RFC 0007, heated cells
  are multiplied by `1 − PLUME_CRATON_WEAKENING`.

S ranks cuts only; it never enters a thinning factor.

**The cut.**
1. Build the plate's node graph: each node's `RIFT_GRAPH_NEIGHBOURS` = 8 nearest
   neighbours, dropping links longer than 2.5 spacings. A link of chord length ℓ between
   nodes i and j costs `ℓ / spacing × (RIFT_COST_FLOOR + (S_i + S_j) / 2)`. The floor makes
   the cut prefer a short path where strength is equal.
2. Terminals: nodes on the flow cut's side A at least `RIFT_CORRIDOR_KM` from it are tied to
   a source, nodes on side B the same distance away to a sink. The corridor bounds how far
   the rift may wander from the flow's line, as the tilt cone did for RFC 0003.
3. Solve the minimum s–t cut with `scipy.sparse.csgraph.maximum_flow` (integer
   capacities), and take the source side as half A
   [[Boykov & Kolmogorov, 2004]](#ref-boykov2004).
4. Admissibility, as for RFC 0003: both halves at least `SPLIT_MIN_NODES`; poles at least
   the relaxed `pole_separation_threshold` apart; RMS fit residual at most
   `1 + RIFT_MAX_RESIDUAL_INCREASE` times the flow cut's; each half's largest connected
   component holds at least `RIFT_MIN_CONNECTED_SHARE` of it. If any test fails, keep the
   flow cut. RFC 0003's area misfit gate is dropped: it only capped how far the cut moved,
   which the corridor now does.

**Using a curved cut.**
- `Plate.split` takes a node mask, through the existing `_plates_from_node_masks`.
- `apply_failed_rift` and `cratons.record_rift_cut` take a per-node distance to the cut
  instead of `|p · n|`: the distance to the nearest node on a cut link, from a KD-tree over
  those nodes. The band and taper are unchanged.
- With RFC 0009, the band of a successful split or failed rift is marked as a belt, so a
  rift that formed once is cheaper to reopen.

**Units and bounds.** Distances in km and chord radians; S and the floor are
dimensionless; capacities are scaled to integers, so the scale changes no cut by more
than rounding. No factor applied to a physical field changes.

**Timestep robustness.** The cut is computed once per split attempt from the current state;
attempts are already at most once per plate per step after the 2 Myr cooldown. No rate is
integrated.

**Runtime.** In the probe, a minimum cut takes at most about 12 ms on the largest plates at
default density, and the whole evaluation (graph, cut, fits) about 30 ms per cut, on
100–200 cuts per Gyr. Acceptance measures it at density 4.

**Save compatibility.** No new state; behaviour only.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| `RIFT_CORRIDOR_KM` | 1,000 km | Design choice: about a craton's width, so a cut can pass either side of one | Sweep 500–2,500 km with G1 and the acceptance runs |
| `RIFT_COST_FLOOR` | 0.05 | Design choice: the probe found 0.05 and 0.25 give nearly the same cuts | Report 0.25 |
| `RIFT_BELT_STRENGTH` | 0 | Design choice: belts are the weakest lithosphere [[Audet & Bürgmann, 2011]](#ref-audet2011) | Sweep 0–0.5 |
| `RIFT_MAX_RESIDUAL_INCREASE` | 0.10 | As RFC 0003: a cut may fit the flow at most 10% worse | Sweep 0.05–0.25 |
| `RIFT_MIN_CONNECTED_SHARE` | 0.99 | Numerical: each half is one piece, allowing a few stray nodes | Fixed |
| `RIFT_GRAPH_NEIGHBOURS` | 8 | Numerical | Fixed |
| `KEEL_FULL_RATIO` | 2.0 | Cratonic keels reach about 200–250 km [[Pearson et al., 2021]](#ref-pearson2021) | Fixed |
| `CRATON_RIFT_RESISTANCE`, `CRATON_STRETCH_RESISTANCE` | 0.75 | #274 | Solution 4: tune to the loss left after 1 |

Tune with `bin/debug/compare_cratons.py` on seeds 7, 11 and 2024 to 1 Gyr, plus seed
456862899 to 3 Gyr, at default density, 5 Myr steps, one run at a time, against the matched
baseline that includes RFC 0009.

## Detailed design (provisional touchpoints)

- `backend/app/merge_split.py`
  - New `_rift_strength(plate)` and `_strength_min_cut(points, strength, cut_normal,
    spacing_rad, min_nodes, pole_separation_threshold, velocities)`, returning the half-A
    mask and the per-node distance to the cut, or `None` to keep the flow cut.
  - `maybe_split_plate`: call it after `_refine_split_cut`; pass the mask to `Plate.split`
    and the distances to `apply_failed_rift` and `cratons.record_rift_cut`.
  - New constants above, with comments citing this RFC.
- `backend/app/plates.py` and `backend/app/sparse_quad_patch.py`: `split` accepts a node
  mask (through `_plates_from_node_masks`); `apply_failed_rift` accepts per-node distances.
- `backend/app/cratons.py`: `record_rift_cut` takes distances.
- Tests:
  - `backend/unit_tests/test_merge_split.py`:
    - on a synthetic plate whose flow cut crosses a craton and a weak belt runs within the
      corridor, the cut follows the belt;
    - with uniform strength, the cut stays within a node spacing of the flow cut;
    - a cut that would leave a half under `SPLIT_MIN_NODES`, or raise the residual over the
      limit, falls back to the flow cut;
    - both halves are connected;
    - the result is identical for the same inputs (deterministic).
  - `backend/unit_tests/test_sparse_quad_patch.py`: `split` by mask preserves every field;
    `apply_failed_rift` with distances matches the normal-based version for a straight cut.
- Docs: `docs/simulation-model.md`, "Merge and split".

## Assumptions and uncertainty

- **The flow chooses the regimes.** The cut assumes the flow's two-regime split is right
  and only the line between them is free, within the corridor.
- **Strength is a proxy.** Craton strength, keel thickness and belt marks stand in for
  temperature, composition, water content and fabric. S only ranks cuts.
- **The kNN graph approximates boundary length.** Link density varies a little with the
  quad layout; the cost floor keeps cuts short but not exactly geodesic.
- **Rifts as whole cuts.** A real rift propagates over Myr and may stall part-way. Solution
  2 models that; this RFC still splits a plate in one step.
- **Without RFC 0009 the cut gains little.** The probe shows the craton blob blocks every
  route; this RFC's value depends on G1.

## Acceptance criteria

- [x] **Gate:** RFC 0009's G3 selects a branch that builds this RFC. It selected the
      minimum cut with resistance fallback (L = 27.6%, B = 0.732; #472, PR #479).
- [ ] Against the matched baseline on seeds 7, 11 and 2024 over 1 Gyr, both fall by at
      least 50%: the share of successful splits whose cut line (the nodes on cut links)
      touches craton, and the craton area in the bands of splits and failed rifts (the
      `split_band_craton_area_m2` and `failed_rift_band_craton_area_m2` counters
      `cratons.record_rift_cut` books).
- [ ] Craton loss to rifting, with solution 4 sized as needed, is at most 2% of live craton
      per Gyr on average across seeds 7, 11 and 2024, and in each of the first three Gyr on
      seed 456862899 (a design target, not an Earth rate).
- [ ] Splits per Gyr and mean plate count change by less than 10%; the mean fit residual
      of post-split plates rises by less than 10%.
- [ ] A cut attempt takes under 100 ms at 10⁵ nodes; step time changes by less than 2% at
      density 4.
- [ ] Ledgers balance; replays are deterministic; unit tests above pass; affected tests,
      ruff and `bin/quality_check.sh` are green.

## Falsification criteria

- **The flow puts regime boundaries through cratons.** If G1 falls in neither building
  branch even on a mosaic (L < 20%, or B > 0.75), the line was never the
  problem. Don't build 1; investigate why the flow's regime boundaries sit under cratons.
- **Curved cuts change the tectonics.** If splits per Gyr, plate counts or residuals move
  past the limits above at the smallest corridor that helps, fall back to solution 4.
- **Belts don't attract rifts.** If, with RFC 0009's belts, fewer than half of cuts that
  could follow a belt within the corridor do so, the strength proxy is too weak a signal;
  revisit `RIFT_BELT_STRENGTH` and the floor before tuning resistance.

## References

- <a id="ref-audet2011"></a>Audet, P., & Bürgmann, R. (2011). Dominant role of tectonic
  inheritance in supercontinent cycles. *Nature Geoscience*, 4, 184–187.
  [doi:10.1038/ngeo1080](https://doi.org/10.1038/ngeo1080)
- <a id="ref-boykov2004"></a>Boykov, Y., & Kolmogorov, V. (2004). An experimental
  comparison of min-cut/max-flow algorithms for energy minimization in vision. *IEEE
  Transactions on Pattern Analysis and Machine Intelligence*, 26(9), 1124–1137.
  [doi:10.1109/TPAMI.2004.60](https://doi.org/10.1109/TPAMI.2004.60)
- <a id="ref-foley2008"></a>Foley, S. F. (2008). Rejuvenation and erosion of the cratonic
  lithosphere. *Nature Geoscience*, 1, 503–510.
  [doi:10.1038/ngeo261](https://doi.org/10.1038/ngeo261)
- <a id="ref-heron2016"></a>Heron, P. J., Pysklywec, R. N., & Stephenson, R. (2016).
  Lasting mantle scars lead to perennial plate tectonics. *Nature Communications*, 7,
  11834. [doi:10.1038/ncomms11834](https://doi.org/10.1038/ncomms11834)
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
- <a id="ref-vauchez1997"></a>Vauchez, A., Barruol, G., & Tommasi, A. (1997). Why do
  continents break-up parallel to ancient orogenic belts? *Terra Nova*, 9, 62–66.
  [doi:10.1111/j.1365-3121.1997.tb00003.x](https://doi.org/10.1111/j.1365-3121.1997.tb00003.x)
- <a id="ref-wilson1966"></a>Wilson, J. T. (1966). Did the Atlantic close and then
  re-open? *Nature*, 211, 676–681. [doi:10.1038/211676a0](https://doi.org/10.1038/211676a0)
