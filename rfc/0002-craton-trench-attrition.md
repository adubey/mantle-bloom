# RFC 0002: Cratons resist sinking: margin-dependent attrition and collision sinks

- **Status:** Draft
- **Date:** 2026-10-13 (revised 2026-10-14 to cover collision sinks, and re-baselined on the
  year-based split pace, #403)
- **Related:** #340 (attrition at trenches), #177 (margin retreat budget), #320 (crust
  transfer), #321 and #331 (terrane docking), #299 (fault shear conservation), #358 (split cooldown in years), #403 (re-baseline); RFCs
  [0001](0001-merge-cap-craton-loss.md), [0005](0005-continental-recycling-balance.md),
  [0008](0008-craton-collision-jamming.md) (buoyant columns resist subduction)

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

The run behind RFCs 0001–0007 shows it. The run is seed 456862899, 19.2 Gyr, saved locally
as `mantle-bloom-seed456862899-19191800000y.mbworld`; the file is on the author's machine
and **not available for download**. Numbers come from `bin/debug/ledger_history.py` on that
save. The save predates #358: it ran post-#340 code (exact commit unknown) with the split
cooldown counted in steps, 200 Myr for most of the run, so it split and merged far less than
the model does now.

- `subducted_m3` was the largest craton sink: 1,128 M km³, **42%** of all craton
  destruction.
- It was steady rather than lumpy: no single 10 Myr sample took more than 0.5% of the live
  craton.
- Per Gyr, it took 6%, 18%, 7%, 12% and 10% of the live craton in the first five Gyr. Later
  windows ran at 20–77%, as the cratons shrank and had more edge per unit volume.

**On the current model** (origin/main 2562f63, 5 Myr steps; `analysis/issue403/`):

| Run | `subducted_m3` per Gyr, % of live craton | Rank among craton sinks |
|---|---|---|
| Seed 7, 0–1 Gyr | 5.9% | 4th, after `delaminated` (54%), `rifted` (13%) and `eroded` (9.7%) |
| Seed 11, 0–1 Gyr | 14% | 3rd, after `eroded` (21%) and `delaminated` (18%) |
| Seed 2024, 0–1 Gyr | 9.3% | 4th, after `delaminated` (31%), `rifted` (16%) and `eroded` (10%) |
| Seed 456862899, Gyr 0 / 1 / 2 | 15% / 22% / 33% | 3rd in Gyr 0 (after `delaminated` 28% and `eroded` 21%); **1st** in Gyr 1 and 2 |

It is still steady: no single step on seed 456862899 takes more than about 0.5% of the
live craton. Over the first Gyr, while continents are large, merges (RFC 0001) and erosion
outweigh it. Once merges stop (1–2 per Gyr on seed 456862899 after the first Gyr), it is the
largest sink, as in the save.

Geologically this is backwards. Archean cratons have survived for billions of years
[[Pearson et al., 2021]](#ref-pearson2021). Craton destruction is rare: the eastern
North China Craton is the textbook case, and it lost its *mantle root* under more than 100
Myr of Pacific subduction while its crust stayed at the surface
[[Zhu et al., 2012]](#ref-zhu2012).

The ledger also lumps three paths into `subducted_m3`:
- `cratons.attrite`, edge wear at trenches;
- the craton on cells `quad_tectonics._retreat` removes down a trench, which should be ~0
  because `cratons.retreat_allowed` keeps craton cells from retreating;
- the lost share at polarized collision fronts (`crust_transfer`, #320).

Until they are split, we can't say for certain how much is attrition.

**The collision sink.** A second account, `collision_reworked_m3`, books craton that "stays
on the surface as ordinary continental crust" (`cratons.py`). For the craton inventory that
is destruction: the crust survives, its craton status doesn't. Three sites book it:
- `crust_transfer._book`: the kept share of a consumed column at a polarized collision
  front (#320), and of a docking terrane (#321);
- `quad_tectonics._accrete_onto_survivors`: suture crust thrust onto the surviving cells;
- the phase audit after `faults.update_faults` (`world.py`), which clips any craton left above
  its column and books it here.

Reading the code narrows this down, but needs the split to confirm:
- The first two take their volume from cells that boundary retreat removed. Retreat never
  removes a cell that still holds craton (`cratons.retreat_allowed`), so they should book
  almost nothing.
- Fault relief (`faults._apply_plate_fault_relief`) moves Hc from donor to receiver cells, and
  moves continental material with it, but not craton. A cratonic donor's craton is left above
  its thinned column, and the audit clips it. That is craton destroyed by faulting with no
  collision involved, booked under a collision name. It also contradicts the rule that
  cratonic crust doesn't flow (`cratons.py`, `orogeny.py`).

On seeds 7, 11 and 2024 to 1 Gyr (origin/main at 1b4f686, reproduced exactly at 2562f63;
RFC 0008's evidence table, `analysis/rfc0008/`), `collision_reworked_m3` took 1.2%, 2.1% and 1.3% of the starting
craton. `subducted_m3` took 5.9%, 14% and 9.3%. The collision sink is small next to
subduction, but it is a conservation error where the fault audit is its source, so it is
fixed whatever its size.

**Cratonic terranes on oceanic plates.** Retreat can't remove their craton cells, and
attrition skips terrane cells (the `trench` mask in `_retreat`), so such a terrane jams at
the trench. It is resolved by docking it as cells (#331). The plate dynamics of the jam,
slab pull from a buoyant column, are RFC 0008's.

The run is a stress test of the model's long-term behaviour, not a calibration interval. It
shows that attrition, applied without regard to margin type, removes cratons on a Gyr
timescale. It does not give a rate to fit, and its late windows, when only small craton
fragments were left, say more about fragment geometry than about attrition.

## Proposed solutions

**What controls subduction erosion.** Margins are either *accretionary*, growing an
accretionary prism from offscraped trench sediment, or *erosive*, losing forearc crust
down the trench. [[von Huene & Scholl, 1991]](#ref-vonhuene1991) and
[[Clift & Vannucchi, 2004]](#ref-clift2004) found the split is set mainly by trench
sediment and convergence rate. Clift & Vannucchi find tectonic erosion favoured where
convergence exceeds 6 cm/yr *and* the trench's sedimentary cover is under 1 km, and
accretion favoured by slow convergence (below 7.6 cm/yr) and/or trench sediment over 1 km.
These are **statistical tendencies across margins, not a rule for any one point**:
- They classify a margin by its *net* behaviour over 10 Myr or more.
- Accretionary margins have short erosive episodes (Nankai, after seamount collisions),
  and erosive margins have short accretionary ones [[Clift & Vannucchi, 2004]](#ref-clift2004).

Erosion happens at the forearc, in the upper plate's crust above the slab. Cratons survive
in part because they sit behind younger, weaker mobile belts that take the deformation
[[Lenardic et al., 2003]](#ref-lenardic2003). Where cratonic lithosphere is
destroyed, the root is removed by refertilisation, heating and delamination
[[Foley, 2008]](#ref-foley2008) [[Lee et al., 2011]](#ref-lee2011)
[[Zhu et al., 2012]](#ref-zhu2012). Wholesale removal of the crust is not how it goes.

1. **Split the ledger first.** Add `attrited_m3` (from `cratons.attrite`) and
   `collision_lost_m3` (crust_transfer's lost and unplaced shares) to the craton ledger.
   `subducted_m3` keeps only what `_retreat` takes down a trench, and an assertion or
   diagnostic flags it as non-zero. Record attrition separately for continental and
   oceanic host plates, since only the continental case is capped by the #177 arc budget.
2. **Margin-dependent attrition, averaged and probabilistic.** Set each trench segment's
   attrition rate from its sediment and convergence, averaged over a margin-scale window in
   space and ~10 Myr in time, as the classification is defined. The rate moves smoothly
   between an accretionary floor and the erosive rate according to how strongly the
   averaged conditions favour erosion. An accretionary margin still wears slowly,
   standing in for its erosive episodes; it never goes to zero. Three variants are tested
   against each other (Design):
   - (a) binary on the averaged inputs;
   - (b) expected rate from a logistic erosive probability, recommended;
   - (c) stochastic regime switching per segment, with a minimum dwell time.
3. **Erode the root, not the crust.** At an erosive contact on a cratonic margin, attrition
   first thins the mantle lithosphere (Hm) of the exposed cells toward the continental
   reference. Only once the keel is gone does the cell lose cratonic status
   (`decratonised_m3`); its crust stays continental. This matches the North China history
   [[Zhu et al., 2012]](#ref-zhu2012) and keeps continental crust that today vanishes.
4. **A plain resistance factor.** Scale attrition by
   `1 − CRATON_ATTRITION_RESISTANCE × strength`, as erosion, rifting and stretching already
   do. It is simple, but has no physical basis for any particular value, so it's listed
   only as the fallback.
5. **Faults move craton with the crust, and cratons don't donate.** In
   `faults._apply_plate_fault_relief`, scale each donor's give by
   `1 − strength(craton)`: cratonic crust doesn't flow, as in `orogeny.relax_orogens`. Carry
   the craton share of whatever a donor does give to its receivers, the way continental
   material already goes (#276). The fault audit then has nothing to clip. What it still
   clips books to a new account, `faulted_m3`, not `collision_reworked_m3`.
6. **Cratonic terranes dock intact.** Docking a terrane by handing its cells to the
   overriding plate (#331) keeps its craton as craton. This RFC depends on #331 rather than
   duplicating it.

**Recommendation:** 1 as a diagnosis gate, then 2(b), with 2(a) and 2(c) as sensitivity
cases, and 5 as a conservation fix whatever the gate finds.
- Solution 1 decides whether attrition is the problem at all. If `cratons.attrite` turns
  out to be a minor part of `subducted_m3`, the RFC re-scopes (see Falsification
  criteria). It also splits `collision_reworked_m3` by site, which sizes solution 5.
- Solution 5 fixes a writer that breaks a stated rule. It's worth doing at any size, and
  first if the fault audit turns out to dominate `collision_reworked_m3`.
- Solution 2 is centred on the literature's thresholds, uses fields the model already
  carries, and gives RFC 0005 the classification it needs.
- Re-measure after 2. If the attrition share is still above target, add 3. It is the more
  physical answer, but it changes how cratons die and needs its own calibration.

**After the re-baseline (#403).** The recommendation stands. `subducted_m3` is no longer the
largest sink in the first Gyr (Problem), so this RFC follows RFC 0001 in the implementation
order. It is still the sink that takes over once merges stop, and it runs at 6–33% per Gyr
against a 3% target, so the case for it holds. `eroded_m3`, 10–21% of the craton in the
first Gyr, isn't addressed by this RFC or any other (#413).

## Design

**Units.**
- Attrition rate: km²/Myr (= m²/yr) of cross-section per unit edge length. A cell of area A
  loses `rate × years / sqrt(A)` metres of craton per step.
- Trench sediment: metres of sediment column on the subducting plate at the contact. The
  model's thicknesses are rock-equivalent (they are shares of Hc), while the literature's
  threshold is a physical, porous sediment column; see *Inputs* for the conversion.
- Convergence: the node's normal closing rate, `rheology.normal_closing_rate_m_per_s`,
  converted to cm/yr for the threshold.

**Phase 1: channel split (gate).** Land the ledger split alone, and run the tuning
seeds. Proceed to phase 2 only if `attrited_m3` is the largest part of what used to be
`subducted_m3`. Otherwise revise the RFC toward whichever channel dominates.

The split also covers the collision sink:
- the fault audit books to new `faulted_m3`;
- `collision_reworked_m3` gains diagnostic sub-accounts by site, outside the balance like
  RFC 0003's rift paths: `reworked_transfer_m3`, `reworked_docking_m3` and
  `reworked_suture_m3`.

If the transfer and suture sites book more than 0.1% of the initial craton per Gyr, craton
is reaching retreating cells by a route `retreat_allowed` doesn't cover; find it before
phase 2.

**Fault conservation (solution 5).** `faults._apply_plate_fault_relief` requests each
donor's debit as a rate times `years_myr` (`debit`, capped at the minimum crust), then
balances donors against receivers' headroom. Solution 5 changes the debit, carries craton
and its provenance, and integrates the resistance over the step.

- *Resisted debit, integrated.* A donor's craton leaves with its crust in proportion, so
  its craton share ρ = craton / Hc is constant while it gives. Its strength is then
  s = min(ρ × Hc / `CRATON_FULL_STRENGTH_HC_M`, 1), and it falls as the donor thins, so
  applying the starting resistance to the whole step depends on the step length. With r
  the unresisted debit rate (`debit / years_myr`, in m/Myr), the donor follows
  `dHc/dt = −r (1 − s)`. For s < 1 that has the exact solution

  `Hc(t) = [1 − (1 − a × Hc₀) × exp(a r t)] / a`, with a = ρ / `CRATON_FULL_STRENGTH_HC_M`

  and the resisted debit over the step is `Hc₀ − Hc(years_myr)`.
  - It is floored at the minimum crust, as `debit` is today.
  - With no craton (a = 0) it reduces to `r × years_myr`, today's debit.
  - A full-strength donor (s = 1) gives nothing, and stays at s = 1.
  - Example: r = 2,000 m/Myr, Hc = 35,000 m, craton = 12,500 m over 10 Myr. The starting
    resistance applied once gives 10,000 m; ten 1 Myr steps give ~11,389 m; the exact
    solution gives ~11,575 m for any partition.
- *Balancing.* The donor–receiver volume balance (`amount`) is then applied as today. When
  receiver headroom binds, the scaled-down debit depends on the partition just as today's
  does. This changes nothing about that.
- *Craton carry.* The craton a donor gives is `ρ × taken`. It is split among receivers in
  proportion to the Hc they get and capped at each receiver's new Hc. Units are metres of
  thickness per cell, balanced in volume with cell areas, as `continental_material_m`
  already is.
- *Provenance.* `craton_formed_years` is a HISTORY field (`surface_fields.py`), and
  `cratons.update` only dates craton newly forming on cells with none, so a carry has to set
  it:
  - a receiver that gains craton takes the oldest date among the donors that gave craton,
    unless its own valid date is older;
  - a donor left with no craton is reset to `CRATON_UNFORMED_YEARS`;
  - `stable_continental_myr` is left alone. The carried craton is already formed, and the
    formation clock only matters for crust not yet cratonic.

**Phase 2: margin-dependent rate.**
- *Segments.* Trench nodes on the overriding (craton-bearing) plate are grouped into
  segments: connected runs of trench cells facing the same neighbour plate, cut into
  pieces no longer than `MARGIN_SEGMENT_KM`.
- *Inputs.* Trench sediment is read at each node's nearest neighbour node through
  `inputs.neighbor_node_index` / `neighbor_plate_id`, as `torque.contact_craton_m` reads
  the neighbour's craton, and only where that cell is oceanic-coded. Continental-coded
  terranes are fragments, not sediment.
  - The model's sediment fields overlap: `mobile_cover_m` is a share of Hc, and
    `mobile_cover_continental_m` is a share of both `mobile_cover_m` and
    `continental_material_m`.
  - So the rock-equivalent column is built from two *disjoint* parts:
    `h_rock = mobile_cover_m + max(continental_material_m − mobile_cover_continental_m, 0)`.
    That is all loose cover, plus the consolidated continental sediment not already
    counted in it. A cell with 600 m of continental loose cover and nothing else reads
    600 m, not 1,200 m.
  - Consolidated non-continental (pelagic) sediment isn't tracked separately, so h_rock
    is a lower bound.
  - The threshold applies to a porous column, so `h = h_rock / (1 − φ̄)`, with φ̄ the mean
    porosity of a trench sediment column (`TRENCH_SEDIMENT_MEAN_POROSITY`).

  Convergence v (cm/yr) is `rheology.normal_closing_rate_m_per_s`. Each input is averaged
  along the segment, length-weighted.
- *Time averaging.* Per node on the overriding plate, the averaged inputs relax toward the
  segment values with time constant τ = `MARGIN_REGIME_AVERAGING_MYR`. Within a step from
  t₀ to t₁, with the segment value x held fixed, the average follows
  `x̄(t) = x + (x̄(t₀) − x) × exp(−(t − t₀)/τ)` exactly.
- *Erosive probability (b).*
  `p = σ((v̄ − v₀)/Δv) × σ((h₀ − h̄)/Δh)`, with σ the logistic function, v₀ = 6 cm/yr and
  h₀ = 1,000 m. p is in [0, 1] by construction.
- *Rate.* `rate = R_acc + (R_ero − R_acc) × p`, in km²/Myr. R_ero = 40 is today's
  `CRATON_ATTRITION_KM2_PER_MYR`; R_acc is the accretionary floor.
- *Integrating the rate over the step.* The rate is nonlinear in the averaged inputs, and
  they change during a step, so `rate(x̄(t₁)) × years` depends on step length.
  - Example: with constant v = 8 cm/yr and h = 0, starting from v̄ = 0, τ = 10 Myr and the
    starting knobs, one 10 Myr step gives ~197 km² of loss per unit edge length, and ten
    1 Myr steps ~148 km², against an exact integral of ~144 km².
  - The loss over a step is therefore `∫ rate(x̄(t)) dt` over [t₀, t₁], along the exact
    trajectory above, by `MARGIN_RATE_QUADRATURE_POINTS`-point Gauss–Legendre quadrature.
    The integrand is smooth (logistics of exponentials), so 8 points are ample for steps
    up to 10 Myr.
  - The loss is then capped by the cell's craton, as `attrite` does today.
- *Variant (a):* `p ∈ {0, 1}` from the averaged inputs.
- *Variant (c):* each segment draws its regime from p, deterministically per
  (seed, segment, epoch), and keeps it for at least `MARGIN_REGIME_MIN_DWELL_MYR`. Epochs
  are absolute time windows (`floor(t / dwell)`), not steps, so the realisation doesn't
  depend on step length. Segments are recomputed each step, so a draw is keyed by the
  segment's plate pair and the cell key of its midpoint on the overriding plate, which
  stays put as long as the segment does.

**Root-first variant (3), if needed.**
- At an erosive cratonic node, Hm thins at `CRATON_KEEL_EROSION_M_PER_MYR` until it reaches
  `REFERENCE_HM_CONTINENTAL_M` (100 km). The mantle removed books to the Hm ledger
  (`hm_ledger`, account `oceanic_and_deep_subduction`).
- Below that, the cell's craton converts to ordinary continental crust, booked to the new
  craton account `decratonised_m3`. Hc and `continental_material_m` are untouched.
- The thinning rate is set so a cratonic margin loses its keel's excess over the reference
  in ~100 Myr of continuous erosive subduction, the North China timescale
  [[Zhu et al., 2012]](#ref-zhu2012). For a 200 km keel that is about 1 km/Myr. The
  model's cratonic Hm is whatever seeding and deformation left, so the first step of this
  variant is to measure its distribution. If cratonic columns carry no keel excess, this
  variant has nothing to act on and is dropped.

**Knobs.**

| Knob | Start | Source | Tuning |
|---|---|---|---|
| h₀, `EROSIVE_MAX_TRENCH_SEDIMENT_M` | 1,000 m | [[Clift & Vannucchi, 2004]](#ref-clift2004) | Centre fixed; check the model's distribution of h at trench contacts first |
| v₀, `EROSIVE_MIN_CONVERGENCE_CM_PER_YR` | 6 cm/yr | [[Clift & Vannucchi, 2004]](#ref-clift2004) | Centre fixed; 7.6 cm/yr (their upper bound for accretion) is the sensitivity case |
| Δh, Δv (logistic widths) | 250 m, 1 cm/yr | Calibration assumption: the source gives tendencies with scatter, not a width | Sensitivity at ½× and 2× |
| φ̄, `TRENCH_SEDIMENT_MEAN_POROSITY` | 0.4 | Calibration assumption: the mean porosity of a ~1 km trench sediment column; no source is cited for it here | Sensitivity at 0.3 and 0.5 |
| `MARGIN_RATE_QUADRATURE_POINTS` | 8 | Numerical: enough that the integrated loss changes by under 0.1% when doubled | Check by doubling |
| `MARGIN_REGIME_AVERAGING_MYR` | 10 Myr | The ≥10 Myr window over which [[Clift & Vannucchi, 2004]](#ref-clift2004) classify margins | Sensitivity at 5 and 20 Myr |
| `MARGIN_SEGMENT_KM` | 500 km | Calibration assumption: the along-strike scale of margin segments in the compilations | Sensitivity at 250 and 1,000 km |
| `MARGIN_REGIME_MIN_DWELL_MYR` (variant c) | 10 Myr | Same window | — |
| R_ero, `CRATON_ATTRITION_KM2_PER_MYR` | 40 (unchanged) | Non-accreting margins [[Scholl & von Huene, 2007]](#ref-scholl2007) | Unchanged |
| R_acc, `CRATON_ATTRITION_ACCRETING_KM2_PER_MYR` | 12 | The accreting-margin value the #340 code comment attributes to [[Scholl & von Huene, 2007]](#ref-scholl2007); not re-verified for this RFC | Verify against the source before use; sensitivity at 0 and 20 |
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

**Runtime.** The inputs read three fields per trench node from a neighbour, with indices
already computed for the step. Segmenting is a connected-components pass over trench
cells. The quadrature evaluates the rate at 8 times per node. All are O(trench nodes),
about 1 ms per plate. Variant 3 adds an O(trench nodes) Hm update. There's no new search.

**Save compatibility.** Two new per-node fields, `trench_sediment_avg_m` and
`trench_convergence_avg_cm_per_yr`, are registered in `surface_fields.py` with blending
remaps. Old saves backfill them with a sentinel, meaning "no history". The first step
initialises them from the current values instead of relaxing toward them. New craton
ledger accounts (`attrited_m3`, `collision_lost_m3`, and `decratonised_m3` for variant 3)
are backfilled at 0.0 by `cratons.ensure_ledger`. A pre-change save keeps its history in
`subducted_m3`.

## Detailed design (provisional touchpoints)

- `backend/app/cratons.py`
  - Add `attrited_m3`, `collision_lost_m3` (and `decratonised_m3`) to
    `CRATON_SINK_ACCOUNTS`.
  - `attrite` books to `attrited_m3` and takes a per-node rate in km²/Myr, in place of
    the constant.
  - New `margin_attrition_loss(...)`: segments, the disjoint sediment column and porosity
    conversion, the averaged-input trajectory, p, and the quadrature-integrated loss;
    variants (a), (b) and (c) behind a module constant.
  - For variant 3, a `wear_keel` helper next to `attrite`.
  - Update the module docstring's Destruction paragraph.
- `backend/app/quad_tectonics.py`, `_retreat`: compute the per-node rate for the `trench`
  set before calling `cratons.attrite`, and update the time-averaged fields; split the
  attrition booking by host plate type. The #177 budget logic is unchanged.
- `backend/app/torque.py`: add a `contact_trench_sediment_m(plate, neighbours, inputs)`
  next to `contact_craton_m`. It reads `mobile_cover_m`, `mobile_cover_continental_m` and
  `continental_material_m` at the nearest neighbour node, on oceanic-coded cells only, and
  returns the disjoint rock-equivalent column. Reused by RFC 0005.
- `backend/app/surface_fields.py`: register `trench_sediment_avg_m` and
  `trench_convergence_avg_cm_per_yr` (new).
- `backend/app/lithosphere_plate.py`: expose the closing rate in the `BoundaryContext` if
  `_retreat` doesn't already have it.
- `backend/app/crust_transfer.py`: book the lost and unplaced craton shares to
  `collision_lost_m3` instead of `subducted_m3`.
- `backend/app/stats.py` and `cratons.diagnostics`: report the new accounts and the
  accretionary fraction of trench length.
- `backend/app/faults.py`, `_apply_plate_fault_relief`: the integrated resisted debit, the
  craton carry and its `craton_formed_years` provenance (solution 5).
- `backend/app/world.py`: the post-fault audit settles to `faulted_m3`.
- `backend/app/cratons.py`: `faulted_m3` in `CRATON_SINK_ACCOUNTS`; the `reworked_*`
  diagnostic sub-accounts, booked from `crust_transfer._book` (by caller: front or terrane)
  and `quad_tectonics._accrete_onto_survivors`.
- `bin/debug/compare_cratons.py`: report the new accounts per Gyr.
- Tests:
  - `backend/unit_tests/test_cratons.py`:
    - p is in [0, 1] and monotone in v and h;
    - a strongly erosive segment wears at 40 km²/Myr × edge, and a strongly accretionary
      one at R_acc;
    - accumulated craton loss (not just the final averaged inputs) agrees within 1%
      between one 10 Myr step and ten 1 Myr steps, under constant inputs, starting far
      from equilibrium (the v = 8 cm/yr, h = 0, v̄ = 0 case above);
    - the ledger balances with the new accounts;
    - an old ledger without them loads.
  - `backend/unit_tests/test_torque.py`: `contact_trench_sediment_m` reads the neighbour's
    nearest node; a cell with 600 m of continental loose cover and no other sediment
    reads 600 m; continental-coded cells read 0.
  - `backend/unit_tests/test_faults.py`:
    - a fault across a full-strength craton moves no Hc out of it;
    - across a half-strength one, the craton it moves arrives on the receivers, live craton is
      unchanged, and the audit clips nothing;
    - the transferred volume is the same, to 1e-9 relative, for one 10 Myr step and ten
      1 Myr steps under constant forcing with no caps binding (the 11,575 m example);
    - an ordinary receiver takes the donor's `craton_formed_years`; a receiver with an older
      date keeps it; an exhausted donor is reset to `CRATON_UNFORMED_YEARS`.
- Docs: `docs/simulation-model.md`, Cratons ("Resistance" and the destruction ledger), and
  the boundary-retreat description in "Quad-surface deformation".

## Assumptions and uncertainty

- **Attrition may not be the main path.** `subducted_m3` lumps three paths. Phase 1
  decides.
- **The classification is a tendency.** The thresholds come from a compilation of modern
  margins, classified by net behaviour over ≥ 10 Myr. The logistic widths, segment length
  and averaging window are calibration assumptions. Sensitivity runs and variants (a) and
  (c) bound how much the result depends on them.
- **The model's trench sediment is untested.** If it is far from Earth's (RFC 0005), the
  thresholds act on the wrong quantity and the classification is skewed. Hence the
  20%-per-regime sanity check. The column also misses consolidated pelagic sediment, and
  the porosity conversion is an assumption.
- **The rates are modern.** R_ero and R_acc are Cenozoic forearc rates. Over Gyr they are
  an assumption, not a measurement.
- **Target versus Earth.** The 3%-per-Gyr acceptance target is a model design target that
  lets cratons last several Gyr. It is not a measured Earth rate: no such rate is
  published for subduction erosion of cratons.

## Acceptance criteria

- [ ] **Gate (phase 1):** the ledger split lands, and on a 1 Gyr replay `subducted_m3` is
      below 0.1% of the initial craton volume, with `attrited_m3` and
      `collision_lost_m3` accounting for the rest. `collision_reworked_m3` is split into
      `faulted_m3` and the `reworked_*` sites. A note in `analysis/` records each channel's
      share before phase 2 starts.
- [ ] With solution 5, `faulted_m3` is below 0.1% of the initial craton per Gyr on seeds 7,
      11 and 2024.
- [ ] Craton loss to attrition is at most **3% of live craton per Gyr** over the first
      3 Gyr on a regenerated seed 456862899, and on average across seeds 7, 11 and 2024,
      for variant (b) at the starting knobs. On origin/main 2562f63, `subducted_m3` (an
      upper bound on attrition until phase 1 splits it) is 15%, 22% and 33% per Gyr on
      seed 456862899, and 9.9% on the tuning-seed mean (5.9%, 14%, 9.3%). The pre-#358
      save ran at 6–18%.
- [ ] The sensitivity runs (variants a and c; logistic widths and averaging window at ½×
      and 2×) are reported, and the 3%-per-Gyr result holds for most of them.
- [ ] Each regime (p < 0.5, p ≥ 0.5) covers at least 20% of trench contact length,
      averaged over a run.
- [ ] Genuine continental area at 1 Gyr is no lower than on `origin/main` for each seed.
- [ ] Ledgers balance (craton within 1e-9 of live volume; continental material closes).
- [ ] Step time on a density-4 world changes by less than 2% (`docs/profiling.md` method).
- [ ] Unit tests above pass; affected tests, ruff and `bin/quality_check.sh` are green.
- [ ] Saves from before the change load.

## Falsification criteria

- **Attrition isn't the sink.** If phase 1 shows `attrited_m3` under half of the old
  `subducted_m3` on the tuning seeds, the premise is wrong. Re-scope to the dominant
  channel; `collision_lost_m3` would point at crust_transfer's partition (#320).
- **Margin type doesn't matter in the model.** If p is near 1 almost everywhere (or near 0)
  on every seed, even at the sensitivity extremes, the model's trench conditions don't
  span Earth's range. The classification can't do its job until RFC 0005 or plate speeds
  are addressed; pause this RFC behind them.
- **Losses just move.** If total craton destruction per Gyr falls by less than half of
  the attrition reduction, the craton is being lost by another route once attrition stops.
  Then the edge exposure itself (no buffer of younger crust) is the problem, and the
  root-first variant or RFC 0005 comes first.
- **The result hangs on knobs.** If the 3%-per-Gyr result holds for some sensitivity cases
  and fails for others by more than a factor of two, the mechanism is underdetermined.
  Don't merge on the starting knobs; revisit with a physically derived rate.

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
