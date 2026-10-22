# Gate 2 per-cell capacity-state audit

## Finding

The replay has enough per-cell state for a **diagnostic, thickness-aware column-capacity map**.
`replay_gate2_capacity.py` now evaluates each event node with its actual `Hc`, `Hm`, crust-type
proxy and Moho-lag scalar, then compares that capacity with the same node's shell `N_ext`. It also
reports by craton-crust share and inherited-belt marker. The replay stores compact per-event
stratum summaries rather than raw node arrays. A full rerun against pinned Gate 1 source reproduced
all 130 candidate time/plate identities and all 380 production cuts.

## Available per-cell fields

| Field | Current replay use | What it supplies and appropriate use/limit |
|---|---|---|
| `crustal_thickness_m` (`Hc`) | Used in each column-capacity solve and shell thickness. | Current crustal column thickness per node. It does not distinguish upper/lower crust or mineral proportions. |
| `mantle_lithosphere_thickness_m` (`Hm`) | Used in each column-capacity solve and shell thickness. | Current modeled mechanical/thermal mantle-lid thickness per node. It is a simulation proxy, not an independently observed lithosphere-asthenosphere depth. |
| `crust_type_code` | Selects continental, active-divergence young-ridge proxy, or mature-oceanic proxy. | Per-node categorical state, including mixed plates. It does not provide stratigraphy or oceanic age. |
| `craton_crust_m` | Used for candidate selection, craton-band coverage, and `craton_crust_m / Hc` strata. | Thickness-equivalent cratonic crust share. It is not keel thickness or mantle composition. |
| `belt_formed_years` | Used as a marked/unmarked inherited-belt stratum. | Whether/when a node was marked as a belt. It does not give fault orientation, damage, or quantitative strength. |
| `moho_thermal_lag_c` | Used in a second, explicitly assumed temperature-profile mode. | Scalar Moho temperature lag. It is not a transient depth-resolved `T(z)` profile. |
| `restite_m` | Available on the event-time plate, not collected by current scorer. | Per-node refractory crustal residue thickness. A composition/thermal proxy only, not a mineral mixture or flow-law parameter. |
| `divergent_age_myr` | Not collected by current scorer; used by tectonic evolution. | Clock accumulated while a node is divergent and reset otherwise. It is not oceanic crust age and cannot calibrate age-dependent cooling. |

The `Plate.age_years` value is plate-wide and also does not provide the age of each oceanic crust cell.

## Missing or incomplete state

- No full depth-temperature profile is stored. `orogeny.moho_temperature_c(Hc, Hm, lag)` returns
  a Moho scalar; the strength diagnostic can use its steady conductive profile, but applying the
  transient lag through depth requires an explicit interpolation assumption.
- No pore-pressure or fluid distribution, mantle water content, local grain size, mineralogy by
  layer, or semibrittle/Peierls rheology state is stored.
- No local strain-rate tensor, frictional fabric orientation, or yield/prestress `K0` field is
  stored. Shell `N_ext` is a diagnostic incremental load; it does not supply the lithosphere's
  background stress state.
- Continental crust is represented by broad type and thickness, not a measured upper/lower crust
  compositional split. Oceanic crust age remains unresolved even though divergent-boundary time is
  tracked.

## What the replay artifact can support

The matched replay covered all 74 source-connected candidates (23/23, 26/26 and 25/25 by seed),
with 56 genuinely disconnected-source candidates excluded from the solver. The 74 accepted cut bands contain 22,896 node samples. The first spatial summary used fixed `P_f/P_lith=0.42` and `phi=30 deg`; it is superseded by the
corrected hydrostatic-water and Byerlee-friction reference in `gate2-capacity-sensitivity.json`.
That re-score places 1,563 nodes (6.83%) above capacity under the hydrostatic profile, and 1,558
(6.80%) with interpolated Moho lag. The pooled node `N_ext/N_cap` ratio has median 0.111 and p95
1.248. These are conditional diagnostic comparisons, not
rupture rates or calibration labels.

| Seed | Connected candidates solved | Band nodes | Over capacity | Median node capacity (TN/m) | Median node `N_ext/N_cap` | p95 node `N_ext/N_cap` |
|---:|---:|---:|---:|---:|---:|---:|
| 7 | 23/23 | 5,942 | 189 (3.18%) | 4.65 | 0.116 | 0.750 |
| 11 | 26/26 | 9,713 | 1,038 (10.69%) | 11.89 | 0.146 | 1.503 |
| 2024 | 25/25 | 7,241 | 336 (4.64%) | 15.64 | 0.058 | 0.956 |

The thermal sensitivity is small for this replay because stored Moho lag is generally near zero;
it changes the pooled over-capacity count by five nodes under the corrected reference. The table is
node-weighted and uses the updated hydrostatic-reference re-score; see the sensitivity report for
all scenarios and source hashes. This does not make the capacity field
well constrained: the local profile still uses one assumed `K0`, pore-pressure ratio, strain rate,
cohesion, friction angle and representative creep laws. The young-oceanic class remains an active-
divergence proxy, not an ocean-age model. The prior split/failed-rift outcome is not used to fit
capacity.

Compact per-seed and pooled summaries are in
[`gate2-spatial-capacity-replay.json`](gate2-spatial-capacity-replay.json); full compact event rows
remain in `/private/tmp/issue487-g2-spatial-g1baseline-seed{7,11,2024}.json`.

## Sources inspected

- `bin/debug/replay_gate2_capacity.py`: reads `Hc`, `Hm`, crust type, craton share, belt marker
  and Moho lag for each candidate; uses `Hc + Hm` as shell thickness and evaluates local column
  capacity against the paired `N_ext` field.
- `backend/app/surface_fields.py`: registers `Hc`, `Hm`, `crust_type_code`, `craton_crust_m`,
  `belt_formed_years`, `moho_thermal_lag_c`, `restite_m`, and `divergent_age_myr` with their
  storage/remap semantics.
- `backend/app/orogeny.py`: computes steady Moho temperature from `Hc/Hm`, applies the scalar
  thermal lag, and documents that this is a proxy.
- `backend/app/lithosphere_plate.py`: advances `divergent_age_myr` only while a node is divergent,
  confirming that it is a boundary-activity clock rather than seafloor age.

## Reproduction / follow-up

The baseline model source is commit `0118752d0ffe75060464204d5b94f81a60fe9a07`. The diagnostic
was run from an archive of that commit with `thin_shell.py`, the replay/mesh scripts and only the
per-node force-field helpers overlaid; the baseline production torque routines were left pinned.
The 100 Myr seed 7 smoke and full event time/plate keys reproduced the prior Gate 1 population. The
full reports were written outside the repository to `/private/tmp`. Keep missing fluid, rheology
and stress-state quantities as scenario uncertainties; do not treat the mapped screen as a
calibrated rupture law.
