# RFC 0019 G3: coarse moisture transport and rain placement

Design and implementation handoff for [RFC 0019](../../rfc/0019-conservative-atmospheric-moisture.md)
and [#898](https://github.com/adubey/mantle-bloom/issues/898). G3 is still open and the
held-out saves remain sealed. The experimental sections record the original protocol. The
[proposed regional checks](#proposed-g3-regional-pattern-checks) reflect the owner's revised
scope. The registered numerical limits remain in force pending approval of the new checks.
See the separate [wind handoff](wind-transport-handoff.md) for possible changes outside this prototype.

## What the candidate must improve

On A and cold, averaging both density-4 maps to 2° reduced the initial best-fit RMS to
0.85–0.90, still well above the 0.50 limit even after giving the candidate a fitted global
scale. On A the initial control delivered too little water to high-latitude land; on cold it
delivered more than legacy while still producing much more dry land. The separate local-lift
treatment improved windward contrast but reduced rain poleward of 60° on both saves and in both
wind modes. The corrected density-4 solver also exceeded the time budget. See the
[#913 report](../issue913/README.md#cold-window-corrected-control).

Test regional delivery first, then place the available rain. A finer-looking map cannot repair
a dry regional budget. A change that improves mean land rain but shifts it into a narrow
coastal strip is also insufficient.

## Data flow and interfaces

1. Sample the native climate inputs once using the existing terrain, ocean mask, temperatures,
   winds and lagged vegetation inputs. Freeze the world and CFD state for each comparison.
   Keep the same mask and lagged cache across arms at a given density. Use
   `_cache_on_grid` in [the G3 harness](../../bin/debug/rfc0019_g3_compare.py) for the
   existing resolution comparison. Reject missing CFD state rather than measuring diagnostic
   fallback as a CFD result.
2. Build an aligned 90×180 transport grid at densities 1 and 4. Density 4 contributes 4×4
   native footprints per coarse cell. Derive cell areas from spherical latitude edges and sum
   source mass rates by tag (ocean, lake/sea, river, vegetation) and receiving reservoir.
   Aggregate the extensive rates, rather than averaging the humidity proxies or treating a
   partial coast as entirely ocean or land. Use the existing binary `is_ocean` from the
   native climate pass; inland waters keep their land-side source tags.
3. Assemble a stationary finite-volume system; solve for column water and calculate its rain
   sinks. Distribute those sink totals to native cells, then use that precipitation for the
   existing biome classification and climate diagnostics. Convert column water to the existing
   near-surface humidity interface separately.
4. The [offline prototype](../../bin/debug/rfc0019_coarse_moisture.py) provides pure array functions
   `aggregate_inputs`, `build_transport_rates`, `solve_budget` and `distribute_rain`.
   Reuse the reference control and metrics in `rfc0019_g3_compare.py`; keep each experimental
   arm independently selectable and retain the old controls for comparison. Reuse
   `transition_rates` for the density-1 zero-mixing parity fixture and
   `solve_stationary` as a small-grid reference, not its residual-to-rain correction as
   a way to pass the new ledger audit.
5. After a passing design and G4 review, factor shared input preparation out of
   [`compute_climate`](../../backend/app/hydroclimate/climate.py) and put the selected
   budget in a **new** `backend/app/hydroclimate/moisture_transport.py`. Use one
   moisture calculation and one final biome classification per update. Keep the
   `skip_moisture` caller in `world._advance_fluid_dynamics` off the moisture path.
   The first prototype may reuse the harness's duplicate sampling; the final timing
   must include real input preparation once, not subtract speculative savings.
   No new persistent atmosphere state is needed. Apply this plan to densities 1 and 4
   first; define exact spherical cell overlaps for nonnested grids or an explicitly
   supported separate grid path before enabling other densities in production.

## Transport and coastal accounting

Start with the corrected #913 control: 2 km effective depth, 9.1 days for producing a source
rate from the humidity proxies, 9.1 days for removal by rain, unit land source and transport
multipliers. Keep the converted sources fixed across experimental arms. These are G3
calibration proxies, not surface reservoir debits. A production water ledger must reconcile
lake/river evaporation with its existing hydrology writer and plant recycling with the
surface water available from the land budget. Avoid paying for one supply in two ledgers or
using the same update's rain as a new source. G3 is a stationary atmosphere screen; G4
and integration must audit the handoff to the surface reservoirs.

For coarse reservoir \(i\), let \(A_i\) be its area in m², \(W_i\) column water in kg/m²,
\(M_i=A_iW_i\) mass in kg, and \(S_i\) its fixed source in kg/s. Each shared face
has a length \(L\) in m and an eastward or northward signed normal wind \(u_n\) in m/s.
For the positive flow direction, the advective transfer is
`L * max(u_n, 0) * W_upwind` kg/s, with the negative direction handled by reversing
the donor and receiver. Use a fixed coordinate convention through aggregation and solve;
the production grid stores rows north to south and the historical prototype reverses them.
Wrap longitude; there are no open polar faces.

Add symmetric mixing with net transfer `G * (W_i - W_j)`, where
`G = K * L / d` (m²/s), `d` is center spacing normal to the face, and
`K` is a horizontal diffusivity in m²/s. Add two nonnegative rates
`G/A_i` and `G/A_j` (1/s) to the same sparse operator. Mixing tends to level
**column water per area**; leveling mass would create a spurious flow on a lat/lon grid.
Assemble shared faces once and audit signed fluxes with equal and opposite owner entries.
Start with isotropic `K = 30,000 m²/s` (design choice: roughly √(2Kτ) = 220 km in
9.1 days, comparable to a coarse meridional cell). Also report numerical spreading from
the no-mixing control; upwind transport already spreads moisture. A down-gradient term is
a bounded proxy for unresolved exchange, and need not reproduce a net poleward storm flux.

At density 4 give a mixed coarse cell two reservoirs, one for each nonempty surface class.
This keeps ocean and land source, rain, and coast transfer distinct. For coarse faces,
split the boundary into native face segments and sum signed normal advective capacities
and mixing conductance by the class pair on either side. Keep opposite-flow capacities
separate. Use native face lengths, coarse normal spacing for mixing, and divide outgoing
capacities by the donor's class area. Within a mixed cell, sum oriented wind capacities
across its native land/ocean edges once. Add symmetric local exchange between touching
classes using `G_coast = A_coarse * f_land * f_ocean / τ_coast` with a proposed
`τ_coast = 2 days` (design choice, fixed in both density-4 arms). This term represents
unresolved exchange inside a coarse cell; there is no second diffusive edge for those same
internal faces. Start with a pure cell per reservoir at density 1; its coastline exchange is
through regular faces. Test this approximation explicitly: disconnected patches of the same
class in a coarse cell share a column, and small coastal fractions can have large advective
rates. Reject cases where shortcuts create obvious coast artifacts or transfer changes enough
water to fail the density comparison. Keep the coastal link fixed when testing between-cell
mixing. Do not repair density disagreement by silently changing masks or reassigning rain.

Assemble `diag(λ + outgoing) - incoming` from vectorized edge lists, where the
initial positive sink is `λ = 1 / 9.1 days`. Start with a direct sparse solve on the
coarse grid. Apply the same operator to each source tag in offline audit runs, reusing the
factorization where useful. Reconstruct a raw cell balance in kg/s from face transfers and
`R_i = λ_i M_i`, independent of the matrix residual. Report maximum and summed absolute
cell residuals relative to total source, global closure, source-tag additivity, and net
ocean-to-land transfer. With no land source, land rain equals net ocean-to-land transfer
to tolerance. Gross inflow must be reported separately so repeated coast crossings cannot
inflate the proof. With zero global source, all mass and fluxes must be zero. Check the raw
solution before any roundoff handling; record any roundoff adjustment and re-audit. A solve
must never hide a bad residual by adding it to rain.

## Rain placement and the restrained legacy preference

First test uniform rain per unit area **within each reservoir**. With regional mass fixed,
compute an optimistic lower bound on fine-map RMS: aggregate candidate and reference
rain over each coarse cell and surface class, then evaluate the area-weighted RMS of
their mean differences using the original fine reference norm as denominator. If that
lower bound is greater than 0.50, no within-reservoir allocation can clear the existing
RMS gate at those totals. Examine upstream supply and belts before tuning fine detail.

In the first placement arm, concentrate some rain on windward terrain and give broad
wet belts a modest preference. Compute signed upslope exposure from a metric terrain gradient
and normalized local wind direction (neutral in calm cells). Smooth the terrain at a fixed
physical scale; use the existing 250 km smooth ascent scale from #913 as a starting point,
but apply it to **placement**, without adding a new rain sink. Use a cheap belt cue built from
the latitude row pattern of smoothed low-level convergence, as in the legacy model's
latitude blending, without running the legacy humidity sweep or using its full rain map.
For example, use `-div(u, v)` on the coarse sphere, blend across longitude by taking
the row mean, smooth over a fixed 5° span (design choice), and interpolate back to native
latitudes. Use bounded scores `t = tanh(exposure / s_t)` and
`b = tanh(belt_cue / s_b)`; set each scale with a documented rule (start with the
global area-weighted RMS of its cue on that snapshot; a zero scale gives a neutral cue).
Use the same rule across modes and worlds.

Candidate starting weights are `w = 1 + 0.35*t + 0.15*b`, so they are in
[0.5, 1.5]. For native cell `j` in reservoir `i`, assign mass rate
`r_j = R_i * A_j*w_j / sum(A_k*w_k)` over that reservoir's native cells. Empty
classes do not get rain; a zero reservoir sink gives zero rain. Publish `r_j/A_j`
as kg/m²/s, converted to mm/yr water equivalent for `ClimateFields.precipitation_mm`.
Use a uniform column water per reservoir initially when projecting back to native cells;
derive diagnostic `humidity_g_kg = 1000*W/(rho_air*h_eff)` with the same density and
depth used for the control. Report humidity range and any downstream sensitivity rather than
silently applying a cap or floor to the conserved water. Audit native sums against each
coarse sink, including separate coastal classes. Inspect block edges; if they are visible,
test smoother placement weights with the same exact totals before revisiting the grid.

The belt cue affects placement inside the available budget. It cannot move rain between
coarse cells or affect a one-cell reservoir at density 1. If the density-1 belts or coarse
regional supply are still wrong, consider a **separate** tuning arm with a bounded broad
rainout preference, `λ_i = λ_0 * (1 + β*b_i)`, starting with `β = 0.15`
(design choice). Keep it in the solve and repeat the cold-land check; no added rain after
the solve. Use the current rain maps as comparators, rather than as hidden source or sink
inputs. A later cold-land or direction-dependent transport hypothesis must earn a new
paired screen. If symmetric mixing erases dry belts before it improves delivery, try
a separately proposed meridionally stronger exchange only if flux and tag diagnostics
show insufficient cross-latitude delivery in both modes.

## Verification and decision order

| Stage | Compare | Decision |
| --- | --- | --- |
| Focused G2 extension | Existing density-1 control at K=0; analytic seam and pole, equal column water at unequal areas, unequal column water, zero wind/source cases, coastal tags, single-island mixed cells, area-conserving allocation and mirror orientation | Reproduce control within floating-point tolerance; transfer and sink audit at ≤1e-8 of total source, nonnegative finite mass and rain, exact zero-source response. Test coastal class totals and synthetic cross-resolution masks. Failure stops the G3 screen. |
| G3 transport screen | A and cold, both wind modes, density 1: K=0 vs 30,000 with identical sources and uniform sink. Controls and candidates are single climate snapshots. | Compare regional rain and import before choosing whether to test a higher resolution. Reject a setting that improves interior moisture only by flattening dry belts. Record cold-land delivery as a tradeoff; under the owner's current tuning priority, a cold-inland decrease alone does not reject a stronger shadow. A mix that only helps one mode is evidence for a mode-specific follow-up, not a pass. |
| G3 placement and coast screen | Same saves/modes at density 4; repeat K=0 and the selected mixing arm with uniform placement on the same coast geometry. Then compare uniform vs bounded placement on the selected arm with the same coarse budget. Record density-1 comparison; check the optimistic RMS bound before iterating on weights. | Require clear improvement in rain patterns with good coast allocation and no loss in dry or cold-land structure. If shape is still wrong at coarse scale, use flux/source tags to decide whether to test the broader sink or transport direction; don't keep tuning the fine texture. |
| G3 freeze and acceptance | If promising, use A, D, E and cold for tuning. A later scale-only comparison can use the existing 0.5/1/2 bounds, selecting at most one constant per wind mode across the tuning set. Then freeze code, constants and metrics and evaluate both registered held-outs, both modes and both densities. | Apply the existing G3 limits: ≤1e-8 closure, ≤5% global rain shift and ≤5 points land/ocean share across densities, registered land/ocean, latitude and Köppen bounds, raw RMS ≤0.50, median density-4 update ≤2× baseline. Clarify the latitude-band percentage-point normalization in the protocol before held-out evaluation. Never fit to held-outs. |

For each tuning comparison, publish common-scale rain and biome maps and a compact JSON table.
Use area-weighted, hemisphere-specific land/ocean band curves (belt peak latitude, width and
wet/dry contrast); windward/leeward pairs defined once from the paired control; land rain
by distance from coast and dry Köppen area (interior supply); and ≥60° land rain and the share
falling at below-freezing temperatures (cold-land supply). Show changes from both the paired
control and current-source climate, per mode. For cold land, distinguish a further shortfall
where the control is already dry from moving an excessive control total toward legacy; check
absolute supply and its geographic placement in both cases. Inspect coast and block seams on
the maps. These are design screens and stop signals, not replacement held-out criteria.

A 2° solve has 16,200 coarse cells before class splitting versus 259,200 density-4 cells.
G2 measured roughly 0.029 s for the coarse synthetic solve section; that does not include
native aggregation, mixed coasts or rain placement. Measure the complete climate update
against the matched current-source call in the same process, warm each path, record at
least three timed repetitions when uncontended, and report medians and peak RSS. Include
any repeat solve or legacy field used by the candidate in its cost. Record candidate/control
commit, save checksums, wind-state availability, source assumptions and commands. No world is
stepped in this screen. G4 is separate; it will measure downstream coupling and test the
registered geological step sizes or document why a second schedule is inapplicable.

G0's historical reproduction and G1's units and broad scale checks stand. The existing G2
result stands for its original operator; the extension above verifies the new operator and
transfer back to the climate grid. [#923](https://github.com/adubey/mantle-bloom/issues/923)
can inform transport hypotheses, but its canceled full-resolution run cannot supply
calibration constants or held-out validation. Review the tuning evidence before G3
acceptance or any production change.

## Focused G2 result

The [synthetic checks](../../bin/debug/check_rfc0019_coarse_moisture.py) passed.
Reproduce the [compact results](g2-coarse-transport.json) from the repository root:

```bash
backend/.venv/bin/python bin/debug/check_rfc0019_coarse_moisture.py \
  --out analysis/rfc0017/g2-coarse-transport.json
```

- With zero mixing on 90×180, the maximum rain difference was
  `5.9e-15` of peak reference rain. The sum of absolute cell residuals was
  `7.9e-16` of total source. The prototype calculates rain directly from the sink.
- With mixing and calm wind, equal column water across unequal cell areas
  produced zero net transfers within roundoff. A localized source exported across
  the longitude seam. The face list closed at both poles, the mirrored wind case
  mirrored its rain, and zero source returned exactly zero mass, rain, and flux.
- A mixed coastal cell with a one-cell island matched the analytic land rain
  from net import. The two internal wind faces matched their geometric capacities.
  A separate case with opposing winds along a coast checked land rain against
  import minus return flow and verified that the source-tag solutions added
  to the combined result.
- Weighted placement preserved each land and ocean sink within roundoff. A coast
  aligned to the coarse grid had matching coarse rain under 1×1 and 4×4 native
  subdivision. The test leaves the sensitivity to unresolved coast geometry on
  real worlds for G3.

These results cover the transport operator and the allocation of specified
weights. The next step is to connect the operator to the G3 harness and compare
zero mixing with the proposed mixing rate on A and cold in both wind modes.

## Paired G3 transport screen

Run [the paired driver](../../bin/debug/rfc0019_g3_mixing.py) once per pinned
[#913 tuning save](../issue913/README.md#recovered-inputs). It verifies each
save's hash before loading and uses the G3 harness to compute the reference
climate, source proxies and rain diagnostics. Within each wind mode the two
arms share one climate snapshot, one source field, a 9.1-day sink, and uniform
rain placement. Only diffusivity changes, from zero to 30,000 m²/s. This is a
90×180 screen of water delivery on frozen terrain.

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_mixing.py <pinned-tuning-save> --out <private-output-directory>
```

Compare the rain and biome maps, hemisphere bands, rain within 500 km of ocean,
rain from 500 to 1,500 km inland, rain at least 1,500 km inland, and land rain
above 60° latitude or below freezing. Use the same windward and leeward masks
within each pair. Check signed coast and 60° crossings to see where extra
rain comes from. Check that any improvement draws the candidate toward the
reference patterns without erasing dry regions. This screen does not measure
density-4 runtime or held-out acceptance.

## G3 result for 30,000 m²/s mixing

The first screen ran both registered saves in both wind modes at 90×180.
The run used base commit `251af936` with the offline prototype in this
worktree. The [compact results](g3-mixing-screen.json) include the save, script,
and result hashes. Full JSON and rain and biome maps are in the private captures
`g3-mixing-a` and `g3-mixing-cold`, alongside verified `SHA256SUMS` files.
Peak RSS was 1.07 GB on A and 1.64 GB on cold. The two probes ran sequentially
after the normal compute-budget check cleared.

| Save and wind | Raw rain RMS relative to current source, control → mix | Land rain ≥60°, mm/yr (current source) | Land rain ≥1,500 km from ocean, mm/yr (current source) | Windward/leeward ratio (current source) |
| --- | ---: | ---: | ---: | ---: |
| A, diagnostic | 1.749 → 1.625 | 114 → 116 (134) | 135 → 137 (200) | 1.40 → 1.30 (4.11) |
| A, CFD | 1.764 → 1.644 | 90 → 145 (353) | 207 → 213 (381) | 1.39 → 1.34 (2.68) |
| Cold, diagnostic | 1.721 → 1.594 | 527 → 618 (292) | 413 → 415 (1,005) | 3.84 → 3.59 (4.50) |
| Cold, CFD | 1.717 → 1.583 | 299 → 465 (396) | 1,248 → 1,215 (1,376) | 2.14 → 1.87 (3.72) |

All eight arms closed globally and locally. The largest summed absolute
cell residual was `6.1e-15` of total source. At the tested setting, raw
rain RMS fell about 7–8% in all four pairs. The rain maps show reduced
concentration along the largest convergence lines, but the narrow bands
persist. Even after fitting one amplitude, relative RMS remains 0.84–0.86.

Mixing raised rain near the ocean in all four pairs. Changes beyond 1,500 km
were small on three pairs and worsened the existing deficit on cold with CFD
wind. Dry Köppen share fell by 0.6–4.5 percentage points but remained well
above the current source in every pair. Mixing pushed high-latitude rain
in the same direction on all four pairs, despite different baseline errors.
It delivered more water to under-supplied high latitudes on A and made the
already wet high latitudes wetter on cold with diagnostic wind. Below-freezing
land rain moved toward current source on both cold arms, which shows that the
latitude and freezing masks reveal different effects.

**Disposition.** Stop this diffusivity setting before the density-4 placement
screen. The lower aggregate error comes with weaker windward/leeward contrast
in all four pairs and little correction to inland deficits. Fine rain
placement cannot change these broad delivery totals. The spatial pattern may
benefit from a different transport closure, but this test supports neither
stronger uniform mixing nor a move to held-out evaluation. A follow-up can
use the signed coast and latitude fluxes to identify which convergence bands
or routes need a different treatment before choosing a new arm. G3 remains open.

## Why the first transport candidate missed

Conservation worked. With a fixed source and a stationary closed globe, global
rain must equal global evaporation. Mixing can only change where that rain falls.
The captures show a low global mean and an unusually high spatial RMS at the
same time. Most cells are too dry while a small area gets very intense rain.

| Save and wind | Mixed mean / reference mean | Mixed median / reference median | Mixed RMS intensity / reference RMS intensity | Mixed land share error |
| --- | ---: | ---: | ---: | ---: |
| A, diagnostic | 0.79 | 0.32 | 1.92 | +0.01 pp |
| A, CFD | 0.76 | 0.23 | 1.93 | −1.04 pp |
| Cold, diagnostic | 0.75 | 0.25 | 1.86 | +7.65 pp |
| Cold, CFD | 0.74 | 0.23 | 1.85 | +6.78 pp |

RMS intensity is the area-weighted RMS of each individual rain map, which
is different from the RMS of their difference. Its ratio can be recovered
from the recorded best-fit scale `a` and fitted relative error `e` as
`sqrt(1-e²)/a`, using the identity for the least-squares projection.
Even the best uniform scale leaves a relative error of 0.84–0.86.
On A the land share is already close while much of the interior is too dry.
On cold the land share is too high by more than the registered five-point
allowance, while the far interior is still too dry. A uniform source change
cannot correct both patterns.

The strongest candidate explanation is that the new budget directly converts
sharp convergence of the near-surface wind into accumulated column water and
then rains it out everywhere on the same timescale. The diagnostic wind
switches meridional direction at the circulation boundaries. The CFD relaxes
toward that same broad wind target. The legacy climate smooths and bounds the
rain contribution of convergence, combines it with a broad humidity baseline,
and separately raises rain on windward slopes. This difference is consistent
with the narrow peaks in both candidate maps and why uniform mixing broadened
them only modestly. The maps and summaries do not quantify how much of the
remaining error comes from those sharp wind transitions, the source map, or
the uniform sink.

## Select the next candidate by measuring what the budget can deliver

1. On A and cold in both wind modes, retain the density-1 source, rain,
   column water and signed face flux arrays. Partition raw rainfall squared
   error by latitude bands, land and ocean, and distance from the coast.
   Report the share of error in the wettest rows and the errors remaining
   after aggregation. Integrate source minus reference rain over each latitude
   cap to estimate the flux it would need. For the flux calculation alone,
   normalize the reference to the same global source so an unmatched total
   cannot appear to require a flux through a pole. Keep the raw reference for
   G3 scoring. Compare the required flux with actual transport across the
   same faces and report coastal and inland deficits separately.
2. Bound what rainout alone could do before choosing another sink law.
   With sources and transition rates fixed, the stationary ledger determines
   `R = S - L M`, where `L M` is net export. Minimize raw area-weighted rain
   squared error over nonnegative mass with
   `M / (18 days) ≤ R ≤ M / (5 days)` and the registered land rain share
   constraint. Each rain rate in the objective is `R / area`. This is an
   offline feasibility check with access to the reference. It gives an
   optimistic lower bound at the given resolution because it can choose each
   sink independently. An infeasible share constraint or a certified bound
   above the target rules out this bracket with those sources and transport.
   If it can fit, inspect what sink pattern it requires and test whether a
   simple rule derived from wind and terrain reproduces the gain. Use the
   fitted sinks only to diagnose what the fixed transport operator can deliver.
3. If the flux comparison confirms excessive convergence into narrow rows
   or inadequate delivery outside them, test one change to the **effective
   column transport**. Start by broadening the abrupt transitions in the
   zonal-mean meridional component over a fixed physical width, while retaining
   the local departures that route coastal and terrain flow. Hold sources,
   sink and fine placement fixed for the paired test. Measure how much
   squared error the change removes from the peak rows and whether it
   retains or improves cold-land delivery, dry interiors and rain shadows.
   Screen the same rule on both wind modes. If the source budget itself
   prevents the required totals or regional delivery, test a source change
   separately rather than fitting more transport parameters to compensate.
4. Advance only a candidate that passes the existing numerical and
   rain-pattern limits as a tuning screen on A, D, E and cold. At density 4,
   compute the lower bound on achievable fine RMS from the coarse land and
   ocean rain totals before trying placement weights. Then measure actual fine
   rain, biomes and the full update time, including all preparation and
   downscaling, against the 2× limit. Confirm the registered latitude-band
   normalization before evaluating the held-outs. Freeze the candidate and
   run the existing held-out protocol to establish a G3 pass.

### Follow-up protocol after Wonder review

Run [the diagnostics](../../bin/debug/rfc0019_g3_diagnose.py) on A and cold
in both saved wind modes. Use one fixed snapshot per save and mode. The
control uses zero added diffusivity, the original tagged sources and a
uniform 9.1-day rainout time. Each intervention uses the same control.
Raw RMS uses native spherical cell areas. For the diagnostic 4° and 8°
scores, average both fields by spherical area and normalize against the
reference at that scale. Longitude tiles the sphere and latitude bins keep
2° caps at both poles. Partition native squared error by signed 15° band
and ocean or land within 0–500, 500–1,500, and at least 1,500 km of the
nearest ocean cell center. Report area and rain mass by region. Sort
latitude rows by reference rain mass and include rows up to and including
the first row that brings their total to 10% of reference rain.

Report signed northward flux at every internal latitude face. Sum modeled
flux across that row of faces. Required northward flux is the rain minus
the source summed over the cap north of the face. Convert reference rain
to kg/s per cell before scaling by total candidate source divided by
total reference rain. Use the scaled rain only for this flux calculation.
Cap flux constrains the total entering the cap. The regional partition
shows which land and ocean regions account for the rain error.

The bounded rainout check uses the zero-diffusivity operator and the control
source. Constrain land rainfall mass to the reference land share ±0.05
times total candidate source, so the constraint is linear. The objective
is squared raw relative RMS, with the 0.50 RMS limit corresponding to 0.25
in that objective. Each tangent to the convex objective underestimates it
over the full feasible set. Minimize the tangent with SciPy HiGHS IPM using
primal and dual feasibility tolerances of `1e-8` in scaled units. Compute
a lower bound from its inequality multipliers and charge any residual
reduced cost against the maximum total mass implied by the 18-day bound.
Subtract an additional numerical margin and report the lower bound
separately from any upper bound. For an upper bound, clamp the inferred
sink rate to the interval, assign the 18-day sink to empty columns, and
independently solve and audit the stationary ledger and land share. Use
up to four tangents with a 15-second limit each. If the bounds straddle
the target or the LP reports infeasibility without a certificate, report
the result as inconclusive.

The transport arms apply a Gaussian to the meridional **cell** velocity
before building face capacities. Use latitude standard deviations of 4°
and 8°, three standard deviations of support, and nearest-row extension
at the poles. Take the row mean with equal longitude weights since cell
areas within a row are equal. Smooth the mean while retaining local
departures. At each width use a matched control that smooths the departures
and keeps the mean. Eastward wind, source, sink and diffusivity use the
control values. Compare both raw and peak-row squared error. Only a gain
over the matched control supports attributing the effect to zonal wind
transitions.

The evaporation arm weights ocean cells by area times wind speed times
the positive humidity deficit. Use saturation humidity at sea surface
temperature minus 0.8 times saturation humidity at model air temperature.
The fixed 0.8 relative humidity is a diagnostic assumption. Rescale the
ocean weights to the control's total ocean source and reuse each land tag.
This arm tests where evaporation occurs. A uniform transfer coefficient
would cancel in the normalization.

Consider a combined arm if both interventions reduce raw RMS across
the four pairs and preserve or improve the cold-land, interior, and
windward/leeward errors against the reference. Report area and total rain
alongside means, especially for the small cold far-interior region. These
are screening decisions. The preregistered G3 metrics and held-out
protocol govern acceptance.

### Follow-up results on A and cold

The four save and wind pairs completed sequentially at 90×180 after the
owner approved the compute override. Each pair ran the control, two
zonal-mean smoothing widths, two matched local-departure controls, and
the equal-total evaporation placement arm. Peak RSS was 1.06 GB on A and
1.68 GB on cold. The largest summed cell residual was `4.6e-15` of source.
The private captures `g3-diagnose-a` and `g3-diagnose-cold` contain the
full arrays and checksum manifests. All capture checksums passed.

| Save and wind | Control RMS at 2° | Control RMS at 8° | 8° zonal smoothing, RMS at 2° | 8° local-departure control, RMS at 2° | Evaporation arm, RMS at 2° | Rainout RMS lower bound |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A, diagnostic | 1.749 | 0.699 | 1.280 | 1.782 | 1.731 | 0.736 |
| A, CFD | 1.764 | 0.708 | 1.219 | 1.779 | 1.768 | 0.839 |
| Cold, diagnostic | 1.721 | 0.790 | 1.367 | 1.774 | 1.702 | 0.689 |
| Cold, CFD | 1.717 | 0.747 | 1.446 | 1.695 | 1.716 | 0.732 |

Aggregating both maps removed much of the apparent error. Broad differences
remain at 8° in every control. The two or three wettest latitude rows
accounted for 33–66% of native control squared error. The 8° zonal-mean
intervention reduced raw error by 16–31% and latitude flux RMSE by 16–18%.
It consistently beat the matched local-departure intervention. It still
left large errors and sometimes worsened regional rain. For example,
below-freezing land rain on cold under CFD fell from 564 to 136 mm/yr
against a reference of 444 mm/yr. With diagnostic winds on cold, ≥60°
land rain rose from 527 to 892 mm/yr against a reference of 292 mm/yr.
The equal-total evaporation arm moved native RMS by at most 1.1% in
either direction and did not consistently help at coarser scales.

HiGHS returned status `Not Set` for all four tangent problems, so the
original LP bounds were inconclusive. A separate relaxation recovered a
lower bound from each archived zero-diffusivity control. With positive
bounded sinks and nonnegative transport, the solutions with uniform 5-day
and 18-day sinks enclose the mass of every solution with a sink inside
that bracket. The mass bounds and the stationary ledger give a range
for rain in each cell. Minimizing error over those ranges while enforcing
global water conservation gives an optimistic lower bound. The relaxation
allows each rain value to vary independently inside its interval and
omits the land share requirement, so it grants more freedom than the
proposed rainout change. It still exceeds 0.50 in every pair.

The archive calculation recovered directed rates on cells with positive
mass from the recorded advective flux and donor mass. Cells with exact
zero mass had no source or incoming flux, and none of the controls had
clipped mass. The recovered operator reproduced the archived stationary
ledger within `4.0e-15` of total source. The [compact bounds
and capture hashes](g3-rainout-bounds.json) can be reproduced with:

```bash
backend/.venv/bin/python bin/debug/rfc0019_rainout_feasibility.py \
  <private-dir>/g3-diagnose-a <private-dir>/g3-diagnose-cold \
  --out analysis/rfc0017/g3-rainout-bounds.json
```

The matched arms support addressing the concentrated zonal transport
first. The tested smoothing widths are insufficient for G3, the tested
evaporation placement does not warrant combining them, and rainout in
the 5–18 day bracket cannot pass with the tested control transport
and source. This evidence leaves other source distributions and
transport that can deliver water differently open for further design.

### Targeted convergence experiment

Use the zonal mean of the signed northward face velocity and multiply by
the circumference of that latitude to get the advective capacity in m²/s.
Work only inside the latitude windows from −68° to −52°, −8° to 8°, and
52° to 68°. These windows bracket the equatorial and subpolar wind reversals
already visible in the control. The drop in capacity from a row's south
face to its north face represents convergence. Redistribute each row's
positive drop inside its window with a latitude Gaussian of 4° standard
deviation normalized to conserve that row's contribution. Leave negative
drops in place and reconstruct the face capacity from the south boundary.
Apply the correction only to the internal latitude faces, retaining
local velocity departures and leaving eastward transport unchanged.
The edge constraint preserves the prescribed **air transport capacity**,
which does not fix water flux because the donor's column water can change.

On A and cold in both wind modes, compare the control, the previous
blanket 8° zonal smoothing arm, and this targeted treatment. Keep source,
rainout and diffusivity at the prior control settings. Score native RMS,
regional rainfall and windward/leeward contrast, and compute the bounded
rainout diagnostic for the targeted operator. Check that narrowing the
intervention preserves or improves the cold-land and inland reference errors
alongside any reduction in rain peaks. If the RMS lower bound remains
above 0.50, changing rainout within the tested bracket cannot rescue it
with these sources. The first screen checks density 1.

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_diagnose.py <pinned-tuning-save> \
  --belt-screen --out <private-output-directory>
```

### Targeted convergence result

Both tuning worlds completed sequentially after the owner approved the
memory override. Each world used its pinned save at 90×180. The controls
and blanket smoothing arms reproduced their archived RMS at all three
scales within `1e-12`. All twelve budgets closed, with a maximum summed
cell residual of `3.9e-15` of total source. Peak RSS was 1.04 GB on A and
1.60 GB on cold.

| Save and wind | Raw RMS, control | Raw RMS, blanket | Raw RMS, targeted | Targeted rainout RMS lower bound |
| --- | ---: | ---: | ---: | ---: |
| A, diagnostic | 1.7490 | 1.2805 | 1.5824 | 0.5080 |
| A, CFD | 1.7639 | 1.2190 | 1.5211 | 0.5666 |
| Cold, diagnostic | 1.7208 | 1.3670 | 1.5235 | 0.5063 |
| Cold, CFD | 1.7173 | 1.4462 | 1.6887 | 0.5856 |

The treatment reduced peak zonal convergence capacity by 58–61% in the
specified windows. Raw rain error fell by 2–14%, less than the 16–31%
reduction from blanket smoothing. The regional effects were uneven. On A
under CFD, land rain poleward of 60° rose from 90 to 266 mm/yr, closer to
the reference of 353 mm/yr. On cold under CFD it rose from 299 to
824 mm/yr against a reference of 396 mm/yr, while below-freezing land rain
fell from 564 to 173 mm/yr against a reference of 444 mm/yr. On cold under
diagnostic wind, rain at least 1,500 km inland fell from 413 to 392 mm/yr
against a reference of 1,005 mm/yr. Preserving air transport capacity at the
window edges still allowed the supply of water to those regions to change.

The rainout relaxation excludes the 0.50 RMS target in all four cases
with this operator, these sources, and 5–18 day sinks. The diagnostic
bounds are close to the target and the CFD bounds have a wider gap.
Combined with the regional regressions, this rules out advancing the
tested setting to density 4. It leaves open a change to the routes
feeding cold land and the interior. G3 remains open.

The private captures `g3-targeted-belts-a` and `g3-targeted-belts-cold`
contain full arrays, script hashes, and verified checksum manifests.
Their `results.json` SHA-256 values are
`7e533079cca05438cc0a50174a3f33dce0ab7ef7fd7835e71ac2fdc2efe91f30`
and `6448837cd580009cdafbd6ff05a9de2b187ff91e52e6635f5ce7b2052e916c1a`,
respectively. The base checkout was `251af936`.

### Moisture and transpiration attribution

The archived control and targeted arrays show accumulation in the water column
as well as concentrated precipitation. Divide `mass_kg` by `area_m2` for
column water, and `rain_kg_s` by `area_m2` for the sink. All eight maps
obey `P = W / 9.1 days` to roundoff. Columns and rainfall therefore share
the same spatial pattern by construction. The reference near-surface humidity
was not recorded in these captures, and would need its own comparison.

| Arms, all four save and wind pairs | Mean column water, kg/m² | Area-weighted median, kg/m² | Area-weighted 99th percentile, kg/m² | Rain falling in the wettest 1% of area |
| --- | ---: | ---: | ---: | ---: |
| Control | 21.9–25.2 | 5.7–8.2 | 411–420 | 18.4–20.1% |
| Targeted belts | 21.9–25.2 | 5.9–8.5 | 376–400 | 17.2–19.5% |

For the last column, sort cells by column water and sum rainfall mass over
the top 1% of area, prorating the boundary cell. The broad global scale is
similar to the short G1 ERA5 sample of 23.69 kg/m². That comparison does
not validate the spatial field. The targeted cold CFD maximum, 3,480 kg/m²,
would exceed the assumed air mass under the prototype's fixed 2 km depth
and 1.2 kg/m³ conversion. Treat it as a failure of that model representation.

The atmospheric transport uses one water column with the near-surface
diagnostic wind or CFD's single-layer wind as its velocity proxy. In a
vertically resolved budget, winds at different heights carry different
amounts of moisture. Varying winds and humidity within the averaging window
also contribute to the mean flux. Those omissions are plausible explanations
for the flawed routing. The paired runs isolate sensitivity to the prescribed
transport and have no vertical information to determine each effect's share.
See the [vertically integrated moisture budget in Wang et al.](https://www.aoml.noaa.gov/phod/docs/jcli-d-12-00284.pdf).

The vegetation term scales prior precipitation by a biome weight and a
recycling fraction before conversion to a source. It accounts for 1.38–1.74%
of total source in the A and cold controls. The sink and transport are linear
at fixed settings, so that same fraction bounds the direct amount of
plant-origin rain globally. Even assigning all of it to the wettest 1% of
area would explain only 7–10% of the rain there. Local recycling and coupled
feedback remain separate uncertainties. Earlier tuning increased all land
sources together with little improvement in rain shape. The separate
[source ownership](../../rfc/0019-conservative-atmospheric-moisture.md#physical-surface-sources-and-ownership)
requirement is to ensure that a surface reservoir pays for any vegetation
flux included in a production budget.

## Proposed G3 regional pattern checks

**Owner direction, 2026-11-17:** For the current tuning run, prioritize recognizable
rain shadows over matching cold-inland rainfall. A quantified cold-inland decrease from
the single-column approximation is an accepted modeling cost and is not by itself a stop
condition. Keep reporting cold-land supply and maps so the consequence remains visible.
This tuning priority does not revise the registered held-out criteria; freeze any replacement
acceptance rules before opening held-out saves.

Use the existing A, D, E and cold tuning cases to settle masks, aggregation,
and minimum meaningful contrasts before evaluating the held-outs. Derive
regions from the matched world's geography, temperatures and wind. Freeze
the rules before comparing candidates and use the same regions across arms.

| Pattern | Check on common-scale maps | Regional measurement |
| --- | --- | --- |
| Rain belts | Broad wet bands and adjacent drier zones are visible where the circulation supports them. | Plot area-weighted rain by latitude in each hemisphere with separate land and ocean curves. Measure peak latitude, belt width, and wet-to-dry contrast. Inspect for spuriously narrow convergence stripes. |
| Rain shadows | Along a prevailing moisture path, exposed slopes are wetter and sheltered slopes are drier beyond major barriers. | Use terrain and prevailing wind to pair windward and leeward regions. Report area-weighted precipitation and the distribution of pairwise contrasts so one very wet range cannot account for the signal. |
| Dry interiors | Large continents retain dry interiors where ocean moisture is depleted along the way. | Compare land rain by distance to coast and dry-biome area. Report sample area and stratify by latitude and exposure so a rain-shadowed coast doesn't define the inland signal. |
| Cold-land precipitation | Cold regions with moisture access receive visible precipitation, while isolated or sheltered cold interiors may remain dry. | Report area-weighted rates and total water equivalent on high-latitude and below-freezing land separately, with a fixed moisture-exposure mask. Check absolute supply and its placement, including the below-freezing share. Use G4 to measure snow and glacier consequences. |

Report precipitation and biome maps on shared scales for each wind mode and
density. Mark missing terrain or insufficient regional area and use the
registered cold or mountain-bearing cases to cover that pattern. Require a
recognizable signal and regional contrast where applicable in both modes.
Set numerical or map-review decision rules from the tuning cases and record
the owner's approval before opening the held-outs. Keep legacy-field RMS,
latitude bands, land/ocean rain and Köppen shares as comparisons. Keep the
registered tolerances for closure, density convergence and runtime in the
proposed gate.

First screen a bounded single-column candidate on A and cold at density 1.
Check maximum and upper-tail column water and the implied near-surface
humidity against local saturation, using the stated effective depth, air
density and temperature. Debit any excess converted to precipitation exactly
once. Hold source amounts fixed in the paired screen and account for water
through transport, rain and any placement step. If water reaching cold land
or the continental interior is insufficient, diagnose regional import before
adjusting fine placement. Advance a promising candidate to density 4 before
freezing it. More detailed vertical transport and plant recycling belong to
the deferred [#944](https://github.com/adubey/mantle-bloom/issues/944) follow-up.

### Saturation-bounded column screen

Compare the unbounded zero-mixing control with one change on A and cold,
in both wind modes at 90×180. Use the same pinned saves, cache, source
mapping, wind and 9.1-day first-order rainout. The local vapor cap is the
model's saturation specific humidity at the cell's annual air temperature,
converted using the same 2 km effective depth and 1.2 kg/m³ air density.
This is a deliberately simple near-surface conversion for the prototype.
It can restrict poleward vapor delivery too aggressively if a cold surface
is an inadequate proxy for the whole column.

In each saturated reservoir, assign excess supply after ordinary rainout
and net export to immediate local precipitation. In unsaturated reservoirs,
the overflow is zero. Solve from the cap downward, releasing cells whose
supplies cannot maintain their cap, and independently audit the face fluxes
and the total rain including overflow. The analytic directed-chain check
covers cascading release, exact overflow, no-source behavior and agreement
with the unbounded solve when the cap is inactive. No world is stepped.

Compare common-scale rain and biome maps, latitude profiles, the existing
strong-ascent/descent masks and adjacent downwind pairs, coast-distance
groups, and below-freezing precipitation within and beyond 500 km of oceans.
The cold-coast split measures proximity to a possible moisture source.
The shadow pairs are a coarse diagnostic of adjacent terrain, with their
sample count reported. If the cap eliminates the pileup but creates narrow
rain bands or strips moisture from cold land and interiors, keep it as a
budget bound and diagnose delivery before advancing to density 4. The
tuning result does not decide the formal G3 gate.

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_bounded_column.py <pinned-tuning-save> \
  --out <private-output-directory>
```

#### Archived transport check while the save runs wait for memory

An inexpensive sensitivity check reused the four 90×180 unmodified-wind
controls from the targeted-convergence captures. It divided each archived
advective vapor flux by the upwind column water to reconstruct the signed
air transport capacity. It recovered the original control mass to `1e-10`
relative tolerance. Where zero donor mass prevented reconstruction, it set
capacity to zero. This covers 0–1,832 faces out of 32,220 per snapshot and
limits interpretation of the changed flow. The cap was uniform at 25 or
50 kg/m² because the archives lack air temperature.

At 50 kg/m², immediate local overflow supplied 57–60% of global precipitation
and rain on land poleward of 60° fell 25–63% relative to the unbounded
controls. At 25 kg/m² the overflow share was 63–67% and polar-land rain fell
28–76%. The cap alone reduces vapor available to travel beyond wind
convergence zones, so it can make regional delivery worse.

The same check projected the signed face capacities toward zero net air
transport per cell, using a metric-weighted graph Laplacian. Removing half
the divergent component increased polar-land rain in all four cases at the
50 kg/m² cap. Full removal kept even the **unbounded** column below
47 kg/m² in all cases, versus peaks of 669–2,036 kg/m² with unchanged
transport. It left some cold-land deficits and reduced the land rain share
on the cold save. The projection changes the air routing significantly and
is a diagnostic of potential leverage. This check has no terrain contrast or
temperature-dependent bound and cannot select a G3 candidate.

| Tuning case | ≥60° land rain with 50 kg/m² cap, mm/yr (original / half / full convergence adjustment) |
| --- | ---: |
| A diagnostic | 54 / 91 / 31 |
| A CFD | 68 / 123 / 75 |
| Cold diagnostic | 195 / 289 / 257 |
| Cold CFD | 178 / 372 / 293 |

The maximum summed cell residual in the archived comparison was below
`4e-15` of total source. Reproduce with the `g3-targeted-belts-a` and
`g3-targeted-belts-cold` directories described above.

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_archive_cap.py <targeted-belt-capture-directory>
```

#### Archived controls with the saved CFD temperature

Read only the saved CFD temperature and wind components from each pinned
save's array archive, using the production resampling to 90×180. Use the
archived sources, land mask and legacy rain. Both recomputed controls
reproduce the archived control masses to `1e-10` relative tolerance.
This measures the local saturation cap and the effect of adjusting
transport using the actual CFD winds, including faces with zero archived
water. The two source worlds are never stepped or fully loaded.

| Case and arm | Peak column, kg/m² | Overflow share | Below-freezing land rain, mm/yr | ≥60° land rain, mm/yr | Legacy relative RMS |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A legacy | — | — | 166 | 353 | — |
| A unbounded | 669 | 0% | 87 | 90 | 1.764 |
| A local cap | 53 | 64% | 65 | 60 | 1.744 |
| A half convergence adjustment and cap | 53 | 48% | 88 | 77 | 1.275 |
| A full convergence adjustment and cap | 43 | 3% | 105 | 76 | 0.554 |
| Cold legacy | — | — | 444 | 396 | — |
| Cold unbounded | 2,036 | 0% | 564 | 299 | 1.717 |
| Cold local cap | 53 | 63% | 75 | 108 | 1.719 |
| Cold half convergence adjustment and cap | 53 | 47% | 141 | 195 | 1.250 |
| Cold full convergence adjustment and cap | 47 | 5% | 259 | 231 | 0.556 |

The local cap prevents excessive column water, but it shifts most rain to
immediate overflow and reduces moisture delivery to below-freezing land.
The full adjustment reduces that shortcut and improves below-freezing
land supply compared with the cap alone. It still leaves both cases with
less cold-land precipitation than the legacy maps and cuts rain more than
1,500 km inland to roughly half of legacy on both cases. The remote region
is 42% of land on A and 3% on cold. All compared budgets close within
`4e-15` of total source. Rain belt and shadow checks, diagnostic-wind
confirmation, and density-4 runtime are separate.

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_archive_cap.py <targeted-belt-capture-directory> \
  --cfd-save <matching-pinned-save>
```

#### Paired local-cap result in both modes

A and cold completed sequentially at 90×180 after the owner's compute
override. Peak RSS was 0.87 GB on A and 1.57 GB on cold. The cap reduced
the four peaks from 669–2,036 kg/m² to about 53.4 kg/m². It converted
62–64% of rain to immediate overflow. Relative RMS changed little
(1.718–1.764 before, 1.714–1.744 after). Common-scale maps show that
narrow convergence stripes remain. All budgets close within `4e-15`.

| Case and wind | Below-freezing land precipitation, legacy / control / cap, mm/yr | Land ≥1,500 km inland, legacy / control / cap, mm/yr | Broad windward/leeward ratio, control / cap |
| --- | ---: | ---: | ---: |
| A diagnostic | 167 / 154 / 134 | 200 / 135 / 85 | 1.40 / 3.89 |
| A CFD | 166 / 87 / 65 | 381 / 207 / 207 | 1.39 / 1.73 |
| Cold diagnostic | 463 / 1,024 / 529 | 1,005 / 413 / 659 | 3.84 / 19.80 |
| Cold CFD | 444 / 564 / 75 | 1,376 / 1,248 / 1,323 | 2.14 / 2.52 |

The distribution of the rain matters. For the cold diagnostic case,
the cap moved an overly wet cold region toward the legacy rate. For the
cold CFD case it undersupplied below-freezing land severely. The larger
ratios for the strongest ascent and descent masks did not improve the
fraction of wetter windward sides in the small adjacent-pair sample
(6–28 pairs depending on the case). Reject the cap by itself as a pattern
candidate. Retain it as a budget bound while screening changes to transport.

Private captures `g3-bounded-column-a` and `g3-bounded-column-cold` contain
maps, captured fields, and verified manifests. The `results.json` SHA-256
values are `03ebd770882f8eb386e65a300bc888247286eef7228adf6381682d5a42a84cf0`
and `816152c7467b3741a5849959d9543245b7570e8ee5546298eb4df7d586a00932`.

#### Convergence adjustment with the local cap in both modes

Use the exact sampled diagnostic and CFD wind, temperature, terrain, and
source fields from the paired runs. Apply the same metric-weighted projection
used in the archived screen to the original shared-face capacities. The
zero-adjustment arm reproduces the captured capped rain to `1e-10` relative
tolerance. Half and full correction retain the same evaporation, cap, and
background sink. Only the moisture transport sees the correction.

| Case and wind | Legacy relative RMS, cap only / half / full | Below-freezing land precipitation, cap only / full, mm/yr | Land ≥1,500 km inland vs. land within 500 km of ocean, full, mm/yr |
| --- | ---: | ---: | ---: |
| A diagnostic | 1.717 / 1.230 / 0.568 | 134 / 122 | 103 vs. 824 |
| A CFD | 1.744 / 1.275 / 0.554 | 65 / 105 | 182 vs. 626 |
| Cold diagnostic | 1.714 / 1.237 / 0.602 | 529 / 429 | 581 vs. 916 |
| Cold CFD | 1.719 / 1.250 / 0.556 | 75 / 259 | 703 vs. 641 |

The full correction removes the thin wet stripes on the maps and cuts the
overflow share to 3–6%, with peak columns at 43–47 kg/m² and summed residual
below `4e-15`. It also weakens the distinction between wet and dry latitude
bands. The peak in the regional rain plot no longer represents an entire
row of excess precipitation. Most inland masks stay drier than the coast.
On cold CFD, the deep interior is wetter than the coast-wide average and is
only 3% of land. The checked lee pairs show only 0.02–3.8 mm/yr median
contrast across the four cases, so the signs alone do not demonstrate
recognizable shadows. Cold-land supply improves sharply in cold CFD and
moves closer to legacy in cold diagnostic. The A high-latitude land remains
undersupplied in both modes.

Advance the convergence-adjusted carrier to a **tuning** screen with a broad
rainout preference based on the original low-level convergence and with
terrain-weighted placement of the available local rain. Test the broad sink
against the adjusted carrier alone first. Check wet-belt width and contrast,
cold-land supply, and inland dryness before adding fine placement. Run the
placement comparison at density 4 to test shadows within each coarse budget
while preserving cold-land supply, mixed-coast totals, and runtime. Keep the
surface wind used for evaporation and exposure untouched. The current arms
do not establish a G3 pass.

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_archive_cap.py <bounded-column-capture-directory> \
  --bounded-capture --maps <private-output-directory>
```

The `projection/results.json` SHA-256 values are
`1317cd367108931bc3dd31f263786cb506e23919e19872630af79de1e80bfa20` (A) and
`9bd17c5e8607e1aba1e2afdd38566be61f9e258e315e2fd00a1cccaba3b7671a` (cold).

#### Broad rainout tuning screen

Compare the fully adjusted carrier and its existing local cap with two
prespecified rainout preferences on the four captured 90×180 fields. Use
the original (unadjusted) signed wind capacity on shared faces. Sum net
inflow by latitude row, divide by row area, and smooth the resulting
convergence rate with a Gaussian of 5° standard deviation. Normalize it by
its area-weighted RMS and bound it with `tanh`. Multiply the background
sink by `1 + β * score`, with the planned `β = 0.15` and a `β = 0.45`
sensitivity. The original map is a comparator. Hold evaporation, the
corrected transport, and the saturation cap fixed across the arms.

Compare area-weighted ocean rain at belt centers (within 8° of 0°, ±60°)
with adjacent waters 8°–20° away. Report center/flank ratios and latitude
profiles alongside the maps. Look for all three broad belts on both wind
modes. A gain of at least 5% in the belt/flank ratio that clears 1.05 is a
useful tuning signal. Reject a belt setting that reduces below-freezing
land precipitation by over 10% from the adjusted control, removes an
existing dry-interior contrast, or restores narrow stripes. These
diagnostic checks guide work on the tuning fields. They are not G3's
held-out acceptance limits.

Neither sink preference produced a useful belt signal. The stronger arm
changed the equatorial center/flank ratio by only 0.9–3.4% and left it at
1.13–1.18. Both subpolar ratios stayed below 1 in all four cases. The
legacy maps also lack a northern wet belt in A by this fixed ocean mask,
while the southern and equatorial belts are distinct. Below-freezing land
precipitation changed by −3% to +7% in the stronger arm. It preserved the
existing inland contrasts and did not fix the cold CFD reversal. All
budgets close below `4e-15`. Do not advance either sink setting to fine
placement. The comparison and maps are in `broad-rainout` under each
bounded-column capture. The `results.json` SHA-256 values are
`0a51361d05858d4dadf965dd37f8f1ccc62d1b24006e9657d3d29e8e89b505b3` (A) and
`71cf25ad911c36ee9aaf750977f6d75fa5eb60a0a0e21f9975434496b3e66dd2` (cold).

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_archive_cap.py <bounded-column-capture-directory> \
  --bounded-capture --belt-screen --maps <private-output-directory>
```

#### Broad moisture delivery tuning screen

Reintroduce only smooth zonal-mean meridional transport on the adjusted
carrier. At each latitude face, divide the sum of original northward
capacities by its latitude-circle length to recover zonal mean speed.
Smooth it by 5° standard deviation, limit the resulting speed to ±2 m/s,
and restore half or all of that broad speed to the corrected carrier's
northward faces. With a 9.1-day background sink, 2 m/s corresponds to
about 1,570 km of travel in the absence of convergence or other rainout.
Retain the uniform rainout, local cap, and source field. Every additional
transfer uses a shared face and debits its donor. This tests whether broad
moisture import can restore rain belts while preserving cold-land supply.
Keep the same belt and inland checks. Measure the fraction of rain made
by saturation overflow and inspect the maps for renewed narrow stripes.

The half-strength transport restored wet equatorial and southern ocean
bands in all four cases and the northern band on the cold world. The
equatorial center/flank ratio rose from 1.11–1.15 to 2.32–2.37, and the
cold CFD below-freezing land rain rose from 259 to 309 mm/yr. The half
arm cut below-freezing land rain on A CFD by 14%, outside the tuning
tolerance. Both strengths brought back narrow saturation-overflow
stripes. Overflow supplied 21–23% of rain at half strength and 37–38%
at full strength. The cold CFD remote interior remained wetter than the
coast, as it is in its legacy comparator. All budgets closed below
`4e-15`. Do not advance either original arm to density 4.

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_archive_cap.py <bounded-column-capture-directory> \
  --bounded-capture --belt-transport --maps <private-output-directory>
```

The `broad-transport/results.json` SHA-256 values are
`6da7a91bd56d4780f7eadd27f9c18938b1617ed9aadfede55343d987347951e1` (A) and
`117e426932bf1cbf8f7adba178b7e8f9cea9568df81f0d0c352904159ba83165` (cold).

#### Spread the saturation overflow before it reaches the surface

The broad air inflow can still end where the meridional velocity reverses,
so immediate overflow deposits much of the rain on one grid row. Test a
bulk approximation of geographic spread for that part of the rain only.
Leave the vapor solve and its local saturation debit unchanged. Assign
each cell's overflow to cells of the same surface class at the same
longitude, using target area and a Gaussian of 4° standard deviation,
truncated at ±8°. Remain within above-freezing or below-freezing land
for each land transfer. Normalize per donor, with the donor cell always
eligible. This books a separate redistribution of condensed water while
preserving the total and land/ocean and freezing-land budgets. The
background rain still falls where it leaves the vapor column.

Compare quarter- and half-strength broad transport with that placement,
holding each transport arm's atmospheric budget fixed for its paired
before/after map. Compare the existing belt ratios and the fraction of
each belt's rain in its wettest single row. Seek broad structure without
the stripes. Keep the previous cold-land and inland checks against the
fully adjusted carrier. The distance is a modeling choice for smoothing
unresolved weather and has no interpretation as the path of a droplet.

The half-strength spread arm reduced the wettest equatorial ocean row's
share of its belt from 29–31% to 17–18% across the four cases, close to the
legacy maps' 15–16%. The equatorial and southern belts remained wet on
both saves and the northern belt was distinct on cold. The quarter arm
had weaker subpolar structure. The half arm's relative RMS was
0.366–0.533, versus 0.594–0.706 with the same transport and local overflow.
Surface class closure was below `4e-16` of total rain. As intended, the
spread changed each freezing-land total by only roundoff compared with
its paired local-overflow arm.

Compared with the fully adjusted carrier, the half-strength arm still
reduces total below-freezing land rain on A CFD by 14%. Almost all of that
decrease occurs within 500 km of the ocean, where rain changes from 830 to
600 mm/yr against legacy's 627 mm/yr. Farther inland, it remains at 57–58
mm/yr against legacy's 136 mm/yr. The inland cold-land deficit stays open.
On cold CFD, cold-land rain rises from 259 to 309 mm/yr, still below
legacy's 444 mm/yr. Keep the half-strength arm as a tuning candidate for
placement. These coarse results do not qualify it for held-out testing.

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_archive_cap.py <bounded-column-capture-directory> \
  --bounded-capture --overflow-placement --maps <private-output-directory>
```

The `overflow-placement/results.json` SHA-256 values are
`876d4657ab208b4e75442e33b6a14f1714310c9e750bd042d51e46b06251293d` (A) and
`8cc756ddb3b11c5c4cdfcdf0e2f2781355a67c5aad5767ef3282e4d45cba0cb8` (cold).

#### Fine terrain placement check using saved elevation

Hold the half-strength coarse rain and the sampled wind fixed. Use the
pinned saves' cached elevation at 0.5° spacing (area-average the finer
cold cache). Repeat matched coarse surface and temperature classes across
their native children so no land, ocean or cold-land allocation can change.
The saved fine ocean masks differ from sampled coarse classifications,
including coarse cells without any matching fine children. This isolated
placement check uses the cached elevation and the coarse class masks.
A fresh density-4 climate comparison is needed to test mixed coasts,
fine temperature, resolution convergence, and total model cost.

Smooth elevation at approximately 250 km, then compute the signed upgrade
along the original wind direction. Use a `tanh` cue scaled by the
area-weighted RMS over land, neutral over ocean and under calm winds.
Compare uniform placement with `1 + 0.35 * cue` inside each coarse budget.
Measure adjacent strong-ascent/descent land pairs with compatible winds.
Report pairs sharing a coarse budget separately from pairs across coarse
boundaries. Each budget must deposit the same amount in both arms.

The terrain cue honored every coarse budget. The maximum residual for any
reservoir was below `1e-18` of total rain. Within a coarse footprint, the
cue made 56–64% of the tested windward cells wetter than their paired
leeward cells. The median differences were below 1 mm/yr in all four
cases, against median baseline precipitation of 67–462 mm/yr. Across
footprints the changes were mixed.

| Case and wind | Within-footprint pairs | Median windward minus lee with terrain, mm/yr | Median uniform rain in those pairs, mm/yr |
| --- | ---: | ---: | ---: |
| A diagnostic | 289 | 0.007 | 67 |
| A CFD | 299 | 0.119 | 143 |
| Cold diagnostic | 123 | 0.560 | 462 |
| Cold CFD | 150 | 0.048 | 118 |

This particular fine placement cue does not produce useful rain shadows.
A single coarse footprint does not span many ridge systems, and the
smoothed cue varies little between adjacent fine cells. Keep the coarse
half-strength transport and spread as the test control for the next arm.
Test windward rainout during passage across a ridge at the transport scale.
Require a matched windward and lee contrast over whole ridge crossings and
audit inland and cold-land delivery before evaluating fine rendering or
opening held-out cases.

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_terrain_placement.py \
  <pinned-tuning-save> <matching-bounded-column-capture>
```

The `terrain-placement/results.json` SHA-256 values are
`49230be7aefdd361a3074b8797a487c7d32e5cea1184a6a94e80ce3fa14be296` (A) and
`95912662dee5181be72c72d7716b3a6a762682d9de24566f0b0994094adabaa9` (cold).

#### Rainout during ridge crossings

On the four captured 2° fields, use the half-strength broad transport,
saturation cap, and overflow spread as the common control. Mark a ridge
crossing where the carrier enters a land crest from a land cell and
continues along the same grid direction into lower land. Require at least
250 m of elevation gain and 250 m of descent over these two face crossings.
Freeze these samples before changing precipitation. Test forward and reverse
transport, including wrapped longitude, with no transfer past the poles.

Intercept a fraction `p * (1 - exp(-rise_m / 1000 m))` of the vapor carried
up to each qualifying crest. Test `p` of 0.10 and 0.25 as design choices.
Debit the full original face transport from the donor, deposit the removed
share at the crest, and credit the crest atmosphere with the remainder.
This changes the steady vapor field and can change how much reaches the
lee. It differs from #913's cellwise sink proportional to standing water
and triggers only on a transport path that crosses a crest. Existing
background rain and overflow placement remain in each arm. Subgrid motion
of cloud water is outside the approximation.

Measure change in the crest-to-lee contrast, including median absolute rain
and the share of crossings where the crest gets wetter and the lee drier.
Compare with the paired control, with cold coast and cold inland rain
reported separately, alongside dry-interior and ocean belt contrasts. Audit
deposition and cellwise closure so rain on the crest gets counted once.
No world is stepped or held-out input opened in this comparison.

At `p = 0.10`, the median gain in crest-minus-lee rain was 14–164 mm/yr,
with a wetter crest and a drier lee in 73–86% of the measured crossings.
The stronger arm increased the contrast farther and removed more water
before it could reach cold inland land on A CFD. That region went from
57.2 mm/yr in the control to 43.2 and 40.3 mm/yr at the two strengths.
On cold CFD the 0.10 arm instead increased below-freezing land rain from
309 to 340 mm/yr. Ocean belts mostly held, except that the already weak
northern band on A CFD weakened farther. Cellwise closure was below
`3e-15` and the spread residual below `5e-16`.

Test one weaker strength, `p = 0.025`, on the same captures. The purpose is
to retain a useful ridge effect with less loss in the already-dry A CFD
cold interior. Hold the crest definition, background transport, cap,
and placement fixed. Compare the median absolute contrast and relative
gain over the matched control, and report both signs and magnitudes.

At `p = 0.025`, the median contrast gain was 4–44 mm/yr across the four
cases. The crest got wetter and the lee drier in 74–86% of crossings. The
weak arm still reduced cold inland rain on A CFD from 57.2 to 49.3 mm/yr
(14%), while preserving or slightly increasing cold inland rain in the
other cases. The A CFD cold coast stayed at 600 mm/yr. Reducing the
global strength alone weakens the shadows in the diagnostic modes while
still affecting the already-dry cold interior on A CFD.

Test one supply-sensitive arm against the original `p = 0.10` result.
Scale the interception fraction by `min(W_donor / 10 kg/m², 1)^2`, using
donor column water from the matched no-ridge solve. The 10 kg/m² break
and square are modeling choices. They reduce the extra rainout when air
already carries little water. They do not introduce a source or relax the
saturation bound. Use the common captured winds and geometry and report
whether this changes the effects on A CFD cold interiors while retaining
measurable rain shadows in both wind modes. Freezing the column values
keeps this tuning comparison linear; the two-solve cost would need to be
included in any density-4 runtime check.

| Case and wind | Ridge crossings | Median contrast gain at `p = 0.10`, mm/yr | At `p = 0.025`, mm/yr | At `p = 0.05`, mm/yr | With the dry taper, mm/yr | Cold inland rain, control / weak / mid / dry taper, mm/yr |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| A diagnostic | 1,019 | 13.9 | 3.8 | 7.3 | 0.03 | 96.0 / 95.6 / 95.2 / 94.8 |
| A CFD | 1,087 | 120.0 | 34.0 | 65.7 | 12.7 | 57.2 / 49.3 / 46.9 / 45.4 |
| Cold diagnostic | 506 | 31.8 | 7.8 | 15.5 | 2.1 | 373.0 / 372.9 / 372.8 / 371.9 |
| Cold CFD | 567 | 164.2 | 44.4 | 84.3 | 84.9 | 409.9 / 414.0 / 417.6 / 410.4 |

The intermediate `p = 0.05` arm retained a measurable shadow: median contrast gain
was 7–84 mm/yr, with crest wetter and lee drier in 73–87% of crossings. It reduced
A CFD cold-inland rain from 57.2 to 46.9 mm/yr (18%) and increased cold CFD
below-freezing land rain from 309.1 to 324.7 mm/yr. The same mode-dependent tradeoff
remains: this setting strengthens the shadow while worsening cold-inland delivery on A.
This tuning arm does not select a setting for density 4 or held-out evaluation.

The dry taper suppressed most of the contrast in the diagnostic modes and
still cut A CFD cold inland precipitation by 21%. The weak uniform arm
was a better tradeoff on these cases, but it still worsened the A deficit.
Keep the face-based rule as evidence that the coarse solver can produce
rain shadows. No tested strength is selected for density 4 or held-out
testing. The budget and direction checks passed, including an analytic
three-cell case (interception plus overflow) and both directions across
the longitude wrap.

The A CFD regional audit traces the loss to net delivery. Net atmospheric
inflow to below-freezing land more than 500 km inland fell from 12.09 to
10.03 million kg/s under the weak arm and to 9.03 million kg/s with the
dry taper. Under the weak arm, gross inflow fell by 5.72 million kg/s
and gross outflow fell by 3.67 million kg/s. Local evaporation was the
same at 2.10 million kg/s, and the contribution of overflow spread across
the regional boundary stayed near 0.615 million kg/s. The weak arm's rain
shortfall therefore comes from lower net moisture delivery. This audit
does not identify which upstream crossings caused the loss. Trace the
affected inflow routes before making another change to the ridge law.

```bash
backend/.venv/bin/python bin/debug/rfc0019_g3_ridge_rainout.py \
  <bounded-column-capture-directory> --maps <private-output-directory>
# Repeat with --weak-sensitivity, --mid-sensitivity or --dry-sensitivity for those arms.
```

The pinned A and cold captures contain `ridge-rainout`, `ridge-weak`,
`ridge-mid`, and `ridge-dry` maps and reports. The `results.json` SHA-256 values are:

| Run | A | Cold |
| --- | --- | --- |
| Ridge strengths | `469dcaaf2d6e0d7730069d0a1178c59b429066b2cb32be378693c5320179acd1` | `4dcc638878043b54555239f355bdc49e193eeb043201806c69909a020de63704` |
| Weak arm | `14d8d0a2eee859b896d7f29e41e226e1eb4411ad48cdb460724441d2b0dc2802` | `f4f5c3112c91b3312e8fef991131c0f7e9937569cc22ef396aacf33b3fc3cc6b` |
| Intermediate 5% arm | `1b98d19f40f7ff6f971686ff406d48d3f9440b94dff0db76f1d729045914e505` | `72be4d970b59d3c041dd4ee074c00b3f08cb7a9c400c237ae72570c10fcaf816` |
| Dry taper | `4a62e849a84ae560b9b179c8430a49b44adda1212b3db6b390026d56e69264ed` | `912afecc4ea9d05a672c2a4eb9149781df833020d8e675592083bb24260b2bca` |

These reports predate the separate cold-inland import audit, which uses
the same solver and inputs. The regional budget can be regenerated with
the same commands.

#### Superseded: four-window ridge-strength screen after reservoir-height correction

**Superseded by the corrected full-density screen at the end of this report.**
The earlier native placement added overflow on top of `budget.rain`, even
though the budget already included that overflow. The stationary budget closure
in those screens was valid, but their placed precipitation totals and maps were
not. Do not use their pattern or density-convergence results.

The first defect found on 2026-11-17 was that the ridge detector indexed
native-cell elevations with reservoir IDs. It was corrected to area-average
terrain elevation into the solver's land/ocean reservoirs. The resulting screen
below was later superseded by the placement correction above. The earlier
A/cold hashes are superseded for provenance. That screen tested all four
strengths on tuning worlds A, D, E, and cold at density 1 in both wind modes. No held-out
save was opened. Each row reports median crest-minus-lee contrast gain and the
fraction of crossings where the crest became wetter and lee drier, for
`p = 0.025 / 0.05 / 0.10 / 0.25`. Gains are mm/yr; cold-inland rates are
control to `p = 0.25`, also in mm/yr.

| Tuning case, wind | Crossings | Contrast gain at 2.5 / 5 / 10 / 25% | Wetter-crest / drier-lee share at 2.5 / 5 / 10 / 25% | Cold-inland rain: control → 25% |
| --- | ---: | ---: | ---: | ---: |
| A, diagnostic | 1,019 | 3.8 / 7.3 / 13.9 / 28.9 | 84 / 83 / 83 / 81% | 96.0 → 91.5 |
| A, CFD | 1,087 | 34.0 / 65.7 / 120.0 / 231.3 | 74 / 73 / 73 / 72% | 57.2 → 40.3 |
| D, diagnostic | 959 | 2.8 / 5.5 / 10.6 / 22.0 | 82 / 82 / 83 / 82% | 82.9 → 76.5 |
| D, CFD | 1,031 | 29.0 / 55.4 / 101.1 / 198.1 | 78 / 78 / 78 / 77% | 48.7 → 40.6 |
| E, diagnostic | 974 | 2.5 / 4.8 / 8.9 / 19.1 | 83 / 83 / 83 / 83% | 86.0 → 80.2 |
| E, CFD | 1,042 | 28.3 / 54.9 / 102.3 / 194.9 | 77 / 77 / 77 / 76% | 38.4 → 30.9 |
| Cold, diagnostic | 506 | 7.8 / 15.5 / 31.8 / 78.6 | 86 / 87 / 86 / 86% | 373.0 → 370.8 |
| Cold, CFD | 567 | 44.4 / 84.3 / 164.2 / 346.3 | 83 / 83 / 83 / 83% | 409.9 → 434.0 |

All tested strengths create the intended local shadow in each tuning case and
wind mode. `p = 0.25` has the largest measured contrast gain. Its cold-inland
rain decreases on A, D, and E under both modes, while it increases on the cold
world under CFD winds; these costs follow the owner's stated tuning priority.
Budget residuals remain below `2.8e-15`. This is a density-1 finite-volume
screen only: it does not establish density convergence, wind-arrow behavior,
the density-4 runtime, or held-out acceptance. The 25% arm is the leading
shadow-priority tuning candidate, not yet a frozen G3 candidate.

Corrected result hashes and artifact manifests are stored in the private
`ridge-reindexed` directories. The `results.json` SHA-256 values are:

| Tuning save | SHA-256 |
| --- | --- |
| A | `42a218527366e7a7b43d102d4392c8b98417dbfefca47f765e3e43bcc6e9c103` |
| D | `6bea4938d3d759cb274ca5d2b31f866f21d023a9a351dc33890650fcdca33a89` |
| E | `9d616a1ee7adf03a9b4a019dceba43e607bed2e4cdc05b3330363e32099ba3c6` |
| Cold | `9340e5746d3fe66d10d653b80cce18903e3d9ca5ed677b823a74acbc8543574f` |

#### Superseded: density-4 A ridge-strength screen

The following density-4 results also double-counted overflow during native
placement. Their ridge detector and stationary closure checks were valid, but
the precipitation totals and maps are superseded by the corrected full-density
screen below.

On 2026-11-17, the A capture was recomputed at density 4 (360×720), then the
ridge-strength screen was evaluated with both solvers and wind modes. For the
finite-volume solver, at `p = 0.025 / 0.05 / 0.10 / 0.25`, median ridge contrast
gains were 0.61 / 1.21 / 2.45 / 6.63 mm/yr for diagnostic winds and 12.54 / 24.94 /
48.84 / 111.86 mm/yr for CFD winds. The wetter-crest and drier-lee fractions were
77.6 / 77.7 / 77.9 / 78.0% and 88.8 / 88.9 / 88.4 / 88.9%, respectively.

At `p = 0.25`, the wind-arrow solver also increased ridge contrast: 4.96 mm/yr
with diagnostic winds and 69.18 mm/yr with CFD winds. Its wetter-crest and
drier-lee shares were 81.4% and 88.4%, respectively. Cold-inland atmospheric
rain changed as follows:

| Solver, wind | Control → 25% cold-inland rain change |
| --- | ---: |
| Finite volume, diagnostic | −2.97% |
| Finite volume, CFD | −1.89% |
| Wind arrows, diagnostic | −2.14% |
| Wind arrows, CFD | −2.88% |

The summed-cell budget residual stayed below `2.7e-15` across the two solvers.

The `p = 0.25` A, D, E, and cold results are:

| Tuning case | Solver, wind | Median contrast gain (mm/yr) | Wetter-crest/drier-lee share | Cold-inland rain change |
| --- | --- | ---: | ---: | ---: |
| A | Finite volume, diagnostic | 6.63 | 78.0% | −2.97% |
| A | Finite volume, CFD | 111.86 | 88.9% | −1.89% |
| A | Wind arrows, diagnostic | 4.96 | 81.4% | −2.14% |
| A | Wind arrows, CFD | 69.18 | 88.4% | −2.88% |
| D | Finite volume, diagnostic | 4.03 | 80.1% | −3.25% |
| D | Finite volume, CFD | 90.16 | 88.9% | −0.46% |
| D | Wind arrows, diagnostic | 3.04 | 81.9% | −3.32% |
| D | Wind arrows, CFD | 59.91 | 89.8% | +0.31% |
| E | Finite volume, diagnostic | 4.24 | 79.2% | −3.11% |
| E | Finite volume, CFD | 103.87 | 88.7% | +0.20% |
| E | Wind arrows, diagnostic | 2.60 | 81.2% | −2.61% |
| E | Wind arrows, CFD | 62.02 | 89.7% | +0.05% |
| Cold | Finite volume, diagnostic | 0.20 | 61.5% | −0.82% |
| Cold | Finite volume, CFD | 194.41 | 79.1% | +0.64% |
| Cold | Wind arrows, diagnostic | 1.62 | 66.8% | −0.45% |
| Cold | Wind arrows, CFD | 93.85 | 80.6% | +0.83% |

These are four tuning worlds at density 4. All four have positive median ridge
contrast gain in both wind modes and solvers. Cold-world diagnostic winds have
the weakest directional share (61.5%). This does not establish the registered
density-convergence limit, runtime acceptance, or held-out acceptance. No held-out
save was opened. The capture, finite-volume comparison,
and two-solver comparison hashes for A are
`439336b38bd8af9e0eeabe2083b3d88f5919f4f627ded03ebe33ed9f5ac69519`,
`88c2fcd1b2c891873ad1e4faa105a3622ae17d9dec81ec567fb459edd45a6946`, and
`af9ab89e6a3ad350d2fcc74da7a32f21546e0d3fe75fe05429e18de42a65d4c6`.
The D capture and two-solver comparison hashes are
`005a510dd2289a820cb5dcbbe405f0a1d3d912795aa206d1c81657db7447591f` and
`9eaf2c4d2c772bb4017a461928e02ee93c1412002b7ea7a1fb3520d21e247343`.
The E capture and two-solver comparison hashes are
`2624234f1d5026f2c55c021790a3a2083344fccf10b47de6523eed8d4476699d` and
`c86f64a81883b7d2adab322f3bbc3d17fd3f7540e426a92f9f02975f5105e7e2`.
The cold capture and two-solver comparison hashes are
`4a6a065dc3f16ad301bee7960114c2b26ab920716f5290fba0d696fa34823ef3` and
`65dc13ab51482209d46e75013fcd5c13bfb925904612053b6817c6c790d72d42`.

#### Corrected full-density ridge screen after overflow-placement fix

On 2026-11-17, I found that `_place_surface_rain` spread `overflow` after
distributing `budget.rain`, although `solve_bounded_budget` had already included
that overflow in `budget.rain`. The placement now distributes
`budget.rain - overflow` normally, then spreads overflow within its surface
class. It checks both surface-class conservation and total placed rain. This
fix invalidates every ridge-placement result above this heading; the stationary
solver's own closure checks remain valid.

The current screen evaluates `p = 0.25` on all four registered tuning worlds,
both wind modes, both solvers, and densities 1 and 4. No held-out save was
opened. Results across all 16 world/mode/density combinations per solver are:

| Solver | Median ridge contrast gain | Wetter-crest and drier-lee share | Cold-inland rain change vs control | Raw RMS to legacy field |
| --- | ---: | ---: | ---: | ---: |
| Finite volume | +8.66 to +346.27 mm/yr | 72.5 to 89.2% | −30.9 to +5.4% | 0.364 to 0.781 |
| Wind arrows | +2.55 to +213.84 mm/yr | 63.4 to 89.9% | −7.6 to +2.3% | 0.732 to 0.945 |

Every case has a positive median contrast gain and a majority of paired
windward/lee crossings move in the intended direction. The finite-volume solver
creates stronger median contrasts; the arrow solver has smaller cold-inland
changes. The old held-out raw-RMS limit of `0.50` is still formal unless the
owner approves replacing it. This tuning screen does not meet it across all
cases, especially for wind arrows.

The registered density-convergence limits pass for `p = 0.25` on both solvers.
Across all worlds and wind modes, density 4 changes global rain by at most
`0.11%` and land/ocean rainfall share by at most `3.83` percentage points from
density 1. Maximum stationary budget residual is below `2.8e-15`; maximum total
placement residual is below `1.9e-14`.

The common-scale maps show an equatorial wet belt on A, D, and E. Its wet-center
to flank ratio spans `2.31–6.54` for finite volume and `8.43–8.76` for wind
arrows across modes and densities. On those same three cases, mean rain more
than 1,500 km inland is `0.13–0.42` of near-coast rain for finite volume and
`0.12–0.51` for wind arrows. The cold-source world remains a separate case: its
cold-inland rain changes are recorded above, without a minimum inland-rain
threshold, per the owner's stated priority. Rain-belt and dry-interior ratios
are tuning diagnostics, not owner-approved held-out limits.

| Density-1 case | SHA-256 | Density-4 case | SHA-256 |
| --- | --- | --- | --- |
| A | `be3064da88911bcd44ba25e43881e99ae8e6992e40935be7ca867eb9180d0d6d` | A | `a600b4b8b5aa46fa9751df6789205082d78751634ecfe05e791c13341be44e07` |
| D | `743ad283c4bf475dd1929d871a5c3cc4fb6070408a38dc7e63ecd804b686c93a` | D | `68c28fc89ddfa0cc1cef025517408c5b1fabce4a194363fd0e29f2265cf8eb88` |
| E | `e114407fde0927280f55ef8390555c7124311b3c778413d0bacff93c069fcf3b` | E | `ae34556c63c33ebdbfae799797802c03d45d95edf868fc2620bd24fbb82aaaa4` |
| Cold | `11e83240404d925de07fdf51a036e17563de29cba3eacce900a5459dddac995c` | Cold | `01efd4c61d719622198d7d3f80a44480148250df8eb53e51c83997c2329c60e4` |

The private corrected map set contains 32 common-scale images across four
tuning cases, two wind modes, two densities, and two solvers. Density-4 runtime,
the proposed pattern decision rules, and held-out acceptance remain open; G3
has not passed.
