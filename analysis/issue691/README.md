# RFC 0014 G8: sediment delivery and outlet-sill placement

Exploration for [#691](https://github.com/adubey/mantle-bloom/issues/691), supporting G8
([#680](https://github.com/adubey/mantle-bloom/issues/680)). No model term was changed and no
coefficient was perturbed. This report does not change RFC 0014's flicker threshold or
authorize an implementation.

## Conclusions

1. **The #680 preflight rates are not sediment rates.** It read `silt_deposited` only. Its
   window A value (2.37 mm/yr) comes from one merge step whose silt supply is 711× the lake's
   external supply and 24× the whole world's. The neighbouring steps give about 0.003 mm/yr.
2. **The preflight measured the wrong connector in A and C.** Its A rim stayed 116 m above the
   water and did not spill. Its C rim lies about 400 m below the anchor's lowest water level and
   2,245 m below the highest level in the component.
3. **"No deposition on the sill" is true of silt only.** In B the spilling high node received
   166–191 m of routed sediment per step (1.7–1.9 mm/yr gross).
4. **Sill rise needs more than a deposition rate.** About 1/6–1/4 of a deposit survives Airy
   compensation, and in B a raised outlet moves to the next-lowest edge.
5. **The registered 0.4–8.8 mm/yr envelope should not be used as a sill range.** It mixes
   quantities and locations, and its source lake is about 90% autochthonous carbonate.

The G8 sediment question needs a redesigned diagnostic, and the silt supply count should be
settled first. See [Recommendation](#recommendation).

## Method

`bin/debug/measure_lake_connector_budget.py` replays each registered G8 checkpoint for five
100 kyr steps with `ice_age_period_years = 0`. The fifth step only supplies the next-step bed
comparison for step four.

Checkpoint SHA-256 values match the #680 registration. `backend/app` and
`backend/constraints.txt` are identical between the registered source
(`5d666ddf8b682a72db8c2189ebbd45a8a861bb20`) and this run's source
(`ea8f0ebb`). Normalized ledger residuals are at or below 1.4e-12.

```bash
NUMBA_CACHE_DIR=/tmp/issue691-numba backend/.venv/bin/python \
  bin/debug/measure_lake_connector_budget.py <window-A-checkpoint.mbworld> \
  --anchors analysis/issue691/anchor-window-a.json \
  --myr 0.5 --step-years 100000 --out /tmp/issue691-window-a.json
```

Run the windows one at a time; peak memory was 6.0–6.7 GB. Compact results are in
[`connector-budget.json`](connector-budget.json). The anchors are copied from the #680 corrected
anchor files, which are not yet on `main`.

The wrappers around `lakes._trace_transfer` and `lakes.solve_lakes_prognostic` only record
values and return the original results. The zero-depth reconciliation below comes from one step
of window A with `bin/debug/measure_lake_ice_phase_loss.py`, which observes
`hydrology.advance_lake_ice` the same way.

## Sediment pathways

| Pathway | Source | Spatial allocation | Limits | Terrain effect |
|---|---|---|---|---|
| Lake silt (`lakes.py`, `silt_deposited`) | `4.1e-6 ×` body supply. Not drawn from any eroded pool, so it is a non-conservative source | Evenly over the final wet area, by wet fraction | Each cell capped at its water depth | Folded into crust by `apply_erosion`, then Airy response |
| Routed fluvial sediment (`route_downstream`) | Eroded material | 15% retained where rivers are big and slow; ends at sinks | Mass conserved | Same |
| Lake spread (`_spread_lake_sediment`) | Routed sediment landing on a lake member | 50% uniform, 50% by depth over the wet component | Mass conserved | Same |
| Glacial, wind, landslide, flattening, marine and leveling | Their own eroded pools | Each pathway's own rule | Hc-cap overflow carried onward | Same |

`ErosionResult.sediment_deposited` holds everything except lake silt. The two are measured
separately here.

**Spill levels** come from the hierarchy's passage heights, which subtract breach notch depth.
Notch and channel depth were zero at every A and B connector recorded here, so that term did not
act in these windows. Deposition fills a notch before a channel, and this was not exercised.

**Airy compensation.** `apply_erosion` books the column change against crust and moves elevation
by the isostatic response. Only about 1/6–1/4 of a deposit survives as surface rise. The measured
net/deposit ratios are 0.165–0.244. Any clearance calculation must include this.

## Silt supply is overcounted on merge steps

`integrated_supply_at_leaf` adds every internal transfer to the receiving leaf, and a final body
sums its leaves. Water that a capped sibling pair passes back and forth is therefore counted as
new supply on every pass. See [Depth-cap recirculation](#depth-cap-recirculation). This
report first attributed the excess to a body's own children; the pass measurement corrected that.

| Window, step | Silt mean (mm/yr) | Component implied supply ÷ component external | Component implied ÷ world external | External-only silt (mm/yr) |
|---|---:|---:|---:|---:|
| A, 221.7 | 0.0032 | 1.0 | 0.02 | 0.0032 |
| **A, 221.8 (event)** | **1.878** | **711** | **24** | **0.0026** |
| A, 221.9 | 0.0027 | 1.8 | 0.03 | 0.0015 |
| B, 222.2 | 0.0061 | 1.0 | 0.21 | 0.0061 |
| B, 222.3 (event) | 0.0086 | 1.7 | 0.50 | 0.0052 |
| C, 201.5 | 1.061 | not measurable | 9.0 | not measurable |

A single lake cannot receive more water than the world supplies, so ratios above 1 in the
fourth column are an accounting artifact, independent of any Earth analogue. Component external
supply reads zero in C because its sinks fall outside the wet component, so C has only the world
bound.

In ordinary steps the component's supply equals its external supply, so the effect is confined
to steps where bodies exchange water repeatedly. The silt is still real to the model: it was about 100 m of the
125 m deposited on A's controlling rim in the event step, and that rim's spill level rose about
12 m. In C at 201.5 Myr, the same step as its silt spike, one dry cell received 5,301 m of
`sediment_deposited` for 899 m of net rise. It received only 12 m of silt, and the link between
the two has not been established.

World lake silt ranges from 2.9e13 to 1.2e15 m³ per step. The upper value is about 12 km³/yr of
non-conservative material.

## Depth-cap recirculation

`bin/debug/measure_lake_recirculation.py` counts the solve's passes and sending pairs over the same
five steps per window. Compact results are in [`recirculation.json`](recirculation.json).

**Mechanism.** A deep basin sits at its tier cap (1,800 or 3,500 m above its floor), below the
saddle it shares with a smaller sibling. The sibling is full to the saddle and spills into the
deep basin. The deep basin spills back and cannot reach the saddle, so the pair never merges. The
surplus moves on every zero-duration pass until the loop ends at one pass per leaf plus one.

- Passes reached the limit in 11 of 15 steps.
- All 42 top-listed pairs with at least 10 sends have a partner at floor + cap.
- Reciprocal pairs carry 94–98% of internal transfer in steps with large cycles. Internal
  transfer reaches 5.8e20 m³ per step against about 1.1e18 m³ of external supply.
- The guard at `lakes.py` lines 790–793 describes this case but does not stop these cycles. The
  branch that lets them through has not been confirmed.

A's event is one of these pairs: at 221.8 Myr, one body is capped at −831 m and the other is full
at −805.6 m. Each sends about 1.6e19 m³ over 989 passes, which explains its silt spike and its
two water levels. These are separate bodies, so a fix that only excludes a body's own children
would not remove this spike.

**Final export.** Surplus still deferred when the passes end is exported at a leaf sink, never the
ocean in these steps. It ranges from 0 to 33% of world external supply per step.

| Window, age (Myr) | Final export ÷ external supply | Ocean overflow (m³) | Sum (m³) |
|---|---:|---:|---:|
| B, 222.2 | 0.268 | 2.23e17 | 5.13e17 |
| B, 222.3 | 0.094 | 5.01e17 | 6.07e17 |
| B, 222.4 | 0.209 | 3.19e17 | 5.57e17 |
| B, 222.5 | 0.000 | 5.74e17 | 5.74e17 |
| C, 201.4 | 0.009 | 4.80e17 | 4.86e17 |
| C, 201.5 | 0.332 | 1.77e17 | 3.85e17 |

In B the sum is steadier than either term, which suggests the cycle may divert water that would
otherwise reach the ocean. A does not show this clearly, and the exported water's destination has
not been traced. Water ledgers still close, so G6 and G7 are unaffected.

## Controlling connectors

A model rim is a graph edge between cells 54–102 km apart. Cell relief reaches ±500 m, so a rim
is not a surveyed sill: A's rim cell has a mean depth of 100 m while its centre stands at the
spill level. Clearances smaller than that sub-grid range are not geomorphically resolved.

| Window | #680 connector | What the water uses |
|---|---|---|
| A | Rim −689 m, an 857 m step between adjacent cells, not spilling, 116 m above the water | Rim −806.5 to −805.6 m, which spills at the event (clearance 0.0). Spill level then −793.0, −792.8 and −794.6 m |
| B | Rim 664.8 m, spilling | Correct edge at the event. The spilling edge then changes each step: nodes 85220, 87566, 85293, 85243 |
| C | Rim −1,469 m, dry | No single level. The component joins cells from −1,302 to +776 m at the event, and the anchor's level runs 278, −1,065, 785, −1,059 and 285 m |

**A** is not sill-limited except at the event step. The surface runs −1,318, −805, −960, −1,346 m,
hundreds of metres below the rim in the other steps.

**B** is the one sill-limited case. Its high node gains 27–30 m net per step, yet the spill level
rises only 3.7 then 1.6 m, because overflow moves to the next-lowest edge. A sill must be
treated as the minimum over alternative routes, not as one node.

**C** is not one lake at one level. At the event, 2 of 58 cells hold 2,968 of its 3,195 km³. The
"component" is a thin film of wet cells around a perched body, so its `ΔV` does not describe a
single water surface and a single connector cannot be assigned to it.

### Sill rise needed for the governing 10% rule

This is a geometric estimate (`ΔV × 10% ÷ event area`) that ignores route switching.

| Window | 10% of event `|ΔV|` (km³) | Sill rise (m) | Deposit at 1/4–1/6 survival (m per 100 kyr) | Equivalent (mm/yr) |
|---|---:|---:|---:|---:|
| A | 49,562 | 79 | 316–473 | 3.2–4.7 |
| B | 13,162 | 4.9 | 20–29 | 0.20–0.29 |
| C | 26,159 | not applicable | not applicable | not applicable |

B's existing outlet deposit already exceeds its value, and its sill still moved less than 4.9 m
per step because of route switching. Other terrain change is the same size as the sill signal: terrain change
outside erosion reached about 9 m per step at the recorded connectors.

## Field evidence

| Source | Value | What it is | Use here |
|---|---|---|---|
| [Rostási et al. 2022](https://doi.org/10.3389/feart.2022.1067105), Lake Balaton | 0.75–0.9 mm/yr | Lake-wide, consolidated, modern | Basin-wide |
| Same | 0.4–0.5 mm/yr | Lake-wide, consolidated, 12–15 kyr mean | Basin-wide, long-term |
| Same | 1.5–1.7 mm/yr | Lake-wide **wet sludge**, about 70% water | Not comparable with consolidated rates |
| Same, citing Pálfi et al. 2021 | 2.1 mm/yr | Mid-lake core, about 150 yr | Local |
| Same, citing Zlinszky et al. 2010 | ~8.8 mm/yr | Seismic, eastern lake, since ~1980 | Local, decadal |
| [Colman et al. 2003](https://pubs.usgs.gov/publication/70121966), Lake Baikal | <0.03 to >0.3 mm/yr | Holocene, ridges to basin floors | Over an order of magnitude of spatial variation in one lake |

Rostási et al. attribute about 89% of Balaton's modern sediment to autochthonous carbonate and
5–6% by mass to streams and the atmosphere, and report sediment as practically absent in the
strait of Tihany. That is the lake's own constriction, so this source does not support aggradation
at a connector. The #680 "lake-wide 1.5–1.7 mm/yr" entry is correct as to source but is a wet
thickness.

Measured rates also fall as the measurement interval grows, the
[Sadler effect](https://en.wikipedia.org/wiki/Sadler_effect), so decadal rates should not be
extended to a 100 kyr step without an interval correction.

For outlets, the sources found describe sills as erosional or bedrock-controlled, for example
[Owens Lake](https://eprints.whiterose.ac.uk/161025/1/Bacon_et_al_2020_Late_Pleistocene_Owens_Lake_overflow_GOA.pdf),
where a hard bedrock sill set the spill level. No source reviewed gives a transferable outlet-sill
aggradation rate. Treat this as absent from what was searched, not as evidence that none exists.

For scale only, 19 Gt/yr ([Beusen et al. 2005](https://www.pbl.nl/en/publications/estimation-of-global-river-transport-of-sediments-and-associated-particulate-c-n-and-p),
95% CI 11–27) over about 30,600 km³/yr ([Cohen et al. 2022](https://meetingorganizer.copernicus.org/EGU22/EGU22-6556.html))
gives about 620 g/m³, or 2.3e-4 by solid volume. `SILT_ACCUMULATION_COEFFICIENT` is about 57×
smaller, as expected for a value calibrated to model behaviour (#363). These figures come from
separate studies, so this is an order-of-magnitude comparison.

## Placement hypotheses

| Hypothesis | Against A/B/C | Discriminating prediction |
|---|---|---|
| Keep the lake-bed blanket | Places nothing on a dry rim, so cannot change a sill by itself | Fixing the supply count leaves A and B spill levels within about the non-event drift |
| Inlet or delta-focused deposition | Top inlets lie 56–2,359 km from the rims, 13 of 42 of them over 1,000 km, far beyond one delta's reach at 54–102 km cells | Rim change stays near zero; response appears as lost volume or shifted shoreline instead |
| Routed deposition at the outlet | Already active in B | Removing routed sediment changes B's spill sequence, and A and C little |
| Sill as a route minimum | Explains B's migration | Spill rise tracks the second-lowest edge, not the deposit on the first |

## Applicable principles

- **Principle 6, account for material through every process:** silt is non-conservative, and its supply
  multi-counts internal transfers.
- **Principle 4, treat each node as an area:** a single rim node is compared with metres of clearance
  while its cell spans up to ±500 m.
- **Principle 3, model fast processes in aggregate:** a one-step spike has been read as a sustained rate.
- **Principle 9, verify before validating:** the Earth-analogue comparison was applied to a quantity that
  had not yet been verified.

## Recommendation

1. **Withdraw the registered sediment preflight and envelope.** The rates, two of its three
   connectors and the 0.4–8.8 mm/yr range do not support the screen. This is a blocked screen,
   not evidence that sediment is unimportant.
2. **Fix the depth-cap cycle first.** Detect it and export the surplus once, then trace where that
   water should go.
3. **Then fix the silt supply count.** Invariant: no body's step supply exceeds world external
   supply. Stopping the cycle may remove most of the excess. This changes model behaviour, so it
   needs the usual replay comparison.
4. **Then run a redesigned G8 diagnostic,** labelled diagnostic-only.

### Draft diagnostic, not yet registered or run

- **Windows:** A and B. Replace C or redefine its event on a single body, because it has no single
  connector. That leaves fewer than the three windows G8 asks for until a replacement is found.
- **Arms:** production; lake silt off; routed lake sediment off; both off. A one-coefficient
  sediment-off control is not one, because routed sediment would remain.
- **Metric:** the existing governing 10% rule, plus the spill-level sequence and the edge chosen
  each step.
- **Accounting:** water ledger residual at or below 0.1% of accepted supply; the routed-sediment
  mass budget closes; silt declared as its own source.
- **Uncertainty:** report the sub-grid relief range beside every clearance, and report node
  density and 10 kyr sidecar sensitivity. A and B share a world, so they are not independent.
- **Not usable:** Earth-analogue calibration of a sill rate, because no transferable source was
  found.

## Remaining uncertainty

- Only five steps per window were measured, on two worlds.
- The per-step supply multiple was measured only for the A and B components; C has the world
  bound only.
- Where final-export water ends up, and whether the cycle affects ocean delivery, is untraced.
- Measurements came from observational wrappers and the hierarchy's reported rims. No arm tested
  what happens with sediment removed.
- Source coverage for outlet aggradation came from a short search and is not a systematic review.

## Related finding outside this issue's scope

Zero-depth rims are real and accounted for. In window A's first step, `advance_lake_ice` reduces
world liquid from 4.52 M to 1.82 M km³, including 2.48 M km³ on 2,356 cells that end at exactly
zero. The ice budget books 2,697,367 km³ of sublimation against 2,697,020 km³ of liquid lost, so the water
is not deleted. It is applied per cell after the solve, so it can break a lake's level surface,
which could contribute to the mixed levels seen in A and C. That is a hypothesis for a separate
G8 mechanism: it has not been tested against the flicker.
