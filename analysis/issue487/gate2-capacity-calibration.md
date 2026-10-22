# Gate 2 strength/capacity calibration audit

## Result

The capacity model is now **scale-calibrated conditionally**, but there is no defensible
single Earth capacity field to assign to every continental node. The matched replay has
candidate load measurements, not observations of which nodes or events ruptured; using the
old algorithm's cuts as yield labels would make the calibration circular. The calibration
therefore checks independently published strength scales against explicit representative
columns and retains the uncertainty range.

The reproducible audit is [`calibrate_gate2_capacity.py`](../../bin/debug/calibrate_gate2_capacity.py)
and its compact output is [`gate2-capacity-calibration.json`](gate2-capacity-calibration.json).
It evaluates the existing 1,728-case continental scenario grid at five total thicknesses,
keeping the model's 35 km crust fixed and assigning the remaining thickness to the mantle lid.
Every percentile below summarizes an equally weighted *designed scenario grid*, not a
probability distribution.

| Total lithosphere | Scenario anchor | Factorial p10 / median / p90 | Grid fraction in published EARS 3–8 TN/m range |
|---:|---:|---:|---:|
| 70 km | 0.91 TN/m | 0.32 / 0.82 / 2.41 TN/m | 5.0% |
| 90 km | 4.09 TN/m | 0.81 / 2.88 / 8.07 TN/m | 37.8% |
| 100 km | **6.28 TN/m** | **1.49 / 4.88 / 12.07 TN/m** | **46.4%** |
| 135 km (current reference thickness) | 15.82 TN/m | 6.17 / 13.98 / 29.62 TN/m | 19.0% |
| 150 km | **20.86 TN/m** | **8.62 / 18.95 / 38.89 TN/m** | 8.0% |

The anchor uses `K0=0.8`, a depth-varying freshwater hydrostatic pressure `P_f=1000 kg/m3 * g * z`,
1330°C lithosphere-base temperature, strain rate `1e-15 s-1`, 20 MPa cohesion, friction coefficient
`mu=0.6` (`phi=atan(mu)`), plagioclase lower crust and dry olivine mantle. The friction coefficient
is the lower end of Byerlee's measured rock-friction range; the upper 0.85 endpoint is included in
the scenario grid. Hydrostatic pore pressure follows the water-column pressure law; unlike the
previous fixed 0.42 fraction, its ratio to lithostatic stress varies with depth and is lower than
0.42 for rock densities used here. K0, cohesion, temperature, strain rate and rheology remain
explicit modeling assumptions. Only the thickness cases and external comparison scale are
literature-motivated. The 100 km reference column's 6.28 TN/m falls within Stamps et al.'s 3–8 TN/m
integrated strength range for the East African Rift's 100 km lithosphere. This is a useful order-of-magnitude check, not a formal equivalence:
Stamps et al. use regional rheology and stress averaging, while this diagnostic measures
additional normal-fault differential capacity from an assumed `K0` state. The factorial median
also falls inside that range, but only 46.4% of the designed grid does. That spread shows the
benchmark does not identify temperature, fluids, rheology, strain rate or `K0` individually.

The 150 km column is a stronger-craton *scenario*, not a calibrated craton value. A 150 km
Tanzania/Congo craton and 70–90 km rift lithosphere are thickness choices in one published
geodynamic model; they motivate distinct columns but do not provide local strength observations.
A separate 2025 rifting model reports peak driving forces of about 22 TN/m for an isotropic
case and at most 10 TN/m with strong inherited anisotropy. Those are driving forces in that
model, not capacity values, so they are only broad force-scale context and are not fitted here.

The 74 connected-source replay candidates have an event-band median-of-medians load of
1.11 TN/m and median event-band p95 of 2.33 TN/m. Those values do not identify capacity:
the events were selected by the legacy algorithm, and there are no independent rupture labels.
The young-oceanic comparison remains uncalibrated because active divergence is only a proxy for
oceanic crust age; `divergent_age_myr` is a boundary-activity clock, not seafloor age. The young
proxy stratum is small (537 cut-band samples).

Anchor capacity changed by at most 1.19% when the integration step doubled from 0.125 to
0.25 km across the tested thicknesses. The 70 km case has the largest depth-grid difference;
the 100 and 150 km anchors differ by 0.43% and 0.12% respectively.

## Calibration choice and limits

Use these as conditional reference anchors for future sensitivity and replay work:

- **Mobile/rift-like continental column:** 100 km total thickness, 6.28 TN/m anchor, and
  1.49–12.07 TN/m p10–p90 scenario bracket.
- **Thick cratonic column:** 150 km total thickness, 20.86 TN/m anchor, and
  8.62–38.89 TN/m p10–p90 scenario bracket. The numeric value is a modeled scenario, not an
  externally calibrated target.
- **Thin rift column:** 70–90 km totals give anchors of 0.91–4.09 TN/m; retain these as a
  separate sensitivity range rather than forcing one rift-wide threshold.

These choices establish a defensible scale range for continued diagnosis. They do **not**
select the production capacity model. The matched replay now maps local columns using per-node
`Hc`, `Hm`, crust type and available proxy fields. A one-factor sensitivity sweep on those exact
node loads shows that the exceedance fraction depends most on `K0`, mantle-water/rheology and pore
pressure assumptions; see [`gate2-capacity-sensitivity.md`](gate2-capacity-sensitivity.md). The
world does not supply depth-resolved temperature, layer composition, pore pressure, mantle water,
local strain rate or inherited-fabric orientation. The 35 km crust held fixed in this audit also
differs from some regional rift profiles. Do not use one anchor for all continental nodes, tune
to the old cut fraction, or close Gate 2 on this result alone. Independently constrain the missing
local profile and rheology state, then complete broader shell load/boundary validation.

## Reproduction

```bash
backend/.venv/bin/python bin/debug/calibrate_gate2_capacity.py \\
  --dz-km 0.125 --output analysis/issue487/gate2-capacity-calibration.json
```

Scale references:

- Byerlee (1978), laboratory friction of rocks, supporting the `mu=0.6-0.85` friction bracket:
  [USGS publication record](https://pubs.usgs.gov/publication/70012542).
- Stamps et al. (2014), East African Rift integrated strength and tectonic stress estimates:
  [Journal of Geophysical Research: Solid Earth](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013JB010717).
- Rajaonarison et al. (2021), regional lithosphere thicknesses used for East African rift and
  craton domains: [Geophysical Research Letters](https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020GL090483).
- Duretz et al. (2025), anisotropy sensitivity of modeled rifting forces:
  [Geochemistry, Geophysics, Geosystems](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2025GC012409).
