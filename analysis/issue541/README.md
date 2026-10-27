# RFC 0014 G5: bounded lake-ice cap diagnosis (#541)

## Updated G5 result: 2 m non-permafrost cap

The G5 diagnostic was rerun with the owner-directed hard 2 m physical-ice ceiling and the
three preregistered #538 diagnosis saves. The diagnostic retains the existing below-cap
seasonal growth/loss proxy and runs all nine sublimation (×0.25/×1/×4) by conduction
(×0.5/×1/×2) settings for 10 Myr. It holds climate, bathymetry, area and water supply fixed.
It assumes non-permafrost lakes; Mantle Bloom has no permafrost or ground-temperature state.

The cap is physical ice thickness `h_i`; its water-equivalent ledger depth is
`H_w = h_i × 917/1000`, bounded by the lake's available water depth. Thus the effective
physical limit is `h_i ≤ min(2 m, lake_depth × 1000/917)`. A cell is bedfast when its cap's
water equivalent reaches lake depth. Lake depths ≤2 m are reported as shallow; the 1% gate
applies only to bedfast area over all lake area deeper than 2 m.

### Baseline setting (1× sublimation, 1× conduction)

| Seed | Physical ice p50 / p90 / p99 / max (m) | Shallow bedfast / wet area | Deep bedfast / deep-water area | Cap volume at 10 Myr (km³) | Sublimation (km³/yr) |
|---:|---:|---:|---:|---:|---:|
| 910211954 | 2.00 / 2.00 / 2.00 / 2.00 | 2.570% | 0% | 20,107.81 | 3,335.50 |
| 635056485 | 2.00 / 2.00 / 2.00 / 2.00 | 0.537% | 0% | 3,331.73 | 734.03 |
| 317739385 | 1.71 / 2.00 / 2.00 / 2.00 | 11.231% | 0% | 7,355.75 | 1,213.42 |

Thickness percentiles include zero-thickness wet cells. For all three saves and all nine
sensitivity settings, physical thickness stays at or below 2 m, deep-water bedfast area is
0%, and cap volume reaches its plateau by 0.1 Myr and remains there through 10 Myr. The
revised deep-water bedfast criterion passes in this fixed-climate G5 diagnostic. Shallow
bedfast area can exceed 1% of total wet area and is reported separately by the owner's
criterion.

### Sensitivity range across nine settings

| Seed | Shallow bedfast / wet area (%) | Deep bedfast / deep area (%) | Cap volume range (km³) |
|---:|---:|---:|---:|
| 910211954 | 1.155–2.679 | 0 | 6,552.77–20,880.56 |
| 635056485 | 0.177–0.537 | 0 | 1,909.81–3,336.73 |
| 317739385 | 1.633–11.271 | 0 | 1,492.23–8,291.71 |

Compact outputs are `rerun-2m-seed910211954.json`, `rerun-2m-seed635056485.json`, and
`rerun-2m-seed317739385.json`. Each preserves the original save checksum and per-run volume
trajectory. Reproduce one save at a time using the command in the diagnostic script's usage
docstring, directing output to `analysis/issue541/rerun-2m-seed<seed>.json`.

This is the RFC 0014 G5 offline prototype, not a production-code replay. It verifies the
2 m physical bound, available-water-equivalent bound, depth-stratified bedfast metric and
fixed-forcing plateau; it does not validate dynamic lake routing or a full production water
ledger.

## Prior result: unbounded cap

The starting D2(c) formulation does **not** pass G5. The baseline sublimation and conduction
settings reach a volume plateau on all three preregistered diagnosis saves, but more than 1% of
wet area freezes to the bed on two saves. Across the requested sensitivity ranges, no single
combination keeps all three saves below 1%: seed 317739385 remains at 1.63% or more. The cap
thickness distribution also changes too much across the knob range to adopt these values without
revising the law.

The saved worlds were those preregistered in #538. Their checksums match that record:

| Seed | Saved age (Myr) | SHA-256 | Freshwater lake nodes | Wet lake area (km²) |
|---:|---:|---|---:|---:|
| 910211954 | 220.4 | `bd4874092efc805e8b174e115fe567c82d35fc2ee606686938b5b6eb96b060d2` | 3,838 | 14,877,406 |
| 635056485 | 452.8 | `ec27d5077bdf43b3a6f0339b5726fa4d110da9eb1f143f96a3004d79f309b045` | 507 | 1,985,204 |
| 317739385 | 200.8 | `6c806f913e4c79bb351b5eb5a143abcfc273df39b3bdf4ce9823fec7d27fb0cd` | 2,270 | 6,355,998 |

The first file's name and #538 record label it 220.2 Myr; the persisted `World.elapsed_years`
in the checksum-matched save is 220.4 Myr. The simulation's persisted age is reported here.
The saves themselves are not committed; the per-save compact JSON is the reproducible output.

### Baseline at 10 Myr

Thickness below is physical ice thickness (water-equivalent `H` converted by 1000/917). The
reported percentiles include zero-thickness wet cells. Frozen-to-bed share and sublimation are
area-weighted over all wet freshwater lake nodes.

| Seed | Thickness p50 / p90 / p99 / max (m) | Frozen to bed (%) | Cap volume (km³) | Sublimation (km³/yr) |
|---:|---:|---:|---:|---:|
| 910211954 | 4.50 / 6.62 / 10.56 / 11.23 | 5.088 | 50,884.51 | 3,335.50 |
| 635056485 | 9.50 / 10.43 / 11.53 / 13.01 | 0.537 | 15,260.98 | 734.03 |
| 317739385 | 1.71 / 3.04 / 4.82 / 6.31 | 18.450 | 8,988.37 | 1,213.42 |

Volume is already unchanged at the 0.1 Myr checkpoint in all runs (and remains unchanged at
each later checkpoint through 10 Myr). Thus the plateau metric passes under fixed climate and
depth. Frozen-to-bed fails the RFC's 1% criterion for two of three baseline saves. The 4.50 m
and 9.50 m median thicknesses on seeds 910211954 and 635056485 exceed the published Lake
Untersee perennial range of 2.2–3.8 m (summer mean 3.0 m); the 1.71 m median on seed 317739385
is near the observed seasonal tundra-lake average seasonal maximum of 1.25 m. The first two
saves' upper tails are also above the perennial range, while 18.45% of seed 317739385's wet
area freezes solid. These comparisons are scale checks, not a one-to-one validation: the model
diagnostic is a static-climate multi-year cap equilibrium, while the seasonal observation is a
within-year maximum.

The 0.25x/1x/4x sublimation by 0.5x/1x/2x conduction runs are in each seed JSON. Across those
nine scenarios per seed, frozen-to-bed area ranges are 1.16–10.27% (910211954), 0.18–0.74%
(635056485), and 1.63–21.50% (317739385). Final cap volume spans 6,553–379,053 km³,
1,910–121,989 km³, and 1,492–57,248 km³ respectively. No RFC 0014 held-out seed was used.

## Method and interpretation

`bin/debug/gate_lake_ice_cap.py` restricts inputs to the three exact SHA-256 values above. It
loads one save, samples the existing climate fields and lake depths on live nodes, excludes
connected ocean, and holds temperature, seasonal amplitude, lake depth and area fixed. For each
wet cell it solves the Stefan conduction ODE against constant melt and sublimation losses, in
closed form, at 0, 0.1, 0.5, 1, 2.5, 5, 7.5 and 10 Myr. This avoids a long integration whose
answer would depend on timestep. The baseline conduction coefficient is derived from bare-ice
Stefan growth, `h = sqrt(2 k FDD / (rho_i L_f))`: using the code's `k=2.1 W m⁻¹ K⁻¹`,
`rho_i=917 kg m⁻³` and `L_f=334 kJ kg⁻¹` gives 3.44 cm per √(°C day). Sensitivity scales `k`
by 0.5 and 2.

For melt, this screen reuses the existing glacier law: 400 m/Myr at 20 annual-mean positive °C,
capped at 3x. Sublimation starts at 0.6 m water equivalent per ice-covered year, multiplied
by each cell's `1 − f_free`, and is tested at 0.15, 0.6 and 2.4 m/yr. The 0.6 m/yr source is an
observed *summer ablation* scale at Lake Untersee; it is not a direct measurement of annual
sublimation. Published observations say the roughly 0.6 m summer loss is compensated by winter
refreezing at the lake bottom, a process this D2(c) prototype omits. This distinction and the
large global fluxes make `s_ice` a weak starting assumption, not a literature-constrained value.

The result does not prove that the physical Stefan law itself is wrong. It shows that combining
fixed air-temperature degree forcing, the glacier melt proxy, a seasonal-fraction sublimation sink, and
a no-basal-heat bed-freeze rule cannot meet the current area threshold robustly across the
diagnosis set. Before implementation, revise D2(c) to resolve the bed-freeze behavior and
reassess whether the Untersee ablation estimate is an appropriate sublimation parameter. Do not
fit against the RFC 0014 held-out seeds.

## Sources

- Leppäranta, M. (1993), “A review of analytical models of sea-ice growth,” *Atmosphere-Ocean*
  31(1), 123–138, [doi:10.1080/07055900.1993.9649465](https://doi.org/10.1080/07055900.1993.9649465).
  The numerical coefficient above is derived from the Stefan equation and the model's ice
  conductivity, density and latent heat; it is the bare-ice value without snow insulation.
- Hirose, T. et al. (2008), “Bottomfast Ice Mapping and the Measurement of Ice Thickness on
  Tundra Lakes Using C-Band Synthetic Aperture Radar Remote Sensing,” *JAWRA* 44(2), 285–292,
  [doi:10.1111/j.1752-1688.2007.00161.x](https://doi.org/10.1111/j.1752-1688.2007.00161.x).
  The study reports 1.25 m as the average maximum seasonal ice thickness in its study region.
- Steel, M. et al. (2015), “Modeling circulation and seasonal fluctuations in perennially
  ice-covered and ice-walled Lake Untersee, Antarctica,” *Limnology and Oceanography* 60(4),
  1231–1245, [doi:10.1002/lno.10086](https://doi.org/10.1002/lno.10086). It reports 2.2–3.8 m
  observed perennial ice thickness, averaging 3.0 m in summer.
- Weisleitner, K. et al. (2019), “Source Environments of the Microbiome in Perennially
  Ice-Covered Lake Untersee, Antarctica,” *Frontiers in Microbiology* 10, 1019,
  [doi:10.3389/fmicb.2019.01019](https://doi.org/10.3389/fmicb.2019.01019). It describes ~60 cm
  summer ice ablation compensated by winter refreezing.

## Files

- `seed910211954.json`
- `seed635056485.json`
- `seed317739385.json`
