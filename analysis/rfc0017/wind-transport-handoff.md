# Wind follow-up for RFC 0019 water delivery

This is a handoff for separate wind work, linked from the
[G3 moisture prototype plan](g3-coarse-transport-plan.md). RFC 0019 owns the column
transport and water ledger. Wind follow-ups own the inputs they provide, and can
run independently. A future wind change that affects delivery should be measured
against the frozen moisture control before it is used for G3 acceptance.

## Start with a comparison that attributes the effect

The [#526 comparison](../issue526/README.md#agreement-with-the-cfd) is useful
screening: one generated world had 62% land-biome agreement and a 0.52 relative
rainfall difference between modes. The older climate swap attributed most of
that precipitation difference to wind, with little effect from temperature
alone. The two worlds in the stepped comparison can drift apart, and the
[#923](https://github.com/adubey/mantle-bloom/issues/923) ERA5 result is a
surface-proxy comparison, not an error score for either model wind.

On identical frozen terrain start with A and cold, both of which have pinned
G3 controls. Capture the two wind fields at the same resolution and compute
the moisture difference with terrain, source rates, rain law, lagged cache,
temperatures and masks fixed. Compare magnitude alone (keep directions),
direction alone (keep speeds), then both. Preserve zeros in the speed-only
arm; exclude zero and weak-vector cells from angle summaries and report the
masked area. Run the full climate comparison once after this attribution,
because surface wind also drives currents and temperature, evaporation and
weathering. Use the official G3 controls for acceptance; these swaps are
diagnostics on tuning cases.

Report hemisphere-specific latitude bands of ocean and land speed, prevailing
angle and signed east/north components; onshore normal wind at coasts; and
coarse convergence with the same physical smoothing length across modes. In
the moisture run report signed water flux across coasts and latitude bands
(especially 60°), rain by distance inland, cold-land precipitation, and
windward/leeward rain. If terrain and temperature are fixed, differences in
source-tag arrival reveal routing changes rather than a changed source. Score
both the open-ocean belts and the downstream maps. Store compact zonal tables,
map comparisons, and the configuration and commit needed to repeat each arm.

Use the checked-in [ERA5 wind climatology](../issue526/era5_wind_climatology.json)
for a cheap Earth scale check on broad ocean latitude bands. Compare model
scalar speed to ERA5 mean scalar speed; compare direction to the direction of
ERA5 mean components when that vector is well defined. The mean speed and
the length of the mean vector are different quantities. An Earth latitude
band is useful for plausibility; an Earth pixel has no corresponding
coordinate on a generated planet.

## Diagnostic wind: direction first; surface roughness separately

**Where to work.** In [`climate.py`](../../backend/app/hydroclimate/climate.py),
`planetary_wind` builds bands and uses `CORIOLIS_DEFLECTION_GAIN`;
`compute_wind` adds a thermal gradient, terrain deflection and wake. The
existing diagnostic ocean mean speed is close to ERA5 after calibration, but
the trades are more meridional than the ERA5 mean flow. Increasing the common
deflection gain would also change wind belts outside the trades, and it would
change the CFD forcing target because that target reuses `planetary_wind`.

**Proposed first experiment.** On the saved snapshots, hold diagnostic speeds
fixed and rotate only the background trade direction toward the observed
ocean-band direction, with a smooth transition out of the trades and the
correct hemisphere signs. First do this as a diagnostic-only snapshot arm
so the CFD target stays fixed. Compare `planetary_wind` output before the
thermal/topography terms with the final wind: if deflection or the thermal term
already dominates the problematic coast, changing the planetary direction may do
little there. Check whether the change moves tropical convergence, redirects
onshore delivery, and shifts rain shadows; abandon it if the maps deteriorate
despite a better zonal wind angle. If promising, test a source-level change with
the same bounded direction adjustment before proposing to share it with CFD.

[#534](https://github.com/adubey/mantle-bloom/issues/534) owns surface roughness:
flat lowland can keep near-ocean speed. A separate roughness-only arm should
measure land/coast wind and moisture import, avoiding a second copy of terrain
drag or wake. If land cover sets roughness, use the lagged cover for the first
coupling test so current rain cannot immediately change current wind. Test whether
that slowdown belongs only in the 10 m field; column moisture can travel above the
surface layer. Do not tune direction and roughness together in the first screen.

## CFD wind: distinguish solver balance from the surface conversion

**Where to work.** In
[`atmosphere_cfd.py`](../../backend/app/hydroclimate/atmosphere_cfd.py),
`_target_wind` projects the same `planetary_wind` target onto the solver
grid. `step_atmosphere_cfd` supplies drag, viscosity and terrain geometry
to `_atmosphere_substep_loop_kernel`. In `climate.surface_wind` the layer
velocity is interpolated back to the climate grid and passed through
unconverted. [#530](https://github.com/adubey/mantle-bloom/issues/530)
records 42–66% of target speed across ocean bands and owns its 10 m contract
and runtime criteria.

**Proposed first experiment.** Hold forcing and terrain fixed and measure
short spin-up/settling on the native CFD grid and after interpolation,
at two fluid resolutions. Record atmospheric day length, sample count,
initial state and signs of approach to steady winds; keep the tectonic world
frozen. Compare an aquaplanet fixture with one terrain snapshot to distinguish
the basic balance from coast and mountains. Inspect target relaxation,
drag, pressure gradients, Coriolis, numerical damping and deflection; report
momentum tendencies and their projection on the local velocity as separate
checks. Coriolis can rotate the flow without directly dissipating energy.
This isolates whether the low speeds appear on the solver grid or are made
worse by remapping. Then test one hypothesis at a time:

- **Surface conversion only:** keep the CFD trajectory and direction fixed;
  convert the diagnostic output using a justified relation tied to a defined
  layer quantity, or clearly label an empirical conversion. Check whether
  speed-only changes improve water delivery and obey #530's band and cost checks.
- **Forcing or drag balance:** if the layer wind itself is too weak or its
  direction misroutes moisture, change one term and measure the equilibrium,
  directional changes, stability and step cost. Inspect the full vector and
  pressure response; matching target speed alone cannot demonstrate correct
  gradients or coastal flow.

Keep the velocity used by currents, evaporation and weathering explicit when
testing a conversion. A CFD surface-speed correction is not automatically the
effective column velocity, just as a diagnostic roughness correction is not a
column drag rule. If routing remains wrong after magnitude calibration,
hand back measured regions and flux directions to the moisture workstream.
The next transport hypothesis can then test a direction correction or a
meridional exchange term against the frozen control instead of inflating the
same scalar.
