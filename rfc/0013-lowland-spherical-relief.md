# Lowland relief from a spherical residual field

- **Status:** Draft
- **Date:** 2026-10-19
- **Related issues:** none yet
- **Related RFCs:** [0009](0009-craton-mosaic.md), [0012](0012-physics-first-simulation.md)

## Problem

Newly generated continental elevation is dominated by broad smooth regions. In four
prototype seeds (7, 42, 2024, 611937962), the fraction of sampled five-point land transects
that changed monotonically in either direction inside a 500 m elevation band was about
78% for the baseline and 69–73% after adding a correlated relief field. A preliminary
sample of ETOPO 2022 gave about 62–64% under the same nominal calculation. These figures
are exploratory: transects overlap, the Earth sample was not area-weighted, and the
prototype and comparison script are not checked into the repository. They suggest a
possible mismatch in small-scale relief, but do not establish its cause or global size.

The affected behavior is the initial topographic texture of continents. This is a model
boundary issue under Design Principle 1 (state boundaries should represent the modeled
system) and Principle 9 (validation should make uncertainty visible): a generic random
field could improve a metric while making terrain less physically interpretable. The
invariant to evaluate is the distribution of within-band transect slopes and spatial
spectra, with ocean and high-elevation terrain excluded from the proposed lowland residual.

## Proposed solutions

Real landscapes contain landforms and relief at a range of spatial scales. Landscape
elevation spectra are not generally featureless random spectra; characteristic scales and
departures from fractal scaling have been observed in measured landscapes
[[Perron et al., 2008]](#ref-perron2008). ETOPO 2022 provides a globally seamless,
validated topography and bathymetry reference suitable for a reproducible comparison,
subject to its resolution and source-data limits [[MacFerrin et al., 2025]](#ref-macferrin2025).
Neither source implies that adding isotropic noise is a physical model of landform
formation.

1. **Keep current relief.** Treat the preliminary difference as a measurement artifact or
   as a limitation better addressed by improving existing relief processes. This avoids
   adding a nonphysical mechanism, but may leave broad smooth lowlands if the difference is
   reproduced.
2. **Add a bounded correlated residual to lowland elevation.** Use a deterministic,
   world-seeded spherical field as unresolved topographic texture, tapering it near sea
   level and at high elevation. This is cheap and reproducible, but remains phenomenological
   and may distort drainage or regional relief.
3. **Condition relief on modeled processes.** Develop relief from plate interiors, rifts,
   faults, erosion, or other represented processes. This is more interpretable but requires
   a validated process model and may not address the same spatial scales.

Recommendation: run the diagnosis gate before choosing an implementation. If it confirms a
robust, resolved lowland mismatch and alternatives cannot explain it, consider solution 2
as an explicitly phenomenological residual. Otherwise prefer solution 1 or 3. The residual
must not be described as tectonic or geological relief.

## Design

### Diagnosis gate

Add a reproducible offline diagnostic in `bin/debug/` that compares baseline generated
terrain, the proposed residual, and Earth terrain. It should:

- Evaluate the monotonic-transect metric in fixed 500 m bands using the same five-point
  windows and four directions for all datasets; report land-area-weighted results and
  uncertainty across spatial blocks, not just the pooled fraction.
- Use ETOPO 2022 15 arc-second data where available, with a documented coarsening procedure
  to match model resolution. Exclude ocean, use the same elevation bands, and report the
  selected geographic coverage and resampling method. Keep downloaded data and large arrays
  out of the repository; record its DOI and checksums/instructions.
- Compare elevation distributions and two-dimensional or radial power spectra by latitude
  and elevation band, so a single transect statistic cannot select the model.
- Run at least 20 fixed terrain seeds at two model node densities. The original four seeds
  (7, 42, 2024, 611937962) are a smoke set, not independent confirmation.

If the seed and block-bootstrap intervals overlap the Earth reference after matching
resolution and land-area weighting, or if spectra show no robust lowland deficit, stop and
retain the baseline. If a deficit persists across seeds, blocks, and resolutions, proceed
to compare the residual against process-conditioned alternatives and update this RFC with
the result before an implementation issue is opened.

### Candidate residual

The prototype is a signed, zero-centered sum of three `SphereNoise` octaves sampled on
unit-sphere positions. `SphereNoise` uses random directional sinusoids, not gradient noise.
At Earth radius, base frequency 40 and octave doubling correspond approximately to
1000 km, 500 km, and 250 km wavelengths along great circles. For each generated world,
the field is seeded once from the world seed and remains spatially continuous and stable
over the simulation; it is not a new independent field per node or turn.

The prototype parameters are calibration choices, not literature-derived values:

| Parameter | Prototype value | Units and role |
|---|---:|---|
| Amplitude scale | 225 | m; multiplies the normalized signed field before tapers (not the resulting RMS) |
| Base frequency | 40 | spherical cycles; approximately 1000 km wavelength |
| Octaves | 3 | wavelengths approximately 1000, 500, 250 km |
| Coast taper | 0–400 | m elevation; suppresses the residual near the land threshold |
| Highland taper | 1,400–5,000 | m elevation; full strength below 1,400 m, fades to zero by 5,000 m |

Apply the coast and highland factors independently and clamp their product to [0, 1]. The
candidate should use `lithosphere.back_elevation_gain_fields` to translate supported
elevation changes into crustal thickness changes, preserving the current isostatic
bookkeeping. A final bound must preserve the land mask and existing elevation limits; the
lower bound must not turn a positive land cell into ocean. The implementation should
confirm these constraints against `initial_terrain.py` before settling the exact clamp.

The field is generated once during initial terrain construction, so its result should be
independent of climate timestep or turn count. No new persistent state is needed if the
field is applied only during generation. Expected runtime is linear in node count and a few
sinusoids per node; benchmark on the largest supported mesh before acceptance. Save
compatibility should be unchanged because generated saves already store elevation and
crustal thickness, not the procedural field.

Tune amplitude and tapers only after the diagnosis gate. Sweep amplitude scale 100–300 m in
50 m increments, coast taper width 200–600 m, and highland full-strength/fade elevations
of 1,000–2,000 m / 4,000–6,000 m. These are sensitivity ranges for calibration, not
geophysical bounds. Compare all candidates on the fixed seed set and retain the smallest
amplitude that meets the accepted metric without degrading the drainage and spectrum
checks.

## Detailed design (provisional touchpoints)

- `backend/app/terrain_noise.py`: existing `SphereNoise`-based relief utilities; likely
  location for a new stateless residual-field helper if the diagnosis selects solution 2.
- `backend/app/initial_terrain.py`: apply the residual during initial terrain construction
  and preserve the land mask.
- `backend/app/lithosphere.py`: existing `back_elevation_gain_fields` for isostatic support;
  avoid duplicating its bookkeeping.
- `backend/unit_tests/test_terrain_noise.py` and a new or existing `initial_terrain` test:
  determinism, zero ocean/highland taper, bounds, preserved land mask, and isostatic support.
- `bin/debug/`: new reproducible terrain-vs-ETOPO diagnostic, if the gate is approved.
- `docs/simulation-model/`: document the residual as a generation-time phenomenological
  field, if implemented.

## Assumptions and uncertainty

- The prototype's four-seed results are exploratory and do not establish an Earth-model
  mismatch. The Earth comparison was not area-weighted and does not account for spatial
  autocorrelation.
- ETOPO grid spacing is not equivalent to model resolution. Resampling, elevation-band
  membership, land masks, and polar coverage can change the metric.
- Correlated sinusoids are a residual parameterization, not a process model. Their apparent
  wavelength and RMS do not map directly to a particular geological mechanism.
- The 225 m RMS and taper elevations came from visual/metric experimentation and need
  sensitivity analysis; they are not constrained by Perron et al. or ETOPO.
- Adding signed relief may affect flow routing, lake formation, erosion, and coastlines even
  if the starting land mask is held fixed. Five-turn climate-only replays with geology
  disabled are useful smoke comparisons, but not long-term validation.

## Acceptance criteria

1. **Diagnosis gate:** the reproducible, area-weighted, block-aware comparison shows a
   lowland mismatch beyond uncertainty at matched resolution across the fixed seeds and
   both node densities. If not, do not implement the residual.
2. The selected solution improves the transect metric toward the Earth reference without
   worsening banded spectra, elevation distribution, or basin/lake coverage beyond agreed
   tolerances; record thresholds and results in the implementation issue before coding.
3. If solution 2 is selected, unit tests establish deterministic output, zero ocean and
   highland residual, valid bounds, unchanged land mask, and isostatically backed elevation
   change.
4. Replays of the four exploratory seeds and the expanded fixed set report initial metrics
   and five-turn climate-only maps with geology disabled. Confirm no timestep dependence
   and no regression in saved-world loading.
5. Runtime increase is below 2% for initial terrain generation on the largest supported
   mesh, or document and accept a different budget before implementation.

## Falsification criteria

- If area weighting, block uncertainty, resolution matching, or spatial spectra remove the
  apparent gap, the premise for a generic residual is falsified; retain current relief or
  investigate the metric.
- If the residual improves monotonicity but yields implausible spectra, basin topology,
  drainage, or strong seed sensitivity, reject this parameterization and revisit
  process-conditioned relief.
- If the residual cannot preserve land/ocean classification while maintaining its target
  amplitude, reject solution 2 rather than relaxing the land-mask constraint.
- If generation-time terrain differs materially with timestep, turn count, or traversal
  order, treat that as an implementation defect and do not accept the design.

## References

- <a id="ref-perron2008"></a>

Perron, J. T., Kirchner, J. W., & Dietrich, W. E. (2008). Spectral signatures of
characteristic spatial scales and nonfractal structure in landscapes. *Journal of
Geophysical Research: Earth Surface*, 113(F4). [doi:10.1029/2007JF000866](https://doi.org/10.1029/2007JF000866)

- <a id="ref-macferrin2025"></a>

MacFerrin, M., Amante, C., Carignan, K., Love, M., & Lim, E. (2025). The Earth Topography
2022 (ETOPO 2022) global DEM dataset. *Earth System Science Data*, 17, 1835–1849.
[doi:10.5194/essd-17-1835-2025](https://doi.org/10.5194/essd-17-1835-2025)
