# RFC 0019 G3 coarse lift arm

Registered before running the lift-aware comparison. This is one diagnostic arm on tuning
worlds, not a new G3 acceptance gate or a production parameter selection.

## Inputs and control

- Tuning saves only: A (`seed910211954-221.6myr.bloom`, SHA-256
  `1f167a8d9e0edf06c3c52b23598dc2f776adc2927f843beb992b5193ad30fd67`) and cold
  (`mantle-bloom-seed317739385-200800000y.bloom`, SHA-256
  `c7e8c4094376b2b57158e36b1066cf529eb475ac0bbf54a6b888b391f294d344`). No held-out save
  is loaded.
- 90×180, finite-volume transport, diagnostic and CFD wind modes.
- The control is the corrected uniform-sink candidate already measured on the same save,
  wind, grid, and solver. Keep its source conversion, mask, effective depth (2 km), source
  time (9.1 days), background rainout (9.1 days), wind multiplier (1), and land-proxy
  multiplier (1) unchanged. The only treatment is added terrain-lift rainout.

## Fixed lift law

1. Smooth elevation on the sphere with a Gaussian of `σ = 250 km`, using a 3σ neighborhood.
   This is a fixed physical scale near one equatorial cell width at density 1.
2. Calculate positive terrain ascent from the same wind used by transport:
   `w_up = max(0, u * ∂h_smooth/∂x + v * ∂h_smooth/∂y)`.
3. Add one bounded local sink:
   `k_lift = min(w_up / 2250 m, 1 / 5 days)`.
   The 2250 m scale is the midpoint of the preregistered 1.5–3.0 km effective-depth range;
   the 5-day cap is the fastest preregistered rainout sensitivity. Add `k_lift` to the
   uniform background sink and debit it through the same stationary water solve.

The law is a deliberately coarse, local lift proxy. It omits cloud-water advection, fallout
delay, and downslope evaporation; this screen cannot validate those omitted processes.

## Predictions and stop

- Terrain-ascent cells receive more rain than their paired descent cells, strengthening the
  area-weighted windward/leeward contrast.
- Smoothing at a fixed physical length broadens the response and should reduce the share of
  error attributable to unresolved fine-scale peaks after 90×180 aggregation.
- Track cold high-latitude land rain and dry-biome share separately; coastal lift may remove
  water before it reaches cold interiors or dry basins.
- Record closure, nonnegativity, source totals, coarse-scale raw and best-fit RMS, rain
  quantiles, hemisphere/surface bands, land rain, Köppen B share, high-latitude land rain,
  windward/leeward pairs, runtime, and RSS for each arm.
- Run no parameter sweep. If this single arm does not improve the predicted terrain contrast
  and coarse-scale pattern on both tuning worlds without a clear cold-land or dry-interior
  regression, stop the local-lift path and revisit transport or the stationary formulation.

These directional checks do not replace or change the approved held-out G3 limits. Closure is
verification; product realism and Earth comparison remain separate validation questions.
