# RFC 0016 diagnosis

## G0 preflight

Status: in progress; no G0 replay or ledger result yet.

RFC 0016 is accepted in PR #724 and tracked by #723. G0 implementation is tracked by
[#725](https://github.com/adubey/mantle-bloom/issues/725). The G0 worktree is based on
`origin/main` at `929d70c4f6b1c01e2c9a1994c17ea86a2d76e272`; the registered RFC 0016 source
commit is `f4438e1c85042de789bd464ce5d2f58c8c0717eb`.

### Current model interfaces

- `hydrology.compute_hydrology` computes annualized source volumes and routed `discharge_m3_s`.
  Its `compute_river_speed` remains a stylized unitless erosion/rendering proxy and is not a
  candidate physical velocity.
- `HydrologyFields` is the per-step, node-aligned output. `erosion.apply_erosion` owns
  persistence and writes selected hydrology outputs back with `Plate.set_fields_on_plate`.
- `surface_fields.SURFACE_FIELDS` registers persistent fields with remap semantics. An
  extensive water-equivalent depth can therefore preserve area-integrated volume through the
  existing refinement, coarsening, and merge paths; the runtime inventory must still be
  converted to and reported in m³ using the node accounting areas.
- The routing graph and discharge are derived each step. A river-water inventory must retain
  state on eligible routed edges even if that edge falls below the rendering `is_river` cutoff.
- Existing `channel_depth` and `channel_width` are geomorphic incision/width outputs; neither
  describes active wetted geometry.

### Source-supported velocity candidates

1. Raymond et al. (2012) report a discharge hydraulic-geometry relation for velocity,
   `U = 0.19 Q^0.29`, with `U` in m/s and `Q` in m³/s. Their compilation contains 1,026 paired
   observations and the velocity fit has `r² = 0.49`; it is a broad empirical estimate, not a
   calibrated law for Mantle-Bloom's synthetic reaches. [Paper](https://doi.org/10.1215/21573689-1597669)
2. Jobson's USGS national regression estimates peak tracer velocity from drainage area, reach
   slope, mean annual discharge, and discharge at measurement. In SI units its predictor is
   `Vp = 0.094 + 0.0143 (Da')^0.919 (Qa')^-0.469 S^0.159 (Q/Da)`, where
   `Da' = Da^1.25 sqrt(g)/Qa` and `Qa' = Q/Qa`. It was developed from about 90 U.S. rivers and
   939 observations; the peak-velocity fit has `r² = 0.70` and RMS error `0.157 m/s`. Because
   peak tracer velocity is not mean advective velocity, this is a sensitivity candidate only
   unless G0 establishes a defensible mapping for residence time. [USGS report](https://pubs.usgs.gov/publication/wri964013)

The candidate comparison still needs an operational drainage-area and mean-annual-discharge
definition for Jobson's equation in this model. Those inputs are not yet registered. No
velocity range or parameterization is frozen, and no river inventory should be interpreted
as accepted until that definition and the paired inventory sensitivity are recorded.

### Checkpoint and accounting readiness

The three registered #680 checkpoint windows, hashes, and the source checkpoint recreation
script are documented in the #680 preregistration. The private B2 save index contains the two
registered source worlds, but not the exact A/B/C checkpoint files. The preregistration says
the checkpoint bytes were produced at `5d666ddf8b682a72db8c2189ebbd45a8a861bb20`; later
production comparisons loaded those exact bytes under `f4438e1c85042de789bd464ce5d2f58c8c0717eb`.

I regenerated A/B from the registered seed-910 source world at both commits. At `5d666dd`, the
checkpoint ages match, but the new file hashes are `9882568e681a0f9f207e3302453b3a92f18865d5640d7baa8d9d8eb09ab4633e`
and `759d83107ae87d34d825b4836addc922d414089774f90739ea00e6a9c9dff5b5`, rather than the
registered `e4b0b4f09b87584aa7f77045f0f079f2309de0132049915a9aeb5589add27340` and
`86fb5b7d2044846ac751906ff996184c7abe5454c93b39033bd82338e8749544`. At `f4438e1`, the hashes are
`a4f60e0d37d4afeea06713a3a78f760d6ef2fbe77f861e304660412bac199ded` and
`44025c97e70e47d870f976721b79f3e53ad089bd2f24757e0838899f0d56daed`. The source-world SHA-256
matches the preregistration in both runs. The v5 checkpoint is a ZIP archive whose manifest
entry receives a current-time ZIP timestamp during serialization, so file hashes can change
even when the serialized state is otherwise identical. The checkpoint hashes alone therefore
do not establish a state mismatch.

The regenerated files are preserved under the private capture directory, grouped by producer
commit. Before using them as documented replacement inputs, replay the 5d666 checkpoints for
the registered four-step A/B targets under production commit f443 and confirm the preregistered
event volumes/connectivity. Keep the original #680 hashes in the report and record the replacement
hashes and the serialization reason before any G0 treatment. The C checkpoint still needs
regeneration and the same check.

RFC 0015 (#544) remains open and has not landed physical precipitation, snow, and glacier
accounting. Therefore the full rain-plus-snow G0 ledger cannot currently pass. RFC 0016 permits
a separately labelled rain-only interim ledger; it does not satisfy the full closure criterion.

### Next G0 actions

1. Freeze the velocity equations, units, valid discharge ranges, and uncertainty comparison,
   including whether Jobson's peak-tracer estimate is admissible for `U`.
2. Implement persistent conservative reach state and per-reach transfer/closure readouts.
3. Regenerate checkpoints from the recorded producer source, then validate their registered
   target events under the pinned production source before using replacement hashes.
4. Report rain-only control closure while marking the full G0 ledger dependent on RFC 0015.

Checkpoint reconstruction replays have run at `5d666dd` and `f4438e1`. No G0 inventory or
ledger replay, paired treatment, or tests have run yet.
