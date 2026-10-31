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
   slope, mean annual discharge, and discharge at measurement. It was developed from about 90
   U.S. rivers and roughly 939 observations; the reported peak-velocity fit has `r² = 0.70` and
   RMS error `0.157 m/s`. Because peak tracer velocity is not mean advective velocity, this is
   a sensitivity candidate only unless G0 establishes a defensible mapping for residence time.
   [USGS report](https://pubs.usgs.gov/publication/wri964013)

The candidate comparison still needs an operational drainage-area and mean-annual-discharge
definition for Jobson's equation in this model. Those inputs are not yet registered. No
velocity range or parameterization is frozen, and no river inventory should be interpreted
as accepted until that definition and the paired inventory sensitivity are recorded.

### Checkpoint and accounting readiness

The three registered #680 checkpoint windows, hashes, and the source checkpoint recreation
script are documented in the #680 preregistration. The private B2 save index contains the two
registered source worlds, but not the exact A/B/C checkpoint files. They can be regenerated
from the source worlds at the registered source commit. Their hashes must be checked before any
replay, and the runs must use the registered source commit rather than this implementation
branch for baseline diagnosis.

RFC 0015 (#544) remains open and has not landed physical precipitation, snow, and glacier
accounting. Therefore the full rain-plus-snow G0 ledger cannot currently pass. RFC 0016 permits
a separately labelled rain-only interim ledger; it does not satisfy the full closure criterion.

### Next G0 actions

1. Freeze the velocity equations, units, valid discharge ranges, and uncertainty comparison,
   including whether Jobson's peak-tracer estimate is admissible for `U`.
2. Implement persistent conservative reach state and per-reach transfer/closure readouts.
3. Regenerate and hash-check the registered checkpoints before a single baseline run.
4. Report rain-only control closure while marking the full G0 ledger dependent on RFC 0015.

No replays or tests have been run from this worktree.
