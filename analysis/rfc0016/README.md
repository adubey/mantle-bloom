# RFC 0016 diagnosis

## G0 preflight

Status: in progress; checkpoint preflight is complete, but G0 river-source selection and the
control ledger remain open.

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

For the checkpoint sensitivity, drainage area is defined as the sum of node accounting areas
upstream of each reach, with directed cycles collapsed before catchment accumulation. The model
has no long-term mean annual discharge per reach, so the Jobson comparison uses the explicit
`Qa=Q` sensitivity assumption; it does not claim that the current 100 kyr interval mean is an
annual climatology. Raymond's equation needs a registered rule for a dry reach that still holds
river water: `Q=0` gives zero velocity and infinite residence time. No velocity range or
parameterization is frozen, and no river inventory should be interpreted as accepted until the
low-flow rule and paired inventory sensitivity are recorded.

### Checkpoint and accounting readiness

The three registered #680 checkpoint windows, hashes, and the source checkpoint recreation
script are documented in the #680 preregistration. The private B2 save index contains the two
registered source worlds, but not the exact A/B/C checkpoint files. The preregistration says
the checkpoint bytes were produced at `5d666ddf8b682a72db8c2189ebbd45a8a861bb20`; later
production comparisons loaded those exact bytes under `f4438e1c85042de789bd464ce5d2f58c8c0717eb`.

I regenerated all three checkpoints from their registered source worlds at the preregistered
producer commit `5d666dd`, then replayed each four-step window under production commit `f4438e1`.
The source-world SHA-256 hashes match #680. The regenerated checkpoint hashes differ from the
registered file hashes, but all three registered target events reproduce under the pinned
production source:

| Window | Registered checkpoint SHA-256 | Recreated checkpoint SHA-256 | Reproduced target ΔV (km³) | Max normalized lake-ledger residual | Peak RSS (MB) |
|---|---|---|---:|---:|---:|
| A | `e4b0b4f09b87584aa7f77045f0f079f2309de0132049915a9aeb5589add27340` | `9882568e681a0f9f207e3302453b3a92f18865d5640d7baa8d9d8eb09ab4633e` | +495,623.986 | 5.82e-13 | 7,465 |
| B | `86fb5b7d2044846ac751906ff996184c7abe5454c93b39033bd82338e8749544` | `759d83107ae87d34d825b4836addc922d414089774f90757e0838899f0d56daed` | +131,625.097 | 1.33e-13 | 7,725 |
| C | `041b7430a6cf8867f2e2c76f18c3328e440add585ee4a2e4b3055e7766c2217e` | `38cc6ff2054cde2db3a65ac0a3e08a09e854175735512b9c8c40ba7a76c6190c` | −261,586.996 | 1.38e-12 | 7,196 |

The v5 checkpoint is a ZIP archive whose manifest entry gets a current-time ZIP timestamp
during serialization. The regenerated files therefore have different byte hashes even when
the checkpoint state reproduces the registered event and ledger behavior. The original
checkpoint files are not in the local/B2 save index; these replacement hashes and their
event validation are recorded here before any G0 treatment. Captures and complete file hashes
are kept in the private `~/mantle-bloom-captures/rfc0016/` directory.

RFC 0015 (#544) remains open and has not landed physical precipitation, snow, and glacier
accounting. Therefore the full rain-plus-snow G0 ledger cannot currently pass. RFC 0016 permits
a separately labelled rain-only interim ledger; it does not satisfy the full closure criterion.

### Checkpoint velocity sensitivity

The read-only `bin/debug/measure_river_velocity_candidates.py` diagnostic was run on all three
replacement checkpoints using the f443 hydrology cache. It counts land reaches outside ocean
and visible lake storage, accumulates catchment node area on the saved flow graph (collapsing
directed cycles), and evaluates the checkpoint's positive physical discharge. The Jobson
sensitivity sets its missing long-term mean discharge `Qa` equal to checkpoint discharge `Q`.
Full summaries are in `velocity-candidates-{a,b,c}.json`.

| Window | Eligible / positive-Q reaches | Raymond inventory (km³) | Jobson peak-tracer sensitivity (km³) | Lake-entering inventory, Raymond / Jobson (km³) |
|---|---:|---:|---:|---:|
| A | 45,246 / 42,630 | 251.435 | 508.782 | 5.933 / 24.550 |
| B | 45,454 / 42,971 | 239.335 | 435.628 | 5.556 / 6.031 |
| C | 39,158 / 34,849 | 463.567 | 799.080 | 4.468 / 4.680 |

The Jobson result changes total steady inventory by 1.7–2.0× relative to Raymond and the
lake-entering subset by 1.05–4.14×. These are initialization sensitivities `sum(Q L/U)`, not
checkpoint-to-checkpoint storage changes, control-ledger closure, or treatment results. The
calculation includes all positive discharges, so both empirical relationships are extrapolated
where model discharge lies outside their source data; several thousand eligible edges also have
zero/nonpositive flow. The RFC still needs a preregistered finite-velocity rule for a dry reach
with retained water.

Raymond provides a mean-flow velocity hydraulic-geometry fit (`U=0.19 Q^0.29`), based on
1,026 paired discharge/width/depth/velocity measurements; its velocity fit has `r²=0.49`.
The source provides a combined empirical fit but does not establish a validity range for
synthetic-world reaches. Jobson provides a second numerical estimate, but it predicts **peak
tracer velocity**, not mean advective velocity; using it to set the linear-reservoir residence
time is unsupported without a peak-to-mean mapping. The old Leopold–Maddock mean-velocity
exponent (`m=0.34`) alone supplies no transferable coefficient. These calculations expose
sensitivity but do not yet satisfy G0's requirement for two admissible, source-supported `U`
parameterizations and a defensible low-flow rule. No river state or G0 treatment has been
implemented; stop before the reservoir implementation gate until that choice is resolved.

### Next G0 actions

1. Resolve a second mean-advective-velocity parameterization (or revise the RFC's two-source
   gate), supported low-flow behavior, and the discharge domain over which each estimate is
   used. The current Jobson peak-tracer sensitivity is not an admissible mean-velocity law.
2. Once that gate is satisfied, implement persistent conservative reach state and per-reach
   transfer/closure readouts.
3. Use the validated A/B/C replacement checkpoints for the G0 control ledger, with their
   producer/production commits and replacement hashes recorded in every run report.
4. Report rain-only control closure while marking the full G0 ledger dependent on RFC 0015.

Checkpoint reconstruction and event validation replays have run at `5d666dd` and `f4438e1`.
Read-only velocity sensitivities have run on A/B/C. No G0 inventory or ledger replay, paired
treatment, or tests have run.
