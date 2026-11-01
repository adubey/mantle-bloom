# RFC 0016 diagnosis

## G0 preflight

Status: **blocked; do not start G1.** The persistent river inventory and all reach transfers
close, and the control reproduces every registered event, but the rain-only land budget closes
only by naming 3.2–5.8% of A and B's input as stranded lake overflow, a defect in the control
(see [G0 control replay](#g0-control-replay)). The interim glacier account closes the full ledger
only by construction, and its balance is too large for its release to support G0's interpretation,
so G0 stays blocked until RFC 0015 (see [Interim glacier account](#interim-glacier-account-732)).

RFC 0016 is accepted in PR #724 and tracked by #723. G0 implementation is tracked by
[#725](https://github.com/adubey/mantle-bloom/issues/725). The G0 branch is based on
`origin/main` at `0f4aac8f`, which has no simulation change from the registered RFC 0016 source
commit `f4438e1c85042de789bd464ce5d2f58c8c0717eb`.

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

### G0 implementation

Each routed reach holds persistent extensive `river_water_m`, updated by the exact constant-inflow
linear reservoir. The central travel time uses Raymond et al.'s `U=0.19Q^0.29` with a one-year
maximum. A save with no river state starts from its saved `Qτ`; lake-covered nodes, sinks and flow
cycles hold no water and pass it on. A steady reach delivers exactly what `route_downstream` does,
source by source, which a unit test holds. Exported lake overflow runs a second reservoir pass,
so an outlet river holds its share, and released river water is its own lake-ledger source.
Eustasy counts the inventory as trapped water.

The first implementation changed the control, and no replay of it exists. It charged dry runoff
the full transit loss, ended water at every visible lake (including floor-draining ones), dropped
wet flow from discharge, and left outlet rivers with no inventory. Those were fixed before the
replays below.

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
annual climatology. Raymond's equation needs a rule for a dry reach that still holds river
water: `Q=0` gives zero velocity and infinite residence time.

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
| B | `86fb5b7d2044846ac751906ff996184c7abe5454c93b39033bd82338e8749544` | `759d83107ae87d34d825b4836addc922d414089774f90739ea00e6a9c9dff5b5` | +131,625.097 | 1.33e-13 | 7,725 |
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
exponent (`m=0.34`) alone supplies no transferable coefficient. Following the owner's direction
that G0 may use estimates, use Raymond as the provisional central estimate and Jobson as a faster
travel-time sensitivity, keeping Jobson's peak-tracer interpretation explicit. Cap residence
time at one year so dry reaches with stored water drain over a finite interval; initialize a
zero-flow reach at zero volume and retain that cap if it later dries. Treat the one-year cap as
an engineering regularization, not a literature result, and include 0.1- and 10-year caps in the
sensitivity report. This is sufficient to proceed with the conservative reach-state
implementation; the comparison remains an uncertainty bound rather than a calibrated river
model.

The window B hash above was garbled in earlier copies of this report and in #725's handoff
comment. The value here is `shasum` of the file.

## G0 control replay

Source commit `9e9c76bc`, on the replacement checkpoints in the table above, four 100 kyr steps
with the ice-age cycle off, as #680 registered them. Arms: Raymond with 0.1-, 1- and 10-year
residence caps, and Jobson with a 1-year cap. Whole-run values are in
[`g0-control-summary.json`](g0-control-summary.json). Replay JSON and per-reach arrays (the central
arm) are outside the repository, under `g0-9e9c76bc/` in the private captures; their `SHA256SUMS`
has SHA-256 `b73811ca9c65639bfd23171dfe1dc82a190419b0b2b3cb1ad7b3b811d4cc263d`.

```bash
backend/.venv/bin/python bin/debug/replay_rfc0016_river_water.py <checkpoint.mbworld> \
  --anchors analysis/rfc0016/target-anchor-window-<w>.json --myr 0.4 --step-years 100000 \
  --checkpoint-myr 0.1 --velocity <raymond|jobson> --residence-cap-years <cap> \
  [--reach-dir <captures>/reaches-<w>] --out <captures>/window-<w>-<velocity>-cap<cap>.json
```

**The control is unchanged.** Target volumes match the registered event, and global lake volume
stays within 1.6 km³ of the `f4438e1` run, out of about 2–4 million km³, in every arm.

| Window | Target V before→after event (km³) | ΔV (km³) | Registered ΔV (km³) | Median step, this / `f4438e1` (s) |
|---|---:|---:|---:|---:|
| A | 383,558.5 (registered)→879,182.5 | +495,624.0 | +495,624 | 7.7 / 7.6 |
| B | 115,361.0→246,986.1 | +131,625.1 | +131,625 | 7.7 / 7.5 |
| C | 264,786.7→3,199.7 | −261,587.0 | −261,587.0 (recreated); −261,592 (registered) | 10.0 / 9.8 |

Step times are indicative only, since other agents' runs overlapped some arms. Peak RSS reached
7.9 GB. A's anchor is dry at the first output, so its starting volume comes from #680.

**River inventory and transfers close.**

| Window | Inventory start→end, central (km³) | Jobson start→end (km³) | Largest reach identity error (m³) | River residual (share of input) | Largest inter-step inventory drift (m³) |
|---|---:|---:|---:|---:|---:|
| A | 229.8→241.5 | 459.6→464.2 | 46 | 9.9e-17 | 4.5e4 |
| B | 237.6→267.7 | 451.4→610.4 | 71 | 2.1e-16 | 3.3e6 |
| C | 448.0→462.3 | 741.2→771.9 | 39 | 1.9e-16 | 6.2e7 |

The reach error is `V₁ − V₀ − (IΔt − O)`, recomputed from the saved arrays, where a reach passes up
to about 1e17 m³ per step. Drift is the inventory change between one hydrology step and the next,
from tectonics and remapping; at most 0.014% of the inventory. The inventory is about 1e-7 of the
water entering rivers in a step, so it is a negligible term in the land budget.

The 0.1- and 10-year caps give the same inventory to the reported precision: only near-dry
reaches reach a cap. Jobson raises the inventory 1.7–2.3×, matching the checkpoint estimate. Neither
moves target-lake supply: new dry runoff to the target is identical across arms to 2 decimals (m³/yr
×1e9), and water released from river storage is about 1e3–1e4 m³/yr against up to 5e11. In window C
the target lake receives effectively no routed supply in these steps.

**The rain-only land ledger does not pass.** After separating the named defect, the residual is
tiny, but the defect is not a sink the RFC's ledger allows.

| Window | Interval input (km³) | Unexplained residual | Stranded overflow, whole run | Stranded overflow, worst step | Lake residual | Scaled melt share |
|---|---:|---:|---:|---:|---:|---:|
| A | 2.54e10 | 2.3e-8 | **5.77%** | 6.10% | 5.6e-13 | 1.6e-5 |
| B | 2.62e10 | 2.0e-5 | **3.16%** | 5.25% | 1.3e-13 | 1.4e-5 |
| C | 3.17e10 | 8.8e-11 | **0.11%** | 0.34% | 1.4e-12 | 1.3e-5 |

**Cause.** The lake solver debits overflow when it leaves the lake system, then
`compute_hydrology` routes it along `flow_target`. Water that ends anywhere but the ocean is credited
to nothing. `bin/debug/measure_rfc0016_unclaimed_settlement.py` traced window A's first step:
all of the 3.33e17 m³ residual is overflow settling on land, against 5.88e17 m³ exported, and the
runoff pass loses none of its supply. 2.60e17 m³ of it lands on one dry basin node at −4,721 m,
on the cycle `8606→8605→8604→8648→8690→8650→8606`, which closes through leaf sink 8690's own spill
edge. `route_downstream` keeps a cycle's remainder like a sink, and nothing reads it. This
predates G0: `f4438e1` routes the same water and drops it the same way, because its lake ledger
closes before this routing. [#696's report](../issue696/README.md#where-exported-water-goes)
measured the same leak as 22.5% of external supply, and this ledger measures it against total
input.

This is design principle 6's failure: a volume leaves one reservoir and reaches none. Rerouting
it sends supply to lakes that don't receive it today and may move the registered events, so it is
not a G0 control change. [#731](https://github.com/adubey/mantle-bloom/issues/731) owns the fix and
its guards, and [#732](https://github.com/adubey/mantle-bloom/issues/732) proposes an interim
frozen-water account for the full ledger, subject to an RFC 0016 revision.

### G0 status

| Check | Result |
|---|---|
| Two preregistered velocity parameterizations, residence cap | Met, with Jobson as a labelled sensitivity (peak tracer velocity, `Qa=Q`) |
| Persistent reach state and conservative update | Met |
| Initial state from saved `Qτ`, recorded as initial state | Met |
| Checkpoints | Replacement hashes, event validated; the registered files are lost |
| Control windows, per-reach records, target supply | Met |
| Rain-only ledger | **Not met**: closes only with stranded overflow named (A 5.77%, B 3.16%, C 0.11%) |
| Full precipitation ledger ≤0.1% | Superseded by the [interim glacier account](#interim-glacier-account-732): holds by construction once #734 closes the rain-only ledger; stability **not met** |

G1 stays closed under RFC 0016's stop rule until both open rows pass.

## Interim glacier account (#732)

Preregistered before any replay. The account is ledger-only, so the expected change to the control
is none.

**Windows.** The G0 checkpoints above, four 100 kyr steps with the ice-age cycle off and the central
river arm (Raymond, one-year cap). Ice-age save: seed 374464662 at 132.6 Myr, SHA-256
`0f8c50c9788684c0ce28050e42a9eb1289f5092ba7fbe3a834c1809ce6a4279c`, the #548 ice-age diagnosis save,
keeping its 8 Myr cycle, for 8 Myr in 100 kyr steps, so one full cycle.

**Arms.** `origin/main` and this branch, run one after the other. A third arm with #734's routing
fix applied measures closure, since that is the condition the ledger criterion names.

| Check | Metric | Threshold |
|---|---|---|
| Control unchanged | per-step lake volume, target ΔV, glacier volume and sea level, branch vs `origin/main` | identical to replay tolerance; expected exactly equal |
| Account identity | `|ΔG − (P_frozen,retained − R_glacier)|` | floating point |
| Partition | land rain + land frozen vs diagnostic precipitation | floating point |
| Budget stability | largest step-to-step `|ΔG|` vs that step's retained input, A, B, C and the ice-age save | reported; the issue asks that swings not exceed inputs |
| Full ledger | residual of `P_diag + M_scaled`, with #734 applied | ≤0.1% of external input, or name the remaining term as a failure |

Expected sizes: model ice divided by `WATER_TO_ICE` makes the opening balance 50,000 times the
model ice volume, so it may well exceed the ocean. Steady ice releases about what the ice retains;
any model-ice swing moves the balance by the same fraction, so large release terms and a failed
stability comparison are possible findings, to be reported rather than tuned (RFC 0016
falsification criteria).

### Results

Source: `origin/main` `cd93d339`, this branch `d0b83742`, and for closure a local merge of this
branch with #734's `425ac919`. Whole-run values are in
[`g0-glacier-account-summary.json`](g0-glacier-account-summary.json). Replay JSON is outside the
repository, under `issue732/` in the private captures.

**The control is unchanged.** On A, B and C, every step's lake volume, target-lake record, glacier
volume, sea level and river inventory are exactly equal between `origin/main` and the branch, a
difference of 0.0. The account has no path into the model.

**The identity and the ledger close.** The account's largest step error is about 5e5 m³ against a
1e21 m³ balance, which is floating point. Land rain plus land frozen precipitation matches
diagnostic precipitation in unit tests.

| Window | Full residual, this branch | Full residual, with #734 | Stranded overflow, branch → with #734 |
|---|---:|---:|---:|
| A | 3.4e-8 | 7.6e-14 | 6.10% → 0 |
| B | 3.3e-5 | 5.2e-14 | 5.25% → 0 |
| C | 1.2e-10 | 8.1e-14 | 0.34% → 0 |

**This table is not evidence for the account.** Closure is true by construction: frozen input
always equals retained plus bare-ground melt, and the account's change plus release cancels
retained, so `full_land_water_residual_m3` reduces to the rain-only `land.residual_m3`. The replays
agree: the two residuals differ by at most 4.6e5 m³, which is floating point. These columns are the
rain-only ledger on a larger denominator, so the improvement with #734 is #734's result. The
branch's residual is under 0.1% only because stranded overflow is a named term, which is not an
allowed sink. The criterion is formally met, but what tests the account is the partition, the
identity, the control comparison and the sizes below. `P_diag` here is land precipitation plus snow
kept by ice on ocean nodes, which is 2–7% of land frozen precipitation.

**Budget stability fails on the mature windows.** The issue asks that the balance not swing more
than its inputs.

| Save | Opening balance (m³) | Largest step change (m³) | Retained input that step (m³) | Change / input | Largest outward / external input |
|---|---:|---:|---:|---:|---:|
| A | 3.79e21 | 6.36e19 | 1.22e18 | **52** | **7.4** |
| B | 3.68e21 | 1.76e19 | 1.01e18 | **17** | **2.3** |
| C | 1.02e21 | 3.16e19 | 1.27e18 | **25** | **3.1** |
| Ice age, 8 Myr | 9.28e19 | 1.42e17 | 4.07e17 | 0.35 | 0.09 |

The ice-age row's ratio is at its largest change; its largest ratio on any step is 0.39, still
under 1.

Two things combine. The opening conversion is the scale: dividing model ice by `WATER_TO_ICE` makes
A's balance about 2,900 times the volume of Earth's ocean (about 1.3e18 m³), and release is a
fraction of that balance, so a 1.7% model-ice loss in A's first step books 6.5e19 m³, seven times
everything that fell in the step. Calving is what moves the fraction.

| Save | Calving share of release, whole run |
|---|---:|
| A | 84% (93% in the first step) |
| B | 64% |
| C | 71% |
| Ice age | 0.5% |

So the mature windows' swings follow sudden model-ice calving, not the melt history. The ice-age
save holds 11–41 times less model ice and barely calves, which is why it stays stable through a full
cycle. Its ledger was not rerun with #734, and it strands up to 13% without it.

**The account has no term for ice changed outside the solver.** `glacier_depth` is remapped by
tectonics, and `step_ice`'s opening volume need not equal the last closing volume, so the account
and model ice can drift apart, which its identity check can't see. The ratio of account to model
ice moved at most 0.06% on A, B and C, but 0.5% over the ice-age cycle, which is comparable to
one step's mean release of about 0.33%. The RFC's whole-world scalar accepts this, and it would
matter over longer runs. If model ice ever reaches zero with a positive balance, that balance
can't be released.

Steady ice would balance, as designed. These windows aren't steady: they lose 1.0–4.0% of the account
in 400 kyr while retained snow adds 0.1–0.5%.

**What this means for G0.** The bookkeeping criterion holds by construction, and the release term
is larger than `P_diag` in 7 of the 12 mature-window steps. The owner's decision is recorded in
[PR #739](https://github.com/adubey/mantle-bloom/pull/739): the account lands as inert bookkeeping
because its identity and control checks pass, the oversized release is a failed G0 interpretation,
and G0 stays blocked, with G1 unavailable, until RFC 0015 resolves the physical glacier-water path.

**A side observation on #734.** Against this branch, its whole-world lake volume differs by up to
84,288 km³ in A and 238,357 km³ in B over the four steps, and glacier volume by up to 77,338 km³ in
B. The target events hold, as #734 reports. G0's rerun on the merged control should record this.
