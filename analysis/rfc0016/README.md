# RFC 0016 diagnosis

## G0 preflight

Status: **in progress; do not start G1.** The owner has bounded G0 and G1 to liquid land water until
RFC 0015 completes. On current source that ledger closes, but only window A survives: #736 removed
window B's event and window C fails the frozen-precipitation screen. G0 needs two replacement
windows, then a rerun. See [Liquid G0 rerun](#liquid-g0-rerun).

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

Source: `origin/main` `cd93d339` against this branch `d0b83742`, before #734 landed. Rerun on
`origin/main` `7d810123`, which has #734, against this branch's merge of it, `4fe51401`. Whole-run values are in
[`g0-glacier-account-summary.json`](g0-glacier-account-summary.json). Replay JSON is outside the
repository, under `issue732/` in the private captures.

**The control is unchanged.** On A, B and C, every step's lake volume, target-lake record, glacier
volume, sea level and river inventory are exactly equal between `origin/main` and the branch, a
difference of 0.0, both before and after #734 landed. The account has no path into the model.

**The identity and the ledger close.** The account's largest step error is about 5e5 m³ against a
1e21 m³ balance, which is floating point. Land rain plus land frozen precipitation matches
diagnostic precipitation in unit tests.

| Window | Full residual, before #734 | Full residual, with #734 | Stranded overflow, before → with #734 |
|---|---:|---:|---:|
| A | 3.4e-8 | 7.6e-14 | 6.10% → 0 |
| B | 3.3e-5 | 5.1e-14 | 5.25% → 0 |
| C | 1.2e-10 | 8.1e-14 | 0.34% → 0 |

**This table is not evidence for the account.** Closure is true by construction: frozen input
always equals retained plus bare-ground melt, and the account's change plus release cancels
retained, so `full_land_water_residual_m3` reduces to the rain-only `land.residual_m3`. The replays
agree: the two residuals differ by at most 4.6e5 m³, which is floating point. These columns are the
rain-only ledger on a larger denominator, so the improvement with #734 is #734's result. The
earlier residual is under 0.1% only because stranded overflow is a named term, which is not an
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
is larger than `P_diag` in 8 of the 12 mature-window steps with #734. The owner's decision is recorded in
[PR #739](https://github.com/adubey/mantle-bloom/pull/739): the account lands as inert bookkeeping
because its identity and control checks pass, the oversized release is a failed G0 interpretation,
and G0 stays blocked, with G1 unavailable, until RFC 0015 resolves the physical glacier-water path.

**A side observation on #734.** Against this branch, its whole-world lake volume differs by up to
84,288 km³ in A and 238,357 km³ in B over the four steps, and glacier volume by up to 77,338 km³ in
B. The target events hold, as #734 reports. G0's rerun on the merged control should record this.

## Why the account releases so much

Asked why release can reach 7.4× input in a balanced ledger. `bin/debug/measure_glacier_release.py`
replays A, B and C on `bf58116f` with the raw `IceLedger`. That commit includes #736's climate noise, so
per-step calving differs from the table above; C's 2.98e19 pulse no longer appears.

**The ledger is correct.** Window A's first step: opening 3.793e21, plus 1.2e18 retained, minus
6.48e19 released, closes at 3.730e21 m³. Release came from opening stock, not that step's snow.

**The stock is on the wrong scale.** Accumulation is scaled to 1/50,000 of real, but stock isn't: sea
level counts A's 7.6e16 m³ of model ice at face value. So the conversion that fits flows gives
about 2,900 oceans of stock, and leaving stock alone makes physical snowfall dwarf any release. No
single conversion fits both. This is what RFC 0015 replaces.

**Model ice is really out of balance.** The account's release fraction copies this ratio.

| Window | Loss / accumulation, by step | Melt only / accumulation |
|---|---:|---:|
| A | 54.8, 9.2, 25.1, 10.8 | 3.4–4.5 |
| B | 20.1, 5.2, 7.4, 11.5 | 3.0–3.8 |
| C | 2.9, 2.5, 2.1, 3.7 | 2.1–3.0 |

The cause of the melt imbalance is unknown. The ice-age cycle isn't it: these saves already have the
cycle off, and a rerun that keeps the save's setting is identical.

**Calving pulses in A and B look like #552.** Stray ice at step start is the largest branch in most
A and B steps, and in those steps most of the ungrounded ice sits on land nodes (78% in A's first
step): columns about 4 km thick, on beds near −4,500 m, with sea level near −740 m. `is_grounded`
assumes seawater over any bed below sea level. Whole columns go at once, so melt tipping one past
the test would explain the pulses, but that isn't shown, and each basin's disconnection from
the ocean is inferred from its land classification. C's calving is on ocean nodes, so #552 doesn't
explain it. Fixing #552 would reduce calving, not the scale mismatch.

## Liquid G0 rerun

Source `50f153ef`, which is simulation-identical to `origin/main` `bf58116f`, and has #734 and #736.
Same 12 arms, commands and replacement checkpoints as [G0 control replay](#g0-control-replay), with
the 10% frozen limit registered in PR #743 before this ran. Whole-run values are in
[`g0-liquid-summary.json`](g0-liquid-summary.json). Replay JSON and per-reach arrays are under
`g0-50f153ef/` in the private captures; their `SHA256SUMS` has SHA-256
`e775bcc1bfae904d6561b0828da22d15d21b4151083094ac304bfcf3b70cd438`.

**The river and liquid gates hold.**

| Window | Inventory, central (km³) | Jobson (km³) | Largest reach error (m³) | Liquid residual, worst step | Stranded overflow |
|---|---:|---:|---:|---:|---:|
| A | 229.5→251.0 | 460.6→543.7 | 44 | 2.0e-13 | 0 |
| B | 239.6→260.2 | 454.9→562.4 | 62 | 2.9e-14 | 0 |
| C | 438.4→467.9 | 719.8→770.5 | 44 | 2.1e-13 | 0 |

The 0.1- and 10-year caps match the central inventory to the reported precision. Jobson raises it
1.6–2.2×. Neither moves the target lake's supply noticeably: river storage releases about 2e4 m³/yr at
most, against 3.5e12 m³/yr of rain on B's lake.

**The windows don't hold.** Target-component volumes, km³. Replays before #736 are the #732 runs with
#734 applied. #736 is the only simulation change between them and this run.

| Window | Checkpoint | Before #736, steps 1–4 | This run, steps 1–4 | Event |
|---|---:|---|---|---|
| A | dry | dry, 879,183, dry, 930,760 | dry, 863,389, 783,206, 856,022 | Fills as before, but no longer dries at step 3 |
| B | 259,257 | 115,361, 246,986, 251,732, 251,451 | 252,737, 257,863, 254,322, 254,560 | **Gone**: no drop and refill |
| C | 137 | 264,787, 3,200, 372,572, 2,950 | 264,787, 3,201, 372,780, 2,966 | Intact |

Whole-world lake volume still swings widely, for example 1.28 to 3.44 million km³ in B's last two
steps, so flicker continues. It just isn't these events.

**Frozen-precipitation screen,** limit 10%, summed over steps with a wet anchor:

| Window | Frozen share of catchment precipitation | Worst step | Result |
|---|---:|---:|---|
| A | 1.0% | 1.0% | Pass |
| B | 1.6% | 1.7% | Pass, but its event is gone |
| C | 18.1% | 22.0% | **Fail** |

C's measure is a lower bound on its uncertainty. No routed supply lands on C's target nodes in any
step, so its 370,000 km³ swings arrive some other way, probably from a neighbouring basin across a
shared spill. That isn't traced, and the catchment follows `flow_target` only, so it may miss
that basin's water. C fails either way.

**G0 conclusion.** River law, reach closure, initial state and liquid ledger pass on this source.
The windows don't: only A is usable, so G1 needs two replacements, chosen and hashed before any
treatment under PR #743. Rerun G0 on them.

## Window rescan

Preregistered before any scan run. PR #743 allows replacing a window whose event is lost or which
fails the frozen-precipitation screen. This scan selects windows on current source and measures
whether #736 changed how often flicker occurs.

**Arms.** Current source `0bda513d`, and `9d6fb5e2`, the last `main` commit before #736, which has
#734. Both replay in series with `bin/debug/measure_lake_step_budget.py`, 100 kyr steps, ice-age
cycle off, for 5 Myr.

**Worlds.** #680's two event worlds, by SHA-256: seed 910211954 at 220.2 Myr (`bd487409…`) and seed
317739385 at 200.8 Myr (`6c806f91…`). Neither is an RFC 0016 held-out save.

**Detector.** #680's unchanged #588 rule: a matched component over 10,000 km², whose successive
volume changes reverse sign, both above 50% of the prior volume. Results from the first 0.4 Myr
after load are excluded in both arms, because the current arm's noise changes at load.

**Flicker comparison, descriptive only.** Reversal counts and comparable-body counts per world and
arm, with the qualifying step ages. Two worlds and one noise history can't separate a real change
from chance, so any difference is reported as an indication, not a result.

**Window selection, current arm only, in order.**

1. Candidates are qualifying reversals. A window is four steps from a checkpoint two steps before
   the reversal, as #680 used. Windows don't overlap.
2. Rank by largest absolute `ΔV`, then take the matched component as the target, with its anchor at
   the wet node nearest its centroid.
3. Save the checkpoint from the source world, and require that a reload reproduces the event.
4. Pass RFC 0016's frozen-precipitation screen: at most 10% across the window.
5. Require nonzero routed supply to the target in the window. G1 tests an upstream store, and window
   C shows a lake with none can't respond through that path. This rule is new, is a selection rule
   rather than a gate, and is disclosed for the owner's review.
6. Take the three highest-ranked windows that pass, with at least one from each world if both
   supply a passing candidate. A registered window keeps no priority: it qualifies only by this rule.

If fewer than three pass, report that and stop; don't relax a threshold. Rerun G0's 12 arms on the
selected windows and record them with hashes before any treatment.


### Rescan results

Source `2e465887`, simulation-identical to `origin/main` `0bda513d`. Whole-run values are in
[`window-rescan-summary.json`](window-rescan-summary.json). Scans, checkpoints and G0 replay JSON are
under `window-rescan/` in the private captures; their `SHA256SUMS` has SHA-256
`ad8ba73bd59e2bc38411d2e5f4f8d4f1f4aaaffdea8667ed72f5b68a8943bc42`.

**Flicker continues after #736.** Steps 5–50 of each scan, so results from the first 0.4 Myr are out.

| World | Arm | Reversals | Comparable bodies | Rate | Steps with a reversal |
|---|---|---:|---:|---:|---:|
| Seed 910211954 | Before #736 | 27 | 271 | 10.0% | 23 of 46 |
| Seed 910211954 | Current | 23 | 239 | 9.6% | 18 of 46 |
| Seed 317739385 | Before #736 | 27 | 106 | 25.5% | 26 of 46 |
| Seed 317739385 | Current | 23 | 110 | 20.9% | 23 of 46 |

The totals are 54 and 46 reversals. As registered, this is an indication only. It doesn't support the
reading that #736 stopped flicker: B's event was lost, but the rate barely moved, and reversals now
fall on different steps. Seed 910211954's save is at 220.4 Myr, not the 220.2 Myr in its file name.

**Candidates.** Eight non-overlapping windows per world, ranked by `|ΔV|`. The top four of seed
910211954 and all eight of seed 317739385 were screened. Every screened reload reproduced the scan's
volume exactly.

| Seed | Checkpoint (Myr) | `|ΔV|` (km³) | Frozen share | Routed supply | Result |
|---|---:|---:|---:|---|---|
| 910211954 | 225.1 | 2,290,848 | 4.3% | Yes | **Window D** |
| 910211954 | 224.4 | 2,288,695 | 51.0% | Yes | Fails frozen screen |
| 910211954 | 224.0 | 2,162,218 | 4.5% | Yes | **Window E** |
| 910211954 | 222.8 | 1,736,163 | 7.0% | Yes | **Window F** |
| 317739385 | All eight | 357,534–383,913 | 19.1–21.6% | None | Fail both |

Seed 317739385 supplies no passing candidate, so the one-per-world clause doesn't apply. Its target
is the lake from window C: no routed supply reaches it in any step. Candidate 224.4's share comes
from one step whose catchment reaches 29 million km².

**Selected windows.** All three start from seed 910211954 and run 0.4 Myr at 100 kyr steps. Each
anchor stays wet on all four steps.

| Window | Checkpoint SHA-256 | Target volume, steps 1–4 (km³) | Anchor file |
|---|---|---|---|
| D | `6a8d6f90…` | 2,293,616, 2,768, 268,996, 788,971 | [`d`](target-anchor-window-d.json) |
| E | `4045b915…` | 106,435, 2,268,652, 6,422, 2,271,028 | [`e`](target-anchor-window-e.json) |
| F | `1f9f2493…` | 1,808,528, 72,365, 694,236, 970,045 | [`f`](target-anchor-window-f.json) |

**G0 on the selected windows passes.** Twelve arms: Raymond at 0.1-, 1- and 10-year caps, and Jobson
at 1 year.

| Window | Inventory, central (km³) | Jobson (km³) | Largest reach error (m³) | Liquid residual, worst step | Stranded overflow | Frozen share, window | Worst step |
|---|---:|---:|---:|---:|---:|---:|---:|
| D | 296.8→284.9 | 296.8→489.4 | 64 | 1.3e-14 | 0 | 4.3% | 6.3% |
| E | 273.1→268.9 | 273.1→519.7 | 32 | 8.1e-13 | 0 | 4.5% | 6.2% |
| F | 245.8→254.1 | 245.8→434.5 | 32 | 2.6e-14 | 0 | 7.0% | 10.8% |

F passes because the screen applies across the window, but one step exceeds 10%. Target volumes match
across all arms of a window to the nearest km³.

**Limits on what these windows can show.**

- All three anchors lie within about 290 km of each other, so "two of three windows" would rest on
  one lake basin in one world.
- The windows replace A as well as B and C. RFC 0016 allows replacement only for an input or
  data-integrity reason, a failed frozen screen or a lost event, and A passes the frozen screen and
  still fills. The rescan's rule that a registered window keeps no priority therefore needs the
  owner's decision, as does the routed-supply rule. A's G0 arms at `50f153ef` are
  simulation-identical to this source, so A can stand in without a rerun. Its anchor is about 1,500 km
  from D's.
- Each event's `|ΔV|` is about 2% of the target's supply over one 100 kyr step, and world river
  inventory is about 0.013% of `|ΔV|`. A river store of this size can't hold the swing, while a change
  of about 0.2% in delivered water could meet the 10% materiality limit. G1's result will turn on
  the flux, not on the storage.
