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
  owner's decision. The owner approved the routed-supply rule. A's G0 arms at `50f153ef` are
  simulation-identical to this source, so A can stand in without a rerun. Its anchor is about 1,500 km
  from D's.
- Each event's `|ΔV|` is about 2% of the target's supply over one 100 kyr step, and world river
  inventory is about 0.013% of `|ΔV|`. A river store of this size can't hold the swing, while a change
  of about 0.2% in delivered water could meet the 10% materiality limit. G1's result will turn on
  the flux, not on the storage.

**Window A's status.** Its registered event survives once, not as a repeating flicker. On current
source, `bin/debug/measure_lake_step_budget.py` from A's checkpoint for 0.8 Myr finds one qualifying
reversal on A's lake, at 221.8 Myr: 424,670→863,389 km³, or `ΔV` +438,719 km³ against the registered
+495,624. The lake then holds at 783,206–856,022 km³ through the window, with no further reversal
on it. The anchor node is dry at step 1 although the matched lake is present, because the anchor is a
location rather than a lake identity. A's checkpoint came from an older source, so its history
differs from the rescan's. Whole-world lake volume at 221.8 Myr is 2.54 million km³ from A's
checkpoint and 1.64 million km³ in the scan. In the scan's history, A's anchor reads dry, dry, 991,432,
dry, dry, 768,249 km³ from 221.6 to 222.1 Myr, so the location flickers, but the detector records no
qualifying reversal there.

### Windows for G1

The owner chose A, D and E on 2026-11-02. This keeps RFC 0016's replacement rule and three-window
criterion unchanged: D replaces B, whose event was lost, and E replaces C, which fails the frozen
screen. A stays, so the registered-window-has-no-priority clause of the rescan's rule is set aside
for it. F isn't used.

| Window | Replaces | Checkpoint SHA-256 | Frozen share | G0 source |
|---|---|---|---:|---|
| A | Stays | `9882568e681a0f9f207e3302453b3a92f18865d5640d7baa8d9d8eb09ab4633e` | 1.0% | `50f153ef` |
| D | B | `6a8d6f90…` | 4.3% | `2e465887` |
| E | C | `4045b915…` | 4.5% | `2e465887` |

G1's pinned source is `0bda513d`, which is simulation-identical to both G0 sources. Anything later,
including #745's starting rock ages, is unmeasured against these windows, so it needs a G0 rerun
before it can replace that pin. Full hashes for D and E are in
[`window-rescan-summary.json`](window-rescan-summary.json).

G0 passes on these three windows: river law, reach closure, initial state, liquid ledger and the
frozen screen. Two limits carry into G1. D and E share a basin, so the two-of-three test has two
basins behind it, and A's checkpoint comes from an older source's history, with a one-way fill.

## G1a preregistration

Frozen by the merge of its PR, before any paired run. Nothing here has run a treatment: the scale
check below only reads each checkpoint's saved climate. Design principles 3, 5 and 8 apply: the
store is a fast process inside a long step, every flux has physical units, and the land law replaces
a proxy with its own physics.

**Why the store is seasonal.** A store of about a metre reaches balance within a few years, far
shorter than a 100 kyr step. With annual-mean forcing the balance is `ET = min(L, Eᵖ)`, whatever the
capacity, so `C` would have no effect at all. Capacity changes the annual water balance through
seasonality [[Milly, 1994]](https://doi.org/10.1029/94WR00586). So the law below runs a monthly bucket
to its repeating year.

**Inputs.** The model has no net radiation, leaf area, diurnal range or cloud, so Penman-type,
Priestley–Taylor and Hargreaves methods can't run without an invented input. It has annual
temperature `Ta`, seasonal amplitude `a` (`biomes.seasonal_temp_amplitude`), axial tilt, the
half-year rain share `s` (`biomes._precip_season`) and the Köppen class.

**Equations, for each dry, unfrozen land node and month `m`.** All depths are mm over the month.

1. Monthly temperature: `Ta_m = Ta + a cos(2π(m − m_w)/12)`, with `m_w` a month after the summer
   solstice.
2. Rain: `L_m = L s/6` in the warm half-year and `L (1 − s)/6` in the cold half. `L` is the node's
   annual liquid rain, so `Σ L_m = L`.
3. Demand: `Eᵖ_m = D · R_m/(λρ) · max(Ta_m + 5, 0)/100 · d_m`, where `R_m` is monthly
   top-of-atmosphere radiation in MJ m⁻² d⁻¹ from the formula in `climate.compute_insolation_w_m2`,
   `λρ` = 2.45 MJ m⁻² mm⁻¹ and `d_m` is days in the month
   [[Oudin et al., 2005]](https://doi.org/10.1016/j.jhydrol.2004.08.026). It needs only these inputs, and
   was proposed for rainfall-runoff modelling after a test of 27 methods on 308 catchments. `D` is the demand factor. The formula replaces `open_water_evaporation_mm` in the
   treatment only.
4. Water limit: `β(S) = min(1, S/((1 − p) C))`, the FAO-56 stress form
   [[Allen et al., 1998]](https://www.fao.org/3/x0490e/x0490e0e.htm), for `S` in `[0, C]`.
5. Store: `dS/dt = L_m/Δ − (Eᵖ_m/Δ) β(S)`, with `Δ` the month length. Within a month this is linear
   above `(1 − p) C` and exponential below it, so it's integrated in closed form with the crossing
   solved exactly. `S` above `C` leaves as saturation excess, the only surface water.
6. Actual ET: `ET_m = ∫ Eᵖ_m β dt`, which can't exceed `S + L_m`. Each flux is debited once.
7. Split: `T_m = f ET_m` and `ET_soil,m = (1 − f) ET_m`. `f` is by Köppen class, below. Wei's ratio
   includes interception, so `ET_soil` covers it too. `f` is a measured actual-flux ratio applied to
   the whole of `ET`, not a stress-dependent split.
8. Interval: step the year from the node's opening `S` until it repeats, to 1e-6 of `C`, booking
   those years explicitly. Then multiply the repeating year by the remaining years. Fluxes therefore
   don't depend on step length once the year repeats, and the interval closes exactly.
9. Opening state: at each checkpoint, `S` starts on the repeating year's January value under the
   checkpoint's climate. That's recorded as initial inventory, not an input, as with `Q τ`.
10. Humidity: on the next climate update, `T` replaces `_vegetation_transpiration_source` on land,
    at `T_mm/PRECIP_MM_PER_G_KG` g/kg, area-binned from nodes to the climate grid.
11. Unchanged: wet-fraction rain, the `is_frozen` mask, `M_scaled`, routing, river reservoirs and
    lakes. Frozen nodes keep control behaviour.

**Parameters.** `C` is uniform across nodes, because the model has no field to justify a pattern.

| Parameter | Low | Central | High | Source, and what's a design choice |
|---|---:|---:|---:|---|
| `C` (mm) | 150 | 300 | 1,000 | Low: [Gupta et al., 2023](https://doi.org/10.1029/2022MS003277), 19,941 km³ in the top 1 m, over about 1.3e8 km² of ice-free land (the area is my assumption). Central: twice that for 2 m, which [Stocker et al., 2023](https://doi.org/10.1038/s41561-023-01125-2) exceed across 37% of vegetated land. High: **design choice**, for deeper water, with no global value behind it. |
| `p` | 0.30 | 0.50 | 0.70 | FAO-56's stated range and common value, for crops; applied to natural land as an assumption. |
| `D` | 0.85 | 1.00 | 1.15 | **Design choice.** No published uncertainty transfers to this formula here. Scaled from the 1.08–1.34 range quoted for the Priestley–Taylor coefficient of 1.26, widened to be symmetric. |
| `f` shift | −0.14 | 0 | +0.18 | [Wei et al., 2017](https://doi.org/10.1002/2016GL072235): global 57.2%, ensemble 42.9–74.9%. Applied to every class, clipped to `[0, 1]`. |

| Köppen group | `f` | Source |
|---|---:|---|
| A | 0.61 | Wei table 1, broad-leaf forest |
| C, D | 0.55 | Wei, needle-leaf and mixed forest |
| BS, ET | 0.65 | Wei, shrubs and grasses, which includes tundra |
| BW | 0.25 | **Design choice.** Wei maps deserts as low and gives no value |
| EF | 0 | No plants |

**Windows and source.** A, D and E on `0bda513d`, with the checkpoints under
[Windows for G1](#windows-for-g1), at 100 kyr steps for 0.4 Myr, and the registered 10 kyr
connectivity trace. No held-out save is used.

**Arms.** Each is paired to the same control, with the same noise history.

| Arm | Change | Role |
|---|---|---|
| Central | Central values | Primary |
| `C` low, high; `p` low, high; `D` low, high; `f` low, high | One at a time | Range: eight arms |
| Partition match | Control, with each node's runoff scaled uniformly by the treatment's land-wide excess ratio for that step, no store, control humidity | Reported, not counted |
| Coupling off | Central, with the control humidity source | Reported, not counted |
| Sham ±0.2% | Control, with runoff scaled by 1.002 and 0.998 | Reported, not counted |

**Readouts.** Those in RFC 0016's G1b, per arm and window, plus each step's land-wide excess ratio and
the humidity source against the control's. The liquid ledger must close to 0.1% in every arm, with
`S` as a reservoir. The materiality rule is unchanged.

### Scale check

From [`g1a-scale-check.json`](g1a-scale-check.json), using `bin/debug/measure_rfc0016_land_demand.py`
on unfrozen land on the climate grid. The private `g1a/SHA256SUMS` has SHA-256
`9505c745a97af1f18b8c7060b72110e386253c93f31988078369f98fbec6de71`. It uses sub-stepped integration,
not the closed form, and it's a scale check: no parameter was chosen from it.

These numbers replace a first version whose warm half-year had seven months, so rain wasn't conserved
(#751). The law's equation 2 was always six; only this script's and the treatment's code were wrong.

| Window | Rain (mm/yr) | Oudin demand (mm/yr) | Control demand (mm/yr) | Treatment ET (mm/yr) | Excess (mm/yr) | Control runoff (mm/yr) | Change |
|---|---:|---:|---:|---:|---:|---:|---:|
| A | 839 | 1,002 | 1,130 | 485 | 354 | 389 | −9.2% |
| D | 785 | 998 | 1,127 | 454 | 331 | 362 | −8.5% |
| E | 758 | 998 | 1,126 | 446 | 311 | 343 | −9.2% |

Central transpiration is 247–269 mm/yr, or 32–33% of rain, within the 39 ± 10% of a compilation of 81
studies ([Schlesinger & Jasechko](https://agu.confex.com/agu/fm14/webprogram/Paper2691.html),
conference abstract), for Earth, which is a sanity check only.

Across the range, excess changes by +4.9% to −21.4% from the control:

| Parameter, low to high | Change in excess across windows |
|---|---:|
| `C` | −5.2% to −11.2% |
| `p` | −7.4% to −10.0% |
| `D` | +4.9% to −21.4% |

**What this means for G1b.**

- **The treatment will probably pass materiality for a reason unrelated to storage.** Each event's
  `|ΔV|` is about 2% of the lake's supply in a step, so a 10% change needs about 0.2% more or less
  delivered water, and the central arm moves land-wide runoff by about 9%. The law's central
  partition differs from Budyko's, and the store can't buffer across steps: at 300 mm, D's 3.3
  million km² catchment holds about 1,000 km³, or 0.04% of its `|ΔV|`. That's why the partition-match
  and sham arms are registered. If they're also material, G1b shows these events respond to any small
  change in water supply, not that storage or plant use caused them.
- **The humidity source roughly triples,** at 2.5–3.2× the control, because the control's proxy is
  tuned to 16% of rain at most. RFC 0016 says to pass the volume, so this registers that as written.
  It's the largest forcing in the package, and coupling off separates it.
- **The step from control to treatment is a shock.** Both arms open from a control checkpoint, so 0.4
  Myr of treatment includes the lake system's response to a sudden change in supply.

**Sources not read in full.** Oudin's equation is from later reproductions of it and FAO-56's `Ks`
from its description of equation 84. The primary pages could not be opened; the method's
catchment counts are from the abstract. Stocker's figure and the Priestley–Taylor range are from
search summaries. Check these before citing them outside this report.

## G1b result

Run on `0bda513d` with the G1a law as frozen, after the warm-half correction in #750, on the three
windows' registered checkpoints. Each window has a control and 13 arms, on two schedules: four
100 kyr steps (production) and forty 10 kyr steps (the connectivity trace). That's 84 runs, none of
which failed. Scores are in [`g1b-summary.json`](g1b-summary.json), from
`bin/debug/score_rfc0016_g1b.py`. The treatment is `bin/debug/rfc0016_land_store.py`, installed by
`replay_rfc0016_river_water.py --treatment`. The run outputs are private. Their `SHA256SUMS` files
have SHA-256 `a158250f…cf1f0` (production) and `95ce6a38…c3ce1` (trace).

```bash
backend/.venv/bin/python bin/debug/replay_rfc0016_river_water.py <checkpoint.mbworld> \
  --anchors analysis/rfc0016/target-anchor-window-<w>.json --myr 0.4 --step-years <100000|10000> \
  --checkpoint-myr 0.1 [--treatment [--capacity-mm C] [--depletion p] [--demand-factor D] \
  [--split-shift f] [--no-humidity] | --runoff-scale 1.002 | --runoff-scale-from <w>-central.json] \
  --out <captures>/<w>-<arm>.json
python3 bin/debug/score_rfc0016_g1b.py --production <100 kyr dir> --trace <10 kyr dir> \
  --out analysis/rfc0016/g1b-summary.json
```

**Accounting closes in every arm.** The worst normalized land residual is 1.1e-12 on production and
6.4e-12 on the trace, against a 0.1% limit. Stranded overflow is zero. The store's identity error
is at most 2.1e-7 mm. Controls reproduce G0's targets exactly. The three events are A +438,718 km³,
D −2,290,848 km³ and E +2,162,218 km³.

**Event readout.** The event's `ΔV` is the target body's volume minus its volume at the previous
step, matched under #588's rule (largest overlap, at least half the smaller body). A's anchor is
dry the step before its event although its lake is present, so the anchor alone gives a different
lake's volume. This readout reproduces the detector's 424,670→863,389 km³ for A. It was added after
the first production matrix, which was then rerun in full; the preregistration names the event,
not how it's read.

### Scores

Event change is `|ΔV|` against the control. The trace column is the longest run of consecutive
10 kyr outputs where the anchor's wet/dry class differs from the control's, out of 40.

| Arm | A event | D event | E event | A trace | D trace | E trace | Windows | Material |
|---|---:|---:|---:|---:|---:|---:|---|---|
| Central | +6.15% | −0.00% | −0.00% | 1 | 1 | **9** | E | No |
| `C` 150 | +0.04% | −0.00% | −0.00% | 1 | 1 | **9** | E | No |
| `C` 1,000 | **+13.69%** | +0.00% | −0.00% | 1 | 1 | **9** | A, E | **Yes** |
| `p` 0.3 | +2.05% | −0.00% | −0.00% | 1 | 1 | **9** | E | No |
| `p` 0.7 | **+10.07%** | −0.00% | −0.00% | 1 | 1 | **6** | A, E | **Yes** |
| `D` 0.85 | +1.28% | −0.00% | −0.00% | 1 | 1 | **8** | E | No |
| `D` 1.15 | **+10.11%** | −0.00% | −0.00% | 1 | 1 | **4** | A, E | **Yes** |
| `f` low | +6.15% | −0.00% | −0.00% | 1 | **3** | **6** | D, E | **Yes** |
| `f` high | +6.15% | −0.00% | −0.00% | 1 | 1 | **9** | E | No |
| Partition match | +0.54% | −0.00% | −0.00% | 1 | 1 | **8** | E | Not counted |
| Coupling off | +6.15% | −0.00% | −0.00% | 1 | 1 | **2** | E | Not counted |
| Sham +0.2% | −0.00% | −0.00% | +0.00% | 1 | 0 | **9** | E | Not counted |
| Sham −0.2% | +0.00% | +0.00% | −0.00% | 1 | 1 | **3** | E | Not counted |

Read as written, the rule makes the package material: four range arms reach two windows.

### Why E can't separate a treatment from noise

**E's trace fails in every arm, including both shams.** The anchor alternates wet and dry every
10 kyr for much of each trace, with stretches where it stays wet. A small change slides one arm's
alternation by one output against the other's, or ends a wet stretch early, and every output after
that differs. A change of 0.2% in runoff, with no store, does this to E as readily as the
treatment does:

```text
E control  W.W.WWWWW.W.W.W.WWWWWWWW.W.WWWWW.W.W.W.W
E sham +   W.W.WWWWW.W.W.WWWW.WWWWWWWWWWWW.W.W.W.W.
E central  W.W.W.WW.W.W.W.W.W.WW.W.W.W.WWW.W.W.W.W.
```

The sham was registered for this question: if a perturbation with no mechanism behind it is
material, the screen measures sensitivity rather than mechanism. On E's trace it is. Without E's
trace, no arm reaches two windows: the four material arms are A's event plus E, or D's trace (a
run of 3, in `f` low only) plus E.

**D's and E's events don't respond to supply.** Every arm leaves both within 0.002%, although
land rain rises 30–50% from the second step (the humidity source, which coupling off removes).
D empties to 2,768 km³ and E fills to 2.27 million km³ whatever reaches them. Their later steps
stay within 0.2%, except D's refill at step 3, which `C` 150 and `f` high raise 11% (on about
269,000 km³). That fits the case [#754](https://github.com/adubey/mantle-bloom/issues/754)
investigates, that these lakes switch on how they connect and spill, not on their supply.

**A's event responds to where the law takes runoff.** The `f` arms and coupling off leave
it at central's +6.15%, because the humidity source acts from the next climate update, so the
event sees only the first step's partition. Partition match takes the same land-wide share as
central, uniformly, and moves it +0.54%. So it's the law's pattern across A's catchment, not the
land-wide total, that moves A, through the treatment's first-step shock. Across arms the response
mostly grows as the first step's runoff falls, from +0.04% (`C` 150, ratio 0.97) to +13.69% (`C`
1,000, 0.84), but not monotonically: `D` 1.15 cuts most (0.75) and gives +10.11%. On production
A's anchor then dries at steps 3–4 in every treatment arm except `C` 150, and in partition match,
where the control and shams stay wet. That's after the registered event, so it isn't scored.

### What G1b shows

- On the rule as written, the package is **material** in four range arms, so it isn't
  deprioritized as a cause.
- Each of those arms depends on E's trace, and E's trace fails under the ±0.2% shams. No arm is
  material in two windows on evidence a sham doesn't also produce.
- Two of the three events (D and E, which share a basin) don't respond to a 30–50% rise in land
  rain. Neither storage nor plant water use can move them through supply.
- A's event responds to the law's runoff pattern, by up to 14%, not to storage: the store holds
  about 0.04% of an event's `|ΔV|`, and `C` 150 moves A under 0.1%.

The materiality rule is the owner's (#680), so reading E's trace as uninformative is a decision for
the owner, not this report. It isn't applied above.

**The coupled arms overstate recycling.** Rain is read from humidity and never debited from it,
and carried moisture decays by a fixed half-life whatever has rained, so transpiration passed to
humidity becomes rain at its cell and again downwind. Land rain rises 38–51% in the coupled arms and
water to the ocean rises with it; coupling off, which keeps the control's source, moves land rain
by under 2%. Until the climate conserves moisture
([#771](https://github.com/adubey/mantle-bloom/issues/771)), only coupling off approximates
recycled water that is redistributed rather than multiplied, so G1b can't evaluate the coupling.

### Assumptions to disclose

- The breach gate uses the control's Budyko runoff in every arm; the treatment changes only the
  water partitioned to rivers.
- The first step's humidity is the control's, because the treatment's transpiration enters at the
  next climate update.
- Each node's opening year starts in month 0, and pelagic cells take `f` 0.55.
- The D and E trace arms were rerun after a script error ran them without their flags, as
  controls. Their first outputs are kept apart from the results.

## Window re-baseline (#791)

Preregistered before any scan run. [G1b](#g1b-result) ran on windows A, D and E at `0bda513d`.
[#754](https://github.com/adubey/mantle-bloom/issues/754) traced D's and E's events to lake-ice
sublimation charged outside the lake solver, and #757 and #758 have since merged. This rescan
asks which target events exist on current source and are clean enough to test a land-water store.
It selects windows only: no treatment arm runs, and G0's 12 river arms rerun once, at the final
pin, before G1b.

**Source.** `origin/main` at `3854792f`, which includes #757 and #758. Every scan runs on this
commit. If `main` moves before the scan ends, the scan stays on `3854792f`; the final pin for G0
and G1b is chosen later and needs its own G0 rerun if it differs.

**Worlds.** The two #680 event worlds, by SHA-256: seed 910211954 at 220.4 Myr (`bd4874092efc805e8b174e115fe567c82d35fc2ee606686938b5b6eb96b060d2`) and
seed 317739385 at 200.8 Myr (`6c806f913e4c79bb351b5eb5a143abcfc273df39b3bdf4ce9823fec7d27fb0cd`), plus window A's
221.6 Myr checkpoint (`9882568e681a0f9f207e3302453b3a92f18865d5640d7baa8d9d8eb09ab4633e`) replayed
for 0.8 Myr. Neither seed is in RFC 0015's held-out set (447276134, 628977350, 295569319), and RFC 0016
hasn't named one. Adding a world is an owner decision, taken before the scan.

**Scan.** `bin/debug/measure_lake_step_budget.py`, 100 kyr steps, ice-age cycle off, 5 Myr per
world. The detector is #680's unchanged #588 rule: a matched component over 10,000 km² whose
successive volume changes reverse sign, both above 50% of the prior volume. Results from the first
0.4 Myr after load are excluded.

**Candidates and ranking.** Qualifying reversals, each a four-step window from a checkpoint two
steps before the reversal, as in the [earlier rescan](#window-rescan). Windows don't overlap. Rank by
largest `|ΔV|`; the matched component is the target and its anchor is the wet node nearest its
centroid.

**Screens, applied in order to ranked candidates.** A candidate that fails one is recorded with the
reason and the next is taken.

1. **Reload.** A checkpoint saved from the source world reproduces the event's volumes exactly.
2. **Frozen precipitation.** At most 10% across the window, over the catchment of every land node
   whose `flow_target` path reaches the target. Same rule as RFC 0016's G0.
3. **Routed supply.** Nonzero routed supply to the target in every step of the window. A
   selection rule, approved by the owner in the earlier rescan.
4. **Artifact screen (new).** The event must not be one the lake solver's accounting produces.
   - The normalized liquid-ledger residual is at most 0.1% in every step.
   - A knockout of lake-ice sublimation, as in #754's K-sub arm, rerun over the window, leaves the
     event at no less than 50% of its `|ΔV|`. This is a driver patch on the replay, not a model
     change. An event that disappears or shrinks below half without sublimation is recorded as a
     sublimation artifact and excluded. The 50% limit is a design choice set before measurement,
     not a literature value.
5. **Basin independence (new).** The three windows' targets share no node set in any step of
   their windows, and no two anchors are within 1,000 km. Both limits are design choices.

**Selection.** The three highest-ranked candidates that pass every screen. Window A competes under
the same screens and keeps no priority; if it passes, it is ranked with the rest. The owner
approved A in the earlier rescan, and this sets that approval aside for the new pin, so the
owner's review applies.

**If fewer than three pass, report which failed which screen and stop.** No threshold is relaxed
after measurement.

**Report.** Reversal counts and comparable-body counts per world on current source, with the
qualifying step ages (descriptive only: two worlds and one noise history can't separate a real
change from chance); whether D's, E's and A's events survive #757 and #758; each selected window's
checkpoint SHA-256, anchor file, target volumes and frozen share. Compact JSON is checked in. Raw
scans, knockout runs and checkpoints stay in the private captures with a `SHA256SUMS` file.

**What this can't show.** Selection ranks by `|ΔV|`, so the windows favour the largest swings, as
before. Two worlds supply the candidates. Passing the artifact screen shows the event survives
removing sublimation, not that supply or storage drives it.

### Re-baseline results

Scans and candidate reloads used pinned source `3854792fad7c54fe97078ee490625f488273bd00`.
The seed 910 save is named `...-220200000y.mbworld`, but its embedded age is 220.4 Myr; its
SHA-256 matches the preregistered input. The first four scan outputs were excluded.

| World | Reversals | Comparable-body observations | Event ages (Myr) |
|---|---:|---:|---|
| Seed 910211954 | 14 | 469 | 221.1, 221.3 (2), 221.7, 221.8, 221.9, 222.6, 222.7 (2), 222.8, 223.9, 224.0, 224.8, 224.9 |
| Seed 317739385 | 5 | 178 | 202.1, 204.4, 205.4 (2), 205.5 |

The separate 221.6 Myr A checkpoint replay had no qualifying reversals in its eight outputs.
The former A event is absent at its registered anchor. The nearest current-source reversals are
4,970 km away at 221.7 Myr and 9,953 km away at 221.8 Myr. D's event is absent around its former
225.3 Myr timing. E's nearest reversal around its former 224.2 Myr timing is 8,514 km from E's
registered anchor. Thus A, D and E do not supply the former target events on this pin.

All eight nonoverlapping ranked candidates reloaded with exact before, after and `ΔV` values.
Candidate ranks below are global, by descending `|ΔV|`. Frozen share is weighted by precipitation
volume over the outputs where the detector-matched component exists; `n/4` shows how many of the
four outputs had a matched target. A missing target fails the every-step routed-supply rule.
Every candidate's whole-lake liquid-ledger residual passed the 0.1% limit.

| Rank | Seed, checkpoint (Myr) | `|ΔV|` (km³) | Frozen share | Routed supply | Result |
|---:|---|---:|---:|---|---|
| 1 | 910211954, 221.1 | 488,100 | 29.59% (4/4) | Every step | Fails frozen screen |
| 2 | 910211954, 222.5 | 460,794 | 5.55% (3/4) | Missing target at step 4 | Fails routed-supply screen |
| 3 | 910211954, 221.5 | 390,911 | 11.55% (3/4) | Missing target at step 3 | Fails frozen and routed-supply screens |
| 4 | 317739385, 201.9 | 44,350 | 0.00% (2/4) | Missing target at steps 3–4 | Fails routed-supply screen |
| 5 | 910211954, 224.6 | 7,644 | 2.48% (4/4) | Every step | Passes individual screens, including K-sub |
| 6 | 317739385, 205.2 | 5,960 | 40.20% (3/4) | Missing target at step 4 | Fails frozen and routed-supply screens |
| 7 | 317739385, 204.2 | 5,310 | 42.93% (2/4) | Missing target at steps 3–4 | Fails frozen and routed-supply screens |
| 8 | 910211954, 223.8 | 2,782 | 86.03% (2/4) | Missing target at steps 3–4 | Fails frozen and routed-supply screens |

**The sole individual passer is not enough to select G1 windows.** Candidate 5 in the global ranking
(seed 910211954, checkpoint 224.6 Myr) has checkpoint SHA-256
`0b0470a7f08ce0cf1a5d06ae172686b2f7d5f0039bd19fd430abed547ee45c19`, target volumes
3,456.502, 11,100.869, 559,447.829 and 1,468.931 km³, and a 2.48% frozen share. Its anchor is
[`target-anchor-window-rescan-2-single-pass.json`](target-anchor-window-rescan-2-single-pass.json).
The maximum normalized liquid-ledger residual is `3.44e-13`. With lake-ice sublimation disabled,
the matched target-body change is 5,076.636 km³, 66.41% of the control `|ΔV|` of 7,644.367 km³,
above the preregistered 50% cutoff.

Only one candidate passes all individual screens, so the three-window basin-independence screen
cannot be satisfied. As preregistered, selection stops here: no G0 river arms or G1 treatment arms
were run, and no threshold was relaxed. The compact values, all anchors, full checkpoint hashes
and per-step screens are in [`window-rescan-2-summary.json`](window-rescan-2-summary.json).
Raw scans, replay outputs and checkpoints remain in the private capture under `rfc0016/window-rescan-2/`
and are covered by its `SHA256SUMS` file.

## Expanded window rescan (#871)

Preregistered on `codex/issue-871-expanded-g1-windows` before the first scan. This expands the
candidate worlds beyond #680's two diagnosis worlds. It selects replacement windows only; it does
not run G0's river matrix or any G1 treatment.

**Source and pool.** Pin source `5912a7d8be999e4148ce956be060932b77241614` (the fetched
`origin/main` at preregistration). The candidate manifest is
[`window-rescan-3-preregistration.json`](window-rescan-3-preregistration.json). It snapshots the
B2 catalog on 2026-11-13, includes every canonical `mantle-bloom-seed*.bloom` save at least 100
Myr old, chooses one save per seed (the greatest saved age; for equal-age duplicates, the name
without a numeric suffix), and excludes all saves from RFC 0015's held-out seeds 447276134,
628977350 and 295569319. Issue-specific checkpoint exports with the `seed<id>-...` naming pattern
are excluded because they are derived diagnosis windows, not independent catalog worlds. No
additional RFC 0016 acceptance save is named in the tracker or current RFC/report.
The two #680 worlds remain in the pool for continuity; additional seeds expand it. Candidate file
SHA-256 values are verified and recorded in the manifest before scanning.

**Primary scan.** Run `bin/debug/measure_lake_step_budget.py` on each listed save, in manifest
order, for 5 Myr at 100 kyr per step with the ice-age cycle disabled. Exclude the first 0.4 Myr
(four outputs) as startup. Keep the source fixed for the entire scan. Report per-world reversal and
comparable-body counts and candidate events using the unchanged #680/#588 detector: a matched
component above 10,000 km² with successive volume changes reversing sign and each change above 50%
of the prior volume.

**Candidate windows and screens.** As in #791, each candidate is a non-overlapping four-step
window starting two steps before a qualifying reversal; rank by descending `|ΔV|`. Use the
detector-matched component as target and the wet node nearest its centroid as anchor. Apply these
screens in order, without altering thresholds: (1) a source-world reload reproduces before, after
and `ΔV`; (2) frozen precipitation is at most 10% across the window over the catchment of each
land node draining to the target; (3) routed supply to the target is nonzero every step; (4) the
normalized liquid-ledger residual is at most 0.1% every step; and (5) a lake-ice-sublimation
knockout retains at least 50% of the control `|ΔV|`. Record every ranked candidate and its first
failing screen. Select at most three highest-ranked individual passers, with no overlapping
windows, disjoint target-node sets and anchors at least 1,000 km apart. Stop if fewer than three
pass; do not relax any screen.

**Step-size sensitivity.** Production remains 100 kyr steps. For every selected target, replay its
four-step 0.4 Myr interval from the same checkpoint at 10 kyr steps and record the same anchor's
lake-volume and connectivity trace. Treat this as a schedule sensitivity, not a replacement for the
production result or a timestep-invariance proof. Interpret the comparison alongside #853's lake
solver slice-length findings and #854's gate reporting requirement. No window is selected or
rejected based on this alternate schedule.

**Decision and artifacts.** Checkpoint hashes, anchor, target `|ΔV|`, frozen share, routed supply
for every step, maximum normalized liquid-ledger residual, knockout response and the 10 kyr
sensitivity go in compact JSON and this report. Save scans, replays and checkpoints only under the
private capture path recorded in the manifest, with `SHA256SUMS` beside them. If three independent
windows do not pass, document the evidence and recommend keeping G1b blocked or requesting an
owner-approved RFC revision. G0 reruns once at the final pin only after three valid windows are
selected.

**Scan progress (updated 2026-11-14).** All completed saves used the registered source and schedule, and
their save hashes matched the manifest. The first four outputs per save are excluded from the
reportable counts. Peaks include the whole scan.

| Seed, checkpoint | Reportable steps | Comparable observations | Reversals | Mean reportable step | Peak RSS | Largest reversal |
|---|---:|---:|---:|---:|---:|---|
| 910211954, 220.2 Myr | 46 | 3,254 | 74 | 45.3 s | 6,873.8 MB | 671,457 km³, step 44 |
| 317739385, 200.8 Myr | 46 | 828 | 67 | 292.5 s | 8,001.7 MB | 1.353 million km³ after 1.343 million km³ growth, steps 9–10 |
| 308161711, 121.5 Myr | 46 | 1,025 | 37 | 19.6 s | 1,206.0 MB | 313,023 km³, step 18 |
| 373936961, 244.4 Myr | 46 | 2,229 | 118 | 50.4 s | 6,422.1 MB | 2.210 million km³, step 37 |
| 374464662, 132.6 Myr | 46 | 611 | 15 | 14.8 s | 1,119.8 MB | 155,352 km³, step 16 |
| 456862899, 19,191.8 Myr | 46 | 1 | 0 | 13.4 s | 2,606.7 MB | none; one comparable observation at step 46 |
| 507012855, 479.0 Myr | 46 | 1,442 | 43 | 48.9 s | 1,236.8 MB | 2.328 million km³, step 49 |
| 544349906, 422.9 Myr | 46 | 1,203 | 33 | 7.5 s | 1,262.5 MB | 1.752 million km³, step 39 |
| 635056485, 452.8 Myr | 46 | 251 | 4 | 172.0 s | 3,867.1 MB | 140,707 km³, step 5 |
| 762517225, 135.2 Myr | 46 | 896 | 20 | 88.9 s | 1,179.6 MB | 606,992 km³, step 12 |
| 997271774, 130.5 Myr | 46 | 449 | 24 | 37.8 s | 1,748.5 MB | 297,049 km³, step 5 |

These are raw detector signals, not passing candidate windows. Full outputs remain in the private
capture, with hashes recorded beside the saves. Scans for seeds 308161711 and 373936961 emitted
repeated glacier maximum-depth warnings with zero meltwater shed. The 19.2 Gyr save showed no
detector reversals despite intermittent lakes reaching 53,216 km³ and several floor-drain
transitions. Seed 507012855 is a lake-rich world with roughly 5.6–11.3 million km³ of lake
volume in reportable steps; its 43 reversals are scan signals only. The 422.9 Myr save also
has large raw lake-volume reversals (33 across 1,203 comparable observations). All pool scans are
complete; candidate screening is in progress. Seed 635056485 repeatedly reached the 5,000 m glacier depth cap
at 0–15 nodes with zero meltwater shed; its 172.0 s mean reportable step made this scan much slower
than the others. It produced four raw reversals across 251 comparable observations. Seed 762517225
also emitted intermittent one-node glacier-depth warnings with zero meltwater shed; its scan found
20 raw reversals across 896 comparable observations.

### Candidate screen progress

Candidates are ranked by descending `|ΔV|` after greedily removing overlapping four-step windows
within each source world. The top three nonoverlapping signals come from distinct seeds:

| Rank | Seed, event step | `|ΔV|` (km³) | First failing screen / status |
|---:|---|---:|---|
| 1 | 507012855, 49 | 2,327,862 | Fails frozen precipitation: 58.08% across the window |
| 2 | 373936961, 37 | 2,209,627 | Fails reload: detector-matched prior volume is not reproduced |
| 3 | 544349906, 39 | 1,751,714 | Fails routed supply: target is absent at window step 4 |
| 4 | 544349906, 46 | 1,703,174 | Passes individual screens; knockout retains 99.91% |
| 5 | 544349906, 28 | 1,498,244 | Passes individual screens; knockout retains 99.999% |
| 6 | 317739385, 10 | 1,353,242 | Fails frozen precipitation: 83.67% across the window |

Candidate 2's checkpoint is at 247.8 Myr (SHA-256
`1be2433e81f66f36b76cdbb89b58335ee5e381305297476a42e6249cf277fbad`), anchor
`[-0.2391, -0.2595, 0.9357]`. The source scan's detector-matched target changed from 430.010 to
2,210,056.907 km³, but the reloaded event output has no matched prior body, so screen 1 fails.
Its after-volume reproduces. Diagnostic catchment values are also recorded: the four-step frozen
share is 67.12%, routed supply is positive each step, and the maximum normalized lake and land
residuals are `2.15e-13` and `7.11e-14`. These later checks do not alter the first-failure result;
the sublimation knockout is not run. The compact metrics are in
[`window-rescan-3-summary.json`](window-rescan-3-summary.json).

Candidate 3's checkpoint is at 426.5 Myr (SHA-256
`276b17601dbcdc23abd256bf27f260b7f6808f6ef0a310b2d70b6e1071393601`), anchor
`[0.9650, -0.2499, 0.0789]`. Reload reproduces the detector-matched 2,655,397.058 to
4,407,111.162 km³ event exactly. Frozen share is 1.65% weighted across the three wet-anchor
outputs; the anchor has no target and no routed supply on step 4. The maximum normalized lake and
land residuals are `7.51e-13` and `4.03e-13`. It fails screen 3, so the sublimation knockout is
not run. See the compact JSON for the per-step values.

Candidate 4's checkpoint is at 427.2 Myr (SHA-256
`82b64049c81c765b221fbc850500f3ae7ee578f65022a32a06387faa0fb148ce`), anchor
`[0.6434, -0.7536, 0.1345]`. Reload reproduces the 2,715,230.583 to 4,418,404.225 km³ event.
Frozen share is 6.59%; routed supply is positive in every step, and the maximum normalized lake
and land residuals are `4.41e-13` and `2.29e-13`. With lake-ice sublimation set to zero, the event
is 1,701,689.422 km³, retaining 99.91% of the control change. It passes all individual screens;
the independence rule applies when additional individual passers are available.

Candidate 5's checkpoint is at 425.4 Myr (SHA-256
`2f59c616610bcd754eab0306e780889f5729c7b09ac8ac0a9a236cdc04418a62`), anchor
`[-0.2405, -0.9630, 0.1217]`. Reload reproduces the 1,548,925.991 to 50,681.820 km³ event.
Frozen share is 1.64%; routed supply is positive in all four steps, and the maximum normalized
lake and land residuals are `1.05e-12` and `5.27e-13`. With sublimation disabled, the event is
−1,498,229.142 km³, retaining 99.999% of the control change. Its pairwise independence check
against candidate 4 passes: anchors are 6,008 km apart and the windows' target-node sets have no
overlap. These are two individual passers; a third is still needed.

Candidate 6's checkpoint is at 201.5 Myr (SHA-256
`5225d37f577110d131bc318396e49f80af93277b49aab67df2f1a2e9e94403fc`), anchor
`[-0.7078, -0.6923, 0.1403]`. Reload reproduces the 1,354,885.552 to 1,643.860 km³ event and
`ΔV = −1,353,241.692 km³`. Frozen share is 83.67%; routed supply is positive in each step, and
maximum normalized lake and land residuals are `6.72e-14` and `1.22e-14`. It fails screen 2, so no
sublimation knockout is run.

Candidate 1's four-step source reload at 100 kyr reproduces the event exactly: the target is
2,378,242.227 km³ before and 4,706,104.689 km³ after, `ΔV = +2,327,862.462 km³`. The weighted
frozen share is 58.08%; its four step shares are 58.06%, 57.92%, 58.24% and 58.09%. Routed target
supply is positive in all four steps. The largest normalized lake and land residuals are
`1.63e-11` and `7.89e-12`, respectively. It fails screen 2, so the sublimation knockout and
independence screens are not run for it. Its checkpoint hash, anchor, per-step values and raw
screen output are recorded in `window-rescan-3-summary.json`; raw files and saves remain in the
private capture.

Candidate 7 (seed 507012855, event step 10) reloads from the checkpoint at 479.7 Myr. The four-step
control reproduces the 780,981.413 to 100,647.589 km³ event (`ΔV = −680,333.823 km³`). Frozen
precipitation exceeds the 10% screen (weighted frozen-to-liquid share 91.35%; per-step shares
83.17%, 82.35%, 146.72% and 133.29%). Routed supply is positive in all four steps; maximum
normalized lake and land residuals are `1.86e-13` and `8.15e-14`. It fails screen 2, so no
sublimation knockout is run.

Candidate 8 (seed 373936961, event step 50) reloads from the after-step-47 checkpoint at 249.1 Myr
(SHA-256 `90a4ea1c070b25ba908251704fb4711f41c91cce5577d9ed04a13499f42b4e81`). The four-step
control reproduces the 681,930.463 to 2,016.227 km³ event (`ΔV = −679,914.236 km³`). Its weighted
frozen-to-liquid precipitation ratio is 25.07 (2,506.58%); all four steps exceed the 10% screen.
Routed supply is positive throughout, and maximum normalized lake and land residuals are
`2.53e-13` and `8.11e-14`. It fails screen 2, so no sublimation knockout is run.

Candidate 9 (seed 910211954, event step 44) reloads from the after-step-41 checkpoint at 224.5 Myr
(SHA-256 `ed562d2ffdfd43d60ca8dca0da880c3f7dfb12c2282c221d1dc1424267a4b9e1`). The four-step
control reproduces the 673,721.222 to 2,264.690 km³ event (`ΔV = −671,456.532 km³`). Its weighted
frozen-to-liquid precipitation ratio is 4.615 (461.50%); all four steps exceed the 10% screen.
Routed supply is positive throughout, and maximum normalized lake and land residuals are
`2.36e-12` and `5.44e-13`. It fails screen 2, so no sublimation knockout is run.

Candidate 10 (seed 373936961, event step 43) reloads from the after-step-40 checkpoint at 248.4 Myr
(SHA-256 `00a6aacf10b8639ffb0b80025701a1fed3509a84e3dfd50166861b18dadf7f91`). The four-step
control reproduces the 670,740.039 to 3,071.312 km³ event (`ΔV = −667,668.726 km³`). Its weighted
frozen-to-liquid precipitation ratio is 3.091 (309.13%); every step exceeds the 10% screen. Routed
target supply is zero in the first three steps. Maximum normalized lake and land residuals are
`1.09e-12` and `3.56e-13`. It fails screen 2, so no sublimation knockout is run.

Candidate 11 (seed 544349906, event step 16) uses the existing after-step-13 checkpoint at 424.2
Myr. The control reproduces the scan's 663,337.060 km³ after-volume, but the anchor's prior body
does not match in the replay (the detector's 488.990 km³ before-volume cannot be reproduced), so
the source-reload screen fails. Its weighted frozen-to-liquid precipitation ratio is 3.01%; routed
supply is positive in all four steps, and maximum normalized lake and land residuals are
`5.87e-14` and `3.16e-14`. Later screens, including the sublimation knockout, are not run.

Candidate 12 (seed 544349906, event step 20) uses the existing after-step-17 checkpoint at 424.6
Myr. The control reproduces the scan's 662,777.491 km³ after-volume, but the detector's 683.698 km³
before-volume is not matched in the replay, so the source-reload screen fails. Its weighted
frozen-to-liquid precipitation ratio is 2.97%; routed supply is positive in all four steps, and
maximum normalized lake and land residuals are `2.77e-13` and `1.38e-13`. Later screens, including
the sublimation knockout, are not run.

Candidate 13 (seed 910211954, event step 30) uses the after-step-27 checkpoint at 223.1 Myr. The
control reproduces the scan's 657,416.858 km³ after-volume, but the detector's 1,177.711 km³ prior
volume is not matched in the replay, so the source-reload screen fails. Frozen share is also above
the 10% limit (weighted frozen-to-liquid ratio 1.332); routed supply is positive in all four steps,
and maximum normalized lake and land residuals are `1.43e-12` and `3.37e-13`. Later screens are not
run.

Candidate 14 (seed 507012855, event step 17) uses the after-step-14 checkpoint at 480.4 Myr. The
control reproduces the scan's 652,630.120 km³ after-volume, but the detector's 343.781 km³ prior
lake is not matched in the replay, so the source-reload screen fails. Frozen share is 0.71%; routed
supply is positive throughout, and maximum normalized lake and land residuals are `7.22e-14` and
`3.35e-14`. Later screens, including the sublimation knockout, are not run.


Candidate 15 (seed 762517225, event step 12) uses the after-step-9 checkpoint at 136.1 Myr. The
control reproduces the 607,472.077 to 480.506 km³ event (`ΔV = −606,991.572 km³`). Frozen share
exceeds the 10% limit (weighted frozen-to-liquid ratio 12.73%; all four steps exceed 10%). Routed
supply is positive throughout, and maximum normalized lake and land residuals are `8.79e-14` and
`3.20e-14`. It fails screen 2, so no sublimation knockout is run.

Candidate 16 (seed 544349906, event step 10) uses the after-step-7 checkpoint at 423.6 Myr (SHA-256
`5e93f0feacf0b54d4016321bade8b0d36f435760ea55c1583a1678c94f2e9b41`). The control exactly
reproduces the 603,549.454 to 410.769 km³ event (`ΔV = −603,138.685 km³`). Weighted frozen share
is 9.15%; routed supply is positive in all four steps, and maximum normalized lake and land
residuals are `6.64e-13` and `3.60e-13`. The sublimation knockout gives `ΔV = −603,079.817 km³`,
retaining 99.99% of the control change.

Candidate 16 is independent of candidate 4: their anchors are 3,514 km apart and their combined
four-step target-node sets do not overlap. It is also independent of candidate 5: their anchors are
2,499 km apart and their target-node sets do not overlap. The three selected windows are ranks 4,
5 and 16.

At 10 kyr resolution, candidate 4's anchor stayed wet at every 0.1 Myr endpoint. Its target-lake
volumes were 3.949, 3.416, 3.698 and 3.609 million km³, compared with 4.528, 2.715, 4.418 and
4.534 million km³ at 100 kyr. Candidate 5 also stayed wet at all endpoints; its target component
had 12 nodes through 0.3 Myr and 10 nodes at 0.4 Myr. Its final endpoint volume was 20,408 km³,
compared with 51,101 km³ at 100 kyr.

At 10 kyr resolution, candidate 16's anchor remained wet at all four endpoints. Its target-lake
volumes were 382, 427, 502 and 583 km³, with 9–10 target nodes, compared with 383, 439, 411 and
380 km³ at 100 kyr. These alternate-cadence results describe volume and connectivity sensitivity;
they do not affect candidate selection.
