# Issue #867: targeted river substeps for persistent valley incision

**Status:** experiment set up, technically verified and peer reviewed (PR #868); the
reference run and the branches have not been run. This page is the preregistration, the technical verification of the three
arms, and the runbook for the agent that runs them. No model change is proposed.

## The question

Lakes and inland seas in Mantle-Bloom are far larger than Earth's. One hypothesis, from
design principle 3 (model fast processes in aggregate over long steps), is a timescale
mismatch: a 100 kyr tectonic step reorganizes drainage, and rivers only erode along the
routing left after the change, so an outlet that incision could have kept open is lost. The
#866 diagnostic measured, over 46,657 outlet changes on seed 544349906, a median required cut
of 81 m on the old outlet path (53.7% at most 100 m) against a median modeled channel growth
of 0 m on that path: inconclusive under its rule, so a causal test was asked for.

The owner's plan ([issue comment](https://github.com/adubey/mantle-bloom/issues/867)) is a
reference trajectory to 25 Myr with loadable checkpoints, then an assumption ladder of three
arms branched from one checkpoint under identical tectonics, from the most permissive bound to
a simulation-driven candidate.

## The arms

All three are module knobs in `backend/app/hydroclimate/river_substeps.py`, off by default
(`river_substeps.ARM`), set per process by the driver. With every knob off the model is
unchanged; the unit tests pin that a one-slice march is bit-identical to the baseline.

| Arm | Spelling | What it assumes | Where it acts |
|---|---|---|---|
| Baseline | `off` | the production model | — |
| Outlet-preservation bound | `outlet-bound` | a *threatened outlet* (a pit whose water reached the ocean last step) may cut its least-climb route at 1 m/yr for the step, stopping at the cut the route needs | `breaching.OutletRelief`: the per-node carving limits and the summed-cost screen are raised for eligible pits; the route stays the ordinary least-climb path; the water gate still applies. The allowance is per node on the route, and at 100 kyr steps it is 100 km, so the binding constraint is the per-node notch ceiling (`BREACH_MAX_NOTCH_M`, 2000 m), not the rate: in effect every threatened pit with a finite-cost route and a passing water gate drains |
| Outlet bound, ceiling | `outlet-bound:all` | every land pit is eligible | same; bounds what *any* incision rule could achieve |
| Uplift-limited bound | `uplift-bound` | the same trigger, but each node's allowance is its rise since the last erosion pass (`channel_reference_elevation_m`): a river keeping pace with uplift | same mechanism, per-node allowance |
| Simulation-driven march | `march:<n>` | the step's river incision and routing run in `n` slices; each slice re-breaches (with the time left in the step as its budget), rebuilds the depression hierarchy and flow directions on the terrain the earlier slices cut, re-routes the step's steady sources, and incises for its share of the step | `erosion.apply_erosion` takes the summed demand through its ordinary caps; `hydro.flow_target`/`discharge_m3_s` become the last slice's |

**Threatened outlet.** Last step's `flow_target` chains are followed to their end
(`river_substeps.exported_to_ocean`); a node whose chain ends in the ocean exported water.
The mask is resampled onto this step's nodes by nearest node, the way the model resamples
its own previous ocean mask. A pit on an exporting node is a threatened outlet: last step
the water here got out. The bootstrap routing at generation counts as a previous step.

**What the march holds fixed.** Climate, the lake solve (levels and spills) and the ice
solve are computed once on the opening routing, as today; the water ledger closes on that
routing. A basin the march opens mid-step keeps its opening water until the next step's
solve drains it (the one-step lag a breach already has). Re-solving lakes per slice is the
broad step-integrated landscape solve of #297, out of this issue's scope.

**Two choices inside the march** (from the PR #868 review) that the reading must allow for.
Each later slice re-breaches with the *time left in the step* as its budget, not its own
slice length, so the `march:100` versus `march:10` comparison is partly convergence of
that choice and not only of routing frequency. And the step's drop cap is widened to the
march's own slice-by-slice cut for the whole step's removal, so rain, weathering, glacier
and seismic removal in a deepened valley can use headroom only the river slices earned
(conservative in practice: each slice's cut was drop-capped on its own terrain).

**Rock and water accounting.** The bounds only change which pits breach; notches are
sub-cell relief as before. The march's cut is a demand that goes through the same
cover-first, headroom-capped, craton-resisted, isostatically compensated write-back as the
baseline's river term, so `removed == deposited + terminal` holds as before (pinned by
`test_the_march_runs_logs_and_keeps_the_ledgers_closed`). The step's drop cap rises to the
march's own slice-by-slice cut where the valley deepened.

## Preregistration

Written before the reference run.

**Input.** Seed 544349906 at node density 4.0 (the #866 world, recreated: #866's saves are
gone), generation defaults, glacier march off as in #866 (`--glacier-march off`), diagnostic
wind, 100 kyr steps. Record the commit, controls and node count (the driver writes them to
`summary.json`).

**Schedule.**

1. Reference: 25 Myr in 250 steps; checkpoint (`.bloom` save, combined and elevation maps)
   at 0, 5, 10, 15, 20 and 25 Myr. Inspect the maps and the Lake Inspector before branching.
2. Branches from the 25 Myr checkpoint (or an earlier one if inspection says the inland
   water is already in place, see below), 5 Myr each in 100 kyr steps, checkpoint every
   1 Myr: `off`, `outlet-bound`, `outlet-bound:all`, `uplift-bound`, `march:10` (10 kyr
   slices). Then the convergence check: `march:100` (1 kyr slices) on the same checkpoint.
   One arm at a time (the compute-budget rule); the driver runs each arm in its own process.

**Primary measures**, per step and compared at matched ages: visible lake area and liquid
volume, sea-tier area and volume, the number of bodies over 10,000 km², endorheic roots and
their area; the exporting land area and old-outlet retention (the share of last step's
exporting area that still exports); channel depth (land mean, p99, max), notch sum, the
arm's own log (eligible, breached and relieved pits; march slices, cost and rerouted nodes).

**Validity.** Every arm must keep the land-water ledger residual at or below 1e-9 of rain
plus meltwater, the lake ledger residual at round-off, and `removed == deposited + terminal`
in the erosion budget; a step that violates one is invalid for inference. The `off` branch's
`state_hash` at its first step must equal the reference's hash at the same age, which shows
the branch is a continuation and not a perturbation (verified on the density-4 probe below).

**Reading the ladder.**

- If `outlet-bound:all` does not reduce lake volume and area materially against `off` over
  5 Myr, incision cannot be the lever, whatever its rate: look at the water balance,
  evaporation, hypsometry and startup water (#856) instead.
- If `outlet-bound:all` helps but `outlet-bound` does not, the lost outlets are not the
  mechanism; basins that never had an outlet are.
- An `outlet-bound` effect is the effect of letting threatened outlets drain, not of a
  1 m/yr rate (see the arms table): do not read it as a rate result. `uplift-bound` is the
  rate-limited arm.
- If `outlet-bound` helps and `uplift-bound` recovers most of it, a rate-limited rule tied to
  uplift is worth an implementation issue.
- If the bounds help and `march:10` recovers a meaningful share, with `march:100` close to
  `march:10`, the simulation-driven route is viable and converging; if `march:100` differs a
  lot from `march:10`, report the interval sensitivity and do not claim convergence (the
  remaining-time breach budget above is one reason the two can differ).
- "Materially" is preregistered as a 10% change in whole-world liquid lake volume at
  matched age, with the direction stated; smaller effects are reported as such.

## Technical verification (done)

Probe on `generate_world(544349906, node_density=4.0)`, 127,958 nodes, glacier march off,
100 kyr steps, this branch, Apple Silicon (16 GB):

| What | Measured |
|---|---|
| Generation | 6.5 s, 0.66 GB peak |
| One 100 kyr step (baseline) | 4.9–5.1 s, about 1.2 GB peak RSS |
| Default-profile save / exact save | 38 MB (1.4 s) / 155 MB (0.09 s) |
| Load default save, step once, compare with the live world's next step | identical `state_hash` (exact save: identical too) |
| `march:4` / `march:10` / `march:25` | +0.7 s / +2.0 s / +5.3 s per step (about 0.21 s per slice) |
| `outlet-bound` / `uplift-bound` | no measurable cost |
| Land-water residual, every arm | at most 8e-13 of rain plus meltwater |

So a 25 Myr reference is about 25 minutes; a 5 Myr branch is 4–5 minutes at `off`,
about 7 at `march:10`, about 22 at `march:100`, and `march:1000` (100 yr slices) would be
about 3 hours. Ten-year slices are out of reach on this path. Each arm runs in its own
process (hydrology's per-process caches otherwise keep the previous arm's world alive). The
first map render is the process's memory peak (about 2.3 GB transient at 1600×800 on a tiny
world); the driver records the peak after each render so it can be told from the step's own.

One-step effects on the probe, from the same 0.3 Myr state (descriptive only; the arms
act on next step's lakes, so a one-step look understates the march and the bounds):

| Arm | Pits | Eligible | Breached | Relieved | Lake volume after the step |
|---|---|---|---|---|---|
| `off` | 555 | — | 29 | — | 2,097,993 km³ |
| `outlet-bound` | 555 | 33 | 43 | 14 | 2,007,124 km³ |
| `outlet-bound:all` | 555 | 555 | 119 | 90 | 1,600,972 km³ |
| `uplift-bound` | 555 | 33 | 29 | 0 | 2,097,993 km³ (identical to `off`) |
| `march:10` | — | — | 74 breached in later slices, 162 nodes rerouted | unchanged this step (lag) |

Two things worth the runner's attention. First, whole-world lake volume is already about
2.1 million km³ at 0.3 Myr on this seed (Earth's lakes hold about 0.18 million km³, the
Caspian 78,000 km³): the excess is present almost from generation, which #856 (startup
water and outlet changes) also points at. The 25 Myr reference will show whether it grows,
holds or drains; if it is in place at 5 Myr, branch from the 5 Myr checkpoint as well.
Second, at this age the uplift bound changes nothing: the rise since the last erosion pass
on the route nodes is below what the closed pits need. That is the arm doing what it says;
it may matter more once mountains are growing.

## Suggested changes to the design

- Add `outlet-bound:all` as a ceiling arm (done above). It costs nothing and answers the
  question that matters most for the underlying problem: can incision of any kind remove
  the inland water at all?
- Candidate intervals: 10 kyr and 1 kyr slices (`march:10`, `march:100`) are affordable and
  give the convergence check the issue asks for; 100 yr is a 3 hour arm; 10 yr is not
  reachable without a local (valley-only) solver, which would be its own issue.
- On the owner's "falling rock, little deposition" point for mountain rivers: the model
  already deposits only where a river is big and slow (`river_speed` under 2 with
  `water_accum_m` over 0.05, `DEPOSITION_FRACTION` 0.15), so steep reaches do not deposit
  today. No change is needed for the arms; the march keeps that rule.
- Keep the glacier march off in all arms, as in #866, so the comparison is about rivers;
  note it as a limitation (glaciers carve channels too).

## Runbook for the executing agent

From the issue worktree, with the compute-budget skill's check first (`--run sim --gb 4`):

```bash
backend/.venv/bin/python bin/debug/river_substep_arms.py reference \
    --seed 544349906 --node-density 4 --myr 25 --checkpoint-myr 5 --glacier-march off \
    --out-dir ~/mantle-bloom-captures/issue867/reference
```

Inspect `reference-*-combined.png` and load a checkpoint in the app (Lake Inspector). Then:

```bash
backend/.venv/bin/python bin/debug/river_substep_arms.py branch \
    ~/mantle-bloom-captures/issue867/reference/reference-25.0Myr.bloom \
    --arm off --arm outlet-bound --arm outlet-bound:all --arm uplift-bound --arm march:10 \
    --myr 5 --checkpoint-myr 1 --glacier-march off \
    --out-dir ~/mantle-bloom-captures/issue867/branch-25
backend/.venv/bin/python bin/debug/river_substep_arms.py branch \
    ~/mantle-bloom-captures/issue867/reference/reference-25.0Myr.bloom \
    --arm march:100 --myr 5 --checkpoint-myr 1 --glacier-march off \
    --out-dir ~/mantle-bloom-captures/issue867/branch-25
```

Outputs: `<label>-steps.json` (one row per step), `<label>-<age>Myr.bloom` and the two maps
per checkpoint, `SHA256SUMS`, and `summary.json`. Check `off-steps.json[0].state_hash`
against the reference row at the same age. Commit only the compact JSON and this page's
findings; saves and images stay in `~/mantle-bloom-captures/issue867/` (use the `b2-saves`
skill to move them). Post the preregistered measures as a table in the issue, per arm at
matched ages, with the ledger residuals and the runtime and peak RSS of each arm, and end
with the recommendation the ladder supports: a small implementation issue, an RFC, or a
redirection to the water balance.
