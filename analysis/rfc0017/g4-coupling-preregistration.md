# RFC 0019 G4: downstream coupling preregistration

**Status: preregistered 2026-11-18, before G4 comparisons.** G3 is formally failed. The
owner-approved exception allows this downstream impact gate to proceed; it does not mark G3
passed. The Round 2 setting `alpha=1, L=1,500 km` is the best safeguard-eligible result and is
used here only as a fixed diagnostic arm. It is not selected for production. No G3 parameter is
retuned from this gate.

**Accounting clarification (2026-11-18, before the paired comparisons were complete):** the
replay exposes `land_sink_settlement` as a named, accounted-for land-water destination where
routed external supply reaches a land sink without an active lake body. It is not stranded or
unlocated water. Report its amount; do not require it to be zero. The separate
`stranded_overflow` field remains required to be zero.

## Question

When the best eligible conservative-moisture diagnostic replaces the current precipitation
field, does RFC 0016's central land-water coupling still close its ledgers, and how does it
change lake supply and volume, ocean export, frozen precipitation, glaciers, erosion, and
biomes? Separate the atmosphere change from the land-store treatment with a paired 2×2 screen.

This gate measures downstream consequences of the current single-column approximation. It
does not establish that the G3 inland/coast pattern passes, validate the humidity proxies as
measured evaporation, or choose a production implementation.

## Registered inputs

Use the three current-source windows selected by RFC 0016's expanded rescan (#871), not its old
A/D/E windows. Those former events do not reproduce on the current source. The selected
windows pass the registered reload, frozen-precipitation, routed-supply, ledger, sublimation
knockout, and basin-independence screens. They are three disjoint target basins in one seed, so
they are not three independent worlds.

All checkpoints come from `mantle-bloom-seed544349906-422900000y.bloom` (SHA-256
`10a127f99d79cf1d2d2272661ea5410613a0aaa6d5c3da6d6639411991a4fbde`). Verify each checkpoint
and anchor file against the private `rfc0016/issue871-window-rescan-3/SHA256SUMS` manifest
before loading. Checkpoints and anchors are stored outside the repository.

| RFC 0016 rank | Checkpoint age | Checkpoint | Checkpoint SHA-256 | Event anchor |
|---:|---:|---|---|---|
| 4 | 427.2 Myr | `checkpoints/seed544349906/after-step-043.bloom` | `82b64049c81c765b221fbc850500f3ae7ee578f65022a32a06387faa0fb148ce` | `anchors/seed544349906-step46.json` |
| 5 | 425.4 Myr | `checkpoints/seed544349906/after-step-025.bloom` | `2f59c616610bcd754eab0306e780889f5729c7b09ac8ac0a9a236cdc04418a62` | `anchors/seed544349906-step28.json` |
| 16 | 423.6 Myr | `checkpoints/seed544349906/after-step-007.bloom` | `5e93f0feacf0b54d4016321bade8b0d36f435760ea55c1583a1678c94f2e9b41` | `anchors/seed544349906-step10.json` |

The registered 100 kyr windows each cover four steps. Repeat each from the same checkpoint at
10 kyr for 40 steps as timestep sensitivity; do not use the alternate schedule to select
parameters or windows. The saves have climate density 4 and diagnostic wind selected. Preserve
those settings in every arm. This is a density-4, diagnostic-wind downstream screen; it does not
replace G3's broader wind and density checks.

## Fixed atmosphere treatment

The candidate arm runs the offline conservative finite-volume solver on the current world's
fields at every climate update. It uses the G3 source proxies, effective column depth, density,
rainout baseline, temperature-derived cap, and native placement code as registered in the
three-round G3 plan. Use coast-focused sink `lambda = lambda0 * (1 + exp(-d/1,500 km))` on land,
`lambda = lambda0` over ocean, with `alpha=1`; retain the G3 ridge-interception strength
`p=0.25`. Use the same saved diagnostic wind mode. The treatment produces precipitation and
column humidity for the ordinary downstream climate, hydrology, ice, erosion, and biome
consumers. It introduces no persisted atmospheric state.

The baseline atmosphere arm uses the current production climate unchanged. The G3 solver's
source mapping remains a proxy: its lake/river and vegetation inputs are not measurements of
actual evaporation or transpiration. This known limitation is part of the interpretation,
not a parameter to tune in G4.

## Paired arms

For each window and step schedule, run these four arms from byte-identical checkpoint inputs,
one process at a time:

| Arm | Atmosphere | RFC 0016 central land-water store |
|---|---|---|
| `legacy-control` | Current production climate | Off |
| `legacy-coupled` | Current production climate | On, default central settings (`C=300 mm`, `p=0.5`, `D=1`, zero split shift, humidity source on) |
| `conservative-control` | Fixed `alpha=1, L=1,500 km, p=0.25` diagnostic | Off |
| `conservative-coupled` | Same fixed diagnostic | On, same central settings |

This separates the candidate atmosphere's downstream shift, the coupling treatment's shift,
and their interaction. All arms use the same saved climate density, wind mode, world state,
ice-age-cycle setting, and deterministic seed/time inputs. No sweep or fallback parameter is
allowed.

## Readouts and validity checks

Record every step, plus the opening state:

- Atmosphere: source, rain, overflow, global and summed-cell residuals, finite/nonnegative
  column water and rainfall, and climate-update time. Require closure at `1e-8` relative to
  total source; do not repair residuals by adding rain.
- Surface water: full land and lake ledgers, runoff, evapotranspiration, lake supply and
  volume, river inventory, and ocean export. Require the existing RFC 0016 normalized land,
  lake, and river ledger residuals to stay at or below `1e-3` each step; stranded overflow must
  be zero. Report named `land_sink_settlement` separately as an accounted-for diagnostic; it
  may be nonzero.
- Frozen water and downstream effects: frozen precipitation volume/share, glacier-water
  ledger and ice volume, rain and river erosion totals, and area-weighted biome shares. Report
  absolute values as well as treatment/control and atmosphere differences.
- Runtime and peak RSS by arm and schedule. Timings are descriptive and only compared between
  arms run without other compute jobs.

Report paired changes for each window separately and as a descriptive median across windows.
Do not call three same-seed windows three independent replications. Do not infer that lake or
glacier effects are caused by plant transpiration alone: atmosphere and surface treatment
interact, and the current G3 source mapping remains approximate.

## Decision rules

- A failure of atmospheric closure, positivity, a required surface ledger, or zero stranded
  overflow makes that arm invalid and blocks G4 completion until the accounting defect is
  resolved.
- If all ledgers close, G4 is an impact report: quantify downstream changes and decide whether
  they are acceptable for a visually oriented, synthetic-world model. No minimum lake, glacier,
  erosion, or biome effect is imposed by this preregistration.
- If the coupled conservative arm has unacceptable downstream changes, keep G3 failed and
  reopen the implementation choice before any production integration. Do not tune this G4
  screen to force a desired downstream response.

## Reproduction and outputs

Implement a debug-only replay wrapper that patches `climate.compute_climate` for the fixed
candidate arm and reuses RFC 0016's existing replay snapshots and treatment. Preserve the
current production code path for the legacy arms. The wrapper must record the loaded checkpoint
SHA-256, source commit, worktree state, parameter fingerprint, and every script hash used.

Store replay JSON, logs, and `SHA256SUMS` under
`~/mantle-bloom-captures/rfc0019/g4/`. Check in the small G4 runner and compact summary only.
Do not check in checkpoints, full replay logs, or per-node arrays. Record the actual outcomes
in `analysis/rfc0017/README.md` and RFC 0019's G4 result entry after all preregistered arms
finish.
