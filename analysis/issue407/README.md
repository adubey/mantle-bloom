# RFC 0007 gate: source attribution on the stripped baseline (#407)

Diagnosis for [RFC 0007](../rfc/0007-mantle-plumes-and-plateaus.md)'s source-attribution
gate, shared with [#406](https://github.com/adubey/mantle-bloom/issues/406). Status:
measurement complete 2026-10-25; the RFC's disjoint-ownership definitions follow in
PR #574. No plume physics was implemented here.

## Reproduce

Run from the repository root after `./bin/setup.sh`. One world at a time
(replays are RAM-heavy, ~65 MB saves, ~5-6 min per 500 Myr seed at density 1):

```sh
backend/.venv/bin/python bin/debug/attribute_reset_sources.py \
  --out analysis/issue407/attribution.json --save-dir /tmp/mantle-bloom-issue407
backend/.venv/bin/python bin/debug/attribute_reset_sources.py --probe-only \
  --out /tmp/probe.json --save-dir /tmp/probe-saves
backend/.venv/bin/python bin/debug/attribute_reset_sources.py --self-test-only \
  --out /tmp/selftest.json --save-dir /tmp/selftest-saves
```

The script refuses to run if `backend/app` or `backend/constraints.txt` differs from
HEAD, so the JSON always describes an unmodified model. Every tap wraps a function,
records, delegates, and is restored afterwards; taps never mutate world state, so a
tapped run steps identically to an untapped one.

## Stripping procedure and controls (shared with #406)

Default `world.generate_world(seed)`, cratons on, every other argument default, no
control overrides — the same inputs as `compare_cratons.run_job(seed, "on", ...)` in
[#403](issue403/README.md), cut to 500 Myr. 5 Myr steps, 100 Myr checkpoints, seeds 7,
11, 2024, serially. Model baseline: `8ed179f`; Python 3.14.6; constraints SHA in
`attribution.json`. Final-world saves stay outside Git (`/tmp/mantle-bloom-issue407`,
SHA recorded in the JSON).

Two taps attribute every melt reset and every continental-ledger booking:

- `continental_ledger.record` wrapper: every booking tagged by caller
  (`module.func:line`) and account. Exact volumes.
- `_erupt_melted_nodes` wrapper (both the `lithosphere_plate` and the
  `quad_tectonics` references): snapshots prior crust code, prior elevation and Hc of
  every melting node, then classifies oceanic→continental recodes and the reset-created
  thickness (`hc_after − hc_before`) the caller books as juvenile on those nodes.
  Volumes use per-node areas on full-plate (`deform_columns`) calls — exact — and the
  plate mean area on donor/cell-subset (`_open_rift`) calls — estimated with the
  max-min spread reported. Tap arithmetic is self-tested (`--self-test-only`: synthetic
  full-plate and subset calls with one exposed recode each; both report exactly
  31,000 m × 1e12 m²).

## Findings

**1. The reset route never fires on the stripped baseline.** Across 3 seeds × 500 Myr
(12,500+ erupt calls, 53,000+ melting nodes), zero melting nodes were oceanic-coded
with prior elevation above sea level and recoded continental, and the melt-site
juvenile bookings that use the post-recode classification (`deform_columns:909`,
`_open_rift:1844`) booked exactly 0 m³ on every seed. Every melt-through on current
main is a sub-sea-level oceanic eruption.

**2. The mechanism exists in code.** The synthetic exposed oceanic column
(`--probe-only`, drives the real `_erupt_melted_nodes`): an oceanic-coded node at
+400 m recodes to continental with Hc 2.5 km → 35 km, and a caller books 32.5 km of
juvenile thickness on it; the drowned (−3000 m) control stays oceanic (→ 11 km, no
juvenile). So the code path is live but unreachable on current main — melt-through
only happens at or below sea level.

**3. Juvenile volumes by caller, gross 0–500 Myr (M km³):**

| Caller | Seed 7 | Seed 11 | Seed 2024 |
|---|---|---|---|
| arc magmatism (`deform_columns:795`) | 626.2 | 636.0 | 676.0 |
| rifted cells (`_open_rift:1891`) | 77.9 | 63.4 | 71.4 |
| magma returns (`magma_transport:334`, `remelted_relaminated_returns_m3`) | 33.3 | 31.4 | 26.6 |
| rift magmatism (`deform_columns:865`) | 20.2 | 29.1 | 16.3 |
| arc-margin seeds (`grow_frontier:1602`) | 5.4 | 9.6 | 11.0 |
| volcano eruptions (`volcanism:79`) | 6.2 | 4.2 | 5.1 |
| melt sites, post-recode (`deform_columns:909`, `_open_rift:1844`) | 0 | 0 | 0 |

Gross `juvenile_additions_m3` totals 736/742/780 M km³ per seed; the ledger closes
each checkpoint with balance error under ~1% of cumulative gross bookings (seed 7:
−6.9 M km³ at 500 Myr). Gross arc production dominates; net continental area is
roughly steady (seed 7: 236→199 M km²; seed 11: 151→145; seed 2024: 123→138) —
no regrowth attributable to resets, which contributed nothing. This baseline serves
#406's stripped-world question with the same numbers.

**4. Breakup background for the plume-stress comparator.** Splits per 500 Myr:
26 (seed 7), 28 (seed 11), 36 (seed 2024); failed rifts 25/24/26 (cumulative event
counts per checkpoint are in the JSON). No plume schedule exists yet, so no
conditional incidence can be computed; when plume code lands, the comparator is:
replay the same scheduled head times and host plates plume-off under identical
controls and step schedule, and report the head-following 20 Myr split fraction
beside the background fraction. The provisional half-of-heads target stays a design
choice, not evidence.

## Reporting standard for plume-code tests (disjoint ownership)

When plume physics is implemented, its tests must report three disjoint quantities,
never only their sum:

- direct plume juvenile: `plume_head`/`plume_tail` emplacement volumes;
- docked plateau conversion: plateau-sourced `converted_m3`;
- reset-route juvenile: melt-reset bookings at recoded nodes, as captured here.

Aggregate `juvenile_additions_m3` growth does not establish plume emplacement,
plateau docking, or preservation. Plateau `converted_m3` is the only figure
comparable to the ~1.1 km³/yr Cenozoic accretion budget.

## What this does not establish

- No oceanic plateau, flood basalt, hotspot track, or weakening exists yet; gross
  emplacement, docking, and conversion are all unmeasured by construction.
- The zero-recode result is for the stripped baseline (default generation, cratons
  on, 500 Myr). Longer runs, other steps, or other densities could in principle find
  an exposed melt-through; the taps above detect it exactly if it happens.
- Cratons were on: the result does not separate craton-suppressed continental
  melt-through from its absence.
