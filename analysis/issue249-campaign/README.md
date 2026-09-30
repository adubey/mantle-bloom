# Issue #249 parity campaign: reference results

Line vs sparse-quad parity runs with the #247 harness (`bin/debug/surface_parity.py`;
gates and tolerances in `docs/surface-parity.md`). They were run on `issue-249-parity-campaign`
at harness commit `f345cc5`. The `repro-lines` set was run at `683247a`, which only fixed a
false positive in the hydrology finite check. Each set's `provenance.json` records the exact
command, commit and platform.

These are the reference numbers a later campaign, for example after fixing #253–#257, should
be compared against.

## Sets

All commands run from `backend/`, with `--out ../analysis/<campaign>/<set>`.

| set | what | command |
|---|---|---|
| `long` | seeds 1–5, density 1, 1 Myr steps, to 400 Myr | `paired --preset long --seeds 1,2,3,4,5 --jobs 5` |
| `repro-equivalent` | the 352.4 Myr reproducer's seed on both surfaces, density 1, to 352 Myr, with renders | `paired --preset long --checkpoints-myr 30,120,240,352 --seeds 804913535 --jobs 2 --render` |
| `stress` | seeds 1–2, density 0.5, to 800 Myr | `paired --preset long --node-density 0.5 --checkpoints-myr 100,200,400,600,800 --audit-every 10 --samples 100000 --seeds 1,2 --jobs 2` |
| `repro-lines` | the real 352.4 Myr line save (`mantle-bloom-seed804913535-352400000y.mbworld`, not in the repo), continued 1 Myr under every audit | `run --preset issue147 --seed 804913535 --surface lines --from-world <save> --step-years 100000 --checkpoints-myr 0.5,1 --audit-every 1 --render` |

Each set directory keeps:

- `seed*-*.json`: deterministic metrics, byte-identical on a rerun of the same code;
- `comparison.json` and `report.md`: gate results;
- `phase-means.json`: per-phase mean seconds per step, taken from the uncommitted per-step
  timings;
- `provenance.json`.

The per-step `*.timings.json` files and logs are not committed. Those runs shared the machine,
so the timings only mean something as quad/line ratios within a set.

Diagnostics behind the focused issues:

| file | script | used in |
|---|---|---|
| `continental-budget-seed1*.json` | `bin/debug/attribute_continental_budget.py` | #253 |
| `quad-retreat-continental-seed1.json` | `bin/debug/probe_quad_retreat_continental.py` | #253 |
| `equal-area-bias-seed7.md` | `bin/debug/measure_equal_area_bias.py` | #257 |

## Comparing a new campaign

1. Rerun the commands above into a new `analysis/<campaign>/` directory.
2. Run `bin/debug/summarize_parity_campaign.py ../analysis/<campaign>`.
3. Diff its `summary.md` against this directory's. Both are re-judged with the gate code
   current at summarize time, so rerun the summary here too if the gates have changed since.

For a single set, `bin/debug/surface_parity.py compare <dir>` re-judges its runs.

## Findings

Full tables are in `summary.md`.

- **Hard invariants (H1–H10) hold in every set.** No folded quads, invalid neighbours, stale
  derived caches or revision drift. Save/load and continuation from a loaded world were
  identical throughout, including on the 158k-node reproducer.
- **The #228 artifacts are gone on quad.** Zero stacked nodes, against 36% on lines at
  400 Myr; anisotropy near zero with no row alignment; no row streaks in the 352 Myr renders.
  Uncovered sphere is 1.3% on quad against 8.0% on lines.
- **Blockers for the #250 cutover:**
  - #253: quad boundary retreat destroys continental crust, through accretion-cap clipping and
    subduction of continental cells riding oceanic plates;
  - #254: the quad deep ocean floor sits about 430 m shallower;
  - #257: whole-world sums assume the nominal node area;
  - #255: a fix or a tolerance decision for about 1.9% double-claimed territory.
- **Not a quad regression:** #256, Hm leaves its caps on both surfaces.
- **Performance (indicative, from shared-machine runs):** quad deform + topology takes 0.28–0.36×
  the line time per step; the whole step takes 0.55–0.66×. The comparison against the #147
  profile on a quiet machine is still to be done.
