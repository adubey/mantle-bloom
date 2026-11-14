# Issue 876: do extensional-load maps line up with RFC 0018's cuts?

Issue [#876](https://github.com/adubey/mantle-bloom/issues/876) explores physics-first rift candidates: computational
candidates (cuts the search evaluates), physical candidates (zones pulled apart and weak enough to deform, which
leave evidence when they fail) and successful rifts. This report is the exploration's step 2. It asks whether a
plate's extensional load relative to rock strength, measured before any eigensolve, matches where RFC 0018 finds
structurally valid cuts and positive breakup probability. RFC 0018's outcomes are the comparison, not ground truth.
No production code, gate or RFC text changed.

## Method

`bin/debug/compare_rift_load_maps.py` uses the 125 split opportunities from [issue 848](../issue848/README.md):
the 20 G0 snapshots (2 Myr) and 8 mature checkpoints from the #848 development seeds (1091931222 at 50 and 100 Myr;
679647835, 229148564 and 1333090219 at 50 and 100 Myr). For each opportunity it builds two maps:

- **Drag map (no solve).** For each node, the largest positive rate at which mantle flow minus the plate's rigid
  rotation pulls apart the ends of its edges. It is divided by rock resistance: RFC 0018's node resistance with
  the flow-residual factor neutralized, so it reflects thickness, cratons, belts, plume heat, faults and thinning.
- **Shell map (one sparse solve per plate).** The thin-shell (#699) incremental extensional resultant under the
  full per-node force field (slab pull, slab drag, ridge push, basal drag, contact drag), divided by the #487
  yield-envelope column capacity, Moho thermal lag included. A ratio of 1 or more means the column yields. Two
  plates could not be meshed.

It then runs the rebuilt mass-balanced generator with the exact repeat shortcut (#848): 1,491 candidates. 39
opportunities have a structurally valid candidate and 20 have positive `P_break`; 117 candidates are valid, 61
of them with positive `p`. Scores are rank AUCs, where 0.5 means no information. Compact results are in
`load-maps.json`.

## Results

| Question | Drag map | Shell map | Node count |
|---|---:|---:|---:|
| Plate has a valid candidate (plate p95) | 0.51 | 0.56 | 0.86 |
| Plate has positive `P_break` (plate p95) | 0.66 | 0.67 | 0.81 |
| Cut is valid (median along cut) | 0.49 | 0.48 | |
| Valid cut has positive `p` (median along cut) | 0.87 | 0.89 | |

1. **Structural validity is unrelated to load.** No map statistic separates valid from invalid cuts, at
   opportunity or cut level (AUC 0.47 to 0.56). Plate size does (0.86). RFC 0018's validity is a geometric
   filter, so its valid cuts are computational candidates. Nothing in the method measures physical candidates.
2. **Breakup probability follows load.** Among valid cuts, those with positive `p` run through more loaded rock.
   The drag result is partly circular, because RFC 0018's opening term and the drag map both come from mantle
   flow. The shell map, which also carries slab pull, is the more independent check.
3. **Cuts are far from yield.** The median shell ratio along a valid cut is 0.068 (0.042 for invalid). Its
   lowest tenth is 0, meaning the line crosses rock in compression, and the median share of cut nodes at or
   past yield is 0. Under one step's incremental load, intact rock carries about 5% of its capacity. A
   breakthrough threshold would therefore need damage or load to accumulate in the rock over many steps.
4. **No useful plate-level skip.** A zero-false-negative threshold on any map skips at most 13 of 125
   opportunities.
5. **Cost.** Drag map 0.18 s for all 125 plates. Shell map 204 s, against 364 s for the RFC 0018 search (with
   the shortcut), so it is not a cheap substitute as built. Reuse across steps, warm starts and a coarser mesh
   are untested.

## Limits

- 28 worlds, and opportunities within a world are not independent.
- The shell stress is an increment from zero load: it ignores prior stress, and the capacity has no inherited-
  weakness term.
- The drag map's divergence is a local proxy, not a stress.

## Reproduce

Run label `issue876-load-maps`. The G0 prototype is `bin/debug/measure_rfc0018_g0.py` on
`origin/codex/issue-797-evidence-rfc0018-g0` (not on `main`); the snapshots are the G0 snapshot set; the
checkpoints come from `bin/debug/generate_issue848_world.py`.

```bash
backend/.venv/bin/python bin/debug/compare_rift_load_maps.py --g0-script <measure_rfc0018_g0.py> \
    --snapshots <rfc0018-g0-snapshots> --seeds <the 20 G0 seeds> --saves <the 8 checkpoints> \
    --out analysis/issue876/load-maps.json --rows <private-dir>/load-map-rows.jsonl
```
