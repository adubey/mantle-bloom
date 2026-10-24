# Issue #468: normalize boundary loads by front length (RFC 0012)

## Principle and invariant

Simulation design principles **4 (treat nodes as finite areas)** and **5 (keep units
consistent)** apply here. The force code sampled an along-front line load at every nearby
plate node, so its total depended on how many rows fell inside the boundary band and on node
density. The repair assigns each distinct nearest-neighbour front node one spacing of line
length and shares that length across rows mapped to it.

The equation-level invariant is: with fixed physical boundary fields, slab pull, sunk-slab
drag, ridge push, and contact drag per represented metre stay fixed as density, band reach,
or repeated rows change. The synthetic probe measured all four ratios between 0.9999999999999997
and 1.0000000000000004 over its 48 cases (densities 1/2/4, reaches 1.5/3/6 spacings, and one
to two repeated input sets). Represented length covered 19.4–91.9% of the geometric span,
depending on side, density, and reach; the force is normalized to the represented front, while
the probe records the coverage gap rather than silently filling it.

The implementation changes no force constants. It shares `ds / rows_for_front_node` across
samples for slab pull, sunk-slab drag, ridge push, and contact drag. Contact traction is
converted through the same weighted line length and each node's accounting area. Basal drag
continues to integrate over plate area.

## Replay comparison

Four fresh worlds ran to 500 Myr at 2 Myr steps on origin/main `f600367` and branch
`3a8a48d` (the implementation commit). These are model-sensitivity observations, not a claim
that four seeds establish a new planetary benchmark. The G1 calibration constants remain
unchanged. Order-dependence checks were disabled; this replay isolates default evolution.

| Seed | Land fraction, main → branch | Plate-rate p50, cm/yr, main → branch | Integrator cap share, main → branch | Runtime, s, main → branch |
|---:|---:|---:|---:|---:|
| 326494039 | 0.229 → 0.161 | 3.79 → 3.30 | 1.82% → 2.14% | 703 → 766 |
| 748935664 | 0.116 → 0.095 | 4.26 → 3.65 | 2.64% → 2.23% | 847 → 1242 |
| 580014012 | 0.279 → 0.264 | 2.84 → 2.21 | 0.51% → 0.27% | 660 → 1226 |
| 111065103 | 0.383 → 0.328 | 3.65 → 2.90 | 2.63% → 1.04% | 716 → 1343 |

The rate p5/p50/p95/p99, clamp floor share, continental Hc/Hm p10/p50/p90/p99, and
continental-node counts are in [comparison.json](comparison.json). Land fraction and median
plate rate decreased in all four sampled branch runs; the size varies by seed. These runs do
not establish that normalization alone should move long-term land or plate-rate distributions
in that direction. Runtime was slower on branch in all four comparisons and highly variable
between seeds; the audit is not a controlled performance benchmark.

## Reproduction

The baseline worktree was detached at `f600367`; the branch was `3a8a48d`.
`bin/debug/audit_force_rules.py` was run serially once per seed in each tree:

```bash
NUMBA_CACHE_DIR=/private/tmp/issue468-numba-cache \
MB_SOURCE_ROOT="$PWD" PYTHONPATH="$PWD/backend" \
backend/.venv/bin/python bin/debug/audit_force_rules.py \
  --seed <seed> --myr 500 --step-myr 2 --out /private/tmp/issue468-<tree>-<seed>.json
```

`bin/debug/probe_boundary_load_normalization.py` ran the 48-case synthetic sweep in
`comparison.json`. Full replay output stayed in `/private/tmp`; only the compact comparison
and inputs used for this report are checked in. Glacier depth clipping diagnostics appeared in
some runs; they were not altered or used as evidence for this issue.
