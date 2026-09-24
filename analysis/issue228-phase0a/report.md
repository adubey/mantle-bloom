# Issue #228 Phase 0a: plate-surface baselines

Measured baselines for the `PlateWithLines` → sparse-quad migration. The metric definitions,
field policy, and proposed tolerances these numbers support are in
[`docs/plate-surface-baseline.md`](../../docs/plate-surface-baseline.md).

## Reference worlds

| World | Config | Source | sha256 (of the .mbworld) |
|---|---|---|---|
| `seed0-000.0myr` … `seed0-060.0myr` | seed 0, node density 4.0, climate 4.0, fluid 2.0, 100 kyr steps (the issue #147 profile config) | `bin/debug/capture_issue228_baselines.py` on `e09020b` | recorded in each `results/<world>.json` → `world.sha256` |
| `mantle-bloom-seed804913535-352400000y` | seed 804913535, density 4.0, 3,524 steps | the #228 reproducer, `~/Downloads/…352400000y.mbworld` | `ae1a03eb…fafb86` |

The seed 0 checkpoints are 115–140 MB each and are not checked in: `worlds/.gitignore`
excludes them. Regenerating all four takes about 90 minutes:

```sh
backend/.venv/bin/python bin/debug/capture_issue228_baselines.py        # -> analysis/issue228-phase0a/worlds/
```

Stepping is deterministic: two separate processes stepping the 10 Myr checkpoint 10 more steps
wrote byte-identical results. So the regenerated files should match the recorded hashes,
provided the code and dependencies haven't changed.

Characterize and summarize:

```sh
W=analysis/issue228-phase0a/worlds R=analysis/issue228-phase0a/results
backend/.venv/bin/python bin/debug/characterize_plate_surface.py --render --out $R \
    $W/seed0-0{00,10,30,60}.0myr.mbworld ~/Downloads/mantle-bloom-seed804913535-352400000y.mbworld
backend/.venv/bin/python bin/debug/characterize_plate_surface.py --advance-steps 10 --out $R \
    $W/seed0-0{10,60}.0myr.mbworld ~/Downloads/mantle-bloom-seed804913535-352400000y.mbworld
python3 bin/debug/summarize_plate_surface.py $R/seed0-*myr.json $R/mantle-bloom-*y.json
```

The diagnostic takes about 5 s per world (plus stepping time with `--advance-steps`). Its JSON
output is byte-identical across reruns.

## Results

| Metric | 0 Myr | 10 Myr | 30 Myr | 60 Myr | 352.4 Myr (reproducer) |
|---|---:|---:|---:|---:|---:|
| Plates | 19 | 29 | 31 | 35 | 35 |
| Nodes | 130,587 | 130,475 | 134,660 | 139,213 | 158,589 |
| Lines | 1,796 | 2,642 | 2,983 | 3,386 | 4,955 |
| One-node lines | 0.3% | 3.8% | 9.7% | 15.0% | 29.3% |
| Rows with overlapping arcs | 0 | 38 | 110 | 186 | 506 |
| Max arc stacking depth | 0 | 3 | 5 | 7 | 15 |
| **Stacked nodes** (< 0.5 s, same plate) | 0.00% | 0.83% | 6.71% | 12.29% | 27.50% |
| Anisotropic nodes (ratio > 2) | 0.01% | 0.08% | 0.22% | 0.73% | 5.38% |
| One-dimensional nodes | 0.00% | 0.00% | 0.00% | 0.06% | 4.24% |
| Anisotropy row alignment (cos 2α)¹ | +0.68 | −0.01 | +0.06 | +0.46 | +0.87 |
| Refinement jump p99 | 1.00 | 1.12 | 1.60 | 1.77 | 1.88 |
| Boundary nodes | 3.5% | 5.0% | 5.6% | 6.3% | 10.9% |
| Thin (tendril) nodes | 0.00% | 0.14% | 0.24% | 0.32% | 4.71% |
| Isolated nodes | 0 | 0 | 0 | 0 | 0 |
| Implied cells: collapsed | 1.75% | 2.39% | 2.67% | 2.89% | 2.81% |
| Implied cells: aspect > 4 | 0.047% | 0.030% | 0.048% | 0.080% | 0.245% |
| Implied cells: skew > 45° | 0.035% | 0.030% | 0.012% | 0.017% | 0.060% |
| Implied cells: folded | 0 | 0.0008% | 0 | 0.0008% | 0.0013% |
| Implied cells: non-conforming | 0.051% | 0.034% | 0.009% | 0.013% | 0.054% |
| Uncovered sphere | 1.33% | 2.44% | 2.74% | 2.95% | 2.33% |
| Void sphere (> 1.5 s from any node) | 0 | 0.044% | 0.011% | 0.151% | 0.013% |
| Multiply covered sphere | 1.60% | 0.62% | 0.13% | 0.14% | 0.48% |
| Nodes inside another plate's outline | 1.77% | 0.80% | 0.41% | 0.55% | 1.38% |
| Plates with holes | 0 | 0 | 2 | 3 | 12 |
| Nominal node area / sphere | 1.000 | 0.999 | 1.031 | 1.066 | 1.215 |
| Voronoi area / nominal, p05 | 1.00 | 0.99 | 0.51 | 0.49 | 0.11 |
| Voronoi area / nominal, p95 | 1.00 | 1.01 | 1.02 | 1.02 | 1.00 |
| Land fraction (node) | 15.3% | 17.2% | 15.9% | 13.9% | 28.5% |
| Sea level (m) | 0 | −366 | −700 | −869 | −826 |
| Σ Hc, nominal area (10⁹ km³) | 8.56 | 8.51 | 8.94 | 9.41 | 16.07 |
| Σ Hc, Voronoi area (10⁹ km³) | 8.56 | 8.48 | 8.69 | 8.90 | 14.05 |
| Σ Hm, nominal area (10⁹ km³) | 19.11 | 25.18 | 28.46 | 28.99 | 38.89 |
| Σ Hm, Voronoi area (10⁹ km³) | 19.10 | 25.07 | 27.80 | 27.59 | 35.07 |
| Continental Hc (effective), nominal | 4.39 | 4.45 | 4.47 | 4.43 | 9.09 |
| Continental Hc (effective), Voronoi | 4.39 | 4.43 | 4.44 | 4.39 | 9.01 |

¹ Averaged over anisotropic nodes only. At 0 and 10 Myr there are too few of them (19 and
~100) for the sign to mean anything. At 352 Myr, 92% of anisotropic neighbourhoods are
elongated along the plate-local row direction.

Seed 0 over time (from `results/seed0-timeseries.jsonl`, written every 10 steps):

| Myr | Plates | Nodes | Lines | One-node lines | Σ Hc nominal (10⁹ km³) | Land (node) | Sea level (m) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 19 | 130,587 | 1,796 | 5 | 8.56 | 15.3% | 0 |
| 10 | 29 | 130,475 | 2,642 | 100 | 8.51 | 17.2% | −366 |
| 20 | 31 | 132,688 | 2,905 | 198 | 8.73 | 16.6% | −582 |
| 30 | 31 | 134,660 | 2,983 | 289 | 8.94 | 15.9% | −700 |
| 40 | 31 | 136,187 | 3,034 | 387 | 9.15 | 15.2% | −753 |
| 50 | 33 | 137,686 | 3,255 | 494 | 9.30 | 14.5% | −811 |
| 60 | 35 | 139,213 | 3,386 | 508 | 9.41 | 13.9% | −869 |

Renders (Eckert IV, 1100×611, 128-colour palette): `results/<world>.{plates,platesDetail,elevation}.png`.
The 352.4 Myr `platesDetail`/`elevation` renders show the row-streak artifacts from #228. The
60 Myr render already shows streaking in young ocean floor and a few void wedges between
plates (the 0.15% void fraction).

## Findings

1. **Fragmentation grows steadily from the start, not just in very long runs.** One-node lines
   go from 0.3% to 15% of all lines in the first 60 Myr (about 8 more per Myr), then to 29% by
   352 Myr. Tendrils, one-dimensional neighbourhoods, and row-aligned anisotropy mostly appear
   late (after 60 Myr), concentrated in a few heavily subducted oceanic plates. Plate 66 in the
   reproducer has a boundary 17× the length of a disc with the same node count, and 71% of its
   nodes are in tendrils.

2. **Same-row arc stacking is the largest representation defect found, and it biases
   conservation accounting.** Arcs sharing a φ overlap in θ, so a plate holds the same patch
   of crust more than once. The stacked fraction goes 0% → 0.8% → 6.7% → 12.3% (0/10/30/60 Myr)
   → 27.5% (352 Myr). The outline code merges overlapping intervals, so this never shows up in
   containment or rendering. It does inflate node count: by 352 Myr there are 21% more nodes
   than the sphere has room for.

   Because every conservation stat is Σ value × nominal node area, stacking reads as crust
   growth. Between 0 and 60 Myr, nominal Σ Hc rises 10% (8.56 → 9.41) but Voronoi-weighted
   Σ Hc rises only 4% (8.56 → 8.90). By 352 Myr the nominal figure overstates Hc by 14% and Hm
   by 11%. Continental Hc is barely affected (≤ 1%); stacking is almost entirely oceanic. This
   bears directly on #216, whose phase budget uses the nominal proxy.

3. **`regularize_line` zero-fills `divergent_age_myr`.** It is the only field-by-field
   `ElevationLine` construction that omits a field (checked by AST scan). Confirmed by direct
   call: a line with divergent age 12 Myr comes back as 0 while `node_created_years` is
   carried. It feeds `relax_young_oceanic_mantle_lithosphere`. See
   docs/plate-surface-baseline.md §4.3.

4. **Implied-cell quality is fine where rows are healthy.** Folded, non-conforming, and
   high-aspect implied cells are all well under 0.3% even at 352 Myr. Collapsed cells (1.7%
   even at generation, from staggered rows of different lengths) are the structural cost of
   the row representation. The migration problem is topology (stubs, stacking, tendrils,
   holes), not the shape of cells in the interior.

5. **Coverage is roughly stable over time.** Uncovered sphere stays at 1.3–3%, true void under
   0.2%, and nodes inside a foreign outline at 0.4–1.8%. Most of the uncovered area is outline
   slack between the staircase polygons, not missing crust. These are the numbers the quad
   representation must match or beat (#228: "at least as strong").

6. **Line/quad parity can only be statistical.** Three RNG streams are keyed by `line_index`
   (`lithosphere_plate.py:396`, `:424`, `volcanism.py:91`).

## Derived-geometry timings

These come from `*.timings.json`, measured on this machine while the capture run was using a
core, so they're indicative only. Rebuilding every plate outline takes 14 ms at 0 Myr, 21 ms at
10–30 Myr, 25 ms at 60 Myr, and 51 ms at 352 Myr. `contains_batch` over 10⁶ sphere samples
across all plates takes 0.18 s at 0 Myr and 0.25 s at 60 Myr. Stepping 10 steps takes 54 s from
10 Myr and 94 s from 60 Myr (with climate). At 352 Myr it's 81 s from the loaded save, which has
more nodes but ran on an unloaded core.

## Using this for Phase 1

`results/*+10steps.json` step each reference world 10 steps (1 Myr) forward on commit
`e09020b` and characterize the result. After the `PlateSurface` refactor, rerunning the
commands above must reproduce every `results/*.json` byte for byte (not `*.timings.json`). A
difference means behaviour changed, whether or not the unit tests notice. If Phase 1 started
from a different `main` commit than `e09020b`, regenerate these files on that commit first.
