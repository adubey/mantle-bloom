# Plate surface baseline: metrics, field policy, and invariants

Issue #228 Phase 0a. This document defines what must stay true while the authoritative plate
terrain moves from `PlateWithLines`/`ElevationLine` to a sparse quad representation. It covers:

- the mesh-quality metrics the new representation is judged by,
- how each persistent field must be treated when territory is remeshed, split, or merged,
- the invariants and proposed tolerances for each migration phase,
- what Phase 0a found that the `PlateSurface` interface (Phase 1) needs to provide.

It deliberately proposes no API. The measurements behind the numbers are in
[`analysis/issue228-phase0a/report.md`](../analysis/issue228-phase0a/report.md). They come from
[`bin/debug/characterize_plate_surface.py`](../bin/debug/characterize_plate_surface.py), which
turns any `.mbworld` into a deterministic JSON document.

## 1. Metrics

All distances are in multiples of the world's line spacing `s = line_spacing_rad(node_density)`
(62.5 km at the default density 4.0). A metric is **neutral** if it needs only node positions
(plus a plate-local frame for orientation), so it can be computed unchanged for a quad plate.
It is **legacy** if it only exists for `ElevationLine`s.

### 1.1 Sample-cloud metrics (neutral)

| Metric | Definition | Healthy lattice |
|---|---|---|
| Nearest-neighbour spacing | Distance to the nearest node of the same plate. | ≈1.0 |
| **Stacked fraction** | Share of nodes with a same-plate neighbour closer than 0.5 s. Two such nodes stand in for the same patch of crust. | 0 |
| Refinement jump | For each node, the largest ratio between its own nearest-neighbour distance and that of any of its 4 nearest neighbours, with stacked neighbours excluded. A 2:1-balanced adaptive mesh stays ≤ 2. | ≈1.0 |
| Anisotropy ratio | √(λmax/λmin) of the tangent-plane covariance of the 8 nearest non-stacked neighbour offsets, capped at 100. | 1.0–1.6 (the staggered rows give ~1.23 median) |
| Anisotropic fraction | Share of nodes with anisotropy ratio > 2. | ≈0 |
| One-dimensional fraction | Share of nodes at the cap: the neighbourhood is effectively a single row. | 0 |
| Row alignment | For anisotropic nodes, cos 2α between the major axis and plate-local east (the row direction). +1 means elongated along rows, −1 across them. This is how the anisotropy in #228 is measured: row-induced anisotropy shows up as a strong positive value. | n/a |
| Boundary node | The node's non-stacked neighbours within 1.5 s leave an empty angular sector > 135°. | |
| Boundary fraction / roughness | Share of boundary nodes; roughness = boundary nodes / 2√(πN), the count a disc of N nodes would have. | roughness ≈1–2 |
| Thin fraction | Boundary nodes with no interior neighbour: tendrils at most two nodes wide. | ≈0 |
| Boundary depth | Graph hops from each node to the nearest boundary node. | |
| Components / isolated nodes | Connected components with nodes linked within `DEFRAG_CONNECT_RADIUS_MULT` (2.5 s); isolated = singleton components. | 1 per plate / 0 |

### 1.2 Cell metrics (neutral definitions, measured on implied cells today)

A quad `(a, b, d, c)` is projected onto the tangent plane at its centroid. Today no real cells
exist, so the diagnostic builds *implied* cells: each pair of consecutive nodes on a row is
joined to the nearest-θ nodes of the next row up. These metrics become exact once cells are
authoritative.

| Metric | Definition |
|---|---|
| Aspect ratio | Longest edge / shortest edge. |
| Skew | max over corners of \|interior angle − 90°\|. Implied cells between staggered rows are inherently sheared (median ~13°), so this is a *floor* for comparison, not a target. |
| Folded | Any corner with a non-positive cross product (a concave or inverted cell). |
| Collapsed | Implied only: both lower nodes map to the same upper node, so the quad degenerates to a triangle. This is the signature of a row stub or a row much sparser than its neighbour. |
| Non-conforming | Any edge longer than 2.5 s: the cell bridges a gap instead of joining neighbours. |
| Area / s² | Cell area relative to a nominal cell. |
| Boundary-cell width | Covered by *thin fraction* and *boundary depth* above. The cell equivalent is the number of cells between opposite boundary edges. |
| Refinement jump | Level difference between edge-adjacent cells (quad). The sample-cloud refinement jump above is the representation-neutral proxy. |

### 1.3 Coverage and area (neutral)

| Metric | Definition |
|---|---|
| Uncovered / multiply-covered sphere | A fixed Fibonacci sample of the sphere (10⁶ points), counted against every plate's `contains_batch`. |
| Void fraction | Uncovered samples farther than 1.5 s from any node: real gaps, not outline slack. |
| Nodes inside another plate's outline | The existing stress-test invariant, evaluated on every node instead of a 20-node sample. |
| Co-located nodes | Nodes within 0.5 s of another plate's node (`plates.compute_node_overlap`). |
| **Voronoi area per node** | Each node's spherical-Voronoi cell over *all* plates' nodes (sums to exactly 4π), divided by the nominal `s²`. This is the area a node actually represents. The current conservation stats assume it is 1.0 everywhere. |
| Nominal area / sphere | `N·s² / 4π`. Above 1 means more nodes than the sphere has room for at nominal spacing, i.e. stacking or overlap. |
| Polygon area / nominal | Per plate: outline area / (node count · s²). |

### 1.4 Legacy line metrics

These have no quad equivalent. They quantify what the migration removes: one-node and two-node
lines, line-length histogram, arcs per row, **rows whose arcs overlap in θ** (and the
duplicated row length), along-row spacing and gaps over `IRREGULARITY_TOLERANCE`, missing rows.

## 2. Field policy

Each persistent value falls into one of these classes. The class decides how the value is
transferred when territory is remeshed (refined, coarsened, regularized, re-latticed), split,
or merged. "Current" is what the line code does today; "Phase 3" is the proposed rule.

| Class | Meaning | Remap rule (Phase 3) |
|---|---|---|
| **Extensive** | An amount per unit area: total = Σ value · area. | Area-weighted conservative transfer: Σ value·area is preserved exactly per operation. Refining splits a parent's amount across its children by area; coarsening sums the children. |
| **Intensive** | A local state or rate, not an amount. | Area-weighted mean on coarsen; interpolate (or copy the parent value) on refine. |
| **Categorical** | A code. | Coarsen: area-weighted majority vote, ties broken deterministically (lowest code). Refine: copy the parent value. Never interpolate. |
| **Boolean provenance** | Once true, it means "this material has had X". | Coarsen: logical OR. Refine: copy the parent value. |
| **Countdown / clock** | Time remaining or elapsed for an ongoing process. | Coarsen: max for "remaining", area-weighted mean for "elapsed". Refine: copy the parent value. |
| **Write-once history** | A timestamp that never changes after being set. | Coarsen: min over children that are set (oldest material wins; sentinel only if all are unset). Refine: copy the parent value. Must survive every operation. |
| **Derived** | Recomputable from other fields. | Never transfer independently; recompute after the fields it derives from have moved. |

### 2.1 Per-node fields (`elevation` + `ElevationLine.OPTIONAL_FIELDS`)

| Field | Class | Current remap (regularize / relattice & merge / fault shear) | Phase 3 rule and notes |
|---|---|---|---|
| `elevation` | Intensive, coupled to Hc/Hm by isostasy | crumple or `np.interp` / nearest / nearest | Move Hc/Hm conservatively, then recompute isostatic elevation and carry the *residual* relief (elevation − isostatic) as an intensive field. Never transfer it independently, or relief and column drift apart. Clipped to [−11 000, 9 000] m. |
| `crustal_thickness_m` (Hc) | **Extensive** (crustal volume) | interp / nearest + global Σ rescale / nearest | The main conserved quantity. Relattice rescales Σ Hc to preserve Σ Hc · *nominal* area. That preserves the proxy but not the actual unique-area volume when nodes are stacked (§4.2). Capped at `MAX_CRUSTAL_THICKNESS_M`; any excess over the cap must go to an explicit overflow sink, not be discarded silently. |
| `mantle_lithosphere_thickness_m` (Hm) | **Extensive** | interp / nearest + scaled with Hc / nearest | As Hc. Capped at `MAX_MANTLE_LITHOSPHERE_THICKNESS_M`. |
| `crust_type_code` | Categorical (0 inherit / 1 oceanic / 2 continental) | nearest | Before voting, resolve `INHERIT` against the owning plate so a merge of two differently typed plates doesn't flip meaning. |
| `lake_depth` | Extensive (water volume) | interp / nearest / nearest | Recomputed by hydrology every climate step, so conservation matters only within a step. |
| `glacier_depth` | Extensive (ice volume) | interp / nearest / nearest | Coupled to the eustatic water budget (`ocean_water_column_m`), so transfers must conserve it. |
| `silt_depth` | Extensive, non-decreasing | interp / nearest / nearest | Only meaningful inside a lake. |
| `soil_depth` | Extensive | interp / nearest / nearest | |
| `coal_deposit_m`, `oil_gas_deposit_m`, `mineral_deposit_m` | Extensive, non-decreasing | interp / nearest / nearest | A remap must not make any node's value drop below what the same material held before, except where shear brings in different material (fault shear already accepts that). |
| `soil_mineral_content`, `soil_organic_content` | Intensive, in [0, 1] | interp / nearest / nearest | Weight by `soil_depth` when coarsening (they are concentrations in the soil column). |
| `channel_depth`, `channel_width` | Intensive, **linear feature** | interp / nearest / nearest | Averaging washes out a one-node-wide river. On coarsen use max depth and keep the width of the deepest child. Better: refine *around* channels rather than coarsen through them. Hydrology rebuilds these from flow, so exact conservation is not required. |
| `divergent_age_myr` | Clock (Myr continuously divergent) | **reset to 0** (bug, §4.3) / nearest / nearest | Mean on coarsen, copy on refine. Drives `relax_young_oceanic_mantle_lithosphere`, so resetting it changes Hm physics. |
| `is_volcano` | Boolean provenance (never reverts) | interp > 0.5 threshold / nearest / nearest | OR on coarsen. The current interp-threshold can *drop* a single-node volcano between two non-volcano samples. |
| `volcano_active_years_remaining` | Countdown | interp / nearest / nearest | Max on coarsen. Interpolating a countdown invents partial lifetimes. |
| `elev_change_reason` | Categorical, diagnostic | nearest | Vote. Prefer structural codes (collision…fault, see `ELEV_CHANGE_STRUCTURAL_OVERRIDE_M_PER_MYR`) over geomorphic ones on ties. |
| `overlap_onset_years` | History, diagnostic (0 = not overlapping) | nearest | Min over non-zero children. Recomputed each step by `merge_split.update_overlap_tracking`. |
| `node_created_years` | **Write-once** (−1 sentinel = predates tracking) | nearest | Min over set children. Must survive every operation; the "Added/Removed Points" view and `GET /world/node_at` depend on it. |

Positions (`phi`, `theta` for lines) are not fields: they are topology. Rigid rotation changes
only `Plate.frame` and must never resample anything.

### 2.2 Per-plate state

`plate_id` (identity), `frame` (local→world rotation), `omega`, `crust_type` (nominal, majority
of effective node types after a split), `age_steps`, `internal_stress`. None of these are
remapped. Split and merge already define them (`_plates_from_node_masks`,
`majority_crust_type`).

### 2.3 World state that refers to plates or nodes

| State | Keyed by | Survives a remesh? |
|---|---|---|
| `pending_magma_parcels` | origin xyz + volume | Yes (geometry-based by design, see magma_transport.py) |
| `faults`, `fault_systems` | `plate_id` + plate-local (φ, θ) trace | Yes for remesh. Split/merge must re-home a fault to the plate that now owns its trace. |
| `earthquakes`, `gap_tracks`, `stranded_basin_tracks`, `removed_points_log`, `corner_notch_log` | world xyz | Yes |
| `collision_progress`, `overlap_progress` | `(plate_id, plate_id)` | Yes |
| `hydrology_cache`, `erosion_cache`, `node_*_cache`, `climate_cache` | global node order at build time | Derived. They must be invalidated whenever topology or node order changes, not just on rotation. |
| RNG streams in `lithosphere_plate._erupt_melted_nodes`, `_ignite_early_rift_volcanoes`, `volcanism.apply_volcanic_activity` | `(seed, year, plate_id, line_index)` | **No.** See §4.4. |
| `GET /world/elevation_point`, `/world/elevation_point_at` | `line_index` into `sorted_nonempty_lines` | UI API. Needs a representation-neutral node address. |

## 3. Invariants and proposed tolerances

Baseline values are from the reference worlds in the report. "Now" gives the fresh-world value
and the value in the 352.4 Myr reproducer.

### 3.1 Hard invariants (every representation, every step)

| Invariant | Tolerance |
|---|---|
| Rigid rotation changes `frame` only: every field array identical, local positions identical | Exact (bitwise) |
| Plate frame orthonormal, det = +1 | 1e−6 (existing stress test) |
| Hc, Hm within caps; elevation within [−11 000, 9 000] m | Exact |
| No NaN/inf in any field | Exact |
| Every field array has the plate's node/cell count | Exact |
| Derived caches (outline, boundary loops, row intervals, KD-trees) are functions of topology revision + frame | Rebuild from scratch must equal the cached value exactly |
| Diagnostic determinism: same world in → same JSON out | Byte-identical (verified) |

### 3.2 Phase 1 (`PlateSurface` over `PlateWithLines`)

Phase 1 must not change behaviour. The check is concrete: running
`characterize_plate_surface.py` on the reference worlds, with and without `--advance-steps`,
must reproduce the checked-in `analysis/issue228-phase0a/results/*.json` **byte for byte**.
Only `*.timings.json` may differ, and it should not get meaningfully slower.

### 3.3 Quad representation (Phases 2–5)

| Metric | Now: fresh → 352.4 Myr | Proposed quad bound |
|---|---|---|
| Stacked fraction | 0 → 27.5% | 0 (hard: authoritative topology cannot hold the same patch twice) |
| One-dimensional fraction | 0 → 4.2% | 0 outside a plate's last surviving strip, ≤ 0.1% |
| Anisotropic fraction (ratio > 2) | 0.01% → 5.4% | ≤ 1% at any age |
| Anisotropic row alignment (cos 2α) | +0.68 → +0.87 | \|mean\| ≤ 0.3: no preferred direction |
| Folded cells | 0 → 1e−5 (implied) | 0 (hard) |
| Collapsed cells | 1.7% → 2.8% (implied) | 0 by construction |
| Aspect > 4 | 0.05% → 0.25% | ≤ 0.1% (refine or repartition above that) |
| Skew > 45° | 0.03% → 0.06% | ≤ 0.1% |
| Non-conforming edges | 0.05% → 0.05% | 0 (neighbours share edges or hanging nodes at a 2:1 level jump) |
| Refinement jump | ≤ 1.1 → p99 1.9, max 4.4 | ≤ 2 (2:1 balance, hard) |
| Thin fraction | 0.003% → 4.7% | ≤ the line baseline at the same age, and ≤ 1% |
| Isolated nodes / extra components after defragment | 0 / 0 | 0 / 0 |
| Uncovered sphere | 1.3% → 2.3% | ≤ line baseline + 0.5 pp |
| Void (> 1.5 s from any node) | 0 → 0.013% | ≤ 0.05% |
| Multiply covered sphere | 1.6% → 0.5% | ≤ line baseline |
| Nodes inside another plate's outline | 1.8% → 1.4% | ≤ line baseline ("at least as strong", #228 acceptance) |
| Voronoi area / nominal, p05–p95 | 1.0–1.0 → 0.11–1.00 | 0.5–1.5. Any wider and extensive fields must carry an explicit per-cell area. |
| Nominal node area / sphere | 1.00 → 1.21 | 1 ± 0.05 |

### 3.4 Conservation tolerances

| Operation | Quantity | Tolerance |
|---|---|---|
| Rigid rotation | all fields | exact |
| Refine / coarsen / regularize | Σ value·area for every extensive field, per plate | 1e−9 relative (float rounding only) |
| Split, defragment | Σ over fragments = Σ of parent, every extensive field | 1e−9 relative |
| Merge | Σ over merged = Σ keep + Σ absorb − the explicitly logged overlap-dedup volume | 1e−9 relative, plus a logged dedup term (never an unlogged loss) |
| Write-once `node_created_years` | min over surviving material | exact |
| `is_volcano` count after coarsen | ≥ number of parent cells containing a volcano | exact |
| Subduction/retreat (deliberate sinks) | removed Hc/Hm | recorded by `phase_budget` (#216), never silent |
| Long-run parity (Phase 5), ensemble of seeds | land_fraction_node, Voronoi-weighted continental Hc, plate count, elevation p05/p50/p95, sea level | quad ensemble mean within 2σ of the line ensemble's seed-to-seed spread at 30/60/120 Myr |

Parity has to be statistical, not bitwise, because of §4.4.

## 4. Findings for the migration

### 4.1 Fragmentation is real and measurable

The 352.4 Myr reproducer matches #228 exactly: 4,955 lines, 1,452 one-node lines (29.3%), and
plate 66 with 459 of 780 lines at one node. The new metrics add that these stubs are
**anisotropic along rows**: 92% of anisotropic neighbourhoods are elongated along the row
direction. Plate 66's boundary is 17× longer than a disc of the same node count, and 71% of its
nodes are in tendrils at most two nodes wide.

### 4.2 Same-row arc stacking: duplicated territory inside a plate

In the reproducer, 506 of 2,632 rows hold two or more arcs (separate `ElevationLine`s at the
same φ) that **overlap in θ**, stacked up to 15 deep. Together they duplicate about 22,900
spacings of row length. 27.5% of all nodes (43,600) sit within 0.5 s of another node of the
*same plate*, almost all of them on a different line at the same φ. Nearly all are on oceanic
plates.

`plates._row_intervals` merges overlapping intervals, so outlines and containment never show
this. The node count does: the world carries 21% more nodes than the sphere has room for at
nominal spacing. Consequences:

- Σ Hc × nominal area overstates the actual (Voronoi-weighted) crustal volume by 14% (16.07 vs
  14.05 × 10⁹ km³). Hm is overstated by 11%, lake volume by 2.8×.
- Every conservation stat and #216's phase budget use the nominal proxy, so node-count churn
  in stacked rows shows up as crust creation or destruction. This is the missing
  "ground-truth unique-area integral" #216 asks for. `extensive_totals_voronoi_weighted` in
  the diagnostic provides one.
- A quad representation removes stacking structurally. The migration's conservation checks
  should be based on Voronoi-weighted or per-cell-area totals, or a correct quad will look as
  if it lost volume.

This looks like a line-mechanics bug (probably interior-subduction arc carving or the
split/defragment partition path) that should be tracked separately from #228.

### 4.3 `regularize_line` resets `divergent_age_myr`

`regularize_line` passes every `OPTIONAL_FIELDS` member to the new `ElevationLine` **except**
`divergent_age_myr`, so the constructor zero-fills it. It's the only field-by-field
construction site that drops a field (checked by AST scan). Every regularize pass therefore
marks a line's nodes as newly divergent, which `rheology.relax_young_oceanic_mantle_lithosphere`
(`divergent_age_myr < 30`) reads. It is a one-line fix plus a regression test, separate from
#228. Quad parity runs must be compared against a line baseline with this fixed, or the two
will differ in Hm for reasons unrelated to the representation.

### 4.4 Capabilities Phase 1 will need

Phase 1 owns the interface. These are requirements that surfaced here, not proposed methods:

1. **Stable, representation-neutral node identity.** Three RNG streams (above) are keyed by
   `line_index`. Any representation change reorders them, so no quad run can reproduce a line
   run bitwise. A neutral key would be stable per step, e.g. `(plate_id, cell id)` or a hash of
   the quantized local position. The UI's `elevation_point` endpoints also address nodes by
   line index.
2. **Per-node area.** Conservation needs an area per node or cell, not the global
   `lithosphere.node_area_m2` constant. §4.2 shows the constant is already wrong by 14% on
   long runs.
3. **Boundary loops with holes, as a derived view.** 12 of 35 plates in the reproducer have
   holes (66 loops in total). Consumers currently get a single keyholed polygon from
   `get_bounding_polygon()`.
4. **Explicit adjacency.** Every neutral metric here had to rebuild adjacency from KD-tree
   radius queries. Hydrology does the same. A surface-level neighbour list would replace both.
5. **Field metadata.** The policy table in §2.1 is currently only in comments and scattered
   `np.interp`/nearest calls. A single registry of field name → class (dtype, default or
   sentinel, remap class) would let Phase 3 implement remapping once and let the contract
   tests check that no field is missing. §4.3 is exactly that kind of omission.
6. **Topology revision distinct from geometry revision.** Derived caches keyed on node order
   (`hydrology_cache`, `erosion_cache`) need to be invalidated on topology change, and
   rotation must not invalidate them.
