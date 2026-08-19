# Simulation Model

## Table of contents

- [Why not a grid](#why-not-a-grid)
- [Plate-local frames](#plate-local-frames)
- [Initial plate generation](#initial-plate-generation)
- [Mantle flow](#mantle-flow)
- [Boundary evolution](#boundary-evolution)
- [Line regularization](#line-regularization)
- [Merge and split](#merge-and-split)
- [Whole-sphere coverage (gap-filling)](#gap-filling)
- [Volcanism](#volcanism)
- [Boundary point reassignment](#reassignment)
- [Projections](#projections)
- [Render image](#render-image)
- [Rotating the view](#rotating-the-view)
- [Plate Inspector](#plate-inspector)
- [Climate](#climate)
  - [Biomes](#biomes)
- [Erosion](#erosion)
- [Bathymetry](#bathymetry)
- [Hydrology (rivers and lakes)](#hydrology)
- [Glaciation](#glaciation)
- [River Inspector](#river-inspector)
- [Coastline](#coastline)
- [Known simplifications](#known-simplifications)

<a id="why-not-a-grid"></a>
## Why not a grid

mantle-bloom does not model the planet as an equirectangular lat/lon grid, and deliberately
so: a fixed lat/lon grid brings a cluster of problems. Pole cells need an artificial
full-clique patch to be mutually adjacent, falloff radii need a
[Dijkstra](https://en.wikipedia.org/wiki/Dijkstra%27s_algorithm) search weighted by real 3D
chord distance to correct for latitude-dependent cell size, and -- the significant one --
elevation/sediment carried by a rotating plate has to be resampled onto the fixed grid every
step via nearest-neighbor
[semi-Lagrangian](https://en.wikipedia.org/wiki/Semi-Lagrangian_scheme) backward-advection,
which doesn't conserve mass where a plate stretches or compresses.

mantle-bloom drops the grid entirely: plates are spherical polygons, and elevation lives on
polylines that rotate *exactly* with their plate (no resampling at all for ordinary motion)
and are only ever touched where a boundary actually creates or destroys crust.

<a id="plate-local-frames"></a>
## Plate-local frames

Each plate owns a rotation matrix `frame` (world = `frame @ local`, see
`geometry.plate_frame_from_seed`) defining its own local spherical coordinate system
`(phi, theta)`, with local `(phi=0, theta=0)` mapping to the plate's original seed point.
Terrain is stored as `ElevationLine`s: each one a fixed plate-local latitude `phi`, holding
elevation samples at plate-local longitude nodes `theta` -- literally a local graticule glued
to the plate.

- **Rotating a plate is exact.** Advancing `frame` by composing the incremental rotation
  from the plate's current Euler pole/rate ([Rodrigues'
  formula](https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula),
  `geometry.rotation_matrix_from_omega`) updates every node's world position automatically,
  with no interpolation and no lost material -- because ordinary rotation never resamples
  anything, there's no semi-Lagrangian resampling step to lose or duplicate mass in the
  first place.
- **Equidistant and parallel, by construction.** Equal `delta-phi` between lines is already
  physically equidistant (meridional spacing on a sphere doesn't depend on latitude), and
  each line's `delta-theta` node spacing is chosen from that line's angular radius
  (`cos(phi)`) to hit `TARGET_LINE_SPACING_KM` (`plates.py`, default 125 km) -- this is what
  avoids the latitude-distortion problem a fixed-resolution lat/lon grid has, where node
  spacing narrows sharply toward the poles.
- **Irregular intervals at boundaries, naturally.** A line only exists for the theta-range
  currently inside the plate's territory; nodes at that cutoff are the ones boundary
  evolution adds or removes (see below), which is exactly where irregular spacing should
  show up.

<a id="initial-plate-generation"></a>
## Initial plate generation (`plates.py`)

Total plate count isn't asked of the caller: `num_plates`, if not given, is drawn from the
seed's own RNG stream (`plates.MIN_AUTO_PLATES = 8` to `plates.MAX_AUTO_PLATES = 20`), so the
world still tiles itself the same deterministic way for a given seed without the caller
having to pick a number. That many seed points are scattered uniformly on the unit sphere
(normalized Gaussian samples), and each gets a plate-local frame built from its own seed
(`geometry.plate_frame_from_seed`).

Three generation choices *are* user-facing -- the UI's "continental plates" and "initial
land" sliders, both 0 to 1 (percent in the UI), defaulting to
`DEFAULT_CONTINENTAL_FRACTION = 0.70` and `DEFAULT_LAND_FRACTION = 0.29`, plus a "point
density" choice (`NODE_DENSITY_CHOICES = (1.0, 4.0)`, see below).

- `continental_fraction`: when given, `round(continental_fraction * num_plates)` plates
  (`rng.choice`, without replacement) are made continental instead of the usual independent
  `CONTINENTAL_FRACTION = 0.4` coin flip per plate, and `num_plates` is bumped up if needed
  so there's still room for at least `MIN_OCEANIC_PLATES` of real ocean floor regardless of
  how high a fraction was requested.
- `land_fraction`: independently controls how much of the *whole sphere* -- not just of
  continental crust -- starts above sea level. This is a genuinely separate knob from
  `continental_fraction`, because continental crust from randomly-seeded Voronoi cells
  doesn't cover a *fixed* area fraction just because a *plate-count* fraction was requested
  (cells vary in size), and because not all continental crust needs to be dry land (compare
  real continental shelves). `_land_noise_threshold` does the work: a one-off whole-sphere
  sweep (independent of any plate's own lattice, at the coarser
  `LAND_FRACTION_SAMPLE_SPACING_KM = 150`, since this only needs to be a statistically
  representative sample) measures the *actual* continental area fraction, then finds the
  noise quantile that would put exactly the right sub-fraction of continental crust above
  sea level to hit the whole-sphere target (capped at 1.0 -- if there isn't enough
  continental area to reach the requested land fraction, every bit of it becomes land and
  that's as close as generation gets). Continental crust's elevation formula becomes
  `amp * (noise - threshold)` instead of the usual `BASE_CONTINENTAL_M + amp * noise` when a
  threshold is available; oceanic crust is untouched (its own base/amplitude make crossing
  sea level implausible regardless). Confirmed directly: at the defaults, measured land
  fraction across several seeds lands within about a percentage point of 29%.

Each plate's elevation lines are populated by `plates.iter_local_lattice`: sweep a full
plate-local `(phi, theta)` lattice at `TARGET_LINE_SPACING_KM` resolution (or a finer one --
see "Elevation point density" below), and for every candidate node, keep it only if this
plate's seed is the *nearest* seed to it (`cKDTree` against all seeds) -- the defining
property of a spherical [Voronoi
diagram](https://en.wikipedia.org/wiki/Voronoi_diagram), computed directly rather than via
an explicit polygon-construction step. Every node ends up owned by exactly one plate, so the
initial tiling has no gaps and no overlaps *by construction* -- there's nothing to
separately verify. Kept nodes get a base elevation by crust type
(`BASE_CONTINENTAL_M = 200`, `BASE_OCEANIC_M = -3800`, continental overridden by
`land_fraction`'s threshold when given, see above) plus a smooth noise texture (`noise.py`,
a small sum of sinusoids with random frequency/phase -- not true gradient noise, just enough
texture to not look perfectly flat).

The same lattice-sweep helper (`plates.build_lines_from_lattice`) is reused by plate merging
(see [Merge and split](#merge-and-split)) -- the only other place a full-footprint sweep is
needed.

**Elevation point density.** The UI's "point density" choice (`node_density`, 1x default or
4x) scales `TARGET_LINE_SPACING_RAD` down via `plates.line_spacing_rad` (halved at 4x --
node count for a fixed area scales with the *square* of resolution, so 4x the nodes needs
half the spacing, not a quarter). Stored on `World.node_density`, set once at generation and
read for that world's entire life, not just at the moment it's generated: every later module
that builds new elevation-line nodes or derives a distance/count threshold from
`TARGET_LINE_SPACING_RAD` -- `line_regrid.py`'s periodic regularization,
`boundary.py`'s per-step growth/merge/reach thresholds, `gaps.py`'s coverage-gap detection
and absorption/spawning, `merge_split.py`'s plate-merge contact distance and split-size
floor, `volcanism.py`'s volcanic-field clustering/coverage -- calls `line_spacing_rad(world.
node_density)` (or scales its own reference constant by the same ratio) instead of reading
the bare module constant. This matters because it's not just a generation-time cosmetic
choice: `line_regrid.py`'s regularize pass in particular runs unconditionally every
`REGULARIZE_INTERVAL_STEPS` steps and, before this threading existed, always resampled a
line back down to the *reference* spacing regardless of what density the world was actually
generated at -- confirmed directly as a real bug during development, a 4x-density world's
own node count reverting to the 1x baseline within the first handful of steps. Every
distance-based threshold derived from `TARGET_LINE_SPACING_RAD` (e.g. `boundary.py`'s
`EXTEND_THRESHOLD_RAD`) scales linearly with the new spacing; every absolute node-*count*
constant tied to a fixed physical area (`gaps.MIN_GAP_POINTS`, `merge_split.SPLIT_MIN_NODES`,
etc.) scales with `node_density` directly, not its square root, since it's already an area
-- see each constant's own comment for the exact reasoning, which predates this option (the
same rescaling used to happen as a one-off hardcoded code change whenever
`TARGET_LINE_SPACING_KM` itself changed; this option just makes it a per-world runtime
choice instead). Genuine fixed physical distances unrelated to sampling resolution (e.g.
`boundary.COLLISION_RANGE_KM`, a real ~400km-wide collision belt) are deliberately *not*
scaled -- only thresholds explicitly defined as multiples of `TARGET_LINE_SPACING_RAD` are.
4x density comes with a real, continuous performance cost, not just a one-time generation
cost -- confirmed directly, roughly a 5x slower per-step time (not just 4x, since flow-
routing/reassignment passes are `O(n log n)` rather than linear) -- which the UI surfaces as
a short note when 4x is selected.

<a id="mantle-flow"></a>
## Mantle flow (`mantle.py`)

A handful of upwelling/downwelling convection centers are placed via a **cubed-sphere
mapping**: a point on a cube face `(u, v, +-1)`, normalized to the unit sphere
(`mantle.cube_to_sphere`) -- "flow in a cube, projected to the sphere," chosen for even
coverage with no pole clustering. Each center contributes a tangential flow vector pointing
away from it (upwelling, positive strength) or toward it (downwelling, negative strength),
with Gaussian falloff by angular distance (`mantle.flow_at`).

Every step, each plate samples this field at its own current footprint (every elevation-line
node plus its boundary loop -- already available, no separate sampling grid needed) and fits
the best-fit rigid rotation via ordinary least squares: minimize
`sum |omega x p_i - v_i|^2`, a linear problem in `omega` solved as one 3x3 system per plate
(`mantle.fit_euler_pole`). `omega`'s direction is the [Euler
pole](https://en.wikipedia.org/wiki/Euler_pole), its magnitude the rotation rate -- this is
the real plate-tectonics formalism, not a hand-set per-plate velocity vector. The new target
is blended with the plate's current `omega` (`VELOCITY_DAMPING = 0.3`, so a plate
accelerates smoothly rather than snapping) and clamped to a plausible speed range
(`MIN/MAX_PLATE_RATE`, equivalent to 0.5-15 cm/yr at `PLANET_RADIUS_KM = 6371`).

<a id="boundary-evolution"></a>
## Boundary evolution (`boundary.py`)

There's no maintained shared-edge structure between plates -- every step, each plate's
elevation-line nodes are matched against a fresh k-d tree of every other plate's current
nodes (`scipy.spatial.cKDTree`). This is self-healing every step rather than requiring an
always-consistent topology, and it's what makes merge/split tractable without a general
spherical polygon-boolean library.

For each node within `MAX_BOUNDARY_EFFECT_RAD` (the widest reach any single effect below
needs -- currently `COLLISION_RANGE_RAD`, 400km) of some other plate's nearest node, the two
plates' relative velocity at that point (from their `omega`s) is decomposed against the
direction toward the neighbor into a **closing rate**: positive means this plate's material
is moving toward the neighbor's (convergent), negative means moving apart (divergent);
`TRANSFORM_RATE_THRESHOLD` (~1 cm/yr equivalent) separates both from transform.

**Convergent boundaries aren't a single effect** -- what happens depends on both plates'
crust type, and how far the effect reaches (and its shape with distance) differs by case:

- **Continent-continent collision** (both plates continental) -> elevation rises
  (`CONVERGENT_MOUNTAIN_RATE_M_PER_MYR`), scaled by an intensity that fades from 1 at zero
  distance to 0 at `COLLISION_RANGE_RAD` (400km) -- a broad crumple zone, matching how wide a
  real collision belt is (e.g. the Himalaya/Tibetan Plateau).
- **Oceanic-under-continental subduction** (continental plate, oceanic neighbor) -> elevation
  rises on the continental side too, but shaped as a **band** (`_band_intensity`): zero right
  at the boundary, peaking at the midpoint of `SUBDUCTION_ARC_INNER_RAD`..
  `SUBDUCTION_ARC_OUTER_RAD` (100-300km), zero again past the outer edge -- a volcanic arc
  forms where the subducting slab has descended deep enough to melt, not at the trench
  itself. This is the one non-monotonic shape here; every other effect peaks at the boundary
  and decays outward.
- **A subducting oceanic plate's own trench** (oceanic plate, any neighbor) -> elevation
  falls (`CONVERGENT_TRENCH_RATE_M_PER_MYR`), fading from 1 at zero distance to 0 at
  `FAR_THRESHOLD_RAD` (~200km) -- unaffected by the neighbor's crust type or by the two cases
  above.
- **Transform** -> elevation rises modestly (`TRANSFORM_UPLIFT_RATE_M_PER_MYR`, a fraction of
  the mountain-building rate -- real strike-slip motion produces at most local pressure-ridge
  relief, not real mountains), fading from 1 at zero distance to 0 at `TRANSFORM_RANGE_RAD`
  (50km).
- **Divergent** -> elevation relaxes exponentially toward a ridge (`oceanic`,
  `DIVERGENT_RIDGE_TARGET_M = -1500`) or rift (`continental`, `DIVERGENT_RIFT_TARGET_M =
  -200`) target -- new crust forming at the boundary. Oceanic ridge spreading is unchanged
  (fades from 1 at zero distance to 0 at `FAR_THRESHOLD_RAD`, ~200km); continental rifting
  reaches much farther (`RIFT_RANGE_RAD`, 300km) -- stretching and thinning the crust
  subsides land well beyond the fault line itself, not just right at it.
- **Structural growth/shrink**, applied independently at each line's two ends (the true
  edge of that line's territory, since lines are contiguous by construction): if divergent
  and the gap has opened past `EXTEND_THRESHOLD_RAD`, insert new nodes at target spacing --
  as many as it actually takes to close the gap (`dist / TARGET_LINE_SPACING_RAD`, capped at
  `MAX_EXTEND_NODES_PER_STEP` as a safety bound, not a normal limit), each given the
  ridge/rift target elevation directly (it's brand new material, not interpolated from
  anything). Inserting a fixed one node per step regardless of gap size used to be the rule;
  at a large `years` step (the UI offers up to 10 Myr) a fast-diverging boundary can open by
  many spacing units in a single step, and one node per step falls further behind every
  step -- that perpetually-reopening leftover gap looked to gap-filling (below) like a
  genuinely new, unclosable gap and kept spawning fresh micro-plates at the same busy
  boundary. If convergent and the gap has closed below `MERGE_THRESHOLD_RAD`, delete the end
  node (crust destroyed/folded away) -- but never the plate's last remaining node. This is
  where mass conservation actually lives: material is only ever created at divergent
  boundaries and destroyed at convergent ones, as literal point insertion/deletion -- there's
  no grid-resampling step that can lose or duplicate it.

<a id="line-regularization"></a>
## Line regularization (`line_regrid.py`)

Per-step boundary evolution only ever touches a line's two ends, so interior spacing stays
regular on its own during ordinary convergent/divergent motion -- but a *transform* boundary
shears nodes along a line without inserting or deleting anything, which can leave spacing
uneven. Every `REGULARIZE_INTERVAL_STEPS` calls to `step_world` (default 5), any line whose
gaps have drifted past `IRREGULARITY_TOLERANCE` (1.5x target spacing, either direction) gets
a fresh evenly-spaced node set across its *existing* extent -- the two endpoints are
preserved exactly, since this never changes where a line's physical edge is, only how
regularly it's sampled -- with elevation re-interpolated onto the new nodes (`np.interp`, 1D
since it's along a single already-ordered curve, not 2D scattered-data interpolation).

<a id="merge-and-split"></a>
## Merge and split (`merge_split.py`)

Unlike rotation and boundary evolution, these are rare, discrete, topology-changing events,
so a one-time resample is an acceptable cost here -- the exact, no-resampling guarantee only
matters for routine per-step motion.

- **Consumption.** A plate whose every elevation node has been deleted (fully subducted), or
  that's been eroded down to a single remaining line -- no real territory left, just a
  sliver along one latitude, regardless of how many nodes are still on that one line -- is
  simply dropped from `world.plates` (`remove_defunct_plates`). Falls directly out of the
  boundary-evolution rule above; no special algorithm needed for either case.
- **Continental collision merge.** If two continental plates have at least
  `MERGE_MIN_CONTACT_NODES` node pairs within `MERGE_CONTACT_DISTANCE_RAD` of each other
  *and* a real closing rate at those points (`boundary.closing_rate`, the same check
  boundary.py uses to classify a convergent boundary), they're fused: keep one plate's
  `frame`, and resample the union footprint from scratch -- a k-d tree over the pre-merge
  combined point cloud, with every candidate lattice node within `MERGE_COVERAGE_RADIUS_RAD`
  of *some* old node kept and given that old node's elevation
  (`plates.build_lines_from_lattice` again). The dropped plate's contribution is gone; the
  boundary between them becomes ordinary interior territory.

  **Why distance alone isn't enough.** plates.py's tiling has no gaps (every point belongs
  to exactly one plate at generation), so *every* pair of neighboring plates is already
  touching along their shared boundary the instant they're generated, regardless of whether
  that boundary is convergent, divergent, or transform. An earlier version checked distance
  only, and for roughly a quarter of random seeds this flagged ordinary neighbors as
  "colliding" before anything had moved at all -- a merge fired on literally the first
  simulation step no matter how small `years` was, because the trigger was a pre-existing
  generation-time condition, not anything that happened during the step. Requiring a genuine
  closing rate, not just proximity, is what actually distinguishes a real collision from any
  other pair of neighbors.

  **Why an instantaneous closing rate still isn't enough.** A curving boundary is commonly
  convergent along one stretch even while the two plates are mostly just sliding past each
  other overall -- a single step can clear the closing-rate check without the plates being
  in anything like a real collision. And real continental collisions play out over tens of
  millions of years, not one step, however small. So a pair only actually merges once
  they've been *continuously* close-and-converging for a sustained, randomized duration
  between `COLLISION_MERGE_MIN_YEARS` and `COLLISION_MERGE_MAX_YEARS` (50-100 Myr,
  randomized per pair via `_collision_threshold_years(seed, pair)` so unrelated collisions
  don't all resolve in lockstep) -- tracked in `World.collision_progress` (pair -> years
  accumulated so far) and advanced by `update_collision_progress` every step. A pair that
  stops being close-and-converging before reaching its threshold has its progress dropped
  entirely, not paused -- the collision didn't sustain. And even once one or more pairs
  *are* ready, `apply_topology_changes` merges at most one per step -- multiple simultaneous
  ready pairs don't all resolve in the same call, matching how every other change in this
  model happens incrementally rather than in a batch. Confirmed directly: scanning random
  seeds at a 1M-year step size, merges now take a realistic 50-100 Myr of sustained
  collision to resolve (verified against `world.events`, see below) rather than firing on
  the first or second step.

  **Event log.** `apply_topology_changes` returns a list of human-readable strings for
  outcomes that actually changed the world this step -- a merge completing (with how many
  million years it took), a plate disappearing (consumption, in either sense above), or a
  split creating a new plate -- which `world.step_world` timestamps with the *post-step*
  `elapsed_years` and appends to `World.events` (capped at `world.MAX_EVENT_LOG_LENGTH`
  entries). `generate_world` logs an initial "world generated" event the same way. A
  collision merely *starting* is deliberately not logged: plates.py's tiling has every
  neighbor pair already touching at generation, so a real fraction of pairs clear the
  proximity-and-closing-rate check in `find_continental_collision_pairs` at some point
  without ever accumulating the sustained duration needed to actually merge -- logging every
  one of those would flood the console with events that don't end up mattering. The API
  returns the full current log on every `/world/generate` and `/world/step` call (see
  api-reference.md) for the frontend's collapsible console.
- **Split.** Each plate's mantle-flow samples are clustered into two groups
  (`scipy.cluster.vq.kmeans2`, k=2) and a separate Euler pole is fit to each. If a single
  rigid rotation fits the whole plate poorly (RMS residual above
  `SPLIT_RMS_RESIDUAL_THRESHOLD`) *and* the two clusters' poles genuinely disagree (more
  than `SPLIT_MIN_POLE_SEPARATION` apart), the plate is cut along the great circle
  equidistant from the two clusters' centroids (`P . (centroid_a - centroid_b) == 0` is
  exactly that circle's plane) and every node partitioned by which side it falls on.

  **Why splitting needs a cooldown.** A single rigid rotation essentially never fits a wide
  footprint's flow samples *exactly* -- any spatially-varying field sampled over a large
  angular extent has some residual, even for a plate with no business splitting. Early
  versions tuned these thresholds against a small test scenario and found, at real scale,
  that ordinary large plates cleared them on every single step; worse, a freshly-split
  daughter plate -- cut from a continuous field, not a genuinely bimodal one -- would often
  still clear the thresholds on the very next step, recursively re-splitting into dozens of
  thin near-parallel slivers within a handful of steps (visually, this looked exactly like
  elevation "banding," since each sliver rotates almost identically to its neighbors).
  `SPLIT_MIN_AGE_STEPS` (a per-plate step counter, reset to 0 on creation by generation,
  split, or merge) requires a plate to exist for a while before it's split-eligible again,
  and `SPLIT_MIN_NODES = 1200` keeps the check off small fragments entirely.

<a id="gap-filling"></a>
## Whole-sphere coverage (gap-filling) (`gaps.py`)

Boundary evolution above only ever grows/shrinks a line's *theta*-direction ends -- a plate
can spread sideways along its existing rows, but it never gains a whole new row toward its
own local pole, and territory a fully-subducted neighbor vacated isn't automatically
reclaimed. Both leave literal gaps: sphere regions no plate currently covers.

Every `line_regrid.REGULARIZE_INTERVAL_STEPS` calls (the same cadence as line
regularization), `gaps.fill_gaps` sweeps a global lattice (`plates.iter_local_lattice` in the
identity frame, reused as a plain lat/lon sweep purely for this one-off detection query),
finds every candidate point farther than `COVERAGE_RADIUS_RAD` from any plate's nearest node,
and clusters the results (`scipy.sparse.csgraph.connected_components` over a k-d-tree radius
graph). Clusters smaller than `MIN_GAP_POINTS` are left alone -- ordinary growth lag that
boundary.py's per-line extension already closes on its own.

Each remaining cluster is resolved one of two ways:

- **A plate's border dominates it** (`DOMINANT_BORDER_FRACTION` of nearby existing nodes),
  **or a young plate has a meaningful share of it** (`age_steps <= YOUNG_PLATE_AGE_STEPS`,
  `YOUNG_PLATE_MIN_BORDER_FRACTION`): absorbed into that plate. Its line set is rebuilt from
  its own local lattice (`plates.build_lines_from_lattice`), preserving elevation wherever
  old data exists (nearest-neighbor lookup, the same technique `merge_split.merge_plates`
  uses) and giving newly-claimed area fresh ridge/rift elevation. This is what actually lets
  a plate grow toward its own pole or reclaim vacated territory.
- **No plate dominates and no young plate qualifies**: the cluster becomes a brand new
  plate -- new crust genuinely forming in open space between separating plates, not
  arbitrarily assigned to one side. `fill_gaps` returns one event message per newly spawned
  plate (`world.step_world` logs each to `World.events`); absorption isn't logged, since it
  grows a plate that already exists rather than adding or removing one.

**Keeping this gradual.** An absorb only claims the part of a gap within `GROWTH_RING_RAD`
of the plate's *existing* nodes (one ring per pass, not an entire possibly-huge cluster at
once -- e.g. an uncovered polar cap nobody's lines reach yet), and only up to
`MAX_ABSORB_NODES_PER_PLATE_PER_CALL` total per plate per call, across every gap it
dominates that call. Without these caps, whichever plate already has the longest border
(i.e. is already the biggest) dominates nearly every nearby gap and can absorb hundreds of
nodes in one call -- confirmed directly during development, where a large plate visibly
ballooned in a single step. Every other change in this model happens incrementally; these
caps keep gap-filling consistent with that instead of being the one place growth can jump.

**Why a young plate gets first claim.** A wide, long-lived, genuinely-shared rift would,
without this exception, spawn a brand new sliver plate at *every* pass -- a fresh spawn
never dominates its own border any better than the last one did, since neither side of a
symmetric spread is "winning." This was also confirmed directly: a busy boundary
fragmenting into a fan of thin, near-parallel micro-plates over a few hundred Myr, each
rotating almost identically to its neighbors. Letting a recently-created plate claim a
meaningful (not necessarily dominant) share of its own neighborhood for a few passes after
it's created breaks that chain -- it gets a chance to consolidate the rift it was born into
before being treated as just another equal competitor.

<a id="volcanism"></a>
## Volcanism (`volcanism.py`)

New continental crust forming where plates are separating -- run on the same
`steps_since_regularize`-gated cadence as gap-filling/line regularization above (detection),
plus every step (eruption and field-lifecycle bookkeeping, alongside erosion.py/bathymetry.py).

**Detection, two passes.** First, every elevation point's own nearest neighbor, whole-world,
unrestricted by plate -- the median of all these distances is "how far apart elevation points
normally sit." Second, only *boundary* points (each plate's own line endpoints,
`Plate.outline_world()`'s own definition of a plate's territory, reused directly) get checked
against each other: for each boundary point, the nearest boundary point on a *different*
plate. If that's more than `GAP_OUTLIER_FACTOR` (3x) the pass-1 median, it fires.

Restricting pass 2 to boundary points, rather than every point whole-world, turned out to be
essential, not a stylistic choice. An earlier version checked every point's own k=4 nearest
neighbors, whole-world, for the *same* point being both the density reference and the outlier
check -- and never fired, across dozens of seeds and step sizes, even with boundary.py's own
line-growth completely disabled. The reason: a point sitting right next to a genuinely wide
inter-plate gap still has plenty of *same-plate* interior neighbors much closer than that gap,
so an unrestricted whole-world nearest-neighbor query never sees the gap at all -- confirmed
directly by a separate cross-plate-specific distance measurement, which found real gaps up to
5x the typical spacing the whole-world query was blind to. Restricting pass 2's search to
boundary-vs-boundary removes the same-plate interior points that were masking the signal,
without changing what "normal spacing" means (pass 1 stays whole-world, since that's a stable
reference regardless of where the boundary happens to be).

Every qualifying pair contributes one new volcano point (the great-circle midpoint between
the two boundary points -- the actual empty space between the separating plates, not either
plate's own territory). All of this pass's new volcano points are then clustered by proximity
(`gaps.cluster_points`, the same connected-components technique gaps.py's own gap-clustering
already uses) and each cluster becomes one brand-new `Plate` with `crust_type="continental"`
("volcanic fields ... result in continental plates," per spec -- the rock has continental
physical properties regardless of whether the two separating plates were themselves oceanic
or continental). Every node of a freshly-spawned field starts as an active volcano.

`VOLCANIC_FIELD_CLUSTER_RADIUS_RAD` is deliberately much wider than gaps.py's own
`CLUSTER_RADIUS_RAD` (~1.5x line spacing, ~187km): gap-clustering there groups points from
one dense, contiguous coverage scan, but pass 2's candidate points here are individual
boundary-point pairs spread out along a whole divergent boundary's length -- at gaps.py's own
radius, a single long rift's many independently-qualifying points never merged into one
cluster at all, each spawning its own tiny field instead (confirmed directly: 11-40+ new
plates from a single clean-up pass on a 10-plate world, the same "boundary fragmenting into a
fan of micro-plates" failure shape gap-filling's own `MIN_GAP_POINTS`/young-plate exception
above were built to prevent). 15x line spacing (~1875km) merges same-rift detections into one
field while still keeping genuinely separate rifts elsewhere on the sphere apart.

A boundary point belonging to a plate that's *currently* a tracked volcanic field
(`World.volcanic_field_plate_ids`) is excluded as a pass-2 *source* candidate -- without this,
a field's own still-forming edge could immediately re-fire against the very neighbor it just
separated from, spawning another field on top of the last one every single clean-up pass.
It's still a valid *target* for some other plate's own check, so a genuinely separate rift on
the field's far side isn't blocked.

**A plate stops being tracked as a volcanic field once fewer than
`VOLCANO_FRACTION_DORMANT_THRESHOLD` (5%) of its own nodes are still `is_volcano`** -- not a
fixed elapsed-time countdown. `is_volcano` never reverts to False once set, so this ratio can
only ever fall, and only by dilution: as the field's own edges grow via ordinary boundary
evolution (or absorb gap territory via gaps.py), each newly-added node starts non-volcanic, so
a field that keeps growing eventually reads as "just an ordinary continental plate that
happens to have a few old volcanoes embedded in it." Checked every step, alongside the
eruption roll below, so the transition is caught within one step of crossing the threshold,
not lagged to the next clean-up interval.

**Eruption, every step.** Each individual volcano point has its own
`volcano_active_years_remaining` (`VOLCANO_ACTIVE_MIN/MAX_YEARS`, 100k-1M years, drawn once at
creation), decremented every step. While active, it rolls a per-step eruption chance
(`1 - exp(-ERUPTION_RATE_PER_MYR * active_years_this_step / 1e6)`, the same
exponential-arrival-rate shape used elsewhere in this codebase, e.g. lake evaporation's own
retention factor -- expected roughly 0.3 to 3 eruption events over a volcano's own full active
life, "occasionally," not every step) and, if it erupts, adds `ERUPTION_ELEVATION_M` (100m) of
new land. Deterministic per `(seed, elapsed_years, plate_id, line_index)`, the same
reproducibility precedent merge_split.py's own per-pair collision threshold sets.
`active_years_this_step` is clamped to the volcano's own *remaining* life, not the full step
size -- a large step (the UI offers up to 10 Myr) shouldn't roll eruption chances for years
past when a short-lived volcano actually went dormant.

**Rendering.** Baked directly into the elevation/plates views' raster the same way lakes and
glaciers are (`VOLCANO_COLOR_RGB`, a hot red-orange distinct from both), drawn after lake but
before glacier so ice still wins where both would apply (a volcano cold enough to glaciate
should read as ice-covered, not lava-red).

<a id="reassignment"></a>
## Boundary point reassignment (`reassign.py`)

Boundary evolution only ever grows or shrinks a line's two *ends* -- it never revisits
whether an *interior* node still actually belongs to the plate carrying it. Enough shearing
along a transform boundary, or a slow rotational drift, can leave a node geometrically
embedded in a neighboring plate's own territory while its data still lives on its original
plate's line. `reassign.reassign_misplaced_points` finds and fixes these.

For every node, it finds that node's `NEIGHBOR_COUNT` (4) nearest neighbors across the whole
world (not just its own plate). If at least `MIN_FOREIGN_NEIGHBORS` (3) of them belong to the
same other plate, the node is misplaced and gets moved:

- **Which line.** Among the *target* plate's lines, the closest one to the node is picked as
  the destination -- but instead of searching every line the target plate owns, it checks
  only the lines its own matching neighbors already sit on. Since those neighbors are, by
  construction, some of the node's closest points in the whole world, the line they're on is
  almost always the actual closest one too; this turns an O(lines-in-target-plate) search
  into an O(1-3) one.
- **Snapping onto it.** An `ElevationLine` is only ever defined at one fixed `phi`, so "on
  line A" and "phi equals A's phi" are the same statement. The node's local phi (in the
  target plate's frame) is overridden to A's phi -- its theta is left as computed, so this is
  a small shift, not a resample, consistent with the node having already been established as
  geometrically close to A.
- **Elevation.** Interpolated as a straight average of the node's own prior elevation and
  whichever existing node on line A sits nearest it in theta.

Every node's fate is decided from one snapshot (a single global k-d tree built once at the
start of the pass) and applied in two batched passes afterward -- all removals from source
lines, then all insertions into destination lines -- rather than mutating plates
mid-decision, so one node's move can't perturb another node's neighbor query within the same
pass. Not logged to `World.events` -- on a busy boundary this can touch many plates and
hundreds of points in a single pass, which flooded the console with little useful signal
even collapsed into one summary line per pass.

**Cadence.** Runs periodically (`REASSIGN_INTERVAL_STEPS`, default 5, its own counter
`World.steps_since_reassign`) but deliberately never the same step as `gaps.fill_gaps` (see
`world.step_world`): both are whole-sphere, cross-plate structural passes, and running them
together would mean each one's "which plate owns what" picture could already be stale by the
time it acts, from the other having just spawned, grown, or reassigned territory out from
under it.

<a id="projections"></a>
## Projections (`projections.py`)

[Behrmann](https://en.wikipedia.org/wiki/Cylindrical_equal-area_projection) (cylindrical
equal-area, standard parallel 30 degrees) and [Eckert
IV](https://en.wikipedia.org/wiki/Eckert_IV_projection) (pseudocylindrical equal-area,
[Snyder 1987](https://pubs.usgs.gov/pp/1395/report.pdf)) both have short closed-form/
Newton-iterated formulas, implemented directly in vectorized numpy rather than pulling in
`pyproj` -- and registered in `projections.PROJECTIONS` so adding a third projection later
is a one-function, one-dict-entry change. Both take/return radians and operate on a
unit-radius sphere; the frontend picks a single uniform pixel scale from the returned data's
bounding box (not independent x/y scales) so the equal-area property actually reads as
equal-area on screen.

<a id="render-image"></a>
## Render image (`render_image.py`)

`GET /world/render` returns a PNG, base64-encoded, rendered entirely server-side --
`render_image.py` computes the same projected coordinates and drawing rules
`MapCanvas.tsx` used to compute client-side from raw JSON, but paints them into a Pillow
image instead of shipping the coordinates over the wire for the browser to draw. The
client's job shrank to "decode this PNG and draw it on a canvas" (see `MapCanvas.tsx`);
every drawing decision -- fill colors, grid-cell sizing, boundary/pole/rotation-arc
overlays, per-node dots -- lives in one place now.

**Why a grid, still.** The elevation-line/gap-filling data is genuinely Lagrangian: nodes
are spaced at `TARGET_LINE_SPACING_KM` *within each plate's own local frame*, not on any
shared, screen-aligned grid. Drawing that raw point cloud directly -- one dot per node --
leaves visible gaps once projected: projected spacing isn't uniform (Behrmann, for
instance, stretches longitude spacing by roughly 50x at high latitudes relative to the
equator, since `x = lon * cos(30deg)` doesn't compensate for the shrinking circumference
the way the latitude term does), so a dot sized to look right at the equator leaves gaps
near the poles, and no single fixed dot size closes the gap everywhere without grossly
overlapping elsewhere. `_render_grid_arrays` fixes this the way the user originally asked:
sweep a uniform lat/lon grid over the whole sphere
(`plates.iter_local_lattice(identity_frame, spacing_rad=GRID_SPACING_RAD)` -- the same
identity-frame trick `gaps.py` uses for its coverage sweep), assign every cell its nearest
elevation node via one `cKDTree` query against every plate's current nodes (no distance
cutoff -- every cell gets *some* value, so there's no gap by construction regardless of how
far the nearest real data happens to be), and project the whole grid the same way
everything else is projected. This is a one-time resample purely for rendering -- it never
touches `world.plates`, so it has no bearing on the mass-conservation properties the rest
of the model is built around (see [Why not a grid](#why-not-a-grid)); a fresh grid is
computed from scratch on every `/world/render` call.

`GRID_SPACING_KM = 100` is a fixed, display-oriented constant, deliberately decoupled from
`plates.TARGET_LINE_SPACING_KM` (the simulation's physics resolution) -- the grid only
needs to look smooth once rasterized, not match whatever resolution the simulation itself
happens to run at. This used to matter for wire-payload size too (when the grid was
serialized as JSON, aliasing it to the physics resolution silently quadrupled the response
size the one time the physics resolution was doubled), but a PNG's size depends on how
compressible its *pixels* are, not on how many samples went into it, so that concern no
longer applies -- `GRID_SPACING_KM` is free to change independently now, and was tightened
from an initial 250km once the "Elevation"/"Plates" views' coastlines read as noticeably
blockier than "Plates (details)" (which draws at the physics resolution): confirmed side by
side that 100km closes most of that gap at a real render size (~220ms at 2200x1222, up from
~117ms), while 250km's cell edges visibly stair-stepped even there.

**Sizing each cell correctly is what actually closes the gaps.** A uniform *sphere* grid
still isn't a uniform *projected* grid -- so each cell is drawn at a size measured from the
projection's own local behavior, not a fixed pixel size: `_row_cell_half_extent` projects
two extra nearby samples per row (one step further in theta, one row further in phi) and
measures the resulting on-screen offset directly, giving that row's cell half-width and
half-height in the same projected units as the point coordinates. theta-only and phi-only
offsets are used deliberately so each measurement isolates one derivative exactly (in both
projections, moving in longitude at fixed latitude never perturbs the projected y
coordinate, and vice versa) rather than approximating a mixed partial derivative.
`render_png` fills each cell as a rectangle of that measured size (with a small
`CELL_OVERLAP_FACTOR` margin so adjacent cells overlap a hair rather than risk a hairline
gap from floating-point rounding), written directly into a numpy pixel array via array
slicing (`_fill_rects`) rather than one Pillow draw call per cell -- cheap even at tens of
thousands of cells, which "Plates (details)" needs for its per-node dots too.

**Resolution is a request parameter, not a fixed constant.** `width`/`height` in the
request become the PNG's exact pixel dimensions; the frontend requests more pixels than its
canvas's displayed CSS size (`App.tsx`'s `RENDER_SCALE`) for a sharper, retina-style render
at the same on-screen footprint. Fixed-pixel-size visual constants (`PADDING_PX`,
`POLE_RADIUS_PX`, boundary/arc line widths, node dot radius) are all defined relative to
`REFERENCE_WIDTH_PX = 1100` and scaled by `width / REFERENCE_WIDTH_PX` when drawing, so a
higher-resolution request doesn't also make those features look proportionally thinner.

**The pole marker is the true Euler pole, colored by plate.** `pole_xyz` is exactly
`plate.omega / |plate.omega|` -- no adjustment toward the plate's own territory. Real
rotation axes are frequently nowhere near the plate they belong to (this is physically
normal for real plate tectonics, not something to correct for), so the marker can land
anywhere on the map, including on top of, or far from, its own plate's fill. Since position
alone can no longer tell you which plate a given pole belongs to, `render_png` fills the
marker with that plate's own color (the same one used for its boundary and rotation arc)
with a white outline for contrast, rather than a single fixed color for every plate.

**The rotation arc replaces the old straight velocity arrow.** A fixed-pixel-radius arc is
drawn around the pole marker itself (`ARC_RADIUS_PX`), swept by an angle between
`ARC_MIN_EXTENT_DEG` and `ARC_MAX_EXTENT_DEG` scaled by how fast the plate is moving
relative to `mantle.MAX_PLATE_RATE`, with an arrowhead at the moving end
(`_draw_rotation_arc`). Because it's centered on an already-projected point and sized
entirely in pixel space, it can't straddle a projection discontinuity near the antimeridian
the way the old arrow (drawn between two separately-projected world points) occasionally
did. Its sweep direction can't be assumed from the sign of `omega` alone -- image y grows
down and projections aren't guaranteed to preserve on-screen handedness -- so it's measured
directly: take a point near the pole, find its true tangential velocity (`omega x point`,
the same formula `boundary.closing_rate` uses elsewhere), project both the point and a small
step along that velocity into pixel space, and see whether the angle around the pole (in
PIL's arc-angle convention) increased or decreased. Confirmed directly that two plates at
the same seed, differing only in the sign of `omega`, render as mirror-image arcs -- the
direction is genuinely sourced from the physics, not a fixed assumption.

**Every view has a legend**, anchored at the map's bottom-left corner -- rendered entirely
client-side (`frontend/src/legendData.ts`/`Legend.tsx`) as an ordinary HTML overlay, not baked
into the PNG the way it used to be (see git history for the old `_draw_legend`/Pillow
version). None of a legend's content is actually data-dependent -- the color stops and symbol
sets are fixed per view, not derived from world state -- so `legendData.ts` is a plain static
table keyed on `MapView`, hand-kept in sync with `render_image.py`'s own color constants (the
same relationship the deleted `elevationColor.ts` used to have, see that history, brought back
for just the legend). Moving it off the server fixed a real limitation: the PNG-baked legend
couldn't update while a rotation drag was in progress (see "Rotating the view" below), since
the live preview only redraws the *previous* rendered frame underneath a wireframe graticule,
never the legend baked into it.

<a id="rotating-the-view"></a>
## Rotating the view

Every map view can be reoriented interactively: press and hold on the map, then drag to pan
it -- simple linear panning (dragging N pixels always pans the view by the same angle,
calibrated so a feature under the cursor at the map's own center stays under the cursor as
you drag), not an arcball trackball (see git history for that earlier version, which allowed
a full 3-DOF roll/twist and had drag sensitivity that varied with distance from the canvas
center). Longitude wraps cyclically and latitude wraps *over* the pole, so dragging far in
one direction just keeps circling the globe rather than clamping at an edge. Release to
commit. This is a pure **view** transform, not a change to the simulated world -- see below
for why that distinction matters and how it's enforced.

**The transform itself.** A single 3x3 rotation matrix, applied to every *real* world
position immediately before it's projected to pixels (`_rotate`, `p @ view_rotation.T`) --
render grid cell positions, plate boundaries/poles/rotation-arc points, climate's whole grid,
wind/current arrows, swell markers, platesDetail's node dots. It travels from the frontend to
`GET /world/render` as a `rotation` query param (9 comma-separated floats, row-major, default
identity -- see `main.py`'s `_parse_view_rotation`) and is never stored on `World`: it's
client-local view state, sent fresh with every render call, exactly like `projection`/`view`
already were, not simulation state. This means it isn't reset by stepping or regenerating,
and different browser tabs can look at the same world from different orientations.

**Why this can never touch climate physics.** `climate.py`'s own grid (`_build_grid`'s
`lat_deg`/`world_xyz`) is the true, fixed planetary frame -- `compute_insolation`'s
zenith-angle law, `coriolis_parameter`, and the latitude-banded wind/current/humidity
structure all key off it directly. The view rotation is applied *only* at the final
projection step, strictly downstream of all of that, so climate simulation results are
completely unaffected by which orientation the user happens to be looking from -- rotating
the view is exactly as consequential as rotating a printed map on a table, never a change to
the planet's actual spin axis or sub-solar point.

**Cell coverage had to become genuinely per-cell.** `_render_grid_arrays` and
`_render_climate_view`'s cell-extent measurement used to take one representative sample per
row (correct at identity rotation, since these projections' scale only depends on latitude,
which is constant along a row) and broadcast it to every cell in that row. Once the view can
rotate arbitrarily, a row of constant *true* latitude can span wildly different *apparent*
positions depending on longitude, so cell size becomes a genuinely per-cell property.
Getting this right took two real fixes, both confirmed by direct visual debugging against a
render that showed pinhole gaps and, separately, huge false smears:

- **Measure by the cell's four corners, not its edge-midpoint neighbors.** A rotation can
  turn an axis-aligned cell into a rotated or skewed quadrilateral, and a rotated square's
  axis-aligned bounding box is set by its corners, not by how far its edge midpoints sit from
  its center -- edge-midpoint measurement underestimates the needed extent by a factor
  approaching sqrt(2) at a 45-degree rotation, which was producing widespread pinhole gaps.
- **Never measure a neighbor by reusing an adjacent cell's already-projected position** (e.g.
  via a `np.roll` index shift, which would otherwise be free). `xyz_to_latlon`'s `atan2` jumps
  from +pi to -pi at the antimeridian, and once the view can rotate arbitrarily that seam can
  land anywhere within a row or grid -- not a rare edge case, since it's a full line across
  the whole sphere that a dense grid sweep is guaranteed to cross somewhere. Reusing an
  arbitrary neighbor's position risks straddling that jump and measuring a wildly wrong
  extent (this was producing the huge false smears). The fix is `_project_offset`: like
  `_project_points`, but for a point known to be a small true angular step from another point
  whose longitude is already known -- it unwraps the new point's longitude to the
  representative within pi of that reference before projecting, which a genuinely small step
  should always permit.

Both fixes together cost a handful of extra full-grid projection passes (corners instead of
one sample per row) -- confirmed via direct benchmarking to add only tens of milliseconds at
real render resolution, comfortably inside the render budget (this is a deliberate mouse-up
action, not a per-frame interactive redraw, so it doesn't need to hit the same bar ordinary
renders do).

**The drag itself is client-side and cheap; the real render only happens on release.**
Re-rendering the actual detailed map (elevation fill, plate boundaries, climate heatmaps) on
every mouse move would be far too slow. `frontend/src/rotation.ts` ports just enough of
`geometry.py`/`projections.py`/`render_image.py` to TypeScript to compute the drag gesture
and preview it live:

- **Pan-to-center**: the drag tracks a target true (lat, lon) -- the point that should sit at
  the display center -- starting from wherever the committed rotation's own center already
  is (`centerOfRotation`) and offsetting it by the drag's pixel delta *from the drag's start*
  (not accumulated incrementally per-event, avoiding drift over a long drag), converted to an
  angle via a pixels-per-radian calibration measured from the projection's own Jacobian at the
  map's center (`getPixelsPerRadian`). `rotationForCenter` then builds the actual rotation
  matrix -- the shortest-arc rotation from the unrotated center to that target
  (`rotationBetween`), chosen specifically because it has no roll/twist about the target axis,
  so repeated small pans never accumulate spurious roll. Longitude wraps mod 2*pi; latitude
  "wraps" by reflecting back down past +-90 degrees with a 180-degree longitude flip, the same
  way walking a full meridian circle returns you to your start with no net longitude change
  (`wrapPanLatLon`).
- **Graticule preview**: while dragging, a wireframe of meridians/parallels every 30 degrees
  is rotated by the in-progress rotation, projected through the same `behrmann`/`eckert4`
  ports (confirmed numerically identical to the Python originals to full float precision),
  and drawn over the last real frame -- `MapCanvas.tsx` redraws the base image
  (`ctx.drawImage`, a cheap raster blit) then the graticule on top, on every mouse move.
- **The projected bounding box is rotation-invariant**, so the frontend's transform
  (scale/offset) can be computed once per `(projection, width, height)` and reused for every
  orientation: a full sphere's projected extent is identical regardless of rotation (rotation
  only permutes which physical point lands at which lat/lon, never changes the *set* of
  lat/lon values a full sphere covers). This is also why the graticule preview aligns
  pixel-perfectly with the still-displayed static frame underneath it -- both use the exact
  same transform.
- **The lat/lon-of-center readout live-updates during the drag** the same way the legend
  itself now does (see above) -- both are ordinary client-side state/DOM, not baked into
  anything the server sends, so neither needs a special-cased exception anymore the way the
  legend readout used to under the old PNG-baked legend.

**Interaction**: press and hold (`LONG_PRESS_MS`, with a small movement tolerance before that
so an accidental brush of the map doesn't start a rotation) shows the graticule at the
currently committed orientation; dragging updates it live; releasing commits the drag's
rotation and triggers the real server render, the same `refresh` path a projection or map
view change already used, just with `rotation` now also part of what's requested.

`MapCanvas.tsx`'s drag gesture itself was later extracted into a shared hook,
`frontend/src/rotationDrag.ts`'s `useRotationDrag`, once the Plate Inspector view (below)
needed the exact same long-press-then-drag-to-pan interaction to rotate a completely different
kind of content (translucent plate ellipses instead of a PNG + graticule). Callbacks are read
through refs updated every render rather than listed in the internal effect's own dependency
array -- deliberately: `onRotationPreview` typically drives a parent state update (the live
lat/lon readout) on every single mousemove, and including callbacks in that dependency array
would tear down and rebuild the drag's window-level listeners on nearly every frame of the
drag itself.

<a id="plate-inspector"></a>
## Plate Inspector

A second interactive map mode, structurally different from every other view: `GET
/world/plates` returns plate outlines and metadata as plain JSON (not a baked PNG), and
`frontend/src/PlateInspector.tsx` renders and drives the interaction entirely client-side --
reusing the same `rotationDrag.ts` gesture as `MapCanvas.tsx` (rotating still works
identically in both views, and the two share one lifted `rotation` state in `App.tsx`, so
switching between them preserves orientation), plus click-to-select and Tab/Shift+Tab to
cycle through plates.

**Every plate is drawn at all times** as a filled, translucent bounding ellipse (color from
`frontend/src/platePalette.ts`, a hand-synced port of `render_image.py`'s `PLATE_PALETTE` so
a plate's color matches its color in the "Plates" view) plus a thin stroke of its true
territorial outline -- confirmed with the user that this replaces having any elevation-image
backdrop at all: overlapping plates are meant to be visible where their ellipses blend. The
selected plate draws last, on top, at higher opacity with a crisp white outline.

**Every plate's individual node points are also drawn**, not just its outline/ellipse --
`GET /world/plates`'s `points` field is every node's own position (`Plate.
all_points_and_elevation()`'s points, not just the outline loop), rounded to 6 decimal places
before serializing (plenty for a unit sphere -- sub-meter at planet scale -- and it matters
now that a payload can carry tens of thousands of points per world, not just the much shorter
outline/ellipse loops). The selected plate's points are drawn last, in a highly visible white
(`SELECTED_POINT_ALPHA`); every other plate's points are drawn dim, in that plate's own color
at low alpha (`POINT_ALPHA`) -- confirmed with the user this dim/bright split is what
distinguishing "visible" from "less visible" colors meant. Unlike the outline/ellipse loops,
individual points don't need `wrapLongitudeNear`'s per-shape seam-unwrapping at all: each
point is an independent dot, not connected to its neighbors by an edge, so there's no
fill/stroke that could bow across the antimeridian -- the same reason `render_image.py`'s
"Plates (details)" view never needed that treatment for its own per-node dots either.
Measured directly: with points enabled, a full-world redraw (all plates, tens of thousands of
points, since every one is re-projected every frame) costs roughly 30-50ms depending on
projection (Eckert IV's forward projection needs Newton iteration per point; Behrmann is
closed-form and somewhat cheaper) -- noticeably heavier than before points existed, though
still comfortably interactive, not chosen to add throttling/level-of-detail for this since it
wasn't asked for.

**The bounding ellipse is a genuine minimum-area enclosing ellipse** (backend
`app/ellipse.py`'s `min_enclosing_ellipse`, Khachiyan's algorithm for the 2D
minimum-volume-enclosing-ellipsoid problem), not a bounding circle rendered as an ellipse by
projection distortion -- confirmed with the user this is what "rotated to fit all points as
closely as possible" meant, since a circle would only ever have one diameter, not two. A
containment-safety post-step (shrinking the fitted shape matrix just enough to cover whichever
input point Khachiyan's convergence tolerance left worst-covered) makes the guarantee
"contains every point" exact rather than merely "within `tol`" -- verified directly against a
closed-form oracle (a rectangle's true MVEE is analytically known) and against every node of
every plate in a real generated world.

Fitting happens in a **local flat-km coordinate system**, not naively in raw (phi, theta) or
directly on the sphere: `geometry.local_tangent_basis` + `azimuthal_equidistant_forward`
project the plate's points around its own `bounding_sphere` centroid into a plane where
radial distance from that center is *exact* real-world km (not exact between two arbitrary
non-center points, but fitting is always relative to that one shared center, so this doesn't
matter here). Fit against *every* node point, not `outline_world()` -- the minimum enclosing
ellipse of a full point set equals that of just its convex hull, and `outline_world()` isn't
a guaranteed hull for a concave plate. The fitted ellipse is sampled back into ~72 world-space
points (`azimuthal_equidistant_inverse`) and sent to the client already in true-frame 3D, same
as the plate's own outline -- the client only ever rotates and projects, never refits.

**Known, deliberately-accepted limitation**: `azimuthal_equidistant_forward` is numerically
singular for a point antipodal (or very close to it) to the projection's own center -- real
Voronoi-seeded plates are angularly compact around their own centroid, so this shouldn't
occur in practice, but a pathological non-convex or hemisphere-spanning point cloud could
trigger it. Not solving general spherical-MVEE for v1; flagged in code and covered by one
deliberately-adversarial test documenting current (imperfect but non-crashing) behavior
rather than leaving it silently untested.

**Click-to-select is a server round trip, not client-side point-in-polygon.** The client
unprojects the click (through whatever view rotation is currently active --
`rotation.ts`'s `unproject`, plus a `matTranspose` since a rotation matrix's inverse is its
own transpose) to a true (lat, lon) and asks `GET /world/plate_at`, which reuses the exact
same nearest-node `cKDTree` lookup the render grid already uses
(`plates.collect_all_points`/`nearest_plate_id`) to find the owning plate. This is more
robust than projected-polygon hit-testing would be (which risks the same antimeridian-seam
problems described above for filled/stroked shapes) and clicks are infrequent enough that one
round trip per click doesn't matter. Tab/Shift+Tab, by contrast, needs no server call at all
-- the full plate list is already in hand client-side, so cycling through it (sorted by
`plate_id`, wrapping both directions) is pure local state.

**Filling a shape that might cross the antimeridian is a genuinely new problem** existing
views didn't have: a *stroke* can simply skip one bad segment (`_stroke_robust_loop`'s
technique, already used for plate boundaries and the graticule), leaving an invisible gap,
but a naive *fill* across a seam-crossing loop would connect two far-apart points and shade a
wildly wrong region. Mitigated per shape: rotate the shape's own center first, use its
longitude as a reference, and unwrap every other point in that same loop to within it
(`rotation.ts`'s `wrapLongitudeNear`, a direct port of `render_image._project_offset`'s
technique) before projecting -- exact as long as the shape's own true angular radius from its
center stays under pi, the same scope boundary as the antipodal-singularity limitation above.

<a id="climate"></a>
## Climate (`climate.py`)

Seven fields -- land temperature, ocean surface temperature, air temperature, wind, ocean
currents, humidity, and precipitation -- computed for five map views (temperature, wind,
ocean currents with swells marked, humidity, precipitation), implementing a climate model
whose mechanisms are richer than a first-principles description would suggest, with the
caveat that mantle-bloom has no vegetation, rivers, or lakes.

**A third, genuinely fixed-shape grid, used only here.** Elevation is Lagrangian (see [Why
not a grid](#why-not-a-grid)); the render grid ([Render image](#render-image)) is a *ragged*
lat/lon sweep, immediately flattened to 1D. Neither supports the array tricks climate leans
on -- `np.roll` wraparound, centered-difference gradients, divergence, land-excluding
neighbor averaging. So climate gets its own equirectangular array, `lat: (H,)` / `lon: (W,)`,
every field `(H, W)`, `GRID_HEIGHT = 90` x `GRID_WIDTH = 180` (2 degrees/cell). It's never
stored on `World` and never touches `world.plates`.

**Fully stateless.** Every field is recomputed from scratch on every render call, from
whatever the *current* plate elevation/crust_type happens to be (the same `cKDTree`
nearest-neighbor sampling `_render_grid_arrays` already uses) -- `step_world` is completely
unaware this module exists; climate only runs when a climate view is actually requested. The
one exception: `World.axial_tilt_deg`, a fixed generation-time property like `seed` (set once
by `generate_world`, read again on every future render), since insolation needs it long after
generation and it isn't something to recompute per call.

**Pipeline order.** Wind needs a temperature field, but the *final* (current-advected) ocean
temperature needs currents, which need wind -- resolved by computing a pre-advection baseline
first and closing the loop only for the final consumer-facing fields:

1. **Insolation** -- `cos(lat)` zenith-angle law (clipped floor), plus axial tilt: with tilt
   0 it's the flat law; otherwise the mean of that same law over
   `AXIAL_TILT_DECLINATION_SAMPLES` declinations swept between `-tilt` and `+tilt` (the
   sub-solar latitude's annual sweep -- this model has no calendar, so it's an annual mean,
   not a season cycle).
2. **Land temperature** = f(insolation) plus elevation-based lapse-rate cooling, kept as part
   of the same base-heating formula rather than a separate causal channel (mountains being
   cold is a consequence of solar heating at altitude).
3. **Ocean temperature baseline** = f(insolation) only, a narrower range plus a freezing
   floor (water's greater thermal inertia) -- pre-advection.
4. **Wind** -- latitude-banded meridional flow (trade winds/westerlies/polar easterlies,
   empirical lookup, the near-surface branch of the real three-cell circulation) plus
   Coriolis zonal deflection (`u = GAIN * sin(lat) * v`), plus an additive term from the real
   local gradient of the pre-advection land/ocean surface temperature (`np.roll`
   centered-difference) -- the empirical banding supplies planetary-scale structure a
   gradient-alone field doesn't produce on its own, the
   gradient term is the genuinely temperature-responsive component. **Mountain
   deflection/Venturi/wake**: smooth elevation
   (Jacobi blur) before differencing into a gradient, cancel wind's into-slope component and
   redirect it tangentially with a speedup factor, ramped so gentle hills don't fully block;
   tangent side chosen by local topology (sample smoothed elevation along both perpendicular
   candidates, pick whichever is lower) -- mesoscale flow-splitting has no Coriolis-preferred
   side the way basin-scale currents do. Wake: walk backward along the (post-deflection) wind
   direction checking for upstream terrain, damping speed near the obstacle and relaxing back
   over a lookback window.
5. **Ocean currents** -- Ekman base (wind rotated by a fixed angle, hemisphere-flipped),
   redirected around coastlines by the same cancel/redirect mechanism as mountain deflection
   but with the tangent side chosen by a *fixed hemisphere sense* instead of local topology
   (real boundary currents have a Coriolis-preferred circulation direction), then smoothed
   along the coast (land-excluding Jacobi averaging, re-deflected each pass, so the effect
   propagates along a coastline rather than staying a single cell deep). **Land swirl**:
   every ocean cell gets a tangential contribution from its *nearest land cell* (`cKDTree`,
   nearest-cell-based rather than landmass-grouped), ramping from 0 at the coast to a peak then
   decaying with distance, direction matching the coastal deflection's hemisphere sense so
   the two agree. **Circumglobal boost**: a speedup on any row with zero land cells anywhere
   along it (a complete ocean ring) -- the stand-in for the Antarctic Circumpolar Current.
   **Wake**: the same backward-walk-and-damp structure as wind's, obstacle test is
   land-instead-of-elevation, plus a per-world-state noise texture (deterministic in `(seed,
   elapsed_years)`, so repeated renders of the same world state don't flicker) standing in
   for turbulent mixing.
6. **Ocean swells** -- convergence (negative divergence, `np.roll` centered differences) of
   the *final* current field, weighted-sampled down to `MAX_OCEAN_SWELLS` points
   (weighted-without-replacement sampling over the convergence field).
7. **Final ocean surface temperature** -- semi-Lagrangian backward advection (single
   fixed-distance backward sample, nearest-cell) of the baseline along the final current
   field: "carried by ocean currents."
8. **Air temperature** -- the land baseline's own solar-heating formula, pulled toward the
   *nearest ocean cell's* final temperature by a distance-based (`exp` e-folding) falloff:
   "moderating effect of oceans," literally, right down to the "nearest ocean and its
   temperature" query being a `cKDTree` chord-distance search -- true 3D distance, so no
   pole/antimeridian special-casing is needed the way a lat/lon-tangent-plane BFS search would
   require. Land temperature itself is never moderated -- only air temperature is.
9. **Humidity** -- an evaporation ceiling over ocean from the local final ocean temperature,
   advected onto land by a wind-driven 2D sweep: a zonal pass and a meridional pass (each a
   sequential flow-direction walk, vectorized across the perpendicular axis, single-column
   evaporation/retention/orographic-dump step per iteration -- pure numpy, no `numba`),
   blended per-cell by each wind component's share of total wind magnitude. The zonal pass's
   sweep *direction* is the same fixed latitude-band lookup wind's meridional structure uses
   (`zonal_direction_for_lat`), not the literal local wind sign. No evapotranspiration term
   (needs vegetation, which doesn't exist here -- an absent input, not a simplified one).
10. **Precipitation** = f(humidity) + an orographic bonus (continuous saturating
    windward-slope moisture dump, from wind blowing up-elevation) -- no zonal
    latitude-climatology baseline (equator/mid-latitude wet bands), cut deliberately.
    Feeds erosion and hydrology (see [Erosion](#erosion) and
    [Hydrology](#hydrology)) but still nothing else (no vegetation).

**Scope, explicitly decided.** Kept out: river outflow feeding currents (no rivers exist),
deep currents, precipitation's zonal climatology baseline. Included, even though richer than
a one-line causal description: axial tilt, wind's mountain deflection/Venturi/wake, ocean
currents' coastal deflection/land swirl/circumglobal boost/wake -- implemented in full, not
simplified down. Dropped outright, not reduced (their inputs don't exist in mantle-bloom):
humidity's evapotranspiration term, river outflow, lake climate influence.

**Rendering.** `render_image.py`'s `CLIMATE_VIEWS` (`temperature`, `wind`, `oceanCurrents`,
`humidity`, `precipitation`, `biome`) route to `_render_climate_view`, a separate path from
the plate-tectonics views since the data source (a real `(H, W)` array, always covering the
whole sphere) is structurally different from the render grid's ragged lattice. Heatmap views
(temperature/humidity/precipitation) reuse the elevation view's color-stop-interpolation
technique with their own stop tables; wind/ocean-currents draw subsampled arrows (numpy-
vectorized projection/direction math, looped only for the unavoidable per-arrow PIL draw
calls), and ocean currents additionally marks each sampled swell point with a small circle.
Temperature/humidity/precipitation additionally draw the current coastline (see
[Coastline](#coastline)) -- a color-scale view carries no land/ocean information on its own,
unlike elevation's own hypsometric coloring. `biome` is categorical, not a heatmap -- see
[Biomes](#biomes) below -- and (like wind/oceanCurrents) skips the separate coastline stroke,
since its own flat Ocean color already reads as a land/ocean boundary on its own.

<a id="biomes"></a>
### Biomes (`biomes.py`)

A pure, stateless classification -- `biomes.classify_biomes(temperature_c, precipitation_mm,
is_ocean)` -- bucketing each climate-grid cell into one of thirteen named biomes (Ocean, Ice,
Tundra, Boreal Forest, five temperate bands from desert to rainforest, and four tropical bands
from desert to rainforest) purely from two axes already computed by this module:
land-surface temperature and annual precipitation. Same two axes the real Whittaker biome
diagram uses, in the same broad cold-to-hot/dry-to-wet relative order, though with this
module's own boundary values -- a simplification in the same spirit as this codebase's other
openly-approximate constants (e.g. erosion.py's `RAIN_EROSION_COEFFICIENT`), not fit against
any specific real-world dataset. `ICE_TEMP_C` reuses `hydrology.GLACIER_ACCUMULATION_TEMP_C`
directly (rather than inventing a second, potentially-inconsistent cold cutoff) so a biome
map's Ice region lines up with where the simulation would actually grow a glacier. `is_ocean`
always wins over temperature/precipitation, since those are land-surface concepts and
`is_ocean` already settles the question for a water cell.

No new per-step state or caching -- like `render_image.py`'s own `temperature_colors`/
`humidity_colors`, this runs fresh from whatever `climate.compute_climate_cached` already
produced, entirely inside `_render_climate_view`'s `"biome"` branch. Implemented with
`np.select` (first-matching condition wins) rather than chained `np.where` overwrites, so
each temperature/precipitation band's cutoffs stay a self-contained, independently checkable
list instead of depending on write order to get boundary cells right.

<a id="erosion"></a>
## Erosion (`erosion.py`)

The other half of the weather<->geology coupling: [Climate](#climate) already has terrain
influencing weather (lapse-rate cooling, mountain wind deflection, orographic rain shadow).
This module is the new direction, weather influencing terrain, implementing an erosion model
cut down to the sources that don't depend on infrastructure mantle-bloom doesn't have.

**Scope cut to four erosion sources.** Coastal-current erosion is dropped (a
distinct source, never implemented here). Weathering's vegetation boost is dropped (no
vegetation field, same reasoning as climate.py's own "deliberately not ported" list).
Rain/sheet erosion, river-channelized erosion, weathering, and glacier erosion all feed into
a downstream deposition pass (see [Hydrology](#hydrology) for the flow-routing graph all of
this depends on, and [Glaciation](#glaciation) for how `glacier_depth` itself is
grown/melted/flowed), so material isn't purely one-way removed anymore: a slow, big river
drops part of its sediment load locally instead of carrying every last grain to the coast.
Glacier-driven **flattening** (broad terrain smoothing under an ice sheet) is a
mantle-bloom-original addition -- see below.

**The mapping problem, and why it's easier in this direction than the reverse.** climate.py
already solves node-cloud -> grid (`_sample_elevation_and_crust`'s cKDTree nearest-neighbor
resample) to build its grid from the current plate state. This module needs the reverse,
grid -> node-cloud: a node's world position converts straight to (lat, lon)
(`geometry.xyz_to_latlon`) and then to a climate-grid (row, col) by direct arithmetic --
mirroring `climate._build_grid`'s own convention exactly -- no tree, no resampling needed at
all, since (unlike the irregular node cloud) the climate grid is already a plain regular
lattice.

**Slope is the one genuinely new piece of math.** climate.py's grid gets slope for free
from neighbor-index differences; an irregular node cloud has no such structure. This reuses
`reassign.py`'s whole-world cKDTree pattern (build once, query `SLOPE_NEIGHBOR_COUNT=4`
nearest neighbors per node) instead: for each node, the elevation drop to the *lowest* of
its nearest neighbors (0 if the node is already a local minimum -- the "slope to lowest
neighbor" definition), divided by the real great-circle distance to that
neighbor. This is a genuine dimensionless rise/run -- elevation drop over real distance, not
elevation drop per grid step (a grid-step measure isn't a true slope at all, since it
silently depends on grid resolution) -- which is why `RAIN_EROSION_COEFFICIENT` was picked
by order-of-magnitude reasoning against `boundary.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR`
(800 m/Myr) and checked against real slope/precipitation distributions from an actual run,
rather than reused from a value tuned against a different scale.
`WEATHERING_COEFFICIENT` needed no such rescaling, since wind speed already uses a directly
comparable scale (`MERIDIONAL_BASE_SPEED = 6.0`) unaffected by the slope re-derivation. River
erosion needs the same re-derivation, one level further removed: it depends on `water_accum`
(see [Hydrology](#hydrology)), itself downstream-accumulated *precipitation* rather than a raw
grid-cell count (a cell-count accumulation would implicitly scale with grid resolution),
so `RIVER_EROSION_COEFFICIENT` was re-derived the same order-of-magnitude way, then checked
against real channel-depth growth over an actual run (see below). Glacier erosion needs no
such re-derivation either -- `glacier_depth` is driven by temperature and precipitation,
which already use consistent physical units (unlike slope), so `GLACIER_EROSION_COEFFICIENT`
needed no rescaling.

**The four formulas**, all computed per-node:

- **Rain/sheet erosion** = `RAIN_EROSION_COEFFICIENT * slope * (precipitation_mm / 1000) *
  dt_myr`.
- **River erosion** = `RIVER_EROSION_COEFFICIENT * channel_boost * water_accum_m^
  RIVER_FLOW_EXPONENT * slope^RIVER_SLOPE_EXPONENT * dt_myr`, where `channel_boost = 1 +
  CHANNEL_EROSION_BOOST * clip(channel_depth / CHANNEL_BOOST_REFERENCE_M, 0, 1)` -- a river
  preferentially re-carves its own established channel, the real mechanism behind
  meandering rivers staying within their valleys rather than cutting a fresh path every
  step. `channel_depth` (persistent, see `plates.ElevationLine`) grows by this same term
  every step (`clip(channel_depth + river_erosion_amount, 0, MAX_CHANNEL_DEPTH_M)`) --
  monotonically non-decreasing, land-only, capped purely as a sanity bound (not a
  physically-derived limit).
- **Weathering** = `WEATHERING_COEFFICIENT * wind_speed * humidity_norm * relief_factor *
  dt_myr`, where `humidity_norm` is humidity clipped to `[0, 1]` against
  `HUMIDITY_REFERENCE`, and `relief_factor = clip(slope / WEATHERING_RELIEF_REFERENCE_SLOPE,
  0, 1)` (same normalize-and-saturate idiom as `humidity_norm`). `relief_factor` is a
  mantle-bloom-original addition -- added after a real bug: unlike rain/river erosion (both already `*
  slope`), weathering originally had no relief dependence at all, so a flat, low-lying coastal
  plain weathered exactly as
  fast as a steep mountain flank -- confirmed directly against a real generated world (a
  single 2 Myr step converted ~4.9% of all land nodes to ocean; at the affected nodes,
  weathering alone -- uniform regardless of relief -- exceeded the node's *entire* elevation,
  well beyond rain/river's own, already slope-gated contribution). `WEATHERING_RELIEF_REFERENCE_SLOPE` (the slope at which weathering reaches full
  strength) was picked against a real world's own slope distribution: median land slope is
  roughly 2.7e-4, far below the reference, so typical flat terrain now gets only a small
  fraction of full weathering, while land already near the 90th percentile of steepness
  (roughly 3e-3) reaches full strength. This roughly halved the single-step land-to-ocean
  flip rate in the same test; the remainder is rain erosion (already slope-gated, but tuned
  with a large coefficient to balance mountain-building rates) eroding low-lying land with
  nonzero slope -- a separate, not-yet-addressed contributor to the same overall trend.
- **Glacier erosion** = `GLACIER_EROSION_COEFFICIENT * slope * ice_factor * dt_myr`, where
  `ice_factor = clip(glacier_depth / GLACIER_EROSION_REFERENCE_DEPTH_M, 0,
  GLACIER_EROSION_MAX_FACTOR)` -- driven by the node's own *actual accumulated ice depth*
  (a real persistent field, see [Glaciation](#glaciation)), not a stateless cold proxy;
  `depth * slope` approximates basal shear stress, a standard real glacial-erosion proxy: a
  flat-bottomed accumulation bowl still correctly erodes near zero regardless of how thick
  the ice sitting in it is. Uses the *previous* step's `glacier_depth` (same one-step lag
  `channel_boost` above already uses for `channel_depth`), since this step's fresh value
  isn't computed until `hydrology.compute_hydrology` runs, just before this term does.
- All four summed, zeroed over ocean nodes (`elevation <= 0`, the sea-level convention used
  everywhere else) -- every source here is a subaerial process. The combined result is
  capped at the node's own drop-to-lowest-neighbor (in meters, not the normalized slope), so
  a single step can't erode a node *past* the valley it drains into and carve a new, lower
  pit.

**Glacier flattening** (`_flatten`, mantle-bloom-original): real
continental ice sheets grind down local relief over broad areas (the Canadian Shield and
Fennoscandia read as glacially smoothed bedrock today, not just eroded lower) -- a genuine
local blur, not a directional erosion/deposition term, so it's applied as a separate signed
elevation delta rather than folded into `erosion_amount`. Each node relaxes toward the mean
elevation of its own `hydrology.py` flow-graph neighbors (reusing that graph rather than a
separate query), scaled by `GLACIER_FLATTEN_RATE_PER_MYR` and the same `ice_factor` glacier
erosion uses -- glacier-free nodes (`ice_factor = 0`) are completely untouched. Confirmed
directly on a real run: near-zero delta at nodes with little local relief, tens to 100+
meters at nodes combining real local relief with thick ice, consistent with "smooths sharp
terrain under ice, leaves already-flat terrain alone."

**Deposition.** The capped, combined erosion amount is routed downstream (see
`hydrology.route_downstream`) with a `retain_fraction` wherever a river qualifies as "big and
slow" -- `river_speed < DEPOSITION_SPEED_THRESHOLD` (a stylized, unitless quantity, faster
where slope is steeper and where more water has accumulated) *and*
`water_accum_m > DEPOSITION_MIN_FLOW_M` (guards against a merely-flat trickle depositing) --
in which case `DEPOSITION_FRACTION` of the material passing through that node settles right
there (a floodplain or delta) instead of continuing on; the rest keeps going, eventually
reaching either another depositing node, an internal sink, or the coast (where it can raise
an *ocean* node's own elevation -- a river delta building outward is real, intended
behavior, not a bug). `line.elevation`, `line.channel_depth` (this module's own), and
`line.lake_depth`/`line.glacier_depth`/`line.silt_depth` (state transitions owned by
`hydrology.py`/`lakes.py`, read directly from `World.hydrology_cache` -- see
[Hydrology](#hydrology)) all get written back
together per line -- no resampling, no topology change, so none of this can interact with
line regularization or point reassignment at all.

**Cadence: every step, no lag on climate -- but a deliberate change from erosion's own
earlier no-hydrology version regarding flow routing.** This module still calls
`climate.compute_climate(world)` fresh every step (no staleness to reason about, same as
before). Flow routing (`hydrology.compute_hydrology`) is comparatively expensive with no JIT
available, so rather than computing it twice per step, this module computes it once and
reuses the result for both erosion and `World.hydrology_cache` (see
[Hydrology](#hydrology)). Runs in `world.step_world`
right after boundary evolution and topology changes, every step (not gated by the periodic
regularize/gap-fill/reassign cadence).

<a id="bathymetry"></a>
## Bathymetry (`bathymetry.py`)

Nothing else pulls submerged continental crust toward any particular depth once it goes
underwater -- generation-time noise can put a continental node a few hundred meters below
sea level, and rifting (see [Boundary evolution](#boundary-evolution)) can push a node deep
underwater without any further constraint on where it settles. This module is a slow
background relaxation (same exponential-toward-a-target style as boundary.py's own divergent
relaxation) pulling every submerged (`elevation <= 0`) continental node toward one of two
targets, chosen by distance to the nearest land node (`elevation > 0`, any plate -- this is a
geographic coastline-proximity question, not a plate-boundary one):

- Within `SHELF_RANGE_RAD` (200km): `SHELF_TARGET_M` (-100m) -- the continental shelf.
- Beyond it: `DEEP_CONTINENTAL_TARGET_M` (-3000m) -- genuinely deep water, though still
  shallower than oceanic crust's own abyssal depth (`plates.BASE_OCEANIC_M = -3800`).

Deliberately continental-only: oceanic crust's average depth already comes from its
generation-time baseline and nothing erodes or otherwise drifts it away on its own
(erosion.py explicitly excludes ocean nodes from both its sources), so it doesn't need a
parallel correction. Relaxes at `BATHYMETRY_RELAX_RATE_PER_MYR` (0.3, slower than
`boundary.DIVERGENT_RELAX_RATE_PER_MYR`'s 0.5 -- a passive equilibration of already-
submerged, non-actively-deforming crust, not an active tectonic process). Runs every step,
right after erosion.

Confirmed directly on a 60 Myr run: submerged continental nodes within 200km of land
averaged -400m (still relaxing toward -100m -- shoreline nodes keep moving as coastlines
shift, so they rarely reach full equilibrium), nodes beyond 200km averaged -2869m (close to
the -3000m target, since deep-water nodes are disturbed far less often); rendered, this
shows up as a visibly lighter shelf band hugging every coastline.

<a id="hydrology"></a>
## Hydrology: rivers and lakes (`hydrology.py`)

This module implements hydrology over mantle-bloom's irregular per-plate node cloud rather
than a fixed grid. Three core algorithms -- steepest-descent flow
direction, priority-flood basin-spill (lake/depression detection), and elevation-ordered
downstream flow accumulation -- all turn out not to actually need a *grid*, only a *graph*:
a regular 8-neighbor grid adjacency is just one convenient substrate for them. This module
builds a graph instead, via a whole-world k-nearest-neighbor query (`FLOW_NEIGHBOR_COUNT =
8`, the same technique `erosion.py`/`reassign.py` already use for their own whole-world
passes), then runs the same three algorithms directly on it.

**Persistence comes for free.** A fixed-grid tectonic model has plates moving relative to
its grid, so a persistent field like `channel_depth` needs deliberate semi-Lagrangian
advection every step just to keep following the crust. mantle-bloom's elevation-line nodes
already rotate exactly with their own plate, so `channel_depth`, `channel_width`, `lake_depth`,
`silt_depth` (see [Lakes](#lakes-are-an-explicit-tree) below), and
`glacier_depth` (see [Glaciation](#glaciation)), stored as ordinary parallel arrays on
`ElevationLine` right alongside `elevation` itself (see
[Why not a grid](#why-not-a-grid)), get that same "just works" persistence for free -- no
advection scheme needed, since rotating a plate never touches those arrays at all. Making
this persistence real required threading every one of these fields through every place an
`ElevationLine` gets rebuilt: preserved unchanged where a rebuild doesn't change node
identity (boundary elevation deltas, bathymetry -- via `dataclasses.replace`, not an
explicit field-by-field reconstruction, specifically so a *future* persistent field is
preserved automatically rather than needing every such call site updated by hand; see
`plates.ElevationLine`'s own docstring for the concrete bug this replaced -- `is_volcano`
silently reset to `False` every step at exactly these two sites, for several steps of actual
development, before it was caught), sliced/concatenated to match where nodes are added or
removed (boundary growth/shrink -- new nodes start at 0, no history to carry; a merge/split's
boolean-mask slice), interpolated alongside elevation where a line gets resampled onto a
fresh spacing (`line_regrid.regularize_line`, since that runs periodically throughout the
simulation and a plain reset would erase rivers/glaciers constantly, not rarely). Two call
sites *do* deliberately reset to 0 rather than preserve: `plates.build_lines_from_lattice`
(generation, gap-fill spawn/absorb, plate merge -- genuinely new or wholesale-resampled
territory has no history to carry) and a reassigned point in `reassign.py` (a rare,
small-scale pass touching only a few boundary-adjacent nodes at a time, unlike
`line_regrid.py`'s near-every-line-every-interval reach).

`flow_target`/`flow_accum`/`river_speed` are deliberately *not* persisted: recomputed fresh
every step, from that step's real climate -- purely this-step derived quantities, cached on
`World.hydrology_cache` only so a later same-turn caller (rendering) doesn't recompute them
again, the same reuse pattern `climate.compute_climate_cached` already established.

**Flow routing is computed once, not twice.** It would be simple enough to compute flow
routing twice per step -- once in `hydrology.py` for the real river_flow/rendering fields,
again inside `erosion.py` for the `water_accum` erosion itself needs -- but flow routing here
has no JIT and is comparatively expensive (confirmed: ~70-100ms on a world with a few
thousand land nodes, a real chunk of a step's total cost), so `erosion.py` computes it once
(`hydrology.compute_hydrology`) and reuses the result for both erosion and
`World.hydrology_cache`, rather than paying for it twice.

**The three algorithms**, all operating on the k-NN graph:

- **Basin-spill** (`_compute_basin_spill`): a multi-source Dijkstra seeded
  from every ocean node, relaxing by *max* (the highest point a path is forced to cross)
  rather than by sum -- a minimax path cost, not a shortest path. Nested sub-basin chains
  collapse to one hop each, cycle-free by construction -- necessary because a naive
  single-neighbor check can't see past more than one nested rim. Still used for flow
  *routing* (`_compute_flow_direction`'s `should_spill`, below) -- lake *detection* itself now
  lives in `lakes.py`'s own, different (nested-basin-aware) algorithm, see
  [Lakes](#lakes-are-an-explicit-tree) below.
- **Flow direction** (`_compute_flow_direction`): among each node's k nearest neighbors
  strictly below its own elevation (a downhill candidate), prefers whichever one already has
  the deepest established channel (`channel_depth > CHANNEL_PREFERENCE_THRESHOLD_M`),
  falling back to plain steepest descent only when *no* downhill candidate has a real channel
  yet -- real rivers meander within, and stay inside, their own existing valley rather than
  recalculating the mathematically steepest path every step. This is a mantle-bloom addition:
  distinct from (and upstream of) `erosion.py`'s own `channel_boost`,
  which only affects how *fast* a node erodes once water is already routed there, not *where*
  it gets routed to begin with. A sink is -1 if there's no downhill candidate at all -- unless
  the sink's current water surface (`elevation + previous step's lake_depth`, a one-step-lagged
  "memory") has already reached its basin's
  true spill point, in which case it redirects to `spill_target` instead of staying a
  dead-end sink forever (turning a filled-past-its-rim lake into a normal through-flowing
  river cell).
- **Downstream accumulation** (`route_downstream`, public -- `erosion.py` reuses it
  directly): a single forward sweep over land nodes in elevation-descending order,
  accumulating a source quantity (precipitation, for `flow_accum`; eroded material, for
  erosion's own deposition pass) along `flow_target` edges, weighted by a `retain_fraction`
  per edge. Correct in one pass because
  every node's target is guaranteed strictly lower, so it's always visited *later* in this
  same order. There is no separate `loss_fraction` term for in-transit river evaporation --
  not needed, and this module already drops temperature-driven effects on
  hydrology to match erosion.py's own "precipitation is enough" simplification.

<a id="lakes-are-an-explicit-tree"></a>
**Lakes are an explicit object tree (`lakes.py`, not `hydrology.py`)**, replacing an earlier
flat per-node flood-fill design (`update_lakes`/`_flood_fill_lake_extent`, since removed).
That design fixed a real bug of its own (letting every previously-wet node, not just the
literal sink, carry its own depth forward and become a flood-fill seed -- see the surrounding
history preserved in `hydrology.py`'s own module docstring), but even with two added rate
limits (a hop budget on growth into new territory, a hard cap on per-node depth change per
step) it only bounded how *fast* a lake's flood-fill could grow, not whether it ever reached a
real equilibrium: a sufficiently large, gently-sloped catchment could keep creeping outward
step after step without ever stopping. `lakes.py` instead computes each basin's actual
floor/rim geometry directly, so a lake's maximum extent is known immediately rather than
discovered incrementally, which is what actually fixes that runaway-growth failure mode.

The core data structure is `lakes.Lake`: an n-ary tree built by `build_lake_hierarchy`, a
two-phase algorithm distinct from `_compute_basin_spill` above (which only finds one
component's own bottleneck to the ocean, with no notion of *nested* sub-basins merging with
each other first) -- a "depression hierarchy" / "watershed by immersion" technique (Barnes et
al.'s fill-spill-merge family), adapted from a regular grid to this codebase's k-NN graph:

- **Phase 1** assigns every node to the catchment of the true local minimum its own *pure*
  steepest descent (not `_compute_flow_direction`'s channel-biased version -- a routing
  concern, orthogonal to basin geometry) eventually drains to, or marks it as draining
  straight to the ocean if that chain never passes through a land local minimum first. Without
  this phase, a naive merge over every raw graph edge would wrongly treat every intermediate
  ridge/rim node passed on the way down to the ocean as its own trivial one-node "lake"
  (confirmed directly during development). Every genuine catchment gets a leaf `Lake` eagerly,
  not lazily on first merge -- a catchment that never merges with anything (a real, fully
  enclosed endorheic basin with no ocean anywhere reachable) would otherwise never get a `Lake`
  object at all, despite being exactly the "no known spill, could fill indefinitely" case this
  needs to represent.
- **Phase 2** walks catchment-boundary edges (weighted by `max(elevation[i], elevation[j])`,
  the saddle/col height a path between them is forced to cross) in ascending order,
  union-finding catchments together -- Kruskal's minimum-spanning-tree construction, except
  the "tree" being built is the *merge tree* of the elevation field itself. A lake's own
  `max_depth` is just the elevation of the very next merge event it takes part in, whether that
  unions it with a sibling lake (creating a parent, whose own `min_depth` is that same
  elevation -- "the point where the basins split") or with the ocean (finalizing it as a root)
  -- both are "reaching the top of the basin," one spilling to another basin, the other to the
  sea.

Each step, `lakes.step_lakes` rebuilds this hierarchy from scratch (`elevation +
prev_silt_depth`, not bare `elevation` -- see below) -- consistent with `flow_target`/
`flow_accum` above being recomputed fresh every step rather than advected -- then resolves
every lake's own water balance top-down: evaporate `prev_level`'s depth (same
`LAKE_EVAPORATION_RATE_PER_MYR`/`LAKE_EVAPORATION_BASELINE_M_PER_MYR` constants and tuning
rationale the old per-node design already used), grow from this step's own inflow (the sum of
`route_downstream`'s `water_deposited` over every one of the lake's own members, since more
than one can be a true sink once several originally separate basins have merged), spread as a
level rise over the lake's own member count (the same area proxy `_lake_component_sizes`/
`LAKE_BREACH_EROSION_COEFFICIENT` already use), clipped to `[floor_elevation, max_depth]`.
Tracking a lake's level as one scalar shared by every member, rather than a per-node depth, is
what removes the old design's threshold-crossing failure mode entirely. Merge/split
transitions fall out of comparing this against each lake's own `min_depth`: a still-separate
child whose own new level reaches its parent's saddle promotes to one merged body (pinned at
exactly the saddle, continuity); an already-merged lake whose new level recedes below its own
saddle demotes back into its children (both starting at exactly the saddle, continuity in
reverse). No cross-step `Lake`-object registry is needed for any of this -- see
`lakes.py`'s own module docstring for why continuity is instead derived from the same
persisted `lake_depth` array (`elevation + prev_lake_depth`, the same "last step's water
surface" the old per-node design itself already used) rather than matching this step's
freshly-rebuilt tree against last step's by member-overlap, which would be fragile to a tree's
exact shape shifting from ordinary terrain churn even when nothing physical about the lake
changed.

**Lakes accumulate silt.** `silt_depth` (a new persistent per-node array, threaded through
`boundary.py`/`line_regrid.py`/`reassign.py`/`merge_split.py` exactly like `lake_depth`) is a
small, ~100x-slower-than-water-growth fraction of the same inflow, settling permanently
(monotonically -- silt never erodes back away) on a lake's own bed. `build_lake_hierarchy`
is given `elevation + silt_depth`, not bare elevation, so a lake's own floor rises as silt
accumulates without touching the real terrain `elevation` other modules read -- a small,
low-inflow lake can plausibly silt in entirely over a long enough run: once its floor has
risen to meet the surrounding rim, it stops registering as a local minimum at all (steepest
descent no longer sees a depression there), and the lake disappears outright, exactly the
"reaches ground level" case a lake should eventually hit.

A river's own `flow_target` can point at a node that's genuinely part of a lake (real inflow),
but a flooded node is never itself classified `is_river` -- checked against this step's
*final* lake extent, after `lakes.step_lakes` runs, not the flat land/ocean split alone -- so a
river's own classification ends at the lake's shore rather than jumping straight across the
water to whatever's on the far side, even though the water itself (`flow_target`/`flow_accum`)
still physically continues through if the lake is spilling.

**A spilling lake erodes its own outlet hard**, proportional to the lake's own surface area,
not just the ordinary precipitation-driven flow passing through that one channel node. Once a
sink's `should_spill` fires (see `_compute_flow_direction`), `compute_hydrology` adds an extra
term into `water_source` right there -- `LAKE_BREACH_EROSION_COEFFICIENT` times the node count
of that lake's own connected component (`_lake_component_sizes`, a union-find grouping over
the same k-NN graph flow routing uses; no separate "physical area" concept exists elsewhere in
this codebase, so node count is the same implicit area proxy `flow_accum` itself already
relies on). That extra term then rides `route_downstream` downstream through the whole outflow
channel exactly like real precipitation would, so it's already fully picked up by erosion.py's
existing river-erosion and channel-growth formulas with no changes needed there -- both
channel_depth and channel_width grow from it automatically. This models a real phenomenon a
plain per-node flow_accum can't: a large standing lake finding an outlet carries far more
erosive force at that point than an equivalent ordinary river of the same instantaneous
discharge, because the whole reservoir sits behind a narrow breach rather than just whatever
rain fell immediately upstream this step -- confirmed directly against a real generated world
that a several-hundred-node lake's breach source (in the low millions) dwarfs the handful of
millimeters of ordinary precipitation reaching that same point, carving its spillway down
hard rather than leaving a shallow, precipitation-scale groove.

**Channel width** (`channel_width`, grown in `erosion.py` alongside `channel_depth`) is a
mantle-bloom addition: standard hydraulic-geometry scaling (width ~
discharge^0.5, the same discharge exponent `RIVER_EROSION_COEFFICIENT`'s own stream-power law
uses) but with no slope term at all -- width is driven by how much water passes through, not
by gradient, so a slow, wide lowland river and a narrow, steep mountain torrent can carry
comparable discharge without needing comparable width. Persistent, monotonically non-
decreasing, capped at `MAX_CHANNEL_WIDTH_M`, the same shape `channel_depth` already has,
rather than a stateless function of the current step's flow alone -- a river that temporarily
dries up doesn't instantly narrow either. Not yet surfaced anywhere in rendering or the River
Inspector's own per-river stats (`hydrology.RiverInfo`) -- tracked, but not yet exposed.

**Rendering.** Both baked directly into the existing elevation/plates/platesDetail views
(no new map-view mode or API surface needed): lakes are always visible, with no separate
overlay toggle needed. Lake cells use
a muddier, less-saturated blue than open ocean (`render_image.LAKE_COLOR_RGB`) via the same
nearest-neighbor grid resample `_render_grid_arrays` already does for elevation/plate_id
(`plates.collect_all_lake_depth`, index-aligned with `plates.collect_all_points`'s own
output so both can share one `cKDTree` query). Rivers are drawn as short line segments, one
per `is_river` node to its own `flow_target` (top `RIVER_FLOW_PERCENTILE = 90.0` of land
`flow_accum`) -- fixed color and
width for every segment regardless of discharge (only *which* segments get drawn varies with
flow magnitude, not how they're drawn).
A second, independent cut applies only here, on top of `is_river`: `render_image.
RIVER_DRAW_MIN_FLOW = 100,000` requires a segment's own `flow_accum` to also clear that
absolute floor before it's drawn on these general-purpose views, so a large world's merely-
top-decile trickles don't clutter every view -- the River Inspector (below) deliberately
keeps listing/drawing every `is_river` network regardless of flow, unaffected by this floor,
since picking a minor tributary out of the full list is exactly what that view is for.

Confirmed live on a real run: river networks render as visibly branching, dendritic
drainage patterns converging toward coasts and lake basins, matching real-world drainage
network shapes; channel_depth grows from 0 to several hundred meters over tens of Myr
without instantly saturating at `MAX_CHANNEL_DEPTH_M`; lakes form, merge, split, and
fluctuate in count/depth across many steps rather than either vanishing or ratcheting
monotonically upward, staying a small fraction of a world's total land nodes even after many
steps rather than creeping to cover a whole continent; `is_river & is_lake` never overlaps on
a real stepped world (a river's own classification genuinely stops at a lake's shore).

<a id="glaciation"></a>
## Glaciation (`hydrology.py`)

A node colder than `GLACIER_ACCUMULATION_TEMP_C` (-10C -- this model has no seasons, so a
mean annual temperature just barely below 0C represents a place with seasonal snow that
would fully melt over a real year, not permanent glaciation) permanently accumulates ice
instead of ever holding liquid water: any precipitation there is treated as fully frozen --
no partial liquid/frozen split, matching `erosion.py`'s existing "use precipitation is
enough" simplification, just gated by the same accumulation threshold rather than dropped
entirely -- and an existing lake sitting there freezes solid into `glacier_depth` this same
step, *before* `flow_target` is (re)computed, so a lake that just froze is correctly treated
as a genuine sink again this step rather than immediately re-filling from this same step's
routed water (see `lakes.step_lakes`'s own docstring). This ordering avoids a specific failure
mode: without the freeze-before-routing ordering and an `is_accumulating` gate on
`lakes.step_lakes`'s own inflow/pin-at-cap logic, a lake formed in a warmer epoch would never
freeze, and a just-frozen basin would re-flood back to its old cap the very same step,
double-counting the same water as both a lake and a glacier at once.

- **Accumulation**: `GLACIER_ACCUMULATION_RATE` converts a step's frozen precipitation into
  meters of ice-depth gain, the same stylized-units-to-meters role `LAKE_FILL_RATE` plays
  for lakes. No cap on `glacier_depth` -- real ice sheets have no basin-capacity analogue.
- **Melt**: once a node warms back above the threshold, `GLACIER_MELT_RATE_M_PER_MYR`
  (scaled by how far above the threshold, capped at `GLACIER_MELT_MAX_FACTOR`) melts ice
  back down, capped so a step can't melt more than actually exists. The melted amount feeds
  directly into that step's water source for `route_downstream` -- real meltwater, feeding
  real river discharge, not a separate accounting bucket.
- **Flow**: a slope-scaled fraction of each node's ice (`GLACIER_FLOW_RATE_PER_MYR`, capped
  at `GLACIER_MAX_FLOW_FRACTION`) moves to its own `flow_target` each step -- the same graph
  water uses, via a direct scatter-add (`np.add.at`) rather than `route_downstream`'s
  elevation-ordered sweep, since glacier flow is a one-hop-per-step process, not a
  full-accumulation-to-terminus one. Ice reaching (or accumulating on) an ocean node is
  discarded rather than piling up -- real sea ice is a different, thinner, seasonal
  phenomenon this model doesn't represent; the same guard reads as calving where a glacier
  reaches a coast.
- **Erosion and flattening**: see [Erosion](#erosion) -- both driven by the *previous* step's
  `glacier_depth` (this step's fresh value isn't ready until this module runs, just before
  those terms are computed), the same one-step lag `channel_boost` already uses.

**Deliberately left out**: no rendering as a distinct color/
layer beyond the same `LAKE_COLOR_RGB`-style baking treatment lakes get (mantle-bloom has no
SNOW biome, so this uses its own `GLACIER_COLOR_RGB`, distinct from both `LAKE_COLOR_RGB` and
`elevation_colors`' own high-peak white/gray stops, applied the same nearest-neighbor-grid-resample way as lakes via
`plates.collect_all_glacier_depth`), no glacial eustatic sea-level coupling (glaciation is
purely local/per-node here), no seasonal accumulation/ablation cycle.

Confirmed live on a real run: glaciers form and grow preferentially at cold (polar)
latitudes, visually distinct from both lakes and high-elevation terrain on the rendered map;
population size and total ice depth fluctuate across steps as plates rotate glaciated nodes
in and out of genuinely cold regions, rather than growing monotonically or vanishing
outright.

<a id="river-inspector"></a>
## River Inspector

A third interactive map mode, same "raw JSON, client renders it" philosophy as
[Plate Inspector](#plate-inspector): `GET /world/rivers` returns every distinct river
network's flow-edge segments plus mouth metadata, plus the current coastline
(`coastline_segments` -- see [Coastline](#coastline); this view has no filled backdrop at
all, so without it there'd be no land/ocean/lake cue whatsoever), as plain JSON, and
`frontend/src/RiverInspector.tsx` renders and drives the interaction entirely client-side --
reusing the same `rotationDrag.ts` gesture (rotation is shared with `MapCanvas`/
`PlateInspector` via one lifted `App.tsx` state, so switching views preserves orientation),
plus click-to-select and Tab/Shift+Tab to cycle through rivers.

**There is no "a river" concept anywhere else in the codebase before this feature.**
`hydrology.py`'s `is_river` is (and remains) a flat per-node boolean mask -- the top decile of
land `flow_accum` -- with no grouping into distinct networks. `hydrology.group_rivers` adds
exactly that grouping, via union-find over `flow_target` edges restricted to nodes that are
`is_river` on *both* ends. This is exact,
not a heuristic: `flow_accum` is monotonically non-decreasing downhill, so once a node clears
the top-decile threshold, every node further downstream in its own chain does too -- meaning
two river nodes joined by a `flow_target` edge always belong to the same real drainage
network, and two nodes in different connected components never do. Grouping is recomputed
fresh from `world.hydrology_cache` on every `/world/rivers` (and `/world/river_at`) call
rather than persisted -- `hydrology_cache` itself is already at most one step stale by design,
and grouping it is cheap, so a `river_id` is only ever meaningful against the most recent
`/world/rivers` response for that reason, not a stable identity across steps (unlike
`plate_id`, which survives a plate's whole lifetime). `App.tsx` resets `selectedRiverId` on
every step for the same reason, not just on generate.

**A network's mouth is always its own max-`flow_accum` member** -- flow only accumulates
downhill, so nothing in the group can out-flow it, which also guarantees the mouth's own
`flow_target` (if any) points *outside* the group: a same-or-higher-`flow_accum` land node
downstream would itself have cleared the threshold and been unioned in already, so the only
things left outside are open ocean or a true dead-end sink. **`mouth_type`** checks
`lake_depth` first (a river can end at a still-spilling/draining lake, which reads more
usefully as `"lake"` than `"ocean"` even though it also has a `flow_target`), then whether
that `flow_target` lands on the ocean, else `"other"` (a dry interior sink).

**`num_tributaries` is an original definition** -- nothing existing to build on. Counts each
member's in-network in-degree (how
many *other* members flow directly into it); a member with in-degree 0 is a headwater, a
separate source stream with nothing upstream of it in this network. A single unbranched
channel has exactly one headwater (itself) and zero tributaries; each additional headwater is
one more distinct stream joining the network somewhere along its course, so
`num_tributaries = headwater_count - 1`.

**Speed** reuses `hydrology.compute_river_speed` (moved there from `erosion.py`, which only
ever consumed it, so `group_rivers` could share it
without a `hydrology.py -> erosion.py` reverse import) evaluated at the mouth, against
`_slope_to_flow_target`'s rise/run to the mouth's own `flow_target` -- the same slope
definition glacier flow already uses, distinct from (but conceptually the same kind of
quantity as) `erosion.py`'s own separate `SLOPE_NEIGHBOR_COUNT`-based slope.

**Rendering**: every river is always drawn as a set of flow-edge line segments (a flat edge
list, not an ordered polyline -- a network can branch), dim for every river and bright,
thicker for the selected one, mirroring Plate Inspector's dim/bright split. Each segment is
projected with the same "unwrap the second point's longitude relative to the first's own"
technique `render_image.py`'s `_draw_rivers` already uses server-side (`_project_offset`),
the two-point analogue of Plate Inspector's `projectLoop` (which unwraps a whole loop relative
to its own center) -- needed because a naive independent projection of both endpoints can bow
a short edge all the way across the map at the antimeridian seam. **The selected river's mouth
draws a ring** (`ctx.arc`, a new drawing primitive Plate Inspector doesn't need at all),
colored by `mouth_type` (`ocean`/`lake`/`other` each a distinct hue) so the endpoint
classification is visible at a glance, not just in the sidebar text.

**Click-to-select** is a server round trip (`GET /world/river_at`), the same nearest-node
`cKDTree` hit-test pattern as `plate_at` (`hydrology.river_at`, searching every network's own
member points at once) -- there's no shape to point-in-polygon test against, only a sparse
node cloud, so this is the natural fit rather than a departure from the pattern. No baked-PNG
elevation backdrop, matching Plate Inspector's own confirmed "the shapes are the whole visual"
choice, rather than adding one just for this view.

<a id="coastline"></a>
## Coastline (`coastline.py`)

The temperature/humidity/precipitation views and the River Inspector all needed a coastline
for the same reason: none of them have any other land/ocean cue at all. The elevation/plates
views show it implicitly through hypsometric coloring, and rivers/lakes draw directly onto
that same backdrop, but a climate value's color scale (or the River Inspector's otherwise
blank canvas) carries no land information on its own -- confirmed live: without this, a river
network on the River Inspector was just a handful of disconnected cyan squiggles floating on
black, with no way to tell where the coast actually was.

**Traced on climate.py's own grid, not the node cloud.** A coastline is inherently about
*land shape*, which needs a "here or not" decision at every point on the sphere -- climate's
`(H, W)` grid already covers the whole sphere densely and uniformly (see
[Climate](#climate)), exactly what's needed, and `world.climate_cache` is already there to
reuse (same one-step-stale caching philosophy as everywhere else that reuses it) rather than
building a second grid just for this.

**Land, ocean, and lake are three separate categories, not two.** `climate.py` has no lake
concept of its own (see its "deliberately not ported" list), so lake-ness is resampled onto
that same grid here, the same nearest-node `cKDTree` technique
`climate._sample_elevation_and_crust` already uses for elevation. The first version of this
folded "land" and "lake" into one combined category contrasted only against ocean -- which
turned out wrong the moment it was tested against an inland lake (one entirely surrounded by
land, the common case): both sides of that boundary read as the same "land-or-lake" value, so
no segment was ever drawn there at all. Fixed by keeping all three categories distinct
(`0=ocean, 1=lake, 2=land`) and tracing an edge wherever *any* pair of adjacent cells
disagrees -- ocean/land (the ordinary coastline), land/lake (an inland lake's own shoreline),
and ocean/lake (a lake that happens to sit right at the coast). A lake mask is also masked to
`~is_ocean` first: `is_lake` and `is_ocean` are two independent nearest-neighbor resamples (of
`hydrology_cache.points` and the plate node cloud respectively), so they can disagree right at
a coastline -- ocean always wins that conflict, rather than a cell reading as both at once.

**A rectilinear trace, not a smoothed contour.** For every pair of horizontally- or
vertically-adjacent grid cells whose category disagrees, the shared edge between them (half a
cell-step in from each of their centers) is one coastline segment -- deliberately matching the
grid's own existing visual language everywhere else it's drawn (render_image.py's
`_fill_rects`, itself rectangular per cell) rather than introducing a marching-squares-style
smoothing pass this codebase has no other use for. Segments are a flat edge list, not an
ordered polyline or closed loop, matching the same `(point_a, point_b)` convention
`hydrology`'s river edges already use -- land/lake shapes aren't simple loops (a lake can sit
inside a landmass that itself isn't one), so there's no ordering to preserve anyway.

**Two consumers, one computation.** `render_image._draw_coastline` projects each segment
(reusing `_project_offset`, same antimeridian-seam handling `_draw_rivers` already needs, for
the same reason: a real, short step between two adjacent grid points can still land on
opposite sides of the seam once the view can rotate arbitrarily) and strokes it into the
temperature/humidity/precipitation PNGs, plus a legend entry. `GET /world/rivers` sends the
same segments as `coastline_segments` (world-space xyz, same shape as `river.segments`) for
the River Inspector to project and draw itself -- see [River Inspector](#river-inspector).
Both call `coastline.compute_coastline_segments` directly rather than each re-deriving the
land/ocean/lake grid their own way, so the coastline a climate view bakes into its PNG and the
one the River Inspector draws client-side are always exactly the same shape.

**A single fixed line color would vanish against parts of the temperature gradient's own
white/black extremes** (see [Render image](#render-image)'s temperature stops), so this is
always drawn as a dark halo pass first, then a lighter line on top -- the same "halo" trick
real maps use for a boundary that has to stay legible over an arbitrary backdrop color, ported
identically to both the server-side PNG stroke and the River Inspector's canvas stroke.

<a id="known-simplifications"></a>
## Known simplifications

Deliberate scoping decisions for v1 (elevation only), each an acceptable line to draw rather
than an oversight:

- **The rendered outline (`Plate.outline_world`) is an envelope, not an exact polygon.**
  It's traced from each line's two current endpoints (ascending phi along the high-theta
  edge, descending back down the low-theta edge), so it's always in sync with the real
  territory -- unlike an earlier version that kept a separate polygon frozen at generation
  and rotated it rigidly, which drifted out of sync after enough stepping and could visibly
  overlap a neighboring plate's stale outline even though the underlying elevation data
  never did. The current approach can't represent a concave notch in a single line's
  extent (see the next bullet) -- a minor visual smoothing, not a data problem, since
  nothing other than rendering reads this outline.
- **Lines are assumed spatially contiguous.** A line's two ends (`theta[0]`, `theta[-1]`)
  are treated as its true territorial edges. A plate that develops a concave notch (or a
  split whose cut crosses one line's span twice) could in principle produce a
  non-contiguous line; this isn't specially detected. In practice a stray gap like that just
  gets treated as a small extra boundary by the next step's adjacency check and self-heals
  the same way any other boundary does.
- **Gap-filling absorption doesn't distinguish "genuinely my own natural growth direction"
  from "I happen to have the longest border nearby."** `MAX_ABSORB_NODES_PER_PLATE_PER_CALL`
  bounds the *rate*, but a plate that's already large still tends to keep winning the
  dominant-border check at its edges over many gap-fill passes, so it can still end up
  larger than a strict "each plate grows toward its own pole independently" rule would give
  -- a coarser approximation than the rest of the model, not a bug, and a possible future
  refinement (e.g. weighting by directional alignment with the plate's own motion rather
  than raw border presence).
- **No biomes/vegetation yet.** Climate (see [Climate](#climate)) feeds erosion, deposition,
  hydrology, and glaciation (see [Erosion](#erosion), [Hydrology](#hydrology), and
  [Glaciation](#glaciation): rain/river/weathering/glacier erosion, downstream deposition,
  rivers, lakes, and glaciers, with glacier flattening as a mantle-bloom-original addition).
  Vegetation's own effects (weathering boost,
  evapotranspiration) remain explicitly out of scope, since their input (a vegetation field)
  doesn't exist here at all -- that work would need to be designed from scratch if revisited
  on this sphere-native foundation.
- **No glacial eustatic sea-level coupling, no seasons.** Glaciation is purely local/per-node
  here -- see [Glaciation](#glaciation).
- **Single in-memory world, no persistence.** See
  [World state](architecture.md#world-state).
