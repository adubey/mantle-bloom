# Simulation Model

## Table of contents

- [Why not a grid](#why-not-a-grid)
- [Plate-local frames](#plate-local-frames)
- [Initial plate generation](#initial-plate-generation)
- [Mantle flow](#mantle-flow)
- [Boundary evolution](#boundary-evolution)
- [Garbage collection](#garbage-collection)
- [Merge and split](#merge-and-split)
- [Whole-sphere coverage (gap-filling)](#gap-filling)
- [Projections](#projections)
- [Render grid](#render-grid)
- [Known simplifications](#known-simplifications)

<a id="why-not-a-grid"></a>
## Why not a grid

mantle-bloom is a from-scratch successor to
[plate-sim](https://github.com/adubey/plate-sim), which modeled the planet as an
equirectangular lat/lon grid. plate-sim's own docs catalog the problems that caused: pole
cells need an artificial full-clique patch to be mutually adjacent, falloff radii need a
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
  anything. This is the direct fix for plate-sim's semi-Lagrangian mass-conservation
  problem.
- **Equidistant and parallel, by construction.** Equal `delta-phi` between lines is already
  physically equidistant (meridional spacing on a sphere doesn't depend on latitude), and
  each line's `delta-theta` node spacing is chosen from that line's angular radius
  (`cos(phi)`) to hit `TARGET_LINE_SPACING_KM` (`plates.py`, default 125 km) -- the direct
  fix for plate-sim's documented "latitude distortion" issue.
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

Continent count *is* user-facing -- the UI's continents slider passes `num_continents`
(`plates.MIN_CONTINENTS = 1` to `plates.MAX_CONTINENTS = 8`). When given, exactly that many
plates (`rng.choice`, without replacement) are made continental instead of the usual
independent `CONTINENTAL_FRACTION = 0.4` coin flip per plate, and `num_plates` is bumped up
if needed so there's still room for at least `MIN_OCEANIC_PLATES` of real ocean floor
regardless of how many continents were requested.

Each plate's elevation lines are populated by `plates.iter_local_lattice`: sweep a full
plate-local `(phi, theta)` lattice at `TARGET_LINE_SPACING_KM` resolution, and for every
candidate node, keep it only if this plate's seed is the *nearest* seed to it (`cKDTree`
against all seeds) -- the defining property of a spherical [Voronoi
diagram](https://en.wikipedia.org/wiki/Voronoi_diagram), computed directly rather than via
an explicit polygon-construction step. Every node ends up owned by exactly one plate, so the
initial tiling has no gaps and no overlaps *by construction* -- there's nothing to
separately verify. Kept nodes get a base elevation by crust type
(`BASE_CONTINENTAL_M = 200`, `BASE_OCEANIC_M = -3800`) plus a smooth noise texture
(`noise.py`, a small sum of sinusoids with random frequency/phase -- not true gradient
noise, just enough texture to not look perfectly flat).

The same lattice-sweep helper (`plates.build_lines_from_lattice`) is reused by plate merging
(see [Merge and split](#merge-and-split)) -- the only other place a full-footprint sweep is
needed.

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

For each node within `FAR_THRESHOLD_RAD` (1.6x target spacing) of some other plate's
nearest node, the two plates' relative velocity at that point (from their `omega`s) is
decomposed against the direction toward the neighbor into a **closing rate**: positive means
this plate's material is moving toward the neighbor's (convergent), negative means moving
apart (divergent); `TRANSFORM_RATE_THRESHOLD` (~1 cm/yr equivalent) separates both from
transform.

- **Convergent + continental** -> elevation rises (`CONVERGENT_MOUNTAIN_RATE_M_PER_MYR`,
  mountain building). **Convergent + oceanic** -> elevation falls
  (`CONVERGENT_TRENCH_RATE_M_PER_MYR`, trench/subduction). Both scaled by an `intensity`
  factor that fades from 1 at zero distance to 0 at `FAR_THRESHOLD_RAD`, so the effect is
  concentrated right at the boundary without a separate distance-falloff pass.
- **Divergent** -> elevation relaxes exponentially toward a ridge (`oceanic`,
  `DIVERGENT_RIDGE_TARGET_M = -1500`) or rift (`continental`,
  `DIVERGENT_RIFT_TARGET_M = -200`) target -- new crust forming at the boundary.
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
  no plate-sim-style resampling step that can lose or duplicate it.

<a id="garbage-collection"></a>
## Garbage collection (`line_regrid.py`)

Per-step boundary evolution only ever touches a line's two ends, so interior spacing stays
regular on its own during ordinary convergent/divergent motion -- but a *transform* boundary
shears nodes along a line without inserting or deleting anything, which can leave spacing
uneven. Every `GC_INTERVAL_STEPS` calls to `step_world` (default 5), any line whose gaps
have drifted past `IRREGULARITY_TOLERANCE` (1.5x target spacing, either direction) gets a
fresh evenly-spaced node set across its *existing* extent -- the two endpoints are preserved
exactly, since GC never changes where a line's physical edge is, only how regularly it's
sampled -- with elevation re-interpolated onto the new nodes (`np.interp`, 1D since it's
along a single already-ordered curve, not 2D scattered-data interpolation).

<a id="merge-and-split"></a>
## Merge and split (`merge_split.py`)

Unlike rotation and boundary evolution, these are rare, discrete, topology-changing events,
so a one-time resample is an acceptable cost here -- the exact, no-resampling guarantee only
matters for routine per-step motion.

- **Consumption.** A plate whose every elevation node has been deleted (fully subducted) is
  simply dropped from `world.plates` -- falls directly out of the boundary-evolution rule
  above, no special algorithm needed.
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
  anything that happened this step -- a new collision starting, a merge completing (with
  how many million years it took), consumption, or a split -- which `world.step_world`
  timestamps with the *post-step* `elapsed_years` and appends to `World.events` (capped at
  `world.MAX_EVENT_LOG_LENGTH` entries). `generate_world` logs an initial "world generated"
  event the same way. The API returns the full current log on every `/world/generate` and
  `/world/step` call (see api-reference.md) for the frontend's collapsible console.
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

Every `line_regrid.GC_INTERVAL_STEPS` calls (the same cadence as garbage collection),
`gaps.fill_gaps` sweeps a global lattice (`plates.iter_local_lattice` in the identity frame,
reused as a plain lat/lon sweep purely for this one-off detection query), finds every
candidate point farther than `COVERAGE_RADIUS_RAD` from any plate's nearest node, and
clusters the results (`scipy.sparse.csgraph.connected_components` over a k-d-tree radius
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
  arbitrarily assigned to one side.

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

<a id="render-grid"></a>
## Render grid (`main.py`)

The elevation-line/gap-filling data is genuinely Lagrangian: nodes are spaced at
`TARGET_LINE_SPACING_KM` *within each plate's own local frame*, not on any shared, screen-
aligned grid. Drawing that raw point cloud directly -- one dot per node -- leaves visible
gaps once projected: projected spacing isn't uniform (Behrmann, for instance, stretches
longitude spacing by roughly 50x at high latitudes relative to the equator, since
`x = lon * cos(30deg)` doesn't compensate for the shrinking circumference the way the
latitude term does), so a dot sized to look right at the equator leaves gaps near the
poles, and no single fixed dot size closes the gap everywhere without grossly overlapping
elsewhere.

`_render_grid` fixes this the way the user asked: sweep a uniform lat/lon grid over the
whole sphere (`plates.iter_local_lattice(identity_frame, spacing_rad=GRID_SPACING_RAD)` --
the same identity-frame trick `gaps.py` uses for its coverage sweep, at
`GRID_SPACING_RAD = TARGET_LINE_SPACING_RAD`), assign every cell its nearest elevation
node via one `cKDTree` query against every plate's current nodes (no distance cutoff --
every cell gets *some* value, so there's no gap by construction regardless of how far the
nearest real data happens to be), and project the whole grid the same way everything else
is projected. This is a one-time resample purely for rendering -- it never touches
`world.plates`, so it has no bearing on the mass-conservation properties the rest of the
model is built around (see [Why not a grid](#why-not-a-grid)); a fresh grid is computed
from scratch on every `/world/render` call.

**Sizing each cell correctly is what actually closes the gaps.** A uniform *sphere* grid
still isn't a uniform *projected* grid -- so each cell is drawn at a size measured from the
projection's own local behavior, not a fixed pixel size: `_row_cell_half_extent` projects
two extra nearby samples per row (one step further in theta, one row further in phi) and
measures the resulting on-screen offset directly, giving that row's `cell_half_width` and
`cell_half_height` in the same projected units as the point coordinates. theta-only and
phi-only offsets are used deliberately so each measurement isolates one derivative exactly
(in both projections, moving in longitude at fixed latitude never perturbs the projected y
coordinate, and vice versa) rather than approximating a mixed partial derivative. The
frontend (`MapCanvas.tsx`) draws each cell as a rectangle of that measured size (with a
small `CELL_OVERLAP_FACTOR` margin so adjacent cells overlap a hair rather than risk a
hairline gap from floating-point rounding) -- confirmed visually to close every gap in both
projections, including the sparsest cells right at the poles.

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
- **No climate, hydrology, erosion, or biomes yet.** Explicitly out of scope for this v1 --
  see plate-sim's own model for the shape that work would take once revisited on this
  sphere-native foundation.
- **Single in-memory world, no persistence.** See
  [World state](architecture.md#world-state).
