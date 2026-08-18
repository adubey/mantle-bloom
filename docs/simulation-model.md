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
- [Boundary point reassignment](#reassignment)
- [Projections](#projections)
- [Render image](#render-image)
- [Rotating the view](#rotating-the-view)
- [Plate Inspector](#plate-inspector)
- [Climate](#climate)
- [Erosion](#erosion)
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

Two generation choices *are* user-facing -- the UI's "continental plates" and "initial land"
sliders, both 0 to 1 (percent in the UI), defaulting to `DEFAULT_CONTINENTAL_FRACTION = 0.70`
and `DEFAULT_LAND_FRACTION = 0.29`.

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
plate-local `(phi, theta)` lattice at `TARGET_LINE_SPACING_KM` resolution, and for every
candidate node, keep it only if this plate's seed is the *nearest* seed to it (`cKDTree`
against all seeds) -- the defining property of a spherical [Voronoi
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
  -200`) target -- new crust forming at the boundary. Unchanged: fades from 1 at zero
  distance to 0 at `FAR_THRESHOLD_RAD` (~200km), same as before.
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

**Every view draws a legend**, anchored at the image's bottom-left corner
(`_draw_legend`), baked into the PNG the same as everything else here rather than a
separate frontend overlay -- one shared panel layout (title, then an optional horizontal
color-gradient bar with tick labels, then any number of symbol rows) covers every case: a
plain color scale reusing that view's own `*_colors` function for the heatmap views
(elevation/temperature/humidity/precipitation, so the legend can never drift out of sync
with the actual fill colors), boundary/pole/rotation-arc symbols for "Plates", a gradient
plus a boundary symbol for "Plates (details)" (its node dots are colored by elevation, the
same scale as "Elevation"), and arrow/ocean/land/swell symbols for wind/oceanCurrents.

<a id="rotating-the-view"></a>
## Rotating the view

Every map view can be reoriented interactively: press and hold on the map, then drag (a full
3-DOF arcball gesture, roll/twist included -- not constrained to north-up panning) to spin a
virtual trackball; release to commit. This is a pure **view** transform, not a change to the
simulated world -- see below for why that distinction matters and how it's enforced.

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

- **Arcball**: a canvas-space point maps to a point on a virtual unit trackball centered on
  the canvas (`screenToArcballVector`, the standard Shoemake 1992 construction, points
  outside the ball's radius clamped to its silhouette rather than left undefined); the
  rotation between the drag's start vector and the current one is recomputed fresh from the
  fixed start point on every move (not accumulated incrementally per-event, avoiding drift
  over a long drag) and composed on top of whatever was already committed from prior drags.
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
- **The legend's lat/lon-of-center readout is the one deliberate exception to "every legend
  element is baked server-side."** It has to live-update during the drag, faster than any
  server round trip could allow, so it's rendered in `App.tsx` from `MapCanvas`'s live
  callback instead.

**Interaction**: press and hold (`LONG_PRESS_MS`, with a small movement tolerance before that
so an accidental brush of the map doesn't start a rotation) shows the graticule at the
currently committed orientation; dragging updates it live; releasing commits the drag's
rotation and triggers the real server render, the same `refresh` path a projection or map
view change already used, just with `rotation` now also part of what's requested.

`MapCanvas.tsx`'s drag gesture itself was later extracted into a shared hook,
`frontend/src/rotationDrag.ts`'s `useRotationDrag`, once the Plate Inspector view (below)
needed the exact same long-press-then-arcball interaction to rotate a completely different
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
ocean currents with swells marked, humidity, precipitation), ported from
[plate-sim](https://github.com/adubey/plate-sim)'s own climate model where its mechanisms are
richer than a first-principles description, with the caveat that mantle-bloom has no
vegetation, rivers, or lakes.

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
   cold is a consequence of solar heating at altitude, matching plate-sim's own
   `compute_temperature`).
3. **Ocean temperature baseline** = f(insolation) only, a narrower range plus a freezing
   floor (water's greater thermal inertia) -- pre-advection.
4. **Wind** -- latitude-banded meridional flow (trade winds/westerlies/polar easterlies,
   empirical lookup, the near-surface branch of the real three-cell circulation) plus
   Coriolis zonal deflection (`u = GAIN * sin(lat) * v`), plus an additive term from the real
   local gradient of the pre-advection land/ocean surface temperature (`np.roll`
   centered-difference) -- the empirical banding supplies planetary-scale structure a
   gradient-alone field doesn't produce (matching plate-sim's own documented finding), the
   gradient term is the genuinely temperature-responsive component. **Mountain
   deflection/Venturi/wake**, ported directly from plate-sim's `wind.py`: smooth elevation
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
   confirmed against plate-sim's actual `_land_swirl_current` to be nearest-cell-based, not
   landmass-grouped, so this is a faithful port), ramping from 0 at the coast to a peak then
   decaying with distance, direction matching the coastal deflection's hemisphere sense so
   the two agree. **Circumglobal boost**: a speedup on any row with zero land cells anywhere
   along it (a complete ocean ring) -- the stand-in for the Antarctic Circumpolar Current.
   **Wake**: the same backward-walk-and-damp structure as wind's, obstacle test is
   land-instead-of-elevation, plus a per-world-state noise texture (deterministic in `(seed,
   elapsed_years)`, so repeated renders of the same world state don't flicker) standing in
   for turbulent mixing.
6. **Ocean swells** -- convergence (negative divergence, `np.roll` centered differences) of
   the *final* current field, weighted-sampled down to `MAX_OCEAN_SWELLS` points (same
   weighted-without-replacement technique `hazards.py`-equivalent code uses elsewhere in
   plate-sim for earthquake/volcano placement).
7. **Final ocean surface temperature** -- semi-Lagrangian backward advection (single
   fixed-distance backward sample, nearest-cell) of the baseline along the final current
   field: "carried by ocean currents."
8. **Air temperature** -- the land baseline's own solar-heating formula, pulled toward the
   *nearest ocean cell's* final temperature by a distance-based (`exp` e-folding) falloff:
   "moderating effect of oceans," literally, right down to the "nearest ocean and its
   temperature" query being a `cKDTree` chord-distance search (true 3D distance, no
   pole/antimeridian special-casing needed) rather than plate-sim's own lat/lon-tangent-plane
   BFS. Land temperature itself is never moderated -- only air temperature is.
9. **Humidity** -- an evaporation ceiling over ocean from the local final ocean temperature,
   advected onto land by a wind-driven 2D sweep: a zonal pass and a meridional pass (each a
   sequential flow-direction walk, vectorized across the perpendicular axis, single-column
   evaporation/retention/orographic-dump step per iteration -- pure numpy, no `numba`),
   blended per-cell by each wind component's share of total wind magnitude. The zonal pass's
   sweep *direction* is the same fixed latitude-band lookup wind's meridional structure uses
   (`zonal_direction_for_lat`), not the literal local wind sign -- matches plate-sim's own
   `compute_humidity` exactly. No evapotranspiration term (needs vegetation, which doesn't
   exist here -- an absent input, not a simplified one).
10. **Precipitation** = f(humidity) + an orographic bonus (continuous saturating
    windward-slope moisture dump, from wind blowing up-elevation) -- no zonal
    latitude-climatology baseline (equator/mid-latitude wet bands), cut deliberately.
    Feeds erosion (see [Erosion](#erosion)) but still nothing else (no rivers, no
    vegetation).

**Scope, explicitly decided.** Kept out: river outflow feeding currents (no rivers exist),
deep currents, precipitation's zonal climatology baseline. Included, even though richer than
a one-line causal description: axial tilt, wind's mountain deflection/Venturi/wake, ocean
currents' coastal deflection/land swirl/circumglobal boost/wake -- all ported from directly
read plate-sim source, not simplified down. Dropped outright, not reduced (their inputs don't
exist in mantle-bloom): humidity's evapotranspiration term, river outflow, lake climate
influence.

**Rendering.** `render_image.py`'s `CLIMATE_VIEWS` (`temperature`, `wind`, `oceanCurrents`,
`humidity`, `precipitation`) route to `_render_climate_view`, a separate path from the
plate-tectonics views since the data source (a real `(H, W)` array, always covering the whole
sphere) is structurally different from the render grid's ragged lattice. Heatmap views
(temperature/humidity/precipitation) reuse the elevation view's color-stop-interpolation
technique with their own stop tables; wind/ocean-currents draw subsampled arrows (numpy-
vectorized projection/direction math, looped only for the unavoidable per-arrow PIL draw
calls), and ocean currents additionally marks each sampled swell point with a small circle.

<a id="erosion"></a>
## Erosion (`erosion.py`)

The other half of the weather<->geology coupling: [Climate](#climate) already has terrain
influencing weather (lapse-rate cooling, mountain wind deflection, orographic rain shadow --
all ported from plate-sim alongside the rest of its climate pipeline). This module is the
new direction, weather influencing terrain, ported from plate-sim's own erosion model but
cut down to the two sources that don't depend on infrastructure mantle-bloom doesn't have.

**Scope cut from plate-sim's five erosion sources.** River/channelized erosion and
deposition both need downhill flow-routing (plate-sim's `flow_target`/D8 accumulation),
which assumes persistent per-cell state on a fixed grid -- building that over a rotating,
irregular per-plate lattice is a separate, harder problem, not attempted here. Glacier
erosion is dropped (no glacier field exists). Weathering's vegetation boost is dropped (no
vegetation field exists, same reasoning as climate.py's own "deliberately not ported"
list). What's left -- rain/sheet erosion and weathering -- are the two sources that need
neither. Erosion here is one-way: material is removed, never redeposited elsewhere, so
there's no persistent sediment field either.

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
its nearest neighbors (0 if the node is already a local minimum -- matching plate-sim's own
"slope to lowest neighbor" definition), divided by the real great-circle distance to that
neighbor. This is a genuine dimensionless rise/run, unlike plate-sim's own slope --
documented there as a known simplification, "elevation drop per grid step, not
drop-over-real-distance" -- which is why `RAIN_EROSION_COEFFICIENT` isn't plate-sim's own
value (tuned against meters-of-drop-per-grid-cell, a very different scale); it was instead
picked by order-of-magnitude reasoning against `boundary.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR`
(800 m/Myr) and checked against real slope/precipitation distributions from an actual run.
`WEATHERING_COEFFICIENT` ports closer to plate-sim's own value, since wind speed uses the
identical scale in both codebases (`MERIDIONAL_BASE_SPEED = 6.0` in each).

**The two formulas**, both computed per-node then subtracted directly from `line.elevation`
(clamped to `MIN_ELEVATION_M`/`MAX_ELEVATION_M`, exactly like `boundary.py`'s own
uplift/trench deltas -- no resampling, no topology change, so this can't interact with line
regularization or point reassignment at all):

- **Rain/sheet erosion** = `RAIN_EROSION_COEFFICIENT * slope * (precipitation_mm / 1000) *
  dt_myr`.
- **Weathering** = `WEATHERING_COEFFICIENT * wind_speed * humidity_norm * dt_myr`, where
  `humidity_norm` is humidity clipped to `[0, 1]` against `HUMIDITY_REFERENCE`.
- Both zeroed over ocean nodes (`elevation <= 0`, the sea-level convention used everywhere
  else) -- rain/weathering erosion is a subaerial process, matching plate-sim's own
  `np.where(is_ocean, 0.0, ...)` masking of these same two sources.
- The *combined* result is capped at the node's own drop-to-lowest-neighbor (in meters, not
  the normalized slope) -- the same cap plate-sim applies, so a single step can't erode a
  node *past* the valley it drains into and carve a new, lower pit.

**Cadence: every step, no lag.** Unlike plate-sim (whose erosion reads the *previous* step's
climate, because its climate is only computed once per step, late in a long pipeline, after
erosion has already run), this module calls `climate.compute_climate(world)` fresh,
immediately before using it -- climate.py is already fully stateless and cheap enough to
call on demand, so there's no staleness to reason about. Runs in `world.step_world` right
after boundary evolution and topology changes, every step (not gated by the periodic
regularize/gap-fill/reassign cadence) -- a deliberate cost/fidelity tradeoff: climate
computation adds real per-step latency (roughly doubling a step's cost), accepted in
exchange for erosion responding smoothly every step rather than arriving in periodic bursts.

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
- **No hydrology (rivers/lakes) or biomes/vegetation yet.** Climate is modeled (see
  [Climate](#climate)) and now feeds erosion (see [Erosion](#erosion): rain/sheet erosion
  and weathering, both ported from plate-sim), but river/lake/vegetation systems -- and the
  channelized erosion, glacier erosion, and deposition that would depend on them -- remain
  explicitly out of scope for this v1. See plate-sim's own model for the shape that work
  would take once revisited on this sphere-native foundation.
- **Single in-memory world, no persistence.** See
  [World state](architecture.md#world-state).
