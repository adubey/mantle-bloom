# Simulation Model

## Table of contents

- [Why not a grid](#why-not-a-grid)
- [Plate-local frames](#plate-local-frames)
- [Initial plate generation](#initial-plate-generation)
- [Mantle flow](#mantle-flow)
- [Plate motion: shift and deform](#boundary-evolution)
- [Line regularization](#line-regularization)
- [Merge and split](#merge-and-split)
- [Whole-sphere coverage (subsumed into deform)](#gap-filling)
- [Volcanism](#volcanism)
- [Boundary point reassignment (subsumed into deform)](#reassignment)
- [Projections](#projections)
- [Render image](#render-image)
- [Rotating the view](#rotating-the-view)
- [Plate Inspector](#plate-inspector)
- [Climate](#climate)
  - [Biomes](#biomes)
- [Ocean/Atmospheric Fluid Dynamics](#ocean-atmospheric-fluid-dynamics)
- [Erosion](#erosion)
- [Bathymetry](#bathymetry)
- [Resources and soil](#resources-and-soil)
- [Hydrology (rivers and lakes)](#hydrology)
- [Glaciation](#glaciation)
- [River Inspector](#river-inspector)
- [Lake Inspector](#lake-inspector)
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
having to pick a number. That many *primary* seed points are scattered uniformly on the unit
sphere (normalized Gaussian samples), and each gets a plate-local frame built from its own
seed (`geometry.plate_frame_from_seed`).

A plate is not a single Voronoi cell, though: alongside the primaries,
`lithosphere_plate.build_plate_tiling` scatters `EXTRA_SITES_PER_PLATE = 2` extra sites per
plate and hands each one to a plate by Prim-style region growing (the still-unassigned site
angularly closest to any already-assigned site joins that site's plate). A plate's territory
is then the *union* of its own sites' Voronoi cells -- still a nearest-site lookup, so still
gap/overlap-free by construction, but the merged cells give lumpier, less convex outlines
than one-cell-per-plate did. `extra_sites_per_plate = 0` recovers the original tiling. Kept
modest deliberately: the more cells a plate fuses, the more concave its outline gets, and
`PlateWithLines`' per-row outline is only an envelope for a genuinely non-convex shape.

Three generation choices *are* user-facing -- the UI's "continental plates" and "initial
land" sliders, both 0 to 1 (percent in the UI), defaulting to
`DEFAULT_CONTINENTAL_FRACTION = 0.70` and `DEFAULT_LAND_FRACTION = 0.29`, plus a "point
density" choice (`NODE_DENSITY_CHOICES = (0.5, 1.0, 2.0, 4.0)`, see below).

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
see "Elevation point density" below), and for every candidate node, keep it only if the
*nearest* site to it (`cKDTree` against all sites, primary and extra) is one this plate owns
-- the defining property of a spherical [Voronoi
diagram](https://en.wikipedia.org/wiki/Voronoi_diagram), computed directly rather than via
an explicit polygon-construction step. Every node ends up owned by exactly one plate, so the
initial tiling has no gaps and no overlaps *by construction* -- there's nothing to
separately verify. Kept nodes get a base elevation by crust type
(`BASE_CONTINENTAL_M = 200`, `BASE_OCEANIC_M = -3800`, continental overridden by
`land_fraction`'s threshold when given, see above) plus a smooth noise texture (`noise.py`,
a small sum of sinusoids with random frequency/phase -- not true gradient noise, just enough
texture to not look perfectly flat), swung by `CONTINENTAL_NOISE_AMPLITUDE_M`/
`OCEANIC_NOISE_AMPLITUDE_M` (2000/900, widened from an original, more timid 1200/500 so a
freshly generated world already shows real relief -- rolling hills and real ocean-basin depth
variation -- rather than a near-flat plain/seafloor waiting for tectonics to draw the first
contours). Confirmed directly across a couple of seeds: land elevation's own standard deviation
rose by roughly 50% (e.g. 306m -> 479m at one seed) and the deepest sampled ocean point got
several hundred meters deeper, without disturbing the land/sea split near sea level itself
(`BASE_CONTINENTAL_M`, or the `land_fraction`-derived threshold when one is given, is
untouched -- only how far the noise texture swings around whichever baseline is already used).

The same lattice-sweep helper (`plates.build_lines_from_lattice`) is reused by plate merging
(see [Merge and split](#merge-and-split)) -- the only other place a full-footprint sweep is
needed.

**Elevation point density.** The UI's "point density" choice (`node_density`, 4x default,
0.5x/1x/2x also available) scales `TARGET_LINE_SPACING_RAD` down via `plates.line_spacing_rad`
(halved at 4x -- node count for a fixed area scales with the *square* of resolution, so 4x the
nodes needs half the spacing, not a quarter; conversely 0.5x, the coarsest option, doubles the
spacing). Stored on `World.node_density`, set once at generation and
read for that world's entire life, not just at the moment it's generated: every later module
that builds new elevation-line nodes or derives a distance/count threshold from
`TARGET_LINE_SPACING_RAD` -- `elevation_lines.py`'s line regularization, `plates.py`'s
`PlateWithLines.deform` (per-turn growth/shrink/claim thresholds -- see [Plate motion: shift
and deform](#boundary-evolution)), `merge_split.py`'s plate-merge contact distance,
split-size floor, and defragmentation connect-radius / fragment-size floor -- calls
`line_spacing_rad(world.node_density)` (or scales its own reference constant by the same
ratio) instead of reading the bare module constant. This matters because
it's not just a generation-time cosmetic choice: `elevation_lines.py`'s regularize pass in
particular runs at the end of every single `deform()` call now (not periodically -- see [Line
regularization](#line-regularization)) and, before this threading existed, always resampled a
line back down to the *reference* spacing regardless of what density the world was actually
generated at -- confirmed directly as a real bug during development, a 4x-density world's
own node count reverting to the 1x baseline within the first handful of steps. Every
distance-based threshold derived from `TARGET_LINE_SPACING_RAD` (e.g. `plates.py`'s
`EXTEND_THRESHOLD_RAD`) scales linearly with the new spacing; every absolute node-*count*
constant tied to a fixed physical area (`merge_split.SPLIT_MIN_NODES`, etc.) scales with
`node_density` directly, not its square root, since it's already an area -- see each
constant's own comment for the exact reasoning, which predates this option (the same
rescaling used to happen as a one-off hardcoded code change whenever `TARGET_LINE_SPACING_KM`
itself changed; this option just makes it a per-world runtime choice instead). Genuine fixed
physical distances unrelated to sampling resolution (e.g. `plates.COLLISION_RANGE_KM`, a real
~400km-wide collision belt) are deliberately *not* scaled -- only thresholds explicitly
defined as multiples of `TARGET_LINE_SPACING_RAD` are. 4x density comes with a real,
continuous performance cost, not just a one-time generation cost -- confirmed directly,
several times slower per-step time (not just 4x, since the polygon-containment classification
`deform()` runs against every near-boundary node is closer to `O(n log n)` than linear) --
which the UI surfaces as a short note when 4x is selected.

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
## Plate motion: shift and deform (`plates.py`'s `Plate.shift`/`Plate.deform`)

Replaces an earlier design (`boundary.py`'s `step_boundaries`) that classified a boundary
node as convergent/divergent/transform from the *velocity* decomposition described further
down -- two plates' relative velocity at a point, decomposed against the direction toward
the neighbor into a signed **closing rate**. That mechanism still exists (`boundary.
closing_rate`) and is still exactly how `merge_split.py` decides whether two continental
plates are actively colliding (see [Merge and split](#merge-and-split)) -- it just no longer
drives ordinary per-step boundary evolution. Classification there is now purely
*geometric*: did a plate's rotated territory end up overlapping a neighbor's, or did it open
up space nobody else claims?

Each plate runs two operations every step, in this order:

1. **`shift(world, years)`** -- refit the plate's Euler pole from the mantle-flow field
   (damped toward the new target, same fit-then-clamp as generation-time, see [Mantle
   flow](#mantle-flow)), then rotate rigidly by `years` at that rate: exact for every node
   the plate carries, no resampling (see [Plate-local frames](#plate-local-frames)). Returns
   `D`, the greatest angular distance any of the plate's own nodes actually moved this
   step -- a real, physically grounded bound on how much crust could plausibly have been
   created, destroyed, or stretched this turn. Runs for every plate; order doesn't matter
   here, since rotating one plate only ever touches its own `frame`.
2. **`deform(world, other_plates, years, D)`** -- reconcile the plate's actual post-`shift`
   footprint against the footprint it's entitled to occupy: the sphere minus every *other*
   currently-live plate's own bounding polygon (`Plate.get_bounding_polygon()`, a live,
   cached outline derived directly from the plate's own current line endpoints -- see [Known
   simplifications](#known-simplifications) for why this is an envelope, not an exact
   polygon). Runs once per plate, in a **freshly randomized order every turn** -- the reason
   is explained below.

A node is **contested** if it's now geometrically contained in some other live plate's
polygon (`geometry.points_in_spherical_polygon`, checked only for nodes a cheap k-d-tree
distance prefilter can't already rule out -- deep-interior nodes are never near enough to
matter). Contested nodes are what used to be "convergent"; nodes that are uncontested and
far from every neighbor are what used to be "divergent" (open, unclaimed territory); nodes
that are uncontested but still close to a neighbor are "transform." This reframing needs no
real polygon union/intersection/subtraction machinery at all: "is this point in the union of
every other plate's territory" is just "is it contained in *any* one of them," a plain
boolean OR over the same containment test used everywhere else.

**Why randomize the processing order.** Two plates can both border the same unclaimed
patch of sphere (a polar cap nobody's reached yet, or ground a subducted neighbor just
vacated) -- if both tried to claim it in the same turn from a shared starting snapshot,
they'd overlap. Instead, each plate's "what am I entitled to" check runs against whatever
state its neighbors are *currently* in: a neighbor already deformed earlier this same turn
reflects this turn's change, one not yet reached doesn't. Whichever plate's turn comes first
claims the space; by the time a later plate checks, that space is already gone from its own
entitled footprint. Randomizing the order every turn (seeded by `(world.seed,
round(world.elapsed_years))`, so a replayed session still reaches the same order given the
same history) is what keeps this fair on average across many turns, rather than always
favoring whichever plate happens to be processed first.

**Convergent effects aren't a single shape** -- what happens depends on both plates'
crust type, and how far the effect reaches (and its shape with distance) differs by case.
The rates/reaches below are unchanged from the model `step_boundaries` used to run; only the
trigger (contested, not a positive closing rate) changed:

- **Continent-continent collision** (both plates continental) -> elevation rises
  (`CONVERGENT_MOUNTAIN_RATE_M_PER_MYR`), scaled by an intensity that fades from 1 at zero
  distance to 0 at `COLLISION_RANGE_RAD` (400km) -- a broad crumple zone, matching how wide a
  real collision belt is (e.g. the Himalaya/Tibetan Plateau). The same collision also adds a
  second, much gentler rise (`FAR_FIELD_MOUNTAIN_RATE_M_PER_MYR`, well under a tenth of the
  near-field rate) far inland: zero out to `FAR_FIELD_COLLISION_INNER_RAD` (1000km, leaving
  the 400-1000km gap where the near-field crumple zone has already faded to nothing
  untouched), ramping to full intensity there and back to zero by
  `FAR_FIELD_COLLISION_OUTER_RAD` (3000km). Real collisions transmit stress this far into the
  continental interior -- the Himalayan-Tibetan and Arabian-Eurasian (Zagros) collisions both
  have deformation reaching comparable distances (Tien Shan/Baikal, Anatolia), the
  Variscan-Appalachian and Uralian orogenies both left belts wider than their core sutures,
  and the Laramide orogeny's basement-cored uplifts sat ~1000-1500km inland of the margin.
  Unlike the divergent cases below, this doesn't relax toward a target -- it adds every step
  it applies, so it can accumulate into a substantial rise over a long-lived collision.

  **Reverse faults: mountain ranges aren't uniformly smooth.** Real shortening in a collision
  belt isn't spread evenly across the whole zone -- fold-thrust belts partition it into
  discrete thrust sheets (fast-rising ridges) separated by footwall synclines/intermontane
  basins that keep rising far more slowly, a real, documented process (the north-south rift
  valleys cutting straight across the Tibetan Plateau's own overall convergent uplift;
  Basin-and-Range-style extension nested inside the Anatolian collision zone). Both the
  near-field collision uplift and the subduction volcanic-arc uplift above are multiplied by a
  `fault_factor` -- 1.0 almost everywhere, dropping to `REVERSE_FAULT_VALLEY_UPLIFT_FACTOR`
  (0.15) on whichever nodes a smooth, deterministic noise field marks as a downthrown block
  (below `REVERSE_FAULT_VALLEY_THRESHOLD`). Those nodes still rise (this is differential uplift
  within an active belt, not literal subsidence), just far slower than their neighbours, so a
  real valley opens up between ranges as the gap widens step after step. The noise is sampled
  in the plate's own *local* frame (`geometry.local_xyz(line.phi, line.theta)`, not world xyz),
  seeded from `(world.seed, plate_id)` only -- so a given downthrown block stays attached to
  the same crust as the plate rotates, a fixed geological feature rather than something that
  reshuffles every step, the same "attached to the crust, not the world" property every other
  persistent field in this codebase already has (see [Why not a grid](#why-not-a-grid)).
  Deliberately not applied to the far-field term, which represents stress transmitted broadly
  into the continental interior, not the belt's own discrete thrust-sheet structure. Confirmed
  directly at a real seed run 20 steps (60 Myr): with this and the seismic-erosion addition
  below, the fraction of land nodes pegged at `MAX_ELEVATION_M` dropped from roughly 9% to
  under 2% versus the same run without either -- mountain ranges keep growing, but no longer
  collapse into a flat plateau at the elevation ceiling.
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
  edge of that line's territory): if an end is uncontested and its gap to the nearest
  neighbor node has opened past `EXTEND_THRESHOLD_RAD`, insert new nodes at target spacing --
  as many as it actually takes to close the gap, capped by both `D` (this step's own physical
  bound -- see `shift()` above) and `MAX_EXTEND_NODES_PER_STEP` as a hard safety ceiling, not
  the normal limit. Each new node gets the ridge/rift target elevation (brand new material,
  not interpolated from anything), *unless* the growth event rolls "overstretched" (see
  below), in which case it comes back as a fresh volcano instead.

  **Torque engine (`lithosphere_plate.py`): new areal crust is always oceanic.** Under the
  Hc/Hm isostasy engine, a growing end (`_grow_or_shrink_line_for_deform`) and a claimed new
  phi row (`_claim_adjacent_territory`) seed the new nodes' lithospheric column via
  `lithosphere_plate.growth_seed_thickness()` -- the *oceanic* reference column regardless of
  the plate's own `crust_type`, because any gap that opens on the sphere is floored by
  sea-floor spreading, not by the neighbouring plate's crust. Seeding a continental plate's
  own reference column there was a real land-area runaway: continental plates continuously
  grow into space subducting oceanic plates vacate, and continental crust never subducts
  back, so every such step permanently converted ocean floor into ~+200 m dry land (measured
  land fraction climbed 0.27 -> 0.48, mean planet elevation rose ~1.7 km over 180 Myr on one
  seed). New oceanic crust on a continental plate lands ~-3.5 km -- a drowned passive margin
  / accreted terrane. Genuine continental rifting is untouched: that thins *existing* crust
  (`rheology.apply_divergent_deformation`), it doesn't grow new nodes here.

  If an end is contested,
  remove however many *consecutive* contested nodes sit there, capped the same way by `D` and
  the safety ceiling -- but never the plate's last remaining node in a line. Growth and
  *ordinary* shrink are end-only: each `ElevationLine` is a single contiguous arc, and
  inserting/deleting anywhere but an end would break that. The one interior case handled is
  **interior subduction** on an oceanic self-plate: a run of at least
  `_INTERIOR_SUBDUCTION_MIN_RUN` contested nodes stranded in a row's *middle*, live nodes on
  both sides (a neighbor -- typically a continental plate, whose own contested nodes never
  subduct -- has rotated bodily over a mid-row patch faster than the end-shrink could
  retreat). That run is carved out and the row's survivors returned as *two* separate
  contiguous `ElevationLine`s at the same `phi`, one per arc; `PlateWithLines.outline_world`
  and `contains_batch` both handle several lines at one `phi`, tracing the gap between the
  arcs as a genuine hole (keyholed out of the plate's polygon, see [Known
  simplifications](#known-simplifications)) rather than claiming it. Total interior deletion
  per call is capped by the safety ceiling, not `D` -- those nodes were overridden over many
  past steps, so clearing them is catch-up cleanup, not this step's own subduction. Before
  this, an interior-only contested patch was left untouched every turn on the reasoning that
  "the neighbor's own growth reaches this row's nearer end before long" -- which never
  happens when the neighbor plants a lobe mid-row and then stops advancing, leaving a frozen
  continental-over-oceanic overlap that didn't heal (seen in a real long-run world:
  `seed 888151728` at 6.9 Myr, plates 9 and 1). This is also where mass conservation lives:
  material is only ever created at open ends and destroyed where genuinely overridden, as
  literal point insertion/deletion -- there's no grid-resampling step that can lose or
  duplicate it.

  **A row never winds past a full revolution.** A line is a circle of plate-local latitude,
  so its theta extent physically can't exceed `2*pi`. Nothing here treats theta as periodic,
  and near a plate's own local pole the "gap to the nearest neighbor" reads as wide open
  forever (the polar cap belongs to nobody), so without a guard a plate that has grown to
  encircle its pole just keeps winding the same ring -- the same ground covered many times
  over (`theta` spans of tens of revolutions were seen on long runs, showing up as
  concentric circles / moire "holes" in the Plate Inspector and as an unbounded contribution
  to plate overlap and node count). End-growth is capped so a row's span never exceeds one
  revolution; once it closes, that end stops. `elevation_lines.regularize_line` also unwinds
  an already-over-wound row (keeping the outermost single revolution) so a world saved before
  this guard heals on its next few steps.
- **Overstretched growth becomes volcanic.** Ordinary ridge/rift fill at a growing end
  instead spawns a fresh volcano (`is_volcano=True`, a random `volcano_active_years_remaining`
  draw, and one *guaranteed* immediate eruption -- see [Volcanism](#volcanism)) with a small
  fixed probability per growth event (`STRETCH_VOLCANO_PROBABILITY`, 0.02), representing "the
  plate has been stretched too thin to keep filling with plain crust." Deliberately
  probabilistic rather than a threshold on some per-call quantity: two threshold designs were
  tried and rejected during development. Checking whether *this line's own existing gap*
  already exceeds target spacing turned out to be dead code -- line regularization (below)
  runs at the end of every `deform()` call and resamples every line back to (within
  tolerance of) exact target spacing, so by the time the *next* call's growth check runs, any
  such gap has already been smoothed away by the *previous* call's own regularize pass.
  Checking whether *this call* needs to insert several nodes at once was confirmed
  empirically unreachable at realistic step sizes and plate rates: sampled 1392 real growth
  events across a running simulation, and 100% of them inserted exactly one node, since
  ordinary per-step divergence rarely outruns a single spacing unit's worth of growth in one
  call regardless of how the threshold was tuned. A small per-event probability sidesteps
  needing any persistent "how long has this been thinning" state (which would have to survive
  regularization, split, and merge) while still producing "occasionally, not constantly"
  volcanic crust at active rifts over a real run.
- **Claiming adjacent territory** (subsumes the old whole-sphere gap-filling pass described
  below): after growing/shrinking each existing line's ends, a plate checks for a whole new
  phi row just beyond its current phi extremes -- the one case ordinary end-growth
  structurally can't reach, since growth only ever extends an *existing* line's own theta
  range, never adds a new line. The new row is only added if it stays at least
  `POLE_CAP_MARGIN_MULT` target spacings clear of the plate's own local pole (`+-pi/2`):
  right at the pole a row's theta step (`spacing / cos(phi)`) blows up and the row degenerates
  into a handful of sub-spacing-circumference rings, so growth *toward* the pole stops short
  and a small permanent polar cap is left unclaimed (harmless -- the render grid's
  nearest-node fill covers it, and `deform()`'s contested test doesn't care). Generation's own
  lattice sweep still fills to the pole, so a plate that legitimately owns its pole *at
  generation* keeps those rings. If that row is open (uncontested) territory, it's added
  directly as a new `ElevationLine` -- not via a full lattice resample (`Plate.grow_into`,
  used elsewhere for rare, merge-scale events): calling that every turn for every plate was
  tried and rejected during development, confirmed to balloon a plate's own node count by
  several times in a single call, since a full resample's own coverage radius around a
  handful of newly-claimed points reconstructs far more lattice area than just those points.
  Reclaiming ground a subducted neighbor vacated *within* an existing row's own theta range
  needs no separate mechanism at all -- the very next time that row's end-growth check runs,
  the vacated neighbor is simply gone from the distance/contested check, and ordinary
  end-growth already extends into it.

**Known limitation: overlap isn't exactly zero, but stays bounded.** `Plate.
get_bounding_polygon()` is an envelope (see [Known simplifications](#known-simplifications)),
and the randomized-order design above means a plate's own "what am I entitled to" check can
be up to one turn stale against a neighbor not yet processed this same turn. Both mean a
node can transiently read as inside a neighbor's polygon without the two plates' actual node
clouds genuinely interpenetrating. Confirmed directly across many stepped turns: sampled
overlap stays low single-digit percent and doesn't grow -- bounded, self-correcting behavior,
not a runaway. (It used to sit in the low teens: the interior-subduction carve-out above
removed the one case -- a neighbor's lobe frozen mid-row on an oceanic plate -- that stayed
put indefinitely rather than self-correcting. Re-running `seed 888151728` from 6.9 Myr, the
plates 9/1 envelope overlap drops from ~15% to a bounded ~2-3% on the first step and stays
there.) A stricter, exactly-zero invariant would need either a self-intersection-safe polygon
construction or a supplementary node-cloud distance guard; not pursued for v1.

**Deep interpenetration is classified too (torque engine, 2026).** The polygon-containment
test used to run only on nodes within `reach_rad` (~3 spacings) of a neighbour node, so a
plate that had slid *deep* over another had its deep-interior overlapping nodes classified
as neither contested nor divergent -- no thickening, no subduction, the overlap just sat
(observed: 15% of one plate on top of another, static). `torque.classify_boundary_nodes` now
also polygon-tests any node inside a neighbour's bounding sphere (cheap triangle-inequality
prefilter). A deep continental overlap then classifies contested -> `rheology` thickens
Hc/Hm -> **mountain uplift**; a deep oceanic overlap classifies contested -> subduction
deletion -> the overlap heals. Every node's onset year is also stamped onto
`ElevationLine.overlap_onset_years` each step (`merge_split.update_overlap_tracking`) and
surfaced as `since_years` in `GET /world/plates` and the `overlapAge` debug render view --
see [debugging.md](debugging.md#overlapage-render-view-plate-overlap-onset).

<a id="line-regularization"></a>
## Line regularization (`elevation_lines.py`)

Per-step `deform()` only ever touches a line's two ends, so interior spacing stays regular
on its own during ordinary convergent/divergent motion -- but a *transform* boundary shears
nodes along a line without inserting or deleting anything, which can leave spacing uneven.
At the end of every single `deform()` call (not periodically -- unlike the earlier
`REGULARIZE_INTERVAL_STEPS`-gated cadence this replaced, since `deform()` itself now needs
every line back at (within tolerance of) exact target spacing before the *next* call's
overstretch check can mean anything -- see [Plate motion: shift and
deform](#boundary-evolution)), any line whose gaps have drifted past
`IRREGULARITY_TOLERANCE` (1.5x target spacing, either direction) gets a fresh evenly-spaced
node set across its *existing* extent -- the two endpoints are preserved exactly, since this
never changes where a line's physical edge is, only how regularly it's sampled -- with
elevation re-interpolated onto the new nodes (`np.interp`, 1D since it's along a single
already-ordered curve, not 2D scattered-data interpolation).

<a id="merge-and-split"></a>
## Merge and split (`merge_split.py`)

Unlike `shift`/`deform`, these are rare, discrete, topology-changing events, so a one-time
resample is an acceptable cost here -- the exact, no-resampling guarantee only matters for
routine per-step motion.

- **Consumption.** A plate with no real remaining territory is simply dropped from
  `world.plates` (`remove_defunct_plates`, run every step). Three shapes count, all falling
  directly out of `deform()`'s own grow/shrink rule -- no special algorithm needed:
  - every elevation node deleted (fully subducted);
  - eroded/subducted down to a single remaining line -- a sliver along one latitude,
    regardless of how many nodes are still on it;
  - a *comb of stubs*: many lines but fewer than two nodes each on average. `deform()`
    shrinks a line only from its ends and never deletes its last node, so a heavily
    subducted oceanic plate decays into a hundred-plus rows of one stranded node apiece --
    a high line count masking that there's no 2D patch left. `has_negligible_territory`'s
    node-to-nonempty-line ratio test catches this; the old "at most one line" check didn't.
    The ratio is scale-free (same at any `node_density`) and sits far below any legitimate
    plate, whose rows carry tens of nodes.
- **Continental collision merge.** If two continental plates have at least
  `MERGE_MIN_CONTACT_NODES` node pairs within `MERGE_CONTACT_DISTANCE_RAD` of each other
  *and* a real closing rate at those points (`boundary.closing_rate` -- the one place this
  velocity-based check still runs; see [Plate motion: shift and
  deform](#boundary-evolution) for why ordinary per-step evolution no longer uses it), they're
  fused: keep one plate's
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
  million years it took), a plate disappearing (consumption, in any of the senses above), a
  split creating a new plate, or a plate fragmenting into disconnected landmasses / shedding
  stranded nodes (see Defragmentation below) -- which `world.step_world` timestamps with the *post-step*
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
  and `SPLIT_MIN_NODES` keeps the check off small fragments entirely.

  **Tuning (2026), for the torque engine.** With `torque.py` driving `omega` a large
  continental plate gets a *good* rigid-rotation fit, so the RMS-residual gate at its
  original 9 cm/yr was never tripped -- ordinary supercontinent-scale plates simply never
  rifted, and with `maybe_split_plate` the only source of new plates, plate count decayed
  monotonically as oceanic plates subducted. The gates were loosened
  (`SPLIT_RMS_RESIDUAL_THRESHOLD` 9->6 cm/yr, `SPLIT_MIN_POLE_SEPARATION` 6->4,
  `SPLIT_MIN_NODES` 1200->700, `SPLIT_SIZE_CERTAIN_RIFT_RAD` pi->2.2 rad,
  `SPLIT_MIN_AGE_STEPS` 15->20), and -- like the collision merge -- `apply_topology_changes`
  now performs **at most one split per step**, so a freshly-generated world (every plate
  clearing its identical cooldown on the same step) staggers its rifts over time instead of
  shattering all at once. A 200-step reproduction goes from 12 plates decaying to 9, to a
  healthy churn oscillating ~18-26.
- **Defragmentation.** `deform()` only ever grows or shrinks a line's *ends*, and never
  deletes its last node -- so subduction or transform shear can carve one plate's node
  cloud into two (or more) fully disconnected landmasses, still carried as a single
  `Plate`, or leave a comb of stranded one-node rows trailing behind it. The split check
  above doesn't catch this: it cuts on mantle-flow *disagreement*, and two lobes of one
  plate are co-moving by definition, so their flow samples never diverge. `deform()`
  doesn't catch it either -- an interior node only reads as contested if it's inside a
  *neighbour's* polygon, and a gap between two lobes of the same plate belongs to nobody.

  So `defragment_plates` (`merge_split.py`, every `DEFRAG_INTERVAL_STEPS` steps -- a
  whole-world k-d-tree pass, cheap but not free, and topology doesn't fragment fast) checks
  it directly. Each plate's nodes are grouped into connected components at
  `DEFRAG_CONNECT_RADIUS_MULT * line_spacing_rad` (`plates.node_components`, via
  `scipy.sparse.csgraph.connected_components` over a k-d-tree radius graph -- ~2.5x the
  world's own line spacing links a genuinely contiguous patch while still separating two
  lobes across a real subduction gap). Every component with at least
  `DEFRAG_FRAGMENT_MIN_NODES` nodes (an area, so it scales with `node_density` directly,
  same as `SPLIT_MIN_NODES`) becomes its own plate; the largest keeps the original plate's
  id, `frame`, `omega`, and age, and the rest are fresh plates carrying a *copy* of that
  `omega` (they were co-moving) and age 0, drawing ids from `World.next_plate_id`. Anything
  smaller is dropped as stranded crust. Nodes are repartitioned by boolean mask through the
  same `ElevationLine.masked` machinery `split` uses, so every persistent per-node field
  survives exactly with no resample.

  A plate that's *all* small components -- no lobe big enough to anchor a plate -- is left
  alone here, not deleted; the comb-of-stubs branch of consumption (above) prunes it on the
  same step. This runs before the collision and split passes so a severed lobe or a ghost
  comb stops polluting neighbour polygons and collision detection first.

<a id="gap-filling"></a>
## Whole-sphere coverage (subsumed into `deform()`)

This used to be a separate periodic pass (`gaps.py`'s `fill_gaps`, since removed): a
whole-sphere lattice sweep every `elevation_lines.REGULARIZE_INTERVAL_STEPS` calls, finding
every point farther than a coverage radius from any plate's nearest node and resolving each
cluster by absorbing it into a dominant (or young) bordering plate, or spawning a brand new
plate if no plate dominated.

In the polygon-based model, "an uncovered gap" and "unclaimed territory" are the same
concept `deform()` already computes every turn for every plate (the sphere minus every other
live plate's own bounding polygon -- see [Plate motion: shift and
deform](#boundary-evolution)), so there's no separate detection pass left to run: a plate
growing toward its own pole, or reclaiming ground a subducted neighbor just vacated, is just
the ordinary "claiming adjacent territory" sub-step of `deform()`, described there.

**One real behavior this doesn't reproduce**: the old pass's "no plate dominates the gap's
border -> spawn a brand new plate" fallback. `deform()`'s claim step only ever grows an
*existing* plate into a whole new row adjacent to its own current phi extremes -- it has no
mechanism to spawn an entirely new plate for a gap that borders no existing plate's reach at
all. In practice this should be rare (the sphere starts fully tiled with no gaps by
construction, see [Initial plate generation](#initial-plate-generation), so a gap can only
open next to whichever plates used to border the space that vacated it), but it's a known,
not-yet-addressed gap in coverage rather than a deliberately preserved behavior.

<a id="volcanism"></a>
## Volcanism (`volcanism.py`, plus `plates.py`'s own `PlateWithLines.deform`)

New continental crust forming where plates are separating. This used to be two halves in one
module: a periodic whole-sphere *detection* pass that spawned brand-new "volcanic field"
plates, plus a per-step *eruption* pass. The detection half is gone -- volcano creation is
now inline in `deform()`'s own overstretched-rift handling (see [Plate motion: shift and
deform](#boundary-evolution)) -- and `volcanism.py` now holds only the per-step eruption
lifecycle, which is unchanged.

**Creation.** When a divergent line end grows, the new nodes come back as a fresh volcano
(rather than plain ridge/rift fill) with a small fixed probability per growth event
(`STRETCH_VOLCANO_PROBABILITY`, 0.02) -- see [Plate motion: shift and
deform](#boundary-evolution) for why this replaced the old whole-sphere gap-outlier scan
(`GAP_OUTLIER_FACTOR`, boundary-point median-spacing comparison) and why it's a probability
roll rather than a deterministic threshold. A volcano node is a node *on the growing plate's
own existing line* -- not a separately spawned `Plate` the way the old detection pass worked.
One consequence: there is no longer any whole-plate "volcanic field" bookkeeping at all --
`World` used to carry a `volcanic_field_plate_ids` set and `volcanism.py` a dormancy check
that relabelled a diluted field as ordinary continental crust, but with nothing creating a
separately-tracked field plate both were dead weight and have been removed. A fresh volcano
gets a random `volcano_active_years_remaining` draw
(`VOLCANO_ACTIVE_MIN/MAX_YEARS`, 100k-1M years) and **one guaranteed immediate eruption**
(`+= ERUPTION_ELEVATION_M`, unconditional -- unlike every later eruption, which is rolled
probabilistically, see below) -- distinguishing a freshly-created rift volcano from an
ordinary volcano's own first, merely-probable eruption.

**Eruption, every step, unchanged.** Each individual volcano point has its own
`volcano_active_years_remaining` (drawn once at creation, whether by `deform()`'s rift
handling above or, previously, the old detection pass), decremented every step. While
active, it rolls a per-step eruption chance
(`1 - exp(-ERUPTION_RATE_PER_MYR * active_years_this_step / 1e6)`, the same
exponential-arrival-rate shape used elsewhere in this codebase, e.g. lake evaporation's own
retention factor -- expected roughly 0.3 to 3 eruption events over a volcano's own full active
life, "occasionally," not every step) and, if it erupts, adds `ERUPTION_ELEVATION_M` (100m) of
new land and grows `mineral_deposit_m` (see [Resources and soil](#resources-and-soil)).
Deterministic per `(seed, elapsed_years, plate_id, line_index)`, the same reproducibility
precedent merge_split.py's own per-pair collision threshold sets. `active_years_this_step` is
clamped to the volcano's own *remaining* life, not the full step size -- a large step (the UI
offers up to 10 Myr) shouldn't roll eruption chances for years past when a short-lived
volcano actually went dormant.

**Rendering.** Baked directly into the elevation/plates views' raster the same way lakes and
glaciers are (`VOLCANO_COLOR_RGB`, a hot red-orange distinct from both), drawn after lake but
before glacier so ice still wins where both would apply (a volcano cold enough to glaciate
should read as ice-covered, not lava-red).

<a id="reassignment"></a>
## Boundary point reassignment (subsumed into `deform()`)

This used to be a separate periodic pass (`reassign.py`, since removed): ordinary per-step
boundary evolution only ever grew or shrunk a line's two *ends*, never revisiting whether an
*interior* node still actually belonged to the plate carrying it, so a node that drifted
into a neighboring plate's own territory (enough shearing along a transform boundary, or
slow rotational drift) could sit there unnoticed until a periodic whole-world scan caught it.

In the polygon-based model this class of bug can't accumulate unnoticed in the first place:
every plate's `deform()` call, every single turn, directly tests every one of its own
near-boundary nodes for containment in a neighbor's current polygon (see [Plate motion:
shift and deform](#boundary-evolution)) -- a node that's drifted into a neighbor's territory
reads as contested the very next turn and gets removed by the ordinary shrink rule, not
waiting for a periodic pass every `REASSIGN_INTERVAL_STEPS` calls. There's no equivalent
"move this node onto the neighbor's own nearest line" step, since a contested node is
deleted (and, on the neighbor's own next turn, its territory is either already covered by
ordinary growth or gets reclaimed via "claiming adjacent territory") rather than transferred
node-for-node -- a coarser mechanism than the old point-relocation logic, but one that keeps
the same underlying invariant (every node ends up owned by whichever plate's polygon
actually contains it) without needing a separate whole-world pass at all.

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
(`plates.iter_local_lattice(identity_frame, spacing_rad=GRID_SPACING_RAD)` -- reusing
`iter_local_lattice` with the identity frame as a plain global lat/lon sweep, rather than any
one plate's own local frame), assign every cell its nearest
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

**The selected-plate panel also reports motion and shape-health diagnostics** (added
2026-08-31, from the long-run plate-geometry investigation in `docs/TODO.md`): the plate's
surface speed in cm/yr and whether it is railed at `mantle.MAX_PLATE_RATE`, its Euler pole,
its age in steps, its median node elevation and submerged fraction (flagged when a
continental plate is mostly under water), the other plates its territory currently overlaps
(and by what share of its own nodes -- one global `cKDTree` pair query in
`main._plate_overlaps`), and any `world.collision_progress` sustained-collision timers it is
part of. These are pure read-outs of existing state (`torque.py`'s `Plate.omega`,
`world.collision_progress`), surfaced because "the plates look wrong on a long run" was
otherwise only diagnosable with an ad-hoc script. See `docs/api-reference.md#get-worldplates`.

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
whose mechanisms are richer than a first-principles description would suggest. Of those
seven, only **wind and air temperature** are read straight off `World.atmosphere_cfd_state` --
a real, continuously time-integrated shallow-water solve (see [Atmospheric wind
solver](#ocean-atmospheric-fluid-dynamics)) -- rather than reconstructed by this module's own
formulas on every call. Land temperature stays a genuine per-call diagnostic (it's also the
CFD state's own equilibrium relaxation target, computed the same way). Ocean currents, ocean
surface temperature, humidity, and precipitation are all per-call diagnostics too: the
shallow-water *ocean* solver had no stable operating point on this grid and was retired
(`compute_ocean_currents` + `advect_ocean_temperature` run here every step instead, fed by the
CFD wind), and the CFD atmosphere's own humidity/precipitation had no orographic-lift,
moisture-convergence, or lake/river/vegetation term and produced a near-zero, uncalibrated
rainfall field that starved erosion and hydrology. The precipitation sweep is at steps 9-10
below. Rivers, lakes, and
vegetation feed back into humidity here (see
step 9 and "Moisture recycling") -- the one place this module reaches outside its own fixed
grid, into the persisted node fields hydrology.py/erosion.py/biomes.py already maintain. The
precipitation field is a steady-state advective sweep, so it already reads as long-term
average rainfall rather than an instantaneous rate.

**A third, genuinely fixed-shape grid, used only here.** Elevation is Lagrangian (see [Why
not a grid](#why-not-a-grid)); the render grid ([Render image](#render-image)) is a *ragged*
lat/lon sweep, immediately flattened to 1D. Neither supports the array tricks climate leans
on -- `np.roll` wraparound, centered-difference gradients, divergence, land-excluding
neighbor averaging. So climate gets its own equirectangular array, `lat: (H,)` / `lon: (W,)`,
every field `(H, W)`, `GRID_HEIGHT = 90` x `GRID_WIDTH = 180` (2 degrees/cell) at the
`climate_density = 1.0` reference -- scaled directly in each dimension by a world's own
`World.climate_density` (the UI's "climate & biome resolution" choice,
`climate.CLIMATE_DENSITY_CHOICES = (0.5, 1.0, 2.0, 4.0)`, `climate.DEFAULT_CLIMATE_DENSITY =
4.0`, see `climate.grid_dimensions`) for a sharper, less pixelated climate map at the cost of
real compute: confirmed directly, `2.0` (4x the cells) costs roughly 2-3x longer to render a
climate/biome-family view and only a modest ~10% longer per simulation step (climate is a
comparatively small share of a step's total cost next to boundary evolution/erosion/hydrology
over the plate node cloud); the default `4.0` (16x the cells) costs roughly 6x longer to
render such a view and roughly 1.5-2x longer per simulation step; `0.5` (a quarter of the
cells) is the corresponding lower-resolution option, for a coarser climate map and slightly
faster steps. Fixed-*degree* offset distances internal to the wind/current
mechanics below (mountain/coastal wake lookback, mountain tangent sampling) stay physically
meaningful at any density via `_REFERENCE_CELL_DEG`, a reference cell size decoupled from the
grid's own actual width -- one exception, noted where it's defined
(`MOUNTAIN_GRADIENT_SMOOTHING_ITERATIONS`): a handful of *iteration*-count-based smoothing
passes blur by a fixed number of grid cells rather than a fixed real distance, so their
real-world smoothing radius shrinks proportionally at a higher density -- a deliberate,
smaller-scope simplification, not rescaled. This grid is never stored on `World` and never
touches `world.plates`.

**Terrain-derived; wind/currents/temperature CFD-sourced, humidity/precipitation diagnostic.**
Elevation/is_ocean are resampled from scratch on every call, from whatever the *current* plate
elevation/crust_type happens to be (the same `cKDTree` nearest-neighbor sampling
`_render_grid_arrays` already uses) -- there's no persistent terrain field of climate's own to
keep in sync, since terrain itself already persists incrementally on the plates. Temperature/
wind/currents are read off the world's own always-on, genuinely prognostic CFD states (see
above); humidity/precipitation are recomputed by this module's own steps 9-10 every call, from
that CFD-sourced wind/ocean temperature. Either way climate is genuinely recomputed every step
regardless of whether a climate view is currently being rendered (see [Erosion](#erosion),
which needs a live climate snapshot for rain/wind/humidity every step), not just on render.
`compute_climate`'s `skip_moisture=True` path (used only by `world._advance_fluid_dynamics`'s
CFD-forcing call, which consumes only elevation/is_ocean/the temperature baselines) skips the
humidity/precipitation sweep and returns those two as zeros.
Two more generation-time-fixed exceptions: `World.axial_tilt_deg` and `World.climate_density`,
both set once by `generate_world` and read again on every future step/render, since insolation
needs the former and the grid's own shape needs the latter long after generation, neither
being something to recompute per call.

**Pipeline order.** Wind needs a temperature field, but the *final* (current-advected) ocean
temperature needs currents, which need wind -- resolved by computing a pre-advection baseline
first and closing the loop only for the final consumer-facing fields:

1. **Insolation** -- `cos(lat)` zenith-angle law (clipped floor) at full sun (`SUNLIGHT =
   1.0`), plus axial tilt: with tilt 0 it's the flat law; otherwise the mean of that same law
   over `AXIAL_TILT_DECLINATION_SAMPLES` declinations swept between `-tilt` and `+tilt` (the
   sub-solar latitude's annual sweep -- this model has no calendar, so it's an annual mean,
   not a season cycle). `World.solar_multiplier` scales it live.
2. **Land temperature** = `LAND_TEMP_MIN_C + LAND_TEMP_RANGE_C * insolation` (`-28 + 55*(·)`,
   calibrated so the equator-to-pole profile lands near `25 / 10 / -1 / -22 C` at lat
   `0 / 45 / 60 / 90` -- an earlier `-60 + 95*(·)` at a dimmed sun put the equator at ~18C
   and everything poleward of ~40 degrees below freezing) minus elevation-based lapse-rate
   cooling, kept as part of the same base-heating formula rather than a separate causal
   channel (mountains being cold is a consequence of solar heating at altitude).
3. **Ocean temperature baseline** = `WATER_TEMP_MIN_C + WATER_TEMP_RANGE_C * insolation`
   (`-2 + 30*(·)`), a narrower range plus a freezing floor (water's greater thermal inertia)
   -- pre-advection.
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
   *plus* a local land-surface moisture source (lake/river evaporation and vegetation
   transpiration -- see "Moisture recycling" below), advected onto land by a wind-driven 2D
   sweep: a zonal pass and a meridional pass (each a sequential flow-direction walk,
   vectorized across the perpendicular axis, single-column evaporation/retention/orographic-
   dump step per iteration -- pure numpy, no `numba`), blended per-cell by each wind
   component's share of total wind magnitude. The zonal pass's sweep *direction* is the same
   fixed latitude-band lookup wind's meridional structure uses (`zonal_direction_for_lat`),
   not the literal local wind sign. Land-cell moisture (advected carry-over plus local source)
   is capped at the same ceiling ocean evaporation saturates at -- without that cap, a long,
   uniformly-sourced land stretch (a continent-spanning rainforest belt) would compound its
   own source additively, cell after cell, toward an asymptote many multiples of any
   physically sensible humidity value; confirmed directly during development as a genuine
   runaway (mean land precipitation still climbing after ten simulated steps, land humidity
   exceeding the ocean's own maximum) before this cap was added.

   **Inland decay is a real per-km half-life, not a per-cell constant.** Each sweep step's
   carried-over moisture is discounted by `climate._retention_factor`, whose flat-land base
   rate (`_zonal_base_retention`/`_meridional_base_retention`) implements the rule of thumb
   that rainfall roughly halves a few hundred km inland (`MOISTURE_HALVING_DISTANCE_KM = 380`
   -- a little longer than the ~300 rule-of-thumb, since the CFD-sourced wind carries moisture
   inland less efficiently than the old diagnostic trade-wind field did, and a strict 300 dried
   continental interiors to near-total desert on the now-warmer planet) -- computed from each step's *actual* physical
   distance (`plates.PLANET_RADIUS_KM`-scaled degrees, `cos(lat)`-corrected for the zonal
   pass, since a step's real longitude distance shrinks toward the poles) rather than a fixed
   per-grid-cell multiplier. This matters because a fixed-per-cell retention silently decays
   inland moisture at a rate that depends on how much real distance one grid step happens to
   span -- which varies with `World.climate_density` (a lower-density world's larger cells
   cover 300km in barely more than one step) and, for the zonal pass, with latitude. An
   earlier fixed `RETENTION_PER_CELL = 0.96` was confirmed directly as the cause of humidity
   and precipitation reading *highest in continental interiors* rather than near coasts: at
   coarser densities a single step barely dented moisture even 300km inland (~95% retained,
   not the physically-expected ~50%), so land-locked interiors sat near the ocean's own
   evaporation ceiling indefinitely, before any local land-surface source was even added on
   top -- unlike real Earth, where roughly 70-80% of precipitation falls over ocean and only
   20-23% over land, concentrated near coasts and windward mountain slopes rather than spread
   evenly (or worse, favoring the interior) across a continent. On-terrain wind slowdown
   (`elevation_factor`, from own-cell elevation and mountain wake) still applies as an
   *additional* multiplicative discount on top of this distance-based base rate, unchanged
   from before -- air moving slower through rough terrain loses proportionally more moisture
   to mixing independent of the distance it covered.

   **Coherent noise breaks the zonal banding.** The ocean evaporation ceiling is a clip of
   the insolation-driven ocean-temperature baseline, so it's very nearly a pure function of
   latitude, and both sweeps inherit those flat parallels -- visible as horizontal banding on
   the humidity map and in the precipitation baseline it scales. The same fractal gaussian
   perturbation the moisture-flux convergence uses (`_coherent_noise` -- a scale-free
   `1/k**MFC_NOISE_SPECTRAL_BETA` spectral field, correlation length
   `MFC_NOISE_CORRELATION_DEG`, seeded per world state on its own RNG stream) is added to the
   blended humidity field at std `HUMIDITY_NOISE_STD` (smaller than `MFC_NOISE_Q_STD` -- it
   perturbs the baseline everywhere, not just a belt edge), re-clipped to
   `[0, MAX_EVAPORATION_CEILING]`. The orographic dump, a separate rain-shadow channel, is
   left untouched.
10. **Precipitation** = f(humidity) + an orographic bonus (continuous saturating
    windward-slope moisture dump, from wind blowing up-elevation) + a **Hadley/Ferrel
    moisture-flux convergence** term (`compute_moisture_flux_convergence`): the
    metric-correct spherical divergence `-div(humidity * wind)` of the CFD wind field,
    smoothed and blended toward its own per-latitude-row mean (`MFC_ZONAL_COHERENCE`) so the
    wet/dry belts read as continuous bands. Converging moisture-bearing wind adds rainfall
    (the ITCZ where the two hemispheres' trades meet, the ~50-65 deg sub-polar front);
    diverging wind multiplies the humidity baseline *down* (the subtropical highs at ~15-35
    deg, the poles under the polar cells). This is the model's zonal precipitation
    climatology -- but emergent from the winds, not a hardcoded function of latitude, so the
    bands shift and distort with the circulation (a supercontinent, an off-centre landmass, a
    monsoon dragging the ITCZ poleward over a summer continent). The convergence bonus is
    capped (`MFC_CONVERGENCE_MAX_Q`) for the same reason the humidity sweep is capped at
    `MAX_EVAPORATION_CEILING`: to keep the precipitation -> rainforest -> transpiration ->
    humidity recycling loop from running away along an equatorial forest belt.
    Because the row-mean blend leaves each latitude nearly flat, the belt edges would
    otherwise fall on hard parallels and show as horizontal banding in the map; a fractal
    gaussian perturbation (`_coherent_noise` -- white noise shaped by a
    `1/k**MFC_NOISE_SPECTRAL_BETA` power spectrum in the spectral domain, so it's isotropic
    and scale-free rather than one dominant bump size, std `MFC_NOISE_Q_STD`, correlation
    length `MFC_NOISE_CORRELATION_DEG`, seeded per world state) added after the blend warps
    and stipples those edges the way real convergence zones meander, without moving the
    wind-derived per-row means.
    Feeds erosion and hydrology (see [Erosion](#erosion) and [Hydrology](#hydrology)), and
    -- one step later, via the moisture-recycling source above -- itself.

**Moisture recycling: rivers, lakes, and vegetation release moisture too** -- the "rain in a
rainforest" effect, where a wet, densely-vegetated region partly sustains its own
precipitation. `compute_humidity` samples `plates.ElevationLine`'s own persisted `lake_depth`/
`channel_depth` fields onto the climate grid (the same nearest-neighbor resample elevation
itself already gets, see `_sample_elevation_and_crust`) to size a lake-evaporation and a
river-evaporation source; a third source, vegetation transpiration, comes from
`biomes.classify_biomes` (`VEGETATION_TRANSPIRATION_BY_BIOME`, a hand-picked weight per Köppen
class peaking at the wet-tropical / oceanic classes, near zero at ice/tundra/desert classes
and at every pelagic ocean class) applied to a *flat* classification with a broadcast latitude
row (continentality via a distance transform on `is_ocean`). Classifying vegetation needs a precipitation value, and
this step's own precipitation is exactly what transpiration itself feeds into -- so, like
every other circular coupling in this module (e.g. ocean current <-> temperature), the loop is
broken with a one-step lag: vegetation is classified from `World.climate_cache`, last step's
already-cached snapshot, not this step's still-unknown one. Zero everywhere on a world's very
first call (`world.climate_cache is None`), self-correcting after one step. A frozen surface
(`air_temperature_c` below `hydrology.FREEZE_POINT_C`) can't evaporate, so lake/river
evaporation -- but not transpiration, already near zero in any biome cold enough to freeze --
is zeroed there.

**Transpiration recycles last step's own rainfall -- it can't manufacture moisture from
nothing.** `_vegetation_transpiration_source` scales `VEGETATION_TRANSPIRATION_BY_BIOME`'s
per-biome weight by `VEGETATION_RECYCLING_FRACTION` (0.16 -- below the real Amazon-basin
quarter-to-half range, tuned down as the climate model got warmer, see the constant's own
comment) *and* by that same cell's own `prev.precipitation_mm` from `World.climate_cache`
(converted back to humidity units) -- a fraction of what actually fell there last turn, not a
flat per-biome constant. An earlier version used a fixed `VEGETATION_TRANSPIRATION_MAX` source
keyed only on biome type, with no real reservoir behind it -- confirmed directly as a real
multi-step runaway: mean land precipitation climbed turn after turn with no equilibrium,
eventually overtaking the ocean's own total. Tying the source to a fraction of real, already-
fallen rain gives the land/ocean split a genuine fixed point instead: stepping two seeds 65-80
turns past generation settles the land precipitation share near ~24% (real Earth is ~20-24%),
rather than climbing without bound toward a land-dominated one.

This also decreases the standing water it comes from, per the same request that asked for the
recycling effect in the first place: lake evaporation already shrinks `lake_depth` (lakes.py's
own water balance, a separate, unrelated mechanism computing the same physical process); river
evaporation is a new loss along `hydrology.route_downstream` (`RIVER_EVAPORATION_*`), reducing
`flow_accum` as water flows through a warm, dry stretch. That evaporation loss can't be the
same one climate.py's own source term reads, since hydrology.py runs *after* climate.py each
step (consuming its precipitation/temperature output) -- climate.py's land-surface term is
necessarily a same-step stand-in sized from the persisted channel_depth instead, not the
literal amount hydrology.py evaporates this same step.

**Scope, explicitly decided.** Kept out: river outflow feeding currents, deep currents.
Included, even though richer than a one-line causal description: axial tilt, wind's mountain
deflection/Venturi/wake, ocean currents' coastal deflection/land swirl/circumglobal
boost/wake, humidity's evapotranspiration term and lake/river climate influence, and
precipitation's zonal wet/dry banding -- the last modelled as the moisture-flux convergence
of the CFD wind field (step 10 above), not the hardcoded latitude climatology an earlier
design had ruled out. All implemented in full, not simplified down. Dropped outright: river
outflow feeding currents, the one river/lake/vegetation-climate coupling still not modeled.

**Rendering.** `render_image.py`'s `CLIMATE_VIEWS` (`temperature`, `wind`, `oceanCurrents`,
`humidity`, `precipitation`, `biome`) route to `_render_climate_view`, a separate path from
the plate-tectonics views since the data source is structurally different from the render
grid's ragged lattice. `temperature`/`wind`/`oceanCurrents` draw the native `(npix,)` CFD
state arrays directly; `humidity`/`precipitation`/`biome` take `compute_climate_cached`'s
diagnostic `(H, W)` fields and resample them nearest-cell onto the HEALPix cells.
Heatmap views (temperature/humidity/precipitation) reuse the elevation view's
color-stop-interpolation technique with their own stop tables; wind/ocean-currents draw
subsampled arrows (numpy-vectorized projection/direction math, looped only for the
unavoidable per-arrow PIL draw calls), and ocean currents additionally marks each sampled
swell point with a small circle. Temperature/humidity/precipitation additionally draw the
current coastline (see [Coastline](#coastline)) -- a color-scale view carries no land/ocean
information on its own, unlike elevation's own hypsometric coloring. `biome` is categorical,
not a heatmap -- see [Biomes](#biomes) below -- and (like wind/oceanCurrents) skips the
separate coastline stroke, since its own flat Ocean color already reads as a land/ocean
boundary on its own.

<a id="biomes"></a>
### Biomes (`biomes.py`)

A pure, stateless classification -- `biomes.classify_biomes(temperature_c, precipitation_mm,
elevation_m, slope, is_ocean, sea_level_m, *, lat_deg, axial_tilt_deg, continentality=None,
dist_to_land_rad=None, has_sea_ice=None)` -- assigning each land cell one of the **31
Köppen-Geiger** 3rd-level classes (`classify_koppen`) and each ocean cell one of **10 pelagic
classes** (`classify_pelagic`), with `is_ocean` settling the land/water split. `BIOME_NAMES` /
`BIOME_COLORS` are the Köppen classes first (descriptive names -- "Humid Subtropical", not
"Cfa"), then the pelagic classes; `OCEAN_IDS` is the id set of the latter.

Real Köppen keys off sub-annual quantities this model never produces -- coldest/warmest-month
temperature and the summer/winter precipitation split -- so `biomes.py` **synthesizes** them:
`_seasonal_temp_amplitude` drives a mean-to-peak seasonal swing from `|lat|`, continentality
(distance inland, 0 at the coast), and `axial_tilt_deg` (a tilt-0 world gets amplitude 0, so
its `s`/`w`/`d` subtypes never occur); `_precip_season` drives a summer precipitation share
and a seasonality concentration from latitude (monsoon belt summer-wet and peaked, a narrow
subtropical dry-summer dip, mid-latitudes even), scaled by the same tilt factor. All of these
constants are visually tuned like the codebase's other openly-approximate values (e.g.
erosion.py's `RAIN_EROSION_COEFFICIENT`), calibrated so an Earth-like world lands its major
zones roughly where Earth's are -- latitude alone can't tell a west coast (Mediterranean) from
an east coast (humid) at the same latitude, so `s` is further gated on maritime
continentality and modest precipitation, and some regimes (the East-Asian-monsoon `Dw` belt)
read as their `f` sibling instead. `ICE_TEMP_C` still reuses
`hydrology.GLACIER_ACCUMULATION_TEMP_C`.

**Pelagic classes** are PPOW's abiotic hierarchy rather than its Earth-geographic province
names (which can't transfer to a different planet): a thermal realm from sea-surface
temperature (polar / cold-temperate / temperate / subtropical / tropical) crossed with a
structural zone -- sea ice (`has_sea_ice` or SST `< -1`), coastal shelf
(`dist_to_land_rad <= SHELF_RANGE_RAD`, matching `geology.SHELF_RANGE_RAD`), equatorial
divergence (`|lat| <= 6` in the tropics), or open ocean / subtropical gyre.

`smooth_biome_field` (the version the map views and `/world/stats` use) computes
`continentality` and coast distance from `is_ocean` itself via a distance transform, plus
sea-ice cover from a passed glacier-depth grid, then runs the same stateless neighbour-vote
cleanup as before on land cells only. `climate.compute_climate` calls it once per computation
and stores `ClimateFields.biome_ids` (see [Climate](#climate)); `stats.py`'s
`biome_land_fraction` (which excludes every `OCEAN_IDS` id) reads that stored field, while the
map-view / hex-export / vegetation-transpiration call sites still classify against their own
grid.

**Wetland and Carboniferous Forest** are no longer displayed climate classes (Köppen has no
such category), but `classify_wetland` and its constants stay: `geology.py`'s per-node coal
formation calls it directly (see [Resources and soil](#resources-and-soil)). Both require
flat, low-lying land (`elevation_m > 0`, `<= WETLAND_MAX_ELEVATION_M`,
`slope <= WETLAND_MAX_SLOPE`); Carboniferous Forest is the warm
(`>= CARBONIFEROUS_MIN_TEMP_C`, aliased to `TROPICAL_TEMP_C`), very wet
(`>= CARBONIFEROUS_MIN_PRECIP_MM`, aliased to `HUMID_MM`) subtype.

**Slope** (`biomes.grid_slope`, still used by `classify_wetland` and `geology.py`) -- real
elevation difference to each cell's steeper of its north/south or east/west neighbor, divided
by that neighbor's real great-circle spacing in meters (longitude narrowed by `cos(lat)`, the
same convention `plates.iter_local_lattice` uses) -- computed on the fine Biome/Combined grid
(`_biome_fields`'s own `elevation_m`, previously unused by these two views) rather than
climate.py's coarser native grid. A different, coarser discretization than
`erosion.compute_slope`'s own node-cloud slope (used by `geology.py`, see below), but both are
genuine rise/run slopes over roughly comparable real distances (~100-125km), so
`WETLAND_MAX_SLOPE` is a shared, visually-tuned cutoff rather than two independently-fit
thresholds -- `stats.py`'s own `biome_land_fraction` reuses `grid_slope` too, against its own
coarser native climate grid, for the same reason.

No new per-step state or caching -- like `render_image.py`'s own `temperature_colors`/
`humidity_colors`, this runs fresh from whatever `climate.compute_climate_cached` already
produced, entirely inside `_render_climate_view`'s `"biome"` branch (and `_render_combined_view`
for "Combined", though Ocean/Intertidal Zone's own biome color is never actually visible there
-- see that view's own legend). Implemented with `np.select` (first-matching condition wins)
rather than chained `np.where` overwrites, so each band's cutoffs stay a self-contained,
independently checkable list instead of depending on write order to get boundary cells right.

**Combined-view shading.** `_render_combined_view` multiplies each land cell's flat Köppen
color by `biomes.biome_relative_shade_factor` -- a *continuous* brightness ramp from
`1 - BIOME_SHADE_AMPLITUDE` at a class's lowest-elevation cell to `1 + BIOME_SHADE_AMPLITUDE`
(±25%) at its highest, linear in that cell's elevation *rank among its own class's cells only*
(so a class pinned to a narrow absolute elevation band still spans the full range). Near real
peaks the result is further blended toward the elevation gradient (`RELIEF_BLEND_MAX`); ocean
cells are the pelagic-province color blended toward the hypsometric depth shade
(`OCEAN_PELAGIC_RELIEF_BLEND`), so deep basins still darken and shelves lighten while every
water cell carries its province's hue. Because that wide color spread means a pixel's *color*
no longer identifies its class, the render is **RGBA** and the per-pixel class/lake/glacier id
is carried in the alpha byte (`alpha = 255 - code`, every classified land and ocean cell now
carrying one; see `render_image.COMBINED_LAKE_ID_CODE`). `frontend/src/legendData.ts` mirrors
that id mapping by hand -- grouping the 41 classes into the ~17 1st/2nd-level legend rows --
and `MapCanvas.tsx` reads it straight off alpha for legend-click-to-highlight, then resets
alpha to opaque before display.

<a id="ocean-atmospheric-fluid-dynamics"></a>
## Atmospheric wind solver (`fluid_dynamics.py`, `atmosphere_cfd.py`)

Most of [Climate](#climate) is *diagnostic* -- recomputed from scratch every call via
latitude-banded heuristics (Ekman currents, coastal deflection, mountain deflection),
not a real fluid solve. **Wind is the exception**: `World.atmosphere_cfd_state` is a genuine,
time-integrated **shallow-water** solve -- real Coriolis/pressure-gradient momentum physics,
explicit numerical time-stepping, real *prognostic* `u`/`v`/`eta`/`temperature_c` that persist
and evolve step to step (rather than being rebuilt every call the way `climate_cache` is), plus
a sustained latitude-banded forcing so the planetary circulation doesn't spin down.

**A shallow-water *ocean* solver (`ocean_cfd.py`) also existed** -- prognostic currents,
surface temperature, and a sediment tracer. It was **retired**: on this HEALPix grid it had no
stable operating point that produced realistic circulation (either an over-damped near-still
interior or a coastal grid-scale instability), and its sediment views had already been dropped.
Ocean currents and the current-advected surface temperature are diagnostic in `climate.py`
every step now (`compute_ocean_currents` + `advect_ocean_temperature`, see below), fed by the
CFD-solved wind.

<a id="mode-toggle"></a>
### Always-on, not a mode

Earlier revisions gated Fluid Dynamics behind a three-way `World.fluid_mode` toggle
(`POST /world/mode`) that froze plate tectonics/climate while it was active. That's gone:
`World.atmosphere_cfd_state` is created once, by `generate_world`
(`atmosphere_cfd.init_atmosphere_cfd`, seeded from the diagnostic `compute_wind` bootstrap),
and never re-initialized again for the rest of that world's life.

**Every `POST /world/step` call advances it**, via `step_world`'s own
`_advance_fluid_dynamics`, right alongside plate tectonics/climate/erosion -- gated on
`World.simulate_climate_biomes` the same way erosion/hydrology already are (not on
`simulate_plate_movement`; wind keeps evolving even with plate movement paused), *and* on
`World.wind_model == "cfd"` -- the diagnostic wind model (see [Wind model](#wind-model)
below) turns `_advance_fluid_dynamics` into a no-op.
Real atmospheric fluid dynamics needs timesteps of hours to days; plate tectonics needs
timesteps of thousands to millions of years -- reconciled by advancing the wind state by a
**fixed real-time increment per tectonics step, regardless of the tectonic `years`
requested**: `atmosphere_cfd.SECONDS_PER_TECTONIC_STEP` (one simulated day). A single
`/world/step` call covering a million tectonic years still only advances the atmosphere by one
simulated day -- an intentional decoupling.

**`refresh_forcing`** (called right before `step_atmosphere_cfd` each tectonics step) keeps the
state's terrain-derived boundary conditions -- `elevation_m`/`is_ocean`/
`equilibrium_temperature_c` -- in sync with the world's evolving plate state, while leaving the
genuinely prognostic `u`/`v`/`eta`/`temperature_c` untouched, so they keep evolving
continuously across tectonics steps rather than resetting. It's fed a `skip_moisture=True`
`compute_climate` snapshot at `World.fluid_density`'s resolution (the humidity/precipitation
sweep is dead weight for forcing that consumes only elevation/is_ocean/temperature baseline).

`climate.py`'s own `compute_wind` diagnostic still exists as the one-time cold-start
bootstrap `init_atmosphere_cfd` falls back to during `generate_world` (before
`World.atmosphere_cfd_state` exists yet) -- *and* as the full wind source whenever the
diagnostic wind model is selected, below.

<a id="wind-model"></a>
### Wind model: `"diagnostic"` (default) vs `"cfd"`

`World.wind_model` (live-adjustable via the Controls window / `POST /world/controls`) picks
which wind field feeds `climate.py`:

- **`"cfd"`** -- the genuine shallow-water solve described above. At the default
  `fluid_density` this is the single most expensive piece of a `/world/step` (hundreds to
  a thousand+ CFL-stable substeps; see [Performance](#fd-performance)).
- **`"diagnostic"`** (default) -- skip that solve entirely. `_advance_fluid_dynamics` becomes a no-op
  (the CFD state is *kept*, just frozen, so switching back resumes from it rather than a
  cold start), and `climate.compute_climate` rebuilds wind from `compute_wind` and air
  temperature from `compute_air_temperature_diagnostic` every call -- the same closed-form
  "ABL" path the cold-start bootstrap uses, but for the whole run. `compute_air_temperature`
  (the maritime-moderation cold-start air temp) is *not* reused here: it's a poor match for
  the CFD's near-radiative-equilibrium temperature field (benchmarks at ~41% land-biome
  agreement vs the CFD, against `compute_air_temperature_diagnostic`'s ~90%).

**Why it's a reasonable trade.** The CFD wind is relaxed hard toward the same
latitude-banded target `compute_wind` builds, and only integrates ~one simulated day per
step, so downstream (ocean currents, humidity advection, precipitation) barely distinguishes
the two wind fields -- keeping the *CFD* air temperature while swapping to diagnostic wind
holds land-biome agreement near 95%. Measured against a ~12-step CFD reference at
`fluid_density=2.0` across three seeds, the fully-diagnostic model reproduces **~84-89% of
the land biome map** and precipitation/temperature **within ~8-11%**, while cutting a
12-step run from ~42 s to ~7 s (**~6x**) at that resolution -- the CFD substep loop itself
is ~15x, diluted by the rest of a step's fixed cost, and the multiple grows with
`fluid_density` (closer to ~1.5x at `0.5`). The residual gap is almost entirely the
air-temperature field's missing advective/diffusive structure -- see `TODO.md` for the
options to close it.

**Rendering.** The `"wind"` and `"temperature"` map views normally draw
`World.atmosphere_cfd_state`'s native HEALPix arrays; under `"diagnostic"` they fall back to
`compute_climate`'s own wind/air-temperature (resampled onto the HEALPix cells like every
other climate field) so the maps reflect the model that's actually driving the sim rather
than a frozen CFD snapshot.

<a id="shallow-water-formulation"></a>
### The shallow-water formulation

Both solvers run on the same fixed equirectangular grid shape `climate.py` already uses (same
`(lat_deg, lon_deg, world_xyz)` shape convention, built via the same `climate.compute_climate`
pipeline), but at `World.fluid_density`'s own resolution rather than `World.climate_density`'s
-- a separate, independently choosable Advanced-settings option (same `climate.
CLIMATE_DENSITY_CHOICES` set, defaulting to match `climate_density`'s own default so an
unchanged world behaves exactly as before), letting a world keep a sharp climate/biome render
grid while running FD mode at a coarser (faster) resolution or vice versa; see
[Performance and grid resolution](#fd-performance) below. The solve is a two-equation shape
(`fluid_dynamics.py` holds every numerical primitive; the retired ocean solver shared them):

- **Momentum** (`u`, `v` -- east/north velocity, m/s): `du/dt = f*v - g'*d(eta)/dx + forcing_x - drag*u + nu*laplacian(u)`, symmetric for `dv/dt`. `f` is the *real* Coriolis parameter (`fluid_dynamics.coriolis_parameter`, `2*OMEGA*sin(lat)`, `OMEGA = 7.292e-5 rad/s` -- a genuine physical constant, unlike `climate.py`'s own simplified `sin(lat)` proxy, since this is an actual momentum equation) and `eta` a geopotential-height anomaly. `forcing_x` also carries the latitude-banded wind-forcing relaxation `WIND_FORCING_RELAXATION_PER_S * (u_target - u)` (see [Atmospheric Fluid Dynamics](#atmosphere-cfd)).
- **Continuity** (`eta`): `d(eta)/dt = -div(H*(u, v))` plus a thermal-relaxation term (see [Atmospheric Fluid Dynamics](#atmosphere-cfd) below), `H` a fixed troposphere-depth constant.
- **Advection-diffusion** of `temperature_c`, via `fluid_dynamics.semi_lagrangian_advect` -- the same backward-trace-and-sample technique `climate.py`'s own `advect_ocean_temperature`/`_sample_at_offset` already use: unconditionally stable regardless of how far `dt * speed` reaches relative to one grid cell, so advection never adds a second, stricter stability constraint on top of the gravity-wave one substep sizing is already built around.

**Reduced gravity.** Real values of `g` and layer depth give a gravity-wave speed near 200
m/s, which at planet-grid resolution forces the CFL-stable substep down near 100 seconds --
correct, but far more substeps than the large-scale wind pattern actually displayed needs to
look right. **Reduced gravity** (`atmosphere_cfd.REDUCED_GRAVITY_M_S2 = 0.5`) is the standard
technique barotropic models use to keep wave speed -- and thus substep count -- tractable
while the solver still *genuinely* integrates Coriolis/forcing/topography every substep.

**Substepping.** A UI "Step" (real seconds -- no relation to tectonic `years`) is subdivided
into as many CFL-stable substeps as the current grid spacing and reduced-gravity wave speed
demand (`fluid_dynamics.cfl_substeps`, CFL safety factor 0.4), computed fresh every call from
the grid's own real spacing (`fluid_dynamics.grid_spacing_m`, `cos(lat)`-corrected the same
way `climate.py` handles meridian convergence throughout) and the current flow speed. `dt`
always evenly divides the requested `seconds` (substep count computed first, `dt` derived
from it), and a hard ceiling (`MAX_SUBSTEPS_PER_STEP`, mirroring `main.MAX_ANIMATION_FRAMES`'s
own "bound worst-case request time" precedent) guards against a single request running
unboundedly long.

<a id="pole-problem"></a>
### The pole problem, and why both solvers need it solved

A fixed lat/lon grid's east-west cell spacing shrinks toward zero at the poles (`cos(lat) ->
0`) -- a real, well-known numerical challenge for exactly this kind of grid ("the pole
problem"), not specific to this codebase. Two consequences, both confirmed directly during
development:

- **CFL sizing.** If substep sizing respected the *true* spacing all the way to the pole, a
  handful of pole-adjacent cells would force an ever-smaller global timestep as grid
  resolution increases, regardless of how coarse the rest of the grid is.
- **Blow-up.** Even with substep sizing fixed, `gradient_m`/`laplacian_m`'s own `1/dx`
  (`1/dx^2` for the viscosity/diffusion terms) still explode at the pole-most rows if `dx`
  is allowed to keep shrinking -- confirmed directly: an early build stayed stable at a
  world's default (coarser) Detail setting but reliably diverged to `inf`/`NaN` within the
  first simulated day at the finest Detail setting, purely from this.

The fix (the same one real lat/lon-grid ocean/atmosphere models use, rather than a much
larger undertaking like a staggered/reduced/icosahedral grid): **`fluid_dynamics.
grid_spacing_m` clamps `dx_m` to its own value at `POLAR_FILTER_START_LAT_DEG` (75 deg) for
every row poleward of it** -- the grid's real spacing never shrinks further than that,
consistently, for every derivative *and* for CFL sizing. **`polar_zonal_filter`** blends each
row poleward of that same latitude toward its own zonal (east-west) mean, ramping to a full
zonal average at the pole-adjacent row, applied to `u`/`v`/`eta` every substep -- this is
what makes the clamp above honest: without also suppressing the *state's* own fine east-west
structure there, the raw physics would keep regenerating exactly the small-scale variation
the clamp assumes doesn't matter.

**The polar cap needs its own extra damping, too.** Removing zonal degrees of freedom near
the pole also removes the small-scale eddies that would normally dissipate momentum there, so
a merely-adequate mid-latitude drag left the zonally-averaged polar band spinning up to a
persistent, unrealistic drift under ordinary forcing (confirmed directly: sustained growth
well past any real wind-driven current, still rising after a 10-day run with nothing else on
the grid unstable). `fluid_dynamics.polar_sponge_drag_per_s` adds an extra damping *rate*
ramping up poleward of the same 75 deg latitude, on top of each solver's own ordinary drag.
Both this and ordinary drag are applied **semi-implicitly** (backward Euler for just the
linear drag term: divide by `(1 + dt*drag)` rather than subtract `dt*drag*u` outright) --
unconditionally stable regardless of how large `dt*drag` gets, which matters specifically
here: an explicit update needs `dt*drag` comfortably under 1, and confirmed directly to blow
up outright (not merely stay too energetic) once the sponge was strong enough to actually fix
the polar over-acceleration it exists to prevent.

**A related, non-polar tuning note.** The solver's own *ordinary* (non-polar) drag
(`atmosphere_cfd.SURFACE_DRAG_BASE_PER_S`) had to be tuned noticeably stronger than an initial
physically-naive guess: at high latitude generally (not
just inside the polar cap), the real Coriolis parameter is largest, which leaves a
weakly-damped system in a lightly-damped, Coriolis-dominated regime -- a persistent
near-inertial oscillation that, at too-weak a drag, hadn't decayed to a physically plausible
speed even after many simulated days. Since drag here is semi-implicit too, raising it is
never a stability risk, only a "how energetic does the flow look" tuning choice -- confirmed
directly by sweeping drag strength until sustained speeds settled into a realistic range.

<a id="ocean-cfd"></a>
### Ocean currents (retired CFD; now diagnostic in `climate.py`)

`ocean_cfd.py` was a shallow-water current solver (Coriolis + pressure-gradient momentum, a
`tau = Cd * |W| * W / MIXED_LAYER_DEPTH_M` wind stress over a near-surface mixed layer, a
temperature tracer the currents advected, and a suspended-sediment dye). On the HEALPix grid
it had no stable operating point: every drag/reduced-gravity/integration-time setting either
over-damped the interior to a near-still ~0.02 m/s or went grid-unstable at coastlines
(rogue cells, surface-height anomalies past 600 m). Fixing that properly means real
numerical-methods work (staggered C-grid, upwind advection, boundary-layer treatment).

It was removed. `climate.py`'s own `compute_ocean_currents` (Ekman base + coastal deflection +
land swirl + circumglobal boost + wake -- a stable, closed-gyre *diagnostic*, the same one that
used to be the CFD's cold-start bootstrap) runs every step now, and
`advect_ocean_temperature` carries the zonal baseline poleward along it. Current *speed* here
is a stylized, unitless quantity (wind-magnitude scale, meaningful only relative to other
cells -- direction and relative magnitude are what the map arrows and the temperature
advection use), not a physical m/s -- unchanged from how the diagnostic model always worked.

<a id="atmosphere-cfd"></a>
### Atmospheric Fluid Dynamics (`atmosphere_cfd.py`)

A shallow-water solve for wind and air temperature. Inputs: **Coriolis force** and
**elevation** (both baked into the momentum equation -- Coriolis directly, elevation via
orographic deflection around high terrain and a lapse-rate-cooled radiative-equilibrium
temperature target, see below).

**Humidity and precipitation are not solved here.** An earlier version carried humidity as a
prognostic wind-advected field with a flat evaporation source and a super-saturation
condensation sink; that produced a near-zero, uncalibrated precipitation field with no
orographic-lift term, which starved erosion and hydrology. Both fields moved back to
`climate.py`'s own diagnostic sweep (orographic lift + lake/river/vegetation moisture
recycling + per-km inland decay), fed by this solver's wind and the ocean solver's surface
temperature -- see [Climate](#climate). This solver keeps only `u`/`v`/`eta`/`temperature_c`.

**Sustained prevailing wind.** The shallow-water momentum equation here has surface drag,
viscosity, and orographic damping but no term that would sustain a planetary circulation on
its own, so the bootstrap winds (`climate.compute_wind`'s latitude bands, ~6 m/s) spun down
to a near-still thermal-gradient balance within a few tectonic steps -- taking the
Ekman-forced ocean currents down with them. `(u, v)` is now additionally relaxed, in the
fast momentum pass, toward a latitude-banded target (`WIND_FORCING_RELAXATION_PER_S` toward
the same trade-wind/westerly/polar-easterly `meridional_dir * MERIDIONAL_BASE_SPEED` +
Coriolis-deflected structure `climate.compute_wind` bootstraps from, *without* its mountain/
temperature-gradient terms -- this solver produces those itself). Steady-state planetary wind
then settles near the target's own ~6-8 m/s scale and holds there step to step, while
Coriolis, the pressure-gradient/thermal term, and orographic deflection still reshape it away
from a pure zonal band.

**Thermal forcing.** `eta` (here, a geopotential-height anomaly rather than a literal sea-
surface height) relaxes toward a target proportional to the local temperature anomaly
(`ETA_TEMPERATURE_COUPLING_M_PER_C`) -- the real atmospheric "hypsometric" relationship (a
warmer air column is thicker, so upper-level geopotential height reads higher over warm
regions), standing in for a full 3D pressure/density solve the same way `climate.py`'s own
wind model already stands in for one. Temperature itself relaxes toward a radiative-
equilibrium target computed once at mode entry from `climate.compute_land_temperature`/
`compute_ocean_temperature_baseline` (frozen elevation/insolation, so this target itself
never changes mid-session) -- diabatic heating from real insolation
(`climate.compute_insolation`, reused directly) plus lapse-rate cooling with elevation,
exactly the mechanism `climate.py`'s own land-temperature formula already uses.

**Orographic deflection.** `_mountain_deflection_tendency` cancels wind's into-slope
component near high terrain and redirects a matching amount tangentially -- the same "cancel
and redirect" *shape* `climate.py`'s own `_mountain_deflection` uses, but expressed here as a
per-second *tendency* (added into `du_dt`/`dv_dt` alongside Coriolis/pressure-gradient/drag)
rather than a direct state overwrite. That distinction is required, not stylistic: an earlier
version applied the cancel-and-redirect transform (including a >1x "Venturi speedup" factor)
directly to `(u, v)` every substep, which is stable as a single diagnostic transform
(`climate.py`'s own version only ever runs once per whole climate computation) but compounds
*geometrically* when the same transform re-applies every substep to an already-redirected
velocity -- confirmed directly as a real bug during development, wind speeds diverging to
absurd values within a single UI "Step" purely from this. The tendency form (energy-neutral:
the damped magnitude is redirected exactly, no speedup factor) only ever damps a *persistent*
into-slope flow, the same stability property ordinary drag already has.

<a id="fd-performance"></a>
### Performance and grid resolution

The wind-solver grid resolution is set by **`World.fluid_density`** (the "Fluid dynamics
resolution" Advanced-settings choice, same `climate.CLIMATE_DENSITY_CHOICES = (0.5, 1.0, 2.0,
4.0)` set `climate_density` itself uses) -- independent of `World.climate_density`, so a world
can keep a sharp climate/biome render grid while running the wind solve at a coarser (faster)
resolution, or the reverse. Defaults to `climate.DEFAULT_CLIMATE_DENSITY = 4.0`, matching `climate_density`'s
own default, so a world generated without touching this setting behaves exactly as it did
before this option existed: a "Step" can take several seconds at the finest setting (many
hundreds of CFL-stable substeps, see `fluid_dynamics.cfl_substeps`, over a ~250k-cell grid),
an accepted trade-off rather than something silently degraded -- lowering `fluid_density`
trades that away deliberately, both by shrinking the cell count itself and, since CFL substep
count scales inversely with grid spacing, by needing fewer substeps to cover the same
requested `seconds` per step. `World.node_density` (plate/elevation-line resolution) is
unrelated and has no bearing on either FD mode.

<a id="fd-render-views"></a>
### Rendering

`_render_climate_view` draws `"wind"` and the land side of `"temperature"` straight off
`World.atmosphere_cfd_state`'s own HEALPix `u`/`v`/`temperature_c`. Every other field -- the
ocean side of `"temperature"`, `"oceanCurrents"` (arrows + swell markers), `"humidity"`,
`"precipitation"`, `"biome"` -- is a `climate.compute_climate_cached` diagnostic (see
[Climate](#climate)) resampled nearest-cell onto the HEALPix cells. `render_image.render_png`
degrades to a plain background-only image if the state is somehow `None`, the same "always
renders *something* standalone" contract every other view in `VIEWS` has. The former
`"oceanCfdSediment"`/`"oceanCfdDeposition"` sediment views were removed with the ocean solver.

<a id="erosion"></a>
## Erosion (`erosion.py`)

The other half of the weather<->geology coupling: [Climate](#climate) already has terrain
influencing weather (lapse-rate cooling, mountain wind deflection, orographic rain shadow).
This module is the new direction, weather influencing terrain, implementing an erosion model
cut down to the sources that don't depend on infrastructure mantle-bloom doesn't have.

**Seven erosion sources: five subaerial, two sea-side.** Weathering's vegetation boost is
dropped (no vegetation field, same reasoning as climate.py's own "deliberately not ported"
list). Rain/sheet erosion, river-channelized erosion, weathering, glacier erosion, and seismic
erosion are the five subaerial sources; they feed into a downstream deposition pass (see
[Hydrology](#hydrology) for the flow-routing graph all of this depends on, and
[Glaciation](#glaciation) for how `glacier_depth` itself is grown/melted/flowed), so material
isn't purely one-way removed anymore: a slow, big river drops part of its sediment load locally
instead of carrying every last grain to the coast. **Submarine erosion** and **coastal
erosion** -- both mantle-bloom-original additions, previously listed as out-of-scope
"coastal-current erosion" -- are the sea floor's and the shoreline's counterparts: they erode
where the five subaerial sources are zeroed, and their debris sheds seaward onto the
surrounding sea floor as marine sediment rather than into any river's flow graph. Submarine
erosion is also what keeps a range built by two *submerged* plates colliding growing far more
slowly than a subaerial one. Glacier-driven **flattening** (broad terrain smoothing under an
ice sheet) and **seismic erosion** (earthquake-triggered landsliding) are also
mantle-bloom-original additions -- see below.

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
the same whole-world cKDTree pattern (build once, query `SLOPE_NEIGHBOR_COUNT=4`
nearest neighbors per node) `PlateWithLines.deform` uses for its own per-plate distance
queries: for each node, the elevation drop to the *lowest* of
its nearest neighbors (0 if the node is already a local minimum -- the "slope to lowest
neighbor" definition), divided by the real great-circle distance to that
neighbor. This is a genuine dimensionless rise/run -- elevation drop over real distance, not
elevation drop per grid step (a grid-step measure isn't a true slope at all, since it
silently depends on grid resolution) -- which is why `RAIN_EROSION_COEFFICIENT` was picked
by order-of-magnitude reasoning against `plates.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR`
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

**The formulas**, all computed per-node:

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
  `GLACIER_EROSION_COEFFICIENT`/`GLACIER_EROSION_MAX_FACTOR` were both raised (0.05 -> 0.09,
  2.0 -> 4.0) from an earlier, more timid starting point: basal shear stress scales with the
  ice's own *weight* (depth), and a thin valley glacier's few hundred meters versus a real
  continental ice sheet's kilometers is wide enough a range that capping the erosive
  multiplier at 2x left genuinely thick, heavy ice under-erosive relative to thin ice.
- **Seismic erosion** = `SEISMIC_EROSION_COEFFICIENT * slope * mountain_height_factor^
  SEISMIC_EROSION_ELEVATION_EXPONENT * dt_myr`, a mantle-bloom-original addition modeling
  earthquake-triggered landsliding -- a real, well-documented contributor in young,
  actively-uplifting ranges (the Himalaya, the Andes) this model previously had no source for
  at all. `mountain_height_factor = clip(elevation / SEISMIC_EROSION_ELEVATION_REFERENCE_M, 0,
  SEISMIC_EROSION_MAX_HEIGHT_FACTOR)`: this model has no explicit fault/stress field, so
  sustained elevation itself stands in for "how tectonically active/seismic this range is,"
  there being no other elevation source that reaches these heights (see
  `plates.CONVERGENT_MOUNTAIN_RATE_M_PER_MYR`). The height term is deliberately superlinear
  (exponent 2) -- real seismicity in young orogens grows faster than linearly with height/
  ongoing convergence, and this is also the model's main mechanism for capping how tall an
  actively colliding range can get before its own erosion catches up with uplift, a real
  "geomorphic ceiling" effect rather than just cosmetic variation. Still needs a slope to fail
  down (an earthquake can't landslide a perfectly flat plateau); `SEISMIC_EROSION_COEFFICIENT`
  was picked the same order-of-magnitude way `RAIN_EROSION_COEFFICIENT` was, against
  `CONVERGENT_MOUNTAIN_RATE_M_PER_MYR` (800 m/Myr). Confirmed directly at a real seed run 20
  steps (60 Myr): together with reverse-fault valleys (see [Plate motion: shift and
  deform](#boundary-evolution)), the fraction of land nodes pegged at `MAX_ELEVATION_M` dropped
  from roughly 9% to under 2% versus the same run without either addition.
- All five summed, zeroed over ocean nodes (`elevation <= 0`, the sea-level convention used
  everywhere else) -- every source above is a subaerial process. The combined result is
  capped at the node's own drop-to-lowest-neighbor (in meters, not the normalized slope), so
  a single step can't erode a node *past* the valley it drains into and carve a new, lower
  pit.
- **Submarine erosion** = `(SUBMARINE_EROSION_COEFFICIENT + SUBMARINE_PRESSURE_COEFFICIENT *
  clip(depth_below_sea_m / SUBMARINE_PRESSURE_REFERENCE_M, 0, 1)) * slope * dt_myr`, applied
  only to ocean nodes. A current-driven baseline (slope alone) plus a term that grows with how
  deep the node still sits (bottom water pressure, gravity-driven slumping of an unbuttressed
  scarp), the whole thing `* slope` so a flat abyssal plain still erodes near zero. Because the
  pressure term fades as a crest climbs toward the surface, the brake on a submerged colliding
  range is strongest while it is deep and eases continuously as it approaches sea level, handing
  off to the subaerial sources the moment it breaches. Coefficients picked the same
  order-of-magnitude way `RAIN_EROSION_COEFFICIENT` was, against
  `CONVERGENT_MOUNTAIN_RATE_M_PER_MYR` (800 m/Myr): at a moderate submarine-ridge slope
  (~0.02) and mid-depth the two terms sum to roughly half the uplift rate -- a sustained drag
  on submarine orogeny, not a hard ceiling.
- **Coastal erosion** = `(COASTAL_EROSION_WAVE_RATE_M_PER_MYR + COASTAL_FROST_MAX_RATE_M_PER_MYR
  * frost_factor) * coastal_proximity * dt_myr`, applied in the near-sea-level band on *both*
  sides of the shoreline (`coastal_proximity = clip(1 - |elevation| / COASTAL_EROSION_BAND_M,
  0, 1)`, so it tapers to zero at ±`COASTAL_EROSION_BAND_M` and also gnaws the crest of a
  mid-ocean range that rises into the band). Wave attack is a flat rate across the band (swell
  energy doesn't depend on the node's own relief -- deliberately *not* relief-gated the way
  weathering is); `frost_factor = exp(-((temperature_c - COASTAL_FROST_PEAK_C) /
  COASTAL_FROST_WIDTH_C)^2)` is a Gaussian peaked just below freezing, since freeze-thaw
  wedging needs the climate to actually cycle through 0 °C -- it falls off toward both
  permanently-frozen and never-freezing climates. Integrated every year like every other rate
  here.
- Submarine + coastal erosion are summed and capped against whatever drop-to-lowest-neighbor
  the subaerial sum didn't already claim (so the joint per-step erosion still can't carve
  below the sea floor a node drains into), then their eroded rock is handed to
  `_spread_marine_sediment`: a submerged source keeps `MARINE_SEDIMENT_LOCAL_FRACTION`
  locally, a subaerial sea-cliff source sheds it all seaward, and the remainder spreads
  inverse-distance-weighted onto the nearest *lower* ocean nodes within `MARINE_SPREAD_RANGE_KM`
  (sediment runs downhill into basins; a source with no lower ocean node in range keeps the
  full amount). Mass is conserved exactly. This marine sediment is folded into the same
  `sediment_deposited` total the `ErosionResult` carries, so [Resources and soil](#resources-and-soil)'s
  shelf oil-and-gas and soil terms see it alongside every other deposition pathway.

**Coastal planation + infill feedback** (`_coastal_openness` / `coastal_planation_amount` /
`_spread_coastal_infill`, mantle-bloom-original): every source above is either purely
subaerial or purely submarine, and none look at coastal *connectivity*, so a
marginally-submerged flat continental shelf sitting right on the waterline is a stable fixed
point that just dithers land↔ocean node-by-node forever (the per-node elevation noise
exceeds the surface's own height above/below sea level). This pass makes a coherent
coastline the stable state instead:

- **Wave exposure.** `_coastal_openness` is a per-node proxy in [0, 1] -- the fraction of
  nodes within `COASTAL_OPENNESS_RANGE_KM` (~150 km, a coarse fetch scale) that are open
  ocean (`hydrology`'s connectivity-aware `is_ocean`, so an inland-lake or enclosed-pit
  shore reads as fully sheltered). Two radius counts
  (`cKDTree.query_ball_point(..., return_length=True)`), density-independent.
- **Planation.** Land within `PLANATION_BAND_M` of sea level is ground down toward a wave-cut
  platform sitting `PLANATION_UNDERCUT_M * exposure` *below* sea level (so a genuinely
  wave-exposed low sheet is cut into open water rather than balanced on the waterline still
  dithering), at `PLANATION_RATE_M_PER_MYR * exposure * proximity * prominence * dt_myr`.
- **Infill.** The planed rock, plus `COASTAL_INFILL_MARINE_FRACTION` of that step's
  submarine + coastal erosion (redirected from `_spread_marine_sediment`'s downhill-to-deep
  spread, which is counterproductive in a shallow embayment), is spread by
  `_spread_coastal_infill` onto shallow ocean within `INFILL_DEPTH_M` of sea level, weighted
  toward the most sheltered (`shelter = 1 - openness / INFILL_SHELTER_REF`), shallowest, and
  most hollow nodes, and toward whatever headroom each still has below its fill cap (so a
  sink stops accreting as it fills -- no per-step overshoot, no iteration). A sheltered sink
  caps `INFILL_MARSH_CREST_M * shelter` *above* sea level, so a silted-up embayment emerges
  as marsh ([Biomes](#biomes)' `classify_wetland`). Exactly mass-conserving via `np.add.at`,
  same as `_spread_marine_sediment`; a source with no sheltered-shallow sink in range keeps
  its full amount locally.
- **Barrier islands** are emergent, not explicitly detected: a shallow sink with land within
  `BARRIER_LANDWARD_KM` that still faces open water (`openness >= BARRIER_MIN_OPENNESS`) gets
  a `BARRIER_PRIORITY` weight boost and a cap raised `BARRIER_CREST_M` above sea level, so a
  shore-parallel bar breaches; the water it then shelters loses open-ocean neighbours, so on
  later steps its own openness falls, its `shelter` rises, and the back-barrier lagoon silts
  up to marsh.
- **Prominence** (`PROMINENCE_REF_M` / `PROMINENCE_MAX`, from each node's height above its
  flow-graph neighbourhood mean) reweights planation up on protrusions and infill up in
  hollows -- what actually collapses a pixel-scale checkerboard, since `openness` measured at
  a ~150 km fetch scale is far too smooth to make a land/ocean call at ~60 km node spacing.
  It only reweights; the rock a prominence-boosted node loses still flows through the
  conservative infill pool.

**Glacier flattening** (`_flatten`, mantle-bloom-original): real
continental ice sheets grind down local relief over broad areas (the Canadian Shield and
Fennoscandia read as glacially smoothed bedrock today, not just eroded lower) -- a genuine
local blur, not a directional erosion/deposition term, so it's applied as a separate signed
elevation delta rather than folded into `erosion_amount`. Each node relaxes toward the mean
elevation of its own `hydrology.py` flow-graph neighbors (reusing that graph rather than a
separate query), scaled by `GLACIER_FLATTEN_RATE_PER_MYR` (0.3, also raised from 0.2 alongside
the glacier-erosion coefficients above, same "heavier ice grinds harder" reasoning) and the
same `ice_factor` glacier erosion uses -- glacier-free nodes (`ice_factor = 0`) are completely
untouched. Confirmed directly on a real run: near-zero delta at nodes with little local relief,
tens to 100+ meters at nodes combining real local relief with thick ice, consistent with
"smooths sharp terrain under ice, leaves already-flat terrain alone."

**Glacial sediment transport pushes material outside the range.** A glacier's scoured load
splits two ways: `GLACIER_TILL_FRACTION` (0.5) settles immediately as subglacial till, right
where the ice picked it up; the rest travels *with the ice itself* rather than joining the
water-routed pool -- `apply_erosion` reuses `hydrology.route_downstream` directly (the same
engine the ordinary river-deposition pool already uses, see [Hydrology](#hydrology)), but along
`hydrology.HydrologyFields.ice_flow_target` (the ice's own real downhill flow path, not water's
-- see [Glaciation](#glaciation) for why those two differ) instead of `flow_target`, retaining
in full the moment a hop reaches a node with less than `GLACIER_VISIBLE_DEPTH_M` of its own ice
-- genuinely outside the glacier -- rather than continuing to travel once there's no more ice
there to carry it. This is a real terminal moraine/outwash deposit built beyond the ice margin,
not debris stranded throughout the glacier's interior. Confirmed directly on a real run: land
nodes sitting right at a glacier's edge (ice-free themselves, with at least one glaciated
neighbor) received roughly 7x the mean sediment deposit of ordinary land elsewhere.

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

**Erosional isostatic compensation.** The net per-step surface change (all erosion minus all
deposition, plus glacial flattening and lake silt) is the rock added to or stripped from a
lithospheric column, so it is applied to `crustal_thickness_m` in full and to `elevation`
only by the resulting Airy response -- `isostatic_elevation(Hc + delta, Hm) -
isostatic_elevation(Hc, Hm)`, the same delta idiom `deform()` uses for tectonic Hc/Hm
changes (see [Isostasy](#isostasy)). An unloaded crustal root rebounds and a sediment pile
subsides under its own weight, so only ~1/6 of subaerial erosion and ~1/4 of submarine
erosion shows up as a surface elevation change; the rest is isostatic. Without this,
coastal + submarine erosion exporting continental crust to the deep ocean planed every
continent flat over a few hundred Myr once orogeny slowed (`docs/TODO.md`); with it, `elevation`
also stays a faithful readout of `isostatic_elevation(Hc, Hm)` between tectonic events rather
than drifting away from it. v1 `PlateWithLines` carries no `Hc` -- those nodes keep the bare
one-for-one response.

**Cadence: every step, no lag on climate -- but a deliberate change from erosion's own
earlier no-hydrology version regarding flow routing.** This module still calls
`climate.compute_climate(world)` fresh every step (no staleness to reason about, same as
before). Flow routing (`hydrology.compute_hydrology`) is comparatively expensive with no JIT
available, so rather than computing it twice per step, this module computes it once and
reuses the result for both erosion and `World.hydrology_cache` (see
[Hydrology](#hydrology)). Runs in `world.step_world`
right after `shift`/`deform` and topology changes, every step (line regularization now runs
inline at the end of every `deform()` call rather than on a periodic cadence -- see [Line
regularization](#line-regularization) -- and the old periodic gap-fill/reassign passes are
gone entirely, see [Whole-sphere coverage](#gap-filling)/[Boundary point
reassignment](#reassignment)).

<a id="bathymetry"></a>
## Bathymetry (`bathymetry.py`)

Submerged crust's own depth is set directly by Airy isostasy on its lithospheric column
(`lithosphere.isostatic_elevation`, see [Isostasy](#isostasy)) -- thin oceanic crust floats
low (abyssal plain), thick continental crust floats near sea level, and rifting that thins a
column subsides it. Nothing relaxes elevation toward a separate shelf/deep-water target
during the simulation any more (an earlier v1 module did; isostasy superseded it).
`SHELF_RANGE_KM`/`SHELF_RANGE_RAD` (200km) -- the shelf width `geology.py` uses to tell
shallow, oil/gas-favorable shelf water from open ocean -- lives here.

Per-plate generation seeds every column purely from its plate's crust type, which leaves
two things a real sea floor doesn't have. **`shape_initial_bathymetry(plates)`** fixes both
**once, at generation** (nothing here touches the ongoing simulation), by thinning `Hc`
(and, in the margin pass, `Hm`) -- via `lithosphere.crustal_thickness_for_submerged_elevation`
where it needs an exact target depth -- so `elevation` stays a faithful isostatic readout:

1. **`_subside_offshore_continental_crust`** -- submerged continental crust otherwise reads
   as one uniform bright shelf however far from land it sits. This drowns it toward the
   abyssal reference depth (`ABYSSAL_REFERENCE_DEPTH_M` ~= -5.2km, the same depth an aged
   oceanic column floats at) by distance to the nearest coastline: shelf within
   `OFFSHORE_SHELF_KM` (200km), full basin depth past `OFFSHORE_ABYSSAL_KM` (1400km),
   smoothstep between. Submerged *continental* nodes only, *downward* only (it never lifts
   crust the generation noise already put deep), and weighted in with depth so coastlines
   hold. Physically: hyper-extended / attenuated continental crust, exactly what real
   drowned continental interiors (Zealandia, submerged plateaus) are.

2. **`_smooth_continental_margins`** -- every continent/ocean plate boundary otherwise
   starts as a vertical cliff. This grades the `Hc`/`Hm` columns across it with an iterated,
   restricted neighbour-average (a Jacobi relaxation of the heat equation) so the seabed
   ramps shelf -> slope -> abyssal over `MARGIN_TRANSITION_KM` (400km):

   - Only *submerged* nodes within that range of opposite-type crust take part; grading
     weight ramps in with depth (`COAST_GUARD_DEPTH_M` -> `MARGIN_FULL_DEPTH_M`) and,
     smoothstepped, with proximity to the contact -- so land, deep-ocean interiors, and
     oceanic/oceanic boundaries (where a ridge or trench genuinely *is* a sharp step) are
     left untouched.
   - Weight-0 nodes are the fixed boundary values the relaxation ramps between. At the
     contact the same-crust and opposite-crust neighbour groups are weighted *equally*
     regardless of node count, so the smoothed field is continuous across the thin seam
     where the two plate lattices meet.

Later tectonics re-sharpen boundaries on their own terms (fresh ocean floor at a rift, old
floor bending into a trench), which is correct. Confirmed directly at the default density:
submerged continental crust now runs from ~-50m on the shelf to ~-5.2km in a large drowned
interior (was a near-flat ~-500m everywhere), and the median elevation step between adjacent
nodes straddling a continent/ocean boundary drops from ~4km to ~2km.

<a id="resources-and-soil"></a>
## Resources and soil (`geology.py`, plus `volcanism.py`'s own eruption roll)

Six new persistent per-node fields on `ElevationLine` (see [Why not a
grid](#why-not-a-grid)/[World state](architecture.md#world-state) for why "persistent" is
free here) -- `soil_depth`/`soil_mineral_content`/`soil_organic_content` (can rise *and*
fall, real soil erodes) and `coal_deposit_m`/`oil_gas_deposit_m`/`mineral_deposit_m`
(monotonically non-decreasing, the same self-reinforcing convention `silt_depth`/
`channel_depth` already establish -- buried peat/hydrocarbons/ore aren't un-buried by a later
climate shift). Threaded through every explicit `ElevationLine` reconstruction site the same
way `channel_depth`/`is_volcano` already are (`PlateWithLines.deform`'s own growth/shrink,
`elevation_lines.py`'s regularize interpolation, `merge_split.py`'s split) -- every other
mutation site already uses `dataclasses.replace`, which copies them automatically (see
`plates.ElevationLine`'s own docstring).

**Minerals** come from real hydrothermal circulation around volcanic activity (porphyry-copper/
VMS-style ore deposits) -- grown directly inside `volcanism.apply_volcanic_activity`'s existing
per-line eruption roll, at the same `erupts` mask that already adds `ERUPTION_ELEVATION_M`, so
no separate detection pass is needed: "an eruption deposits mineral-rich material" is exactly
what that mask already means. This only ever fires where this codebase's own volcanism model
already places `is_volcano` nodes (rift-spawned volcanic fields, see
[Volcanism](#volcanism)) -- a fair collapse of real geology's two ore-forming settings
(mid-ocean-ridge VMS deposits, arc porphyry deposits) into the one volcanism mechanism this
simulator has.

**`erosion.apply_erosion` now returns an `ErosionResult`** (`points`/`elevation`/`slope`/
`rain`/`river`/`weathering`/`sediment_deposited`/`net_elevation_change_m`/`temperature_c`/
`precipitation_mm`, all this step's already-computed per-node terms, `None` for an empty world)
instead of nothing, so
`geology.apply_resource_formation` -- called right after it (and after `bathymetry.py`/
`volcanism.py`, so `mineral_deposit_m` is already this step's fresh value before soil reads it)
from `world.step_world` -- can reuse them directly rather than re-deriving them a second time.
Two of its own private helpers became public for the same reason (mirroring
`hydrology.compute_river_speed`'s own earlier relocation): `climate_grid_indices`,
`compute_slope`. `geology.py`'s own node ordering comes directly from
`World.hydrology_cache.line_refs` -- the exact `(plate, line_index, start, end)` list that
produced `hydro`'s own flat arrays -- rather than a second independent gather, so there's
nothing to keep in sync by hand.

**Coal** forms from peat accumulating in flat, low-lying, waterlogged land -- Wetland or
Carboniferous Forest (see [Biomes](#biomes)'s `classify_wetland`, shared directly so the map's
Carboniferous Forest region always lines up with where coal is actually forming fastest).
Carboniferous Forest (warm, tropical swamp) accumulates several times faster than plain
Wetland (a cooler bog/marsh) -- real Carboniferous/Permian coal is predominantly of that
tropical-swamp origin.

**Oil and gas** form the same way in shallow marine shelf water instead of on land --
organic-rich (mostly planktonic) sediment settling in a quiet, sunlit, nutrient-fed shelf sea,
real petroleum source rocks. Gated to shelf water via the same shelf-distance concept
`bathymetry.py` already established (`bathymetry.SHELF_RANGE_RAD`, a fresh land-distance
`cKDTree` query mirroring that module's own technique), boosted near a river mouth by reusing
`hydrology.py`'s own `water_deposited` directly (real petroleum provinces cluster near deltas
-- the Gulf of Mexico, Niger Delta -- for exactly this reason; an ocean node a river empties
into already tells this module "nutrient runoff reaches the sea here," no separate river-mouth
search needed).

**Soil** forms from weathered rock (a fraction of `erosion.py`'s own weathering term becomes
fresh regolith in place) plus organic matter from a warmth/moisture-driven `productivity` term
(this codebase has no vegetation field to draw on directly), is stripped by fast rain/river
erosion (real topsoil loss), and gets an extra deposit wherever a slow, big river drops its
sediment load -- `erosion.py`'s own floodplain deposition term (`sediment_deposited`), reused
directly rather than re-derived: soil "carried down river ... and deposited in flood plains"
is literally the same mechanism erosion.py already models for sediment generally.
`soil_organic_content` relaxes toward that same `productivity` term (same
exponential-toward-target style `bathymetry.py`'s own shelf relaxation already uses);
`soil_mineral_content` relaxes toward a target driven by the node's own accumulated
`mineral_deposit_m` (weathered/hydrothermal rock feeding the soil above it) plus a small
baseline from ordinary rock weathering. Both are zeroed wherever `soil_depth` reads as zero
(bare rock has nothing to hold either) or over ocean. The richest soil needs *both* high
mineral and high organic content at once -- the Soil Quality map view scores fertility as
`sqrt(mineral_content * organic_content)`, a geometric mean that rewards having both far more
than either alone, rather than a plain average.

**Generation-time "initial soil maturity."** `geology.seed_initial_soil`, called once from
`world.generate_world` right after `plates.generate_plates` (not stored on `World` afterward --
a one-time seed, the same treatment `continental_fraction`/`land_fraction` already get), seeds
`soil_depth`/`soil_organic_content`/`soil_mineral_content` on land nodes scaled by the UI's
"initial soil maturity" slider (0 to 1, default 0). At 0, every land node starts at exactly
zero soil -- `ElevationLine`'s own zero defaults already give this, so the function is a no-op
rather than special-casing it; the planet is barren by default, the same way `channel_depth`/
`glacier_depth` start empty at generation. Deliberately *not* climate-informed at seed time (no
biome differentiation yet at generation, unlike the organic-content relaxation
`apply_resource_formation` drives every step thereafter) -- just a uniform-ish starting
maturity with a `SphereNoise` texture for visual variation, the same role `plates.py`'s own
initial elevation noise plays, on a fresh RNG stream (`seed + 2`) so it doesn't disturb
`plates.generate_plates`' own stream or `world.py`'s `seed + 1` mantle-center stream.

**Rendering.** Two new node-cloud-derived views (`render_image.py`'s `RESOURCE_VIEWS`), sharing
one fine-grid resample (`_resource_fields`, structured like `_biome_fields` but narrower --
these views need none of `climate.py`'s own fields) that -- unlike Biome/Combined -- reads only
always-defined persistent fields, so both render sensibly even before the first step (all-zero/
barren, same as any other freshly generated world's persistent fields). "Resources" is a
categorical-ish overlay: a muted, low-saturation land/ocean backdrop, then coal (land) or oil &
gas (ocean -- the two never spatially overlap, one being strictly land-only and the other
strictly ocean-only) blended in by richness fraction, then minerals blended on top last (can
co-occur with either, since volcanism isn't restricted by crust type). "Soil Quality" is a
continuous fertility heatmap (the same color-stop-interpolation technique
temperature/humidity/precipitation already use), plus the coastline overlay (a continuous
color scale carries no land/ocean cue on its own).

A third node-cloud-derived view, **"geomorph"** (a debug view -- Map View dropdown's
**Debug > Erosion & Deposition**; see [debugging.md](debugging.md)), colours every node by
`ErosionResult.net_elevation_change_m` -- this step's post-erosion elevation minus its
pre-erosion elevation, i.e. erosion minus every deposition pathway plus the small
flatten/lake-siltation terms, but *not* tectonic deform, isostasy, or volcanism. It is
retained on `World.erosion_cache` (a one-step-stale cache, same tolerance as
`climate_cache`/`hydrology_cache`, not persisted) purely so this view can nearest-node
resample it onto the same fine grid `_resource_fields` uses. The scale is diverging (warm
brown/orange where the step net-lowered a node, cool blue where it net-raised one, neutral
grey in the +-few-metre band, clamped past +-60 m/step) with the coastline overlaid. Its
purpose is to make the per-step lumpiness of near-sea-level deposition -- a +200 m spike on
one node, ~0 on its neighbour, the coastal-checkerboard mechanism -- visible, since it shows
up in no other view. A flat neutral field before the world has been stepped.

<a id="hydrology"></a>
## Hydrology: rivers and lakes (`hydrology.py`)

This module implements hydrology over mantle-bloom's irregular per-plate node cloud rather
than a fixed grid. Three core algorithms -- steepest-descent flow
direction, priority-flood basin-spill (lake/depression detection), and elevation-ordered
downstream flow accumulation -- all turn out not to actually need a *grid*, only a *graph*:
a regular 8-neighbor grid adjacency is just one convenient substrate for them. This module
builds a graph instead, via a whole-world k-nearest-neighbor query (`FLOW_NEIGHBOR_COUNT =
8`, the same technique `erosion.py` already uses for its own whole-world slope pass), then
runs the same three algorithms directly on it.

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
identity (`deform()`'s own elevation deltas, bathymetry -- via `dataclasses.replace`, not an
explicit field-by-field reconstruction, specifically so a *future* persistent field is
preserved automatically rather than needing every such call site updated by hand; see
`plates.ElevationLine`'s own docstring for the concrete bug this replaced -- `is_volcano`
silently reset to `False` every step at exactly these two sites, for several steps of actual
development, before it was caught), sliced/concatenated to match where nodes are added or
removed (`deform()`'s own growth/shrink -- new nodes start at 0, no history to carry; a
merge/split's boolean-mask slice), interpolated alongside elevation where a line gets
resampled onto a fresh spacing (`elevation_lines.regularize_line`, since that now runs at the
end of every `deform()` call and a plain reset would erase rivers/glaciers constantly, not
rarely). One call site *does* deliberately reset to 0 rather than preserve:
`plates.build_lines_from_lattice` (generation, plate merge -- genuinely new or
wholesale-resampled territory has no history to carry).

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
  river cell). A node below `FREEZE_POINT_C` (see [Glaciation](#glaciation) for why this is a
  separate, warmer threshold than the one permanent glaciers use) is unconditionally forced to
  -1 *after* the should_spill check above, rather than folded into the same check: an
  ordinary, unobstructed downhill node has `filled_elevation` exactly equal to its own
  elevation, so forcing it to look like a sink first would make it trivially satisfy
  `water_surface >= filled_elevation` and get redirected right back onto its own normal
  neighbor via `spill_target` -- silently undoing the freeze for most ordinary terrain (caught
  directly by this module's own test coverage). Applying the override to the final
  `flow_target` instead has no such escape hatch.
- **Downstream accumulation** (`route_downstream`, public -- `erosion.py` reuses it
  directly): a single forward sweep over land nodes in elevation-descending order,
  accumulating a source quantity (precipitation, for `flow_accum`; eroded material, for
  erosion's own deposition pass) along `flow_target` edges, weighted by a `retain_fraction`
  per edge (deposited locally, e.g. sediment) and, as of the moisture-recycling addition, a
  separate `loss_fraction` per edge (in-transit river evaporation, `RIVER_EVAPORATION_*`) --
  unlike a retained share, a lost share simply vanishes rather than being added to `deposited`,
  since it left as atmospheric moisture rather than staying in the channel. Correct in one
  pass because every node's target is guaranteed strictly lower, so it's always visited
  *later* in this same order.

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
`PlateWithLines.deform`/`elevation_lines.py`/`merge_split.py` exactly like `lake_depth`) is a
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
output so both can share one `cKDTree` query). Rivers are drawn as short line segments, one per drawn `is_river` node to its own
`flow_target`, in a fixed color. Which networks and how wide is `render_image._rivers_to_draw`:
group `is_river` nodes into connected drainage networks (`hydrology.group_rivers`), keep the
strongest `river_draw_max_networks(world)` by mouth `flow_accum`
(`RIVER_DRAW_MAX_NETWORKS_BY_NODE_DENSITY`, ~6-10 by `node_density`), and drop any whose
mouth can't clear the small absolute trickle floor `river_draw_min_flow(world)`
(`RIVER_DRAW_MIN_FLOW_BY_NODE_DENSITY`, also used to trim sub-threshold headwater stubs).
This is a *world-relative* cut because mouth `flow_accum` -- a physical upstream water total
-- spans orders of magnitude between an arid and a soaked world, so no single absolute floor
generalizes. Line width is a 1/2/3 px tier off each segment's `flow_accum` as a fraction of
its own network's mouth flow (`RIVER_WIDTH_TIER_FRACTIONS`), then capped by the network's
size rank (`RIVER_WIDTH_CAP_BY_RANK` = `(3,2,2)`, tail 1) so only the single largest river is
ever drawn 3 px wide. The River Inspector (below) deliberately keeps listing/drawing every
`is_river` network regardless of flow, unaffected by any of this, since picking a minor
tributary out of the full list is exactly what that view is for.

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

**Two separate cold thresholds, not one.** `FREEZE_POINT_C` (0C, the real phase-change point)
governs whether precipitation falls as snow/ice rather than rain, and whether standing lake
water or flowing river water freezes solid *this same step* -- a node below it never holds
liquid water at all: any precipitation there is treated as fully frozen (no partial
liquid/frozen split, matching `erosion.py`'s existing "use precipitation is enough"
simplification), an existing lake sitting there freezes solid into `glacier_depth`, and a
river reaching it stops flowing entirely this step. `GLACIER_ACCUMULATION_TEMP_C` (-10C) stays
a separate, colder reference used only by the melt-rate formula below (where melt bottoms out
at zero) -- this model has no seasons, so a mean annual temperature only slightly below
freezing represents a place with seasonal snow/ice that still melts back down every step, not
permanent glaciation; only the much colder zone keeps enough of its own snowfall through every
step's melt term to actually build a permanent ice sheet. Confirmed directly this stays
self-correcting: a node at, say, -3C genuinely freezes its lake/river solid every step (the
warmer threshold), but the ice it forms also melts back at a real, if reduced, rate the very
same step, converging to a near-zero *net* accumulation there rather than a spreading glacier.

Freezing (whichever threshold) happens *before* `flow_target` is (re)computed, so a lake that
just froze is correctly treated as a genuine sink again this step rather than immediately
re-filling from this same step's routed water (see `lakes.step_lakes`'s own docstring). This
ordering avoids a specific failure mode: without the freeze-before-routing ordering and an
`is_frozen` gate on `lakes.step_lakes`'s own inflow/pin-at-cap logic, a lake formed in a
warmer epoch would never freeze, and a just-frozen basin would re-flood back to its old cap
the very same step, double-counting the same water as both a lake and a glacier at once.

**A frozen river doesn't just vanish.** `_compute_flow_direction` forces a frozen node's
`flow_target` to -1 unconditionally (see [Hydrology](#hydrology) for why this override can't
reuse the ordinary should-spill escape valve a genuine closed-basin lake gets) -- water
arriving there this step has nowhere to go, so `route_downstream` deposits it in place
(`water_deposited`), and that deposited amount is folded straight into `glacier_depth` right
after routing finishes (`_update_glaciers` itself already ran by then, so this newly-frozen
river ice starts flowing/melting from next step onward -- an accepted one-step lag, the same
tolerance this codebase generally has for that kind of staleness). A river's flow also loses a
real, evaporating fraction of itself every step it's *not* frozen (`RIVER_EVAPORATION_*`,
applied as `route_downstream`'s new `loss_fraction`, distinct from `retain_fraction`'s
locally-deposited share -- an evaporated share simply vanishes rather than settling anywhere)
-- warmer nodes lose more, capped well short of 1.0 so a single very large step can't
evaporate an entire river in one hop. Excluded at an already-frozen node (already blocked
above) and at a spilling lake's own breach node (`should_spill`, see [Hydrology](#hydrology)):
that term already models a lake's concentrated overflow surge sized purely from its own
surface area, and shouldn't also lose a further fraction to evaporation in the very same step
it's cut loose. This evaporated water is genuinely lost from the river system, and is exactly
half of the "rivers, lakes, and vegetation release moisture" ask (see [Climate](#climate)'s own
"Moisture recycling" section for the other half -- that lost moisture actually reappearing as
atmospheric humidity, which for architectural reasons can't be the literal same-step amount
this module evaporates, since climate.py runs *before* this module each step).

- **Accumulation**: `GLACIER_ACCUMULATION_RATE` converts a step's frozen precipitation into
  meters of ice-depth gain, the same stylized-units-to-meters role `LAKE_FILL_RATE` plays
  for lakes. The *same* rate also converts the water a river deposits when it freezes solid
  at a cold sink (`water_deposited`, an accumulated flux -- not added 1:1 the way a standing
  frozen lake's own `lake_depth`, already a depth, is).
- **Melt**: `GLACIER_MELT_RATE_M_PER_MYR`, scaled by how far the node's temperature sits above
  `GLACIER_ACCUMULATION_TEMP_C` specifically (capped at `GLACIER_MELT_MAX_FACTOR`) -- not
  `FREEZE_POINT_C` -- melts ice back down every step regardless of whether that step's water
  is currently freezing solid or not (this is exactly what keeps the 0C-to--10C band from
  building a permanent ice sheet, see above), capped so a step can't melt more than actually
  exists. On top of that, a **depth-squared basal melt/sublimation** term
  (`GLACIER_BASAL_MELT_M_PER_MYR * (depth / GLACIER_BASAL_MELT_REFERENCE_DEPTH_M) ** 2`) runs
  *unconditionally*, even where the surface melt factor is zero: without an always-on,
  thickness-scaled sink, ice converging on a flat-floored interior sink (glacier flow scales
  with bed slope, ~0 there) accumulates without bound -- a slow cosmetic quirk when
  precipitation was near-zero, a runaway once precipitation became real. `GLACIER_MAX_DEPTH_M`
  is a hard backstop beneath that, shedding any one-step pile past ~5 km into meltwater. All
  melt -- surface, basal, and overflow -- feeds directly into that step's water source for
  `route_downstream` (real meltwater, real river discharge).
- **Flow**: a slope-scaled fraction of each node's ice (`GLACIER_FLOW_RATE_PER_MYR`, capped
  at `GLACIER_MAX_FLOW_FRACTION`) moves to its own flow target each step, via a direct
  scatter-add (`np.add.at`) rather than `route_downstream`'s elevation-ordered sweep, since
  glacier flow is a one-hop-per-step process, not a full-accumulation-to-terminus one. Ice
  reaching (or accumulating on) an ocean node is discarded rather than piling up -- real sea
  ice is a different, thinner, seasonal phenomenon this model doesn't represent; the same
  guard reads as calving where a glacier reaches a coast.

  **Ice flows on its own target, not water's.** Ice moves downhill under its own weight
  regardless of whether a given step happens to be below freezing -- a real glacier's flow is
  driven by gravity/internal ice deformation, not by the day's temperature. Liquid water's own
  `flow_target` (`_compute_flow_direction`, [Hydrology](#hydrology)) is deliberately forced
  shut at any frozen node (rivers freeze over); reusing that same target for glacier flow would
  have meant *no* glacier ever advances, since a genuinely glaciated node (cold enough to be
  accumulating ice) is almost always frozen. `HydrologyFields.ice_flow_target` is computed by
  the same steepest-descent/channel-preference/spill-redirect logic as water's own target
  (`_compute_flow_direction`'s new `apply_freeze=False` call), just without that final freeze
  override, and it's what glacier flow (and erosion.py's own glacial sediment transport, see
  [Erosion](#erosion)) actually uses. Confirmed directly this was a real, previously-silent gap
  rather than a hypothetical one: at a real seed run, ice_flow_target differed from water's own
  flow_target at roughly 10% of all nodes -- glacier flow was reaching real ground water's own
  routing never could, exactly the population of nodes this fix was meant to unlock. The
  spill-redirect test that `ice_flow_target` inherits adds the node's *own ice depth* to its
  water-surface term (`elevation + prev lake depth + prev glacier depth`), so a deep enough
  ice mass in a closed interior basin overtops its confining rim and drains toward the ocean
  via `spill_target`, the way real ice sheets feed outlet glaciers.
- **Erosion and flattening**: see [Erosion](#erosion) -- both driven by the *previous* step's
  `glacier_depth` (this step's fresh value isn't ready until this module runs, just before
  those terms are computed), the same one-step lag `channel_boost` already uses.

**Deliberately left out**: no rendering as a distinct color/
layer beyond the same `LAKE_COLOR_RGB`-style baking treatment lakes get (mantle-bloom has no
SNOW biome, so this uses its own `GLACIER_COLOR_RGB`, distinct from both `LAKE_COLOR_RGB` and
`elevation_colors`' own high-peak white/gray stops, applied the same nearest-neighbor-grid-resample way as lakes via
`plates.collect_all_glacier_depth`), no glacial eustatic sea-level coupling (glaciation is
purely local/per-node here), no seasonal accumulation/ablation cycle. Interior/convergence
ice depth is bounded three ways -- the ice-surface spill redirect (drains overfull closed
basins toward the ocean), the depth-squared basal melt (a soft equilibrium in the
low-thousands of metres), and `GLACIER_MAX_DEPTH_M` (a hard ~5 km backstop, excess shed to
meltwater) -- so a landlocked accumulation center reaches a finite equilibrium rather than
growing without bound. Over a long run an ice age still *builds* gradually as the climate
settles (glaciated node count climbs across tens of Myr), but per-node depths stay physical.

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
land `flow_accum`, plus one addition: a land node not itself under visible ice but with at
least one neighbor that is (`GLACIER_VISIBLE_DEPTH_M`), and that's actually carrying real
outflow this step (`flow_accum > 0`), is also marked `is_river` directly. Rivers commonly begin
right where a glacier's meltwater first emerges from the ice, but a fresh headwater has, by
definition, no upstream tributaries yet, so it can easily fail to clear the ordinary top-decile
threshold on its own even though it's a perfectly real river source -- this gives such a stream
a visible source right at the glacier's edge instead of only appearing once downstream
tributaries happen to push it over the percentile cut. With no grouping into distinct networks
otherwise. `hydrology.group_rivers` adds
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

<a id="lake-inspector"></a>
## Lake Inspector

A fourth interactive map mode, same "raw JSON, client renders it" philosophy as
[River Inspector](#river-inspector): `GET /world/lakes` returns every currently-visible lake's
basin info as plain JSON, and `frontend/src/LakeInspector.tsx` renders and drives the
interaction entirely client-side -- reusing the same rotation gesture and Tab/Shift+Tab
cycling pattern as the other inspectors, plus a click that always resolves to *something*
useful (the user's own spec: click a lake and see its inflows and where its water would spill
out; click ordinary dry land and see "the basin nonetheless," including its own lowest point).

**A "lake" here means one connected component of `lake_depth >
hydrology.LAKE_MIN_VISIBLE_DEPTH_M` nodes** -- exactly the mask `render_image.py`/`coastline.py`
already use to draw standing water, grouped via `hydrology.lake_components` (an undirected
union-find over the k-NN flow graph, factored out of what used to be `_lake_component_sizes`'s
own private tally so `main.py` could reuse the real grouping, not just a count). This is a
different unit than `lakes.Lake` itself: a `Lake` is a node in the depression *hierarchy* --
one basin's own catchment, dry or wet, geometric regardless of whether water currently reaches
it -- while a Lake Inspector "lake" is specifically the *currently-manifested* connected body
of standing water, which can span more than one `Lake` once several small basins have merged.
`main._smallest_lake_containing` bridges the two: given a wet component, it walks this step's
freshly-rebuilt hierarchy (`lakes.build_lake_hierarchy`, same "recomputed fresh, not persisted"
philosophy `group_rivers` already follows) top-down for the smallest tree node whose own
`members` fully contains that component -- the basin exactly as currently manifested, so its
own `outlet_node_idx`/`max_depth` is the *current* rim, not an interior saddle that's already
flooded by an earlier merge.

**A lake's displayed extent is its wet members only, not its whole geometric catchment.** A
`Lake`'s own `members` routinely includes dry higher ground on the way down to the basin's
floor (see `lakes.py`'s own docstring) -- drawing all of it as "lake" would both misrepresent
what's actually underwater and mean a click near the lake's own edge could land on a dry
catchment member and read back as a plain basin instead (confirmed directly while building
this: the very first version did exactly that on a small, newly-formed lake). `main.
_lake_basin_summary` filters to `fields.lake_depth`-wet members before reporting `member_xyz`/
`member_count`/`floor_xyz` whenever `is_lake` is true; a dry basin has no such distinction to
make, so it keeps reporting its whole catchment.

**A land click always resolves to a basin, not just a lake hit.** `GET /world/lake_at` finds
the node nearest the click, then classifies it into one of four `kind`s: `"lake"` (currently
flooded -- same info shape as a `/world/lakes` entry), `"basin"` (dry land that's still part of
a real catchment -- `main._leaf_lakes_by_node` maps every catchment-owning node straight to its
own leaf `Lake`, no hierarchy walk needed since a dry click doesn't need to align with any
currently-merged water body), `"no_basin"` (an ordinary hillslope whose steepest-descent chain
reaches the ocean without ever passing a local minimum -- `lakes.py`'s own `_OCEAN_CATCHMENT`,
never part of any basin at all), or `"ocean"`. Every non-ocean, non-no_basin case reports the
same shape: `floor_xyz`/`floor_elevation_m` (the basin's own lowest point) and
`outlet_xyz`/`outlet_elevation_m` -- "the lowest point of the edge of the basin," the saddle a
river out of it would source from, `null` for an unresolved closed/endorheic basin with no
known spill (a legitimate state, not missing data), plus every `RiverInfo` (see
[River Inspector](#river-inspector)) whose own mouth lands somewhere in the basin, regardless
of that river's own `mouth_type` label -- both `"lake"` and `"other"` read as "a river that
ends here" from the basin's own point of view.

**Rendering** is a point cloud, not traced shapes: each lake's own wet members (and, for a
selected dry basin, its whole catchment) are drawn as dots, dim for every lake and bright for
the selected basin, mirroring Plate Inspector's dim/bright split -- there's no natural polygon
outline for a scattered k-NN node group the way a plate's own boundary line already gives one
for free. The selected basin's floor and outlet (and every inflowing river's own mouth) each
draw a ring, the same drawing primitive River Inspector's mouth marker uses, so "where would
this basin drain out" is visible at a glance rather than only in the sidebar text. The same
`coastline_segments` River Inspector already fetches is reused rather than a second copy, since
`/world/lakes` computes the identical boundary -- both endpoints expose it independently
(matching `/world/rivers`' own precedent) but the frontend only needs to fetch it once.

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

- **The bounding polygon (`Plate.outline_world`/`get_bounding_polygon`) is an envelope, not
  an exact polygon -- and now a load-bearing one, not just a rendering convenience.** It's
  traced as a *staircase* from each line's current endpoints, stepping at the midpoint phi
  between adjacent rows so a straight diagonal never cuts across a concave notch between two
  rows with very different theta extents (an earlier, smoother scanline version did exactly
  that, and -- since `PlateWithLines.deform` now uses this same polygon to decide
  contested/open territory every turn, see [Plate motion: shift and
  deform](#boundary-evolution) -- confirmed directly to cause real over-claiming, not just a
  cosmetic smoothing). Since the split-row work below, the outline is traced as the exact
  boundary of the union of every row's theta-interval(s) (`_plate_outline_loops`), a general
  rectilinear-union trace that also handles a plate momentarily in two disconnected pieces,
  or with an interior hole, without the old single-staircase's diagonal over-claiming. It's
  still one `(n, 3)` vertex array -- holes and disjoint pieces are joined by zero-width
  keyhole seams (`_stitch_loops`) that the winding-number test cancels through -- so every
  `get_bounding_polygon()` consumer is unchanged. Always in sync with the real territory
  (read live from the same line data, never a separately-tracked, driftable copy). Not
  guaranteed self-intersection-safe for an arbitrarily large lateral shift between rows (e.g.
  heavy transform shearing) -- a residual source of the bounded-but-nonzero overlap noted in
  [Plate motion: shift and deform](#boundary-evolution).
- **Each `ElevationLine` is one contiguous arc; a row may carry several of them.** A line's
  two ends (`theta[0]`, `theta[-1]`) are its true territorial edges, and `deform()`'s
  ordinary grow/shrink only ever touches those ends. What it *can't* fix that way is a run of
  overridden nodes stranded in a row's interior with live nodes either side -- so `deform()`
  carves that run out (oceanic self-plate only) and hands the row back as two separate
  contiguous `ElevationLine`s at the same `phi`, with the gap between them a real hole in the
  plate's territory (see [Plate motion: shift and deform](#boundary-evolution)). `outline_
  world` / `contains_batch` / `_RowLookup` all take several lines per `phi`; `split` and
  defragmentation likewise now keep every arc of a partition-severed row (`split_into_
  contiguous_runs`) rather than dropping all but the largest. An earlier version instead left
  an interior contested patch untouched, punching no hole but never resolving the overlap
  either -- the frozen `seed 888151728` plates 9/1 case.
- **Vegetation is a derived climate classification, not a persisted field.** Climate (see
  [Climate](#climate)) feeds erosion, deposition, hydrology, and glaciation (see
  [Erosion](#erosion), [Hydrology](#hydrology), and [Glaciation](#glaciation):
  rain/river/weathering/glacier erosion, downstream deposition, rivers, lakes, and glaciers,
  with glacier flattening as a mantle-bloom-original addition) and, since the moisture-
  recycling addition (see [Climate](#climate)'s own "Moisture recycling" section), humidity's
  evapotranspiration term feeds back the other way too, sourced from `biomes.classify_biomes`.
  Vegetation still isn't a *persisted* per-node field anywhere, though (`biomes.py`'s
  classification, and climate.py's own transpiration weight derived from it, are both
  recomputed fresh from climate every call) -- so weathering's own vegetation boost
  (erosion.py) and soil formation's vegetation input (geology.py, which still falls back to a
  warmth/moisture `productivity` proxy) remain out of scope, since neither has a persisted
  field to read like hydrology.py's `lake_depth`/`channel_depth`.
- **No glacial eustatic sea-level coupling, no seasons.** Glaciation is purely local/per-node
  here -- see [Glaciation](#glaciation).
- **Single in-memory world, no persistence.** See
  [World state](architecture.md#world-state).
