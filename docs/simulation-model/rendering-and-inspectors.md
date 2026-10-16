# Rendering and inspectors

[Back to the Simulation Model index](../simulation-model.md).

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

**Why a grid, still.** The plate terrain is genuinely Lagrangian: nodes are spaced at
`TARGET_LINE_SPACING_KM` *within each plate's own local frame*, not on any
shared, screen-aligned grid. Drawing that raw point cloud directly -- one dot per node --
leaves visible gaps once projected: projected spacing isn't uniform (Behrmann, for
instance, stretches longitude spacing by roughly 50x at high latitudes relative to the
equator, since `x = lon * cos(30deg)` doesn't compensate for the shrinking circumference
the way the latitude term does), so a dot sized to look right at the equator leaves gaps
near the poles, and no single fixed dot size closes the gap everywhere without grossly
overlapping elsewhere. `_render_grid_arrays` fixes this the way the user originally asked:
sweep a uniform lat/lon grid over the whole sphere
(`elevation_lines.iter_local_lattice(identity_frame, spacing_rad=GRID_SPACING_RAD)` -- reusing
`iter_local_lattice` with the identity frame as a plain global lat/lon sweep, rather than any
one plate's own local frame), assign every cell its nearest
elevation node via one `cKDTree` query against every plate's current nodes (no distance
cutoff -- every cell gets *some* value, so there's no gap by construction regardless of how
far the nearest real data happens to be), and project the whole grid the same way
everything else is projected. This is a one-time resample purely for rendering -- it never
touches `world.plates`, so it has no bearing on the mass-conservation properties the rest
of the model is built around (see [Why not a grid](foundations.md#why-not-a-grid)); a fresh grid is
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
projection's own local behavior, not a fixed pixel size: `_compute_projected_climate_grid`
projects each cell's four corners (half a cell away in latitude and longitude, through the
current view rotation) and takes the largest on-screen offset in x and in y as that cell's
half-width and half-height. This is per cell, not per row, because a rotated view no longer
keeps a latitude row at one projected scale (see [Rotating the view](#rotating-the-view)).
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
2026-08-31, from the long-run plate-geometry investigation in [GitHub issue #119](https://github.com/adubey/mantle-bloom/issues/119)): the plate's
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

<a id="fault-inspector"></a>
## Plates & Faults view

An interactive map mode, same "raw JSON, client renders it" philosophy as
[River Inspector](hydrology.md#river-inspector) / [Lake Inspector](hydrology.md#lake-inspector). It merges the old
PNG "Plates" view and the standalone "Fault lines" inspector into one canvas
(`frontend/src/PlatesAndFaults.tsx`): plate outlines / bounding ellipses / node points from
`GET /world/plates` with click-to-select + Tab/Shift+Tab cycling and a metadata sidebar
(speed, Euler pole, ...), the intraplate faults + fault systems from `GET /world/faults`
drawn over them (display-only -- selecting a plate emphasises its own strands), and an
earthquake (`GET /world/earthquakes`) + volcano (`GET /world/volcanoes`) activity overlay
toggled together from the sidebar. The same `rotationDrag.ts` gesture, rotation shared with
`MapCanvas` via one lifted `App.tsx` state. The dedicated "Plate Inspector" map mode
(plates only, no faults/overlay) is still available separately.

**`fault_id` is stable across a step**, unlike `river_id` / `lake_id`: faults carry a
persistent identity for their whole lifetime (a monotonic `World.next_fault_id` counter).

**`distance_from_boundary_km`** is recomputed live per call (nearest cross-plate node via one
`cKDTree` per plate, built once per `/world/faults` call), separate from the stored
`birth_distance_from_boundary_km` -- so the sidebar can show a fault drifting away from (or
toward) a boundary as its plate moves.

**Rendering**: each trace is an ordered polyline (not a branch-capable edge list like a
river), colored by regime -- reverse magenta-red, normal yellow-green, strike-slip warm
yellow, matching `render_image.py`'s `_ELEV_REASON_RGB` fault stops (codes 15/16/17) and the
[Last elevation change](#render-image) view. Active faults draw solid; locked-up scars draw
dashed and dim. Every fault also gets a midpoint mark -- filled for active, a hollow ring for
a scar -- because a real fault is only tens of km long and its whole trace can project to a
few pixels at world scale, so the line alone wouldn't register. The selected fault draws
bright and thick with both endpoints ringed (the same `ctx.arc` primitive the River / Lake
inspectors use). Same per-edge antimeridian
longitude-unwrap as the River Inspector (`ScreenProjector.segment`). Reuses the `coastline_segments` already
fetched for the River / Lake inspectors -- this view has no other land/ocean cue.

**Fault systems** (`fault_systems` in the same `/world/faults` payload -- see
[Faults](tectonics.md#faults)) draw *under* the strand family: a broad translucent belt (the zone the
strands scatter across) plus a thin dashed centerline for the master lineament, regime-
colored. Selecting a strand highlights its parent system's belt and shows the system's
regime / length / age in the sidebar -- there's no separate system hit-test; the strands are
what you click.

**Click-to-select** is a server round trip (`GET /world/fault_at`), nearest-trace within a
120 km threshold -- there's no shape to point-in-polygon test, only sparse polylines, so a
min-distance hit-test is the natural fit (as with `river_at` / `lake_at`).
