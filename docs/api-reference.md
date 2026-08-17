# API Reference

Routes over a single in-memory world (see
[architecture.md#world-state](architecture.md#world-state)) -- no world id, one world at a
time.

## `POST /world/generate`

Request body:

```json
{
  "seed": 1,
  "num_plates": null,
  "continental_fraction": 0.7,
  "land_fraction": 0.29,
  "axial_tilt_deg": 23.5,
  "num_mantle_centers": 8
}
```

All fields optional. `num_plates` defaults to `null` -- omit it (the frontend always does)
to let the world tile itself into a plausible plate count from `seed` alone
(`plates.MIN_AUTO_PLATES` to `plates.MAX_AUTO_PLATES`, see
[simulation-model.md#initial-plate-generation](simulation-model.md#initial-plate-generation)),
or pass an explicit count to override it. `continental_fraction` and `land_fraction` are the
UI's two generation sliders (0 to 1, defaulting in the frontend to
`plates.DEFAULT_CONTINENTAL_FRACTION = 0.70` and `plates.DEFAULT_LAND_FRACTION = 0.29`);
both optional, each falling back to its own default behavior when omitted (an independent
per-plate coin flip for crust type, a fixed elevation floor for land) -- see
[simulation-model.md#initial-plate-generation](simulation-model.md#initial-plate-generation).
`axial_tilt_deg` is the UI's third generation slider (degrees), defaulting to
`world.DEFAULT_AXIAL_TILT_DEG = 23.5` (Earth's real tilt) -- doesn't affect plate generation,
only `climate.py`'s insolation on future renders (see
[simulation-model.md#climate](simulation-model.md#climate)), which is why it's stored on
`World` rather than consumed once here. `num_mantle_centers` defaults to
`world.DEFAULT_MANTLE_CENTERS = 8`. Replaces whatever world previously existed.

Response: a summary --

```json
{
  "seed": 1,
  "elapsed_years": 0.0,
  "num_plates": 13,
  "events": [{ "elapsed_years": 0.0, "message": "World generated with 13 plates (4 continental)." }]
}
```

`events` is the *entire* current event log (capped at `world.MAX_EVENT_LOG_LENGTH = 200`
entries, oldest dropped first), not just what changed this call -- simplest for the frontend,
which just replaces its displayed console with it on every response, and small enough not to
matter on the wire. See
[simulation-model.md#merge-and-split](simulation-model.md#merge-and-split) for what generates
an event (merges, plates disappearing, splits -- a collision merely starting is deliberately
not logged) and how each one is worded.

## `POST /world/step`

Request body:

```json
{ "years": 2000000 }
```

Advances the current world by `years` (see
[simulation-model.md](simulation-model.md) for what a step actually does). Returns the same
summary shape as `/world/generate`, with `events` reflecting anything logged up through this
step. `404` if no world has been generated yet.

## `GET /world/render?projection=behrmann|eckert4&view=elevation|plates|platesDetail|temperature|wind|oceanCurrents|humidity|precipitation&width=1100&height=611&rotation=1,0,0,0,1,0,0,0,1`

Renders the current world as a PNG, base64-encoded. All drawing (elevation fill, plate-color
fill, boundary outlines, pole markers, rotation arcs, per-node dots) happens server-side
(see [simulation-model.md#render-image](simulation-model.md#render-image)) -- the client
just decodes and paints the image, it never sees raw coordinate data. `400` for an
unrecognized projection/view name, a width/height outside `[1, main.MAX_RENDER_DIMENSION_PX]`
(4000), or a malformed `rotation`, `404` if no world has been generated yet.

```json
{
  "projection": "behrmann",
  "elapsed_years": 2000000,
  "image_base64": "iVBORw0KGgoAAAANSUhEUgAA..."
}
```

- `view` selects what gets drawn: `"elevation"` (colored by height/depth), `"plates"`
  (colored by owning plate, plus boundary outlines/pole markers/rotation arcs), or
  `"platesDetail"` (each plate's raw elevation-line nodes as dots, colored by elevation,
  plus boundary outlines) -- the frontend's Map View dropdown picks this directly. Five more
  views come from `climate.py` (see
  [simulation-model.md#climate](simulation-model.md#climate)): `"temperature"`,
  `"humidity"`, and `"precipitation"` are heatmaps; `"wind"` and `"oceanCurrents"` draw
  subsampled direction/magnitude arrows, and `"oceanCurrents"` additionally marks detected
  ocean swells with small circles.
- `rotation` is the map's current view orientation (see
  [simulation-model.md#rotating-the-view](simulation-model.md#rotating-the-view)): a
  row-major 3x3 rotation matrix as 9 comma-separated floats, applied to every real-world
  position immediately before it's projected. Optional, defaults to identity (today's
  behavior, center at lat=0/lon=0). Purely a render-time transform -- it's never stored on
  `World` and has no bearing on climate simulation results, which key off the true,
  un-rotated planetary frame regardless of what orientation is currently being viewed.
- `width`/`height` are the returned image's exact pixel dimensions. The frontend requests
  more than its canvas's displayed CSS size (`RENDER_SCALE` in `App.tsx`) for a sharper,
  retina-style render at the same on-screen footprint; line widths, dot/pole radii, and
  padding all scale with the requested width (see `render_image.py`'s `pixel_scale`) so a
  higher-resolution request doesn't also make those features look thinner.
- `image_base64` decodes to a PNG. The frontend builds `data:image/png;base64,<this>` as an
  `<img>` source and draws it onto the canvas with `drawImage` -- see `MapCanvas.tsx`.

## `GET /world/plates`

The "Plate Inspector" map mode's data source (see
[simulation-model.md#plate-inspector](simulation-model.md#plate-inspector)) -- unlike
`/world/render`, this returns plain JSON, not a baked PNG; `frontend/src/PlateInspector.tsx`
renders and drives the interaction itself. Un-rotated/true-frame throughout (no `rotation`
param -- the client applies its current view rotation only at draw time). `404` if no world
has been generated yet.

```json
{
  "elapsed_years": 2000000,
  "plates": [
    {
      "plate_id": 0,
      "crust_type": "continental",
      "num_rows": 85,
      "num_points": 3437,
      "outline": [[0.98, 0.12, -0.05], ["..."]],
      "bounding_ellipse": {
        "center_xyz": [0.97, 0.15, -0.04],
        "diameter_a_km": 13478.2,
        "diameter_b_km": 7148.9,
        "outline": [[0.91, 0.28, -0.04], ["..."]]
      }
    }
  ]
}
```

- `num_rows` is the count of `ElevationLine`s with at least one node (a plate can carry
  zero-node placeholder rows -- see `merge_split.py`'s own consumption checks -- excluded
  here to match `outline_world()`'s own filtering). `num_points` is the total node count
  across every row.
- `outline` traces the plate's live territorial boundary (`Plate.outline_world()`) as a
  closed loop of world-space unit vectors.
- `bounding_ellipse` is `null` for an empty plate (`num_points == 0`), otherwise the minimum-
  area ellipse enclosing every one of the plate's nodes (see
  [simulation-model.md#plate-inspector](simulation-model.md#plate-inspector) for the fitting
  method) -- `diameter_a_km`/`diameter_b_km` are the major/minor diameters in real km,
  `outline` is ~72 sampled world-space points around its perimeter for drawing.

## `GET /world/plate_at?lat_deg=0&lon_deg=0`

The Plate Inspector's click hit-test: which plate owns the node nearest `(lat_deg, lon_deg)`
-- both in the *true* (un-rotated) frame; the client unprojects its click through whatever
view rotation is active first (see `rotation.ts`'s `unproject`/`matTranspose`). `400` for
non-finite input, `404` if no world has been generated yet.

```json
{ "plate_id": 3 }
```

`plate_id` is `null` only if every plate is empty (shouldn't happen via the API, but a
freshly-constructed `World` with no plates at all has no node to be nearest to).
