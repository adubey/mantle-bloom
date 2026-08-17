# API Reference

Three routes over a single in-memory world (see
[architecture.md#world-state](architecture.md#world-state)) -- no world id, one world at a
time.

## `POST /world/generate`

Request body:

```json
{ "seed": 1, "num_plates": null, "num_continents": 4, "num_mantle_centers": 8 }
```

All fields optional. `num_plates` defaults to `null` -- omit it (the frontend always does)
to let the world tile itself into a plausible plate count from `seed` alone
(`plates.MIN_AUTO_PLATES` to `plates.MAX_AUTO_PLATES`, see
[simulation-model.md#initial-plate-generation](simulation-model.md#initial-plate-generation)),
or pass an explicit count to override it. `num_continents` is the UI's continents slider
(`plates.MIN_CONTINENTS` to `plates.MAX_CONTINENTS`); also optional, falling back to an
independent per-plate coin flip when omitted. `num_mantle_centers` defaults to
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
an event (collisions starting/merging, consumption, splits) and how each one is worded.

## `POST /world/step`

Request body:

```json
{ "years": 2000000 }
```

Advances the current world by `years` (see
[simulation-model.md](simulation-model.md) for what a step actually does). Returns the same
summary shape as `/world/generate`, with `events` reflecting anything logged up through this
step. `404` if no world has been generated yet.

## `GET /world/render?projection=behrmann|eckert4&view=elevation|plates|platesDetail&width=1100&height=611`

Renders the current world as a PNG, base64-encoded. All drawing (elevation fill, plate-color
fill, boundary outlines, pole markers, velocity arrows, per-node dots) happens server-side
(see [simulation-model.md#render-image](simulation-model.md#render-image)) -- the client
just decodes and paints the image, it never sees raw coordinate data. `400` for an
unrecognized projection/view name or a width/height outside `[1, main.MAX_RENDER_DIMENSION_PX]`
(4000), `404` if no world has been generated yet.

```json
{
  "projection": "behrmann",
  "elapsed_years": 2000000,
  "image_base64": "iVBORw0KGgoAAAANSUhEUgAA..."
}
```

- `view` selects what gets drawn: `"elevation"` (colored by height/depth), `"plates"`
  (colored by owning plate, plus boundary outlines/pole markers/velocity arrows), or
  `"platesDetail"` (each plate's raw elevation-line nodes as dots, colored by elevation,
  plus boundary outlines) -- the frontend's Map View dropdown picks this directly.
- `width`/`height` are the returned image's exact pixel dimensions. The frontend requests
  more than its canvas's displayed CSS size (`RENDER_SCALE` in `App.tsx`) for a sharper,
  retina-style render at the same on-screen footprint; line widths, dot/pole radii, and
  padding all scale with the requested width (see `render_image.py`'s `pixel_scale`) so a
  higher-resolution request doesn't also make those features look thinner.
- `image_base64` decodes to a PNG. The frontend builds `data:image/png;base64,<this>` as an
  `<img>` source and draws it onto the canvas with `drawImage` -- see `MapCanvas.tsx`.
