# API Reference

Three routes over a single in-memory world (see
[architecture.md#world-state](architecture.md#world-state)) -- no world id, one world at a
time.

## `POST /world/generate`

Request body:

```json
{ "seed": 1, "num_plates": 12, "num_mantle_centers": 8 }
```

All fields optional (defaults: `world.DEFAULT_NUM_PLATES = 12`,
`world.DEFAULT_MANTLE_CENTERS = 8`). Replaces whatever world previously existed.

Response: a summary --

```json
{ "seed": 1, "elapsed_years": 0.0, "num_plates": 12 }
```

## `POST /world/step`

Request body:

```json
{ "years": 2000000 }
```

Advances the current world by `years` (see
[simulation-model.md](simulation-model.md) for what a step actually does). Returns the same
summary shape as `/world/generate`. `404` if no world has been generated yet.

## `GET /world/render?projection=behrmann|eckert4`

Every plate's elevation-line nodes, projected to 2D. `400` for an unrecognized projection
name, `404` if no world has been generated yet.

```json
{
  "projection": "behrmann",
  "elapsed_years": 2000000,
  "plates": [
    {
      "plate_id": 0,
      "crust_type": "continental",
      "lines": [
        { "points": [[x, y], ...], "elevation": [m, ...] }
      ]
    }
  ]
}
```

`points` and `elevation` are index-aligned and the same length. Coordinates are in the
projection's own planar units for a unit-radius sphere (see
[simulation-model.md#projections](simulation-model.md#projections)); the frontend picks its
own pixels-per-unit scale from the data's bounding box.
