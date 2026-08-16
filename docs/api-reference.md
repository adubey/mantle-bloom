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

Every plate's elevation-line nodes plus everything the "Plates" map view needs (pole,
velocity arrow, boundary outline), all projected to 2D. `400` for an unrecognized
projection name, `404` if no world has been generated yet.

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
      ],
      "pole": [x, y],
      "rotation_rate_deg_per_myr": 4.2,
      "velocity_arrow": { "start": [x, y], "end": [x, y] },
      "boundary": [[x, y], ...]
    }
  ]
}
```

`lines[].points`/`lines[].elevation` are index-aligned and the same length. All coordinates
(`lines[].points`, `pole`, `velocity_arrow.start`/`.end`, `boundary`) are in the projection's
own planar units for a unit-radius sphere (see
[simulation-model.md#projections](simulation-model.md#projections)); the frontend picks a
single pixels-per-unit scale from the union of all of them so switching map views never
rescales or re-centers.

- `pole` -- the plate's current Euler pole direction (`omega` normalized), or `null` if the
  plate isn't moving (`omega` is ~0).
- `velocity_arrow` -- a short great-circle arc from the plate's seed point in its current
  velocity direction, length scaled by how fast it's moving relative to
  `mantle.MAX_PLATE_RATE`; `null` under the same condition as `pole`.
- `boundary` -- the plate's rough boundary outline. Cosmetic only -- see
  [simulation-model.md#known-simplifications](simulation-model.md#known-simplifications).
