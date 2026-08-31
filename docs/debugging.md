# Debugging & diagnostic views

Working through a degradation (plate geometry going bad on long runs, coastlines dithering
pixel-by-pixel) usually needs a number the ordinary map views don't surface. This file
documents the debug-only views and endpoints that exist for that, and how to read them.
`docs/TODO.md` ("Diagnostic views & debug output") tracks the ones still worth building.

Debug views live in the frontend's **Map View → "Debug >"** dropdown group and, like every
other render view, are just `GET /world/render?view=...` PNGs (see
[api-reference.md](api-reference.md)'s `/world/render`). They are billed the same as any
render and carry no simulation side effects.

---

## Speckle (coastal-dither) overlay

**View:** `?view=speckle` · Map View → "Coastal dither (speckle)" · `render_image._render_speckle_view`

### What it's for

The "speckled low-relief coastlines" problem (see `docs/TODO.md`): a marginally-submerged
flat shelf whose per-node elevation noise is larger than its own height above/below sea
level, so neighbouring nodes flip land↔ocean and the coast renders as a checkerboard instead
of a shoreline. On the Elevation / Biome / Combined views that just looks like a fuzzy,
slightly-noisy coast; there was no way to see *where* the coast is a genuine checkerboard vs.
a clean line, or to make a legible before/after for a coastal-feedback change without an
ad-hoc script. This view is that script, checked in.

### The metric

For every elevation node in the raw plate node cloud (`plates.collect_all_points`), with
`sea_level = World.sea_level_m`:

- **class** = `elevation > sea_level` (land) vs. `<= sea_level` (ocean). Raw elevation only —
  no `hydrology` connectivity filter, deliberately, so this matches what the investigation
  scripts computed and so an enclosed sub-sea-level lake shore still shows up.
- **near** = `|elevation - sea_level| < SPECKLE_NEAR_BAND_M` (120 m). Only these nodes are
  drawn; everything else is just backdrop.
- **coastal-dither fraction** = of a near node's `SPECKLE_NEIGHBOR_K` (8) nearest neighbours,
  the share that are the *opposite* class from the node itself. `0.0` = the whole
  neighbourhood agrees (a coherent shoreline); higher = more disagreement.

`coastal_dither_fraction(points, elevation, sea_level_m) -> (fraction, near)` is a plain
module-level function — call it directly from a probe script against a loaded `.mbworld`.

**Reading the numbers.** The metric's natural scale is set by taking *k* nearest neighbours
on an irregular 2D node cloud:

| situation | fraction |
|---|---|
| coherent shoreline (monotonic ramp across sea level) | `< ~0.35` |
| a perfect land/ocean checkerboard | `~0.5` (the 4 orthogonal neighbours flip, the 4 diagonals don't) |
| random per-node dither | `~0.5` |
| a genuinely isolated speck (one land node ringed entirely by ocean, or vice versa) | `→ 1.0` |

So `SPECKLE_FLAG_FRACTION` (0.75) flags **isolated specks** — the single-pixel islands and
ponds — not the mixed zone. The colour ramp puts the `~0.5` checkerboard band firmly in
"hot" (orange) territory below the flag threshold.

### Reading the render

- **Backdrop:** muted olive land (`SPECKLE_LAND_BACKDROP_RGB`) / dark blue ocean
  (`SPECKLE_OCEAN_BACKDROP_RGB`), split at raw sea level.
- **Near-sea-level nodes:** a dot per node, coloured by fraction — green (clean) → yellow →
  orange (`~checkerboard`) → red (approaching isolated). `speckle_colors()` /
  `_SPECKLE_STOP_F` / `_SPECKLE_STOP_RGB`.
- **Flagged nodes** (fraction ≥ 0.75): an oversized **magenta** square (`SPECKLE_FLAG_RGB`),
  so isolated specks stand out over the ramp even at a glance.

A clean coast reads as a thin green thread one node wide. A dithering drowned shelf reads as
a broad orange/red smear with magenta flecks. Inland lake shores also light up (the metric
doesn't know they aren't ocean) — usually useful, occasionally noise.

### Doing a before/after

```python
from pathlib import Path
from app import render_image
from app.persistence import load_world_bytes    # or world.generate_world + step_world

world = load_world_bytes(Path("~/Downloads/mantle-bloom-seed888151728-85000000y.mbworld").expanduser().read_bytes())
frac, near = render_image.coastal_dither_fraction(
    *render_image.plates.collect_all_points(world.plates)[:2], world.sea_level_m
)
print(f"near={near.sum()}  flagged={(frac >= 0.75).sum()}  mean_frac={frac[near].mean():.3f}")
open("/tmp/speckle.png", "wb").write(
    render_image.render_png(world, "behrmann", "speckle", 1400, 770)
)
```

Step the world N times each way (feedback change on vs. off) and compare `flagged`,
`frac[near].mean()`, and the two PNGs. The investigation's headline metric — "315 of 768
band nodes flip land↔ocean every step" — is this same `near` band; a fix should drop
`flagged` and `mean_frac` and visibly thin the smear.

### Constants (`render_image.py`)

| constant | default | meaning |
|---|---|---|
| `SPECKLE_NEAR_BAND_M` | 120.0 | half-width of the sea-level band the overlay draws |
| `SPECKLE_NEIGHBOR_K` | 8 | nearest neighbours averaged for the fraction |
| `SPECKLE_FLAG_FRACTION` | 0.75 | fraction at/above which a node gets the magenta marker |

---

## Plate Inspector motion / shape / overlap fields

**Endpoint:** `GET /world/plates` (`main._plate_summary`) · Map View → "Plate Inspector"

Per plate, alongside the geometry: `speed_cm_per_yr` + `at_max_rate` (railed at
`mantle.MAX_PLATE_RATE`, shown red), `euler_pole` (lat/lon), `age_steps`,
`median_elevation_m` + `submerged_fraction` (red when a continental plate is >50% under
water), `overlaps` (which other plates this one's territory sits on top of, and by what
fraction of its own nodes — `main._plate_overlaps`), and `collisions`
(`world.collision_progress` timers involving the plate). See `docs/TODO.md` ("Plate geometry
degrades on long runs") for what these numbers turned up.

---

## Still worth building

See `docs/TODO.md` → "Diagnostic views & debug output" for the current list: a per-node
geomorph-rate (`sediment_deposited` / net `dElev`) diverging map, a stranded sub-sea-level
basin report, lake-churn event-log dedup, and a standalone
`python -m app.<something> <save.mbworld>` plate-diagnostics dump.
