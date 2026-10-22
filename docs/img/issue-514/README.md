# Issue #514 documentation figures

These starter SVGs are conceptual diagrams, not output from a particular seed or saved world.
They introduce the plate-local surface cells and the Hc/Hm column vocabulary. Keep them
editable, and replace or complement them with simulation renders when a documented behavior
needs evidence from a real run.

## Conceptual landscape illustrations

The JPG files are the unannotated source artwork. Each matching `*-annotated.svg` embeds
that artwork and adds editable SVG text and arrows. Keep the JPG sources checked in so the
art can be replaced without rebuilding the annotations.

- `simulation-model-overview-annotated.svg` — planet systems overview.
- `climate-landscape-annotated.svg` — wind, rainfall, rain shadow, ocean currents, and polar climate.
- `hydrology-glacier-lakes-river-annotated.svg` — glacier meltwater, connected lake basins, and river outlet.
- `globe-to-map-annotated.svg` — globe, coordinate grid, and map projection.
- `surface-processes-mountain-to-ocean-annotated.svg` — uplift, erosion, sediment transport, and crustal root.

## Figures created

- `plate-territory.svg` — oblique view of one plate's cell footprint and its independent
  neighbor lattice; used in [Plate-local frames](../../simulation-model/foundations.md#plate-local-frames).
- `plate-relief.svg` — the plate lattice carrying elevation relief and surface fields; also
  used in [Plate-local frames](../../simulation-model/foundations.md#plate-local-frames).
- `lithosphere-cross-section.svg` — conceptual convergent mountain belt with ground surface,
  Moho, crust, thickened mantle lithosphere, lithosphere base and asthenosphere; used in the
  [tectonics and lithosphere](../../simulation-model/tectonics.md) page.

## Useful next figures

1. **Seeded plate render:** show a generated world's actual plate outlines and elevation from
   the same saved world, with seed, elapsed time, projection and view stated in the caption.
   This can sit beside the conceptual pair and show what the model currently produces.
2. **Cross-section from a saved world:** capture the app's mantle cross-section along a
   stated slice, with its topography, crust/lithosphere columns and optional boundary markers.
   Record the save, slice orientation, position and time; this would replace the schematic for
   readers who want model output.
3. **Boundary change sequence:** show the same plate edge before and after a step, labeling
   newly opened cells, contested cells and cells that remain. This would make `deform()`'s
   geometric territory reconciliation easier to follow.
4. **Convergent versus divergent columns:** use paired conceptual or saved-world sections to
   compare crust and mantle-lithosphere thickening in a collision with thinning at a rift, and
   mark how the surface responds through isostasy.
5. **Inspector screenshot:** show the Plate Inspector or Plates & Faults view with one plate
   selected, plus the matching globe/map context. This would help users connect the model
   diagrams to controls in the app.

Simulation-render captions should identify the world/save, seed, elapsed time, projection or
slice, and any display settings that affect what is visible. Keep conceptual figures explicitly
marked as such.
