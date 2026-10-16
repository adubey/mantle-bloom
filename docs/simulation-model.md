# Simulation Model

This page is the map to the model documentation. Each topic page contains the full
description and links to related implementation details.

## Topics

- [Foundations and generation](simulation-model/foundations.md)
- [Tectonics and lithosphere](simulation-model/tectonics.md)
- [Climate and atmospheric flow](simulation-model/climate.md)
- [Surface processes and resources](simulation-model/surface-processes.md)
- [Hydrology and glaciation](simulation-model/hydrology.md)
- [Rendering and inspectors](simulation-model/rendering-and-inspectors.md)

## Section links

These links preserve the section anchors used by existing references.

<a id="why-not-a-grid"></a>
[Why not a grid](simulation-model/foundations.md#why-not-a-grid)
<a id="plate-local-frames"></a>
[Plate-local frames](simulation-model/foundations.md#plate-local-frames)
<a id="initial-plate-generation"></a>
[Initial plate generation (plates.py)](simulation-model/foundations.md#initial-plate-generation)
<a id="worldsketch"></a>
[Human-made worlds (sketch-driven generation) (worldsketch.py)](simulation-model/foundations.md#worldsketch)
<a id="mantle-flow"></a>
[Mantle flow (mantle.py)](simulation-model/foundations.md#mantle-flow)
<a id="boundary-evolution"></a>
[Plate motion: shift and deform (plates.py's Plate.shift/Plate.deform)](simulation-model/tectonics.md#boundary-evolution)
<a id="quad-deformation"></a>
[Quad-surface deformation (quad_tectonics.py)](simulation-model/tectonics.md#quad-deformation)
<a id="orogenic-relief"></a>
[Orogenic relief before delamination (orogeny.py)](simulation-model/tectonics.md#orogenic-relief)
<a id="merge-and-split"></a>
[Merge and split (merge_split.py)](simulation-model/tectonics.md#merge-and-split)
<a id="gap-filling"></a>
[Whole-sphere coverage: local thinning-then-melting, plus a whole-sphere fallback (gaps.py)](simulation-model/tectonics.md#gap-filling)
<a id="per-node-crust-type"></a>
[Per-node crust type](simulation-model/tectonics.md#per-node-crust-type)
<a id="cratons"></a>
[Cratons (cratons.py)](simulation-model/tectonics.md#cratons)
<a id="volcanism"></a>
[Volcanism (volcanism.py, plus lithosphere_plate.deform_columns)](simulation-model/tectonics.md#volcanism)
<a id="faults"></a>
[Faults (intraplate) (faults.py)](simulation-model/tectonics.md#faults)
<a id="projections"></a>
[Projections (projections.py)](simulation-model/rendering-and-inspectors.md#projections)
<a id="render-image"></a>
[Render image (render_image.py)](simulation-model/rendering-and-inspectors.md#render-image)
<a id="rotating-the-view"></a>
[Rotating the view](simulation-model/rendering-and-inspectors.md#rotating-the-view)
<a id="plate-inspector"></a>
[Plate Inspector](simulation-model/rendering-and-inspectors.md#plate-inspector)
<a id="climate"></a>
[Climate (climate.py)](simulation-model/climate.md#climate)
<a id="biomes"></a>
[Biomes (biomes.py)](simulation-model/climate.md#biomes)
<a id="ocean-atmospheric-fluid-dynamics"></a>
[Atmospheric wind solver (fluid_dynamics.py, atmosphere_cfd.py)](simulation-model/climate.md#ocean-atmospheric-fluid-dynamics)
<a id="mode-toggle"></a>
[Always-on, not a mode](simulation-model/climate.md#mode-toggle)
<a id="wind-model"></a>
[Wind model: "diagnostic" (default) vs "cfd"](simulation-model/climate.md#wind-model)
<a id="shallow-water-formulation"></a>
[The shallow-water formulation](simulation-model/climate.md#shallow-water-formulation)
<a id="pole-problem"></a>
[The pole problem, and why both solvers need it solved](simulation-model/climate.md#pole-problem)
<a id="ocean-cfd"></a>
[Ocean currents (retired CFD; now diagnostic in climate.py)](simulation-model/climate.md#ocean-cfd)
<a id="atmosphere-cfd"></a>
[Atmospheric Fluid Dynamics (atmosphere_cfd.py)](simulation-model/climate.md#atmosphere-cfd)
<a id="fd-performance"></a>
[Performance and grid resolution](simulation-model/climate.md#fd-performance)
<a id="fd-render-views"></a>
[Rendering](simulation-model/climate.md#fd-render-views)
<a id="erosion"></a>
[Erosion (erosion.py)](simulation-model/surface-processes.md#erosion)
<a id="bathymetry"></a>
[Bathymetry (bathymetry.py)](simulation-model/surface-processes.md#bathymetry)
<a id="eustatic-sea-level"></a>
[Eustatic sea level (eustasy.py)](simulation-model/surface-processes.md#eustatic-sea-level)
<a id="tuning-knobs"></a>
[Live tuning knobs (world.TUNING_MULTIPLIER_FIELDS)](simulation-model/surface-processes.md#tuning-knobs)
<a id="resources-and-soil"></a>
[Resources and soil (geology.py, plus volcanism.py's own eruption roll)](simulation-model/surface-processes.md#resources-and-soil)
<a id="hydrology"></a>
[Hydrology: rivers and lakes (hydrology.py)](simulation-model/hydrology.md#hydrology)
<a id="lakes-are-an-explicit-tree"></a>
[Lakes are an explicit object tree (lakes.py, not hydrology.py)](simulation-model/hydrology.md#lakes-are-an-explicit-tree)
<a id="lake-vs-sea-tiers"></a>
[Lake vs. sea tiers: a depth ceiling by size (2026-09-10)](simulation-model/hydrology.md#lake-vs-sea-tiers)
<a id="glaciation"></a>
[Glaciation (hydrology.py)](simulation-model/hydrology.md#glaciation)
<a id="river-inspector"></a>
[River Inspector](simulation-model/hydrology.md#river-inspector)
<a id="lake-inspector"></a>
[Lake Inspector](simulation-model/hydrology.md#lake-inspector)
<a id="fault-inspector"></a>
[Plates & Faults view](simulation-model/rendering-and-inspectors.md#fault-inspector)
<a id="coastline"></a>
[Coastline (coastline.py)](simulation-model/hydrology.md#coastline)
<a id="known-simplifications"></a>
[Known simplifications](simulation-model/foundations.md#known-simplifications)
<a id="initial-plate-generation-platespy"></a>
[Initial plate generation (plates.py)](simulation-model/foundations.md#initial-plate-generation-platespy)
<a id="human-made-worlds-sketch-driven-generation-worldsketchpy"></a>
[Human-made worlds (sketch-driven generation) (worldsketch.py)](simulation-model/foundations.md#human-made-worlds-sketch-driven-generation-worldsketchpy)
<a id="mantle-flow-mantlepy"></a>
[Mantle flow (mantle.py)](simulation-model/foundations.md#mantle-flow-mantlepy)
<a id="plate-motion-shift-and-deform-platespys-plateshiftplatedeform"></a>
[Plate motion: shift and deform (plates.py's Plate.shift/Plate.deform)](simulation-model/tectonics.md#plate-motion-shift-and-deform-platespys-plateshiftplatedeform)
<a id="quad-surface-deformation-quadtectonicspy"></a>
[Quad-surface deformation (quad_tectonics.py)](simulation-model/tectonics.md#quad-surface-deformation-quadtectonicspy)
<a id="orogenic-relief-before-delamination-orogenypy"></a>
[Orogenic relief before delamination (orogeny.py)](simulation-model/tectonics.md#orogenic-relief-before-delamination-orogenypy)
<a id="merge-and-split-mergesplitpy"></a>
[Merge and split (merge_split.py)](simulation-model/tectonics.md#merge-and-split-mergesplitpy)
<a id="whole-sphere-coverage-local-thinning-then-melting-plus-a-whole-sphere-fallback-gapspy"></a>
[Whole-sphere coverage: local thinning-then-melting, plus a whole-sphere fallback (gaps.py)](simulation-model/tectonics.md#whole-sphere-coverage-local-thinning-then-melting-plus-a-whole-sphere-fallback-gapspy)
<a id="cratons-cratonspy"></a>
[Cratons (cratons.py)](simulation-model/tectonics.md#cratons-cratonspy)
<a id="volcanism-volcanismpy-plus-lithosphereplatedeformcolumns"></a>
[Volcanism (volcanism.py, plus lithosphere_plate.deform_columns)](simulation-model/tectonics.md#volcanism-volcanismpy-plus-lithosphereplatedeformcolumns)
<a id="faults-intraplate-faultspy"></a>
[Faults (intraplate) (faults.py)](simulation-model/tectonics.md#faults-intraplate-faultspy)
<a id="projections-projectionspy"></a>
[Projections (projections.py)](simulation-model/rendering-and-inspectors.md#projections-projectionspy)
<a id="render-image-renderimagepy"></a>
[Render image (render_image.py)](simulation-model/rendering-and-inspectors.md#render-image-renderimagepy)
<a id="climate-climatepy"></a>
[Climate (climate.py)](simulation-model/climate.md#climate-climatepy)
<a id="biomes-biomespy"></a>
[Biomes (biomes.py)](simulation-model/climate.md#biomes-biomespy)
<a id="atmospheric-wind-solver-fluiddynamicspy-atmospherecfdpy"></a>
[Atmospheric wind solver (fluid_dynamics.py, atmosphere_cfd.py)](simulation-model/climate.md#atmospheric-wind-solver-fluiddynamicspy-atmospherecfdpy)
<a id="always-on-not-a-mode"></a>
[Always-on, not a mode](simulation-model/climate.md#always-on-not-a-mode)
<a id="wind-model-diagnostic-default-vs-cfd"></a>
[Wind model: "diagnostic" (default) vs "cfd"](simulation-model/climate.md#wind-model-diagnostic-default-vs-cfd)
<a id="the-shallow-water-formulation"></a>
[The shallow-water formulation](simulation-model/climate.md#the-shallow-water-formulation)
<a id="the-pole-problem-and-why-both-solvers-need-it-solved"></a>
[The pole problem, and why both solvers need it solved](simulation-model/climate.md#the-pole-problem-and-why-both-solvers-need-it-solved)
<a id="ocean-currents-retired-cfd-now-diagnostic-in-climatepy"></a>
[Ocean currents (retired CFD; now diagnostic in climate.py)](simulation-model/climate.md#ocean-currents-retired-cfd-now-diagnostic-in-climatepy)
<a id="atmospheric-fluid-dynamics-atmospherecfdpy"></a>
[Atmospheric Fluid Dynamics (atmosphere_cfd.py)](simulation-model/climate.md#atmospheric-fluid-dynamics-atmospherecfdpy)
<a id="performance-and-grid-resolution"></a>
[Performance and grid resolution](simulation-model/climate.md#performance-and-grid-resolution)
<a id="rendering"></a>
[Rendering](simulation-model/climate.md#rendering)
<a id="erosion-erosionpy"></a>
[Erosion (erosion.py)](simulation-model/surface-processes.md#erosion-erosionpy)
<a id="bathymetry-bathymetrypy"></a>
[Bathymetry (bathymetry.py)](simulation-model/surface-processes.md#bathymetry-bathymetrypy)
<a id="eustatic-sea-level-eustasypy"></a>
[Eustatic sea level (eustasy.py)](simulation-model/surface-processes.md#eustatic-sea-level-eustasypy)
<a id="live-tuning-knobs-worldtuningmultiplierfields"></a>
[Live tuning knobs (world.TUNING_MULTIPLIER_FIELDS)](simulation-model/surface-processes.md#live-tuning-knobs-worldtuningmultiplierfields)
<a id="effect-sizes"></a>
[Effect sizes](simulation-model/surface-processes.md#effect-sizes)
<a id="resources-and-soil-geologypy-plus-volcanismpys-own-eruption-roll"></a>
[Resources and soil (geology.py, plus volcanism.py's own eruption roll)](simulation-model/surface-processes.md#resources-and-soil-geologypy-plus-volcanismpys-own-eruption-roll)
<a id="hydrology-rivers-and-lakes-hydrologypy"></a>
[Hydrology: rivers and lakes (hydrology.py)](simulation-model/hydrology.md#hydrology-rivers-and-lakes-hydrologypy)
<a id="lake-vs-sea-tiers-a-depth-ceiling-by-size-2026-09-10"></a>
[Lake vs. sea tiers: a depth ceiling by size (2026-09-10)](simulation-model/hydrology.md#lake-vs-sea-tiers-a-depth-ceiling-by-size-2026-09-10)
<a id="glaciation-hydrologypy"></a>
[Glaciation (hydrology.py)](simulation-model/hydrology.md#glaciation-hydrologypy)
<a id="plates-faults-view"></a>
[Plates & Faults view](simulation-model/rendering-and-inspectors.md#plates-faults-view)
<a id="coastline-coastlinepy"></a>
[Coastline (coastline.py)](simulation-model/hydrology.md#coastline-coastlinepy)
