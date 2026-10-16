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
<a id="known-simplifications"></a>
[Known simplifications](simulation-model/foundations.md#known-simplifications)
<a id="hydrology"></a>
<a id="lakes-are-an-explicit-tree"></a>
[Lakes are an explicit object tree](simulation-model/hydrology.md#lakes-are-an-explicit-tree)
[Hydrology: rivers and lakes (hydrology.py)](simulation-model/hydrology.md#hydrology)
<a id="lake-vs-sea-tiers"></a>
[Lake vs. sea tiers: a depth ceiling by size (2026-09-10)](simulation-model/hydrology.md#lake-vs-sea-tiers)
<a id="glaciation"></a>
[Glaciation (hydrology.py)](simulation-model/hydrology.md#glaciation)
<a id="river-inspector"></a>
[River Inspector](simulation-model/hydrology.md#river-inspector)
<a id="lake-inspector"></a>
[Lake Inspector](simulation-model/hydrology.md#lake-inspector)
<a id="coastline"></a>
[Coastline (coastline.py)](simulation-model/hydrology.md#coastline)
<a id="projections"></a>
[Projections (projections.py)](simulation-model/rendering-and-inspectors.md#projections)
<a id="render-image"></a>
[Render image (render_image.py)](simulation-model/rendering-and-inspectors.md#render-image)
<a id="rotating-the-view"></a>
[Rotating the view](simulation-model/rendering-and-inspectors.md#rotating-the-view)
<a id="plate-inspector"></a>
[Plate Inspector](simulation-model/rendering-and-inspectors.md#plate-inspector)
<a id="fault-inspector"></a>
[Plates & Faults view](simulation-model/rendering-and-inspectors.md#fault-inspector)
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
