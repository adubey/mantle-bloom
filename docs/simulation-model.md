# Simulation Model

This guide explains how Mantle Bloom builds and changes a planet. It is written for readers
who know basic geology and programming, but may not know the specialist terms used in
tectonics, numerical models, or this codebase. Each topic starts with a plain-language map
of the main ideas, then gives the precise rules, units, and code locations for readers who
need implementation detail.

The simulation represents terrain as cells attached to moving plates. Climate, water, ice,
and erosion act on that terrain; map views turn the model state into pictures. The diagram
shows the main links between topics, not every calculation or feedback loop.

```mermaid
flowchart LR
    F[Planet setup and plate cells] --> T[Plate motion and crust]
    T --> S[Land, sea floor, and erosion]
    S --> H[Rivers, lakes, and ice]
    H --> C[Climate and rainfall]
    C --> S
    T --> R[Maps and inspectors]
    S --> R
    H --> R
    C --> R
```

Use the topic pages below for the full explanation. The topic list is also a code map: it
points to the main Python modules that implement each part. When model behavior changes,
update the matching topic page so readers and developers have the same description of the
current model.

## Topics and implementation map

| Topic | Start here for | Main implementation areas |
|---|---|---|
| [Foundations and generation](simulation-model/foundations.md) | Sphere geometry, plate-local coordinates, initial and sketch-driven worlds, mantle flow | `geometry.py`, `plates.py` generation, `continent_layout.py`, `plate_partition.py`, `layout_lattice.py`, `synthetic_relief.py`, `worldsketch.py`, `mantle.py` |
| [Tectonics and lithosphere](simulation-model/tectonics.md) | Plate motion, boundary deformation, crust and mantle columns, faults, volcanism | `plates.py`, `quad_tectonics.py`, `lithosphere_plate.py`, `rheology.py`, `orogeny.py`, `merge_split.py`, `gaps.py`, `cratons.py`, `faults.py`, `volcanism.py` |
| [Climate and atmospheric flow](simulation-model/climate.md) | Climate fields, biomes, atmospheric solver and grid transfers | `hydroclimate/climate.py`, `fluid_dynamics.py`, `atmosphere_cfd.py`, `biomes.py` |
| [Surface processes and resources](simulation-model/surface-processes.md) | Erosion, isostasy, bathymetry, sea level, tuning and geology | `hydroclimate/erosion.py`, `bathymetry.py`, `eustasy.py`, `geology.py`, `lithosphere.py` |
| [Hydrology and glaciation](simulation-model/hydrology.md) | Flow routing, rivers, lakes, lake hierarchy, glaciers and coastlines | `hydroclimate/hydrology.py`, `lakes.py`, `ice_sheet.py`, `coastline.py` |
| [Rendering and inspectors](simulation-model/rendering-and-inspectors.md) | Map projections, image generation, view rotation, plate and fault inspection | `projections.py`, `render_image.py`, backend routes in `main.py`, frontend map and inspector components |

## Older section links

These anchors keep older links working; use the topic links above for new references.

<a id="why-not-a-grid"></a>
[Why not a grid](simulation-model/foundations.md#why-not-a-grid)
<a id="plate-local-frames"></a>
[Plate-local frames](simulation-model/foundations.md#plate-local-frames)
<a id="initial-plate-generation"></a>
[Initial plate generation (plates.py)](simulation-model/foundations.md#initial-plate-generation)
<a id="cratons-first-layout"></a>
[Cratons-first layout (continent_layout.py, plate_partition.py)](simulation-model/foundations.md#cratons-first-layout)
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
