# Education applications

*Living document — exploratory ideas for using Mantle-Bloom in teaching and learning, not a committed roadmap.*

Mantle-Bloom simulates plate tectonics, mantle flow, atmospheric circulation, hydrology, erosion, glaciation, and biomes in a single coupled system (see [the simulation model overview](../docs/simulation-model.md)). Because it runs planet history forward in time and exposes internal state through map views and inspectors, it fits several educational workflows beyond passive observation.

## 1. Classroom and learning use cases

### Watching the simulation play out (time-scale intuition)
- **Slow processes made visible:** Students rarely get to see continental drift, mountain building, river capture, or ice-sheet advance happen continuously. Stepping or playing a world builds intuition for how slow tectonic rates reshape continents and climates over millions of years.
- **Demonstrations:** Works well as a live visual aid during lectures on plate boundaries, orogeny, rain shadows, or sea-level change.

### Guided mini-experiments & sensitivity labs
- **Controlled A/B experiments:** Starting from the same saved `.mbworld`, students can adjust live tuning multipliers (`world.TUNING_MULTIPLIER_FIELDS` — e.g., erosion, volcanism, sea level) and compare the forked outcome against an unchanged control world.
- **Hypothesis testing:** Encourages scientific method habits: state a prediction, change one parameter, run for N steps, and inspect how hypsometry, coastlines, river networks, or biomes diverge.

### Forensic "Mystery Planet" labs
- **Reconstructing history from evidence:** Instead of generating a fresh world, an instructor provides a mature `.mbworld` save and asks students to act as planetary geologists using the built-in inspectors (**Plate Inspector**, **Plates & Faults view**, **River Inspector**, **Lake Inspector**).
- **Example prompts:**
  - *Where did two continents collide in the past?* (look for craton sutures, orogenic relief, and fault patterns)
  - *Why is this coastal or near-coastal region arid?* (trace prevailing winds and mountain rain shadows)
  - *Which direction is this mantle plume's plate moving, and how will the coastline change next?*

### Systems thinking & coupled feedback loops
- **Cross-system cause and effect:** Most Earth-science tools isolate tectonics, weather, and hydrology. Mantle-Bloom couples them each step, letting students trace multi-hop feedback chains:
  - *Tectonics → Climate → Erosion → Isostasy:* Uplift redirects atmospheric flow and creates a rain shadow, concentrating precipitation and erosion on the windward side, which alters crustal loading and isostatic rebound.
  - *Glaciation & Eustasy:* Ice-sheet growth or retreat shifts eustatic sea level, flooding river valleys or splitting/merging lakes.

### "What-if" counterfactuals with sketch-driven worlds
- **Intentional geography:** Using sketch-driven generation (`worldsketch.py`), students can design specific continental layouts instead of relying on random seeds.
- **Counterfactual experiments:** Compare a single supercontinent against scattered island arcs, or place a major mountain barrier at the equator vs. mid-latitudes and observe how atmospheric circulation and biomes respond.

### Map projections & spatial literacy
- **Interactive distortion:** Mantle-Bloom includes multiple map projections and live view rotation. Students can rotate an evolving world and see firsthand how different projections distort polar ice sheets, plate areas, and wind belts—or test how a polar continent behaves when rotated toward the equator.

### Computational science & numerical modeling
- **Readable modeling trade-offs:** [docs/simulation-model.md](../docs/simulation-model.md) documents *why* specific numerical simplifications were chosen (e.g., plate-local frames instead of a global Eulerian grid, the atmospheric pole problem, diagnostic vs. shallow-water CFD wind solvers, and numerical guards like lake vs. sea tiers).
- **Code-to-phenomenon lab:** Useful for high school CS or undergraduate Earth-science/numerical-methods courses where students read a model rule, modify a parameter or solver choice, and observe the resulting planetary behavior.

### Human geography & worldbuilding crossovers
- **Resource & settlement reasoning:** Because geology, volcanism, hydrology, and biomes generate soil, resources, rivers, lakes, and coastlines, students can use a generated world to reason about where agriculture, trade routes, or resource conflicts would emerge—and how volcanic events or sea-level shifts would disrupt them.

---

## 2. Existing capabilities that support education

| Capability | Educational role | Where it lives |
|---|---|---|
| Live simulation playback (`Play` / `Step`) | Time-scale intuition, lecture demos | Frontend viewer + backend step loop |
| Save/load (`.mbworld`) | Distributing "mystery worlds" and A/B lab baselines | `docs/save-compatibility.md` |
| Sketch-driven generation | Testing "what-if" continental layouts | `worldsketch.py` |
| Live tuning multipliers | Controlled sensitivity experiments | `world.TUNING_MULTIPLIER_FIELDS` |
| Plate, Fault, River & Lake inspectors | Forensic inquiry and data reading | `docs/simulation-model/rendering-and-inspectors.md`, `hydrology.md` |
| Map projections & view rotation | Cartography and projection distortion labs | `projections.py` |
| Standalone desktop binary | Zero-setup classroom/student distribution (no Python/Node needed) | `docs/packaging.md` |

---

## 3. Future ideas to explore (not committed)

- **Curated `.mbworld` scenario packs:** Pre-built starter worlds (e.g., "Supercontinent breakup", "Ice-age peak", "Active subduction arc") paired with short lab worksheets.
- **Classroom-friendly "Compare" workflow:** Side-by-side or diff view for two runs from the same baseline save to make A/B tuning labs effortless.
- **Guided / Mystery mode:** Ability to lock tuning knobs or hide plate-history overlays until students submit their diagnosis.
- **Exportable lab artifacts:** One-click timelapse GIFs or summary CSV/JSON metrics (land area, mean elevation, biome breakdown) for student lab reports.

## 4. Open questions

- **Target audience:** Is the sweet spot middle/high school Earth science, undergraduate geophysics/CS, or informal science outreach?
- **Distribution:** Single-process desktop binary (`bin/package.sh`) vs. multi-user hosted classrooms (noting the single-world limit in [docs/hosting.md](../docs/hosting.md)).
