// Legend content for every /world/render map view -- rendered client-side (see Legend.tsx)
// rather than baked into the PNG server-side the way it used to be (see backend
// app/render_image.py's own "No server-side legend" comment). None of this is actually
// data-dependent (the color stops/symbols are fixed per view, not derived from world state),
// so a plain static table keyed on MapView is enough -- no new API endpoint needed. Kept in
// sync by hand with render_image.py's own color constants, the same relationship
// elevationColor.ts used to have before it was deleted in favor of server-side rendering
// (see that module's docstring) -- this is the one place that duplication is back, and only
// for the legend's own display, not for the actual map coloring.
//
// The elevation gradient's ticks are in meters *relative to the current sea level*, not an
// absolute figure -- render_image.py's elevation_colors shifts its stops by
// World.sea_level_m (see main.py's /world/controls), so "0" on this static bar always means
// "at whatever sea level currently is," accurate regardless of that live-adjustable value
// without this legend needing to know it.

import type { MapView } from "./api";

export type SwatchKind = "line" | "square" | "circle" | "ring" | "arrow" | "arc";

export interface LegendSymbol {
  kind: SwatchKind;
  color: string;
  outline?: string;
  label: string;
}

export interface GradientStop {
  value: number;
  color: string;
}

export interface LegendTick {
  value: number;
  label: string;
}

export interface LegendGradient {
  min: number;
  max: number;
  stops: GradientStop[];
  ticks: LegendTick[];
}

export interface LegendSpec {
  title: string;
  gradient?: LegendGradient;
  symbols: LegendSymbol[];
}

const rgb = (r: number, g: number, b: number) => `rgb(${r}, ${g}, ${b})`;

// --- Colors, ported from render_image.py (see that module for the authoritative values). ---

const RIVER_COLOR = rgb(77, 216, 230);
const LAKE_COLOR = rgb(58, 92, 122);
const GLACIER_COLOR = rgb(255, 255, 255);
const COASTLINE_COLOR = rgb(235, 235, 235);
const SYMBOL_COLOR = rgb(200, 205, 215);
const WIND_ARROW_COLOR = rgb(230, 230, 255);
const CURRENT_ARROW_COLOR = rgb(140, 210, 255);
const OCEAN_BACKDROP_COLOR = rgb(18, 28, 55);
const LAND_BACKDROP_COLOR = rgb(40, 46, 34);

// _ELEVATION_STOP_E / _ELEVATION_STOP_RGB.
const ELEVATION_GRADIENT: LegendGradient = {
  min: -11000,
  max: 9000,
  stops: [
    { value: -11000, color: rgb(10, 10, 40) },
    { value: -4000, color: rgb(15, 40, 110) },
    { value: -1500, color: rgb(40, 110, 190) },
    { value: -200, color: rgb(110, 170, 210) },
    { value: 0, color: rgb(200, 210, 150) },
    { value: 200, color: rgb(90, 150, 60) },
    { value: 1200, color: rgb(170, 160, 90) },
    { value: 3000, color: rgb(120, 90, 60) },
    { value: 6000, color: rgb(195, 188, 178) },
    { value: 9000, color: rgb(222, 217, 210) },
  ],
  ticks: [
    { value: -8000, label: "-8k" },
    { value: -4000, label: "-4k" },
    { value: 0, label: "0" },
    { value: 4000, label: "4k" },
    { value: 8000, label: "8k" },
  ],
};

// _TEMPERATURE_STOP_C / _TEMPERATURE_STOP_RGB.
const TEMPERATURE_GRADIENT: LegendGradient = {
  min: -60,
  max: 30,
  stops: [
    { value: -60, color: rgb(255, 255, 255) },
    { value: -45, color: rgb(148, 0, 211) },
    { value: -30, color: rgb(75, 0, 130) },
    { value: -15, color: rgb(0, 0, 255) },
    { value: 0, color: rgb(0, 200, 0) },
    { value: 10, color: rgb(255, 255, 0) },
    { value: 15, color: rgb(255, 140, 0) },
    { value: 20, color: rgb(255, 0, 0) },
    { value: 30, color: rgb(0, 0, 0) },
  ],
  ticks: [
    { value: -60, label: "-60°" },
    { value: 0, label: "0°" },
    { value: 10, label: "10°" },
    { value: 20, label: "20°" },
    { value: 30, label: "30°" },
  ],
};

// _HUMIDITY_STOP_V / _HUMIDITY_STOP_RGB.
const HUMIDITY_GRADIENT: LegendGradient = {
  min: 0,
  max: 1.4,
  stops: [
    { value: 0, color: rgb(120, 100, 60) },
    { value: 0.3, color: rgb(170, 150, 90) },
    { value: 0.6, color: rgb(140, 170, 110) },
    { value: 0.9, color: rgb(70, 150, 140) },
    { value: 1.2, color: rgb(30, 100, 140) },
    { value: 1.4, color: rgb(15, 60, 110) },
  ],
  ticks: [
    { value: 0, label: "Dry" },
    { value: 1.4, label: "Humid" },
  ],
};

// _PRECIPITATION_STOP_MM / _PRECIPITATION_STOP_RGB.
const PRECIPITATION_GRADIENT: LegendGradient = {
  min: 0,
  max: 6000,
  stops: [
    { value: 0, color: rgb(180, 160, 100) },
    { value: 500, color: rgb(200, 190, 110) },
    { value: 1200, color: rgb(140, 180, 90) },
    { value: 2400, color: rgb(60, 140, 90) },
    { value: 4000, color: rgb(40, 90, 150) },
    { value: 6000, color: rgb(20, 40, 120) },
  ],
  ticks: [
    { value: 0, label: "0" },
    { value: 3000, label: "3000" },
    { value: 6000, label: "6000" },
  ],
};

// biomes.BIOME_NAMES / biomes.BIOME_COLORS, in the same order (the 31 Köppen-Geiger land
// classes -- descriptive names, internal 3rd-level detail -- then the 10 pelagic ocean
// classes). Kept as raw [r, g, b] tuples (not just CSS strings) and exported, so
// MapCanvas.tsx's legend-click-to-highlight can compare them directly against decoded canvas
// pixel data: the Biome view draws every pixel as exactly one of these colors (see backend
// app/render_image.py's _render_biome_view -- no coastline/graticule overlay), so an exact
// RGB match picks out a clicked group's cells with no server round-trip. Hand-synced with
// biomes.py, same precedent this file always had with render_image.py's color constants.
export const BIOME_RGB_ENTRIES: [string, [number, number, number]][] = [
  ["Tropical Rainforest", [26, 51, 21]],
  ["Tropical Monsoon", [34, 66, 29]],
  ["Tropical Savanna", [74, 92, 44]],
  ["Tropical Savanna (Dry Summer)", [92, 103, 52]],
  ["Hot Desert", [176, 141, 99]],
  ["Cold Desert", [150, 136, 104]],
  ["Hot Semi-Arid", [156, 134, 78]],
  ["Cold Semi-Arid", [124, 128, 82]],
  ["Hot-Summer Mediterranean", [150, 138, 66]],
  ["Warm-Summer Mediterranean", [128, 136, 76]],
  ["Cold-Summer Mediterranean", [112, 126, 86]],
  ["Humid Subtropical (Dry Winter)", [94, 116, 54]],
  ["Subtropical Highland", [100, 122, 66]],
  ["Cold Subtropical Highland", [106, 124, 80]],
  ["Humid Subtropical", [58, 96, 44]],
  ["Oceanic", [72, 104, 56]],
  ["Subpolar Oceanic", [92, 112, 76]],
  ["Mediterranean Continental (Hot Summer)", [132, 126, 80]],
  ["Mediterranean Continental (Warm Summer)", [124, 124, 84]],
  ["Mediterranean Subarctic", [120, 122, 92]],
  ["Extremely Cold Mediterranean Subarctic", [128, 126, 104]],
  ["Monsoon Continental (Hot Summer)", [96, 110, 56]],
  ["Monsoon Continental (Warm Summer)", [104, 116, 68]],
  ["Monsoon Subarctic", [116, 118, 84]],
  ["Extremely Cold Monsoon Subarctic", [128, 122, 100]],
  ["Hot-Summer Humid Continental", [84, 104, 52]],
  ["Warm-Summer Humid Continental", [100, 114, 66]],
  ["Subarctic (Boreal)", [122, 112, 80]],
  ["Extremely Cold Subarctic", [140, 128, 100]],
  ["Tundra", [150, 152, 132]],
  ["Ice Cap", [234, 238, 240]],
  ["Tropical Open Ocean", [7, 17, 43]],
  ["Subtropical Gyre", [4, 12, 36]],
  ["Equatorial Divergence", [14, 40, 78]],
  ["Tropical Coastal Waters", [34, 82, 122]],
  ["Temperate Open Ocean", [10, 26, 58]],
  ["Temperate Shelf", [30, 70, 106]],
  ["Cold-Temperate Open Ocean", [18, 42, 72]],
  ["Cold-Temperate Shelf", [38, 74, 104]],
  ["Polar Ocean", [32, 62, 94]],
  ["Polar Sea Ice", [226, 233, 237]],
];

const RGB_BY_LABEL = new Map(BIOME_RGB_ENTRIES.map(([label, tuple]) => [label, tuple]));

// The legend shows 1st-level Köppen groups split by selected 2nd-level distinctions (and the
// pelagic classes by thermal realm), not all 41 internal classes. Each group lists the class
// labels (from BIOME_RGB_ENTRIES) it covers; a click highlights every member. `swatch` is the
// representative member whose color the legend row shows.
interface BiomeGroup {
  label: string;
  swatch: string;
  members: string[];
}

const LAND_GROUPS: BiomeGroup[] = [
  { label: "Tropical Rainforest", swatch: "Tropical Rainforest", members: ["Tropical Rainforest", "Tropical Monsoon"] },
  { label: "Tropical Savanna", swatch: "Tropical Savanna", members: ["Tropical Savanna", "Tropical Savanna (Dry Summer)"] },
  { label: "Desert", swatch: "Hot Desert", members: ["Hot Desert", "Cold Desert"] },
  { label: "Semi-Arid Steppe", swatch: "Hot Semi-Arid", members: ["Hot Semi-Arid", "Cold Semi-Arid"] },
  {
    label: "Mediterranean",
    swatch: "Hot-Summer Mediterranean",
    members: ["Hot-Summer Mediterranean", "Warm-Summer Mediterranean", "Cold-Summer Mediterranean"],
  },
  {
    label: "Humid Subtropical",
    swatch: "Humid Subtropical",
    members: ["Humid Subtropical", "Humid Subtropical (Dry Winter)"],
  },
  {
    label: "Subtropical Highland",
    swatch: "Subtropical Highland",
    members: ["Subtropical Highland", "Cold Subtropical Highland"],
  },
  { label: "Oceanic", swatch: "Oceanic", members: ["Oceanic", "Subpolar Oceanic"] },
  {
    label: "Humid Continental",
    swatch: "Warm-Summer Humid Continental",
    members: [
      "Hot-Summer Humid Continental",
      "Warm-Summer Humid Continental",
      "Mediterranean Continental (Hot Summer)",
      "Mediterranean Continental (Warm Summer)",
      "Monsoon Continental (Hot Summer)",
      "Monsoon Continental (Warm Summer)",
    ],
  },
  {
    label: "Subarctic (Boreal)",
    swatch: "Subarctic (Boreal)",
    members: [
      "Subarctic (Boreal)",
      "Extremely Cold Subarctic",
      "Mediterranean Subarctic",
      "Extremely Cold Mediterranean Subarctic",
      "Monsoon Subarctic",
      "Extremely Cold Monsoon Subarctic",
    ],
  },
  { label: "Tundra", swatch: "Tundra", members: ["Tundra"] },
  { label: "Ice Cap", swatch: "Ice Cap", members: ["Ice Cap"] },
];

// Coldest-first, so "Sea Ice" and "Polar Ocean" sit right after land's "Ice Cap" in
// BIOME_GROUPS -- the two frozen biomes (Ice Cap, Sea Ice) end up adjacent, and the whole
// ocean block stays contiguous next to them (see BIOME_GROUPS).
const OCEAN_GROUPS: BiomeGroup[] = [
  { label: "Sea Ice", swatch: "Polar Sea Ice", members: ["Polar Sea Ice"] },
  { label: "Polar Ocean", swatch: "Polar Ocean", members: ["Polar Ocean"] },
  {
    label: "Cold-Temperate Seas",
    swatch: "Cold-Temperate Open Ocean",
    members: ["Cold-Temperate Open Ocean", "Cold-Temperate Shelf"],
  },
  { label: "Temperate Seas", swatch: "Temperate Open Ocean", members: ["Temperate Open Ocean", "Temperate Shelf"] },
  {
    label: "Tropical Seas",
    swatch: "Tropical Open Ocean",
    members: ["Tropical Open Ocean", "Subtropical Gyre", "Equatorial Divergence", "Tropical Coastal Waters"],
  },
];

const BIOME_GROUPS: BiomeGroup[] = [...LAND_GROUPS, ...OCEAN_GROUPS];

const groupSwatchColor = (g: BiomeGroup): string => {
  const t = RGB_BY_LABEL.get(g.swatch)!;
  return rgb(t[0], t[1], t[2]);
};

const COASTLINE_SYMBOL: LegendSymbol = { kind: "line", color: COASTLINE_COLOR, label: "Coastline" };

// Resources view -- see backend app/render_image.py's RESOURCE_LAND_BACKDROP_RGB/
// RESOURCE_OCEAN_BACKDROP_RGB/COAL_COLOR_RGB/OIL_GAS_COLOR_RGB/MINERAL_COLOR_RGB.
const RESOURCE_LAND_COLOR = rgb(82, 76, 68);
const RESOURCE_OCEAN_COLOR = rgb(22, 32, 48);
const COAL_COLOR = rgb(28, 24, 22);
const OIL_GAS_COLOR = rgb(92, 56, 14);
const MINERAL_COLOR = rgb(175, 62, 205);

// Speckle / coastal-dither debug view -- see backend app/render_image.py's SPECKLE_* constants
// (_SPECKLE_STOP_F / _SPECKLE_STOP_RGB, SPECKLE_*_BACKDROP_RGB, SPECKLE_FLAG_RGB). The
// gradient runs over the coastal-dither fraction (share of a near-sea-level node's nearest
// neighbours on the opposite side of the waterline), 0 = clean coast, ~0.5 = checkerboard.
const SPECKLE_LAND_COLOR = rgb(54, 58, 46);
const SPECKLE_OCEAN_COLOR = rgb(24, 34, 52);
const SPECKLE_FLAG_COLOR = rgb(255, 0, 200);
const SPECKLE_GRADIENT: LegendGradient = {
  min: 0,
  max: 1,
  stops: [
    { value: 0.0, color: rgb(60, 130, 90) },
    { value: 0.2, color: rgb(150, 180, 70) },
    { value: 0.35, color: rgb(240, 205, 70) },
    { value: 0.5, color: rgb(240, 140, 50) },
    { value: 1.0, color: rgb(230, 50, 50) },
  ],
  ticks: [
    { value: 0, label: "Clean" },
    { value: 0.5, label: "Dither" },
    { value: 1, label: "Isolated" },
  ],
};

// geology.py's own soil_fertility_colors stops (_SOIL_STOP_V / _SOIL_STOP_RGB) -- fertility
// is already a [0, 1] score (sqrt(mineral_content * organic_content)), not a physical unit.
const SOIL_QUALITY_GRADIENT: LegendGradient = {
  min: 0,
  max: 1,
  stops: [
    { value: 0.0, color: rgb(168, 150, 120) },
    { value: 0.15, color: rgb(150, 120, 80) },
    { value: 0.35, color: rgb(120, 90, 55) },
    { value: 0.6, color: rgb(80, 58, 36) },
    { value: 1.0, color: rgb(35, 26, 18) },
  ],
  ticks: [
    { value: 0, label: "Barren" },
    { value: 1, label: "Rich" },
  ],
};

// render_image.py's geomorph_colors stops (_GEOMORPH_STOP_M / _GEOMORPH_STOP_RGB) -- this
// step's net per-node elevation change in metres, a diverging scale centred on 0 (warm =
// the step net-lowered a node, cool = it net-raised one), clamped past +-60 m/step.
const GEOMORPH_GRADIENT: LegendGradient = {
  min: -60,
  max: 60,
  stops: [
    { value: -60, color: rgb(120, 42, 20) },
    { value: -20, color: rgb(206, 96, 44) },
    { value: -4, color: rgb(224, 200, 170) },
    { value: 0, color: rgb(232, 232, 232) },
    { value: 4, color: rgb(168, 206, 214) },
    { value: 20, color: rgb(52, 132, 184) },
    { value: 60, color: rgb(18, 52, 112) },
  ],
  ticks: [
    { value: -60, label: "-60 m" },
    { value: 0, label: "0" },
    { value: 60, label: "+60 m" },
  ],
};

// Combined mode's per-pixel class id, carried in the render's alpha channel (see backend
// app/render_image.py's COMBINED_LAKE_ID_CODE comment): alpha = 255 - code, where code is a
// class's index in BIOME_RGB_ENTRIES + 1 (that order matches backend biomes.BIOME_NAMES, same
// hand-sync precedent as BIOME_RGB_ENTRIES itself -- every classified cell now carries one,
// land Köppen and ocean pelagic alike), or one of the two overlay codes below, or 0 for a gap
// between cells. MapCanvas.tsx reads this straight off each pixel's alpha byte -- exact, no
// RGB match -- so Combined's colors are free to span a wide shaded-relief range.
const COMBINED_LAKE_ID_CODE = BIOME_RGB_ENTRIES.length + 1;
const COMBINED_GLACIER_ID_CODE = BIOME_RGB_ENTRIES.length + 2;

function combinedIdCode(label: string): number | null {
  const idx = BIOME_RGB_ENTRIES.findIndex(([l]) => l === label);
  return idx < 0 ? null : idx + 1;
}

function groupFor(label: string): BiomeGroup | undefined {
  return BIOME_GROUPS.find((g) => g.label === label);
}

export interface PaletteEntry {
  label: string;
  colors: [number, number, number][];
}

// Two ways to pick out a clicked swatch's pixels (see MapCanvas.tsx's applyHighlight):
//   - `palette` + `tolerance`: nearest-neighbor RGB match across every classifiable label
//     (Biome view -- every pixel is exactly one of BIOME_RGB_ENTRIES, tolerance 0).
//   - `idCodes`: exact match on the per-pixel id carried in the alpha channel (Combined view
//     -- see combinedIdCode above). When set, `palette`/`tolerance` are unused.
export interface HighlightTarget {
  selected: string;
  palette: PaletteEntry[];
  tolerance: number;
  idCodes?: number[];
}

// Resolves a clicked legend group label to a HighlightTarget for the two clickable views (see
// Legend.tsx). A legend row is a 1st/2nd-level group covering several internal classes, so a
// click highlights every member: Biome matches on the members' exact fixed pixel colors,
// Combined on the members' alpha-channel id codes. The "Sea Ice" group also folds in the
// glacier overlay code (glaciated ocean paints over Polar-Sea-Ice cells with nearly the same
// colour).
export function highlightTargetFor(view: MapView, label: string): HighlightTarget | null {
  if (view === "combined" && label === "Lake") {
    return { selected: label, palette: [], tolerance: 0, idCodes: [COMBINED_LAKE_ID_CODE] };
  }
  if (view === "combined" && label === "Glacier (ice cover)") {
    return { selected: label, palette: [], tolerance: 0, idCodes: [COMBINED_GLACIER_ID_CODE] };
  }

  const group = groupFor(label);
  if (!group) return null;

  if (view === "biome") {
    // Nearest-neighbour match across every group, so one group's highlight can't bleed into a
    // neighbour whose colour lands nearby -- `selected` picks out this group's members.
    const palette: PaletteEntry[] = BIOME_GROUPS.map((g) => ({
      label: g.label,
      colors: g.members.map((m) => RGB_BY_LABEL.get(m)!),
    }));
    return { selected: label, palette, tolerance: 0 };
  }
  if (view === "combined") {
    const idCodes = group.members
      .map((m) => combinedIdCode(m))
      .filter((c): c is number => c !== null);
    if (label === "Sea Ice") idCodes.push(COMBINED_GLACIER_ID_CODE);
    if (idCodes.length === 0) return null;
    return { selected: label, palette: [], tolerance: 0, idCodes };
  }
  return null;
}

export function legendFor(view: MapView): LegendSpec | null {
  switch (view) {
    case "elevation":
      return {
        title: "Elevation (m)",
        gradient: ELEVATION_GRADIENT,
        symbols: [
          { kind: "line", color: RIVER_COLOR, label: "River" },
          { kind: "square", color: LAKE_COLOR, label: "Lake" },
          { kind: "square", color: GLACIER_COLOR, label: "Glacier (ice cover)" },
        ],
      };
    case "plates":
      return {
        title: "Plates",
        symbols: [
          { kind: "line", color: SYMBOL_COLOR, label: "Plate boundary" },
          { kind: "circle", color: SYMBOL_COLOR, outline: "#ffffff", label: "Euler pole" },
          { kind: "arc", color: SYMBOL_COLOR, label: "Rotation (rate & direction)" },
        ],
      };
    case "platesDetail":
      return {
        title: "Elevation (m)",
        gradient: ELEVATION_GRADIENT,
        symbols: [{ kind: "line", color: SYMBOL_COLOR, label: "Plate boundary" }],
      };
    case "temperature":
      return { title: "Temperature (°C)", gradient: TEMPERATURE_GRADIENT, symbols: [COASTLINE_SYMBOL] };
    case "humidity":
      return { title: "Humidity", gradient: HUMIDITY_GRADIENT, symbols: [COASTLINE_SYMBOL] };
    case "precipitation":
      return { title: "Precipitation (mm/yr)", gradient: PRECIPITATION_GRADIENT, symbols: [COASTLINE_SYMBOL] };
    case "wind":
      return {
        title: "Wind",
        symbols: [
          { kind: "arrow", color: WIND_ARROW_COLOR, label: "Speed (arrow length)" },
          { kind: "square", color: OCEAN_BACKDROP_COLOR, label: "Ocean" },
          { kind: "square", color: LAND_BACKDROP_COLOR, label: "Land" },
        ],
      };
    case "oceanCurrents":
      return {
        title: "Ocean currents",
        symbols: [
          { kind: "arrow", color: CURRENT_ARROW_COLOR, label: "Speed (arrow length)" },
          { kind: "square", color: OCEAN_BACKDROP_COLOR, label: "Ocean" },
          { kind: "square", color: LAND_BACKDROP_COLOR, label: "Land" },
          { kind: "ring", color: "#ffffff", label: "Ocean swell" },
        ],
      };
    case "biome":
      return {
        // 1st-level Köppen groups split by selected 2nd-level distinctions, then the pelagic
        // ocean classes by thermal realm (see BIOME_GROUPS). Each row covers several internal
        // 3rd-level classes; a click highlights every member (see highlightTargetFor).
        title: "Köppen Climate",
        symbols: BIOME_GROUPS.map((g) => ({ kind: "square", color: groupSwatchColor(g), label: g.label }) as LegendSymbol),
      };
    case "combined": {
      // Same grouped Köppen/pelagic list as the Biome view (see BIOME_GROUPS) for
      // consistency between the two, plus the river/lake/ice overlays drawn on top. The
      // hypsometric elevation gradient isn't shown as its own bar: land brightness already
      // varies visibly with elevation and ocean with depth (see _render_combined_view), so
      // the grouped swatches stay the single click target.
      //
      // The three overlays sit where their subject matter does rather than all bunched at
      // the front: "Glacier (ice cover)" between land's "Ice Cap" and ocean's "Sea Ice" with
      // the other frozen swatches, and "Lake"/"River" right after the ocean block.
      const groupSymbol = (g: BiomeGroup): LegendSymbol => ({
        kind: "square",
        color: groupSwatchColor(g),
        label: g.label,
      });
      return {
        title: "Elevation & Köppen Climate",
        symbols: [
          ...LAND_GROUPS.map(groupSymbol),
          { kind: "square", color: GLACIER_COLOR, label: "Glacier (ice cover)" },
          ...OCEAN_GROUPS.map(groupSymbol),
          { kind: "square", color: LAKE_COLOR, label: "Lake" },
          { kind: "line", color: RIVER_COLOR, label: "River" },
        ],
      };
    }
    case "resources":
      return {
        title: "Resources",
        symbols: [
          { kind: "square", color: RESOURCE_LAND_COLOR, label: "Land" },
          { kind: "square", color: RESOURCE_OCEAN_COLOR, label: "Ocean" },
          { kind: "square", color: COAL_COLOR, label: "Coal" },
          { kind: "square", color: OIL_GAS_COLOR, label: "Oil & Gas" },
          { kind: "square", color: MINERAL_COLOR, label: "Minerals" },
        ],
      };
    case "soilQuality":
      return { title: "Soil Quality", gradient: SOIL_QUALITY_GRADIENT, symbols: [COASTLINE_SYMBOL] };
    case "speckle":
      return {
        title: "Coastal dither (opposite-class neighbour fraction)",
        gradient: SPECKLE_GRADIENT,
        symbols: [
          { kind: "square", color: SPECKLE_FLAG_COLOR, label: "Flagged (fraction ≥ 0.75)" },
          { kind: "square", color: SPECKLE_LAND_COLOR, label: "Land backdrop" },
          { kind: "square", color: SPECKLE_OCEAN_COLOR, label: "Ocean backdrop" },
        ],
      };
    case "geomorph":
      return { title: "Elevation change / step", gradient: GEOMORPH_GRADIENT, symbols: [COASTLINE_SYMBOL] };
    default:
      return null; // plateInspector/riverInspector never had a server-drawn legend either
  }
}
