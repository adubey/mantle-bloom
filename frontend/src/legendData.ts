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

// biomes.BIOME_NAMES / biomes.BIOME_COLORS. Kept as raw [r, g, b] tuples (not just CSS
// strings) and exported, so MapCanvas.tsx's legend-click-to-highlight feature can compare
// them directly against decoded canvas pixel data -- the Biome view draws every pixel as
// one of exactly these colors (see backend app/render_image.py's _render_biome_view: no
// coastline/graticule overlay drawn on top), so an exact RGB match is enough to pick out a
// clicked biome's cells with no server round-trip.
export const BIOME_RGB_ENTRIES: [string, [number, number, number]][] = [
  ["Ocean", [18, 28, 55]],
  ["Ice", [223, 235, 240]],
  ["Tundra", [156, 171, 158]],
  ["Boreal Forest", [61, 96, 74]],
  ["Temperate Desert", [176, 152, 116]],
  ["Temperate Grassland", [168, 178, 107]],
  ["Woodland/Shrubland", [126, 143, 90]],
  ["Temperate Seasonal Forest", [79, 121, 66]],
  ["Temperate Rainforest", [42, 94, 68]],
  ["Subtropical Desert", [214, 178, 115]],
  ["Savanna", [196, 178, 92]],
  ["Tropical Seasonal Forest", [58, 122, 66]],
  ["Tropical Rainforest", [26, 84, 46]],
  ["Wetland", [101, 111, 66]],
  ["Carboniferous Forest", [20, 66, 40]],
  ["Intertidal Zone", [70, 120, 130]],
];

const BIOME_ENTRIES: [string, string][] = BIOME_RGB_ENTRIES.map(([label, [r, g, b]]) => [label, rgb(r, g, b)]);

const COASTLINE_SYMBOL: LegendSymbol = { kind: "line", color: COASTLINE_COLOR, label: "Coastline" };

// Resources view -- see backend app/render_image.py's RESOURCE_LAND_BACKDROP_RGB/
// RESOURCE_OCEAN_BACKDROP_RGB/COAL_COLOR_RGB/OIL_GAS_COLOR_RGB/MINERAL_COLOR_RGB.
const RESOURCE_LAND_COLOR = rgb(82, 76, 68);
const RESOURCE_OCEAN_COLOR = rgb(22, 32, 48);
const COAL_COLOR = rgb(28, 24, 22);
const OIL_GAS_COLOR = rgb(92, 56, 14);
const MINERAL_COLOR = rgb(175, 62, 205);

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

// The actual fixed pixel colors Combined mode paints for these two overlays (see backend
// app/render_image.py's LAKE_COLOR_RGB/GLACIER_COLOR_RGB) -- kept separate from the
// LAKE_COLOR/GLACIER_COLOR swatch colors above since Glacier's swatch is a deliberately pure
// white for legend legibility and doesn't match the actual rendered pixel color; matching
// (see highlightTargetFor below) needs the real one.
const LAKE_RENDER_RGB: [number, number, number] = [58, 92, 122];
const GLACIER_RENDER_RGB: [number, number, number] = [221, 240, 245];

// Mirrors backend app/render_image.py's _LAND_SHADE_STOP_E/_LAND_SHADE_STOP_FACTOR -- the
// same elevation-keyed brightness multiplier Combined mode applies to a biome's flat color
// (see that module's _land_shade_factor). Combined's legend-click-to-highlight (see
// highlightTargetFor/MapCanvas.tsx) needs the actual range of shaded RGBs a biome's color can
// appear as on the rendered map, since -- unlike the plain Biome view, whose pixels are each
// exactly one of BIOME_RGB_ENTRIES with nothing else drawn on top -- Combined never paints a
// biome's flat, unshaded color for any land cell.
const LAND_SHADE_FACTORS = [0.78, 0.92, 1.05, 1.2, 1.38];

function shadedVariants([r, g, b]: [number, number, number]): [number, number, number][] {
  return LAND_SHADE_FACTORS.map((f) => [
    Math.min(255, Math.round(r * f)),
    Math.min(255, Math.round(g * f)),
    Math.min(255, Math.round(b * f)),
  ]);
}

// How close (Euclidean RGB distance) a Combined-view pixel has to be to one of a biome's
// shaded variants to still count as "that biome" for highlighting -- needs slack because
// Combined additionally blends land color toward the elevation gradient above
// RELIEF_ELEVATION_RANGE_M (see render_image.py's RELIEF_BLEND_MAX), which this deliberately
// doesn't try to reverse: a real mountain peak inside a forest biome reading as "too rocky to
// highlight" once its blend gets large is an acceptable, even fitting, edge case -- the whole
// point of that blend is to visually stop reading as flat forest color up there.
const COMBINED_MATCH_TOLERANCE = 28;

export interface PaletteEntry {
  label: string;
  colors: [number, number, number][];
}

// `palette` carries every other classifiable label's own colors alongside `selected`'s (not
// just the clicked one) so MapCanvas.tsx's applyHighlight can pick each pixel's *nearest*
// label across the whole set rather than merely testing "close enough to the selected one" in
// isolation -- see that function's own comment for why a nearest-neighbor classification is
// required to keep one biome's highlight from bleeding into another whose shaded variants
// happen to land nearby in RGB space.
export interface HighlightTarget {
  selected: string;
  palette: PaletteEntry[];
  tolerance: number;
}

// Resolves a clicked legend swatch label to a HighlightTarget (see MapCanvas.tsx's
// applyHighlight) for the two views whose legends are clickable (see Legend.tsx). Biome's own
// pixels are an exact, single fixed color per swatch (tolerance 0, and since every pixel is
// exactly one of these colors to begin with, nearest-neighbor classification never actually
// has to break a close call). Combined's land swatches can appear as any of several
// elevation-shaded variants (see shadedVariants above) within COMBINED_MATCH_TOLERANCE, while
// its Lake/Glacier swatches are still exact fixed colors, same as Biome's. Ocean/Intertidal
// Zone are excluded from Combined's palette entirely -- both are is_ocean cells that Combined
// always paints with the elevation gradient instead of any biome color (see backend
// app/render_image.py's _render_combined_view), so neither ever appears as a distinguishable
// color a click could highlight; they still appear as swatches in the legend for list parity
// with Biome (see Legend.tsx), just as non-functional ones.
export function highlightTargetFor(view: MapView, label: string): HighlightTarget | null {
  if (view === "biome") {
    if (!BIOME_RGB_ENTRIES.some(([l]) => l === label)) return null;
    const palette: PaletteEntry[] = BIOME_RGB_ENTRIES.map(([l, rgbTuple]) => ({ label: l, colors: [rgbTuple] }));
    return { selected: label, palette, tolerance: 0 };
  }
  if (view === "combined") {
    if (label === "Ocean" || label === "Intertidal Zone") return null;
    const palette: PaletteEntry[] = [
      { label: "Lake", colors: [LAKE_RENDER_RGB] },
      { label: "Glacier (ice cover)", colors: [GLACIER_RENDER_RGB] },
      ...BIOME_RGB_ENTRIES.filter(([l]) => l !== "Ocean" && l !== "Intertidal Zone").map(([l, rgbTuple]) => ({
        label: l,
        colors: shadedVariants(rgbTuple),
      })),
    ];
    if (!palette.some((p) => p.label === label)) return null;
    return { selected: label, palette, tolerance: COMBINED_MATCH_TOLERANCE };
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
        title: "Biome",
        symbols: BIOME_ENTRIES.map(([label, color]) => ({ kind: "square", color, label }) as LegendSymbol),
      };
    case "combined":
      return {
        // Ocean/land relief both follow the elevation gradient (see backend
        // app/render_image.py's _render_combined_view), and land is additionally tinted by
        // biome -- but unlike the Elevation/Plates Detail views, that gradient isn't shown as
        // its own bar here: land's per-cell brightness already visibly varies with elevation
        // (see shadedVariants/highlightTargetFor above, mirroring backend
        // _land_shade_factor), so the swatches below stay the single reference a click can
        // target, without a separate scale implying elevation is the primary thing being
        // shown. The biome list itself matches Biome's legend swatch-for-swatch (see
        // BIOME_ENTRIES) for consistency between the two views, Ocean and Intertidal Zone
        // included even though neither's own biome color is ever actually visible in Combined
        // (see highlightTargetFor's own comment) -- both always render via the elevation
        // gradient instead.
        title: "Combined (true color)",
        symbols: [
          { kind: "line", color: RIVER_COLOR, label: "River" },
          { kind: "square", color: LAKE_COLOR, label: "Lake" },
          { kind: "square", color: GLACIER_COLOR, label: "Glacier (ice cover)" },
          ...BIOME_ENTRIES.map(([label, color]) => ({ kind: "square", color, label }) as LegendSymbol),
        ],
      };
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
    default:
      return null; // plateInspector/riverInspector never had a server-drawn legend either
  }
}
