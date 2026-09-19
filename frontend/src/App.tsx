import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { KeyboardEvent } from "react";
import "./index.css";
import {
  animateWorld, fetchCornerNotchLog, fetchDebugScenarios, fetchEarthquakes, fetchElevationPoint, fetchElevationPointAt, fetchFaults, fetchLakes, fetchNodeAt, fetchPlates, fetchPointSample, fetchRivers, fetchStats, fetchStatsHistory, fetchVolcanoes, fetchWorldSummary, generateDebugWorld, generateWorld, renderWorld, stepWorld, stopAnimation, updateControls,
  TUNING_MULTIPLIER_KEYS,
} from "./api";
import type {
  AnimateResponse, CornerNotchLogEntry, DebugScenario, EarthquakeSummary, ElevationPointResponse, FaultSummary, FaultSystemSummary, LakeAtResponse, LakeSummary, MapView, NodeAtResponse, PlateSummary, PointSample, Projection, RenderResponse, RiverSummary, Segment, TuningKey, TuningMultipliers, VolcanoSummary, WorldStats, WorldSummary,
} from "./api";
import MapCanvas from "./MapCanvas";
import SketchEditor from "./SketchEditor";
import PlateInspector from "./PlateInspector";
import RiverInspector from "./RiverInspector";
import LakeInspector from "./LakeInspector";
import PlatesAndFaults from "./PlatesAndFaults";
import EventConsole from "./EventConsole";
import CornerNotchLogPanel from "./CornerNotchLogPanel";
import StatsModal from "./StatsModal";
import ControlsModal from "./ControlsModal";
import AdvancedSettingsModal from "./AdvancedSettingsModal";
import FileModal from "./FileModal";
import AnimationModal from "./AnimationModal";
import SaveAnimationModal from "./SaveAnimationModal";
import Legend from "./Legend";
import MeasureOverlay from "./MeasureOverlay";
import ProgressBar from "./ProgressBar";
import { PREMADE_WORLDS } from "./premadeWorlds";
import { faultKindForLegendLabel, highlightTargetFor } from "./legendData";
import { centerOfRotation, IDENTITY_ROTATION, rotationForCenter } from "./rotation";
import type { Mat3 } from "./rotation";
import { getCookie, setCookie } from "./cookies";

// The map is displayed at DISPLAY_WIDTH x DISPLAY_HEIGHT (CSS pixels, unchanged from
// before), but the image requested from the server is RENDER_SCALE times bigger -- a
// sharper, retina-style render at the same on-screen size, rather than a bigger map.
// Whether the page was loaded with ?deb in the URL -- gates a cluster of developer-only UI:
// the Generate World dialog's "Debugging Worlds" tab (see backend debug_worlds.py), the Event
// Console and corner-notch log side panels, and the "Debug >" group of Map View choices.
// Without it, none of that is present at all (not just disabled), so an ordinarily-generated
// world's sidebar/dialogs stay uncluttered by default. Read once at module load (this is a
// browser-only SPA bundle, so `window` is always available) since it's a one-time flag for the
// whole page load, not something that changes while the app is open. Works identically on the
// packaged prod server and the dev/debug server -- it only ever looks at the URL the browser
// actually loaded, never which backend is serving it.
const DEBUG_UI = new URLSearchParams(window.location.search).has("deb");

const DISPLAY_WIDTH = 1100;
const DISPLAY_HEIGHT = 611;
const RENDER_SCALE = 2;
const RENDER_WIDTH = DISPLAY_WIDTH * RENDER_SCALE;
const RENDER_HEIGHT = DISPLAY_HEIGHT * RENDER_SCALE;
const STEP_YEARS_OPTIONS = [10_000, 100_000, 1_000_000, 10_000_000];
const PLAY_INTERVAL_MS = 400;
// Percent, matching backend app/plates.py's DEFAULT_CONTINENTAL_FRACTION/DEFAULT_LAND_FRACTION.
const DEFAULT_CONTINENTAL_PERCENT = 70;
const DEFAULT_LAND_PERCENT = 29;
// Degrees, matching backend app/world.py's DEFAULT_AXIAL_TILT_DEG (Earth's real tilt).
const DEFAULT_AXIAL_TILT_DEG = 23.5;
// Percent, matching backend app/geology.py's own default (0, a fully barren starting world --
// see geology.seed_initial_soil).
const DEFAULT_INITIAL_SOIL_MATURITY_PERCENT = 0;
// The Generate dialog's single "Detail" control, driving both backend app/plates.py's
// NODE_DENSITY_CHOICES/DEFAULT_NODE_DENSITY (elevation point density) and app/climate.py's
// CLIMATE_DENSITY_CHOICES/DEFAULT_CLIMATE_DENSITY (climate & biome resolution) together --
// the two share the same discrete multiplier set, so one dial covers both rather than
// asking the user to reason about two separately. A discrete set, not a free-form slider,
// since there's no continuous unit for "how many points/cells," only "how many times as
// many." "Very Low" is the coarsest, fastest option; the default ("Standard") is one step
// short of the finest. "High" is the exception to the "one dial, one multiplier" rule above:
// it trades more cost for sharper detail unevenly -- 1.5x Standard's point density but 2x its
// climate-grid resolution -- so it carries its own `climateDensity` override (see
// climateDensityForDetail below) instead of reusing `value` as both.
const DETAIL_CHOICES: { value: number; label: string; climateDensity?: number }[] = [
  { value: 6, label: "High", climateDensity: 8 },
  { value: 4, label: "Standard" },
  { value: 2, label: "Medium" },
  { value: 1, label: "Low" },
  { value: 0.5, label: "Very Low" },
];
const DEFAULT_DETAIL = 4;
// node_density and climate_density (see DETAIL_CHOICES' own comment) agree for every choice
// except "High", which asks for more climate resolution than point density -- look up its
// override, falling back to `detail` itself (node_density == climate_density) everywhere else.
function climateDensityForDetail(detail: number): number {
  return DETAIL_CHOICES.find((d) => d.value === detail)?.climateDensity ?? detail;
}
// The Advanced-settings dialog's own "Fluid dynamics resolution" choice (see backend app/
// world.py's World.fluid_density) -- same shape as DETAIL_CHOICES but capped at its own top
// choice, "High" (value 2, unrelated to -- and unchanged by -- the "High" now atop
// DETAIL_CHOICES above; the two dials' labels aren't on a shared scale): the atmospheric wind
// solver runs every step (see docs/simulation-model.md#ocean-atmospheric-fluid-dynamics), so
// there's no "only pay for Standard (or finer) when you opt in" case left to justify offering
// this dial anything past that, matching backend app/climate.py's own FLUID_DENSITY_CHOICES.
const FLUID_DETAIL_CHOICES: { value: number; label: string }[] = [
  { value: 2, label: "High" },
  { value: 1, label: "Medium" },
  { value: 0.5, label: "Low" },
];
const DEFAULT_FLUID_DETAIL = 2;
// MIN_PLATES matches backend app/plates.py's MIN_AUTO_PLATES; the manual slider deliberately
// runs past the auto (seed-based) MAX_AUTO_PLATES (20) up to 40, so a hand-picked value can be
// denser than "Auto" would ever choose. The backend puts no upper clamp on an explicit
// num_plates; ids past 20 just wrap the 20-colour PLATE_PALETTE (cosmetic only).
const MIN_PLATES = 8;
const MAX_PLATES = 40;
const DEFAULT_PLATES = 14;
// The Advanced-settings "Voronoi points" slider -- the total number of Voronoi seed points the
// plate tiling scatters before merging cells down to the chosen plate count (see backend
// lithosphere_plate.generate_plates' voronoi_points param). Higher makes plate outlines
// lumpier/less convex (and, for a sketch-driven world, makes the continental/oceanic boundary
// hug the drawn coastline more tightly -- see DEFAULT_VORONOI_POINTS_SKETCH below); lower makes
// them smoother/coarser. Max was 10,000 (down to 2,000, see issue #128): VoronoiPreview.tsx's
// client-side preview does a brute-force O(240x120xpoints) scan that got noticeably slow near
// the old max.
const MIN_VORONOI_POINTS = 8;
const MAX_VORONOI_POINTS = 2000;
// Default for "random" (procedural, no sketch) worlds.
const DEFAULT_VORONOI_POINTS_RANDOM = 500;
// Default for sketch-driven worlds ("Human-made" and "Premade worlds", including Pangaea) --
// the new max. A sketch's land/sea comes from the drawing itself, but only for nodes whose
// plate is already continental/oceanic to begin with (see lithosphere_plate.generate_plates);
// that continental/oceanic boundary is a Voronoi cell boundary between the sketch-placed sites,
// not the sketch's own outline, so more points (smaller cells) makes generated coastlines
// resemble the actual drawing much more closely instead of clipping/filling past a coarse,
// blobby plate boundary.
const DEFAULT_VORONOI_POINTS_SKETCH = MAX_VORONOI_POINTS;
// Matching backend app/world.py's World.sea_level_m/World.solar_multiplier defaults.
const DEFAULT_SEA_LEVEL_M = 0;
const DEFAULT_SOLAR_MULTIPLIER = 1;
// Matching backend app/world.py's World.ice_age_period_years default -- 0 == no ice-age cycle.
const DEFAULT_ICE_AGE_PERIOD_YEARS = 0;
// Every geomorphic-budget tuning knob defaults to 1.0 (untuned) -- matches every
// *_multiplier default on backend app/world.py's World.
const DEFAULT_TUNING: TuningMultipliers = Object.fromEntries(
  TUNING_MULTIPLIER_KEYS.map((k) => [k, 1]),
) as TuningMultipliers;
// Matching backend app/world.py's World.simulate_plate_movement/World.simulate_climate_biomes
// defaults -- both on, i.e. a normal full simulation.
const DEFAULT_SIMULATE_PLATE_MOVEMENT = true;
const DEFAULT_SIMULATE_CLIMATE_BIOMES = true;
// Matching backend app/world.py's World.wind_model default -- the fast closed-form diagnostic
// (ABL) wind: it reproduces ~85-90% of the CFD biome map for a fraction of the per-step cost,
// so it's the better starting point; switch to "cfd" in Controls for the full shallow-water solve.
const DEFAULT_WIND_MODEL = "diagnostic";
// Matching backend app/world.py's World.fault_deformation_mode default. "fault" localises
// plate-boundary deformation onto fault lines (faults spawn boundary-hugging, so the
// collision zone still deforms -- as fault-tracking ridges rather than one smooth swell);
// "boundary" is the pre-faults-rework behaviour (smooth uplift/rift bands at the polygon
// edge); "both" runs the boundary bands plus the scaled-up fault relief. See faults.py /
// LithospherePlate.deform.
const DEFAULT_FAULT_DEFORMATION_MODE = "fault";
// Off by default for an ordinarily-generated/loaded world -- see World.debug_diagnostics.
const DEFAULT_DEBUG_DIAGNOSTICS = false;

function randomSeed(): number {
  return Math.floor(Math.random() * 1_000_000_000);
}

function formatLatLon(latDeg: number, lonDeg: number): string {
  const latDir = latDeg >= 0 ? "N" : "S";
  const lonDir = lonDeg >= 0 ? "E" : "W";
  return `${Math.abs(latDeg).toFixed(1)}°${latDir}, ${Math.abs(lonDeg).toFixed(1)}°${lonDir}`;
}

function isIdentityRotation(rotation: Mat3): boolean {
  return rotation.every((v, i) => v === IDENTITY_ROTATION[i]);
}

// The map view to show after switching to `mode` -- keeps the current view if it's already
// valid there (e.g. switching Tectonics & Climate -> Ocean Fluid Dynamics -> back doesn't
// force "elevation" back on if the user had "temperature" selected), and falls back to that
// Persists the map's view state (projection/mapView/rotation) across a browser refresh --
// these three are otherwise pure client-local React state (see the `rotation` field's own
// comment above), so without this a refresh would silently reset the view to its defaults
// even when /world/summary below finds the same world still sitting in server memory. Kept
// as its own small cookie rather than folded into anything server-side since it's display
// state, not simulation state -- same reasoning `rotation` itself already gets.
const VIEW_COOKIE_NAME = "mantle-bloom-view";
const MAP_VIEW_CHOICES = new Set<MapView>([

  "elevation", "platesDetail", "speckle", "temperature", "wind", "oceanCurrents", "humidity", "precipitation", "biome", "combined",
  "resources", "soilQuality", "geomorph", "elevReason", "overlapAge", "nodeAge", "plateInspector", "riverInspector", "lakeInspector", "platesAndFaults",
]);
const PROJECTION_CHOICES = new Set<Projection>(["behrmann", "eckert4"]);
// The Map View select's "Debug >" optgroup (see the DEBUG_UI comment above) -- kept as its own
// set so a restored VIEW_COOKIE_NAME cookie pointing at one of these can be sanitized back to
// the default when the page loads without ?deb, rather than landing on a view whose picker
// option is no longer there to switch away from.
const DEBUG_MAP_VIEWS = new Set<MapView>([
  "platesDetail", "speckle", "platesAndFaults", "geomorph", "elevReason", "overlapAge", "nodeAge",
  "plateInspector", "riverInspector", "lakeInspector",
]);

interface ViewCookie {
  projection: Projection;
  mapView: MapView;
  rotation: Mat3;
}

function loadViewCookie(): ViewCookie | null {
  const raw = getCookie(VIEW_COOKIE_NAME);
  if (!raw) return null;
  try {
    const parsed = JSON.parse(raw);
    if (!PROJECTION_CHOICES.has(parsed.projection) || !MAP_VIEW_CHOICES.has(parsed.mapView)) return null;
    if (!Array.isArray(parsed.rotation) || parsed.rotation.length !== 9 || !parsed.rotation.every((v: unknown) => typeof v === "number" && Number.isFinite(v))) {
      return null;
    }
    return { projection: parsed.projection, mapView: parsed.mapView, rotation: parsed.rotation };
  } catch {
    return null;
  }
}

// Read once at module load (this is a browser-only SPA bundle, so `document` is always
// available) rather than inside the component -- every useState initializer below just reads
// this same snapshot, so a page that never had a cookie set falls back to the same defaults
// it always used.
const initialView = loadViewCookie();

export default function App() {
  const [showGenerateDialog, setShowGenerateDialog] = useState(false);
  // "Random" (today's noise-driven generation, default) vs. "Human-made" (a drawn or loaded
  // coastline) vs. "Premade worlds" (one of PREMADE_WORLDS' built-in coastline sketches) vs.
  // "Debugging Worlds" (see below) -- see SketchEditor.tsx and backend app/worldsketch.py.
  // sketchImageDataUrl is the captured drawing/image/premade pick as a full data URL (so it
  // can double as an <img> preview source); handleGenerate strips its
  // `data:image/png;base64,` prefix before sending. "Premade worlds" rides the exact same
  // sketchImageDataUrl/generation path "Human-made" does -- picking one just pre-fills it with
  // a built-in image instead of a drawn/loaded one. None of "random"/"human"/"premade" resets
  // on a successful Generate, so re-opening the dialog to tweak the seed/detail keeps the same
  // drawing/pick around rather than discarding it.
  const [generateMode, setGenerateMode] = useState<"random" | "human" | "premade" | "debug">("random");
  const [sketchImageDataUrl, setSketchImageDataUrl] = useState<string | null>(null);
  // The selected Premade world's backendId (see premadeWorlds.ts) -- rides alongside
  // sketchImageDataUrl/seed to generateWorld's own premadeWorldId param so the backend can
  // ground plate placement/mantle convection in real data for Earth/Pangaea; null for Z&D
  // (and for every other tab, where it's simply unused).
  const [premadeWorldId, setPremadeWorldId] = useState<"earth" | "pangaea" | "got" | null>(null);
  const [showSketchEditor, setShowSketchEditor] = useState(false);
  // "Debugging Worlds" tab (see backend debug_worlds.py) -- tiny scripted plate scenarios for
  // fast iteration on the gap-filling problem. The scenario list is fetched once (static,
  // server-authoritative so the frontend never hardcodes it) and the picker defaults to
  // whichever comes first.
  const [debugScenarios, setDebugScenarios] = useState<DebugScenario[]>([]);
  const [debugScenario, setDebugScenario] = useState<string>("");
  useEffect(() => {
    if (!DEBUG_UI) return; // "Debugging Worlds" tab is hidden without ?deb -- nothing to populate it with
    fetchDebugScenarios()
      .then((r) => {
        setDebugScenarios(r.scenarios);
        setDebugScenario((cur) => cur || r.scenarios[0]?.name || "");
      })
      .catch(() => {
        // Best-effort -- if this fails the "Debugging Worlds" tab just shows an empty picker;
        // every other Generate World tab is unaffected.
      });
  }, []);
  // "Load an image" (Human-made tab) -- a hidden file input triggered by a plain button, read
  // via FileReader straight to a data URL, same shape SketchEditor's own "Done" produces.
  const loadImageInputRef = useRef<HTMLInputElement>(null);
  const handleLoadSketchFile = useCallback((e: React.ChangeEvent<HTMLInputElement>) => {
    const file = e.target.files?.[0];
    e.target.value = ""; // so picking the same file again still fires onChange
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      if (typeof reader.result === "string") setSketchImageDataUrl(reader.result);
    };
    reader.readAsDataURL(file);
  }, []);
  const [showAdvancedSettings, setShowAdvancedSettings] = useState(false);
  const [seed, setSeed] = useState(randomSeed());
  const [continentalPercent, setContinentalPercent] = useState(DEFAULT_CONTINENTAL_PERCENT);
  const [landPercent, setLandPercent] = useState(DEFAULT_LAND_PERCENT);
  const [axialTiltDeg, setAxialTiltDeg] = useState(DEFAULT_AXIAL_TILT_DEG);
  const [detail, setDetail] = useState(DEFAULT_DETAIL);
  // The Advanced-settings dialog's own "Fluid dynamics resolution" choice -- a separate dial
  // from `detail` above: unlike node_density/climate_density (both driven off `detail`, modulo
  // "High"'s own climateDensity override -- see climateDensityForDetail), this only affects
  // Ocean/Atmospheric Fluid Dynamics's own grid (see backend app/world.py's World.fluid_density),
  // not plate/climate resolution, so it's worth letting the user pick independently rather than
  // folding it into `detail` too. Defaults to DEFAULT_FLUID_DETAIL ("High", *this* dial's own
  // top choice), matching the backend's own FLUID_DENSITY_CHOICES cap -- see that constant's
  // own comment for why it's lower than DETAIL_CHOICES' own top choice.
  const [fluidDensity, setFluidDensity] = useState(DEFAULT_FLUID_DETAIL);
  const [initialSoilMaturityPercent, setInitialSoilMaturityPercent] = useState(DEFAULT_INITIAL_SOIL_MATURITY_PERCENT);
  const [autoPlates, setAutoPlates] = useState(true);
  const [numPlates, setNumPlates] = useState(DEFAULT_PLATES);
  const [voronoiPoints, setVoronoiPoints] = useState(DEFAULT_VORONOI_POINTS_RANDOM);

  const [stepYears, setStepYears] = useState(STEP_YEARS_OPTIONS[1]);
  const [projection, setProjection] = useState<Projection>(initialView?.projection ?? "eckert4");
  const [mapView, setMapView] = useState<MapView>(
    initialView && (DEBUG_UI || !DEBUG_MAP_VIEWS.has(initialView.mapView)) ? initialView.mapView : "combined",
  );
  // handleStep's own closure over `mapView` is captured when the step *starts*, so if the map
  // mode changes while that step is still in flight, its post-step refresh would otherwise
  // render the world's new (stepped) state back in the old, now-stale map mode -- overwriting
  // whatever the mode change's own refresh (the effect below) just showed. This ref always
  // holds the current value, so the post-step refresh renders whatever mode is selected by the
  // time the step actually finishes, not whatever was selected when it started.
  const mapViewRef = useRef(mapView);
  mapViewRef.current = mapView;
  // The map's current view orientation (see rotation.ts and
  // docs/simulation-model.md#rotating-the-view) -- lives entirely here, like
  // projection/mapView, sent fresh with every render call rather than stored server-side
  // (it's client-local view state, not simulation state). centerLatLon is derived display
  // state: it tracks `rotation` normally, but during an active drag MapCanvas overrides it
  // continuously via onRotationPreview, since the real legend is baked server-side into the
  // PNG and can't update mid-drag.
  const [rotation, setRotation] = useState<Mat3>(initialView?.rotation ?? IDENTITY_ROTATION);
  const [centerLatLon, setCenterLatLon] = useState(() => {
    if (!initialView) return { lat: 0, lon: 0 };
    const [latRad, lonRad] = centerOfRotation(initialView.rotation);
    return { lat: (latRad * 180) / Math.PI, lon: (lonRad * 180) / Math.PI };
  });
  // Press-and-hold on the "Center: ..." readout (see the render below) opens these editable
  // lat/lon fields so the view can be recentered on a typed coordinate instead of only by
  // dragging. centerEditValue holds the raw text (not numbers) while typing, so an in-progress
  // entry like "-" or "12." doesn't get clobbered by parsing/reformatting on every keystroke;
  // it's only parsed and validated on submit (see handleCommitCenterEdit).
  const [editingCenter, setEditingCenter] = useState(false);
  const [centerEditValue, setCenterEditValue] = useState({ lat: "", lon: "" });
  const [centerEditError, setCenterEditError] = useState<string | null>(null);
  // Legend-click-to-highlight (see Legend.tsx/MapCanvas.tsx) -- only ever meaningful on the
  // views whose legend swatches are clickable (Biome, Combined, "Last elevation change" and
  // "Plates & Faults"), so it's cleared any time the view changes away from all of them rather
  // than silently carrying a stale selection into a view whose legend can't reflect or clear it.
  const [highlightedBiome, setHighlightedBiome] = useState<string | null>(null);
  useEffect(() => {
    if (
      mapView !== "biome" && mapView !== "combined" && mapView !== "elevReason" && mapView !== "platesAndFaults"
    ) {
      setHighlightedBiome(null);
    }
  }, [mapView]);
  // Memoized so its identity only changes with the selection itself, not on every render --
  // MapCanvas.tsx's highlight-toggle effect is keyed directly on this object's identity.
  const highlightTarget = useMemo(
    () => (highlightedBiome ? highlightTargetFor(mapView, highlightedBiome) : null),
    [mapView, highlightedBiome],
  );
  // Elevation view's own "Mountains" / "Plains & Plateaus" legend toggles (see Legend.tsx) --
  // baked server-side into the rendered PNG (render_image.py's terrain-relief overlay), not a
  // client-side highlight like highlightedBiome above, so they need their own refresh() call
  // rather than a pure client-side redraw. Same "reset on leaving the view, don't persist
  // across a browser refresh" precedent as highlightedBiome.
  const [showMountains, setShowMountains] = useState(false);
  const [showPlainsPlateaus, setShowPlainsPlateaus] = useState(false);
  useEffect(() => {
    if (mapView !== "elevation") {
      setShowMountains(false);
      setShowPlainsPlateaus(false);
    }
  }, [mapView]);
  // refresh() reads these at call time (see its own definition below) rather than taking them
  // as arguments -- avoids threading two more params through every one of refresh()'s several
  // call sites, the same `mapViewRef` pattern used just above for the same reason.
  const terrainOverlayRef = useRef({ showMountains, showPlainsPlateaus });
  terrainOverlayRef.current = { showMountains, showPlainsPlateaus };
  const [summary, setSummary] = useState<WorldSummary | null>(null);
  const [renderData, setRenderData] = useState<RenderResponse | null>(null);
  // Plate Inspector's own data (see PlateInspector.tsx) -- true-frame/rotation-independent,
  // so unlike renderData it's only refetched on an actual world-state change (generate/step),
  // never on a projection/rotation change. selectedPlateId only makes sense within one
  // world's lifetime (plate ids aren't stable across a regenerate), so it resets there too.
  const [platesData, setPlatesData] = useState<PlateSummary[]>([]);
  const [selectedPlateId, setSelectedPlateId] = useState<number | null>(null);
  // River Inspector's own data (see RiverInspector.tsx) -- same true-frame/world-state-only
  // refresh pattern as platesData. river_id is only meaningful against the most recent
  // /world/rivers response (rivers are regrouped fresh every call, no persistent identity --
  // see backend app/hydrology.py's group_rivers), so it's reset on every generate *and* step,
  // not just generate like selectedPlateId.
  const [riversData, setRiversData] = useState<RiverSummary[]>([]);
  const [selectedRiverId, setSelectedRiverId] = useState<number | null>(null);
  // The land/lake-vs-ocean boundary (see backend app/coastline.py), fetched alongside rivers
  // from the same /world/rivers response -- the River Inspector has no other land/ocean cue.
  const [coastlineSegments, setCoastlineSegments] = useState<Segment[]>([]);
  // Lake Inspector's own data (see LakeInspector.tsx) -- same true-frame/world-state-only
  // refresh pattern as riversData, and reuses the same coastlineSegments fetched above rather
  // than a second copy (both /world/rivers and /world/lakes compute the identical boundary).
  // `selectedBasin` is whatever's currently displayed -- either one of `lakesData`
  // (`selectedBasin.is_lake`) or a dry basin/no-basin/ocean result from a land click, which
  // isn't itself a member of `lakesData` (see api.ts's LakeAtResponse) -- `selectedKind` carries
  // which of those a `null` `selectedBasin` actually means (nothing clicked yet vs. an ocean or
  // no-basin click). Reset on every generate *and* step, same as selectedRiverId: lake_id is
  // only meaningful against the most recent /world/lakes response.
  const [lakesData, setLakesData] = useState<LakeSummary[]>([]);
  const [selectedBasin, setSelectedBasin] = useState<LakeSummary | null>(null);
  const [selectedBasinKind, setSelectedBasinKind] = useState<LakeAtResponse["kind"] | null>(null);
  // Intraplate fault + fault-system + activity data for the "Plates & Faults" view -- same
  // true-frame/world-state-only refresh pattern as platesData. Faults aren't individually
  // selectable there (plate selection is), so there's no fault-id selection state.
  const [faultsData, setFaultsData] = useState<FaultSummary[]>([]);
  const [faultSystemsData, setFaultSystemsData] = useState<FaultSystemSummary[]>([]);
  // Recent earthquakes + current volcano vents for the "Plates & Faults" view's activity
  // overlay (see faults.Earthquake / GET /world/earthquakes, volcanism.py / GET
  // /world/volcanoes). Refreshed alongside faultsData; both are pure overlays, nothing to
  // select. `showQuakesVolcanoes` is the sidebar checkbox that toggles the whole overlay.
  const [earthquakesData, setEarthquakesData] = useState<EarthquakeSummary[]>([]);
  const [volcanoesData, setVolcanoesData] = useState<VolcanoSummary[]>([]);
  const [showQuakesVolcanoes, setShowQuakesVolcanoes] = useState(true);
  // The Elevation & Biome / Elevation / Biome views' click-to-inspect popup (see
  // MapCanvas.tsx's onProbe and the popup JSX below). `displayX`/`displayY` place it over the
  // map in CSS pixels; `sample` fills in once GET /world/sample_at resolves. Cleared on any
  // view/projection/rotation change and on every render refresh (a step moves everything
  // that was under it), so a stale popup never lingers over ground it no longer describes.
  const [probe, setProbe] = useState<
    | { displayX: number; displayY: number; latDeg: number; lonDeg: number; status: "loading" | "ok" | "error"; sample: PointSample | null }
    | null
  >(null);
  // Tags each sample_at request so a slow earlier response (or a dismissal) can't overwrite a
  // newer click's -- same monotonic-id guard as renderRequestIdRef below.
  const probeRequestIdRef = useRef(0);
  const handleProbe = useCallback(
    (next: { displayX: number; displayY: number; latDeg: number; lonDeg: number } | null) => {
      const requestId = ++probeRequestIdRef.current;
      if (!next) {
        setProbe(null);
        return;
      }
      setProbe({ ...next, status: "loading", sample: null });
      fetchPointSample(next.latDeg, next.lonDeg)
        .then((sample) => {
          if (requestId === probeRequestIdRef.current) setProbe({ ...next, status: "ok", sample });
        })
        .catch(() => {
          // A click that races a generate/step (no world yet), or a dropped request -- show
          // the failure in place rather than a silently empty popup.
          if (requestId === probeRequestIdRef.current) setProbe({ ...next, status: "error", sample: null });
        });
    },
    [],
  );
  // Any of these means the popup's anchor no longer maps to the same ground -- rotation/
  // projection/view change moves the map under it, and a new renderData means a step or a
  // Controls edit just changed what's there -- so drop it rather than leave it floating.
  useEffect(() => {
    probeRequestIdRef.current++;
    setProbe(null);
  }, [mapView, projection, rotation, renderData]);
  // The "Added/Removed Points" (nodeAge) view's own click-to-inspect popup -- same shape and
  // lifecycle as `probe` above, kept as a separate state (rather than widening PointSample's
  // union) since the two views' popups show entirely different fields.
  const [nodeProbe, setNodeProbe] = useState<
    | { displayX: number; displayY: number; latDeg: number; lonDeg: number; status: "loading" | "ok" | "error"; result: NodeAtResponse | null }
    | null
  >(null);
  const nodeProbeRequestIdRef = useRef(0);
  const handleNodeProbe = useCallback(
    (next: { displayX: number; displayY: number; latDeg: number; lonDeg: number } | null) => {
      const requestId = ++nodeProbeRequestIdRef.current;
      if (!next) {
        setNodeProbe(null);
        return;
      }
      setNodeProbe({ ...next, status: "loading", result: null });
      fetchNodeAt(next.latDeg, next.lonDeg)
        .then((result) => {
          if (requestId === nodeProbeRequestIdRef.current) setNodeProbe({ ...next, status: "ok", result });
        })
        .catch(() => {
          if (requestId === nodeProbeRequestIdRef.current) setNodeProbe({ ...next, status: "error", result: null });
        });
    },
    [],
  );
  useEffect(() => {
    nodeProbeRequestIdRef.current++;
    setNodeProbe(null);
  }, [mapView, projection, rotation, renderData]);
  // The "Points" (platesDetail) debug view's click-to-inspect + arrow-key navigation (see
  // MapCanvas.tsx's highlightLine prop and fetchElevationPointAt/fetchElevationPoint in
  // api.ts). Kept as its own probe, same reasoning as nodeProbe above (a different view's
  // popup shows entirely different fields) -- but unlike every other probe here, this one
  // also drives a keyboard-navigable selection rather than being purely mouse-driven, so the
  // full ElevationPointResponse (not just the fields the popup renders) is kept in state:
  // arrow-key stepping below reads plate_id/line_index/point_index straight back out of it.
  const [pointProbe, setPointProbe] = useState<
    | {
        displayX: number; displayY: number; latDeg: number; lonDeg: number;
        status: "loading" | "ok" | "error"; result: ElevationPointResponse | null;
      }
    | null
  >(null);
  const pointProbeRequestIdRef = useRef(0);
  // Focused right after a successful click (see handlePointProbe) so ArrowLeft/Right and
  // Shift+ArrowLeft/Right work immediately without an extra click on the map first --
  // same "focus on selection" pattern PlateInspector's own containerRef uses for Tab/Shift+Tab.
  const pointContainerRef = useRef<HTMLDivElement>(null);
  const handlePointProbe = useCallback(
    (next: { displayX: number; displayY: number; latDeg: number; lonDeg: number } | null) => {
      const requestId = ++pointProbeRequestIdRef.current;
      if (!next) {
        setPointProbe(null);
        return;
      }
      setPointProbe({ ...next, status: "loading", result: null });
      fetchElevationPointAt(next.latDeg, next.lonDeg)
        .then((result) => {
          if (requestId === pointProbeRequestIdRef.current) {
            setPointProbe({ ...next, status: "ok", result });
            pointContainerRef.current?.focus();
          }
        })
        .catch(() => {
          if (requestId === pointProbeRequestIdRef.current) setPointProbe({ ...next, status: "error", result: null });
        });
    },
    [],
  );
  useEffect(() => {
    pointProbeRequestIdRef.current++;
    setPointProbe(null);
  }, [mapView, projection, rotation, renderData]);
  // ArrowLeft/Right steps to the previous/next point on the selected line (wrapping);
  // Shift+ArrowLeft/Right steps to the previous/next line on the same plate instead, ordered
  // by ascending phi (also wrapping) -- see backend plates.sorted_nonempty_lines -- keeping
  // the same point_index (the server clamps it into the new line's own range, see
  // GET /world/elevation_point). No-op with nothing selected, or outside the Points view.
  const handlePointKeyDown = useCallback(
    (e: KeyboardEvent<HTMLDivElement>) => {
      if (mapView !== "platesDetail") return;
      const result = pointProbe?.result;
      if (!result) return;
      const dir = e.key === "ArrowRight" ? 1 : e.key === "ArrowLeft" ? -1 : 0;
      if (dir === 0) return;
      e.preventDefault();
      const requestId = ++pointProbeRequestIdRef.current;
      const lineIndex = e.shiftKey
        ? (result.line.line_index + dir + result.line.num_lines) % result.line.num_lines
        : result.line.line_index;
      const pointIndex = e.shiftKey
        ? result.point.index
        : (result.point.index + dir + result.line.num_points) % result.line.num_points;
      fetchElevationPoint(result.plate_id, lineIndex, pointIndex).then((next) => {
        if (requestId === pointProbeRequestIdRef.current) {
          setPointProbe((cur) => (cur ? { ...cur, status: "ok", result: next } : cur));
        }
      });
    },
    [mapView, pointProbe],
  );
  // Stats panel data (see StatsModal.tsx) -- `stats` is the latest snapshot, `statsHistory`
  // accumulates one entry per generate/step (deduped by elapsed_years) for the panel's graph
  // tabs. Recorded continuously, not just while the modal is open, so opening it later still
  // shows the full history since the world was generated. The backend now keeps its own copy
  // too (World.stats_history, persisted in every saved .mbworld -- see backend app/world.py),
  // which handleWorldReplaced below restores this from after a Load -- otherwise a loaded
  // save would resume with an empty chart despite however much history it was saved with.
  const [stats, setStats] = useState<WorldStats | null>(null);
  const [statsHistory, setStatsHistory] = useState<WorldStats[]>([]);
  const [showStatsModal, setShowStatsModal] = useState(false);
  // Live world controls (see ControlsModal.tsx and backend app/world.py's World.sea_level_m/
  // World.solar_multiplier) -- reset to their defaults on every fresh Generate, same as
  // every other generation-time value here, even though these two are adjustable afterward
  // too (unlike the others) via /world/controls.
  const [seaLevelM, setSeaLevelM] = useState(DEFAULT_SEA_LEVEL_M);
  const [solarMultiplier, setSolarMultiplier] = useState(DEFAULT_SOLAR_MULTIPLIER);
  // Full period (years) of the glacial<->interglacial cycle; 0 == disabled. Live-adjustable
  // via Controls like seaLevelM/solarMultiplier. See backend app/world.py's World.ice_age_period_years.
  const [iceAgePeriodYears, setIceAgePeriodYears] = useState(DEFAULT_ICE_AGE_PERIOD_YEARS);
  // Same live-adjustable-via-Controls pattern as seaLevelM/solarMultiplier above -- lets the
  // user run plate tectonics only, climate & biomes only, or (the default) both together. See
  // backend app/world.py's World.simulate_plate_movement/World.simulate_climate_biomes.
  const [simulatePlateMovement, setSimulatePlateMovement] = useState(DEFAULT_SIMULATE_PLATE_MOVEMENT);
  const [simulateClimateBiomes, setSimulateClimateBiomes] = useState(DEFAULT_SIMULATE_CLIMATE_BIOMES);
  // "cfd" (shallow-water solve) or "diagnostic" (fast closed-form ABL wind) -- see backend
  // app/world.py's World.wind_model. Live-adjustable via Controls like the toggles above.
  const [windModel, setWindModel] = useState(DEFAULT_WIND_MODEL);
  // "boundary" / "fault" / "both" -- see backend app/world.py's World.fault_deformation_mode.
  // Live-adjustable via Controls like windModel.
  const [faultDeformationMode, setFaultDeformationMode] = useState(DEFAULT_FAULT_DEFORMATION_MODE);
  // Gate for the verbose _fill_corner_notch decision log -- see backend World.debug_diagnostics
  // / GET /world/corner_notch_log. Live-adjustable via Controls like windModel.
  const [debugDiagnostics, setDebugDiagnostics] = useState(DEFAULT_DEBUG_DIAGNOSTICS);
  // The corner-notch log itself, refreshed alongside faultsData (see refreshCornerNotchLog) --
  // populated regardless of debugDiagnostics' current value (cheap to fetch; stays empty when
  // off), so toggling the flag on mid-session shows entries from the very next step.
  const [cornerNotchLog, setCornerNotchLog] = useState<CornerNotchLogEntry[]>([]);
  // Geomorphic-budget tuning knobs (see DEFAULT_TUNING / backend World's *_multiplier
  // group) -- one object of dimensionless multipliers, live-adjustable via Controls, reset
  // to all-1.0 on a fresh Generate and synced from the loaded world on Load.
  const [tuning, setTuning] = useState<TuningMultipliers>(DEFAULT_TUNING);
  const [showControlsModal, setShowControlsModal] = useState(false);
  const [showFileModal, setShowFileModal] = useState(false);
  // The Record toolbar button's dialog (see AnimationModal.tsx) -- separate from showFileModal
  // now that recording isn't nested inside "File...".
  const [showAnimationModal, setShowAnimationModal] = useState(false);
  // The div wrapping whichever map view component is currently mounted -- see
  // FileModal.tsx's own "Save Image" comment for why it reads the live <canvas> straight
  // out of this DOM node rather than a ref threaded through four separate components.
  const mapWrapperRef = useRef<HTMLDivElement>(null);
  const [playing, setPlaying] = useState(false);
  const [busy, setBusy] = useState(false);
  const [stepping, setStepping] = useState(false);
  const [error, setError] = useState<string | null>(null);
  // Background animation (recording) -- see AnimationModal.tsx and handleStartAnimation
  // below. `animation` is non-null only while a run is in flight -- it holds the live frame
  // count for the sidebar progress display, but that display only makes sense for a
  // fixed-frame-count run (`unbounded: false`): "frame 12 of 20" is meaningful, "frame 12 of
  // (whenever I press Stop)" isn't, so the unbounded case relies on the main map's live
  // preview plus the ticking "elapsed" readout instead (see the sidebar render below).
  // `animationResult` holds the finished MP4 (+ world summary) until the user saves or
  // discards it via SaveAnimationModal, which pops up as soon as it's set -- reached whether
  // the run finished on its own or the toolbar's Stop button ended it early, so ending a
  // recording always leads straight to that save-or-discard choice rather than a sidebar
  // panel to go looking for. The run streams each frame's PNG straight onto the main map, and
  // holds the server's world lock throughout, so every world-mutating / rendering control is
  // disabled while `animation` is set (see `animating`). The toolbar's "Stop" button goes
  // through api.ts's stopAnimation() rather than aborting the request, so the stream ends
  // cleanly and still returns the video (see handleStopAnimation below).
  const [animation, setAnimation] = useState<{ frame: number; total: number; unbounded: boolean } | null>(null);
  const [animationResult, setAnimationResult] = useState<AnimateResponse | null>(null);
  const animating = animation !== null;

  // Two refresh() calls can be in flight at once -- e.g. changing map mode while a step is
  // in flight fires one from the mode-change effect below (for the pre-step world) and
  // another from handleStep's own post-step Promise.all (for the post-step world) -- and
  // nothing guarantees they resolve in the order they were issued (server thread-pool
  // scheduling, retries, network jitter). Without a guard, whichever response happens to
  // land last wins even if it's the stale one, which is what let the map mode visibly revert
  // when a step's own render landed after a mode-change render that was issued earlier but
  // resolved later. This ref tags every call with a monotonic id and only ever commits the
  // response from the most recently *issued* call, so a stale response is silently dropped
  // instead of overwriting a newer one.
  const renderRequestIdRef = useRef(0);
  const refresh = useCallback(async (proj: Projection, view: MapView, viewRotation: Mat3) => {
    if (view === "plateInspector" || view === "riverInspector" || view === "lakeInspector" || view === "platesAndFaults") return; // none of these use renderData -- see below
    const requestId = ++renderRequestIdRef.current;
    const { showMountains: mountains, showPlainsPlateaus: plainsPlateaus } = terrainOverlayRef.current;
    try {
      const data = await renderWorld(proj, view, RENDER_WIDTH, RENDER_HEIGHT, viewRotation, mountains, plainsPlateaus);
      if (requestId === renderRequestIdRef.current) setRenderData(data);
    } catch (e) {
      if (requestId === renderRequestIdRef.current) setError(String(e));
    }
  }, []);

  const refreshPlates = useCallback(async () => {
    try {
      const data = await fetchPlates();
      setPlatesData(data.plates);
    } catch (e) {
      setError(String(e));
    }
  }, []);

  const refreshRivers = useCallback(async () => {
    try {
      const data = await fetchRivers();
      setRiversData(data.rivers);
      setCoastlineSegments(data.coastline_segments);
    } catch (e) {
      setError(String(e));
    }
  }, []);

  const refreshLakes = useCallback(async () => {
    try {
      const data = await fetchLakes();
      setLakesData(data.lakes);
    } catch (e) {
      setError(String(e));
    }
  }, []);

  const refreshFaults = useCallback(async () => {
    try {
      const [faults, quakes, volcanoes] = await Promise.all([fetchFaults(), fetchEarthquakes(), fetchVolcanoes()]);
      setFaultsData(faults.faults);
      setFaultSystemsData(faults.fault_systems);
      setEarthquakesData(quakes.earthquakes);
      setVolcanoesData(volcanoes.volcanoes);
    } catch (e) {
      setError(String(e));
    }
  }, []);

  // Best-effort, same spirit as recordStats -- a failed fetch here shouldn't surface as the
  // main error line or block generate/step, since this is a debug-only side panel.
  const refreshCornerNotchLog = useCallback(async () => {
    try {
      const log = await fetchCornerNotchLog();
      setCornerNotchLog(log.entries);
    } catch {
      // ignored -- see comment above
    }
  }, []);

  // Shared by recordStats below and the animation progress handler (see handleStartAnimation)
  // -- both just land a WorldStats snapshot, one fetched, one riding along an animate() stream
  // line, so both dedupe against the history's last entry the same way.
  const applyStats = useCallback((s: WorldStats) => {
    setStats(s);
    setStatsHistory((prev) =>
      prev.length > 0 && prev[prev.length - 1].elapsed_years === s.elapsed_years ? prev : [...prev, s],
    );
  }, []);

  // Stats are a secondary/best-effort feature -- a failed fetch here (e.g. a transient
  // network blip) shouldn't surface as the main error line or block generate/step, unlike
  // refresh/refreshPlates above which are core to the map actually updating.
  const recordStats = useCallback(async () => {
    try {
      const s = await fetchStats();
      applyStats(s);
    } catch {
      // ignored -- see comment above
    }
  }, [applyStats]);

  const handleGenerate = useCallback(async () => {
    setBusy(true);
    setError(null);
    try {
      // The "Human-made"/"Premade worlds" tabs' sketch rides along as a bare base64 payload --
      // strip the data URL's `data:image/png;base64,` prefix generateWorld's caller
      // (SketchEditor's toDataURL / the file-picker's FileReader / premadeWorlds.ts) all produce.
      const sketchBase64 =
        (generateMode === "human" || generateMode === "premade") && sketchImageDataUrl
          ? sketchImageDataUrl.split(",", 2)[1] ?? null
          : null;
      const s =
        generateMode === "debug"
          ? await generateDebugWorld(debugScenario, seed)
          : await generateWorld(
              seed, continentalPercent / 100, landPercent / 100, axialTiltDeg, detail, initialSoilMaturityPercent / 100,
              climateDensityForDetail(detail), fluidDensity, autoPlates ? null : numPlates, voronoiPoints, sketchBase64,
              generateMode === "premade" ? premadeWorldId : null,
            );
      setSummary(s);
      setSelectedPlateId(null);
      setSelectedRiverId(null);
      setSelectedBasin(null);
      setSelectedBasinKind(null);
      setShowGenerateDialog(false);
      setStatsHistory([]); // plate ids and elapsed_years both reset with a fresh world
      setSeaLevelM(DEFAULT_SEA_LEVEL_M); // live controls reset with a fresh world too
      setSolarMultiplier(DEFAULT_SOLAR_MULTIPLIER);
      setIceAgePeriodYears(DEFAULT_ICE_AGE_PERIOD_YEARS);
      setSimulatePlateMovement(DEFAULT_SIMULATE_PLATE_MOVEMENT);
      setSimulateClimateBiomes(DEFAULT_SIMULATE_CLIMATE_BIOMES);
      setWindModel(DEFAULT_WIND_MODEL);
      setFaultDeformationMode(DEFAULT_FAULT_DEFORMATION_MODE);
      // A debug world starts with diagnostics already on server-side (see debug_worlds.py) --
      // match that here rather than resetting to the ordinary default.
      setDebugDiagnostics(generateMode === "debug" ? true : DEFAULT_DEBUG_DIAGNOSTICS);
      setTuning(DEFAULT_TUNING);
      await Promise.all([refresh(projection, mapView, rotation), refreshPlates(), refreshRivers(), refreshLakes(), refreshFaults(), refreshCornerNotchLog(), recordStats()]);
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }, [
    seed, continentalPercent, landPercent, axialTiltDeg, detail, fluidDensity, initialSoilMaturityPercent, autoPlates, numPlates, voronoiPoints,
    generateMode, sketchImageDataUrl, premadeWorldId, debugScenario, projection, mapView, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, refreshFaults, refreshCornerNotchLog, recordStats,
  ]);


  // Debounced so dragging a Controls slider doesn't fire a network request (and force a
  // climate recompute, see main.py's /world/controls) on every single pixel of movement --
  // only once movement has paused briefly. Local slider state (seaLevelM/solarMultiplier)
  // still updates immediately on every change, so the slider itself never feels laggy.
  const controlsDebounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  // Accumulates every touched control between debounce fires (and deep-merges `tuning`), so
  // adjusting two different knobs inside the 150ms window doesn't drop the first one.
  const pendingControlsRef = useRef<{
    seaLevelM?: number;
    solarMultiplier?: number;
    iceAgePeriodYears?: number;
    simulatePlateMovement?: boolean;
    simulateClimateBiomes?: boolean;
    windModel?: string;
    faultDeformationMode?: string;
    debugDiagnostics?: boolean;
    tuning?: Partial<TuningMultipliers>;
  }>({});
  const pushControls = useCallback((next: typeof pendingControlsRef.current) => {
    const pending = pendingControlsRef.current;
    pendingControlsRef.current = { ...pending, ...next, tuning: { ...pending.tuning, ...next.tuning } };
    if (controlsDebounceRef.current) clearTimeout(controlsDebounceRef.current);
    controlsDebounceRef.current = setTimeout(async () => {
      const payload = pendingControlsRef.current;
      pendingControlsRef.current = {};
      try {
        await updateControls(payload);
        await Promise.all([refresh(projection, mapViewRef.current, rotation), recordStats()]);
      } catch (e) {
        setError(String(e));
      }
    }, 150);
  }, [projection, rotation, refresh, recordStats]);

  const handleSeaLevelChange = useCallback((v: number) => {
    setSeaLevelM(v);
    pushControls({ seaLevelM: v });
  }, [pushControls]);

  const handleTuningChange = useCallback((key: TuningKey, v: number) => {
    setTuning((prev) => ({ ...prev, [key]: v }));
    pushControls({ tuning: { [key]: v } });
  }, [pushControls]);

  const handleTuningReset = useCallback(() => {
    setTuning(DEFAULT_TUNING);
    pushControls({ tuning: DEFAULT_TUNING });
  }, [pushControls]);

  // Checkboxes, not sliders, so no dragging concern -- still routed through the same
  // (harmlessly short) debounce as the sliders above, for one shared code path.
  const handleSimulatePlateMovementChange = useCallback((v: boolean) => {
    setSimulatePlateMovement(v);
    pushControls({ simulatePlateMovement: v });
  }, [pushControls]);

  const handleSimulateClimateBiomesChange = useCallback((v: boolean) => {
    setSimulateClimateBiomes(v);
    pushControls({ simulateClimateBiomes: v });
  }, [pushControls]);

  const handleSolarMultiplierChange = useCallback((v: number) => {
    setSolarMultiplier(v);
    pushControls({ solarMultiplier: v });
  }, [pushControls]);

  const handleIceAgePeriodChange = useCallback((v: number) => {
    setIceAgePeriodYears(v);
    pushControls({ iceAgePeriodYears: v });
  }, [pushControls]);

  const handleWindModelChange = useCallback((v: string) => {
    setWindModel(v);
    pushControls({ windModel: v });
  }, [pushControls]);

  const handleFaultDeformationModeChange = useCallback((v: string) => {
    setFaultDeformationMode(v);
    pushControls({ faultDeformationMode: v });
  }, [pushControls]);

  const handleDebugDiagnosticsChange = useCallback((v: boolean) => {
    setDebugDiagnostics(v);
    pushControls({ debugDiagnostics: v });
  }, [pushControls]);

  const handleStep = useCallback(async () => {
    if (!summary) return;
    setStepping(true);
    setError(null);
    try {
      const s = await stepWorld(stepYears);
      setSummary(s);
      setSelectedRiverId(null); // rivers are regrouped fresh every step -- a stale id could point at an unrelated network
      setSelectedBasin(null); // lakes are regrouped fresh every step too -- same reasoning
      setSelectedBasinKind(null);
      // mapViewRef.current, not mapView -- see the ref's own comment above.
      await Promise.all([
        refresh(projection, mapViewRef.current, rotation), refreshPlates(), refreshRivers(), refreshLakes(), refreshFaults(), refreshCornerNotchLog(), recordStats(),
      ]);
    } catch (e) {
      setError(String(e));
      setPlaying(false);
    } finally {
      setStepping(false);
    }
  }, [summary, stepYears, projection, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, refreshFaults, refreshCornerNotchLog, recordStats]);

  // FileModal's "Load World" -- a loaded world fully replaces the current one, same as a
  // fresh Generate (see handleGenerate above), plus syncing every live Controls value
  // (seaLevelM/solarMultiplier/simulatePlateMovement/simulateClimateBiomes) to the *loaded*
  // world's own real values: calling updateControls with no fields set changes nothing but
  // still returns the current world's current values (see api.ts's updateControls/backend
  // app/main.py's /world/controls), which is simpler than adding a new endpoint just to read
  // them back. statsHistory is restored the same way, from GET /world/stats_history -- the
  // loaded save's own recorded series (see World.stats_history), not reset to empty the way
  // it used to be before the backend persisted this itself.
  const handleWorldReplaced = useCallback(async (s: WorldSummary) => {
    setBusy(true);
    setError(null);
    try {
      setSummary(s);
      setSelectedPlateId(null);
      setSelectedRiverId(null);
      setSelectedBasin(null);
      setSelectedBasinKind(null);
      // Awaited (not fire-and-forget) so it lands before the Promise.all below's own
      // recordStats() call appends this instant's reading on top of it -- otherwise the two
      // setStatsHistory calls could resolve in either order and recordStats()'s plain-array
      // append would risk being clobbered by this one landing second.
      try {
        const { history } = await fetchStatsHistory();
        setStatsHistory(history);
      } catch {
        // best-effort, same spirit as recordStats below -- see its own comment
      }
      const controls = await updateControls({});
      setSeaLevelM(controls.sea_level_m);
      setSolarMultiplier(controls.solar_multiplier);
      setIceAgePeriodYears(controls.ice_age_period_years);
      setSimulatePlateMovement(controls.simulate_plate_movement);
      setSimulateClimateBiomes(controls.simulate_climate_biomes);
      setWindModel(controls.wind_model);
      setFaultDeformationMode(controls.fault_deformation_mode);
      setDebugDiagnostics(controls.debug_diagnostics);
      setTuning(Object.fromEntries(TUNING_MULTIPLIER_KEYS.map((k) => [k, controls[k]])) as TuningMultipliers);
      await Promise.all([refresh(projection, mapView, rotation), refreshPlates(), refreshRivers(), refreshLakes(), refreshFaults(), refreshCornerNotchLog(), recordStats()]);
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }, [projection, mapView, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, refreshFaults, refreshCornerNotchLog, recordStats]);

  // Post-animation refresh -- an animation run already advanced the world for real (see
  // api.ts's animateWorld), so this just runs the same post-step refresh handleStep does.
  const handleWorldAdvanced = useCallback(async (s: WorldSummary) => {
    setSummary(s);
    setSelectedRiverId(null);
    setSelectedBasin(null);
    setSelectedBasinKind(null);
    await Promise.all([refresh(projection, mapView, rotation), refreshPlates(), refreshRivers(), refreshLakes(), refreshFaults(), refreshCornerNotchLog(), recordStats()]);
  }, [projection, mapView, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, refreshFaults, refreshCornerNotchLog, recordStats]);

  // The Record toolbar button's AnimationModal -- run the whole recording in the background.
  // Each streamed frame's PNG is painted straight onto the main map (the run holds the server
  // world lock, so a normal render would 503 -- see api.ts/main.py), so the main display
  // doubles as the animation preview; every world-mutating control is disabled until it
  // finishes or the toolbar's Stop button ends it (see `animating`, handleStopAnimation).
  const handleStartAnimation = useCallback(
    async ({ numFrames, stepsPerFrame, unbounded }: { numFrames: number; stepsPerFrame: number; unbounded: boolean }) => {
      setShowAnimationModal(false);
      setError(null);
      setAnimationResult(null);
      setAnimation({ frame: 0, total: numFrames, unbounded });
      const view = mapViewRef.current;
      try {
        const result = await animateWorld(
          projection, view, RENDER_WIDTH, RENDER_HEIGHT, rotation, stepYears, stepsPerFrame, numFrames,
          (p) => {
            setAnimation({ frame: p.frame, total: p.total, unbounded });
            // Keep the sidebar's "elapsed" readout climbing frame by frame instead of it
            // sitting frozen at the pre-recording value until the whole run resolves.
            setSummary((prev) => (prev ? { ...prev, elapsed_years: p.elapsedYears } : prev));
            // Likewise for the Stats panel -- GET /world/stats would just block on the
            // world lock for the whole run, so this stream-carried snapshot (see api.ts's
            // AnimateProgress) is what lets it update live instead (GitHub issue #155).
            if (p.stats) applyStats(p.stats);
            if (p.imageBase64) {
              // Paint the frame onto the main map. renderRequestIdRef is bumped so any
              // pre-animation refresh() still in flight can't clobber it afterward.
              renderRequestIdRef.current++;
              setRenderData({ projection, elapsed_years: p.elapsedYears, image_base64: p.imageBase64 });
            }
          },
        );
        // Reached whether the run finished on its own or the toolbar's Stop button ended it
        // early (see api.ts's stopAnimation -- the stream still finishes cleanly and returns
        // the video either way, so both cases resolve here rather than throwing).
        setAnimationResult(result);
        await handleWorldAdvanced(result);
      } catch (e) {
        setError(String(e));
      } finally {
        setAnimation(null);
      }
    },
    [projection, rotation, stepYears, handleWorldAdvanced, applyStats],
  );

  // The toolbar's "Stop" button while recording -- asks the server to end the animation
  // cleanly after its current frame (see api.ts's stopAnimation) so it still finishes with a
  // complete, saveable video rather than being abandoned mid-run. Fire-and-forget: the
  // animateWorld() call above is what actually resolves once the stream's final "done"
  // message arrives; a failure to even deliver the stop request surfaces as `error` so it
  // isn't silently swallowed.
  const handleStopAnimation = useCallback(() => {
    stopAnimation().catch((e) => setError(String(e)));
  }, []);

  // The save-or-discard dialog's "Save" button (see SaveAnimationModal.tsx) -- downloads the
  // finished MP4, then closes the dialog same as Discard would.
  const handleSaveAnimation = useCallback(async () => {
    if (!animationResult) return;
    const { mime, videoBase64, elapsed_years, seed: animSeed } = animationResult;
    const blob = await (await fetch(`data:${mime};base64,${videoBase64}`)).blob();
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = `mantle-bloom-seed${animSeed}-animation-${Math.round(elapsed_years)}y.mp4`;
    a.click();
    URL.revokeObjectURL(url);
    setAnimationResult(null);
  }, [animationResult]);

  // The save-or-discard dialog's "Discard" button -- throws away the finished video. Doesn't
  // touch the world itself, which already permanently advanced when the recording ran.
  const handleDiscardAnimation = useCallback(() => {
    setAnimationResult(null);
  }, []);

  // Resets the view orientation back to the default (lat=0/lon=0, see rotation.ts) -- just
  // updates state, same as a completed drag; the effect below does the actual re-fetch.
  const handleRecenter = useCallback(() => {
    setRotation(IDENTITY_ROTATION);
    setCenterLatLon({ lat: 0, lon: 0 });
  }, []);

  // Click on the "Center: ..." readout opens the editable lat/lon popup (see the editingCenter
  // state above).
  const handleOpenCenterEdit = useCallback(() => {
    if (busy || !summary || animating || editingCenter) return;
    setCenterEditValue({ lat: centerLatLon.lat.toFixed(4), lon: centerLatLon.lon.toFixed(4) });
    setCenterEditError(null);
    setEditingCenter(true);
  }, [busy, summary, animating, editingCenter, centerLatLon]);

  const handleCancelCenterEdit = useCallback(() => {
    setEditingCenter(false);
    setCenterEditError(null);
  }, []);

  // Validates the typed lat/lon before committing it as the new view center -- same
  // rotationForCenter math MapCanvas's drag-to-pan gesture uses to commit a completed drag
  // (see rotationDrag.ts), just fed a typed target instead of a dragged one.
  const handleCommitCenterEdit = useCallback(() => {
    const latText = centerEditValue.lat.trim();
    const lonText = centerEditValue.lon.trim();
    const lat = Number(latText);
    const lon = Number(lonText);
    if (
      latText === "" || lonText === "" || !Number.isFinite(lat) || !Number.isFinite(lon) ||
      lat < -90 || lat > 90 || lon < -180 || lon > 180
    ) {
      setCenterEditError("Latitude must be -90 to 90 and longitude -180 to 180.");
      return;
    }
    setRotation(rotationForCenter((lat * Math.PI) / 180, (lon * Math.PI) / 180));
    setCenterLatLon({ lat, lon });
    setEditingCenter(false);
    setCenterEditError(null);
  }, [centerEditValue]);

  const selectedPlate = platesData.find((p) => p.plate_id === selectedPlateId) ?? null;
  const selectedRiver = riversData.find((r) => r.river_id === selectedRiverId) ?? null;

  // Re-render with the current world whenever the projection, map view, view rotation, or the
  // Elevation view's terrain-relief toggles change -- all baked server-side into the returned
  // image (see api.ts's renderWorld). A completed drag (see MapCanvas's onRotationCommitted
  // below) just updates `rotation` state; this effect is what actually re-fetches, the same
  // pattern projection/mapView already used before rotation existed.
  useEffect(() => {
    if (summary) {
      refresh(projection, mapView, rotation);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [projection, mapView, rotation, showMountains, showPlainsPlateaus]);

  // Persists projection/mapView/rotation to VIEW_COOKIE_NAME on every change, so the next
  // page load's `initialView` (see loadViewCookie above) picks up wherever the user left off.
  useEffect(() => {
    setCookie(VIEW_COOKIE_NAME, JSON.stringify({ projection, mapView, rotation }));
  }, [projection, mapView, rotation]);

  // Restores the map after a browser refresh: the world itself lives entirely in server
  // memory (see backend main.py's module docstring), so it's often still there even though
  // this component's own state (summary, platesData, ...) always starts blank on a fresh
  // mount. If /world/summary finds one, this replaces the blank starting state with it (via
  // the same "replace the current world" path FileModal's Load World uses) using the
  // projection/mapView/rotation the cookie above just restored; a 404 (no world generated
  // yet this server session) just leaves the normal blank/Generate-dialog state in place.
  useEffect(() => {
    let cancelled = false;
    (async () => {
      try {
        const s = await fetchWorldSummary();
        if (!cancelled) await handleWorldReplaced(s);
      } catch {
        // no world in server memory -- nothing to restore
      }
    })();
    return () => {
      cancelled = true;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const stepRef = useRef(handleStep);
  stepRef.current = handleStep;

  // Self-scheduling rather than setInterval: each step must finish before the next is
  // scheduled, so a slow step (larger worlds take longer -- a step that also lands on a
  // regularize/gap-fill or reassignment interval, or triggers a merge/split, can take
  // noticeably longer than a routine one) can never overlap with the next tick. A
  // fixed-cadence setInterval would keep firing regardless of whether the previous request
  // had returned.
  useEffect(() => {
    if (!playing) return;
    let cancelled = false;
    let timeoutId: ReturnType<typeof setTimeout>;
    const tick = () => {
      timeoutId = setTimeout(async () => {
        if (cancelled) return;
        await stepRef.current();
        if (!cancelled) tick();
      }, PLAY_INTERVAL_MS);
    };
    tick();
    return () => {
      cancelled = true;
      clearTimeout(timeoutId);
    };
  }, [playing]);

  return (
    <div style={{ fontFamily: "system-ui, sans-serif", color: "#e6e8ef", padding: 24 }}>
      <div style={{ display: "flex", gap: 16, alignItems: "flex-start" }}>
        <div style={{ display: "flex", flexDirection: "column", gap: 8, width: 170, flexShrink: 0 }}>
          <h1 style={{ fontSize: 20, marginTop: 0, marginBottom: 4, textAlign: "center" }}>Mantle Bloom</h1>
          <p style={{ opacity: 0.7, marginTop: 0, marginBottom: 16 }}>
            Physical World Builder
          </p>

          <button
            onClick={() => setShowGenerateDialog(true)}
            disabled={busy || animating}
            className={!summary ? "generate-flashing" : undefined}
            style={{ fontSize: 12 }}
          >
            Generate World
          </button>
          {busy && <ProgressBar label="Generating world" />}

          <button onClick={() => setShowStatsModal(true)} disabled={!summary} style={{ fontSize: 12 }}>
            📊 Stats
          </button>

          <button
            onClick={() => setShowControlsModal(true)}
            disabled={!summary || animating}
            style={{ fontSize: 12 }}
          >
            🎛️ Controls
          </button>

          <button onClick={() => setShowFileModal(true)} disabled={animating} style={{ fontSize: 12 }}>
            📁 File...
          </button>

          <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
            <legend style={{ fontSize: 11 }}>Time</legend>
            <label style={{ display: "block", marginBottom: 6 }}>
              Years per step
              <select
                value={stepYears}
                onChange={(e) => setStepYears(Number(e.target.value))}
                style={{ width: "100%", fontSize: 12 }}
              >
                {STEP_YEARS_OPTIONS.map((y) => (
                  <option key={y} value={y}>
                    {y.toLocaleString()}
                  </option>
                ))}
              </select>
            </label>
            {/* Play / Stop / Record -- the classic tape-deck order. Play advances the world
                continuously (one step every PLAY_INTERVAL_MS, sized by "Years per step"
                above) until Stop is pressed; Record opens AnimationModal to configure and
                start a video recording, which Stop also ends (cleanly -- see
                handleStopAnimation). Stop is only enabled while one of the two is running.
                Ending a recording, whether by Stop or by running its course, resolves
                handleStartAnimation's animateWorld() call and pops up SaveAnimationModal
                (below) to save or discard the finished MP4 -- there's no separate "ready"
                panel to go looking for elsewhere in the sidebar for that part. A fixed-frame
                recording does still get its own sidebar panel *while it runs* (frame X of Y,
                below); a "keep going until Stop" one (the AnimationModal default) instead
                relies on the map's live preview and the ticking "elapsed" readout for
                progress, since it has no fixed total to show a fraction of. */}
            <div style={{ display: "flex", gap: 6 }}>
              <button
                onClick={() => setPlaying(true)}
                disabled={busy || stepping || !summary || animating || playing}
                title="Play"
                aria-label="Play"
                // Green while clickable, same "color signals affordance" treatment as Record's
                // red (below) -- both go back to the default (undefined) color once disabled.
                style={{
                  flex: 1, fontSize: 15, lineHeight: 1, padding: "4px 0",
                  color: busy || stepping || !summary || animating || playing ? undefined : "#27c93f",
                }}
              >
                ▶
              </button>
              <button
                onClick={() => (animating ? handleStopAnimation() : setPlaying(false))}
                disabled={!playing && !animating}
                title="Stop"
                aria-label="Stop"
                style={{ flex: 1, fontSize: 15, lineHeight: 1, padding: "4px 0" }}
              >
                ⏹
              </button>
              <button
                onClick={() => setShowAnimationModal(true)}
                disabled={busy || !summary || animating || playing}
                title="Record animation"
                aria-label="Record animation"
                // Flashes red/grey while a recording is in flight (see index.css's
                // .record-flashing, which overrides the static `color` below for as long as
                // it's applied) so it reads as "recording" the way a camcorder's tally light
                // does, distinct from just being disabled like the other controls.
                className={animating ? "record-flashing" : undefined}
                style={{
                  flex: 1, fontSize: 15, lineHeight: 1, padding: "4px 0",
                  color: busy || !summary || animating || playing ? undefined : "#ff5f56",
                }}
              >
                ⏺
              </button>
            </div>
            {stepping && <ProgressBar label="Stepping world" style={{ marginTop: 6 }} />}
          </fieldset>

          <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
            <legend style={{ fontSize: 11 }}>Map View</legend>
            <select
              value={mapView}
              onChange={(e) => setMapView(e.target.value as MapView)}
              disabled={animating || !summary}
              style={{ width: "100%", marginBottom: 6, fontSize: 12 }}
            >
              <optgroup label="Maps">
                <option value="combined">Elevation &amp; Biome</option>
                <option value="elevation">Elevation</option>
                <option value="biome">Biome</option>
                <option value="temperature">Temperature</option>
                <option value="wind">Wind</option>
                <option value="oceanCurrents">Ocean currents</option>
                <option value="humidity">Humidity</option>
                <option value="precipitation">Precipitation</option>
                <option value="resources">Resources</option>
                <option value="soilQuality">Soil Quality</option>
              </optgroup>
              {DEBUG_UI && (
                <optgroup label="Debug &gt;">
                  <option value="platesDetail">Points</option>
                  <option value="speckle">Coastal dither (speckle)</option>
                  <option value="platesAndFaults">Plates &amp; Faults</option>
                  <option value="geomorph">Erosion &amp; Deposition</option>
                  <option value="elevReason">Last elevation change</option>
                  <option value="overlapAge">Plate overlap age</option>
                  <option value="nodeAge">Added/removed points</option>
                  <option value="plateInspector">Plate Inspector</option>
                  <option value="riverInspector">Rivers</option>
                  <option value="lakeInspector">Lake Inspector</option>
                </optgroup>
              )}
            </select>
            <select
              value={projection}
              onChange={(e) => setProjection(e.target.value as Projection)}
              disabled={animating || !summary}
              style={{ width: "100%", fontSize: 12 }}
            >
              <option value="behrmann">Behrmann (cylindrical equal-area)</option>
              <option value="eckert4">Eckert IV (pseudocylindrical equal-area)</option>
            </select>
            {/* position: relative anchors the popup below to this row specifically (not the
                whole fieldset), so it floats over the Re-center button rather than pushing it
                down -- same floating-card treatment the map's click-to-inspect probe popup uses
                (see the `probe &&` block above). */}
            <div style={{ position: "relative", marginTop: 6 }}>
              <div
                onClick={handleOpenCenterEdit}
                title={!summary || animating || editingCenter ? undefined : "Click to type a center coordinate"}
                style={{ opacity: 0.8, cursor: !summary || animating || editingCenter ? undefined : "pointer", userSelect: "none" }}
              >
                Center: {formatLatLon(centerLatLon.lat, centerLatLon.lon)}
              </div>
              {editingCenter && (
                <form
                  onSubmit={(e) => {
                    e.preventDefault();
                    handleCommitCenterEdit();
                  }}
                  onKeyDown={(e) => {
                    if (e.key === "Escape") handleCancelCenterEdit();
                  }}
                  style={{
                    position: "absolute",
                    top: "100%",
                    left: 0,
                    right: 0,
                    marginTop: 4,
                    zIndex: 30,
                    background: "#151a2e",
                    border: "1px solid #333",
                    borderRadius: 6,
                    padding: 8,
                    boxShadow: "0 2px 10px rgba(0, 0, 0, 0.5)",
                  }}
                >
                  <div style={{ display: "flex", gap: 4 }}>
                    <input
                      type="number"
                      inputMode="decimal"
                      step="any"
                      min={-90}
                      max={90}
                      value={centerEditValue.lat}
                      onChange={(e) => setCenterEditValue((v) => ({ ...v, lat: e.target.value }))}
                      placeholder="lat (-90 to 90)"
                      aria-label="Center latitude"
                      autoFocus
                      style={{ width: "50%", fontSize: 12 }}
                    />
                    <input
                      type="number"
                      inputMode="decimal"
                      step="any"
                      min={-180}
                      max={180}
                      value={centerEditValue.lon}
                      onChange={(e) => setCenterEditValue((v) => ({ ...v, lon: e.target.value }))}
                      placeholder="lon (-180 to 180)"
                      aria-label="Center longitude"
                      style={{ width: "50%", fontSize: 12 }}
                    />
                  </div>
                  {centerEditError && (
                    <div style={{ color: "#ff8080", fontSize: 10, marginTop: 2 }}>{centerEditError}</div>
                  )}
                  <div style={{ display: "flex", gap: 4, marginTop: 4 }}>
                    <button type="submit" style={{ flex: 1, fontSize: 12 }}>
                      Set
                    </button>
                    <button type="button" onClick={handleCancelCenterEdit} style={{ flex: 1, fontSize: 12 }}>
                      Cancel
                    </button>
                  </div>
                </form>
              )}
            </div>
            <button
              onClick={handleRecenter}
              disabled={busy || !summary || animating || isIdentityRotation(rotation) || editingCenter}
              style={{ width: "100%", marginTop: 6, fontSize: 12 }}
            >
              Re-center
            </button>
          </fieldset>

          {/* Only for a fixed-frame-count recording -- "frame 12 of 20" is meaningful
              progress, but a "keep going until Stop is pressed" run has no total to show a
              fraction of, so it relies on the main map's live preview and the ticking
              "elapsed" readout above instead (see the `animation` state comment). */}
          {animation && !animation.unbounded && (
            <fieldset style={{ border: "1px solid #3a4d8f", borderRadius: 6, padding: 8, fontSize: 12 }}>
              <legend style={{ fontSize: 11 }}>Animation</legend>
              <div style={{ height: 6, borderRadius: 3, background: "#2a3050", overflow: "hidden" }}>
                <div
                  style={{
                    height: "100%",
                    width: `${(animation.frame / animation.total) * 100}%`,
                    background: "#5b8cff",
                    transition: "width 0.2s linear",
                  }}
                />
              </div>
              <div style={{ fontSize: 11, opacity: 0.7, margin: "4px 0 0" }}>
                Rendering frame {animation.frame} of {animation.total}… the map is stepping in
                the background. Press ⏹ Stop above to finish early.
              </div>
            </fieldset>
          )}

          {(mapView === "plateInspector" || mapView === "platesAndFaults") && (
            <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
              <legend style={{ fontSize: 11 }}>Selected plate</legend>
              {selectedPlate ? (
                <div style={{ opacity: 0.9 }}>
                  <div>id: {selectedPlate.plate_id}</div>
                  <div>crust: {selectedPlate.crust_type}</div>
                  {selectedPlate.num_rows != null && <div>rows: {selectedPlate.num_rows}</div>}
                  <div>points: {selectedPlate.num_points}</div>
                  <div>age: {selectedPlate.age_steps} steps</div>
                  <div style={{ color: selectedPlate.at_max_rate ? "#e06c4b" : undefined }}>
                    speed: {selectedPlate.speed_cm_per_yr.toFixed(2)} cm/yr{selectedPlate.at_max_rate ? " (railed at MAX)" : ""}
                  </div>
                  {selectedPlate.euler_pole && (
                    <div>
                      euler pole: {selectedPlate.euler_pole.lat_deg.toFixed(0)}&deg;, {selectedPlate.euler_pole.lon_deg.toFixed(0)}&deg;
                    </div>
                  )}
                  {selectedPlate.median_elevation_m != null && (
                    <div
                      style={{
                        color:
                          selectedPlate.crust_type === "continental" && selectedPlate.submerged_fraction > 0.5
                            ? "#e06c4b"
                            : undefined,
                      }}
                    >
                      median elev: {selectedPlate.median_elevation_m.toFixed(0)} m &middot;{" "}
                      {(selectedPlate.submerged_fraction * 100).toFixed(0)}% submerged
                    </div>
                  )}
                  {selectedPlate.bounding_ellipse && (
                    <>
                      <div>diameter A: {selectedPlate.bounding_ellipse.diameter_a_km.toFixed(0)} km</div>
                      <div>diameter B: {selectedPlate.bounding_ellipse.diameter_b_km.toFixed(0)} km</div>
                    </>
                  )}
                  {selectedPlate.overlaps.length > 0 && (
                    <div style={{ marginTop: 4 }}>
                      overlaps:{" "}
                      {selectedPlate.overlaps
                        .filter((o) => o.fraction >= 0.01)
                        .map(
                          (o) =>
                            `#${o.plate_id} (${(o.fraction * 100).toFixed(0)}%` +
                            (o.since_years != null
                              ? `, since ${(o.since_years / 1e6).toFixed(0)} My`
                              : "") +
                            ")",
                        )
                        .join(", ") || "<1% only"}
                    </div>
                  )}
                  {selectedPlate.collisions.length > 0 && (
                    <div>
                      colliding:{" "}
                      {selectedPlate.collisions
                        .map((c) => `#${c.plate_id} (${(c.years / 1e6).toFixed(1)} My)`)
                        .join(", ")}
                    </div>
                  )}
                </div>
              ) : (
                <div style={{ opacity: 0.6 }}>Click a plate, or press Tab.</div>
              )}
            </fieldset>
          )}

          {mapView === "platesAndFaults" && (
            <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
              <legend style={{ fontSize: 11 }}>Overlays</legend>
              <label style={{ display: "flex", gap: 6, alignItems: "center" }}>
                <input
                  type="checkbox"
                  checked={showQuakesVolcanoes}
                  onChange={(e) => setShowQuakesVolcanoes(e.target.checked)}
                />
                Earthquakes &amp; volcanoes
              </label>
            </fieldset>
          )}

          {mapView === "riverInspector" && (
            <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
              <legend style={{ fontSize: 11 }}>Selected river</legend>
              {selectedRiver ? (
                <div style={{ opacity: 0.9 }}>
                  <div>id: {selectedRiver.river_id}</div>
                  <div>ends at: {selectedRiver.mouth_type}</div>
                  <div>flow rate: {selectedRiver.flow_rate.toFixed(1)}</div>
                  <div>speed: {selectedRiver.speed.toFixed(2)}</div>
                  <div>tributaries: {selectedRiver.num_tributaries}</div>
                  <div>nodes: {selectedRiver.num_nodes}</div>
                </div>
              ) : (
                <div style={{ opacity: 0.6 }}>
                  {riversData.length > 0 ? "Click a river, or press Tab." : "No rivers yet -- step the world forward."}
                </div>
              )}
            </fieldset>
          )}

          {mapView === "lakeInspector" && (
            <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
              <legend style={{ fontSize: 11 }}>Selected basin</legend>
              {selectedBasin ? (
                <div style={{ opacity: 0.9 }}>
                  <div>{selectedBasin.is_lake ? `lake id: ${selectedBasin.lake_id}` : "dry basin (no standing water)"}</div>
                  <div>nodes: {selectedBasin.member_count}</div>
                  <div>lowest point: {selectedBasin.floor_elevation_m.toFixed(0)} m</div>
                  {selectedBasin.is_lake && (
                    <>
                      <div>water level: {selectedBasin.water_elevation_m?.toFixed(0)} m</div>
                      <div>{selectedBasin.is_spilling ? "spilling over its outlet" : "not currently spilling"}</div>
                    </>
                  )}
                  <div>
                    exit (edge low point):{" "}
                    {selectedBasin.outlet_elevation_m != null ? `${selectedBasin.outlet_elevation_m.toFixed(0)} m` : "none known -- closed basin"}
                  </div>
                  <div style={{ marginTop: 6 }}>
                    inflowing rivers: {selectedBasin.inflow_rivers.length === 0 ? "none" : ""}
                  </div>
                  {selectedBasin.inflow_rivers.map((r, i) => (
                    <div key={i} style={{ paddingLeft: 8, opacity: 0.8 }}>
                      - flow {r.flow_rate.toFixed(1)}, {r.num_nodes} nodes
                    </div>
                  ))}
                  <div style={{ marginTop: 6 }}>
                    outflowing river:{" "}
                    {selectedBasin.outflow_river
                      ? ""
                      : selectedBasin.is_spilling
                        ? "not yet a real channel"
                        : "none -- not currently spilling"}
                  </div>
                  {selectedBasin.outflow_river && (
                    <div style={{ paddingLeft: 8, opacity: 0.8 }}>
                      - ends at {selectedBasin.outflow_river.mouth_type}, flow{" "}
                      {selectedBasin.outflow_river.flow_rate.toFixed(1)}, {selectedBasin.outflow_river.num_nodes} nodes
                    </div>
                  )}
                </div>
              ) : selectedBasinKind === "ocean" ? (
                <div style={{ opacity: 0.6 }}>That's open ocean.</div>
              ) : selectedBasinKind === "no_basin" ? (
                <div style={{ opacity: 0.6 }}>This point drains straight to the ocean -- no enclosed basin here.</div>
              ) : (
                <div style={{ opacity: 0.6 }}>
                  {lakesData.length > 0 ? "Click a lake or basin, or press Tab." : "Click any point on land to inspect its basin."}
                </div>
              )}
            </fieldset>
          )}

          {summary && (
            <div style={{ fontSize: 11, opacity: 0.8 }}>
              <div>seed: {summary.seed}</div>
              <div>plates: {summary.num_plates}</div>
              <div>elapsed: {(summary.elapsed_years / 1e6).toFixed(stepYears === 10_000 ? 2 : 1)} Myr</div>
              {animation && <div>{animation.frame} frames recorded</div>}
            </div>
          )}
          {error && <div style={{ color: "#ff8080", fontSize: 11 }}>{error}</div>}

          {DEBUG_UI && (
            <>
              <EventConsole events={summary?.events ?? []} />
              <div style={{ marginTop: 8 }}>
                <CornerNotchLogPanel entries={cornerNotchLog} enabled={debugDiagnostics} />
              </div>
            </>
          )}
        </div>

        <div>
          <div ref={mapWrapperRef} style={{ position: "relative", width: DISPLAY_WIDTH, height: DISPLAY_HEIGHT }}>
          {mapView === "plateInspector" ? (
            <PlateInspector
              plates={platesData}
              width={RENDER_WIDTH}
              height={RENDER_HEIGHT}
              displayWidth={DISPLAY_WIDTH}
              displayHeight={DISPLAY_HEIGHT}
              projection={projection}
              rotation={rotation}
              selectedPlateId={selectedPlateId}
              onSelectPlate={setSelectedPlateId}
              onRotationPreview={(latDeg, lonDeg) => setCenterLatLon({ lat: latDeg, lon: lonDeg })}
              onRotationCommitted={(newRotation) => setRotation(newRotation)}
              interactionDisabled={animating}
            />
          ) : mapView === "riverInspector" ? (
            <RiverInspector
              rivers={riversData}
              coastlineSegments={coastlineSegments}
              width={RENDER_WIDTH}
              height={RENDER_HEIGHT}
              displayWidth={DISPLAY_WIDTH}
              displayHeight={DISPLAY_HEIGHT}
              projection={projection}
              rotation={rotation}
              selectedRiverId={selectedRiverId}
              onSelectRiver={setSelectedRiverId}
              onRotationPreview={(latDeg, lonDeg) => setCenterLatLon({ lat: latDeg, lon: lonDeg })}
              onRotationCommitted={(newRotation) => setRotation(newRotation)}
              interactionDisabled={animating}
            />
          ) : mapView === "lakeInspector" ? (
            <LakeInspector
              lakes={lakesData}
              coastlineSegments={coastlineSegments}
              width={RENDER_WIDTH}
              height={RENDER_HEIGHT}
              displayWidth={DISPLAY_WIDTH}
              displayHeight={DISPLAY_HEIGHT}
              projection={projection}
              rotation={rotation}
              selectedBasin={selectedBasin}
              onSelect={(kind, basin) => {
                setSelectedBasinKind(kind);
                setSelectedBasin(basin);
              }}
              onRotationPreview={(latDeg, lonDeg) => setCenterLatLon({ lat: latDeg, lon: lonDeg })}
              onRotationCommitted={(newRotation) => setRotation(newRotation)}
              interactionDisabled={animating}
            />
          ) : mapView === "platesAndFaults" ? (
            <PlatesAndFaults
              plates={platesData}
              faults={faultsData}
              faultSystems={faultSystemsData}
              earthquakes={earthquakesData}
              volcanoes={volcanoesData}
              coastlineSegments={coastlineSegments}
              showQuakesVolcanoes={showQuakesVolcanoes}
              highlightedFaultKind={highlightedBiome ? faultKindForLegendLabel(highlightedBiome) : null}
              width={RENDER_WIDTH}
              height={RENDER_HEIGHT}
              displayWidth={DISPLAY_WIDTH}
              displayHeight={DISPLAY_HEIGHT}
              projection={projection}
              rotation={rotation}
              selectedPlateId={selectedPlateId}
              onSelectPlate={setSelectedPlateId}
              onRotationPreview={(latDeg, lonDeg) => setCenterLatLon({ lat: latDeg, lon: lonDeg })}
              onRotationCommitted={(newRotation) => setRotation(newRotation)}
              interactionDisabled={animating}
            />
          ) : (
            // tabIndex + onKeyDown: only the Points (platesDetail) view's ArrowLeft/Right and
            // Shift+ArrowLeft/Right do anything (see handlePointKeyDown), but the container is
            // always present/focusable, same as PlateInspector's own wrapper -- harmless on
            // every other view since handlePointKeyDown no-ops when mapView isn't platesDetail.
            <div ref={pointContainerRef} tabIndex={0} onKeyDown={handlePointKeyDown} style={{ outline: "none", display: "inline-block" }}>
              <MapCanvas
                imageBase64={renderData?.image_base64 ?? null}
                width={RENDER_WIDTH}
                height={RENDER_HEIGHT}
                displayWidth={DISPLAY_WIDTH}
                displayHeight={DISPLAY_HEIGHT}
                projection={projection}
                rotation={rotation}
                onRotationPreview={(latDeg, lonDeg) => setCenterLatLon({ lat: latDeg, lon: lonDeg })}
                onRotationCommitted={(newRotation) => setRotation(newRotation)}
                highlightTarget={highlightTarget}
                onProbe={
                  mapView === "nodeAge"
                    ? handleNodeProbe
                    : mapView === "platesDetail"
                      ? handlePointProbe
                      : mapView === "combined" || mapView === "elevation" || mapView === "biome"
                        ? handleProbe
                        : undefined
                }
                highlightLine={
                  mapView === "platesDetail" && pointProbe?.result
                    ? { pointsXyz: pointProbe.result.line_points_xyz, selectedIndex: pointProbe.result.point.index }
                    : null
                }
                alphaEncodedIds={mapView === "combined" || mapView === "biome"}
                interactionDisabled={animating}
                coastlineSegments={coastlineSegments}
              />
            </div>
          )}
          {probe && (
            <div
              style={{
                position: "absolute",
                left: Math.max(4, Math.min(probe.displayX + 12, DISPLAY_WIDTH - 186)),
                top: Math.max(4, Math.min(probe.displayY + 12, DISPLAY_HEIGHT - 150)),
                width: 174,
                background: "#151a2e",
                border: "1px solid #333",
                borderRadius: 6,
                padding: "8px 10px",
                fontSize: 11,
                lineHeight: 1.6,
                boxShadow: "0 2px 10px rgba(0, 0, 0, 0.5)",
                zIndex: 10,
              }}
            >
              <button
                type="button"
                title="Close"
                onClick={() => handleProbe(null)}
                style={{
                  position: "absolute", top: 3, right: 4, width: 18, height: 18, padding: 0,
                  border: "none", background: "transparent", color: "#999", cursor: "pointer",
                  fontSize: 14, lineHeight: "18px",
                }}
              >
                ×
              </button>
              <div style={{ opacity: 0.6, marginBottom: 4 }}>{formatLatLon(probe.latDeg, probe.lonDeg)}</div>
              {probe.status === "loading" && <div style={{ opacity: 0.7 }}>Sampling…</div>}
              {probe.status === "error" && <div style={{ color: "#ff8080" }}>Couldn’t sample this point.</div>}
              {probe.status === "ok" && probe.sample && (
                <div style={{ display: "grid", gridTemplateColumns: "auto 1fr", columnGap: 8, rowGap: 1 }}>
                  <span style={{ opacity: 0.55 }}>Elevation</span>
                  <span>{Math.round(probe.sample.elevation_m).toLocaleString()} m</span>
                  <span style={{ opacity: 0.55 }}>Biome</span>
                  <span>{probe.sample.biome}</span>
                  <span style={{ opacity: 0.55 }}>Precip.</span>
                  <span>{Math.round(probe.sample.precipitation_mm).toLocaleString()} mm</span>
                  <span style={{ opacity: 0.55 }}>Temp.</span>
                  <span>{probe.sample.temperature_c.toFixed(1)} °C</span>
                  <span style={{ opacity: 0.55 }}>Plate</span>
                  <span>{probe.sample.plate_id ?? "—"}</span>
                </div>
              )}
            </div>
          )}
          {nodeProbe && (
            <div
              style={{
                position: "absolute",
                left: Math.max(4, Math.min(nodeProbe.displayX + 12, DISPLAY_WIDTH - 200)),
                top: Math.max(4, Math.min(nodeProbe.displayY + 12, DISPLAY_HEIGHT - 170)),
                width: 188,
                background: "#151a2e",
                border: "1px solid #333",
                borderRadius: 6,
                padding: "8px 10px",
                fontSize: 11,
                lineHeight: 1.6,
                boxShadow: "0 2px 10px rgba(0, 0, 0, 0.5)",
                zIndex: 10,
              }}
            >
              <button
                type="button"
                title="Close"
                onClick={() => handleNodeProbe(null)}
                style={{
                  position: "absolute", top: 3, right: 4, width: 18, height: 18, padding: 0,
                  border: "none", background: "transparent", color: "#999", cursor: "pointer",
                  fontSize: 14, lineHeight: "18px",
                }}
              >
                ×
              </button>
              <div style={{ opacity: 0.6, marginBottom: 4 }}>{formatLatLon(nodeProbe.latDeg, nodeProbe.lonDeg)}</div>
              {nodeProbe.status === "loading" && <div style={{ opacity: 0.7 }}>Sampling…</div>}
              {nodeProbe.status === "error" && <div style={{ color: "#ff8080" }}>Couldn’t sample this point.</div>}
              {nodeProbe.status === "ok" && nodeProbe.result && (
                <div>
                  {nodeProbe.result.node && (
                    <div style={{ display: "grid", gridTemplateColumns: "auto 1fr", columnGap: 8, rowGap: 1 }}>
                      <span style={{ opacity: 0.55 }}>Plate</span>
                      <span>{nodeProbe.result.node.plate_id}</span>
                      <span style={{ opacity: 0.55 }}>Phi</span>
                      <span>{nodeProbe.result.node.phi.toFixed(4)} rad</span>
                      <span style={{ opacity: 0.55 }}>Theta</span>
                      <span>{nodeProbe.result.node.theta.toFixed(4)} rad</span>
                      <span style={{ opacity: 0.55 }}>Elevation</span>
                      <span>{Math.round(nodeProbe.result.node.elevation_m).toLocaleString()} m</span>
                      <span style={{ opacity: 0.55 }}>Created</span>
                      <span>
                        {nodeProbe.result.node.node_created_years < 0
                          ? "predates tracking"
                          : `${(nodeProbe.result.node.node_created_years / 1e6).toFixed(1)} My`}
                      </span>
                    </div>
                  )}
                  {nodeProbe.result.removed && (
                    <div
                      style={{
                        marginTop: nodeProbe.result.node ? 6 : 0,
                        paddingTop: nodeProbe.result.node ? 6 : 0,
                        borderTop: nodeProbe.result.node ? "1px solid #333" : "none",
                        display: "grid",
                        gridTemplateColumns: "auto 1fr",
                        columnGap: 8,
                        rowGap: 1,
                      }}
                    >
                      <span style={{ opacity: 0.55 }}>Removed</span>
                      <span>{(nodeProbe.result.removed.removed_years / 1e6).toFixed(1)} My ago</span>
                      <span style={{ opacity: 0.55 }}>From plate</span>
                      <span>{nodeProbe.result.removed.plate_id}</span>
                    </div>
                  )}
                  {!nodeProbe.result.node && !nodeProbe.result.removed && (
                    <div style={{ opacity: 0.7 }}>No recent activity here.</div>
                  )}
                </div>
              )}
            </div>
          )}
          {pointProbe && (
            <div
              style={{
                position: "absolute",
                left: Math.max(4, Math.min(pointProbe.displayX + 12, DISPLAY_WIDTH - 200)),
                top: Math.max(4, Math.min(pointProbe.displayY + 12, DISPLAY_HEIGHT - 190)),
                width: 188,
                background: "#151a2e",
                border: "1px solid #333",
                borderRadius: 6,
                padding: "8px 10px",
                fontSize: 11,
                lineHeight: 1.6,
                boxShadow: "0 2px 10px rgba(0, 0, 0, 0.5)",
                zIndex: 10,
              }}
            >
              <button
                type="button"
                title="Close"
                onClick={() => handlePointProbe(null)}
                style={{
                  position: "absolute", top: 3, right: 4, width: 18, height: 18, padding: 0,
                  border: "none", background: "transparent", color: "#999", cursor: "pointer",
                  fontSize: 14, lineHeight: "18px",
                }}
              >
                ×
              </button>
              <div style={{ opacity: 0.6, marginBottom: 4 }}>{formatLatLon(pointProbe.latDeg, pointProbe.lonDeg)}</div>
              {pointProbe.status === "loading" && <div style={{ opacity: 0.7 }}>Sampling…</div>}
              {pointProbe.status === "error" && <div style={{ color: "#ff8080" }}>Couldn’t sample this point.</div>}
              {pointProbe.status === "ok" && pointProbe.result && (
                <div>
                  <div style={{ display: "grid", gridTemplateColumns: "auto 1fr", columnGap: 8, rowGap: 1 }}>
                    <span style={{ opacity: 0.55 }}>Plate</span>
                    <span>{pointProbe.result.plate_id}</span>
                    <span style={{ opacity: 0.55 }}>Phi</span>
                    <span>{pointProbe.result.point.phi.toFixed(4)} rad</span>
                    <span style={{ opacity: 0.55 }}>Theta</span>
                    <span>{pointProbe.result.point.theta.toFixed(4)} rad</span>
                    <span style={{ opacity: 0.55 }}>Elevation</span>
                    <span>{Math.round(pointProbe.result.point.elevation_m).toLocaleString()} m</span>
                  </div>
                  <div
                    style={{
                      marginTop: 6, paddingTop: 6, borderTop: "1px solid #333",
                      display: "grid", gridTemplateColumns: "auto 1fr", columnGap: 8, rowGap: 1,
                    }}
                  >
                    <span style={{ opacity: 0.55 }}>Point</span>
                    <span>{pointProbe.result.point.index + 1} of {pointProbe.result.line.num_points}</span>
                    <span style={{ opacity: 0.55 }}>Line</span>
                    <span>{pointProbe.result.line.line_index + 1} of {pointProbe.result.line.num_lines}</span>
                  </div>
                </div>
              )}
            </div>
          )}
          <MeasureOverlay
            projection={projection}
            width={RENDER_WIDTH}
            height={RENDER_HEIGHT}
            displayWidth={DISPLAY_WIDTH}
            displayHeight={DISPLAY_HEIGHT}
            disabled={animating}
          />
          </div>
          <Legend
            mapView={mapView}
            projection={projection}
            width={RENDER_WIDTH}
            height={RENDER_HEIGHT}
            displayWidth={DISPLAY_WIDTH}
            displayHeight={DISPLAY_HEIGHT}
            mapWrapperRef={mapWrapperRef}
            hasWorld={!!summary}
            highlightedBiome={highlightedBiome}
            onBiomeClick={(label) => setHighlightedBiome((cur) => (cur === label ? null : label))}
            showMountains={showMountains}
            onToggleMountains={() => setShowMountains((cur) => !cur)}
            showPlainsPlateaus={showPlainsPlateaus}
            onTogglePlainsPlateaus={() => setShowPlainsPlateaus((cur) => !cur)}
          />
          <p style={{ fontSize: 11, opacity: 0.6, marginTop: 4 }}>
            {mapView === "plateInspector"
              ? "Click a plate to select it. Tab / Shift+Tab cycles plates. Press and hold, then drag to rotate."
              : mapView === "riverInspector"
                ? "Click a river to select it. Tab / Shift+Tab cycles rivers. Press and hold, then drag to rotate."
                : mapView === "lakeInspector"
                  ? "Click a lake or any point on land to inspect its basin. Tab / Shift+Tab cycles lakes. Press and hold, then drag to rotate."
                  : mapView === "platesAndFaults"
                  ? "Click a plate to select it (its fault strands emphasise, and its Euler pole + a speed-scaled motion arc appear). Tab / Shift+Tab cycles plates. Click a fault type in the legend to isolate that regime. Toggle the earthquake & volcano overlay in the sidebar. Press and hold, then drag to rotate."
                  : mapView === "combined" || mapView === "elevation" || mapView === "biome"
                    ? "Click any point for its elevation, biome, precipitation, temperature, and plate. Press and hold, then drag to rotate."
                    : mapView === "platesDetail"
                      ? "Click a point for its phi/theta and its ElevationLine. Left/Right steps along the line, Shift+Left/Right steps to the next/previous line. Press and hold, then drag to rotate."
                      : "Press and hold, then drag the map to rotate it."}
          </p>
        </div>
      </div>

      {showGenerateDialog && (
        <div
          style={{
            position: "fixed",
            inset: 0,
            background: "rgba(0, 0, 0, 0.6)",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
          }}
        >
          <div
            style={{
              background: "#151a2e",
              border: "1px solid #333",
              borderRadius: 8,
              padding: 20,
              minWidth: 280,
              width: 380,
              maxWidth: "66vw",
              maxHeight: "85vh",
              overflowY: "auto",
            }}
          >
            <h2 style={{ fontSize: 16, marginTop: 0, marginBottom: 12 }}>Generate World</h2>

            <div style={{ display: "flex", marginBottom: 16, borderBottom: "1px solid #333" }}>
              {(DEBUG_UI
                ? (["random", "human", "premade", "debug"] as const)
                : (["random", "human", "premade"] as const)
              ).map((mode) => (
                <button
                  key={mode}
                  type="button"
                  onClick={() => {
                    // Only reset on an actual mode change -- re-clicking the already-active tab
                    // must not discard a manual slider adjustment. "debug" is excluded from the
                    // reset entirely: generateDebugWorld never reads voronoi_points at all (see
                    // handleGenerate), so resetting it there would just be a misleading no-op
                    // value shown in Advanced settings.
                    if (mode !== generateMode) {
                      setGenerateMode(mode);
                      if (mode !== "debug") {
                        // Reset to the mode-appropriate default (see DEFAULT_VORONOI_POINTS_RANDOM/
                        // _SKETCH's own comments) so switching tabs doesn't leave a value picked
                        // for a different mode's fidelity needs -- a manual adjustment within a
                        // mode is kept until the next tab switch.
                        setVoronoiPoints(mode === "random" ? DEFAULT_VORONOI_POINTS_RANDOM : DEFAULT_VORONOI_POINTS_SKETCH);
                      }
                    }
                  }}
                  style={{
                    flex: 1,
                    padding: "6px 0",
                    fontSize: 12,
                    border: "none",
                    background: "none",
                    cursor: "pointer",
                    color: generateMode === mode ? "#e6e8ef" : "#8b8fa3",
                    borderBottom: generateMode === mode ? "2px solid #5b8cff" : "2px solid transparent",
                  }}
                >
                  {mode === "random"
                    ? "Random"
                    : mode === "human"
                      ? "Human-made"
                      : mode === "premade"
                        ? "Premade worlds"
                        : "Debugging Worlds"}
                </button>
              ))}
            </div>

            {generateMode === "human" && (
              <div style={{ marginBottom: 16 }}>
                {sketchImageDataUrl ? (
                  <div>
                    <img
                      src={sketchImageDataUrl}
                      alt="Drawn/loaded coastline preview"
                      style={{ width: "100%", borderRadius: 4, border: "1px solid #333", display: "block" }}
                    />
                    <div style={{ display: "flex", gap: 6, marginTop: 6 }}>
                      <button type="button" onClick={() => setShowSketchEditor(true)} disabled={busy} style={{ flex: 1, fontSize: 12 }}>
                        Redraw
                      </button>
                      <button type="button" onClick={() => loadImageInputRef.current?.click()} disabled={busy} style={{ flex: 1, fontSize: 12 }}>
                        Load a different image
                      </button>
                      <button type="button" onClick={() => setSketchImageDataUrl(null)} disabled={busy} style={{ fontSize: 12 }}>
                        Clear
                      </button>
                    </div>
                  </div>
                ) : (
                  <div style={{ display: "flex", gap: 8 }}>
                    <button type="button" onClick={() => setShowSketchEditor(true)} disabled={busy} style={{ flex: 1, fontSize: 12 }}>
                      ✏️ Draw a map
                    </button>
                    <button type="button" onClick={() => loadImageInputRef.current?.click()} disabled={busy} style={{ flex: 1, fontSize: 12 }}>
                      📂 Load an image
                    </button>
                  </div>
                )}
                <input
                  ref={loadImageInputRef}
                  type="file"
                  accept="image/*"
                  style={{ display: "none" }}
                  onChange={handleLoadSketchFile}
                />
                <div style={{ fontSize: 11, color: "#999", marginTop: 6 }}>
                  Land and sea come from this coastline instead of "Initial land" below
                  (Advanced settings) -- draw, or load, a mostly-white image with the coastline
                  as a dark outline; the map's corners are assumed ocean. "Continental plates"
                  still controls how many plates a large landmass gets split across, the way
                  real continents straddle more than one.
                </div>
              </div>
            )}

            {generateMode === "debug" && (
              <div style={{ marginBottom: 16 }}>
                <label style={{ display: "block" }}>
                  Scenario
                  <select
                    value={debugScenario}
                    onChange={(e) => setDebugScenario(e.target.value)}
                    style={{ width: "100%", marginTop: 4 }}
                  >
                    {debugScenarios.map((s) => (
                      <option key={s.name} value={s.name}>
                        {s.label}
                      </option>
                    ))}
                  </select>
                </label>
                <div style={{ fontSize: 11, color: "#999", marginTop: 6 }}>
                  Tiny, low-resolution plate configurations with pinned (scripted, not
                  torque-driven) motion, for fast iteration on the plate-boundary gap-filling
                  problem -- see docs/debugging.md. Debug diagnostics (the corner-notch decision
                  log) start on automatically.
                </div>
              </div>
            )}

            {generateMode === "premade" && (
              <div style={{ marginBottom: 16 }}>
                <div style={{ display: "flex", gap: 8 }}>
                  {PREMADE_WORLDS.map((world) => (
                    <button
                      key={world.id}
                      type="button"
                      onClick={() => {
                        setSketchImageDataUrl(world.dataUrl);
                        setSeed(world.seed);
                        setPremadeWorldId(world.backendId);
                      }}
                      disabled={busy}
                      style={{
                        flex: 1,
                        padding: 4,
                        display: "flex",
                        flexDirection: "column",
                        alignItems: "center",
                        gap: 4,
                        cursor: "pointer",
                        background: "none",
                        border: sketchImageDataUrl === world.dataUrl ? "2px solid #5b8cff" : "2px solid #333",
                        borderRadius: 4,
                      }}
                    >
                      <img
                        src={world.dataUrl}
                        alt={`${world.label} coastline preview`}
                        style={{ width: "100%", borderRadius: 2, background: "#fff", display: "block" }}
                      />
                      <span style={{ fontSize: 11, color: "#e6e8ef", textAlign: "center" }}>{world.label}</span>
                    </button>
                  ))}
                </div>
                <div style={{ fontSize: 11, color: "#999", marginTop: 6 }}>
                  Same coastline-driven generation as "Human-made" above, just starting from one
                  of these built-in maps instead of a drawn or loaded one. "Present-day Earth"
                  also seeds real-world mountain ranges and rivers, not just coastlines.
                  "Pangaea" reassembles real continent outlines into the classic supercontinent
                  fit; "Dragons &amp; Zombie World" is an original interpretation of a
                  well-known fantasy world's geography, not a reproduction of any map.
                </div>
              </div>
            )}

            <label style={{ display: "block", marginBottom: 16 }}>
              Seed
              <div style={{ display: "flex", gap: 6 }}>
                <input
                  type="number"
                  value={seed}
                  onChange={(e) => setSeed(Number(e.target.value))}
                  disabled={generateMode === "premade"}
                  style={{ flex: 1 }}
                />
                <button
                  type="button"
                  title="Randomize seed"
                  aria-label="Randomize seed"
                  onClick={() => setSeed(randomSeed())}
                  disabled={generateMode === "premade"}
                >
                  🎲
                </button>
              </div>
              {generateMode === "premade" && (
                <div style={{ fontSize: 11, color: "#999", marginTop: 4 }}>
                  Fixed for premade worlds -- each one is tuned to this seed so it always
                  generates looking like its preview above, not a random plate layout.
                </div>
              )}
            </label>

            {generateMode !== "debug" && (
              <label style={{ display: "block", marginBottom: 16 }}>
                Detail
                <select
                  value={detail}
                  onChange={(e) => setDetail(Number(e.target.value))}
                  style={{ width: "100%", marginTop: 4 }}
                >
                  {DETAIL_CHOICES.map((d) => (
                    <option key={d.value} value={d.value}>
                      {d.label}
                    </option>
                  ))}
                </select>
                <div style={{ fontSize: 11, color: "#999", marginTop: 4 }}>
                  Elevation point density and climate & biome resolution together. Higher is
                  sharper -- less pixelated Temperature/Wind/Currents/Humidity/Precipitation/
                  Biome/Elevation &amp; Biome/Resources/Soil Quality maps and more elevation-line nodes -- but
                  simulation steps and rendering both run slower. Lower runs faster but coarser.
                </div>
              </label>
            )}

            <div style={{ display: "flex", gap: 8, justifyContent: "space-between" }}>
              {generateMode !== "debug" && (
                <button type="button" onClick={() => setShowAdvancedSettings(true)} disabled={busy}>
                  Advanced settings
                </button>
              )}
              <div style={{ display: "flex", gap: 8, marginLeft: generateMode === "debug" ? "auto" : undefined }}>
                <button onClick={() => setShowGenerateDialog(false)} disabled={busy}>
                  Cancel
                </button>
                <button
                  onClick={handleGenerate}
                  disabled={
                    busy ||
                    ((generateMode === "human" || generateMode === "premade") && !sketchImageDataUrl) ||
                    (generateMode === "debug" && !debugScenario)
                  }
                >
                  Generate
                </button>
              </div>
            </div>
          </div>
        </div>
      )}

      {showAdvancedSettings && (
        <AdvancedSettingsModal
          landPercent={landPercent}
          continentalPercent={continentalPercent}
          autoPlates={autoPlates}
          numPlates={numPlates}
          minPlates={MIN_PLATES}
          maxPlates={MAX_PLATES}
          voronoiPoints={voronoiPoints}
          minVoronoiPoints={MIN_VORONOI_POINTS}
          maxVoronoiPoints={MAX_VORONOI_POINTS}
          effectivePlateCount={autoPlates ? DEFAULT_PLATES : numPlates}
          generateMode={generateMode}
          sketchImageDataUrl={sketchImageDataUrl}
          seed={seed}
          axialTiltDeg={axialTiltDeg}
          initialSoilMaturityPercent={initialSoilMaturityPercent}
          fluidDensity={fluidDensity}
          fluidDensityChoices={FLUID_DETAIL_CHOICES}
          onLandPercentChange={setLandPercent}
          onContinentalPercentChange={setContinentalPercent}
          onAutoPlatesChange={setAutoPlates}
          onNumPlatesChange={setNumPlates}
          onVoronoiPointsChange={setVoronoiPoints}
          onAxialTiltDegChange={setAxialTiltDeg}
          onInitialSoilMaturityPercentChange={setInitialSoilMaturityPercent}
          onFluidDensityChange={setFluidDensity}
          onClose={() => setShowAdvancedSettings(false)}
        />
      )}

      {showStatsModal && <StatsModal stats={stats} history={statsHistory} onClose={() => setShowStatsModal(false)} />}

      {showControlsModal && (
        <ControlsModal
          seaLevelM={seaLevelM}
          solarMultiplier={solarMultiplier}
          iceAgePeriodYears={iceAgePeriodYears}
          simulatePlateMovement={simulatePlateMovement}
          simulateClimateBiomes={simulateClimateBiomes}
          windModel={windModel}
          faultDeformationMode={faultDeformationMode}
          debugDiagnostics={debugDiagnostics}
          tuning={tuning}
          onSeaLevelChange={handleSeaLevelChange}
          onSolarMultiplierChange={handleSolarMultiplierChange}
          onIceAgePeriodChange={handleIceAgePeriodChange}
          onSimulatePlateMovementChange={handleSimulatePlateMovementChange}
          onSimulateClimateBiomesChange={handleSimulateClimateBiomesChange}
          onWindModelChange={handleWindModelChange}
          onFaultDeformationModeChange={handleFaultDeformationModeChange}
          onDebugDiagnosticsChange={handleDebugDiagnosticsChange}
          onTuningChange={handleTuningChange}
          onTuningReset={handleTuningReset}
          onClose={() => setShowControlsModal(false)}
        />
      )}

      {showSketchEditor && (
        <SketchEditor
          initialImageDataUrl={sketchImageDataUrl}
          onDone={(dataUrl) => {
            setSketchImageDataUrl(dataUrl);
            setShowSketchEditor(false);
          }}
          onCancel={() => setShowSketchEditor(false)}
        />
      )}

      {showFileModal && (
        <FileModal
          hasWorld={!!summary}
          seed={summary?.seed ?? null}
          elapsedYears={summary?.elapsed_years ?? null}
          mapView={mapView}
          mapWrapperRef={mapWrapperRef}
          onClose={() => setShowFileModal(false)}
          onWorldReplaced={handleWorldReplaced}
        />
      )}

      {showAnimationModal && (
        <AnimationModal
          hasWorld={!!summary}
          stepYears={stepYears}
          mapView={mapView}
          onClose={() => setShowAnimationModal(false)}
          onStartAnimation={handleStartAnimation}
        />
      )}

      {animationResult && (
        <SaveAnimationModal
          result={animationResult}
          onSave={handleSaveAnimation}
          onDiscard={handleDiscardAnimation}
        />
      )}
    </div>
  );
}
