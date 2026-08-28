import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import "./index.css";
import {
  fetchLakes, fetchPlates, fetchRivers, fetchStats, fetchWorldSummary, generateWorld, renderWorld, stepWorld, updateControls,
} from "./api";
import type {
  LakeAtResponse, LakeSummary, MapView, PlateSummary, Projection, RenderResponse, RiverSummary, Segment, WorldStats, WorldSummary,
} from "./api";
import MapCanvas from "./MapCanvas";
import PlateInspector from "./PlateInspector";
import RiverInspector from "./RiverInspector";
import LakeInspector from "./LakeInspector";
import EventConsole from "./EventConsole";
import StatsModal from "./StatsModal";
import ControlsModal from "./ControlsModal";
import AdvancedSettingsModal from "./AdvancedSettingsModal";
import FileModal from "./FileModal";
import Legend from "./Legend";
import { highlightTargetFor } from "./legendData";
import { centerOfRotation, IDENTITY_ROTATION } from "./rotation";
import type { Mat3 } from "./rotation";
import { getCookie, setCookie } from "./cookies";

// The map is displayed at DISPLAY_WIDTH x DISPLAY_HEIGHT (CSS pixels, unchanged from
// before), but the image requested from the server is RENDER_SCALE times bigger -- a
// sharper, retina-style render at the same on-screen size, rather than a bigger map.
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
// many." "Low" is the coarsest, fastest option; the default ("Very High") is the finest.
const DETAIL_CHOICES: { value: number; label: string }[] = [
  { value: 4, label: "Very High" },
  { value: 2, label: "High" },
  { value: 1, label: "Medium" },
  { value: 0.5, label: "Low" },
];
const DEFAULT_DETAIL = 4;
// The Advanced-settings dialog's own "Fluid dynamics resolution" choice (see backend app/
// world.py's World.fluid_density) -- same shape as DETAIL_CHOICES but capped at "High": the
// atmospheric wind solver runs every step (see docs/simulation-model.md#ocean-atmospheric-
// fluid-dynamics), so there's no "only pay for Very High when you opt in" case left to justify
// offering it, matching backend app/climate.py's own FLUID_DENSITY_CHOICES.
const FLUID_DETAIL_CHOICES: { value: number; label: string }[] = [
  { value: 2, label: "High" },
  { value: 1, label: "Medium" },
  { value: 0.5, label: "Low" },
];
const DEFAULT_FLUID_DETAIL = 2;
// Matching backend app/plates.py's MIN_AUTO_PLATES/MAX_AUTO_PLATES -- the same range the
// world's own "Auto" (seed-based) plate count is drawn from, so an explicit slider value
// always lands somewhere the auto behavior could plausibly have picked too.
const MIN_PLATES = 8;
const MAX_PLATES = 20;
const DEFAULT_PLATES = 14;
// Matching backend app/world.py's World.sea_level_m/World.solar_multiplier defaults.
const DEFAULT_SEA_LEVEL_M = 0;
const DEFAULT_SOLAR_MULTIPLIER = 1;
// Matching backend app/world.py's World.simulate_plate_movement/World.simulate_climate_biomes
// defaults -- both on, i.e. a normal full simulation.
const DEFAULT_SIMULATE_PLATE_MOVEMENT = true;
const DEFAULT_SIMULATE_CLIMATE_BIOMES = true;
// Matching backend app/world.py's World.wind_model default -- the shallow-water CFD solve.
const DEFAULT_WIND_MODEL = "cfd";

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
  "elevation", "plates", "platesDetail", "temperature", "wind", "oceanCurrents", "humidity", "precipitation", "biome", "combined",
  "resources", "soilQuality", "plateInspector", "riverInspector", "lakeInspector",
]);
const PROJECTION_CHOICES = new Set<Projection>(["behrmann", "eckert4"]);

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
  const [showAdvancedSettings, setShowAdvancedSettings] = useState(false);
  const [seed, setSeed] = useState(randomSeed());
  const [continentalPercent, setContinentalPercent] = useState(DEFAULT_CONTINENTAL_PERCENT);
  const [landPercent, setLandPercent] = useState(DEFAULT_LAND_PERCENT);
  const [axialTiltDeg, setAxialTiltDeg] = useState(DEFAULT_AXIAL_TILT_DEG);
  const [detail, setDetail] = useState(DEFAULT_DETAIL);
  // The Advanced-settings dialog's own "Fluid dynamics resolution" choice -- same
  // DETAIL_CHOICES set as `detail` above, but a separate dial: unlike node_density/
  // climate_density (merged into `detail`), this only affects Ocean/Atmospheric Fluid
  // Dynamics's own grid (see backend app/world.py's World.fluid_density), not plate/climate
  // resolution, so it's worth letting the user pick independently rather than folding it into
  // `detail` too. Defaults to DEFAULT_FLUID_DETAIL ("High"), matching the backend's own
  // FLUID_DENSITY_CHOICES cap -- see that constant's own comment for why it's lower than
  // DETAIL_CHOICES' own "Very High".
  const [fluidDensity, setFluidDensity] = useState(DEFAULT_FLUID_DETAIL);
  const [initialSoilMaturityPercent, setInitialSoilMaturityPercent] = useState(DEFAULT_INITIAL_SOIL_MATURITY_PERCENT);
  const [autoPlates, setAutoPlates] = useState(true);
  const [numPlates, setNumPlates] = useState(DEFAULT_PLATES);

  const [stepYears, setStepYears] = useState(STEP_YEARS_OPTIONS[1]);
  const [projection, setProjection] = useState<Projection>(initialView?.projection ?? "eckert4");
  const [mapView, setMapView] = useState<MapView>(initialView?.mapView ?? "elevation");
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
  // Legend-click-to-highlight (see Legend.tsx/MapCanvas.tsx) -- only ever meaningful on the
  // Biome and Combined views (the only legends whose swatches are clickable), so it's cleared
  // any time the view changes away from both rather than silently carrying a stale selection
  // into a view whose legend can't reflect or clear it.
  const [highlightedBiome, setHighlightedBiome] = useState<string | null>(null);
  useEffect(() => {
    if (mapView !== "biome" && mapView !== "combined") setHighlightedBiome(null);
  }, [mapView]);
  // Memoized so its identity only changes with the selection itself, not on every render --
  // MapCanvas.tsx's highlight-toggle effect is keyed directly on this object's identity.
  const highlightTarget = useMemo(
    () => (highlightedBiome ? highlightTargetFor(mapView, highlightedBiome) : null),
    [mapView, highlightedBiome],
  );
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
  // Stats panel data (see StatsModal.tsx) -- `stats` is the latest snapshot, `statsHistory`
  // accumulates one entry per generate/step (deduped by elapsed_years) for the panel's graph
  // tabs, built entirely client-side since the backend endpoint itself is stateless (see
  // backend app/stats.py). Recorded continuously, not just while the modal is open, so
  // opening it later still shows the full history since the world was generated.
  const [stats, setStats] = useState<WorldStats | null>(null);
  const [statsHistory, setStatsHistory] = useState<WorldStats[]>([]);
  const [showStatsModal, setShowStatsModal] = useState(false);
  // Live world controls (see ControlsModal.tsx and backend app/world.py's World.sea_level_m/
  // World.solar_multiplier) -- reset to their defaults on every fresh Generate, same as
  // every other generation-time value here, even though these two are adjustable afterward
  // too (unlike the others) via /world/controls.
  const [seaLevelM, setSeaLevelM] = useState(DEFAULT_SEA_LEVEL_M);
  const [solarMultiplier, setSolarMultiplier] = useState(DEFAULT_SOLAR_MULTIPLIER);
  // Same live-adjustable-via-Controls pattern as seaLevelM/solarMultiplier above -- lets the
  // user run plate tectonics only, climate & biomes only, or (the default) both together. See
  // backend app/world.py's World.simulate_plate_movement/World.simulate_climate_biomes.
  const [simulatePlateMovement, setSimulatePlateMovement] = useState(DEFAULT_SIMULATE_PLATE_MOVEMENT);
  const [simulateClimateBiomes, setSimulateClimateBiomes] = useState(DEFAULT_SIMULATE_CLIMATE_BIOMES);
  // "cfd" (shallow-water solve) or "diagnostic" (fast closed-form ABL wind) -- see backend
  // app/world.py's World.wind_model. Live-adjustable via Controls like the toggles above.
  const [windModel, setWindModel] = useState(DEFAULT_WIND_MODEL);
  const [showControlsModal, setShowControlsModal] = useState(false);
  const [showFileModal, setShowFileModal] = useState(false);
  // The div wrapping whichever map view component is currently mounted -- see
  // FileModal.tsx's own "Save Image" comment for why it reads the live <canvas> straight
  // out of this DOM node rather than a ref threaded through four separate components.
  const mapWrapperRef = useRef<HTMLDivElement>(null);
  const [playing, setPlaying] = useState(false);
  const [busy, setBusy] = useState(false);
  const [stepping, setStepping] = useState(false);
  const [error, setError] = useState<string | null>(null);

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
    if (view === "plateInspector" || view === "riverInspector" || view === "lakeInspector") return; // none of these use renderData -- see below
    const requestId = ++renderRequestIdRef.current;
    try {
      const data = await renderWorld(proj, view, RENDER_WIDTH, RENDER_HEIGHT, viewRotation);
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

  // Stats are a secondary/best-effort feature -- a failed fetch here (e.g. a transient
  // network blip) shouldn't surface as the main error line or block generate/step, unlike
  // refresh/refreshPlates above which are core to the map actually updating.
  const recordStats = useCallback(async () => {
    try {
      const s = await fetchStats();
      setStats(s);
      setStatsHistory((prev) =>
        prev.length > 0 && prev[prev.length - 1].elapsed_years === s.elapsed_years ? prev : [...prev, s],
      );
    } catch {
      // ignored -- see comment above
    }
  }, []);

  const handleGenerate = useCallback(async () => {
    setBusy(true);
    setError(null);
    try {
      const s = await generateWorld(
        seed, continentalPercent / 100, landPercent / 100, axialTiltDeg, detail, initialSoilMaturityPercent / 100,
        detail, fluidDensity, autoPlates ? null : numPlates,
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
      setSimulatePlateMovement(DEFAULT_SIMULATE_PLATE_MOVEMENT);
      setSimulateClimateBiomes(DEFAULT_SIMULATE_CLIMATE_BIOMES);
      setWindModel(DEFAULT_WIND_MODEL);
      await Promise.all([refresh(projection, mapView, rotation), refreshPlates(), refreshRivers(), refreshLakes(), recordStats()]);
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }, [
    seed, continentalPercent, landPercent, axialTiltDeg, detail, fluidDensity, initialSoilMaturityPercent, autoPlates, numPlates,
    projection, mapView, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, recordStats,
  ]);


  // Debounced so dragging a Controls slider doesn't fire a network request (and force a
  // climate recompute, see main.py's /world/controls) on every single pixel of movement --
  // only once movement has paused briefly. Local slider state (seaLevelM/solarMultiplier)
  // still updates immediately on every change, so the slider itself never feels laggy.
  const controlsDebounceRef = useRef<ReturnType<typeof setTimeout> | null>(null);
  const pushControls = useCallback((next: {
    seaLevelM?: number;
    solarMultiplier?: number;
    simulatePlateMovement?: boolean;
    simulateClimateBiomes?: boolean;
    windModel?: string;
  }) => {
    if (controlsDebounceRef.current) clearTimeout(controlsDebounceRef.current);
    controlsDebounceRef.current = setTimeout(async () => {
      try {
        await updateControls(next);
        await Promise.all([refresh(projection, mapView, rotation), recordStats()]);
      } catch (e) {
        setError(String(e));
      }
    }, 150);
  }, [projection, mapView, rotation, refresh, recordStats]);

  const handleSeaLevelChange = useCallback((v: number) => {
    setSeaLevelM(v);
    pushControls({ seaLevelM: v });
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

  const handleWindModelChange = useCallback((v: string) => {
    setWindModel(v);
    pushControls({ windModel: v });
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
        refresh(projection, mapViewRef.current, rotation), refreshPlates(), refreshRivers(), refreshLakes(), recordStats(),
      ]);
    } catch (e) {
      setError(String(e));
      setPlaying(false);
    } finally {
      setStepping(false);
    }
  }, [summary, stepYears, projection, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, recordStats]);

  // FileModal's "Load World" -- a loaded world fully replaces the current one, same as a
  // fresh Generate (see handleGenerate above), plus syncing every live Controls value
  // (seaLevelM/solarMultiplier/simulatePlateMovement/simulateClimateBiomes) to the *loaded*
  // world's own real values: calling updateControls with no fields set changes nothing but
  // still returns the current world's current values (see api.ts's updateControls/backend
  // app/main.py's /world/controls), which is simpler than adding a new endpoint just to read
  // them back.
  const handleWorldReplaced = useCallback(async (s: WorldSummary) => {
    setBusy(true);
    setError(null);
    try {
      setSummary(s);
      setSelectedPlateId(null);
      setSelectedRiverId(null);
      setSelectedBasin(null);
      setSelectedBasinKind(null);
      setStatsHistory([]);
      const controls = await updateControls({});
      setSeaLevelM(controls.sea_level_m);
      setSolarMultiplier(controls.solar_multiplier);
      setSimulatePlateMovement(controls.simulate_plate_movement);
      setSimulateClimateBiomes(controls.simulate_climate_biomes);
      setWindModel(controls.wind_model);
      await Promise.all([refresh(projection, mapView, rotation), refreshPlates(), refreshRivers(), refreshLakes(), recordStats()]);
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }, [projection, mapView, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, recordStats]);

  // FileModal's "Make Animation" -- it already advanced the world for real (see
  // api.ts's animateWorld), so this just runs the same post-step refresh handleStep does.
  const handleWorldAdvanced = useCallback(async (s: WorldSummary) => {
    setSummary(s);
    setSelectedRiverId(null);
    setSelectedBasin(null);
    setSelectedBasinKind(null);
    await Promise.all([refresh(projection, mapView, rotation), refreshPlates(), refreshRivers(), refreshLakes(), recordStats()]);
  }, [projection, mapView, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, recordStats]);

  // Resets the view orientation back to the default (lat=0/lon=0, see rotation.ts) -- just
  // updates state, same as a completed drag; the effect below does the actual re-fetch.
  const handleRecenter = useCallback(() => {
    setRotation(IDENTITY_ROTATION);
    setCenterLatLon({ lat: 0, lon: 0 });
  }, []);

  const selectedPlate = platesData.find((p) => p.plate_id === selectedPlateId) ?? null;
  const selectedRiver = riversData.find((r) => r.river_id === selectedRiverId) ?? null;

  // Re-render with the current world whenever the projection, map view, or view rotation
  // changes -- all three are baked server-side into the returned image (see api.ts's
  // renderWorld). A completed drag (see MapCanvas's onRotationCommitted below) just updates
  // `rotation` state; this effect is what actually re-fetches, the same pattern
  // projection/mapView already used before rotation existed.
  useEffect(() => {
    if (summary) {
      refresh(projection, mapView, rotation);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [projection, mapView, rotation]);

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

          <button onClick={() => setShowGenerateDialog(true)} disabled={busy} style={{ fontSize: 12 }}>
            Generate World
          </button>

          <button onClick={() => setShowStatsModal(true)} disabled={!summary} style={{ fontSize: 12 }}>
            📊 Stats
          </button>

          <button
            onClick={() => setShowControlsModal(true)}
            disabled={!summary}
            style={{ fontSize: 12 }}
          >
            🎛️ Controls
          </button>

          <button onClick={() => setShowFileModal(true)} style={{ fontSize: 12 }}>
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
            <div style={{ display: "flex", gap: 6 }}>
              <button onClick={handleStep} disabled={busy || stepping || !summary} style={{ flex: 1, fontSize: 12 }}>
                Step
              </button>
              <button onClick={() => setPlaying((p) => !p)} disabled={busy || !summary} style={{ flex: 1, fontSize: 12 }}>
                {playing ? "Pause" : "Play"}
              </button>
            </div>
          </fieldset>

          <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
            <legend style={{ fontSize: 11 }}>Map View</legend>
            <select
              value={mapView}
              onChange={(e) => setMapView(e.target.value as MapView)}
              style={{ width: "100%", marginBottom: 6, fontSize: 12 }}
            >
              <optgroup label="Tectonics & Climate">
                <option value="plates">Plates</option>
                <option value="platesDetail">Plates (details)</option>
                <option value="plateInspector">Plate Inspector</option>
                <option value="riverInspector">River Inspector</option>
                <option value="lakeInspector">Lake Inspector</option>
                <option value="elevation">Elevation</option>
                <option value="temperature">Temperature</option>
                <option value="wind">Wind</option>
                <option value="oceanCurrents">Ocean currents</option>
                <option value="humidity">Humidity</option>
                <option value="precipitation">Precipitation</option>
                <option value="biome">Biome</option>
                <option value="combined">Combined</option>
                <option value="resources">Resources</option>
                <option value="soilQuality">Soil Quality</option>
              </optgroup>
            </select>
            <select
              value={projection}
              onChange={(e) => setProjection(e.target.value as Projection)}
              style={{ width: "100%", fontSize: 12 }}
            >
              <option value="behrmann">Behrmann (cylindrical equal-area)</option>
              <option value="eckert4">Eckert IV (pseudocylindrical equal-area)</option>
            </select>
            <div style={{ marginTop: 6, opacity: 0.8 }}>
              Center: {formatLatLon(centerLatLon.lat, centerLatLon.lon)}
            </div>
            <button
              onClick={handleRecenter}
              disabled={busy || !summary || isIdentityRotation(rotation)}
              style={{ width: "100%", marginTop: 6, fontSize: 12 }}
            >
              Re-center
            </button>
          </fieldset>

          {mapView === "plateInspector" && (
            <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
              <legend style={{ fontSize: 11 }}>Selected plate</legend>
              {selectedPlate ? (
                <div style={{ opacity: 0.9 }}>
                  <div>id: {selectedPlate.plate_id}</div>
                  <div>crust: {selectedPlate.crust_type}</div>
                  {selectedPlate.num_rows != null && <div>rows: {selectedPlate.num_rows}</div>}
                  <div>points: {selectedPlate.num_points}</div>
                  {selectedPlate.bounding_ellipse && (
                    <>
                      <div>diameter A: {selectedPlate.bounding_ellipse.diameter_a_km.toFixed(0)} km</div>
                      <div>diameter B: {selectedPlate.bounding_ellipse.diameter_b_km.toFixed(0)} km</div>
                    </>
                  )}
                </div>
              ) : (
                <div style={{ opacity: 0.6 }}>Click a plate, or press Tab.</div>
              )}
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
              <div>elapsed: {(summary.elapsed_years / 1e6).toFixed(1)} Myr</div>
            </div>
          )}
          {error && <div style={{ color: "#ff8080", fontSize: 11 }}>{error}</div>}

          <EventConsole events={summary?.events ?? []} />
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
            />
          ) : (
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
              alphaEncodedIds={mapView === "combined"}
            />
          )}
          </div>
          <Legend
            mapView={mapView}
            highlightedBiome={highlightedBiome}
            onBiomeClick={(label) => setHighlightedBiome((cur) => (cur === label ? null : label))}
          />
          <p style={{ fontSize: 11, opacity: 0.6, marginTop: 4 }}>
            {mapView === "plateInspector"
              ? "Click a plate to select it. Tab / Shift+Tab cycles plates. Press and hold, then drag to rotate."
              : mapView === "riverInspector"
                ? "Click a river to select it. Tab / Shift+Tab cycles rivers. Press and hold, then drag to rotate."
                : mapView === "lakeInspector"
                  ? "Click a lake or any point on land to inspect its basin. Tab / Shift+Tab cycles lakes. Press and hold, then drag to rotate."
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
            <h2 style={{ fontSize: 16, marginTop: 0 }}>Generate World</h2>

            <label style={{ display: "block", marginBottom: 16 }}>
              Seed
              <div style={{ display: "flex", gap: 6 }}>
                <input
                  type="number"
                  value={seed}
                  onChange={(e) => setSeed(Number(e.target.value))}
                  style={{ flex: 1 }}
                />
                <button
                  type="button"
                  title="Randomize seed"
                  aria-label="Randomize seed"
                  onClick={() => setSeed(randomSeed())}
                >
                  🎲
                </button>
              </div>
            </label>

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
                Biome/Combined/Resources/Soil Quality maps and more elevation-line nodes -- but
                simulation steps and rendering both run slower. Lower runs faster but coarser.
              </div>
            </label>

            <div style={{ display: "flex", gap: 8, justifyContent: "space-between" }}>
              <button type="button" onClick={() => setShowAdvancedSettings(true)} disabled={busy}>
                Advanced settings
              </button>
              <div style={{ display: "flex", gap: 8 }}>
                <button onClick={() => setShowGenerateDialog(false)} disabled={busy}>
                  Cancel
                </button>
                <button onClick={handleGenerate} disabled={busy}>
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
          axialTiltDeg={axialTiltDeg}
          initialSoilMaturityPercent={initialSoilMaturityPercent}
          fluidDensity={fluidDensity}
          fluidDensityChoices={FLUID_DETAIL_CHOICES}
          onLandPercentChange={setLandPercent}
          onContinentalPercentChange={setContinentalPercent}
          onAutoPlatesChange={setAutoPlates}
          onNumPlatesChange={setNumPlates}
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
          simulatePlateMovement={simulatePlateMovement}
          simulateClimateBiomes={simulateClimateBiomes}
          windModel={windModel}
          onSeaLevelChange={handleSeaLevelChange}
          onSolarMultiplierChange={handleSolarMultiplierChange}
          onSimulatePlateMovementChange={handleSimulatePlateMovementChange}
          onSimulateClimateBiomesChange={handleSimulateClimateBiomesChange}
          onWindModelChange={handleWindModelChange}
          onClose={() => setShowControlsModal(false)}
        />
      )}

      {showFileModal && (
        <FileModal
          hasWorld={!!summary}
          seed={summary?.seed ?? null}
          elapsedYears={summary?.elapsed_years ?? null}
          projection={projection}
          mapView={mapView}
          rotation={rotation}
          renderWidth={RENDER_WIDTH}
          renderHeight={RENDER_HEIGHT}
          mapWrapperRef={mapWrapperRef}
          onClose={() => setShowFileModal(false)}
          onWorldReplaced={handleWorldReplaced}
          onWorldAdvanced={handleWorldAdvanced}
        />
      )}
    </div>
  );
}
