import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import "./index.css";
import {
  animateWorld, fetchEarthquakes, fetchFaults, fetchLakes, fetchPlates, fetchPointSample, fetchRivers, fetchStats, fetchVolcanoes, fetchWorldSummary, generateWorld, renderWorld, stepWorld, updateControls,
  TUNING_MULTIPLIER_KEYS,
} from "./api";
import type {
  AnimateResponse, EarthquakeSummary, FaultSummary, FaultSystemSummary, LakeAtResponse, LakeSummary, MapView, PlateSummary, PointSample, Projection, RenderResponse, RiverSummary, Segment, TuningKey, TuningMultipliers, VolcanoSummary, WorldStats, WorldSummary,
} from "./api";
import MapCanvas from "./MapCanvas";
import PlateInspector from "./PlateInspector";
import RiverInspector from "./RiverInspector";
import LakeInspector from "./LakeInspector";
import PlatesAndFaults from "./PlatesAndFaults";
import EventConsole from "./EventConsole";
import StatsModal from "./StatsModal";
import ControlsModal from "./ControlsModal";
import AdvancedSettingsModal from "./AdvancedSettingsModal";
import FileModal from "./FileModal";
import Legend from "./Legend";
import { faultKindForLegendLabel, highlightTargetFor } from "./legendData";
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
  "resources", "soilQuality", "geomorph", "elevReason", "overlapAge", "plateInspector", "riverInspector", "lakeInspector", "platesAndFaults",
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
  const [mapView, setMapView] = useState<MapView>(initialView?.mapView ?? "combined");
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
  // Geomorphic-budget tuning knobs (see DEFAULT_TUNING / backend World's *_multiplier
  // group) -- one object of dimensionless multipliers, live-adjustable via Controls, reset
  // to all-1.0 on a fresh Generate and synced from the loaded world on Load.
  const [tuning, setTuning] = useState<TuningMultipliers>(DEFAULT_TUNING);
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
  // Background animation (see FileModal's "Start Animation" and handleStartAnimation below).
  // `animation` is non-null only while a run is in flight -- it holds the live frame count for
  // the sidebar progress bar; `animationResult` holds the finished MP4 (+ world summary) until
  // the user saves or dismisses it. The run streams each frame's PNG straight onto the main
  // map, and holds the server's world lock throughout, so every world-mutating / rendering
  // control is disabled while `animation` is set (see `animating`). `animCancelRef` carries
  // the AbortController that the sidebar "Cancel" button trips.
  const [animation, setAnimation] = useState<{ frame: number; total: number } | null>(null);
  const [animationResult, setAnimationResult] = useState<AnimateResponse | null>(null);
  const animCancelRef = useRef<AbortController | null>(null);
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
      setIceAgePeriodYears(DEFAULT_ICE_AGE_PERIOD_YEARS);
      setSimulatePlateMovement(DEFAULT_SIMULATE_PLATE_MOVEMENT);
      setSimulateClimateBiomes(DEFAULT_SIMULATE_CLIMATE_BIOMES);
      setWindModel(DEFAULT_WIND_MODEL);
      setFaultDeformationMode(DEFAULT_FAULT_DEFORMATION_MODE);
      setTuning(DEFAULT_TUNING);
      await Promise.all([refresh(projection, mapView, rotation), refreshPlates(), refreshRivers(), refreshLakes(), refreshFaults(), recordStats()]);
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }, [
    seed, continentalPercent, landPercent, axialTiltDeg, detail, fluidDensity, initialSoilMaturityPercent, autoPlates, numPlates,
    projection, mapView, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, refreshFaults, recordStats,
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
        refresh(projection, mapViewRef.current, rotation), refreshPlates(), refreshRivers(), refreshLakes(), refreshFaults(), recordStats(),
      ]);
    } catch (e) {
      setError(String(e));
      setPlaying(false);
    } finally {
      setStepping(false);
    }
  }, [summary, stepYears, projection, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, refreshFaults, recordStats]);

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
      setIceAgePeriodYears(controls.ice_age_period_years);
      setSimulatePlateMovement(controls.simulate_plate_movement);
      setSimulateClimateBiomes(controls.simulate_climate_biomes);
      setWindModel(controls.wind_model);
      setFaultDeformationMode(controls.fault_deformation_mode);
      setTuning(Object.fromEntries(TUNING_MULTIPLIER_KEYS.map((k) => [k, controls[k]])) as TuningMultipliers);
      await Promise.all([refresh(projection, mapView, rotation), refreshPlates(), refreshRivers(), refreshLakes(), refreshFaults(), recordStats()]);
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }, [projection, mapView, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, refreshFaults, recordStats]);

  // Post-animation refresh -- an animation run already advanced the world for real (see
  // api.ts's animateWorld), so this just runs the same post-step refresh handleStep does.
  const handleWorldAdvanced = useCallback(async (s: WorldSummary) => {
    setSummary(s);
    setSelectedRiverId(null);
    setSelectedBasin(null);
    setSelectedBasinKind(null);
    await Promise.all([refresh(projection, mapView, rotation), refreshPlates(), refreshRivers(), refreshLakes(), refreshFaults(), recordStats()]);
  }, [projection, mapView, rotation, refresh, refreshPlates, refreshRivers, refreshLakes, refreshFaults, recordStats]);

  // FileModal's "Start Animation" -- close the File dialog and run the whole animation in the
  // background. Each streamed frame's PNG is painted straight onto the main map (the run holds
  // the server world lock, so a normal render would 503 -- see api.ts/main.py), so the main
  // display doubles as the animation preview; the sidebar shows progress + a Cancel button,
  // and every world-mutating control is disabled until it finishes (see `animating`).
  const handleStartAnimation = useCallback(
    async ({ numFrames, yearsPerFrame }: { numFrames: number; yearsPerFrame: number }) => {
      setShowFileModal(false);
      setError(null);
      setAnimationResult(null);
      setAnimation({ frame: 0, total: numFrames });
      const ctrl = new AbortController();
      animCancelRef.current = ctrl;
      const view = mapViewRef.current;
      try {
        const result = await animateWorld(
          projection, view, RENDER_WIDTH, RENDER_HEIGHT, rotation, yearsPerFrame, numFrames,
          (p) => {
            setAnimation({ frame: p.frame, total: p.total });
            if (p.imageBase64) {
              // Paint the frame onto the main map. renderRequestIdRef is bumped so any
              // pre-animation refresh() still in flight can't clobber it afterward.
              renderRequestIdRef.current++;
              setRenderData({ projection, elapsed_years: 0, image_base64: p.imageBase64 });
            }
          },
          ctrl.signal,
        );
        setAnimationResult(result);
        await handleWorldAdvanced(result);
      } catch (e) {
        if (ctrl.signal.aborted) {
          // User pressed Cancel -- the world is left wherever the last completed frame put it.
          await handleWorldAdvanced(await fetchWorldSummary().catch(() => null) ?? summary!);
        } else {
          setError(String(e));
        }
      } finally {
        setAnimation(null);
        animCancelRef.current = null;
      }
    },
    [projection, rotation, handleWorldAdvanced, summary],
  );

  const handleCancelAnimation = useCallback(() => {
    animCancelRef.current?.abort(new DOMException("animation cancelled", "AbortError"));
  }, []);

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
  }, [animationResult]);

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

          <button onClick={() => setShowGenerateDialog(true)} disabled={busy || animating} style={{ fontSize: 12 }}>
            Generate World
          </button>

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
            <div style={{ display: "flex", gap: 6 }}>
              <button onClick={handleStep} disabled={busy || stepping || !summary || animating} style={{ flex: 1, fontSize: 12 }}>
                Step
              </button>
              <button onClick={() => setPlaying((p) => !p)} disabled={busy || !summary || animating} style={{ flex: 1, fontSize: 12 }}>
                {playing ? "Pause" : "Play"}
              </button>
            </div>
          </fieldset>

          <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
            <legend style={{ fontSize: 11 }}>Map View</legend>
            <select
              value={mapView}
              onChange={(e) => setMapView(e.target.value as MapView)}
              disabled={animating}
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
              <optgroup label="Debug &gt;">
                <option value="platesDetail">Points</option>
                <option value="speckle">Coastal dither (speckle)</option>
                <option value="platesAndFaults">Plates &amp; Faults</option>
                <option value="geomorph">Erosion &amp; Deposition</option>
                <option value="elevReason">Last elevation change</option>
                <option value="overlapAge">Plate overlap age</option>
                <option value="plateInspector">Plate Inspector</option>
                <option value="riverInspector">Rivers</option>
                <option value="lakeInspector">Lakes</option>
              </optgroup>
            </select>
            <select
              value={projection}
              onChange={(e) => setProjection(e.target.value as Projection)}
              disabled={animating}
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
              disabled={busy || !summary || animating || isIdentityRotation(rotation)}
              style={{ width: "100%", marginTop: 6, fontSize: 12 }}
            >
              Re-center
            </button>
          </fieldset>

          {(animation || animationResult) && (
            <fieldset style={{ border: "1px solid #3a4d8f", borderRadius: 6, padding: 8, fontSize: 12 }}>
              <legend style={{ fontSize: 11 }}>Animation</legend>
              {animation ? (
                <>
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
                  <div style={{ fontSize: 11, opacity: 0.7, margin: "4px 0 6px" }}>
                    Rendering frame {animation.frame} of {animation.total}… the map is stepping in the background.
                  </div>
                  <button onClick={handleCancelAnimation} style={{ width: "100%", fontSize: 12 }}>
                    Cancel
                  </button>
                </>
              ) : (
                <>
                  <div style={{ fontSize: 11, opacity: 0.7, marginBottom: 6 }}>Animation ready.</div>
                  <div style={{ display: "flex", gap: 6 }}>
                    <button onClick={handleSaveAnimation} style={{ flex: 1, fontSize: 12 }}>
                      Save MP4
                    </button>
                    <button onClick={() => setAnimationResult(null)} style={{ flex: 1, fontSize: 12 }}>
                      Dismiss
                    </button>
                  </div>
                </>
              )}
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
              onProbe={mapView === "combined" || mapView === "elevation" || mapView === "biome" ? handleProbe : undefined}
              alphaEncodedIds={mapView === "combined"}
              interactionDisabled={animating}
            />
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
                  : mapView === "platesAndFaults"
                  ? "Click a plate to select it (its fault strands emphasise, and its Euler pole + a speed-scaled motion arc appear). Tab / Shift+Tab cycles plates. Click a fault type in the legend to isolate that regime. Toggle the earthquake & volcano overlay in the sidebar. Press and hold, then drag to rotate."
                  : mapView === "combined" || mapView === "elevation" || mapView === "biome"
                    ? "Click any point for its elevation, biome, precipitation, temperature, and plate. Press and hold, then drag to rotate."
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
                Biome/Elevation &amp; Biome/Resources/Soil Quality maps and more elevation-line nodes -- but
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
          iceAgePeriodYears={iceAgePeriodYears}
          simulatePlateMovement={simulatePlateMovement}
          simulateClimateBiomes={simulateClimateBiomes}
          windModel={windModel}
          faultDeformationMode={faultDeformationMode}
          tuning={tuning}
          onSeaLevelChange={handleSeaLevelChange}
          onSolarMultiplierChange={handleSolarMultiplierChange}
          onIceAgePeriodChange={handleIceAgePeriodChange}
          onSimulatePlateMovementChange={handleSimulatePlateMovementChange}
          onSimulateClimateBiomesChange={handleSimulateClimateBiomesChange}
          onWindModelChange={handleWindModelChange}
          onFaultDeformationModeChange={handleFaultDeformationModeChange}
          onTuningChange={handleTuningChange}
          onTuningReset={handleTuningReset}
          onClose={() => setShowControlsModal(false)}
        />
      )}

      {showFileModal && (
        <FileModal
          hasWorld={!!summary}
          seed={summary?.seed ?? null}
          elapsedYears={summary?.elapsed_years ?? null}
          stepYears={stepYears}
          mapView={mapView}
          mapWrapperRef={mapWrapperRef}
          onClose={() => setShowFileModal(false)}
          onWorldReplaced={handleWorldReplaced}
          onStartAnimation={handleStartAnimation}
        />
      )}
    </div>
  );
}
