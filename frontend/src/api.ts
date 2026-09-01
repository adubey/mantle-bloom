// Overridable at build/dev time via VITE_API_BASE (see bin/restart.sh's --backend-port), so a
// non-default backend port stays wired up correctly instead of silently pointing at :8000.
const API_ROOT = import.meta.env.VITE_API_BASE ?? "http://localhost:8000";

export const API_BASE = API_ROOT;

export type Projection = "behrmann" | "eckert4";

// The "kind of map" requested from /world/render -- the server now bakes this choice
// directly into the returned image (fill color rule, whether boundaries/poles/velocity
// arrows or raw per-plate node dots are drawn), rather than the frontend deciding how to
// draw raw coordinate data it fetched.
export type MapView =
  | "elevation"
  | "plates"
  | "platesDetail"
  | "speckle"
  | "temperature"
  | "wind"
  | "oceanCurrents"
  | "humidity"
  | "precipitation"
  | "biome"
  | "combined"
  | "resources"
  | "soilQuality"
  | "geomorph"
  | "plateInspector"
  | "riverInspector"
  | "lakeInspector";

export interface WorldEvent {
  elapsed_years: number;
  message: string;
}

export interface WorldSummary {
  seed: number;
  elapsed_years: number;
  num_plates: number;
  events: WorldEvent[];
}

export interface RenderResponse {
  projection: Projection;
  elapsed_years: number;
  // Raw PNG bytes, base64-encoded -- decode via `data:image/png;base64,${image_base64}`.
  // See backend app/render_image.py for how it's drawn.
  image_base64: string;
}

// Plate Inspector data -- unlike RenderResponse, this is raw JSON the client renders itself
// (see PlateInspector.tsx), not a baked PNG. Un-rotated/true-frame throughout: the client
// applies its current view rotation only at draw time (see rotation.ts).
export interface BoundingEllipse {
  center_xyz: [number, number, number];
  diameter_a_km: number; // major
  diameter_b_km: number; // minor
  outline: [number, number, number][];
}

export interface PlateOverlap {
  plate_id: number;
  fraction: number; // share of THIS plate's nodes sitting on top of plate_id
}

export interface PlateSummary {
  plate_id: number;
  crust_type: "continental" | "oceanic";
  // null for a Plate representation with no row concept (see backend app/main.py's
  // _plate_summary).
  num_rows: number | null;
  num_points: number;
  // Plate motion + shape health (see _plate_summary). `at_max_rate` true for most plates at
  // once is a long-run pathology; a continental plate with a high `submerged_fraction` is
  // over-stretched. `overlaps`/`collisions` expose stalled territory conflicts.
  speed_cm_per_yr: number;
  at_max_rate: boolean;
  euler_pole: { lat_deg: number; lon_deg: number } | null;
  age_steps: number;
  median_elevation_m: number | null;
  submerged_fraction: number;
  overlaps: PlateOverlap[];
  collisions: { plate_id: number; years: number }[];
  outline: [number, number, number][];
  // Every node's own position (not just the outline loop) -- see PlateInspector.tsx, which
  // plots these individually, bright for the selected plate and dim for every other one.
  points: [number, number, number][];
  bounding_ellipse: BoundingEllipse | null;
}

export interface PlatesResponse {
  elapsed_years: number;
  plates: PlateSummary[];
}

// River Inspector data -- same "raw JSON, client renders itself" philosophy as
// PlatesResponse (see RiverInspector.tsx). `river_id` is only meaningful against the most
// recent /world/rivers response -- rivers are regrouped fresh on every call, not given a
// persistent identity across steps (see backend app/hydrology.py's group_rivers).
export interface RiverSummary {
  river_id: number;
  num_nodes: number;
  // Each entry is one flow edge (a node and its downstream flow_target) as a pair of world-
  // space points -- a flat edge list, not an ordered polyline, since a river network can
  // branch (see PlateInspector.tsx's segment-based river rendering for the same shape used
  // server-side in render_image.py's _draw_rivers).
  segments: Segment[];
  mouth_xyz: [number, number, number];
  mouth_type: "ocean" | "lake" | "other";
  flow_rate: number;
  speed: number;
  num_tributaries: number;
}

// A world-space segment pair -- shared shape for both RiverSummary.segments and
// coastline_segments below.
export type Segment = [[number, number, number], [number, number, number]];

export interface RiversResponse {
  elapsed_years: number;
  rivers: RiverSummary[];
  // The land/lake-vs-ocean boundary (see backend app/coastline.py) -- included here because
  // the River Inspector draws no filled backdrop at all, so without this there's no
  // land/ocean cue in this view whatsoever.
  coastline_segments: Segment[];
}

// A river feeding into a lake/basin (see LakeSummary.inflow_rivers below) -- a small inline
// summary, not a full RiverSummary, since the Lake Inspector only needs to say "a river ends
// here, this big" rather than draw that river's own full path (that's what the River
// Inspector's own map mode is for).
export interface InflowRiver {
  mouth_xyz: [number, number, number];
  flow_rate: number;
  num_nodes: number;
}

// Lake Inspector data -- same "raw JSON, client renders itself" philosophy as PlatesResponse/
// RiversResponse (see LakeInspector.tsx). `lake_id` is only meaningful against the most recent
// /world/lakes response -- lakes are regrouped fresh every call, same non-persistent-identity
// contract as `river_id` (see backend app/hydrology.py's group_rivers/lake_components).
// This same shape is also embedded (with `lake_id: null`) in a `LakeAtResponse` for a click
// that lands on a dry, never-flooded basin -- see backend app/main.py's `_lake_basin_summary`.
export interface LakeSummary {
  lake_id: number | null;
  is_lake: boolean;
  member_count: number;
  // Every member node's own position -- what the Lake Inspector plots per node, mirroring
  // PlateSummary.points (bright for the selected basin, dim for every other one).
  member_xyz: [number, number, number][];
  floor_xyz: [number, number, number];
  floor_elevation_m: number;
  // "The lowest point of the edge of the basin" -- the saddle a river out of it would source
  // from. `null` for an unresolved closed/endorheic basin with no known spill (a legitimate
  // state, not missing data -- see backend app/lakes.py's own docstring).
  outlet_xyz: [number, number, number] | null;
  outlet_elevation_m: number | null;
  // `null` unless `is_lake` (a dry basin has no water surface at all).
  water_elevation_m: number | null;
  is_spilling: boolean;
  inflow_rivers: InflowRiver[];
}

export interface LakesResponse {
  elapsed_years: number;
  lakes: LakeSummary[];
  coastline_segments: Segment[];
}

// The Lake Inspector's click hit-test result -- see backend app/main.py's `lake_at` for what
// each `kind` means. `basin` is `null` for "ocean"/"no_basin", set otherwise.
export interface LakeAtResponse {
  kind: "lake" | "basin" | "no_basin" | "ocean";
  basin: LakeSummary | null;
}

// Stats panel data (see backend app/stats.py) -- a stateless snapshot of the *current* world;
// the client is what accumulates a history over time (see App.tsx's recordStats), same
// division of responsibility /world/render already has with renderData. Land/ocean
// temperature/elevation/depth fields are `null` when their domain has no grid cells at all
// (e.g. an all-ocean world has no land cells, so elevation_*/land_*/air_* are all null; a
// world with no ocean at all would similarly null out ocean_depth_*/ocean_temperature_*) --
// StatsModal must handle that, not assume every field is always present.
export interface WorldStats {
  elapsed_years: number;
  // Single running totals, not a spatial snapshot distribution like every field below --
  // the Simulation tab is what turns a run of these into a min/max/mean/std-dev over time.
  plate_count: number;
  elevation_point_count: number;
  land_fraction: number;
  ocean_fraction: number;
  // Land only (height above the current sea level) -- see stats.py's module docstring.
  elevation_min_m: number | null;
  elevation_max_m: number | null;
  elevation_mean_m: number | null;
  elevation_std_m: number | null;
  // Ocean only, positive meters below the current sea level.
  ocean_depth_min_m: number | null;
  ocean_depth_max_m: number | null;
  ocean_depth_mean_m: number | null;
  ocean_depth_std_m: number | null;
  land_temperature_min_c: number | null;
  land_temperature_max_c: number | null;
  land_temperature_mean_c: number | null;
  land_temperature_std_c: number | null;
  air_temperature_min_c: number | null;
  air_temperature_max_c: number | null;
  air_temperature_mean_c: number | null;
  air_temperature_std_c: number | null;
  ocean_temperature_min_c: number | null;
  ocean_temperature_max_c: number | null;
  ocean_temperature_mean_c: number | null;
  ocean_temperature_std_c: number | null;
  precipitation_min_mm: number;
  precipitation_max_mm: number;
  precipitation_mean_mm: number;
  precipitation_std_mm: number;
  // Percent of *land* cells (0 to 1 each) in each named biome (see backend app/biomes.py) --
  // "Ocean" is never a key here (it's always 0% of land, by construction) and the dict is
  // `{}` entirely when there are no land cells at all.
  biome_land_fraction: Record<string, number>;
}

async function asJson<T>(resp: Response): Promise<T> {
  if (!resp.ok) {
    const detail = await resp.text();
    throw new Error(`${resp.status} ${resp.statusText}: ${detail}`);
  }
  return resp.json() as Promise<T>;
}

async function asBlob(resp: Response): Promise<Blob> {
  if (!resp.ok) {
    const detail = await resp.text();
    throw new Error(`${resp.status} ${resp.statusText}: ${detail}`);
  }
  return resp.blob();
}

// `numPlates`: the "Number of plates" slider (see backend app/plates.py's MIN_AUTO_PLATES/
// MAX_AUTO_PLATES) -- `null`/omitted means "Auto," the world's original behavior of tiling
// itself into a plausible plate count drawn from the seed alone (see plates.generate_plates).
// continentalFraction/landFraction are the dialog's two sliders (0 to 1): the fraction of
// plates made continental, and the fraction of the whole sphere that starts above sea level
// (independent of the first -- see plates.py's _land_noise_threshold for how the two
// combine). axialTiltDeg is the dialog's fourth slider (degrees) -- doesn't affect plate
// generation, only climate.py's insolation at render time (see world.py's
// World.axial_tilt_deg). nodeDensity is the dialog's "point density" choice (0.5, 1, 2, or 4,
// see plates.NODE_DENSITY_CHOICES) -- how many elevation-line nodes each plate starts with, and
// stays scaled to for the rest of that world's life (see world.py's World.node_density).
// initialSoilMaturity is the dialog's fifth slider (0 to 1) -- how much soil the world starts
// with (0 = fully barren, matching every other persistent field's own zero-start default; see
// backend app/geology.py's seed_initial_soil), a one-time generation-time seed, not stored on
// World the way nodeDensity is (see world.generate_world's own docstring for why).
// climateDensity is the dialog's "climate & biome resolution" choice (0.5, 1, 2, or 4, see
// backend app/climate.py's CLIMATE_DENSITY_CHOICES) -- how finely climate.py's own grid (and the
// Biome/Combined/Resources/Soil-Quality views' own render grid, scaled the same way) resolves
// temperature/wind/humidity/precipitation; stored on World for the rest of that world's life,
// same reasoning nodeDensity's own storage gives (see world.py's World.climate_density).
// fluidDensity is the Advanced-settings dialog's "Fluid dynamics resolution" choice (same
// 0.5/1/2/4 set) -- how finely the continuous Ocean/Atmospheric CFD's own grid resolves
// currents/wind, independent of climateDensity (see world.py's World.fluid_density).
export function generateWorld(
  seed: number,
  continentalFraction: number,
  landFraction: number,
  axialTiltDeg: number,
  nodeDensity: number,
  initialSoilMaturity: number,
  climateDensity: number,
  fluidDensity: number,
  numPlates: number | null,
): Promise<WorldSummary> {
  return fetch(`${API_BASE}/world/generate`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      seed,
      num_plates: numPlates,
      continental_fraction: continentalFraction,
      land_fraction: landFraction,
      axial_tilt_deg: axialTiltDeg,
      node_density: nodeDensity,
      initial_soil_maturity: initialSoilMaturity,
      climate_density: climateDensity,
      fluid_density: fluidDensity,
    }),
  }).then(asJson<WorldSummary>);
}

// Live-adjusts sea level, solar heat, and/or the plate-movement/climate-biomes step toggles
// on the current world (the "Controls" window -- see backend app/world.py's
// World.sea_level_m/World.solar_multiplier/World.simulate_plate_movement/
// World.simulate_climate_biomes and main.py's /world/controls) -- any argument may be
// omitted to leave that value untouched. Forces an immediate server-side climate recompute,
// so the caller should re-fetch /world/render and /world/stats right after this resolves to
// see the effect.
export interface ControlsState {
  sea_level_m: number;
  solar_multiplier: number;
  simulate_plate_movement: boolean;
  simulate_climate_biomes: boolean;
  wind_model: string;
}

export function updateControls(controls: {
  seaLevelM?: number;
  solarMultiplier?: number;
  simulatePlateMovement?: boolean;
  simulateClimateBiomes?: boolean;
  windModel?: string;
}): Promise<ControlsState> {
  return fetch(`${API_BASE}/world/controls`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      sea_level_m: controls.seaLevelM,
      solar_multiplier: controls.solarMultiplier,
      simulate_plate_movement: controls.simulatePlateMovement,
      simulate_climate_biomes: controls.simulateClimateBiomes,
      wind_model: controls.windModel,
    }),
  }).then(asJson<ControlsState>);
}

// The backend rejects an overlapping /world/step with 503 (see backend app/main.py's
// _step_lock) rather than queueing it, since a step already in flight can take seconds on a
// large world. Retry on a short delay until the in-flight step finishes and the lock frees up
// -- any other error (network failure, 404 for no world, etc.) still propagates immediately.
const STEP_RETRY_DELAY_MS = 300;

export async function stepWorld(years: number): Promise<WorldSummary> {
  for (;;) {
    const resp = await fetch(`${API_BASE}/world/step`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ years }),
    });
    if (resp.status === 503) {
      await new Promise((resolve) => setTimeout(resolve, STEP_RETRY_DELAY_MS));
      continue;
    }
    return asJson<WorldSummary>(resp);
  }
}

// width/height are the actual pixel dimensions of the returned image -- pass a higher
// resolution than the canvas's displayed CSS size for a sharper (retina-style) render; the
// server scales line widths/dot/pole sizes to match (see render_image.py's pixel_scale).
// `rotation` is the map's current view orientation (see rotation.ts and
// docs/simulation-model.md#rotating-the-view) -- a flattened row-major 3x3 matrix, default
// identity (center at lat=0/lon=0) when omitted, matching backend main.py's own default.
export function renderWorld(
  projection: Projection,
  view: MapView,
  width: number,
  height: number,
  rotation?: number[],
): Promise<RenderResponse> {
  const params = new URLSearchParams({ projection, view, width: String(width), height: String(height) });
  if (rotation) params.set("rotation", rotation.join(","));
  return fetch(`${API_BASE}/world/render?${params}`).then(asJson<RenderResponse>);
}

export function fetchPlates(): Promise<PlatesResponse> {
  return fetch(`${API_BASE}/world/plates`).then(asJson<PlatesResponse>);
}

// The Plate Inspector's click hit-test -- (latDeg, lonDeg) must already be in the *true*
// (un-rotated) frame; the caller unprojects the click through whatever view rotation is
// currently active first (see rotation.ts's unproject/matTranspose).
export function fetchPlateAt(latDeg: number, lonDeg: number): Promise<{ plate_id: number | null }> {
  const params = new URLSearchParams({ lat_deg: String(latDeg), lon_deg: String(lonDeg) });
  return fetch(`${API_BASE}/world/plate_at?${params}`).then(asJson<{ plate_id: number | null }>);
}

export interface PointSample {
  lat_deg: number;
  lon_deg: number;
  elevation_m: number;
  is_ocean: boolean;
  biome_id: number;
  biome: string;
  temperature_c: number;
  precipitation_mm: number;
  plate_id: number | null;
}

// The Elevation & Biome / Elevation / Biome views' click-to-inspect popup (see App.tsx) --
// like fetchPlateAt, (latDeg, lonDeg) must already be in the *true* (un-rotated) frame; the
// caller unprojects the click through the active view rotation first.
export function fetchPointSample(latDeg: number, lonDeg: number): Promise<PointSample> {
  const params = new URLSearchParams({ lat_deg: String(latDeg), lon_deg: String(lonDeg) });
  return fetch(`${API_BASE}/world/sample_at?${params}`).then(asJson<PointSample>);
}

export function fetchStats(): Promise<WorldStats> {
  return fetch(`${API_BASE}/world/stats`).then(asJson<WorldStats>);
}

// Re-syncs with whatever world is already sitting in server memory -- used on page load to
// restore the map after a browser refresh (see App.tsx) without generating a new world.
// Rejects (404) if no world has been generated yet in this server session.
export function fetchWorldSummary(): Promise<WorldSummary> {
  return fetch(`${API_BASE}/world/summary`).then(asJson<WorldSummary>);
}

export function fetchRivers(): Promise<RiversResponse> {
  return fetch(`${API_BASE}/world/rivers`).then(asJson<RiversResponse>);
}

// The River Inspector's click hit-test -- same true-frame contract as fetchPlateAt.
export function fetchRiverAt(latDeg: number, lonDeg: number): Promise<{ river_id: number | null }> {
  const params = new URLSearchParams({ lat_deg: String(latDeg), lon_deg: String(lonDeg) });
  return fetch(`${API_BASE}/world/river_at?${params}`).then(asJson<{ river_id: number | null }>);
}

export function fetchLakes(): Promise<LakesResponse> {
  return fetch(`${API_BASE}/world/lakes`).then(asJson<LakesResponse>);
}

// "File > Save World" -- the entire current world state as an opaque binary blob (see
// backend app/persistence.py -- deliberately pickle, no cross-version compatibility
// promise). The caller downloads this via the usual SPA blob-download pattern
// (URL.createObjectURL + a synthetic <a download>), not a plain navigation, so the browser
// doesn't need a same-origin/CORS-friendly direct link.
export function saveWorld(): Promise<Blob> {
  return fetch(`${API_BASE}/world/save`).then(asBlob);
}

// "File > Load World" -- `file` is the raw bytes of a file saveWorld() previously produced
// (an <input type="file">'s File object, or any Blob). Replaces whatever world previously
// existed, same as generateWorld -- the caller should run the same full refresh sequence it
// runs after generateWorld.
export function loadWorld(file: Blob): Promise<WorldSummary> {
  return fetch(`${API_BASE}/world/load`, {
    method: "POST",
    headers: { "Content-Type": "application/octet-stream" },
    body: file,
  }).then(asJson<WorldSummary>);
}

export interface AnimateResponse extends WorldSummary {
  // The rendered animation, base64-encoded -- decode via `data:${mime};base64,${videoBase64}`.
  videoBase64: string;
  // The video's MIME type (currently always "video/mp4" -- H.264, see backend
  // app/render_image.py's stream_animation_mp4).
  mime: string;
}

// Progress during a Make Animation run -- `frame` of `total` frames rendered so far.
export interface AnimateProgress {
  frame: number;
  total: number;
}

// Each animation frame is a full step_world + render (see backend app/main.py's
// /world/animate), and up to MAX_ANIMATION_FRAMES=240 of those run back-to-back server-side
// before the video is complete -- by far the slowest request the app makes, and one that's
// only gotten slower as the simulation itself has picked up more per-step work (erosion,
// sediment redistribution, coastline stabilization, ...). A big/slow world's 240-frame worst
// case genuinely can run past any fixed whole-request deadline (the old 15-min total cap was
// aborting healthy runs mid-render), so instead of bounding the whole run we bound only the
// gap *between* stream lines: the endpoint sends a `{type: "progress"}` line after every
// single frame, and no one frame's step_world + render comes anywhere near this long, so a
// silence past it means the server or connection actually died rather than the render merely
// being slow.
const ANIMATE_IDLE_TIMEOUT_MS = 3 * 60 * 1000;

// "File > Make Animation" -- renders `numFrames` frames of `view`/`projection`'s progress,
// starting from the world's current state (frame 0) and stepping it forward by
// `yearsPerFrame` real years between each subsequent frame, encoded as an H.264/MP4 video.
// **This permanently advances the world** by `(numFrames - 1) * yearsPerFrame` years, the
// same as calling stepWorld that many times -- not a side-effect-free preview (see backend
// app/main.py's /world/animate). The caller should run the same post-step refresh sequence
// it runs after stepWorld.
//
// The endpoint streams newline-delimited JSON: one `{type: "progress", frame, total}` per
// frame (surfaced via `onProgress`, for a progress bar), then a final `{type: "done", ...}`
// carrying the video. A mid-stream `{type: "error"}` line (rendering blew up after the
// response already 200'd) is re-thrown here like any other failure. The run is aborted only
// if the stream goes silent for ANIMATE_IDLE_TIMEOUT_MS -- an actively-progressing render,
// however long in total, keeps resetting that window.
export async function animateWorld(
  projection: Projection,
  view: MapView,
  width: number,
  height: number,
  rotation: number[] | undefined,
  yearsPerFrame: number,
  numFrames: number,
  onProgress?: (progress: AnimateProgress) => void,
): Promise<AnimateResponse> {
  // Idle watchdog: abort the fetch if no bytes arrive for ANIMATE_IDLE_TIMEOUT_MS. Rearmed
  // on every chunk (see the read loop below), so an actively-streaming run never trips it.
  const idleAbort = new AbortController();
  let idleTimer: ReturnType<typeof setTimeout> | undefined;
  const rearmIdleTimer = () => {
    clearTimeout(idleTimer);
    idleTimer = setTimeout(
      () => idleAbort.abort(new DOMException("animation stream stalled -- no data from server", "TimeoutError")),
      ANIMATE_IDLE_TIMEOUT_MS,
    );
  };
  rearmIdleTimer();

  try {
    const resp = await fetch(`${API_BASE}/world/animate`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        projection,
        view,
        width,
        height,
        rotation: rotation ? rotation.join(",") : null,
        years_per_frame: yearsPerFrame,
        num_frames: numFrames,
      }),
      signal: idleAbort.signal,
    });
    if (!resp.ok || !resp.body) {
      const detail = await resp.text();
      throw new Error(`${resp.status} ${resp.statusText}: ${detail}`);
    }

    const reader = resp.body.getReader();
    const decoder = new TextDecoder();
    let buffered = "";
    let done: AnimateResponse | null = null;

    const handleLine = (line: string) => {
      if (!line.trim()) return;
      const msg = JSON.parse(line) as Record<string, unknown>;
      if (msg.type === "progress") {
        onProgress?.({ frame: msg.frame as number, total: msg.total as number });
      } else if (msg.type === "error") {
        throw new Error(String(msg.detail));
      } else if (msg.type === "done") {
        done = {
          seed: msg.seed as number,
          elapsed_years: msg.elapsed_years as number,
          num_plates: msg.num_plates as number,
          events: msg.events as WorldEvent[],
          videoBase64: msg.video_base64 as string,
          mime: msg.mime as string,
        };
      }
    };

    for (;;) {
      const { value, done: streamDone } = await reader.read();
      if (streamDone) break;
      rearmIdleTimer();
      buffered += decoder.decode(value, { stream: true });
      const lines = buffered.split("\n");
      buffered = lines.pop() ?? "";
      for (const line of lines) handleLine(line);
    }
    if (buffered) handleLine(buffered);

    if (!done) throw new Error("animation stream ended without a result");
    return done;
  } finally {
    clearTimeout(idleTimer);
  }
}

// "File > Export Hex Grid" data -- a geodesic-icosahedron hex/pentagon tiling of the sphere
// (see backend app/geodesic.py), independent of the plate simulation's own node cloud, with
// the current world's elevation/biome sampled onto each tile. `vertices` is a shared,
// globally-indexed corner-point list -- `tiles[i].corner_vertex_ids`/`neighbor_ids` are
// index-aligned (neighbor_ids[k] is the tile sharing the edge between corner_vertex_ids[k]
// and corner_vertex_ids[(k+1) % length]) -- see docs/hex-export-format.md for the full
// shape and neighbor-finding pseudocode for a client application.
export interface HexTile {
  id: number;
  is_pentagon: boolean;
  center_lat_deg: number;
  center_lon_deg: number;
  center_xyz: [number, number, number];
  corner_vertex_ids: number[];
  elevation_m: number;
  is_ocean: boolean;
  biome: string;
  neighbor_ids: number[];
}

export interface HexGridExport {
  planet_radius_km: number;
  frequency: number;
  num_tiles: number;
  vertices: [number, number, number][];
  tiles: HexTile[];
}

// `frequency` is one of a small preset set (see backend app/geodesic.py's
// FREQUENCY_CHOICES) -- the UI offers it as a size dropdown, not a free-form input, same
// "discrete choices" pattern as node_density/climate_density in generateWorld.
export function exportHexGrid(frequency: number): Promise<HexGridExport> {
  return fetch(`${API_BASE}/world/export_hexgrid`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ frequency }),
  }).then(asJson<HexGridExport>);
}

// The Lake Inspector's click hit-test -- same true-frame contract as fetchPlateAt/fetchRiverAt,
// but (unlike those) always resolves to *something* informative for any land click, not just a
// hit on the thing being inspected -- see LakeAtResponse's own comment.
export function fetchLakeAt(latDeg: number, lonDeg: number): Promise<LakeAtResponse> {
  const params = new URLSearchParams({ lat_deg: String(latDeg), lon_deg: String(lonDeg) });
  return fetch(`${API_BASE}/world/lake_at?${params}`).then(asJson<LakeAtResponse>);
}
