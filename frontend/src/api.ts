export const API_BASE = "http://localhost:8000";

export type Projection = "behrmann" | "eckert4";

// The "kind of map" requested from /world/render -- the server now bakes this choice
// directly into the returned image (fill color rule, whether boundaries/poles/velocity
// arrows or raw per-plate node dots are drawn), rather than the frontend deciding how to
// draw raw coordinate data it fetched.
export type MapView =
  | "elevation"
  | "plates"
  | "platesDetail"
  | "temperature"
  | "wind"
  | "oceanCurrents"
  | "humidity"
  | "precipitation"
  | "plateInspector";

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

export interface PlateSummary {
  plate_id: number;
  crust_type: "continental" | "oceanic";
  num_rows: number;
  num_points: number;
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

// Stats panel data (see backend app/stats.py) -- a stateless snapshot of the *current* world;
// the client is what accumulates a history over time (see App.tsx's recordStats), same
// division of responsibility /world/render already has with renderData. Land/air/ocean
// temperature and precipitation fields are `null` when their domain has no grid cells at all
// (e.g. an all-ocean world has no land cells, so every land_*/air_* field is null) --
// StatsModal must handle that, not assume every field is always present.
export interface WorldStats {
  elapsed_years: number;
  land_fraction: number;
  ocean_fraction: number;
  elevation_min_m: number;
  elevation_max_m: number;
  elevation_mean_m: number;
  land_temperature_min_c: number | null;
  land_temperature_max_c: number | null;
  land_temperature_mean_c: number | null;
  air_temperature_min_c: number | null;
  air_temperature_max_c: number | null;
  air_temperature_mean_c: number | null;
  ocean_temperature_min_c: number | null;
  ocean_temperature_max_c: number | null;
  ocean_temperature_mean_c: number | null;
  precipitation_min_mm: number;
  precipitation_max_mm: number;
  precipitation_mean_mm: number;
}

async function asJson<T>(resp: Response): Promise<T> {
  if (!resp.ok) {
    const detail = await resp.text();
    throw new Error(`${resp.status} ${resp.statusText}: ${detail}`);
  }
  return resp.json() as Promise<T>;
}

// No total plate count here -- the world tiles itself into a plausible number from the seed
// alone (see backend app/plates.py's generate_plates). continentalFraction/landFraction are
// the dialog's two sliders (0 to 1): the fraction of plates made continental, and the
// fraction of the whole sphere that starts above sea level (independent of the first --
// see plates.py's _land_noise_threshold for how the two combine). axialTiltDeg is the
// dialog's third slider (degrees) -- doesn't affect plate generation, only climate.py's
// insolation at render time (see world.py's World.axial_tilt_deg).
export function generateWorld(
  seed: number,
  continentalFraction: number,
  landFraction: number,
  axialTiltDeg: number,
): Promise<WorldSummary> {
  return fetch(`${API_BASE}/world/generate`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({
      seed,
      continental_fraction: continentalFraction,
      land_fraction: landFraction,
      axial_tilt_deg: axialTiltDeg,
    }),
  }).then(asJson<WorldSummary>);
}

export function stepWorld(years: number): Promise<WorldSummary> {
  return fetch(`${API_BASE}/world/step`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ years }),
  }).then(asJson<WorldSummary>);
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

export function fetchStats(): Promise<WorldStats> {
  return fetch(`${API_BASE}/world/stats`).then(asJson<WorldStats>);
}
