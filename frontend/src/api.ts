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
  | "precipitation";

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
export function renderWorld(
  projection: Projection,
  view: MapView,
  width: number,
  height: number,
): Promise<RenderResponse> {
  const params = new URLSearchParams({ projection, view, width: String(width), height: String(height) });
  return fetch(`${API_BASE}/world/render?${params}`).then(asJson<RenderResponse>);
}
