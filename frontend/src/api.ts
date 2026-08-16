export const API_BASE = "http://localhost:8000";

export type Projection = "behrmann" | "eckert4";

export interface WorldSummary {
  seed: number;
  elapsed_years: number;
  num_plates: number;
}

export interface RenderLine {
  points: [number, number][];
  elevation: number[];
}

export interface VelocityArrow {
  start: [number, number];
  end: [number, number];
}

export interface RenderPlate {
  plate_id: number;
  crust_type: "continental" | "oceanic";
  lines: RenderLine[];
  pole: [number, number] | null;
  rotation_rate_deg_per_myr: number;
  velocity_arrow: VelocityArrow | null;
  boundary: [number, number][];
}

export interface RenderResponse {
  projection: Projection;
  elapsed_years: number;
  plates: RenderPlate[];
}

async function asJson<T>(resp: Response): Promise<T> {
  if (!resp.ok) {
    const detail = await resp.text();
    throw new Error(`${resp.status} ${resp.statusText}: ${detail}`);
  }
  return resp.json() as Promise<T>;
}

// No plate count here -- the world tiles itself into a plausible number of plates from the
// seed alone (see backend app/plates.py's generate_plates).
export function generateWorld(seed: number): Promise<WorldSummary> {
  return fetch(`${API_BASE}/world/generate`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ seed }),
  }).then(asJson<WorldSummary>);
}

export function stepWorld(years: number): Promise<WorldSummary> {
  return fetch(`${API_BASE}/world/step`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ years }),
  }).then(asJson<WorldSummary>);
}

export function renderWorld(projection: Projection): Promise<RenderResponse> {
  return fetch(`${API_BASE}/world/render?projection=${projection}`).then(asJson<RenderResponse>);
}
