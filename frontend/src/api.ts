export const API_BASE = "http://localhost:8000";

export type Projection = "behrmann" | "eckert4";

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

// No total plate count here -- the world tiles itself into a plausible number from the seed
// alone (see backend app/plates.py's generate_plates). numContinents is the continents
// slider -- exactly that many plates are made continental.
export function generateWorld(seed: number, numContinents: number): Promise<WorldSummary> {
  return fetch(`${API_BASE}/world/generate`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ seed, num_continents: numContinents }),
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
