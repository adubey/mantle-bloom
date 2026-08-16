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

export interface RenderPlate {
  plate_id: number;
  crust_type: "continental" | "oceanic";
  lines: RenderLine[];
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

export function generateWorld(seed: number, numPlates: number): Promise<WorldSummary> {
  return fetch(`${API_BASE}/world/generate`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ seed, num_plates: numPlates }),
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
