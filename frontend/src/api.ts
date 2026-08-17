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

// A uniform lat/lon grid covering the whole sphere, each cell assigned its nearest
// elevation node -- gives the map full, gap-free coverage regardless of how sparse the
// underlying per-plate line data looks once projected. This is what MapCanvas actually
// paints as the filled map; RenderPlate.lines is the raw physical data (kept for
// reference/debugging), not what's drawn as the territory fill.
export interface RenderGrid {
  points: [number, number][];
  elevation: number[];
  plate_id: number[];
  crust_type: ("continental" | "oceanic")[];
  // Each cell's half-width/half-height in the *same projected units* as points -- sized
  // from the projection's actual local behavior at that latitude (see backend
  // main.py's _row_cell_half_extent), not a fixed screen size. Projected spacing isn't
  // uniform (e.g. Behrmann stretches longitude spacing near the poles by ~50x relative to
  // the equator), so drawing every cell the same pixel size leaves gaps somewhere no
  // matter what size is picked.
  cell_half_width: number[];
  cell_half_height: number[];
}

export interface RenderResponse {
  projection: Projection;
  elapsed_years: number;
  plates: RenderPlate[];
  grid: RenderGrid;
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

// includeLines defaults to (and should stay) false for every view except "Plates (details)"
// -- it's raw per-plate node data (position + elevation for every simulation node), ~2MB of
// JSON that the other views never read (see backend main.py's render() docstring).
export function renderWorld(projection: Projection, includeLines: boolean = false): Promise<RenderResponse> {
  return fetch(`${API_BASE}/world/render?projection=${projection}&include_lines=${includeLines}`).then(
    asJson<RenderResponse>,
  );
}
