import { API_BASE } from "./api";
import { asJson } from "./apiResponse";

// The logged-out viewer's routes (backend app/main.py's /public/*, app/public_view.py, #373).
// The allowance rides on an HttpOnly cookie, which a same-origin fetch sends by default.

export type PublicView = "elevation" | "biome";

export interface PublicWorld {
  revision: string;
  seed: number;
  elapsed_years: number;
}

export interface PublicState {
  world: PublicWorld | null;
  // null: unlimited (the owner, viewing /view).
  steps_remaining: number | null;
  step_limit: number;
  step_years: number;
}

export interface PublicStepResult {
  world: PublicWorld;
  // False when the world had already moved on, or another viewer's step was running, and this
  // request was served that result without computing (or charging) anything.
  computed: boolean;
  steps_remaining: number | null;
}

// A step the server refused with a status the viewer explains: 429 (allowance spent) or 503
// (the owner is using the world).
export class PublicStepRefused extends Error {
  readonly status: number;
  constructor(status: number) {
    super(`step refused (${status})`);
    this.status = status;
  }
}

export function fetchPublicState(): Promise<PublicState> {
  return fetch(`${API_BASE}/public/world`).then(asJson<PublicState>);
}

export function renderPublicWorld(view: PublicView, width: number, height: number): Promise<PublicWorld & { image_base64: string }> {
  const params = new URLSearchParams({ view, width: String(width), height: String(height) });
  return fetch(`${API_BASE}/public/render?${params}`).then(asJson<PublicWorld & { image_base64: string }>);
}

export async function stepPublicWorld(revision: string): Promise<PublicStepResult> {
  const resp = await fetch(`${API_BASE}/public/step`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ revision }),
  });
  if (resp.status === 429 || resp.status === 503) throw new PublicStepRefused(resp.status);
  return asJson<PublicStepResult>(resp);
}
