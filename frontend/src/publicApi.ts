import { API_BASE, type WorldId } from "./api";
import { asJson } from "./apiResponse";
import { publicWorldBase } from "./worldRoutes";

// The logged-out viewer's routes (backend app/main.py's /public/worlds/*,
// app/public_view.py, #373, #563). Each addresses one world the owner marked
// public. The allowance rides on an HttpOnly cookie, which a same-origin fetch
// sends by default; it is per browser across all public worlds.

export type PublicView = "elevation" | "biome" | "combined";

export interface PublicWorld {
  revision: string;
  seed: number;
  elapsed_years: number;
}

export interface PublicWorldListing extends PublicWorld {
  world_id: WorldId;
  name: string;
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

export function fetchPublicWorlds(): Promise<{ worlds: PublicWorldListing[] }> {
  return fetch(`${API_BASE}/public/worlds`).then(asJson<{ worlds: PublicWorldListing[] }>);
}

// With `after` (the revision on screen), the server holds the request until the
// addressed world moves past it or ~25 s pass, so every viewer follows other
// viewers' steps (#374).
export function fetchPublicState(worldId: WorldId, after?: string, signal?: AbortSignal, guestDemo = false): Promise<PublicState> {
  const params = new URLSearchParams();
  if (after !== undefined) params.set("after", after);
  if (guestDemo) params.set("guest_demo", "true");
  const encoded = params.toString();
  const query = encoded ? `?${encoded}` : "";
  return fetch(`${publicWorldBase(API_BASE, worldId)}/world${query}`, { signal }).then(asJson<PublicState>);
}

export function renderPublicWorld(worldId: WorldId, view: PublicView, width: number, height: number): Promise<PublicWorld & { image_base64: string }> {
  const params = new URLSearchParams({ view, width: String(width), height: String(height) });
  return fetch(`${publicWorldBase(API_BASE, worldId)}/render?${params}`).then(asJson<PublicWorld & { image_base64: string }>);
}

export async function stepPublicWorld(worldId: WorldId, revision: string, guestDemo = false): Promise<PublicStepResult> {
  const query = guestDemo ? "?guest_demo=true" : "";
  const resp = await fetch(`${publicWorldBase(API_BASE, worldId)}/step${query}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ revision }),
  });
  if (resp.status === 429 || resp.status === 503) throw new PublicStepRefused(resp.status);
  return asJson<PublicStepResult>(resp);
}
