// World-ID route helpers (#563): every world operation addresses its world
// explicitly. Pure functions so the frontend test suite can cover them.

export type WorldId = string;

/** The selected world, creating one on demand when there is none yet (#592).
 * Returns the current selection unchanged, so opening the Generate World
 * dialog stays click-to-open; only confirming Generate mints a world, never a
 * page load. */
export async function ensureWorldId(
  current: WorldId | null,
  create: () => Promise<WorldId>,
): Promise<WorldId> {
  if (current) return current;
  return create();
}

/** The most recently updated world's id, or null when the list is empty.
 * Centralizes the "adopt the latest world" rule (#592) shared by the initial
 * selection effect and the on-demand ensure path, so both pick the same world. */
export function latestWorldId(worlds: { world_id: WorldId; updated_at: number }[]): WorldId | null {
  if (worlds.length === 0) return null;
  return [...worlds].sort((a, b) => b.updated_at - a.updated_at)[0].world_id;
}

/** Base path of one world's admin routes, e.g. `/worlds/abc123`. */
export function worldBase(apiBase: string, worldId: WorldId): string {
  return `${apiBase}/worlds/${encodeURIComponent(worldId)}`;
}

/** Base path of one public world's routes, e.g. `/public/worlds/abc123`. */
export function publicWorldBase(apiBase: string, worldId: WorldId): string {
  return `${apiBase}/public/worlds/${encodeURIComponent(worldId)}`;
}

/** The `?world=` selection from a page URL, or null when absent/empty. */
export function worldIdFromSearch(search: string): WorldId | null {
  const id = new URLSearchParams(search).get("world");
  return id ? id : null;
}

/** `url` with its `world` parameter set, preserving every other parameter
 * (notably `deb`) and the hash. */
export function withWorldId(url: string, worldId: WorldId): string {
  const parsed = new URL(url, "http://placeholder.invalid");
  parsed.searchParams.set("world", worldId);
  return `${parsed.pathname}${parsed.search}${parsed.hash}`;
}
