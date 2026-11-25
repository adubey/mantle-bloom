// Persists whether the user collapsed the animation RunOverlay to the sidebar
// ProgressBar across page refreshes (#1033). Scoped to sessionStorage so a
// refresh in the same tab restores the user's chosen view, while a brand-new
// tab defaults to the expanded RunOverlay.

const KEY_PREFIX = "mantle-bloom:animation-overlay-collapsed:";

function defaultStorage(): Storage | null {
  try {
    return typeof window !== "undefined" ? window.sessionStorage : null;
  } catch {
    return null;
  }
}

export function readAnimationOverlayCollapsed(worldId: string, storage: Storage | null = defaultStorage()): boolean {
  if (!storage) return false;
  try {
    return storage.getItem(`${KEY_PREFIX}${worldId}`) === "true";
  } catch {
    return false;
  }
}

export function writeAnimationOverlayCollapsed(
  worldId: string,
  collapsed: boolean,
  storage: Storage | null = defaultStorage(),
): void {
  if (!storage) return;
  try {
    storage.setItem(`${KEY_PREFIX}${worldId}`, collapsed ? "true" : "false");
  } catch {
    // best-effort: storage may be unavailable in restricted contexts
  }
}

export function clearAnimationOverlayCollapsed(worldId: string, storage: Storage | null = defaultStorage()): void {
  if (!storage) return;
  try {
    storage.removeItem(`${KEY_PREFIX}${worldId}`);
  } catch {
    // best-effort
  }
}

// Format a rolling step duration for the animation progress UI: seconds/step,
// switching to minutes/step once a step takes more than 120 seconds (#1033).
export function formatStepDuration(secondsPerStep: number): string {
  if (secondsPerStep > 120) {
    return `${(secondsPerStep / 60).toFixed(1)} min/step`;
  }
  if (secondsPerStep >= 10) {
    return `${secondsPerStep.toFixed(1)} s/step`;
  }
  return `${secondsPerStep.toFixed(2)} s/step`;
}
