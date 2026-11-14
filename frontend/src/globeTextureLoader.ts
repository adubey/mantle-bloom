// Scheduling and retry controller for GlobeSelector's elevation texture (#501).
// Kept pure (no DOM/React imports) so the state machine can be unit-tested directly.

export interface GlobeTextureTarget {
  worldId: string;
  worldRevision: number;
  visible: boolean;
}

export interface GlobeTextureLoaderOptions<T> {
  loadPixels: (worldId: string, worldRevision: number) => Promise<T>;
  onLoaded: (pixels: T) => void;
  onCleared?: () => void;
  debounceMs?: number;
  retryBaseMs?: number;
  retryMaxMs?: number;
  setTimer?: (fn: () => void, ms: number) => ReturnType<typeof setTimeout>;
  clearTimer?: (id: ReturnType<typeof setTimeout>) => void;
}

const DEFAULT_DEBOUNCE_MS = 250;
const DEFAULT_RETRY_BASE_MS = 1000;
const DEFAULT_RETRY_MAX_MS = 15000;

export function createGlobeTextureLoader<T>({
  loadPixels,
  onLoaded,
  onCleared,
  debounceMs = DEFAULT_DEBOUNCE_MS,
  retryBaseMs = DEFAULT_RETRY_BASE_MS,
  retryMaxMs = DEFAULT_RETRY_MAX_MS,
  setTimer = setTimeout,
  clearTimer = clearTimeout,
}: GlobeTextureLoaderOptions<T>) {
  let current: GlobeTextureTarget = { worldId: "", worldRevision: 0, visible: false };
  let loaded: { worldId: string; worldRevision: number } | null = null;
  let inFlight: { worldId: string; worldRevision: number } | null = null;
  let timer: ReturnType<typeof setTimeout> | null = null;
  let retryAttempt = 0;
  let generation = 0;
  let disposed = false;

  const cancelTimer = () => {
    if (timer !== null) {
      clearTimer(timer);
      timer = null;
    }
  };

  const isTargetValid = (target: GlobeTextureTarget) =>
    Boolean(target.worldId) && target.worldRevision > 0 && target.visible;

  const isUpToDate = () =>
    loaded !== null &&
    loaded.worldId === current.worldId &&
    loaded.worldRevision === current.worldRevision;

  const scheduleFetch = (delayMs: number) => {
    cancelTimer();
    if (disposed || !isTargetValid(current) || isUpToDate()) return;
    if (delayMs <= 0) {
      startFetch();
      return;
    }
    timer = setTimer(() => {
      timer = null;
      startFetch();
    }, delayMs);
  };

  const startFetch = () => {
    cancelTimer();
    if (disposed || !isTargetValid(current) || isUpToDate()) return;
    if (
      inFlight !== null &&
      inFlight.worldId === current.worldId &&
      inFlight.worldRevision === current.worldRevision
    ) {
      return;
    }

    const { worldId, worldRevision } = current;
    const reqGen = ++generation;
    inFlight = { worldId, worldRevision };

    loadPixels(worldId, worldRevision).then(
      (pixels) => {
        if (disposed || reqGen !== generation) return;
        inFlight = null;
        if (current.worldId !== worldId) return;
        loaded = { worldId, worldRevision };
        retryAttempt = 0;
        onLoaded(pixels);
        if (current.worldRevision !== worldRevision && isTargetValid(current)) {
          scheduleFetch(debounceMs);
        }
      },
      () => {
        if (disposed || reqGen !== generation) return;
        inFlight = null;
        if (!isTargetValid(current)) return;
        if (current.worldId !== worldId || current.worldRevision !== worldRevision) {
          scheduleFetch(loaded === null ? 0 : debounceMs);
          return;
        }
        const delay = Math.min(retryMaxMs, retryBaseMs * 2 ** retryAttempt);
        retryAttempt += 1;
        timer = setTimer(() => {
          timer = null;
          startFetch();
        }, delay);
      },
    );
  };

  return {
    update(next: GlobeTextureTarget) {
      if (disposed) return;
      const worldChanged = next.worldId !== current.worldId;
      const revisionChanged = next.worldRevision !== current.worldRevision;

      if (worldChanged) {
        loaded = null;
        inFlight = null;
        retryAttempt = 0;
        generation += 1;
        cancelTimer();
        onCleared?.();
      } else if (revisionChanged) {
        retryAttempt = 0;
        cancelTimer();
      }

      if (!next.visible) {
        retryAttempt = 0;
        cancelTimer();
      }

      current = next;
      if (!isTargetValid(current) || isUpToDate()) {
        cancelTimer();
        return;
      }

      if (
        inFlight !== null &&
        inFlight.worldId === current.worldId &&
        inFlight.worldRevision === current.worldRevision
      ) {
        return;
      }

      const immediate = loaded === null && inFlight === null;
      scheduleFetch(immediate ? 0 : debounceMs);
    },

    dispose() {
      disposed = true;
      cancelTimer();
      generation += 1;
      inFlight = null;
    },
  };
}
