import { phaseLinks } from "./phaseInfo";

// Dwell and lag-bounded queue for Viewer's "Learn more" phase links.
// Kept pure (no DOM/React imports) so the state machine can be unit-tested directly.

export const LEARN_MORE_DWELL_MS = 1000;
// Maximum age of a queued update before it is dropped from the front of the queue;
// displayed links can trail the live phase by up to maxLagMs + dwellMs.
export const LEARN_MORE_MAX_LAG_MS = 2000;

export interface PhaseDwellOptions {
  onPhaseChange: (phase: string | null) => void;
  hasLinks?: (phase: string) => boolean;
  phaseKey?: (phase: string) => string;
  dwellMs?: number;
  maxLagMs?: number;
  now?: () => number;
  setTimer?: (fn: () => void, ms: number) => ReturnType<typeof setTimeout>;
  clearTimer?: (id: ReturnType<typeof setTimeout>) => void;
}

interface QueuedPhase {
  phase: string;
  enqueuedAt: number;
}

function defaultPhaseKey(phase: string): string {
  const links = phaseLinks(phase);
  return links.length > 0 ? links.map((link) => link.url).join("\n") : phase;
}

export function createPhaseDwellQueue({
  onPhaseChange,
  hasLinks = (phase) => phaseLinks(phase).length > 0,
  phaseKey = defaultPhaseKey,
  dwellMs = LEARN_MORE_DWELL_MS,
  maxLagMs = LEARN_MORE_MAX_LAG_MS,
  now = () => Date.now(),
  setTimer = setTimeout,
  clearTimer = clearTimeout,
}: PhaseDwellOptions) {
  let displayed: string | null = null;
  let displayedKey: string | null = null;
  let actualPhase: string | null = null;
  const queue: QueuedPhase[] = [];
  let timer: ReturnType<typeof setTimeout> | null = null;
  let disposed = false;

  const cancelTimer = () => {
    if (timer !== null) {
      clearTimer(timer);
      timer = null;
    }
  };

  const expungeStale = (t: number) => {
    while (queue.length > 0 && t - queue[0].enqueuedAt > maxLagMs) {
      queue.shift();
    }
  };

  const show = (phase: string | null) => {
    cancelTimer();
    if (phase !== displayed) {
      displayed = phase;
      displayedKey = phase !== null ? phaseKey(phase) : null;
      onPhaseChange(phase);
    }
    if (phase !== null) {
      timer = setTimer(() => {
        timer = null;
        advance();
      }, dwellMs);
    }
  };

  const advance = () => {
    if (disposed) return;
    const t = now();
    expungeStale(t);
    while (queue.length > 0 && phaseKey(queue[0].phase) === displayedKey) {
      queue.shift();
    }
    if (queue.length > 0) {
      const next = queue.shift()!;
      show(next.phase);
      return;
    }
    if (actualPhase === null) {
      show(null);
    } else if (phaseKey(actualPhase) !== displayedKey) {
      show(actualPhase);
    }
  };

  const update = (rawPhase: string | null | undefined) => {
    if (disposed) return;
    const t = now();
    const nextPhase = rawPhase && hasLinks(rawPhase) ? rawPhase : null;
    actualPhase = nextPhase;

    if (nextPhase === null) {
      if (displayed === null || timer !== null) return;
      show(null);
      return;
    }

    if (displayed === null) {
      show(nextPhase);
      return;
    }

    const nextKey = phaseKey(nextPhase);
    if (timer === null) {
      if (nextKey !== displayedKey) show(nextPhase);
      return;
    }

    expungeStale(t);
    const tailKey = queue.length > 0 ? phaseKey(queue[queue.length - 1].phase) : displayedKey;
    if (nextKey === tailKey) return;

    queue.push({ phase: nextPhase, enqueuedAt: t });
  };

  const dispose = () => {
    disposed = true;
    cancelTimer();
    queue.length = 0;
  };

  return { update, dispose };
}
