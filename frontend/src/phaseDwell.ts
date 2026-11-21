import { phaseLinks } from "./phaseInfo";

// Dwell and lag-bounded queue for Viewer's "Learn more" phase links.
// Kept pure (no DOM/React imports) so the state machine can be unit-tested directly.

export const LEARN_MORE_DWELL_MS = 1000;
export const LEARN_MORE_MAX_LAG_MS = 2000;

export interface PhaseDwellOptions {
  onPhaseChange: (phase: string | null) => void;
  hasLinks?: (phase: string) => boolean;
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

export function createPhaseDwellQueue({
  onPhaseChange,
  hasLinks = (phase) => phaseLinks(phase).length > 0,
  dwellMs = LEARN_MORE_DWELL_MS,
  maxLagMs = LEARN_MORE_MAX_LAG_MS,
  now = () => Date.now(),
  setTimer = setTimeout,
  clearTimer = clearTimeout,
}: PhaseDwellOptions) {
  let displayed: string | null = null;
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
    if (queue.length > 0) {
      const next = queue.shift()!;
      show(next.phase);
      return;
    }
    if (actualPhase === null) {
      show(null);
    } else if (actualPhase !== displayed) {
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

    if (timer === null) {
      if (nextPhase !== displayed) show(nextPhase);
      return;
    }

    const tailPhase = queue.length > 0 ? queue[queue.length - 1].phase : displayed;
    if (nextPhase === tailPhase) return;

    expungeStale(t);
    queue.push({ phase: nextPhase, enqueuedAt: t });
  };

  const reset = () => {
    cancelTimer();
    queue.length = 0;
    actualPhase = null;
    if (displayed !== null) {
      displayed = null;
      onPhaseChange(null);
    }
  };

  const dispose = () => {
    disposed = true;
    cancelTimer();
    queue.length = 0;
  };

  return { update, reset, dispose };
}
