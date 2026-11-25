export const LAKE_STEP_YEARS_MIN = 50;
export const LAKE_STEP_YEARS_MAX = 1_000;
export const GLACIER_STEP_YEARS_MIN = 100;
export const GLACIER_STEP_YEARS_MAX = 1_000;
export const PLATE_STEP_YEARS_MIN = 1_000;
export const PLATE_STEP_YEARS_MAX = 100_000;

const PLATE_LOG_MIN = Math.log10(PLATE_STEP_YEARS_MIN);
const PLATE_LOG_MAX = Math.log10(PLATE_STEP_YEARS_MAX);
const PLATE_LOG_SAMPLES = 200;

function roundPlateYears(raw: number): number {
  const step = raw < 10_000 ? 100 : 1_000;
  const rounded = Math.round(raw / step) * step;
  return Math.min(PLATE_STEP_YEARS_MAX, Math.max(PLATE_STEP_YEARS_MIN, rounded));
}

// Strictly increasing, deduplicated log-spaced ticks from 1,000 to 100,000 yr so every
// keyboard arrow press on the range input advances to a distinct year value (PR #1000).
export const PLATE_STEP_YEARS_TICKS: readonly number[] = (() => {
  const seen = new Set<number>([PLATE_STEP_YEARS_MIN, 10_000, PLATE_STEP_YEARS_MAX]);
  for (let i = 0; i <= PLATE_LOG_SAMPLES; i++) {
    const logVal = PLATE_LOG_MIN + (i / PLATE_LOG_SAMPLES) * (PLATE_LOG_MAX - PLATE_LOG_MIN);
    seen.add(roundPlateYears(10 ** logVal));
  }
  return Array.from(seen).sort((a, b) => a - b);
})();

export function plateYearsToTickIndex(years: number): number {
  const clamped = Math.min(PLATE_STEP_YEARS_MAX, Math.max(PLATE_STEP_YEARS_MIN, years));
  const targetLog = Math.log10(clamped);
  let bestIndex = 0;
  let bestDist = Infinity;
  for (let i = 0; i < PLATE_STEP_YEARS_TICKS.length; i++) {
    const dist = Math.abs(Math.log10(PLATE_STEP_YEARS_TICKS[i]) - targetLog);
    if (dist < bestDist) {
      bestDist = dist;
      bestIndex = i;
    }
  }
  return bestIndex;
}

export function maxSubstepsBound(stepYears: number, substepYears: number): number {
  return Math.max(1, Math.ceil(stepYears / Math.max(substepYears, 1)));
}

export function clampMaxSubsteps(raw: number, maxBound: number): number {
  if (!Number.isFinite(raw)) return 1;
  return Math.min(maxBound, Math.max(1, Math.round(raw)));
}
