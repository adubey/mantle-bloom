export const LAKE_STEP_YEARS_MIN = 50;
export const LAKE_STEP_YEARS_MAX = 1_000;
export const GLACIER_STEP_YEARS_MIN = 100;
export const GLACIER_STEP_YEARS_MAX = 1_000;
export const PLATE_STEP_YEARS_MIN = 1_000;
export const PLATE_STEP_YEARS_MAX = 100_000;

export const PLATE_LOG_MIN = Math.log10(PLATE_STEP_YEARS_MIN);
export const PLATE_LOG_MAX = Math.log10(PLATE_STEP_YEARS_MAX);

export function plateYearsToLog(years: number): number {
  const clamped = Math.min(PLATE_STEP_YEARS_MAX, Math.max(PLATE_STEP_YEARS_MIN, years));
  return Math.log10(clamped);
}

export function logToPlateYears(logValue: number): number {
  const raw = 10 ** Math.min(PLATE_LOG_MAX, Math.max(PLATE_LOG_MIN, logValue));
  const step = raw < 10_000 ? 100 : 1_000;
  const rounded = Math.round(raw / step) * step;
  return Math.min(PLATE_STEP_YEARS_MAX, Math.max(PLATE_STEP_YEARS_MIN, rounded));
}
