import { TUNING_MULTIPLIER_KEYS } from "./api";
import type { SplitMethod, StopMode, TuningMultipliers } from "./api";

// Standard UI defaults for a fresh session and every freshly generated UI
// world (#1007). These are UI choices only: the backend model defaults in
// backend/app/world.py are unchanged, so API/replay worlds keep the model
// values. Where a default below diverges from the backend world default,
// App's fresh-Generate flow pushes it to the new world through the existing
// controls endpoint (see freshWorldControls below); loading an existing
// world still adopts that world's own stored values instead.

// (1) Generate -> Advanced: keep initial land outside the polar caps.
export const DEFAULT_KEEP_LAND_OFF_POLES = true;
// (2) Tectonic plate step: 10,000 yr/step (Time dialog's Plates slider).
export const DEFAULT_STEP_YEARS = 10_000;
// (5) Lake substeps: 50 yr, stopping once lakes resolve.
export const DEFAULT_LAKE_STEP_YEARS = 50;
export const DEFAULT_LAKE_STOP_MODE: StopMode = "resolve";
// (3) Glacier march on, and (4) 250 yr substeps stopping once ice resolves.
// With the 10,000 yr plate step above a march is 40 substeps, not the
// minutes-long Myr-step march the old off default guarded against.
export const DEFAULT_SIMULATE_GLACIER_MARCH = true;
export const DEFAULT_GLACIER_STEP_YEARS = 250;
export const DEFAULT_GLACIER_STOP_MODE: StopMode = "resolve";
// (6) RFC 0018 rift-line plate splitting.
export const DEFAULT_SPLIT_METHOD: SplitMethod = "rift_lines";
// (7) Collision braking off: both brake multipliers at 0x, every other
// geomorphic-budget knob at the model's untuned 1x.
export const DEFAULT_TUNING: TuningMultipliers = {
  ...(Object.fromEntries(TUNING_MULTIPLIER_KEYS.map((k) => [k, 1])) as TuningMultipliers),
  contact_drag_multiplier: 0,
  craton_drag_multiplier: 0,
};

// The controls payload a freshly generated world needs so its stored values
// match the UI defaults above. Plate step years travel on step/Play, and
// keep-land-off-poles / glacier march travel on the generate request itself,
// so none of those belong here.
export function freshWorldControls(): {
  lakeStepYears: number;
  lakeStopMode: StopMode;
  glacierStepYears: number;
  glacierStopMode: StopMode;
  splitMethod: SplitMethod;
  tuning: TuningMultipliers;
} {
  return {
    lakeStepYears: DEFAULT_LAKE_STEP_YEARS,
    lakeStopMode: DEFAULT_LAKE_STOP_MODE,
    glacierStepYears: DEFAULT_GLACIER_STEP_YEARS,
    glacierStopMode: DEFAULT_GLACIER_STOP_MODE,
    splitMethod: DEFAULT_SPLIT_METHOD,
    tuning: DEFAULT_TUNING,
  };
}
