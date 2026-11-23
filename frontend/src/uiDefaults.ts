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
// (7) The #998 plate-motion setting: collision braking at 0.5x (craton component left at 1x) and the static mantle
// flow at 4x, every other geomorphic-budget knob at the model's untuned 1x. With the force forms below it is the one
// setting that reaches Earth's plate speeds and neighbour coherence together (analysis/issue998/README.md, sections
// 13 and 17: plate speed 3-5 cm/yr, contact median near 3 cm/yr, opposition 0.51, 10 of 20 development seeds accepted
// at 50 Myr). The 4x flow is a diagnostic amplitude, not a physical claim, and the rift-line parameters those runs
// used were RFC 0018's calibrated set, which the UI does not push.
export const DEFAULT_TUNING: TuningMultipliers = {
  ...(Object.fromEntries(TUNING_MULTIPLIER_KEYS.map((k) => [k, 1])) as TuningMultipliers),
  contact_drag_multiplier: 0.5,
  mantle_flow_multiplier: 4,
};
// (8) RFC 0021's force forms the #998 setting was run with (B.1 ridge push, B.2 sunk-slab drag, B.3 age force column).
export const DEFAULT_RIDGE_PUSH_FORM = "potential_energy";
export const DEFAULT_SLAB_DRAG_FORM = "trench_normal";
export const DEFAULT_AGE_FORCE_COLUMN = true;

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
  ridgePushForm: string;
  slabDragForm: string;
  ageForceColumn: boolean;
  tuning: TuningMultipliers;
} {
  return {
    lakeStepYears: DEFAULT_LAKE_STEP_YEARS,
    lakeStopMode: DEFAULT_LAKE_STOP_MODE,
    glacierStepYears: DEFAULT_GLACIER_STEP_YEARS,
    glacierStopMode: DEFAULT_GLACIER_STOP_MODE,
    splitMethod: DEFAULT_SPLIT_METHOD,
    ridgePushForm: DEFAULT_RIDGE_PUSH_FORM,
    slabDragForm: DEFAULT_SLAB_DRAG_FORM,
    ageForceColumn: DEFAULT_AGE_FORCE_COLUMN,
    tuning: DEFAULT_TUNING,
  };
}
