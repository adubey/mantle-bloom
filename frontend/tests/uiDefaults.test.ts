import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { TUNING_MULTIPLIER_KEYS } from "../src/api.js";
import {
  DEFAULT_GLACIER_STEP_YEARS,
  DEFAULT_GLACIER_STOP_MODE,
  DEFAULT_KEEP_LAND_OFF_POLES,
  DEFAULT_LAKE_STEP_YEARS,
  DEFAULT_LAKE_STOP_MODE,
  DEFAULT_SIMULATE_GLACIER_MARCH,
  DEFAULT_SPLIT_METHOD,
  DEFAULT_STEP_YEARS,
  DEFAULT_TUNING,
  freshWorldControls,
} from "../src/uiDefaults.js";
import { PLATE_STEP_YEARS_TICKS } from "../src/timeSettings.js";

// The standard UI defaults from #1007. Backend model defaults are unchanged;
// these pin only what a fresh UI session / fresh UI world starts from.
describe("uiDefaults (#1007)", () => {
  it("sets the standard defaults", () => {
    assert.equal(DEFAULT_KEEP_LAND_OFF_POLES, true);
    assert.equal(DEFAULT_STEP_YEARS, 10_000);
    assert.equal(DEFAULT_SIMULATE_GLACIER_MARCH, true);
    assert.equal(DEFAULT_GLACIER_STEP_YEARS, 250);
    assert.equal(DEFAULT_GLACIER_STOP_MODE, "resolve");
    assert.equal(DEFAULT_LAKE_STEP_YEARS, 50);
    assert.equal(DEFAULT_LAKE_STOP_MODE, "resolve");
    assert.equal(DEFAULT_SPLIT_METHOD, "rift_lines");
  });

  it("sets the #998 plate-motion setting and leaves every other knob at 1x", () => {
    assert.equal(DEFAULT_TUNING.contact_drag_multiplier, 0.5);
    assert.equal(DEFAULT_TUNING.craton_drag_multiplier, 1);
    assert.equal(DEFAULT_TUNING.mantle_flow_multiplier, 4);
    for (const key of TUNING_MULTIPLIER_KEYS) {
      if (key === "contact_drag_multiplier" || key === "mantle_flow_multiplier") continue;
      assert.equal(DEFAULT_TUNING[key], 1, key);
    }
  });

  it("pushes exactly the divergent defaults to a fresh world", () => {
    assert.deepEqual(freshWorldControls(), {
      lakeStepYears: 50,
      lakeStopMode: "resolve",
      glacierStepYears: 250,
      glacierStopMode: "resolve",
      splitMethod: "rift_lines",
      ridgePushForm: "potential_energy",
      slabDragForm: "trench_normal",
      ageForceColumn: true,
      tuning: DEFAULT_TUNING,
    });
  });

  it("keeps the plate step on the Time slider's tick grid", () => {
    assert.ok(PLATE_STEP_YEARS_TICKS.includes(DEFAULT_STEP_YEARS));
  });
});
