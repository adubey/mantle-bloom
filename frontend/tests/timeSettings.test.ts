import assert from "node:assert/strict";
import { describe, it } from "node:test";
import {
  clampMaxSubsteps,
  maxSubstepsBound,
  PLATE_STEP_YEARS_MAX,
  PLATE_STEP_YEARS_MIN,
  PLATE_STEP_YEARS_TICKS,
  plateYearsToTickIndex,
} from "../src/timeSettings.js";

describe("timeSettings log slider", () => {
  it("spans 1,000 to 100,000 yr with strictly increasing distinct ticks", () => {
    assert.equal(PLATE_STEP_YEARS_TICKS[0], PLATE_STEP_YEARS_MIN);
    assert.equal(PLATE_STEP_YEARS_TICKS.at(-1), PLATE_STEP_YEARS_MAX);
    assert.ok(PLATE_STEP_YEARS_TICKS.includes(10_000));
    for (let i = 1; i < PLATE_STEP_YEARS_TICKS.length; i++) {
      assert.ok(PLATE_STEP_YEARS_TICKS[i] > PLATE_STEP_YEARS_TICKS[i - 1]);
    }
  });

  it("round-trips every tick index so keyboard arrows never snap back", () => {
    for (let i = 0; i < PLATE_STEP_YEARS_TICKS.length; i++) {
      assert.equal(plateYearsToTickIndex(PLATE_STEP_YEARS_TICKS[i]), i);
    }
    assert.equal(plateYearsToTickIndex(100), 0);
    assert.equal(plateYearsToTickIndex(1_000_000), PLATE_STEP_YEARS_TICKS.length - 1);
  });

  it("computes and clamps max substeps bounds from step and substep years", () => {
    assert.equal(maxSubstepsBound(10_000, 250), 40);
    assert.equal(maxSubstepsBound(10_000, 300), 34);
    assert.equal(maxSubstepsBound(1_000, 2_000), 1);
    assert.equal(clampMaxSubsteps(0, 40), 1);
    assert.equal(clampMaxSubsteps(100, 20), 20);
    assert.equal(clampMaxSubsteps(Number.NaN, 20), 1);
  });
});
