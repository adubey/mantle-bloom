import assert from "node:assert/strict";
import { describe, it } from "node:test";
import {
  logToPlateYears,
  PLATE_LOG_MAX,
  PLATE_LOG_MIN,
  PLATE_STEP_YEARS_MAX,
  PLATE_STEP_YEARS_MIN,
  plateYearsToLog,
} from "../src/timeSettings.js";

describe("timeSettings log slider", () => {
  it("maps endpoints and decade midpoints cleanly", () => {
    assert.equal(logToPlateYears(PLATE_LOG_MIN), PLATE_STEP_YEARS_MIN);
    assert.equal(logToPlateYears(4), 10_000);
    assert.equal(logToPlateYears(PLATE_LOG_MAX), PLATE_STEP_YEARS_MAX);
  });

  it("clamps out-of-range inputs and round-trips decade ticks", () => {
    assert.equal(plateYearsToLog(100), PLATE_LOG_MIN);
    assert.equal(plateYearsToLog(1_000_000), PLATE_LOG_MAX);
    assert.equal(logToPlateYears(plateYearsToLog(1_000)), 1_000);
    assert.equal(logToPlateYears(plateYearsToLog(10_000)), 10_000);
    assert.equal(logToPlateYears(plateYearsToLog(100_000)), 100_000);
  });
});
