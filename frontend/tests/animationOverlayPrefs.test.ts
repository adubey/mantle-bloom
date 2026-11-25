import assert from "node:assert/strict";
import { describe, it } from "node:test";
import {
  clearAnimationOverlayCollapsed,
  formatStepDuration,
  readAnimationOverlayCollapsed,
  writeAnimationOverlayCollapsed,
} from "../src/animationOverlayPrefs.ts";

function makeMemoryStorage(): Storage {
  const map = new Map<string, string>();
  return {
    get length() {
      return map.size;
    },
    clear() {
      map.clear();
    },
    getItem(key: string) {
      return map.has(key) ? map.get(key)! : null;
    },
    key(index: number) {
      return Array.from(map.keys())[index] ?? null;
    },
    removeItem(key: string) {
      map.delete(key);
    },
    setItem(key: string, value: string) {
      map.set(key, String(value));
    },
  };
}

describe("animationOverlayPrefs (#1033)", () => {
  it("defaults to expanded RunOverlay (false) when no prior preference exists", () => {
    const storage = makeMemoryStorage();
    assert.equal(readAnimationOverlayCollapsed("world-1", storage), false);
    assert.equal(readAnimationOverlayCollapsed("world-1", null), false);
  });

  it("persists collapsed ProgressBar and expanded RunOverlay per world across refreshes", () => {
    const storage = makeMemoryStorage();
    writeAnimationOverlayCollapsed("world-1", true, storage);
    writeAnimationOverlayCollapsed("world-2", false, storage);
    assert.equal(readAnimationOverlayCollapsed("world-1", storage), true);
    assert.equal(readAnimationOverlayCollapsed("world-2", storage), false);

    writeAnimationOverlayCollapsed("world-1", false, storage);
    assert.equal(readAnimationOverlayCollapsed("world-1", storage), false);
  });

  it("clears the stored preference when the animation ends", () => {
    const storage = makeMemoryStorage();
    writeAnimationOverlayCollapsed("world-1", true, storage);
    clearAnimationOverlayCollapsed("world-1", storage);
    assert.equal(readAnimationOverlayCollapsed("world-1", storage), false);
  });

  it("formats step durations in seconds/step up to 120s and minutes/step above 120s", () => {
    assert.equal(formatStepDuration(0.426), "0.43 s/step");
    assert.equal(formatStepDuration(12.34), "12.3 s/step");
    assert.equal(formatStepDuration(120), "120.0 s/step");
    assert.equal(formatStepDuration(126), "2.1 min/step");
    assert.equal(formatStepDuration(210), "3.5 min/step");
  });
});
