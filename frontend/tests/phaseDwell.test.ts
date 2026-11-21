import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { createPhaseDwellQueue } from "../src/phaseDwell";

function createTestHarness() {
  let currentMs = 0;
  let nextId = 1;
  const timers = new Map<number, { fn: () => void; at: number }>();
  const history: Array<{ phase: string | null; at: number }> = [];

  const controller = createPhaseDwellQueue({
    onPhaseChange: (phase) => history.push({ phase, at: currentMs }),
    hasLinks: (phase) => phase !== "Unknown phase",
    dwellMs: 1000,
    maxLagMs: 2000,
    now: () => currentMs,
    setTimer: (fn, ms) => {
      const id = nextId++;
      timers.set(id, { fn, at: currentMs + ms });
      return id as unknown as ReturnType<typeof setTimeout>;
    },
    clearTimer: (id) => {
      timers.delete(id as unknown as number);
    },
  });

  const advanceTo = (targetMs: number) => {
    while (true) {
      let earliestId: number | null = null;
      let earliestAt = Infinity;
      for (const [id, entry] of timers) {
        if (entry.at <= targetMs && entry.at < earliestAt) {
          earliestAt = entry.at;
          earliestId = id;
        }
      }
      if (earliestId === null) break;
      const entry = timers.get(earliestId)!;
      timers.delete(earliestId);
      currentMs = entry.at;
      entry.fn();
    }
    currentMs = targetMs;
  };

  return { controller, history, advanceTo, setTime: (ms: number) => { currentMs = ms; } };
}

describe("createPhaseDwellQueue", () => {
  it("shows the first phase immediately and keeps each queued phase visible for ~1s", () => {
    const { controller, history, advanceTo } = createTestHarness();

    controller.update("Moving plates");
    assert.deepEqual(history, [{ phase: "Moving plates", at: 0 }]);

    advanceTo(200);
    controller.update("Reshaping plates");
    advanceTo(400);
    controller.update("Reshaping plates"); // duplicate fraction update
    advanceTo(500);
    controller.update(null); // step finished early

    assert.deepEqual(history, [{ phase: "Moving plates", at: 0 }]);

    advanceTo(1000);
    assert.deepEqual(history, [
      { phase: "Moving plates", at: 0 },
      { phase: "Reshaping plates", at: 1000 },
    ]);

    advanceTo(2000);
    assert.deepEqual(history, [
      { phase: "Moving plates", at: 0 },
      { phase: "Reshaping plates", at: 1000 },
      { phase: null, at: 2000 },
    ]);
  });

  it("clears immediately on null once the 1s dwell has already elapsed", () => {
    const { controller, history, advanceTo } = createTestHarness();

    controller.update("Moving plates");
    advanceTo(1500);
    assert.deepEqual(history, [{ phase: "Moving plates", at: 0 }]);

    controller.update(null);
    assert.deepEqual(history, [
      { phase: "Moving plates", at: 0 },
      { phase: null, at: 1500 },
    ]);
  });

  it("expunges stale items from the front of the queue when > 2s behind, then catches up to actual state", () => {
    const { controller, history, advanceTo } = createTestHarness();

    controller.update("P0");
    advanceTo(100);
    controller.update("P1");
    advanceTo(200);
    controller.update("P2");
    advanceTo(300);
    controller.update("P3");
    advanceTo(400);
    controller.update("P4");

    advanceTo(1000);
    assert.equal(history.at(-1)?.phase, "P1");

    advanceTo(2000);
    assert.equal(history.at(-1)?.phase, "P2");

    // At 3000ms, P3 (enqueued at 300ms) and P4 (400ms) are >2s old and expunged from the front,
    // but actualPhase is still P4 so it catches up to P4 directly.
    advanceTo(3000);
    assert.deepEqual(history, [
      { phase: "P0", at: 0 },
      { phase: "P1", at: 1000 },
      { phase: "P2", at: 2000 },
      { phase: "P4", at: 3000 },
    ]);
  });

  it("ignores phases without links and resets cleanly on world switch", () => {
    const { controller, history, advanceTo } = createTestHarness();

    controller.update("Unknown phase");
    assert.deepEqual(history, []);

    controller.update("Moving plates");
    advanceTo(200);
    controller.update("Unknown phase");
    controller.reset();

    assert.deepEqual(history, [
      { phase: "Moving plates", at: 0 },
      { phase: null, at: 200 },
    ]);
  });
});
