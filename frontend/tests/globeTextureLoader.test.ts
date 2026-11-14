import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { createGlobeTextureLoader } from "../src/globeTextureLoader.js";

function createFakeTimers() {
  let nextId = 1;
  const tasks = new Map<number, { at: number; fn: () => void }>();
  let now = 0;

  return {
    setTimer(fn: () => void, ms: number): ReturnType<typeof setTimeout> {
      const id = nextId++;
      tasks.set(id, { at: now + ms, fn });
      return id as unknown as ReturnType<typeof setTimeout>;
    },
    clearTimer(id: ReturnType<typeof setTimeout>) {
      tasks.delete(id as unknown as number);
    },
    async advanceBy(ms: number) {
      const target = now + ms;
      while (true) {
        let earliest: { id: number; at: number; fn: () => void } | null = null;
        for (const [id, task] of tasks) {
          if (task.at <= target && (!earliest || task.at < earliest.at)) {
            earliest = { id, ...task };
          }
        }
        if (!earliest) break;
        tasks.delete(earliest.id);
        now = earliest.at;
        earliest.fn();
        await Promise.resolve();
      }
      now = target;
    },
    pendingCount() {
      return tasks.size;
    },
  };
}

function flushMicrotasks() {
  return Promise.resolve().then(() => Promise.resolve());
}

describe("createGlobeTextureLoader", () => {
  it("skips revision-0 and missing worldId, then fetches once on first world load", async () => {
    const timers = createFakeTimers();
    const calls: Array<[string, number]> = [];
    const loaded: string[] = [];

    const loader = createGlobeTextureLoader<string>({
      loadPixels: async (worldId, rev) => {
        calls.push([worldId, rev]);
        return `${worldId}:${rev}`;
      },
      onLoaded: (pixels) => loaded.push(pixels),
      setTimer: timers.setTimer,
      clearTimer: timers.clearTimer,
    });

    loader.update({ worldId: "", worldRevision: 0, visible: true });
    loader.update({ worldId: "w1", worldRevision: 0, visible: true });
    await flushMicrotasks();
    assert.deepEqual(calls, []);

    loader.update({ worldId: "w1", worldRevision: 1, visible: true });
    await flushMicrotasks();
    assert.deepEqual(calls, [["w1", 1]]);
    assert.deepEqual(loaded, ["w1:1"]);
    loader.dispose();
  });

  it("issues no fetches while hidden and catches up lazily when shown", async () => {
    const timers = createFakeTimers();
    const calls: Array<[string, number]> = [];
    const loaded: string[] = [];

    const loader = createGlobeTextureLoader<string>({
      loadPixels: async (worldId, rev) => {
        calls.push([worldId, rev]);
        return `${worldId}:${rev}`;
      },
      onLoaded: (pixels) => loaded.push(pixels),
      debounceMs: 250,
      setTimer: timers.setTimer,
      clearTimer: timers.clearTimer,
    });

    loader.update({ worldId: "w1", worldRevision: 1, visible: true });
    await flushMicrotasks();
    assert.deepEqual(calls, [["w1", 1]]);

    // Step repeatedly while hidden: no renders should fire.
    loader.update({ worldId: "w1", worldRevision: 2, visible: false });
    loader.update({ worldId: "w1", worldRevision: 3, visible: false });
    loader.update({ worldId: "w1", worldRevision: 4, visible: false });
    await timers.advanceBy(1000);
    assert.deepEqual(calls, [["w1", 1]]);

    // Becoming visible again fetches the latest revision once.
    loader.update({ worldId: "w1", worldRevision: 4, visible: true });
    await timers.advanceBy(250);
    await flushMicrotasks();
    assert.deepEqual(calls, [
      ["w1", 1],
      ["w1", 4],
    ]);
    assert.deepEqual(loaded, ["w1:1", "w1:4"]);
    loader.dispose();
  });

  it("debounces rapid revision bumps while visible", async () => {
    const timers = createFakeTimers();
    const calls: Array<[string, number]> = [];

    const loader = createGlobeTextureLoader<string>({
      loadPixels: async (worldId, rev) => {
        calls.push([worldId, rev]);
        return `${worldId}:${rev}`;
      },
      onLoaded: () => {},
      debounceMs: 250,
      setTimer: timers.setTimer,
      clearTimer: timers.clearTimer,
    });

    loader.update({ worldId: "w1", worldRevision: 1, visible: true });
    await flushMicrotasks();

    loader.update({ worldId: "w1", worldRevision: 2, visible: true });
    await timers.advanceBy(100);
    loader.update({ worldId: "w1", worldRevision: 3, visible: true });
    await timers.advanceBy(100);
    loader.update({ worldId: "w1", worldRevision: 4, visible: true });
    await timers.advanceBy(249);
    assert.deepEqual(calls, [["w1", 1]]);

    await timers.advanceBy(1);
    await flushMicrotasks();
    assert.deepEqual(calls, [
      ["w1", 1],
      ["w1", 4],
    ]);
    loader.dispose();
  });

  it("retries with backoff after a failed fetch without a world change", async () => {
    const timers = createFakeTimers();
    const calls: Array<[string, number]> = [];
    const loaded: string[] = [];
    let failCount = 2;

    const loader = createGlobeTextureLoader<string>({
      loadPixels: async (worldId, rev) => {
        calls.push([worldId, rev]);
        if (failCount > 0) {
          failCount -= 1;
          throw new Error("503");
        }
        return `${worldId}:${rev}`;
      },
      onLoaded: (pixels) => loaded.push(pixels),
      retryBaseMs: 1000,
      retryMaxMs: 4000,
      setTimer: timers.setTimer,
      clearTimer: timers.clearTimer,
    });

    loader.update({ worldId: "w1", worldRevision: 1, visible: true });
    await flushMicrotasks();
    assert.deepEqual(calls, [["w1", 1]]);
    assert.deepEqual(loaded, []);

    // First retry after 1000ms (still fails).
    await timers.advanceBy(1000);
    await flushMicrotasks();
    assert.equal(calls.length, 2);
    assert.deepEqual(loaded, []);

    // Second retry after 2000ms (succeeds).
    await timers.advanceBy(1999);
    assert.equal(calls.length, 2);
    await timers.advanceBy(1);
    await flushMicrotasks();
    assert.equal(calls.length, 3);
    assert.deepEqual(loaded, ["w1:1"]);
    assert.equal(timers.pendingCount(), 0);
    loader.dispose();
  });

  it("stops retrying after maxRetries and clears only when leaving a non-empty world", async () => {
    const timers = createFakeTimers();
    let calls = 0;
    let cleared = 0;
    const loaded: string[] = [];

    const loader = createGlobeTextureLoader<string>({
      loadPixels: async (worldId, rev) => {
        calls += 1;
        if (worldId === "w1") throw new Error("404");
        return `${worldId}:${rev}`;
      },
      onLoaded: (pixels) => loaded.push(pixels),
      onCleared: () => {
        cleared += 1;
      },
      retryBaseMs: 1000,
      retryMaxMs: 4000,
      maxRetries: 2,
      setTimer: timers.setTimer,
      clearTimer: timers.clearTimer,
    });

    loader.update({ worldId: "w1", worldRevision: 1, visible: true });
    await flushMicrotasks();
    assert.equal(cleared, 0); // "" -> "w1" does not fire onCleared
    assert.equal(calls, 1);

    await timers.advanceBy(1000);
    await flushMicrotasks();
    assert.equal(calls, 2);

    await timers.advanceBy(2000);
    await flushMicrotasks();
    assert.equal(calls, 3);

    // No further retries scheduled after maxRetries (2).
    assert.equal(timers.pendingCount(), 0);
    await timers.advanceBy(10000);
    assert.equal(calls, 3);

    loader.update({ worldId: "w2", worldRevision: 1, visible: true });
    await flushMicrotasks();
    assert.equal(cleared, 1); // "w1" -> "w2" fires onCleared once
    assert.deepEqual(loaded, ["w2:1"]);
    loader.dispose();
  });
});
