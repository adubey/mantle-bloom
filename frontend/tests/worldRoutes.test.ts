import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { ensureWorldId, latestWorldId, publicWorldBase, withWorldId, worldBase, worldIdFromSearch } from "../src/worldRoutes.js";

describe("worldRoutes", () => {
  it("builds admin and public bases with encoded IDs", () => {
    assert.equal(worldBase("", "abc 123"), "/worlds/abc%20123");
    assert.equal(publicWorldBase("http://x", "abc 123"), "http://x/public/worlds/abc%20123");
  });

  it("reads the ?world= selection from the page URL", () => {
    assert.equal(worldIdFromSearch("?world=abc&deb"), "abc");
    assert.equal(worldIdFromSearch("?deb"), null);
    assert.equal(worldIdFromSearch(""), null);
  });

  it("sets ?world= while preserving other params and the hash", () => {
    // URLSearchParams serializes the bare `deb` flag as `deb=`; the app's
    // `has("deb")` gate accepts both spellings.
    assert.equal(withWorldId("http://x/app?deb#frag", "abc"), "/app?deb=&world=abc#frag");
    assert.equal(withWorldId("http://x/app", "abc"), "/app?world=abc");
    assert.equal(withWorldId("http://x/app?world=old", "abc"), "/app?world=abc");
  });
});

describe("ensureWorldId", () => {
  it("keeps the current selection without creating a world", async () => {
    let created = 0;
    const id = await ensureWorldId("abc", async () => {
      created += 1;
      return "new-id";
    });
    assert.equal(id, "abc");
    assert.equal(created, 0);
  });

  it("creates the first world when nothing is selected yet (#592)", async () => {
    const id = await ensureWorldId(null, async () => "new-id");
    assert.equal(id, "new-id");
  });
});

describe("latestWorldId", () => {
  it("returns null for an empty world list", () => {
    assert.equal(latestWorldId([]), null);
  });

  it("picks the most recently updated world without reordering the list", () => {
    const worlds = [
      { world_id: "old", updated_at: 1 },
      { world_id: "new", updated_at: 3 },
      { world_id: "mid", updated_at: 2 },
    ];
    assert.equal(latestWorldId(worlds), "new");
    assert.deepEqual(
      worlds.map((w) => w.world_id),
      ["old", "new", "mid"],
    );
  });
});
