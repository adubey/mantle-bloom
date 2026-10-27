import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { publicWorldBase, withWorldId, worldBase, worldIdFromSearch } from "../src/worldRoutes.js";

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
