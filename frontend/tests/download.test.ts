import { resetDom } from "./domSetup";
import assert from "node:assert/strict";
import { afterEach, beforeEach, describe, it, mock } from "node:test";
import { saveRecordedVideo } from "../src/download";

// A few bytes that aren't valid UTF-8, so a text round trip would corrupt them.
const VIDEO_BYTES = new Uint8Array([0x00, 0x00, 0x00, 0x20, 0x66, 0x74, 0x79, 0x70, 0xff, 0xfe, 0x80]);
const VIDEO_BASE64 = Buffer.from(VIDEO_BYTES).toString("base64");

describe("saveRecordedVideo", () => {
  let created: Blob[];
  let revoked: string[];
  let clicks: { href: string; download: string }[];
  let onClick: (event: Event) => void;

  beforeEach(() => {
    resetDom();
    created = [];
    revoked = [];
    clicks = [];
    mock.timers.enable({ apis: ["setTimeout"] });
    // The served app's `connect-src 'self'` refuses a fetch of a data: URL (#961).
    mock.method(globalThis, "fetch", async (input: unknown) => {
      if (String(input).startsWith("data:")) throw new TypeError("Failed to fetch");
      throw new Error(`unexpected fetch of ${String(input)}`);
    });
    mock.method(URL, "createObjectURL", (blob: Blob) => {
      created.push(blob);
      return `blob:http://localhost/${created.length}`;
    });
    mock.method(URL, "revokeObjectURL", (url: string) => {
      revoked.push(url);
    });
    onClick = (event: Event) => {
      const anchor = event.target as HTMLAnchorElement;
      clicks.push({ href: anchor.href, download: anchor.download });
      event.preventDefault();
    };
    document.addEventListener("click", onClick);
  });

  afterEach(() => {
    document.removeEventListener("click", onClick);
    mock.timers.reset();
    mock.restoreAll();
  });

  it("downloads the decoded video without fetching a data: URL", async () => {
    saveRecordedVideo({ videoBase64: VIDEO_BASE64, mime: "video/mp4", seed: 42, elapsed_years: 2_500_000.4 });

    assert.equal(created.length, 1);
    assert.equal(created[0].type, "video/mp4");
    assert.deepEqual(new Uint8Array(await created[0].arrayBuffer()), VIDEO_BYTES);
    assert.deepEqual(clicks, [{ href: "blob:http://localhost/1", download: "mantle-bloom-seed42-animation-2500000y.mp4" }]);
    assert.equal(document.querySelectorAll("a").length, 0);
  });

  it("keeps the object URL alive until the browser has started the download", () => {
    saveRecordedVideo({ videoBase64: VIDEO_BASE64, mime: "video/mp4", seed: 1, elapsed_years: 0 });

    assert.deepEqual(revoked, []);
    mock.timers.tick(60_000);
    assert.deepEqual(revoked, ["blob:http://localhost/1"]);
  });

  it("throws on a payload it can't decode, so the caller keeps the video", () => {
    assert.throws(() => saveRecordedVideo({ videoBase64: "not base64!", mime: "video/mp4", seed: 1, elapsed_years: 0 }));
    assert.deepEqual(clicks, []);
  });
});
