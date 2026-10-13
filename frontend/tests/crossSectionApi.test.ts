import assert from "node:assert/strict";
import { test } from "node:test";
import { requestCrossSection } from "../src/crossSectionApi.js";

test("requests the selected cross-section frame with its dimensions and abort signal", async () => {
  const controller = new AbortController();
  let requestedUrl = "";
  let requestedSignal: AbortSignal | null | undefined;

  const response = await requestCrossSection(
    "http://localhost:8123",
    73,
    119,
    1100,
    611,
    controller.signal,
    async (input, init) => {
      requestedUrl = String(input);
      requestedSignal = init?.signal;
      return new Response(JSON.stringify({ image_base64: "cG5n" }), { status: 200 });
    },
  );

  const url = new URL(requestedUrl);
  assert.equal(url.origin, "http://localhost:8123");
  assert.equal(url.pathname, "/world/cross-section");
  assert.deepEqual(Object.fromEntries(url.searchParams), {
    plane_angle_deg: "73",
    position_index: "119",
    width: "1100",
    height: "611",
  });
  assert.equal(requestedSignal, controller.signal);
  assert.deepEqual(response, { image_base64: "cG5n" });
});

test("preserves the backend status when section rendering is unavailable", async () => {
  await assert.rejects(
    requestCrossSection("", 0, 0, 1100, 611, undefined, async () =>
      new Response("section rendering is not implemented", { status: 501, statusText: "Not Implemented" })),
    /501 Not Implemented/,
  );
});
