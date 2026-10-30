import assert from "node:assert/strict";
import { test } from "node:test";
import { crossSectionUrl, requestCrossSection } from "../src/crossSectionApi.js";
import { advanceCrossSectionPosition, nearbyCrossSectionPositions } from "../src/crossSectionFrames.js";
import { projectGreatCirclePoint } from "../src/crossSectionLocator.js";

test("requests the selected cross-section frame with its dimensions and abort signal", async () => {
  const controller = new AbortController();
  let requestedUrl = "";
  let requestedSignal: AbortSignal | null | undefined;

  const response = await requestCrossSection(
    "http://localhost:8123/worlds/abc123",
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
  assert.equal(url.pathname, "/worlds/abc123/cross-section");
  assert.deepEqual(Object.fromEntries(url.searchParams), {
    plane_angle_deg: "73",
    position_index: "119",
    width: "1100",
    height: "611",
  });
  assert.equal(requestedSignal, controller.signal);
  assert.deepEqual(response, { image_base64: "cG5n" });
});

test("encodes the plumes flag only when enabled", () => {
  const without = crossSectionUrl("http://x", 10, 2, 320, 200);
  assert.equal(without.includes("show_plumes"), false);

  const withPlumes = crossSectionUrl("http://x", 10, 2, 320, 200, false, false, true);
  assert.equal(withPlumes.includes("show_plumes=true"), true);
});

test("preserves the backend status when section rendering is unavailable", async () => {
  await assert.rejects(
    requestCrossSection("", 0, 0, 1100, 611, undefined, async () =>
      new Response("section rendering is not implemented", { status: 501, statusText: "Not Implemented" })),
    /501 Not Implemented/,
  );
});

test("selects only nearby prefetch frames and wraps around the great circle", () => {
  assert.deepEqual(nearbyCrossSectionPositions(119, 120, 8), [118, 0, 117, 1, 116, 2, 115, 3]);
  assert.deepEqual(nearbyCrossSectionPositions(4, 8, 20), [3, 5, 2, 6, 1, 7, 0]);
  const remainingFrames = nearbyCrossSectionPositions(37, 120, 119);
  assert.equal(remainingFrames.length, 119);
  assert.equal(new Set(remainingFrames).size, 119);
  assert.equal(remainingFrames.includes(37), false);
});

test("fast arc keyboard steps wrap across two positions", () => {
  assert.equal(advanceCrossSectionPosition(0, 120, -1, 2), 118);
  assert.equal(advanceCrossSectionPosition(119, 120, 1, 2), 1);
});

test("locator shows the equator horizontally and oblique planes tilted", () => {
  const equator = projectGreatCirclePoint(0, Math.PI / 2, 130, 130, 106);
  assert.ok(Math.abs(equator.x - 130) < 1e-10);
  assert.ok(equator.y > 130);

  const meridian = projectGreatCirclePoint(90, 0, 130, 130, 106);
  assert.ok(Math.abs(meridian.x - 130) < 1e-10);
  assert.ok(meridian.y < 130);

  const oblique = projectGreatCirclePoint(45, 0, 130, 130, 106);
  assert.ok(oblique.x > 130);
  assert.ok(oblique.y < 130);
});
