import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { test } from "node:test";
import { asJson, expandColumnarStatsHistory, type ColumnarStatsHistory } from "../src/apiResponse.js";

const RESPONSE_SIZE_BYTES = 15_000_000;

test("expands scalar columns and nested map columns into stats snapshots", () => {
  const payload: ColumnarStatsHistory = {
    format: "columnar-v1",
    length: 3,
    columns: {
      elapsed_years: [0, 10, 20],
      biome_land_fraction: {
        keys: ["forest", "desert"],
        columns: { forest: [0.8, null, 0.5], desert: [0.2, null, 0.5] },
        states: [2, 1, 2],
      },
      optional_map: {
        keys: [],
        columns: {},
        states: [0, 1, 0],
      },
    },
  };

  assert.deepEqual(expandColumnarStatsHistory(payload), [
    { elapsed_years: 0, biome_land_fraction: { forest: 0.8, desert: 0.2 }, optional_map: null },
    { elapsed_years: 10, biome_land_fraction: {}, optional_map: {} },
    { elapsed_years: 20, biome_land_fraction: { forest: 0.5, desert: 0.5 }, optional_map: null },
  ]);
});

test("parses and expands a 15 MB columnar JSON response", async (t) => {
  const rowCount = 150;
  const padding = "x".repeat(100_000);
  const payload: ColumnarStatsHistory = {
    format: "columnar-v1",
    length: rowCount,
    columns: {
      elapsed_years: Array.from({ length: rowCount }, (_, i) => i * 1_000_000),
      probe: Array(rowCount).fill(padding),
    },
  };
  const body = JSON.stringify(payload);
  const responseBytes = Buffer.byteLength(body);
  assert.ok(responseBytes >= RESPONSE_SIZE_BYTES);

  const startedAt = performance.now();
  const parsed = await asJson<ColumnarStatsHistory>(
    new Response(body, { status: 200, headers: { "Content-Type": "application/json" } }),
  );
  const history = expandColumnarStatsHistory(parsed);
  const elapsedMs = performance.now() - startedAt;

  assert.equal(history.length, rowCount);
  assert.equal(history[0].elapsed_years, 0);
  assert.equal((history.at(-1)?.probe as string).length, padding.length);
  t.diagnostic(`Parsed and expanded ${responseBytes.toLocaleString()} bytes in ${elapsedMs.toFixed(1)} ms`);
});

if (process.env.STATS_HISTORY_FIXTURE) {
  test("expands a captured stats-history response", async (t) => {
    const body = await readFile(process.env.STATS_HISTORY_FIXTURE!, "utf8");
    const parsed = await asJson<ColumnarStatsHistory>(
      new Response(body, { status: 200, headers: { "Content-Type": "application/json" } }),
    );
    const startedAt = performance.now();
    const history = expandColumnarStatsHistory(parsed);
    const elapsedMs = performance.now() - startedAt;

    assert.equal(history.length, parsed.length);
    assert.equal(history[0].elapsed_years, 0);
    assert.equal(history.at(-1)?.elapsed_years, 246_400_000);
    assert.ok(Object.keys(history[0].biome_land_fraction as object).length > 0);
    t.diagnostic(
      `Expanded captured ${Buffer.byteLength(body).toLocaleString()}-byte response with ${history.length.toLocaleString()} records in ${elapsedMs.toFixed(1)} ms`,
    );
  });
}
