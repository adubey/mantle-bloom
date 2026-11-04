import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { existsSync, readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { groupName, PLAYING_FALLBACK_PHASE, phaseLinks, playingPhaseView } from "../src/phaseInfo.js";

// Every phase label the backend streams (backend app/world.py #629) must resolve to at
// least one working learn link, all of them https. Unknown phases get no box rather
// than a wrong link.
const BACKEND_PHASES = [
  "Creating continents",
  "Creating cratons",
  "Starting mantle convection",
  "Seeding climate and oceans",
  "Moving plates",
  "Reshaping plates",
  "Merging and splitting plates",
  "Filling seafloor gaps",
  "Moving magma",
  "Raising mantle plumes",
  "Resolving glaciers",
  "Raising volcanoes",
  "Forming natural resources",
  "Balancing sea level",
  // The user-visible spans (backend app/trace.py, #801).
  "Laying out continents and mountain belts",
  "Dividing the sphere into plates",
  "Building pseudo-history",
  "Shaping initial bathymetry",
  "Dating the initial rocks",
  "Seeding rivers and glaciers",
  "Solving mantle flow and plate motion",
  "Detecting collisions",
  "Updating faults",
  "Solving wind and ocean currents",
  "Resolving climate",
  "Resolving rivers and lakes",
  "Marching glaciers",
  "Eroding and depositing sediment",
];

describe("phaseLinks", () => {
  it("resolves every backend phase to at least one https link", () => {
    for (const phase of BACKEND_PHASES) {
      const links = phaseLinks(phase);
      assert.ok(links.length > 0, `${phase} has no learn links`);
      for (const link of links) {
        assert.match(link.url, /^https:\/\//, `${phase} link is not https: ${link.url}`);
      }
    }
  });

  it("returns no links for unknown or missing phases", () => {
    assert.deepEqual(phaseLinks("Reticulating splines"), []);
    assert.deepEqual(phaseLinks(null), []);
    assert.deepEqual(phaseLinks(undefined), []);
  });
});

describe("groupName", () => {
  it("drops the bar from a span group and leaves a bare group alone", () => {
    assert.equal(groupName("Stepping/Hydrometeorology (climate, water & ice)"), "Hydrometeorology (climate, water & ice)");
    assert.equal(groupName("Moving plates"), "Moving plates");
  });
});

describe("playingPhaseView", () => {
  it("shows the runner's phase while a step is under way", () => {
    assert.deepEqual(playingPhaseView({ label: "Updating faults", group: "Stepping/Moving plates", fraction: 0.22 }), {
      phase: "Updating faults",
      group: "Stepping/Moving plates",
      fraction: 0.22,
    });
  });

  it("falls back to the generic line between steps, not the last phase", () => {
    // Regression for #809 review: a null phase used to keep "Balancing sea level" at ~97%.
    assert.deepEqual(playingPhaseView(null), { phase: PLAYING_FALLBACK_PHASE, group: null, fraction: undefined });
  });
});

// GitHub's heading anchor: lowercase, drop everything but letters, digits, spaces, hyphens and
// underscores, then spaces become hyphens.
function githubAnchor(heading: string): string {
  return heading
    .trim()
    .toLowerCase()
    .replace(/[^a-z0-9 _-]/g, "")
    .replace(/ /g, "-");
}

describe("docs links", () => {
  // phaseInfo.ts links into docs/simulation-model on GitHub; check each file and anchor against the
  // local docs, so a renamed heading fails here instead of showing a dead link in the overlay.
  const repoRoot = join(dirname(fileURLToPath(import.meta.url)), "..", "..", "..");
  const prefix = "https://github.com/adubey/mantle-bloom/blob/main/";

  it("points every docs link at a file and heading that exist", () => {
    const urls = new Set<string>();
    for (const phase of BACKEND_PHASES) for (const link of phaseLinks(phase)) if (link.url.startsWith(prefix)) urls.add(link.url);
    assert.ok(urls.size > 10, "found no docs links to check");
    for (const url of urls) {
      const [path, anchor] = url.slice(prefix.length).split("#");
      const file = join(repoRoot, path);
      assert.ok(existsSync(file), `${url}: ${path} does not exist`);
      if (anchor === undefined) continue;
      const anchors = readFileSync(file, "utf8")
        .split("\n")
        .filter((line) => /^#{1,6} /.test(line))
        .map((line) => githubAnchor(line.replace(/^#{1,6} /, "")));
      assert.ok(anchors.includes(anchor), `${url}: no heading with anchor #${anchor}`);
    }
  });
});
