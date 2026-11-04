import assert from "node:assert/strict";
import { describe, it } from "node:test";
import { groupName, phaseLinks } from "../src/phaseInfo.js";

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
