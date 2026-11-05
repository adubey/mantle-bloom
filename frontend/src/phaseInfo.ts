// Learn-more links for the run overlay's phase text (#629): one entry per phase label
// the backend streams can send (see backend app/world.py's generate_world_progress and
// step_world_progress, and the user-visible spans in app/trace.py). Wikipedia pages were checked
// to return 200 when added (all 23, 2026-11-04); tests/phaseInfo.test.ts checks every docs link's
// file and heading offline. Unknown phases get no box rather than a wrong link.

export interface PhaseLink {
  label: string;
  url: string;
}

const DOCS = "https://github.com/adubey/mantle-bloom/blob/main/docs/simulation-model";
const WIKI = "https://en.wikipedia.org/wiki";

const PHASE_LINKS: Record<string, PhaseLink[]> = {
  "Creating continents": [
    { label: "Continent — Wikipedia", url: `${WIKI}/Continent` },
    { label: "How worlds are generated", url: `${DOCS}/foundations.md#initial-plate-generation-platespy` },
  ],
  "Creating cratons": [
    { label: "Craton — Wikipedia", url: `${WIKI}/Craton` },
    { label: "Cratons in the model", url: `${DOCS}/tectonics.md#cratons` },
  ],
  "Starting mantle convection": [
    { label: "Mantle convection — Wikipedia", url: `${WIKI}/Mantle_convection` },
    { label: "Mantle flow in the model", url: `${DOCS}/foundations.md#mantle-flow-mantlepy` },
  ],
  "Seeding climate and oceans": [
    { label: "Climate — Wikipedia", url: `${WIKI}/Climate` },
    { label: "Climate in the model", url: `${DOCS}/climate.md#climate-climatepy` },
  ],
  "Moving plates": [
    { label: "Plate tectonics — Wikipedia", url: `${WIKI}/Plate_tectonics` },
    { label: "How a plate changes each step", url: `${DOCS}/tectonics.md#how-a-plate-changes-during-each-time-step` },
  ],
  "Reshaping plates": [
    { label: "Plate tectonics — Wikipedia", url: `${WIKI}/Plate_tectonics` },
    { label: "How a plate changes each step", url: `${DOCS}/tectonics.md#how-a-plate-changes-during-each-time-step` },
  ],
  "Merging and splitting plates": [
    { label: "Plate tectonics — Wikipedia", url: `${WIKI}/Plate_tectonics` },
    { label: "Merge and split in the model", url: `${DOCS}/tectonics.md#merge-and-split` },
  ],
  "Filling seafloor gaps": [
    { label: "Seafloor spreading — Wikipedia", url: `${WIKI}/Seafloor_spreading` },
    { label: "Gap filling in the model", url: `${DOCS}/tectonics.md#gap-filling` },
  ],
  "Moving magma": [{ label: "Magma — Wikipedia", url: `${WIKI}/Magma` }],
  "Raising mantle plumes": [
    { label: "Mantle plume — Wikipedia", url: `${WIKI}/Mantle_plume` },
    { label: "Mantle plumes in the model", url: `${DOCS}/tectonics.md#mantle-plumes` },
  ],
  "Resolving glaciers": [
    { label: "Glacier — Wikipedia", url: `${WIKI}/Glacier` },
    { label: "Glaciation in the model", url: `${DOCS}/hydrology.md#glaciation-hydrologypy` },
  ],
  "Raising volcanoes": [
    { label: "Volcano — Wikipedia", url: `${WIKI}/Volcano` },
    { label: "Volcanism in the model", url: `${DOCS}/tectonics.md#volcanism` },
  ],
  "Forming natural resources": [
    { label: "Ore — Wikipedia", url: `${WIKI}/Ore` },
    { label: "Surface processes and resources", url: `${DOCS}/surface-processes.md` },
  ],
  // Span-driven labels (#801): the same names as backend app/trace.py's user-visible spans.
  "Laying out continents and mountain belts": [
    { label: "Continent — Wikipedia", url: `${WIKI}/Continent` },
    { label: "Orogeny — Wikipedia", url: `${WIKI}/Orogeny` },
    { label: "How worlds are generated", url: `${DOCS}/foundations.md#initial-plate-generation-platespy` },
  ],
  "Dividing the sphere into plates": [
    { label: "Voronoi diagram — Wikipedia", url: `${WIKI}/Voronoi_diagram` },
    { label: "How worlds are generated", url: `${DOCS}/foundations.md#initial-plate-generation-platespy` },
  ],
  "Building pseudo-history": [
    { label: "Relief from a synthetic tectonic history", url: `${DOCS}/foundations.md#relief-from-a-synthetic-tectonic-history-synthetic_reliefpy` },
  ],
  "Shaping initial bathymetry": [{ label: "Bathymetry — Wikipedia", url: `${WIKI}/Bathymetry` }],
  "Dating the initial rocks": [{ label: "Radiometric dating — Wikipedia", url: `${WIKI}/Radiometric_dating` }],
  "Seeding rivers and glaciers": [
    { label: "River — Wikipedia", url: `${WIKI}/River` },
    { label: "Glacier — Wikipedia", url: `${WIKI}/Glacier` },
    { label: "Rivers and lakes in the model", url: `${DOCS}/hydrology.md#hydrology-rivers-and-lakes-hydrologypy` },
  ],
  "Solving mantle flow and plate motion": [
    { label: "Mantle convection — Wikipedia", url: `${WIKI}/Mantle_convection` },
    { label: "Plate tectonics — Wikipedia", url: `${WIKI}/Plate_tectonics` },
    { label: "Mantle flow in the model", url: `${DOCS}/foundations.md#mantle-flow-mantlepy` },
  ],
  "Detecting collisions": [
    { label: "Convergent boundary — Wikipedia", url: `${WIKI}/Convergent_boundary` },
    { label: "How a plate changes each step", url: `${DOCS}/tectonics.md#how-a-plate-changes-during-each-time-step` },
  ],
  "Updating faults": [
    { label: "Fault — Wikipedia", url: `${WIKI}/Fault_(geology)` },
    { label: "Faults in the model", url: `${DOCS}/tectonics.md#faults` },
  ],
  "Solving wind and ocean currents": [
    { label: "Atmospheric circulation — Wikipedia", url: `${WIKI}/Atmospheric_circulation` },
    { label: "Ocean current — Wikipedia", url: `${WIKI}/Ocean_current` },
    { label: "Atmospheric fluid dynamics in the model", url: `${DOCS}/climate.md#atmospheric-fluid-dynamics-atmosphere_cfdpy` },
  ],
  "Resolving climate": [
    { label: "Climate — Wikipedia", url: `${WIKI}/Climate` },
    { label: "Climate in the model", url: `${DOCS}/climate.md#climate-climatepy` },
  ],
  "Resolving rivers and lakes": [
    { label: "River — Wikipedia", url: `${WIKI}/River` },
    { label: "Lake — Wikipedia", url: `${WIKI}/Lake` },
    { label: "Rivers and lakes in the model", url: `${DOCS}/hydrology.md#hydrology-rivers-and-lakes-hydrologypy` },
  ],
  "Marching glaciers": [
    { label: "Glacier — Wikipedia", url: `${WIKI}/Glacier` },
    { label: "Glaciation in the model", url: `${DOCS}/hydrology.md#glaciation-hydrologypy` },
  ],
  "Eroding and depositing sediment": [
    { label: "Erosion — Wikipedia", url: `${WIKI}/Erosion` },
    { label: "Erosion in the model", url: `${DOCS}/surface-processes.md#erosion-erosionpy` },
  ],
  "Balancing sea level": [
    { label: "Sea level — Wikipedia", url: `${WIKI}/Sea_level` },
    { label: "Eustatic sea level in the model", url: `${DOCS}/surface-processes.md#eustatic-sea-level-eustasypy` },
  ],
};

// The text of a "<bar>/<group>" label from the backend's spans (#801): the group alone, since the
// overlay's headline already names the bar.
export function groupName(group: string): string {
  const slash = group.indexOf("/");
  return slash === -1 ? group : group.slice(slash + 1);
}

// What the overlay shows for a playing world, from the runner's published phase (#801). Between
// steps the phase is null; showing the last one would leave the overlay naming finished work near
// 100%, so fall back to the generic line with no group and an indeterminate bar.
export const PLAYING_FALLBACK_PHASE = "Advancing the world";

export interface PlayingPhaseView {
  phase: string;
  group: string | null;
  fraction: number | undefined;
}

export function playingPhaseView(phase: { label: string; group: string | null; fraction: number } | null): PlayingPhaseView {
  if (!phase) return { phase: PLAYING_FALLBACK_PHASE, group: null, fraction: undefined };
  return { phase: phase.label, group: phase.group, fraction: phase.fraction };
}

export function phaseLinks(phase: string | null | undefined): PhaseLink[] {
  if (!phase) return [];
  return PHASE_LINKS[phase] ?? [];
}
