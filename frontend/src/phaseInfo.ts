// Learn-more links for the run overlay's phase text (#629): one entry per phase label
// the backend streams can send (see backend app/world.py's generate_world_progress and
// step_world_progress). Every URL here was verified to resolve when added; unknown
// phases get no box rather than a wrong link.

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
  "Balancing sea level": [
    { label: "Sea level — Wikipedia", url: `${WIKI}/Sea_level` },
    { label: "Eustatic sea level in the model", url: `${DOCS}/surface-processes.md#eustatic-sea-level-eustasypy` },
  ],
};

export function phaseLinks(phase: string | null | undefined): PhaseLink[] {
  if (!phase) return [];
  return PHASE_LINKS[phase] ?? [];
}
