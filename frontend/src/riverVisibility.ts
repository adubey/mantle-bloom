import type { RiverSummary } from "./api";

// Suppress short streams that add visual noise. This preserves the former 30-node cutoff at
// default node_density=4, but applies it in meters so it doesn't move with node density.
const MIN_RIVER_LENGTH_M = 1_875_000;

export const isVisibleRiver = (river: RiverSummary): boolean =>
  river.length_m >= MIN_RIVER_LENGTH_M;
