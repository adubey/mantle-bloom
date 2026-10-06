import type { RiverSummary } from "./api";

// Suppress tiny, unbranched streams that add visual noise without representing a river
// network. Branched rivers remain visible regardless of length.
export const isVisibleRiver = (river: RiverSummary): boolean =>
  river.num_tributaries > 0 || river.num_nodes >= 30;
