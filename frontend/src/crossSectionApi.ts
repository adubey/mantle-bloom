import { asJson } from "./apiResponse.js";

export interface CrossSectionResponse {
  // Raw PNG bytes, base64-encoded, following the /worlds/{world_id}/render image response convention.
  image_base64: string;
}

export function crossSectionUrl(
  apiBase: string,
  planeAngleDeg: number,
  positionIndex: number,
  width: number,
  height: number,
  showBoundaries = false,
  showFaults = false,
): string {
  const params = new URLSearchParams({
    plane_angle_deg: String(planeAngleDeg),
    position_index: String(positionIndex),
    width: String(width),
    height: String(height),
  });
  if (showBoundaries) params.set("show_boundaries", "true");
  if (showFaults) params.set("show_faults", "true");
  return `${apiBase}/cross-section?${params}`;
}

export function requestCrossSection(
  apiBase: string,
  planeAngleDeg: number,
  positionIndex: number,
  width: number,
  height: number,
  signal?: AbortSignal,
  fetcher: typeof fetch = fetch,
  showBoundaries = false,
  showFaults = false,
): Promise<CrossSectionResponse> {
  return fetcher(crossSectionUrl(apiBase, planeAngleDeg, positionIndex, width, height, showBoundaries, showFaults), { signal })
    .then(asJson<CrossSectionResponse>);
}
