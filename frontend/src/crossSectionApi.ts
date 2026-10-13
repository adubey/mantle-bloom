import { asJson } from "./apiResponse.js";

export interface CrossSectionResponse {
  // Raw PNG bytes, base64-encoded, following the /world/render image response convention.
  image_base64: string;
}

export function crossSectionUrl(
  apiBase: string,
  planeAngleDeg: number,
  positionIndex: number,
  width: number,
  height: number,
): string {
  const params = new URLSearchParams({
    plane_angle_deg: String(planeAngleDeg),
    position_index: String(positionIndex),
    width: String(width),
    height: String(height),
  });
  return `${apiBase}/world/cross-section?${params}`;
}

export function requestCrossSection(
  apiBase: string,
  planeAngleDeg: number,
  positionIndex: number,
  width: number,
  height: number,
  signal?: AbortSignal,
  fetcher: typeof fetch = fetch,
): Promise<CrossSectionResponse> {
  return fetcher(crossSectionUrl(apiBase, planeAngleDeg, positionIndex, width, height), { signal })
    .then(asJson<CrossSectionResponse>);
}
