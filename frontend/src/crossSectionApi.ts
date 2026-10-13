import { asJson } from "./apiResponse.js";

export interface CrossSectionResponse {
  // Raw PNG bytes, base64-encoded, following the /world/render image response convention.
  image_base64: string;
}

export function crossSectionUrl(
  apiBase: string,
  latitudeDeg: number,
  longitudeIndex: number,
  width: number,
  height: number,
): string {
  const params = new URLSearchParams({
    latitude_deg: String(latitudeDeg),
    longitude_index: String(longitudeIndex),
    width: String(width),
    height: String(height),
  });
  return `${apiBase}/world/cross-section?${params}`;
}

export function requestCrossSection(
  apiBase: string,
  latitudeDeg: number,
  longitudeIndex: number,
  width: number,
  height: number,
  signal?: AbortSignal,
  fetcher: typeof fetch = fetch,
): Promise<CrossSectionResponse> {
  return fetcher(crossSectionUrl(apiBase, latitudeDeg, longitudeIndex, width, height), { signal })
    .then(asJson<CrossSectionResponse>);
}
