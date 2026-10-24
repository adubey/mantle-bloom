import { renderWorld } from "./api";
import type { WorldId } from "./worldRoutes";
import type { Vec3 } from "./rotation";
import { getRenderTransform, project, toPixels } from "./rotation";

// The small textured globes beside the map: the cross-section's locator and the sidebar's
// globe selector (#497). Both paint an orthographic disc client-side from one small elevation
// map, so dragging them repaints on every pointer move without a server round trip.

export interface MapPixels {
  data: Uint8ClampedArray;
  width: number;
  height: number;
}

// [center, right, up]: the true-frame unit vectors shown at the disc's centre, toward its
// right edge, and toward its top edge.
export type GlobeBasis = [Vec3, Vec3, Vec3];

const MAP_WIDTH = 720;
const MAP_HEIGHT = 360;

let cachedMap: { worldId: WorldId; worldRevision: number; pixels: Promise<MapPixels> } | null = null;

function decodeMap(dataUrl: string): Promise<MapPixels> {
  return new Promise((resolve, reject) => {
    const source = new Image();
    source.onload = () => {
      const canvas = document.createElement("canvas");
      canvas.width = source.width;
      canvas.height = source.height;
      const context = canvas.getContext("2d", { willReadFrequently: true });
      if (!context) {
        reject(new Error("2D canvas unavailable"));
        return;
      }
      context.drawImage(source, 0, 0);
      resolve({ data: context.getImageData(0, 0, source.width, source.height).data, width: source.width, height: source.height });
    };
    source.onerror = () => reject(new Error("Could not decode the globe map"));
    source.src = dataUrl;
  });
}

// The unrotated Behrmann elevation map, decoded once per (world, revision) and shared by
// both globes.
export function loadWorldMapPixels(worldId: WorldId, worldRevision: number): Promise<MapPixels> {
  if (cachedMap?.worldId === worldId && cachedMap?.worldRevision === worldRevision) return cachedMap.pixels;
  const pixels = renderWorld(worldId, "behrmann", "elevation", MAP_WIDTH, MAP_HEIGHT)
    .then((map) => decodeMap(`data:image/png;base64,${map.image_base64}`));
  const entry = { worldId, worldRevision, pixels };
  cachedMap = entry;
  pixels.catch(() => {
    if (cachedMap === entry) cachedMap = null;
  });
  return pixels;
}

// Fills `target` with the disc of radius `radius` centred on (cx, cy), lit toward the viewer;
// pixels off the disc are left untouched (transparent in a fresh ImageData).
export function paintGlobe(map: MapPixels, [center, right, up]: GlobeBasis, target: ImageData, cx: number, cy: number, radius: number): void {
  const transform = getRenderTransform("behrmann", map.width, map.height);
  const { data } = target;
  for (let y = 0; y < target.height; y++) {
    for (let x = 0; x < target.width; x++) {
      const sx = (x + 0.5 - cx) / radius;
      const sy = (cy - y - 0.5) / radius;
      const squared = sx * sx + sy * sy;
      if (squared > 1) continue;
      const depth = Math.sqrt(1 - squared);
      const vx = right[0] * sx + up[0] * sy + center[0] * depth;
      const vy = right[1] * sx + up[1] * sy + center[1] * depth;
      const vz = right[2] * sx + up[2] * sy + center[2] * depth;
      const lat = Math.asin(Math.max(-1, Math.min(1, vz)));
      const lon = Math.atan2(vy, vx);
      const [mapX, mapY] = project("behrmann", lat, lon);
      const [imageX, imageY] = toPixels(transform, mapX, mapY);
      const column = Math.max(0, Math.min(map.width - 1, Math.round(imageX)));
      const row = Math.max(0, Math.min(map.height - 1, Math.round(imageY)));
      const sourceIndex = (row * map.width + column) * 4;
      const targetIndex = (y * target.width + x) * 4;
      const shade = 0.72 + 0.28 * depth;
      data[targetIndex] = map.data[sourceIndex] * shade;
      data[targetIndex + 1] = map.data[sourceIndex + 1] * shade;
      data[targetIndex + 2] = map.data[sourceIndex + 2] * shade;
      data[targetIndex + 3] = 255;
    }
  }
}
