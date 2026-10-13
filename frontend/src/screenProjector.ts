import type { Projection } from "./api";
import type { Mat3, RenderTransform, Vec3 } from "./rotation";
import {
  getRenderTransform, latLonToXyz, matApply, matTranspose, project, toPixels, unproject, wrapLongitudeNear, xyzToLatLon,
} from "./rotation";

export type Pixel = [number, number];

// Same "skip segments longer than N times the median" seam-break technique as backend
// render_image._stroke_robust_loop.
const SEGMENT_BREAK_FACTOR = 6;

// Projects *true* world points onto the canvas under one (possibly in-progress) view rotation:
// the projection steps every client-drawn view (MapCanvas and the inspectors) repeats for each
// frame it draws.
export class ScreenProjector {
  readonly transform: RenderTransform;
  private readonly projection: Projection;
  private readonly rotation: Mat3;

  constructor(projection: Projection, rotation: Mat3, width: number, height: number) {
    this.projection = projection;
    this.rotation = rotation;
    this.transform = getRenderTransform(projection, width, height);
  }

  private toScreen(lat: number, lon: number): Pixel {
    const [x, y] = project(this.projection, lat, lon);
    return toPixels(this.transform, x, y);
  }

  // A lone point. Fine for independent dots, but two endpoints projected this way can bow a
  // short edge clear across the map at the antimeridian seam -- use segment()/loop() for
  // anything stroked or filled.
  point(p: Vec3): Pixel {
    const [lat, lon] = xyzToLatLon(matApply(this.rotation, p));
    return this.toScreen(lat, lon);
  }

  // One edge (two true world points a short 3D hop apart), unwrapping the second point's
  // longitude relative to the first's -- the same technique backend render_image.py's
  // _draw_rivers uses server-side (_project_offset).
  segment(a: Vec3, b: Vec3): [Pixel, Pixel] {
    const [latA, lonA] = xyzToLatLon(matApply(this.rotation, a));
    const [latB, rawLonB] = xyzToLatLon(matApply(this.rotation, b));
    return [this.toScreen(latA, lonA), this.toScreen(latB, wrapLongitudeNear(rawLonB, lonA))];
  }

  // A closed loop, unwrapping every point's longitude relative to the loop's own (rotated)
  // center -- see rotation.ts's wrapLongitudeNear for why a naive projection of a rotated
  // loop can bow a filled shape across the whole map.
  loop(points: Vec3[], center: Vec3): Pixel[] {
    const [, centerLon] = xyzToLatLon(matApply(this.rotation, center));
    return points.map((p) => {
      const [lat, lon] = xyzToLatLon(matApply(this.rotation, p));
      return this.toScreen(lat, wrapLongitudeNear(lon, centerLon));
    });
  }
}

// The true (un-rotated) [latDeg, lonDeg] under a click at backing-store pixel (x, y), or null
// for a click off the sphere's projected silhouette.
export function trueLatLonDegAt(
  projection: Projection, rotation: Mat3, width: number, height: number, backingX: number, backingY: number,
): [number, number] | null {
  const transform = getRenderTransform(projection, width, height);
  const x = (backingX - transform.offsetX) / transform.scale;
  const y = -(backingY - transform.offsetY) / transform.scale;
  const latLon = unproject(projection, x, y);
  if (!latLon) return null;
  const displayXyz = latLonToXyz(latLon[0], latLon[1]); // in the *display* (rotated) frame
  const trueXyz = matApply(matTranspose(rotation), displayXyz); // transpose == inverse for a rotation matrix
  const [trueLat, trueLon] = xyzToLatLon(trueXyz);
  return [(trueLat * 180) / Math.PI, (trueLon * 180) / Math.PI];
}

// Tab/Shift+Tab selection cycling: the id after (or before) `currentId` in ascending order,
// wrapping; the first (or last) id when nothing in `ids` is selected; null when `ids` is empty.
export function cycleId(ids: number[], currentId: number | null, backwards: boolean): number | null {
  if (ids.length === 0) return null;
  const sorted = [...ids].sort((a, b) => a - b);
  const currentIndex = currentId == null ? -1 : sorted.indexOf(currentId);
  const delta = backwards ? -1 : 1;
  const nextIndex = currentIndex === -1 ? (backwards ? sorted.length - 1 : 0) : (currentIndex + delta + sorted.length) % sorted.length;
  return sorted[nextIndex];
}

// Strokes a closed loop of projected pixels, skipping any segment much longer than the loop's
// median: a shape whose true angular extent exceeds what one reference-longitude unwrap can
// resolve still crosses the antimeridian seam, rarely, even after ScreenProjector.loop.
export function strokeRobustLoop(
  ctx: CanvasRenderingContext2D, pixels: Pixel[], color: string, lineWidth: number, canvasWidth: number,
) {
  const n = pixels.length;
  if (n < 2) return;
  const segLengths: number[] = [];
  for (let i = 0; i < n; i++) {
    const next = pixels[(i + 1) % n];
    segLengths.push(Math.hypot(next[0] - pixels[i][0], next[1] - pixels[i][1]));
  }
  const sorted = [...segLengths].sort((a, b) => a - b);
  const median = sorted[Math.floor(sorted.length / 2)] || 0;
  const breakThreshold = Math.max(median * SEGMENT_BREAK_FACTOR, 20 * (canvasWidth / 1100));
  ctx.strokeStyle = color;
  ctx.lineWidth = lineWidth;
  for (let i = 0; i < n; i++) {
    if (segLengths[i] > breakThreshold) continue;
    const next = pixels[(i + 1) % n];
    ctx.beginPath();
    ctx.moveTo(pixels[i][0], pixels[i][1]);
    ctx.lineTo(next[0], next[1]);
    ctx.stroke();
  }
}
