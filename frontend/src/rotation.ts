// Client-side math for the "rotate the planet" feature: just enough of geometry.py/
// projections.py/render_image.py ported to TypeScript to (a) drive the arcball drag gesture
// and (b) draw a cheap wireframe graticule preview while dragging, without needing a round
// trip to the server for every mouse move (the real detailed map -- elevation fill, plate
// boundaries, climate heatmaps -- is far too slow to re-render at that rate; see
// docs/simulation-model.md#rotating-the-view). The final committed rotation is sent to
// `/world/render`'s `rotation` query param (see api.ts's renderWorld) for the real render.

import type { Projection } from "./api";

export type Vec3 = [number, number, number];
// Row-major 3x3, flattened -- matches the wire format render.py's _parse_view_rotation
// expects (9 comma-separated floats).
export type Mat3 = number[];

export const IDENTITY_ROTATION: Mat3 = [1, 0, 0, 0, 1, 0, 0, 0, 1];

export function latLonToXyz(latRad: number, lonRad: number): Vec3 {
  const cosLat = Math.cos(latRad);
  return [cosLat * Math.cos(lonRad), cosLat * Math.sin(lonRad), Math.sin(latRad)];
}

export function xyzToLatLon([x, y, z]: Vec3): [number, number] {
  const lat = Math.asin(Math.min(1, Math.max(-1, z)));
  const lon = Math.atan2(y, x);
  return [lat, lon];
}

export function vecLength([x, y, z]: Vec3): number {
  return Math.hypot(x, y, z);
}

export function vecNormalize(v: Vec3): Vec3 {
  const len = vecLength(v);
  if (len < 1e-12) return [0, 0, 1];
  return [v[0] / len, v[1] / len, v[2] / len];
}

export function vecDot(a: Vec3, b: Vec3): number {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

export function vecCross(a: Vec3, b: Vec3): Vec3 {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}

export function matIdentity(): Mat3 {
  return [...IDENTITY_ROTATION];
}

// b applied first, then a -- i.e. the combined matrix represents "rotate by b, then by a",
// matching standard matrix-multiplication composition order (a @ b).
export function matMultiply(a: Mat3, b: Mat3): Mat3 {
  const out = new Array(9).fill(0);
  for (let r = 0; r < 3; r++) {
    for (let c = 0; c < 3; c++) {
      let sum = 0;
      for (let k = 0; k < 3; k++) sum += a[r * 3 + k] * b[k * 3 + c];
      out[r * 3 + c] = sum;
    }
  }
  return out;
}

export function matApply(m: Mat3, [x, y, z]: Vec3): Vec3 {
  return [
    m[0] * x + m[1] * y + m[2] * z,
    m[3] * x + m[4] * y + m[5] * z,
    m[6] * x + m[7] * y + m[8] * z,
  ];
}

// Rotation matrices are orthonormal, so transpose == inverse -- used to go from a click's
// world position in the *display* (rotated) frame back to the *true* frame (see
// PlateInspector.tsx's click handling).
export function matTranspose(m: Mat3): Mat3 {
  return [m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]];
}

// Rodrigues' rotation formula -- same construction as backend geometry.py's rotation_matrix,
// ported directly so the two stay in exact agreement.
export function rotationMatrix(axis: Vec3, angleRad: number): Mat3 {
  const n = vecLength(axis);
  if (n < 1e-12 || Math.abs(angleRad) < 1e-12) return matIdentity();
  const [ax, ay, az] = [axis[0] / n, axis[1] / n, axis[2] / n];
  const sin = Math.sin(angleRad);
  const cos = Math.cos(angleRad);
  const t = 1 - cos;
  return [
    t * ax * ax + cos, t * ax * ay - sin * az, t * ax * az + sin * ay,
    t * ax * ay + sin * az, t * ay * ay + cos, t * ay * az - sin * ax,
    t * ax * az - sin * ay, t * ay * az + sin * ax, t * az * az + cos,
  ];
}

// The rotation that carries unit vector `from` onto unit vector `to`, via the shortest arc.
export function rotationBetween(from: Vec3, to: Vec3): Mat3 {
  const dot = Math.min(1, Math.max(-1, vecDot(from, to)));
  if (dot > 1 - 1e-12) return matIdentity(); // already (almost) the same point
  if (dot < -1 + 1e-12) {
    // Antipodal: any axis perpendicular to `from` gives a valid 180-degree rotation.
    const helper: Vec3 = Math.abs(from[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0];
    const axis = vecNormalize(vecCross(from, helper));
    return rotationMatrix(axis, Math.PI);
  }
  const axis = vecNormalize(vecCross(from, to));
  return rotationMatrix(axis, Math.acos(dot));
}

// --- Arcball: maps a canvas-space point to a point on a virtual unit trackball centered on
// the canvas, the standard Shoemake (1992) construction. Points outside the ball's radius
// are clamped to its silhouette (z = 0) rather than left undefined, so a drag that leaves the
// nominal circle still produces a sensible (increasingly edge-on, twisting) rotation instead
// of breaking.
export function screenToArcballVector(px: number, py: number, centerX: number, centerY: number, radius: number): Vec3 {
  let x = (px - centerX) / radius;
  let y = -(py - centerY) / radius; // screen y grows down; the trackball's y should grow up
  const lenSq = x * x + y * y;
  if (lenSq > 1) {
    const scale = 1 / Math.sqrt(lenSq);
    x *= scale;
    y *= scale;
    return [x, y, 0];
  }
  return [x, y, Math.sqrt(1 - lenSq)];
}

// --- Projections: direct ports of backend/app/projections.py, same formulas, same units
// (radians in, planar units out for a unit-radius sphere).

const BEHRMANN_STANDARD_PARALLEL = (30 * Math.PI) / 180;
const BEHRMANN_COS_PHI0 = Math.cos(BEHRMANN_STANDARD_PARALLEL);

function behrmann(lat: number, lon: number): [number, number] {
  return [lon * BEHRMANN_COS_PHI0, Math.sin(lat) / BEHRMANN_COS_PHI0];
}

const ECKERT4_C = 2 / Math.sqrt(Math.PI * (4 + Math.PI));
const ECKERT4_RHS_SCALE = 2 + Math.PI / 2;
const ECKERT4_NEWTON_ITERS = 30;

function eckert4Theta(lat: number): number {
  const target = ECKERT4_RHS_SCALE * Math.sin(lat);
  if (Math.abs(Math.abs(lat) - Math.PI / 2) < 1e-9) return Math.sign(lat) * (Math.PI / 2);
  let theta = lat / 2;
  for (let i = 0; i < ECKERT4_NEWTON_ITERS; i++) {
    const sinT = Math.sin(theta);
    const cosT = Math.cos(theta);
    const f = theta + sinT * cosT + 2 * sinT - target;
    const fPrime = 2 * cosT * (cosT + 1);
    theta -= f / (Math.abs(fPrime) < 1e-9 ? 1e-9 : fPrime);
  }
  return theta;
}

function eckert4(lat: number, lon: number): [number, number] {
  const theta = eckert4Theta(lat);
  return [ECKERT4_C * lon * (1 + Math.cos(theta)), ECKERT4_C * Math.PI * Math.sin(theta)];
}

export function project(projection: Projection, lat: number, lon: number): [number, number] {
  return projection === "behrmann" ? behrmann(lat, lon) : eckert4(lat, lon);
}

// --- Inverse projections: for turning a click's planar (x, y) back into (lat, lon), so
// PlateInspector can unproject a click and ask the server which plate is there
// (/world/plate_at). Both return null outside their valid domain (e.g. a click landed in the
// canvas's background padding, off the sphere's own projected silhouette) -- callers must
// skip the plate_at call on null, not coerce it to NaN.

function inverseBehrmann(x: number, y: number): [number, number] | null {
  const sinLat = y * BEHRMANN_COS_PHI0;
  if (sinLat < -1 || sinLat > 1) return null;
  const lon = x / BEHRMANN_COS_PHI0;
  if (lon < -Math.PI || lon > Math.PI) return null;
  return [Math.asin(sinLat), lon];
}

// Closed-form in this direction even though the forward direction needs Newton iteration:
// y = C*pi*sin(theta) gives theta directly via one asin; plugging that theta into the same
// defining equation (theta + sin(theta)cos(theta) + 2sin(theta) = (2+pi/2)*sin(lat)) then
// solves lat via one more asin, since lat only ever appears there inside a sin.
function inverseEckert4(x: number, y: number): [number, number] | null {
  const sinTheta = y / (ECKERT4_C * Math.PI);
  if (sinTheta < -1 || sinTheta > 1) return null;
  const theta = Math.asin(sinTheta);
  const sinLat = (theta + Math.sin(theta) * Math.cos(theta) + 2 * Math.sin(theta)) / ECKERT4_RHS_SCALE;
  if (sinLat < -1 || sinLat > 1) return null;
  const denom = ECKERT4_C * (1 + Math.cos(theta));
  if (Math.abs(denom) < 1e-9) return null; // exactly at a pole (theta=+-pi/2): lon undefined
  return [Math.asin(sinLat), x / denom];
}

export function unproject(projection: Projection, x: number, y: number): [number, number] | null {
  return projection === "behrmann" ? inverseBehrmann(x, y) : inverseEckert4(x, y);
}

// --- The render transform (scale/offset from a bounding-box fit) -- see render_image.py's
// render_png/_render_climate_view. Deliberately kept in exact sync (REFERENCE_WIDTH_PX,
// PADDING_PX) since this needs to reproduce the server's own transform pixel-for-pixel for
// the graticule preview to align with the last real rendered frame underneath it.
const REFERENCE_WIDTH_PX = 1100;
const PADDING_PX = 20;

export interface RenderTransform {
  scale: number;
  offsetX: number;
  offsetY: number;
}

// The whole sphere's projected bounding box is identical regardless of view rotation (see
// docs/simulation-model.md#rotating-the-view) -- rotation only permutes which physical point
// lands at which lat/lon, never changes the set of lat/lon values a full sphere covers. So
// this only needs `projection`/`width`/`height` as inputs, and is cheap enough (a few
// thousand points, computed once, memoized below) to just compute directly rather than
// hardcode per-projection constants.
function computeExtent(projection: Projection, width: number, height: number): RenderTransform {
  let minX = Infinity;
  let maxX = -Infinity;
  let minY = Infinity;
  let maxY = -Infinity;
  const LAT_SAMPLES = 181;
  const LON_SAMPLES = 361;
  for (let i = 0; i < LAT_SAMPLES; i++) {
    const lat = (-90 + (180 * i) / (LAT_SAMPLES - 1)) * (Math.PI / 180);
    for (let j = 0; j < LON_SAMPLES; j++) {
      const lon = (-180 + (360 * j) / (LON_SAMPLES - 1)) * (Math.PI / 180);
      const [x, y] = project(projection, lat, lon);
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    }
  }

  const pixelScale = width / REFERENCE_WIDTH_PX;
  const paddingPx = PADDING_PX * pixelScale;
  const dataW = Math.max(maxX - minX, 1e-9);
  const dataH = Math.max(maxY - minY, 1e-9);
  const scale = Math.min((width - 2 * paddingPx) / dataW, (height - 2 * paddingPx) / dataH);
  const offsetX = width / 2 - (scale * (minX + maxX)) / 2;
  const offsetY = height / 2 + (scale * (minY + maxY)) / 2;
  return { scale, offsetX, offsetY };
}

const extentCache = new Map<string, RenderTransform>();

export function getRenderTransform(projection: Projection, width: number, height: number): RenderTransform {
  const key = `${projection}:${width}:${height}`;
  let cached = extentCache.get(key);
  if (!cached) {
    cached = computeExtent(projection, width, height);
    extentCache.set(key, cached);
  }
  return cached;
}

export function toPixels(transform: RenderTransform, x: number, y: number): [number, number] {
  return [transform.scale * x + transform.offsetX, -transform.scale * y + transform.offsetY];
}

// Unwraps `lon` to the representative within +-pi of `referenceLon` -- a direct port of the
// technique backend render_image._project_offset uses server-side. atan2's own +pi/-pi jump
// at the antimeridian means two points a tiny true angular distance apart can still land on
// opposite sides of it; a *known-nearby* reference longitude (e.g. a shape's own rotated
// center) lets every other point in that same shape be unwrapped safely before projecting,
// so a filled shape's edge never bows across the whole map (see PlateInspector.tsx).
export function wrapLongitudeNear(lon: number, referenceLon: number): number {
  const twoPi = 2 * Math.PI;
  return referenceLon + (((((lon - referenceLon + Math.PI) % twoPi) + twoPi) % twoPi) - Math.PI);
}

// --- Graticule: meridians and parallels for the live drag preview.

const GRATICULE_STEP_DEG = 30;
const MERIDIAN_SAMPLES = 37;
const PARALLEL_SAMPLES = 73;

export interface GraticuleLine {
  points: Vec3[]; // world-space unit vectors, un-rotated -- caller applies the current rotation
}

let cachedGraticule: GraticuleLine[] | null = null;

export function getGraticule(): GraticuleLine[] {
  if (cachedGraticule) return cachedGraticule;
  const lines: GraticuleLine[] = [];

  for (let lonDeg = -180; lonDeg < 180; lonDeg += GRATICULE_STEP_DEG) {
    const lon = (lonDeg * Math.PI) / 180;
    const points: Vec3[] = [];
    for (let i = 0; i < MERIDIAN_SAMPLES; i++) {
      const lat = (-90 + (180 * i) / (MERIDIAN_SAMPLES - 1)) * (Math.PI / 180);
      points.push(latLonToXyz(lat, lon));
    }
    lines.push({ points });
  }

  for (let latDeg = -60; latDeg <= 60; latDeg += GRATICULE_STEP_DEG) {
    const lat = (latDeg * Math.PI) / 180;
    const points: Vec3[] = [];
    for (let i = 0; i < PARALLEL_SAMPLES; i++) {
      const lon = (-180 + (360 * i) / (PARALLEL_SAMPLES - 1)) * (Math.PI / 180);
      points.push(latLonToXyz(lat, lon));
    }
    lines.push({ points });
  }

  cachedGraticule = lines;
  return lines;
}
