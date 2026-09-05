import { useEffect, useRef } from "react";
import { plateColor } from "./platePalette";

// A rough, client-side-only preview of how the plate tiling would carve up a Human-made
// sketch, shown live under the "Voronoi points" slider in AdvancedSettingsModal. It is a
// deliberate approximation of backend lithosphere_plate.build_plate_tiling / worldsketch.
// sketch_plate_sites -- enough to see how boundary density and land/ocean bias change with the
// slider, not a faithful port. The real generation runs on the sphere with angular
// region-growing; here we work on a small equirectangular grid with planar (longitude-wrapped)
// distance.

// Working grid (2:1 equirectangular, matching the sketch canvas aspect). Small on purpose --
// the per-cell nearest-point scan below is O(cells x points).
const GW = 240;
const GH = 120;
// On-screen size the small grid is scaled up to (smoothing disabled, so cells stay crisp).
const VIEW_W = 336;
const VIEW_H = 168;
// Debounce so dragging the slider doesn't recompute on every intermediate value.
const RECOMPUTE_DEBOUNCE_MS = 120;
// Luminance below this (0-255) counts as coastline/ink -- a flood-fill barrier and drawn dark
// on top of the plate colours. White paper is 255; the sketch's coast/river/mountain inks are
// all well below this (see SketchEditor.tsx's TOOL_COLORS).
const INK_LUMA_THRESHOLD = 150;

interface Props {
  // Full data URL (`data:image/png;base64,...`), same shape App.tsx stores the sketch as.
  sketchImageDataUrl: string;
  seed: number;
  // The "Voronoi points" slider value -- total seed points.
  numPoints: number;
  // Plate count to merge down to (App passes the explicit count, or its DEFAULT_PLATES
  // stand-in when "Auto" is on).
  plateCount: number;
  continentalPercent: number;
}

// Small deterministic PRNG (mulberry32) -- keeps the preview stable for a given seed the way
// the real seed-driven generation is.
function mulberry32(a: number): () => number {
  return () => {
    a |= 0;
    a = (a + 0x6d2b79f5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function wrappedDx(ax: number, bx: number): number {
  let dx = Math.abs(ax - bx);
  if (dx > GW / 2) dx = GW - dx;
  return dx;
}

interface Cell {
  ocean: Uint8Array;
  coast: Uint8Array;
}

// Classify each grid cell as ink/ocean/land: ocean is flood-filled from the four map corners
// through non-ink cells (x wraps, y does not), matching worldsketch.py's corner-seeded fill.
function classifyGrid(img: ImageData): Cell {
  const n = GW * GH;
  const coast = new Uint8Array(n);
  for (let i = 0; i < n; i++) {
    const r = img.data[i * 4];
    const g = img.data[i * 4 + 1];
    const b = img.data[i * 4 + 2];
    const luma = 0.299 * r + 0.587 * g + 0.114 * b;
    coast[i] = luma < INK_LUMA_THRESHOLD ? 1 : 0;
  }
  const ocean = new Uint8Array(n);
  const stack: number[] = [];
  const seedCell = (x: number, y: number) => {
    const i = y * GW + x;
    if (!coast[i] && !ocean[i]) {
      ocean[i] = 1;
      stack.push(i);
    }
  };
  seedCell(0, 0);
  seedCell(GW - 1, 0);
  seedCell(0, GH - 1);
  seedCell(GW - 1, GH - 1);
  while (stack.length) {
    const i = stack.pop()!;
    const x = i % GW;
    const y = (i / GW) | 0;
    const neighbors = [
      [(x + 1) % GW, y],
      [(x - 1 + GW) % GW, y],
      [x, y - 1],
      [x, y + 1],
    ];
    for (const [nx, ny] of neighbors) {
      if (ny < 0 || ny >= GH) continue;
      const ni = ny * GW + nx;
      if (!coast[ni] && !ocean[ni]) {
        ocean[ni] = 1;
        stack.push(ni);
      }
    }
  }
  return { ocean, coast };
}

interface Point {
  x: number;
  y: number;
  plate: number;
}

function buildPoints(cells: Cell, seed: number, numPoints: number, plateCount: number, continentalPercent: number): Point[] {
  const rand = mulberry32((seed >>> 0) || 1);
  const land: number[] = [];
  const ocean: number[] = [];
  for (let i = 0; i < GW * GH; i++) {
    if (cells.coast[i]) continue;
    if (cells.ocean[i]) ocean.push(i);
    else land.push(i);
  }
  const pick = (pool: number[], fallback: number[]): number => {
    const src = pool.length ? pool : fallback.length ? fallback : null;
    if (!src) return ((rand() * GH) | 0) * GW + ((rand() * GW) | 0);
    return src[(rand() * src.length) | 0];
  };

  const plates = Math.max(1, plateCount);
  const nCont = Math.min(plates, Math.max(0, Math.round((continentalPercent / 100) * plates)));
  const primaries: Point[] = [];
  for (let p = 0; p < plates; p++) {
    const cell = p < nCont ? pick(land, ocean) : pick(ocean, land);
    primaries.push({ x: cell % GW, y: (cell / GW) | 0, plate: p });
  }

  const points: Point[] = [...primaries];
  const extras = Math.max(0, numPoints - plates);
  for (let e = 0; e < extras; e++) {
    const x = rand() * GW;
    const y = rand() * GH;
    // Merge each extra point into the nearest primary's plate (approximates the real
    // angular region-growing).
    let best = 0;
    let bestD = Infinity;
    for (let p = 0; p < primaries.length; p++) {
      const dx = wrappedDx(x, primaries[p].x);
      const dy = y - primaries[p].y;
      const d = dx * dx + dy * dy;
      if (d < bestD) {
        bestD = d;
        best = primaries[p].plate;
      }
    }
    points.push({ x, y, plate: best });
  }
  return points;
}

function render(canvas: HTMLCanvasElement, points: Point[], cells: Cell) {
  const out = new ImageData(GW, GH);
  for (let y = 0; y < GH; y++) {
    for (let x = 0; x < GW; x++) {
      const i = y * GW + x;
      let best = 0;
      let bestD = Infinity;
      for (let p = 0; p < points.length; p++) {
        const dx = wrappedDx(x + 0.5, points[p].x);
        const dy = y + 0.5 - points[p].y;
        const d = dx * dx + dy * dy;
        if (d < bestD) {
          bestD = d;
          best = points[p].plate;
        }
      }
      const [r, g, b] = plateColor(best);
      const o = i * 4;
      if (cells.coast[i]) {
        // Drawn coastline/ink stays visible on top of the plate fill.
        out.data[o] = 26;
        out.data[o + 1] = 26;
        out.data[o + 2] = 26;
      } else {
        out.data[o] = r;
        out.data[o + 1] = g;
        out.data[o + 2] = b;
      }
      out.data[o + 3] = 255;
    }
  }
  const src = document.createElement("canvas");
  src.width = GW;
  src.height = GH;
  src.getContext("2d")!.putImageData(out, 0, 0);
  const ctx = canvas.getContext("2d")!;
  ctx.imageSmoothingEnabled = false;
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.drawImage(src, 0, 0, canvas.width, canvas.height);
}

export default function VoronoiPreview({ sketchImageDataUrl, seed, numPoints, plateCount, continentalPercent }: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    let cancelled = false;
    const timer = window.setTimeout(() => {
      const img = new Image();
      img.onload = () => {
        if (cancelled || !canvasRef.current) return;
        const grid = document.createElement("canvas");
        grid.width = GW;
        grid.height = GH;
        const gctx = grid.getContext("2d");
        if (!gctx) return;
        gctx.drawImage(img, 0, 0, GW, GH);
        const data = gctx.getImageData(0, 0, GW, GH);
        const cells = classifyGrid(data);
        const points = buildPoints(cells, seed, numPoints, plateCount, continentalPercent);
        render(canvasRef.current, points, cells);
      };
      img.src = sketchImageDataUrl;
    }, RECOMPUTE_DEBOUNCE_MS);
    return () => {
      cancelled = true;
      window.clearTimeout(timer);
    };
  }, [sketchImageDataUrl, seed, numPoints, plateCount, continentalPercent]);

  return (
    <canvas
      ref={canvasRef}
      width={VIEW_W}
      height={VIEW_H}
      style={{
        width: "100%",
        maxWidth: VIEW_W,
        aspectRatio: "2 / 1",
        borderRadius: 4,
        border: "1px solid #333",
        display: "block",
      }}
    />
  );
}
