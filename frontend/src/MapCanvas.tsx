import { useEffect, useRef } from "react";
import type { RenderResponse } from "./api";
import { elevationColor } from "./elevationColor";

export type MapView = "elevation" | "plates" | "platesDetail";

interface Props {
  data: RenderResponse | null;
  view: MapView;
  width: number;
  height: number;
}

const PADDING_PX = 20;
const POLE_RADIUS_PX = 5;
const ARROWHEAD_LENGTH_PX = 7;
const NODE_DOT_RADIUS_PX = 1.6;
// Cells are drawn slightly larger than their measured half-extent so adjacent cells
// overlap a hair rather than risk a hairline gap from floating-point rounding.
const CELL_OVERLAP_FACTOR = 1.15;
// Velocity arrows are always meant to be short, local indicators (backend
// ARROW_BASE_ANGULAR_LENGTH_RAD caps them at ~0.15 rad on the sphere). If a plate's seed
// happens to sit near the map's antimeridian or a pole, projecting its (geodesically short)
// arrow endpoint can land it on the "other side" of the projection -- a tiny step in world
// space, but a huge jump in projected coordinates -- and drawing a straight line between
// the two would paint a stray line across most of the map. Any arrow this long on screen is
// that artifact, not a real velocity indicator, so it's skipped rather than drawn.
const MAX_ARROW_FRACTION_OF_CANVAS = 0.15;

// A fixed categorical palette so each plate reads as a distinct region across
// generate/step calls (plate_id is stable within one world's lifetime).
const PLATE_PALETTE = [
  "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
  "#42d4f4", "#f032e6", "#bcf60c", "#fabebe", "#469990",
  "#e6beff", "#9a6324", "#fffac8", "#800000", "#aaffc3",
  "#808000", "#ffd8b1", "#000075", "#a9a9a9", "#ffffff",
];

function plateColor(plateId: number): string {
  return PLATE_PALETTE[plateId % PLATE_PALETTE.length];
}

interface Transform {
  scale: number;
  offsetX: number;
  offsetY: number;
}

function toPixel(t: Transform, x: number, y: number): [number, number] {
  return [t.scale * x + t.offsetX, -t.scale * y + t.offsetY];
}

// A plate's live outline (see backend Plate.outline_world) is built from many short
// consecutive hops, but a plate whose local frame happens to put its own pole -- or the
// map projection's antimeridian -- inside its territory can still produce one or two
// segments that are wildly longer than the rest once projected. Rather than draw those as
// a straight line straight across the map, break the path there: find the typical
// (median) segment length for this loop and skip drawing anything far longer than that.
function strokeRobustLoop(ctx: CanvasRenderingContext2D, points: [number, number][], t: Transform, color: string) {
  if (points.length < 2) return;
  const pixels = points.map(([x, y]) => toPixel(t, x, y));
  const n = pixels.length;
  const lengths = pixels.map(([x0, y0], i) => {
    const [x1, y1] = pixels[(i + 1) % n];
    return Math.hypot(x1 - x0, y1 - y0);
  });
  const sorted = [...lengths].sort((a, b) => a - b);
  const median = sorted[Math.floor(sorted.length / 2)] || 0;
  const breakThreshold = Math.max(median * 6, 20);

  ctx.strokeStyle = color;
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  let penDown = false;
  for (let i = 0; i < n; i++) {
    const [x0, y0] = pixels[i];
    const [x1, y1] = pixels[(i + 1) % n];
    if (!penDown) {
      ctx.moveTo(x0, y0);
      penDown = true;
    }
    if (lengths[i] > breakThreshold) {
      penDown = false;
    } else {
      ctx.lineTo(x1, y1);
    }
  }
  ctx.stroke();
}

function drawArrow(ctx: CanvasRenderingContext2D, x0: number, y0: number, x1: number, y1: number, color: string) {
  ctx.strokeStyle = color;
  ctx.fillStyle = color;
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.moveTo(x0, y0);
  ctx.lineTo(x1, y1);
  ctx.stroke();

  const angle = Math.atan2(y1 - y0, x1 - x0);
  ctx.beginPath();
  ctx.moveTo(x1, y1);
  ctx.lineTo(
    x1 - ARROWHEAD_LENGTH_PX * Math.cos(angle - Math.PI / 6),
    y1 - ARROWHEAD_LENGTH_PX * Math.sin(angle - Math.PI / 6),
  );
  ctx.lineTo(
    x1 - ARROWHEAD_LENGTH_PX * Math.cos(angle + Math.PI / 6),
    y1 - ARROWHEAD_LENGTH_PX * Math.sin(angle + Math.PI / 6),
  );
  ctx.closePath();
  ctx.fill();
}

export default function MapCanvas({ data, view, width, height }: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    ctx.fillStyle = "#0b1020";
    ctx.fillRect(0, 0, width, height);

    if (!data || data.plates.length === 0) return;

    // Bounding box over every coordinate the payload carries (the render grid, plate
    // boundaries, poles, velocity arrows) regardless of the active view, so switching
    // views never rescales or re-centers the map.
    let minX = Infinity;
    let maxX = -Infinity;
    let minY = Infinity;
    let maxY = -Infinity;
    const consider = (x: number, y: number) => {
      if (x < minX) minX = x;
      if (x > maxX) maxX = x;
      if (y < minY) minY = y;
      if (y > maxY) maxY = y;
    };
    for (const [x, y] of data.grid.points) consider(x, y);
    for (const plate of data.plates) {
      for (const [x, y] of plate.boundary) consider(x, y);
      if (plate.pole) consider(plate.pole[0], plate.pole[1]);
      if (plate.velocity_arrow) {
        consider(...plate.velocity_arrow.start);
        consider(...plate.velocity_arrow.end);
      }
    }
    if (!isFinite(minX)) return;

    const dataW = Math.max(maxX - minX, 1e-9);
    const dataH = Math.max(maxY - minY, 1e-9);
    const scale = Math.min((width - 2 * PADDING_PX) / dataW, (height - 2 * PADDING_PX) / dataH);
    const t: Transform = {
      scale,
      offsetX: width / 2 - (scale * (minX + maxX)) / 2,
      offsetY: height / 2 + (scale * (minY + maxY)) / 2,
    };

    // The render grid: a uniform lat/lon sweep with every cell already assigned its
    // nearest data point (see backend main.py's _render_grid), each drawn at the size its
    // own row actually needs (cell_half_width/height) so the filled map has no gaps
    // regardless of projection distortion. Skipped for "platesDetail", which shows the raw
    // per-plate line data (with its real gaps) instead of this gap-free reprojection of it.
    if (view === "elevation" || view === "plates") {
      const { grid } = data;
      for (let i = 0; i < grid.points.length; i++) {
        const [x, y] = grid.points[i];
        const [px, py] = toPixel(t, x, y);
        const w = grid.cell_half_width[i] * t.scale * 2 * CELL_OVERLAP_FACTOR;
        const h = grid.cell_half_height[i] * t.scale * 2 * CELL_OVERLAP_FACTOR;

        if (view === "elevation") {
          ctx.fillStyle = elevationColor(grid.elevation[i]);
        } else {
          ctx.fillStyle = plateColor(grid.plate_id[i]);
        }
        ctx.fillRect(px - w / 2, py - h / 2, w, h);
      }
    }

    if (view === "platesDetail") {
      // The actual line-based data structure: each plate's elevation nodes, drawn as dots
      // colored by elevation rather than resampled onto a gap-free grid -- see this dot's
      // own point density (and the gaps between lines) for what the simulation really
      // tracks per plate.
      for (const plate of data.plates) {
        for (const line of plate.lines) {
          for (let i = 0; i < line.points.length; i++) {
            const [x, y] = line.points[i];
            const [px, py] = toPixel(t, x, y);
            ctx.fillStyle = elevationColor(line.elevation[i]);
            ctx.beginPath();
            ctx.arc(px, py, NODE_DOT_RADIUS_PX, 0, 2 * Math.PI);
            ctx.fill();
          }
        }
      }
    }

    if (view === "plates" || view === "platesDetail") {
      for (const plate of data.plates) {
        const color = plateColor(plate.plate_id);
        strokeRobustLoop(ctx, plate.boundary, t, color);

        if (view === "plates" && plate.velocity_arrow) {
          const [sx, sy] = toPixel(t, ...plate.velocity_arrow.start);
          const [ex, ey] = toPixel(t, ...plate.velocity_arrow.end);
          const maxArrowPx = Math.min(width, height) * MAX_ARROW_FRACTION_OF_CANVAS;
          if (Math.hypot(ex - sx, ey - sy) <= maxArrowPx) {
            drawArrow(ctx, sx, sy, ex, ey, "#ffffff");
          }
        }

        if (view === "plates" && plate.pole) {
          const [px, py] = toPixel(t, plate.pole[0], plate.pole[1]);
          ctx.fillStyle = "#ff2d55";
          ctx.strokeStyle = "#ffffff";
          ctx.lineWidth = 1;
          ctx.beginPath();
          ctx.arc(px, py, POLE_RADIUS_PX, 0, 2 * Math.PI);
          ctx.fill();
          ctx.stroke();
        }
      }
    }
  }, [data, view, width, height]);

  return <canvas ref={canvasRef} width={width} height={height} style={{ borderRadius: 8 }} />;
}
