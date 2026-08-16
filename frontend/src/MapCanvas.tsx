import { useEffect, useRef } from "react";
import type { RenderResponse } from "./api";
import { elevationColor } from "./elevationColor";

export type MapView = "elevation" | "plates";

interface Props {
  data: RenderResponse | null;
  view: MapView;
  width: number;
  height: number;
}

const PADDING_PX = 20;
const POINT_SIZE_PX = 2;
const POLE_RADIUS_PX = 5;
const ARROWHEAD_LENGTH_PX = 7;

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

    // Bounding box over every coordinate the payload carries (lines, boundary, pole,
    // velocity arrow) regardless of the active view, so switching views never rescales
    // or re-centers the map.
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
    for (const plate of data.plates) {
      for (const line of plate.lines) {
        for (const [x, y] of line.points) consider(x, y);
      }
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

    if (view === "elevation") {
      for (const plate of data.plates) {
        for (const line of plate.lines) {
          for (let i = 0; i < line.points.length; i++) {
            const [px, py] = toPixel(t, line.points[i][0], line.points[i][1]);
            ctx.fillStyle = elevationColor(line.elevation[i]);
            ctx.fillRect(px - POINT_SIZE_PX / 2, py - POINT_SIZE_PX / 2, POINT_SIZE_PX, POINT_SIZE_PX);
          }
        }
      }
      return;
    }

    // "plates" view: territory (faint, categorical) + boundary outline + pole marker +
    // velocity arrow, per plate.
    for (const plate of data.plates) {
      const color = plateColor(plate.plate_id);

      ctx.fillStyle = color;
      ctx.globalAlpha = 0.35;
      for (const line of plate.lines) {
        for (const [x, y] of line.points) {
          const [px, py] = toPixel(t, x, y);
          ctx.fillRect(px - POINT_SIZE_PX / 2, py - POINT_SIZE_PX / 2, POINT_SIZE_PX, POINT_SIZE_PX);
        }
      }
      ctx.globalAlpha = 1;

      if (plate.boundary.length > 1) {
        ctx.strokeStyle = color;
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        plate.boundary.forEach(([x, y], i) => {
          const [px, py] = toPixel(t, x, y);
          if (i === 0) ctx.moveTo(px, py);
          else ctx.lineTo(px, py);
        });
        ctx.closePath();
        ctx.stroke();
      }

      if (plate.velocity_arrow) {
        const [sx, sy] = toPixel(t, ...plate.velocity_arrow.start);
        const [ex, ey] = toPixel(t, ...plate.velocity_arrow.end);
        drawArrow(ctx, sx, sy, ex, ey, "#ffffff");
      }

      if (plate.pole) {
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
  }, [data, view, width, height]);

  return <canvas ref={canvasRef} width={width} height={height} style={{ borderRadius: 8 }} />;
}
