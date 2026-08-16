import { useEffect, useRef } from "react";
import type { RenderResponse } from "./api";
import { elevationColor } from "./elevationColor";

interface Props {
  data: RenderResponse | null;
  width: number;
  height: number;
}

const PADDING_PX = 20;
const POINT_SIZE_PX = 2;

export default function MapCanvas({ data, width, height }: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    ctx.fillStyle = "#0b1020";
    ctx.fillRect(0, 0, width, height);

    if (!data || data.plates.length === 0) return;

    let minX = Infinity;
    let maxX = -Infinity;
    let minY = Infinity;
    let maxY = -Infinity;
    for (const plate of data.plates) {
      for (const line of plate.lines) {
        for (const [x, y] of line.points) {
          if (x < minX) minX = x;
          if (x > maxX) maxX = x;
          if (y < minY) minY = y;
          if (y > maxY) maxY = y;
        }
      }
    }
    if (!isFinite(minX)) return;

    const dataW = Math.max(maxX - minX, 1e-9);
    const dataH = Math.max(maxY - minY, 1e-9);
    const scale = Math.min((width - 2 * PADDING_PX) / dataW, (height - 2 * PADDING_PX) / dataH);
    const offsetX = width / 2 - (scale * (minX + maxX)) / 2;
    const offsetY = height / 2 + (scale * (minY + maxY)) / 2;

    for (const plate of data.plates) {
      for (const line of plate.lines) {
        for (let i = 0; i < line.points.length; i++) {
          const [x, y] = line.points[i];
          const px = scale * x + offsetX;
          const py = -scale * y + offsetY;
          ctx.fillStyle = elevationColor(line.elevation[i]);
          ctx.fillRect(px - POINT_SIZE_PX / 2, py - POINT_SIZE_PX / 2, POINT_SIZE_PX, POINT_SIZE_PX);
        }
      }
    }
  }, [data, width, height]);

  return <canvas ref={canvasRef} width={width} height={height} style={{ borderRadius: 8 }} />;
}
