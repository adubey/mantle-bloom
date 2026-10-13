import type { PlateSummary, Projection } from "./api";
import { fetchPlateAt } from "./api";
import type { Mat3, Vec3 } from "./rotation";
import { ScreenProjector, cycleId } from "./screenProjector";

export const INSPECTOR_BACKGROUND = "#0b1020";
// Matches backend render_image.py's COASTLINE_COLOR_RGB/COASTLINE_HALO_RGB: a dark halo under a
// light line, so the coastline reads over any backdrop -- the inspectors have no filled
// backdrop, so without it there'd be no land/ocean cue at all.
export const COASTLINE_RGB = "235, 235, 235";
export const COASTLINE_HALO_RGB = "15, 15, 15";

export interface InspectorFrame {
  ctx: CanvasRenderingContext2D;
  projector: ScreenProjector;
  pixelScale: number; // relative to the 1100 px reference width
  lineWidth: number;
}

// Clears `canvas` to the inspector background and sets up projection for one frame at
// `previewRotation`. Null until the canvas is mounted.
export function beginInspectorFrame(
  canvas: HTMLCanvasElement | null, projection: Projection, previewRotation: Mat3, width: number, height: number,
): InspectorFrame | null {
  const ctx = canvas?.getContext("2d");
  if (!canvas || !ctx) return null;
  ctx.fillStyle = INSPECTOR_BACKGROUND;
  ctx.fillRect(0, 0, width, height);
  const pixelScale = width / 1100;
  return {
    ctx, projector: new ScreenProjector(projection, previewRotation, width, height), pixelScale, lineWidth: Math.max(1, pixelScale),
  };
}

// Strokes each edge on its own (a flat edge list, not a polyline -- a river network branches),
// projected with ScreenProjector.segment so no edge bows across the antimeridian seam.
export function strokeSegments(
  ctx: CanvasRenderingContext2D, projector: ScreenProjector, segments: [Vec3, Vec3][], color: string, width: number,
  dash: number[] = [],
) {
  ctx.strokeStyle = color;
  ctx.lineWidth = width;
  ctx.setLineDash(dash);
  for (const [a, b] of segments) {
    const [[x1, y1], [x2, y2]] = projector.segment(a, b);
    ctx.beginPath();
    ctx.moveTo(x1, y1);
    ctx.lineTo(x2, y2);
    ctx.stroke();
  }
  ctx.setLineDash([]);
}

// The halo, then the lighter line on top.
export function strokeCoastline(
  ctx: CanvasRenderingContext2D, projector: ScreenProjector, segments: [Vec3, Vec3][], lineWidth: number,
) {
  strokeSegments(ctx, projector, segments, `rgba(${COASTLINE_HALO_RGB}, 1.0)`, lineWidth * 2.6);
  strokeSegments(ctx, projector, segments, `rgba(${COASTLINE_RGB}, 1.0)`, lineWidth * 1.1);
}

// Click-to-select and Tab cycling for the two plate views (PlateInspector, PlatesAndFaults),
// as InspectorCanvas's onPick/onCycle.
export function plateSelection(
  plates: PlateSummary[], selectedPlateId: number | null, onSelectPlate: (id: number | null) => void,
) {
  return {
    onPick: (latDeg: number, lonDeg: number) =>
      fetchPlateAt(latDeg, lonDeg).then(({ plate_id }) => onSelectPlate(plate_id)),
    onCycle: (backwards: boolean) => {
      const next = cycleId(plates.map((p) => p.plate_id), selectedPlateId, backwards);
      if (next != null) onSelectPlate(next);
    },
  };
}
