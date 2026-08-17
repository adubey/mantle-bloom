import { useEffect, useRef } from "react";
import type { KeyboardEvent } from "react";
import type { PlateSummary, Projection } from "./api";
import { fetchPlateAt } from "./api";
import type { Mat3, RenderTransform, Vec3 } from "./rotation";
import {
  getRenderTransform, latLonToXyz, matApply, matTranspose, project, toPixels, unproject, wrapLongitudeNear, xyzToLatLon,
} from "./rotation";
import { useRotationDrag } from "./rotationDrag";
import { plateColor } from "./platePalette";

interface Props {
  plates: PlateSummary[];
  width: number;
  height: number;
  displayWidth: number;
  displayHeight: number;
  projection: Projection;
  // The currently *committed* view rotation -- same rotation state MapCanvas shares, so
  // switching between the two views preserves orientation.
  rotation: Mat3;
  selectedPlateId: number | null;
  onSelectPlate: (id: number | null) => void;
  onRotationPreview: (latDeg: number, lonDeg: number) => void;
  onRotationCommitted: (rotation: Mat3) => void;
}

const BACKGROUND = "#0b1020";
const ELLIPSE_FILL_ALPHA = 0.28;
const ELLIPSE_FILL_ALPHA_SELECTED = 0.6;
const OUTLINE_ALPHA = 0.35;
const SELECTED_OUTLINE_RGB = "255, 255, 255";
// Matches backend render_image.py's NODE_DOT_RADIUS_PX -- same visual size as the "Plates
// (details)" view's own per-node dots.
const POINT_RADIUS_PX = 1.6;
const POINT_ALPHA = 0.18; // other (non-selected) plates -- "less visible"
const SELECTED_POINT_RGB = "255, 255, 255"; // the selected plate -- "visible"
const SELECTED_POINT_ALPHA = 0.9;
// Same "skip segments longer than N times the median" seam-break technique MapCanvas's
// graticule uses (itself a port of backend render_image._stroke_robust_loop) -- necessary
// for *stroking* a shape that crosses the antimeridian even after per-shape unwrapping (see
// projectLoop below), since a shape's true angular extent can still, rarely, exceed what a
// single reference-longitude unwrap can fully resolve.
const SEGMENT_BREAK_FACTOR = 6;

// Renders every plate's outline + bounding ellipse + individual node points (see backend
// app/plates.py's plate_bounding_ellipse and GET /world/plates) as an interactive display --
// no server-baked PNG for this view. Every plate is always drawn as a filled, translucent
// ellipse (so overlapping plates visibly blend) plus its own points, dim; the selected plate
// is drawn last, brighter, with a crisp white true-boundary outline and its points in a
// highly visible white. Reuses the same long-press-drag rotate gesture as MapCanvas (see
// rotationDrag.ts) plus click-to-select and Tab/Shift+Tab to cycle plates.
export default function PlateInspector({
  plates, width, height, displayWidth, displayHeight, projection, rotation,
  selectedPlateId, onSelectPlate, onRotationPreview, onRotationCommitted,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);

  // Projects a closed loop of *true* world points, unwrapping every point's longitude
  // relative to the loop's own (rotated) center first -- see rotation.ts's
  // wrapLongitudeNear's docstring for why a naive projection of a rotated loop can bow a
  // filled shape across the whole map at the antimeridian seam.
  const projectLoop = (worldPts: Vec3[], centerXyz: Vec3, previewRotation: Mat3, transform: RenderTransform): [number, number][] => {
    const [, centerLon] = xyzToLatLon(matApply(previewRotation, centerXyz));
    return worldPts.map((p) => {
      const r = matApply(previewRotation, p);
      const lat = Math.asin(Math.min(1, Math.max(-1, r[2])));
      const lon = wrapLongitudeNear(Math.atan2(r[1], r[0]), centerLon);
      const [x, y] = project(projection, lat, lon);
      return toPixels(transform, x, y);
    });
  };

  const strokeRobustLoop = (
    ctx: CanvasRenderingContext2D, pixels: [number, number][], color: string, lineWidth: number,
  ) => {
    const n = pixels.length;
    if (n < 2) return;
    const segLengths: number[] = [];
    for (let i = 0; i < n; i++) {
      const next = pixels[(i + 1) % n];
      segLengths.push(Math.hypot(next[0] - pixels[i][0], next[1] - pixels[i][1]));
    }
    const sorted = [...segLengths].sort((a, b) => a - b);
    const median = sorted[Math.floor(sorted.length / 2)] || 0;
    const breakThreshold = Math.max(median * SEGMENT_BREAK_FACTOR, 20 * (width / 1100));
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
  };

  const drawPlates = (previewRotation: Mat3) => {
    const canvas = canvasRef.current;
    const ctx = canvas?.getContext("2d");
    if (!canvas || !ctx) return;
    ctx.fillStyle = BACKGROUND;
    ctx.fillRect(0, 0, width, height);

    const transform = getRenderTransform(projection, width, height);
    const pixelScale = width / 1100;
    const lineWidth = Math.max(1, pixelScale);
    const pointRadius = POINT_RADIUS_PX * pixelScale;
    const pointSize = pointRadius * 2;

    // Individual node points -- unlike projectLoop, these are independent dots (no polygon
    // edge connects one to the next), so there's no antimeridian-seam-crossing fill/stroke
    // to worry about and no wrapLongitudeNear step is needed, matching how backend
    // render_image.py's "Plates (details)" view projects its own per-node dots.
    const drawPoints = (points: Vec3[], color: string) => {
      ctx.fillStyle = color;
      for (const p of points) {
        const r = matApply(previewRotation, p);
        const lat = Math.asin(Math.min(1, Math.max(-1, r[2])));
        const lon = Math.atan2(r[1], r[0]);
        const [x, y] = project(projection, lat, lon);
        const [px, py] = toPixels(transform, x, y);
        ctx.fillRect(px - pointRadius, py - pointRadius, pointSize, pointSize);
      }
    };

    const drawOne = (plate: PlateSummary, fillAlpha: number, selected: boolean) => {
      const [r, g, b] = plateColor(plate.plate_id);
      const ellipse = plate.bounding_ellipse;
      const centerXyz: Vec3 = ellipse ? ellipse.center_xyz : plate.outline[0];
      if (!centerXyz) return;

      if (ellipse) {
        const pixels = projectLoop(ellipse.outline, centerXyz, previewRotation, transform);
        ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${fillAlpha})`;
        ctx.beginPath();
        pixels.forEach(([px, py], i) => (i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py)));
        ctx.closePath();
        ctx.fill();
      }

      if (plate.outline.length > 0) {
        const pixels = projectLoop(plate.outline, centerXyz, previewRotation, transform);
        const color = selected ? `rgba(${SELECTED_OUTLINE_RGB}, 1)` : `rgba(${r}, ${g}, ${b}, ${OUTLINE_ALPHA})`;
        strokeRobustLoop(ctx, pixels, color, selected ? lineWidth * 2 : lineWidth);
      }

      const pointColor = selected
        ? `rgba(${SELECTED_POINT_RGB}, ${SELECTED_POINT_ALPHA})`
        : `rgba(${r}, ${g}, ${b}, ${POINT_ALPHA})`;
      drawPoints(plate.points, pointColor);
    };

    for (const plate of plates) {
      if (plate.plate_id === selectedPlateId) continue;
      drawOne(plate, ELLIPSE_FILL_ALPHA, false);
    }
    const selected = plates.find((p) => p.plate_id === selectedPlateId);
    if (selected) drawOne(selected, ELLIPSE_FILL_ALPHA_SELECTED, true);
  };

  useEffect(() => {
    drawPlates(rotation);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [plates, selectedPlateId, projection, rotation, width, height]);

  const handleClick = (backingX: number, backingY: number) => {
    const transform = getRenderTransform(projection, width, height);
    const x = (backingX - transform.offsetX) / transform.scale;
    const y = -(backingY - transform.offsetY) / transform.scale;
    const latLon = unproject(projection, x, y);
    if (!latLon) return; // clicked outside the sphere's projected silhouette
    const displayXyz = latLonToXyz(latLon[0], latLon[1]); // in the *display* (rotated) frame
    const trueXyz = matApply(matTranspose(rotation), displayXyz); // transpose == inverse for a rotation matrix
    const [trueLat, trueLon] = xyzToLatLon(trueXyz);
    fetchPlateAt((trueLat * 180) / Math.PI, (trueLon * 180) / Math.PI)
      .then(({ plate_id }) => {
        onSelectPlate(plate_id);
        containerRef.current?.focus();
      })
      .catch(() => {
        // A click that races a generate/step (no world, or a stale plate list) -- ignored,
        // matching how a stray click elsewhere in the app fails silently rather than
        // surfacing a toast for a purely exploratory gesture.
      });
  };

  useRotationDrag({
    elementRef: canvasRef,
    width, height, displayWidth, displayHeight, rotation,
    onFrame: drawPlates,
    onRotationPreview,
    onRotationCommitted,
    onClick: handleClick,
  });

  useEffect(() => {
    containerRef.current?.focus();
  }, []);

  const handleKeyDown = (e: KeyboardEvent<HTMLDivElement>) => {
    if (e.key !== "Tab") return;
    e.preventDefault();
    if (plates.length === 0) return;
    const ids = plates.map((p) => p.plate_id).sort((a, b) => a - b);
    const currentIndex = selectedPlateId == null ? -1 : ids.indexOf(selectedPlateId);
    const delta = e.shiftKey ? -1 : 1;
    const nextIndex = currentIndex === -1 ? (e.shiftKey ? ids.length - 1 : 0) : (currentIndex + delta + ids.length) % ids.length;
    onSelectPlate(ids[nextIndex]);
  };

  return (
    <div ref={containerRef} tabIndex={0} onKeyDown={handleKeyDown} style={{ outline: "none", display: "inline-block" }}>
      <canvas
        ref={canvasRef}
        width={width}
        height={height}
        style={{ borderRadius: 8, width: displayWidth, height: displayHeight }}
      />
    </div>
  );
}
