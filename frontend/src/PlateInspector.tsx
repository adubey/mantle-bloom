import { useEffect, useRef } from "react";
import type { PlateSummary } from "./api";
import InspectorCanvas from "./InspectorCanvas";
import type { InspectorViewProps } from "./InspectorCanvas";
import type { Mat3, Vec3 } from "./rotation";
import { plateColor } from "./platePalette";
import { beginInspectorFrame, plateSelection } from "./inspectorDrawing";
import { strokeRobustLoop } from "./screenProjector";

interface Props extends InspectorViewProps {
  worldId: string;
  plates: PlateSummary[];
  selectedPlateId: number | null;
  onSelectPlate: (id: number | null) => void;
}

const ELLIPSE_FILL_ALPHA = 0.28;
const ELLIPSE_FILL_ALPHA_SELECTED = 0.6;
// Outlines are always colored by plate identity (see plateColor) -- selected just gets full
// opacity and a thicker stroke instead of the muted default, rather than switching to a
// different color, so the outline's hue keeps meaning "which plate" in both states.
const OUTLINE_ALPHA = 0.35;
const OUTLINE_ALPHA_SELECTED = 1.0;
// Points are always this same neutral white/off-white -- deliberately *never* tinted by
// plateColor, so a dot never reads as the same color as that plate's own outline: the two
// are a different visual category (individual node vs. traced boundary) as well as usually a
// different brightness (selection state). Matches backend render_image.py's
// NODE_DOT_RADIUS_PX for the dot's on-screen size.
const POINT_RGB = "255, 255, 255";
const POINT_RADIUS_PX = 1.6;
const POINT_ALPHA = 0.18; // other (non-selected) plates -- "less visible"
const SELECTED_POINT_ALPHA = 0.9; // the selected plate -- "visible"

// Renders every plate's outline + bounding ellipse + individual node points (see backend
// app/plates.py's plate_bounding_ellipse and GET /worlds/{world_id}/plates) as an interactive display --
// no server-baked PNG for this view. Every plate is always drawn as a filled, translucent
// ellipse (so overlapping plates visibly blend) plus a plate-colored outline stroke plus its
// own points (always neutral white/off-white, on top of everything else, deliberately never
// the outline's own color, so it's always clear exactly where each point sits relative to
// the traced boundary); the selected plate is drawn last, brighter, with a full-opacity
// outline and highly visible points. Reuses the same long-press-drag rotate gesture as
// MapCanvas (see rotationDrag.ts) plus click-to-select and Tab/Shift+Tab to cycle plates.
export default function PlateInspector({ worldId, plates, selectedPlateId, onSelectPlate, ...view }: Props) {
  const { width, height, projection, rotation } = view;
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const drawPlates = (previewRotation: Mat3) => {
    const frame = beginInspectorFrame(canvasRef.current, projection, previewRotation, width, height);
    if (!frame) return;
    const { ctx, projector, pixelScale, lineWidth } = frame;
    const pointRadius = POINT_RADIUS_PX * pixelScale;
    const pointSize = pointRadius * 2;

    // Individual node points -- independent dots (no polygon edge connects one to the next),
    // so there's no antimeridian-seam-crossing fill/stroke to worry about, matching how
    // backend render_image.py's "Plates (details)" view projects its own per-node dots.
    const drawPoints = (points: Vec3[], color: string) => {
      ctx.fillStyle = color;
      for (const p of points) {
        const [px, py] = projector.point(p);
        ctx.fillRect(px - pointRadius, py - pointRadius, pointSize, pointSize);
      }
    };

    const drawOne = (plate: PlateSummary, fillAlpha: number, selected: boolean) => {
      const [r, g, b] = plateColor(plate.plate_id);
      const ellipse = plate.bounding_ellipse;
      const centerXyz: Vec3 = ellipse ? ellipse.center_xyz : plate.outline[0];
      if (!centerXyz) return;

      if (ellipse) {
        const pixels = projector.loop(ellipse.outline, centerXyz);
        ctx.fillStyle = `rgba(${r}, ${g}, ${b}, ${fillAlpha})`;
        ctx.beginPath();
        pixels.forEach(([px, py], i) => (i === 0 ? ctx.moveTo(px, py) : ctx.lineTo(px, py)));
        ctx.closePath();
        ctx.fill();
      }

      if (plate.outline.length > 0) {
        const pixels = projector.loop(plate.outline, centerXyz);
        const outlineAlpha = selected ? OUTLINE_ALPHA_SELECTED : OUTLINE_ALPHA;
        strokeRobustLoop(ctx, pixels, `rgba(${r}, ${g}, ${b}, ${outlineAlpha})`, selected ? lineWidth * 2 : lineWidth, width);
      }

      // Drawn last (on top of both the ellipse fill and the outline stroke above), in a
      // color the outline never uses, so it's always possible to tell exactly where each
      // point sits relative to the traced boundary.
      const pointAlpha = selected ? SELECTED_POINT_ALPHA : POINT_ALPHA;
      drawPoints(plate.points, `rgba(${POINT_RGB}, ${pointAlpha})`);
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

  return <InspectorCanvas {...view} canvasRef={canvasRef} onFrame={drawPlates} {...plateSelection(worldId, plates, selectedPlateId, onSelectPlate)} />;
}
