import { useEffect, useRef } from "react";
import type { RiverSummary, Segment } from "./api";
import { fetchRiverAt } from "./api";
import InspectorCanvas from "./InspectorCanvas";
import type { InspectorViewProps } from "./InspectorCanvas";
import type { Mat3 } from "./rotation";
import { isVisibleRiver } from "./riverVisibility";
import { beginInspectorFrame, strokeCoastline, strokeSegments } from "./inspectorDrawing";
import { cycleId } from "./screenProjector";

interface Props extends InspectorViewProps {
  rivers: RiverSummary[];
  coastlineSegments: Segment[];
  selectedRiverId: number | null;
  onSelectRiver: (id: number | null) => void;
}

// Matches backend render_image.py's RIVER_COLOR_RGB, so a river reads the same color here as
// it does baked into the elevation view's rendering.
const RIVER_RGB = "77, 216, 230";
const RIVER_ALPHA = 0.3; // other (non-selected) rivers -- "less visible"
const SELECTED_RIVER_ALPHA = 1.0;

// The endpoint circle's color carries the same mouth-type classification the sidebar text
// does -- a quick visual read of "where does this river end" without needing to read the
// panel, without inventing a whole new legend (these three colors are close enough to the
// existing lake/ocean/land palette elsewhere in the app to read intuitively).
const MOUTH_RING_COLOR: Record<RiverSummary["mouth_type"], string> = {
  ocean: "rgba(90, 170, 255, 1.0)",
  lake: "rgba(120, 220, 190, 1.0)",
  other: "rgba(230, 170, 70, 1.0)",
};
const MOUTH_RING_RADIUS_PX = 7;


// Renders every river network's flow-edge segments (see backend app/hydrology.py's
// group_rivers and GET /world/rivers) as an interactive display -- no server-baked PNG for
// this view, same "raw JSON, client draws it" philosophy as PlateInspector.tsx. Every river
// is always drawn, dim; the selected one is drawn last, bright, plus a ring at its own mouth
// (colored by endpoint type). Reuses the same long-press-drag rotate gesture as MapCanvas/
// PlateInspector (see rotationDrag.ts) plus click-to-select and Tab/Shift+Tab to cycle rivers.
export default function RiverInspector({ rivers, coastlineSegments, selectedRiverId, onSelectRiver, ...view }: Props) {
  const { width, height, projection, rotation } = view;
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const drawRivers = (previewRotation: Mat3) => {
    const frame = beginInspectorFrame(canvasRef.current, projection, previewRotation, width, height);
    if (!frame) return;
    const { ctx, projector, pixelScale, lineWidth } = frame;

    strokeCoastline(ctx, projector, coastlineSegments, lineWidth);

    const drawOne = (river: RiverSummary, selected: boolean) => {
      const alpha = selected ? SELECTED_RIVER_ALPHA : RIVER_ALPHA;
      strokeSegments(ctx, projector, river.segments, `rgba(${RIVER_RGB}, ${alpha})`, selected ? lineWidth * 2.5 : lineWidth);

      if (!selected) return;
      const [px, py] = projector.point(river.mouth_xyz);
      const ringRadius = MOUTH_RING_RADIUS_PX * pixelScale;
      ctx.strokeStyle = MOUTH_RING_COLOR[river.mouth_type];
      ctx.lineWidth = lineWidth * 2;
      ctx.beginPath();
      ctx.arc(px, py, ringRadius, 0, 2 * Math.PI);
      ctx.stroke();
      ctx.fillStyle = MOUTH_RING_COLOR[river.mouth_type];
      ctx.beginPath();
      ctx.arc(px, py, ringRadius * 0.3, 0, 2 * Math.PI);
      ctx.fill();
    };

    const visibleRivers = rivers.filter(isVisibleRiver);
    for (const river of visibleRivers) {
      if (river.river_id === selectedRiverId) continue;
      drawOne(river, false);
    }
    const selected = visibleRivers.find((r) => r.river_id === selectedRiverId);
    if (selected) drawOne(selected, true);
  };

  useEffect(() => {
    drawRivers(rotation);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [rivers, coastlineSegments, selectedRiverId, projection, rotation, width, height]);

  // A click on a river too minor to draw selects nothing, rather than an invisible river.
  const pickRiver = (latDeg: number, lonDeg: number) =>
    fetchRiverAt(latDeg, lonDeg).then(({ river_id }) => {
      const visibleRiverId = rivers.some(
        (river) => river.river_id === river_id && isVisibleRiver(river),
      ) ? river_id : null;
      onSelectRiver(visibleRiverId);
    });

  const cycleRiver = (backwards: boolean) => {
    const next = cycleId(rivers.filter(isVisibleRiver).map((r) => r.river_id), selectedRiverId, backwards);
    if (next != null) onSelectRiver(next);
  };

  return <InspectorCanvas {...view} canvasRef={canvasRef} onFrame={drawRivers} onPick={pickRiver} onCycle={cycleRiver} />;
}
