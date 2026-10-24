import { useEffect, useRef } from "react";
import type { LakeAtResponse, LakeSummary, Segment } from "./api";
import { fetchLakeAt } from "./api";
import InspectorCanvas from "./InspectorCanvas";
import type { InspectorViewProps } from "./InspectorCanvas";
import type { Mat3, Vec3 } from "./rotation";
import { beginInspectorFrame, strokeCoastline, strokeSegments } from "./inspectorDrawing";
import { cycleId } from "./screenProjector";

interface Props extends InspectorViewProps {
  worldId: string;
  lakes: LakeSummary[];
  coastlineSegments: Segment[];
  // Whatever's currently displayed -- either one of `lakes` (selectedBasin.is_lake) or a dry
  // basin returned by a land click (not itself a member of `lakes` -- see api.ts's
  // LakeAtResponse). `null` means nothing selected yet.
  selectedBasin: LakeSummary | null;
  onSelect: (kind: LakeAtResponse["kind"], basin: LakeSummary | null) => void;
}

// Lake water -- close to render_image.py's own LAKE_COLOR_RGB in hue, brightened for
// visibility as scattered dots against INSPECTOR_BACKGROUND rather than a filled hypsometric wash.
const LAKE_RGB = "90, 170, 255";
const LAKE_ALPHA = 0.16; // other (non-selected) lakes -- "less visible"
const SELECTED_LAKE_ALPHA = 0.9;
// A selected *dry* basin (no visible water) -- same amber RiverInspector's own
// MOUTH_RING_COLOR uses for a river's "other" (dry interior sink) mouth type, so the color
// already reads as "dry land, not water" from that view.
const BASIN_RGB = "230, 170, 70";
const SELECTED_BASIN_ALPHA = 0.9;
const POINT_RADIUS_PX = 1.6;

// The basin's own lowest point -- always drawn for whatever's selected, lake or dry basin.
const FLOOR_RING_COLOR = "rgba(230, 235, 245, 1.0)";
// "The lowest point of the edge of the basin" -- where a river out of it would source from.
const OUTLET_RING_COLOR = "rgba(255, 205, 70, 1.0)";
const INFLOW_RING_COLOR = "rgba(120, 220, 190, 1.0)"; // matches RiverInspector's own "lake" mouth-ring color
const MARKER_RADIUS_PX = 6;

// Matches RiverInspector.tsx's own RIVER_RGB (in turn render_image.py's RIVER_COLOR_RGB) --
// the selected basin's outflow reads as "a river," the same color everywhere else in this
// codebase draws one, rather than inventing a fourth meaning for a color already spoken for.
const OUTFLOW_RIVER_RGB = "77, 216, 230";

// Renders every currently-visible lake (see backend app/main.py's GET /worlds/{world_id}/lakes) as an
// interactive point-cloud display -- no server-baked PNG, same "raw JSON, client draws it"
// philosophy as PlateInspector/RiverInspector. Every lake is always drawn, dim; the selected
// basin (a lake from `lakes`, or a dry basin returned by a land click -- see LakeAtResponse) is
// drawn last, bright, plus its own floor and outlet markers, a ring at every inflowing
// river's mouth, and (when it's actively spilling and that spill has grown into a real river
// network -- see api.ts's OutflowRiver) its own live outflow drawn as a full river path, not
// just another ring. Reuses the same long-press-drag rotate gesture as the other inspectors, plus
// click-to-select and Tab/Shift+Tab to cycle through `lakes` specifically (a dry basin has no
// enumerable identity to cycle through -- see this module's own Tab handler).
export default function LakeInspector({ worldId, lakes, coastlineSegments, selectedBasin, onSelect, ...view }: Props) {
  const { width, height, projection, rotation } = view;
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const drawLakes = (previewRotation: Mat3) => {
    const frame = beginInspectorFrame(canvasRef.current, projection, previewRotation, width, height);
    if (!frame) return;
    const { ctx, projector, pixelScale, lineWidth } = frame;
    const pointRadius = POINT_RADIUS_PX * pixelScale;
    const pointSize = pointRadius * 2;

    strokeCoastline(ctx, projector, coastlineSegments, lineWidth);

    const drawPoints = (points: Vec3[], color: string) => {
      ctx.fillStyle = color;
      for (const p of points) {
        const [px, py] = projector.point(p);
        ctx.fillRect(px - pointRadius, py - pointRadius, pointSize, pointSize);
      }
    };

    const drawRing = (p: Vec3, color: string) => {
      const [px, py] = projector.point(p);
      const r = MARKER_RADIUS_PX * pixelScale;
      ctx.strokeStyle = color;
      ctx.lineWidth = lineWidth * 2;
      ctx.beginPath();
      ctx.arc(px, py, r, 0, 2 * Math.PI);
      ctx.stroke();
      ctx.fillStyle = color;
      ctx.beginPath();
      ctx.arc(px, py, r * 0.3, 0, 2 * Math.PI);
      ctx.fill();
    };

    for (const lake of lakes) {
      if (selectedBasin && selectedBasin.is_lake && lake.lake_id === selectedBasin.lake_id) continue;
      drawPoints(lake.member_xyz, `rgba(${LAKE_RGB}, ${LAKE_ALPHA})`);
    }

    if (selectedBasin) {
      const [rgb, alpha] = selectedBasin.is_lake ? [LAKE_RGB, SELECTED_LAKE_ALPHA] : [BASIN_RGB, SELECTED_BASIN_ALPHA];
      drawPoints(selectedBasin.member_xyz, `rgba(${rgb}, ${alpha})`);
      // Drawn before the rings below (so the outlet marker still reads clearly on top of it)
      // and bold -- the same "selected" line width RiverInspector.tsx's own drawOne uses -- so
      // the basin's actual live outflow reads as a real river leaving the lake, not just an
      // outlet marker with no channel drawn at all.
      if (selectedBasin.outflow_river) {
        strokeSegments(ctx, projector, selectedBasin.outflow_river.segments, `rgba(${OUTFLOW_RIVER_RGB}, 0.95)`, lineWidth * 2.5);
      }
      for (const river of selectedBasin.inflow_rivers) {
        drawRing(river.mouth_xyz, INFLOW_RING_COLOR);
      }
      if (selectedBasin.outlet_xyz) drawRing(selectedBasin.outlet_xyz, OUTLET_RING_COLOR);
      drawRing(selectedBasin.floor_xyz, FLOOR_RING_COLOR);
    }
  };

  useEffect(() => {
    drawLakes(rotation);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [lakes, coastlineSegments, selectedBasin, projection, rotation, width, height]);

  const pickBasin = (latDeg: number, lonDeg: number) =>
    fetchLakeAt(worldId, latDeg, lonDeg).then(({ kind, basin }) => onSelect(kind, basin));

  // Cycles through `lakes` only: a dry basin has no enumerable identity to cycle through.
  const cycleLake = (backwards: boolean) => {
    const currentId = selectedBasin && selectedBasin.is_lake ? selectedBasin.lake_id : null;
    const next = cycleId(lakes.map((l) => l.lake_id as number), currentId, backwards);
    if (next == null) return;
    onSelect("lake", lakes.find((l) => l.lake_id === next) ?? null);
  };

  return <InspectorCanvas {...view} canvasRef={canvasRef} onFrame={drawLakes} onPick={pickBasin} onCycle={cycleLake} />;
}
