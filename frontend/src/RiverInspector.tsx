import { useEffect, useRef } from "react";
import type { KeyboardEvent } from "react";
import type { Projection, RiverSummary, Segment } from "./api";
import { fetchRiverAt } from "./api";
import type { Mat3, RenderTransform, Vec3 } from "./rotation";
import { getRenderTransform, latLonToXyz, matApply, matTranspose, project, toPixels, unproject, wrapLongitudeNear, xyzToLatLon } from "./rotation";
import { useRotationDrag } from "./rotationDrag";

interface Props {
  rivers: RiverSummary[];
  coastlineSegments: Segment[];
  width: number;
  height: number;
  displayWidth: number;
  displayHeight: number;
  projection: Projection;
  // The currently *committed* view rotation -- same rotation state MapCanvas/PlateInspector
  // share, so switching between views preserves orientation.
  rotation: Mat3;
  selectedRiverId: number | null;
  onSelectRiver: (id: number | null) => void;
  onRotationPreview: (latDeg: number, lonDeg: number) => void;
  onRotationCommitted: (rotation: Mat3) => void;
  // Inert rotate-drag + click while a background animation holds the world lock (see App.tsx).
  interactionDisabled?: boolean;
}

const BACKGROUND = "#0b1020";
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

// Matches backend render_image.py's COASTLINE_COLOR_RGB/COASTLINE_HALO_RGB -- same
// dark-halo-plus-light-line stroke, drawn there for the same reason it's needed here: this
// view has no filled backdrop at all, so without a coastline there's no land/ocean cue
// whatsoever, unlike the elevation/plates views' own hypsometric coloring.
const COASTLINE_RGB = "235, 235, 235";
const COASTLINE_HALO_RGB = "15, 15, 15";

// Renders every river network's flow-edge segments (see backend app/hydrology.py's
// group_rivers and GET /world/rivers) as an interactive display -- no server-baked PNG for
// this view, same "raw JSON, client draws it" philosophy as PlateInspector.tsx. Every river
// is always drawn, dim; the selected one is drawn last, bright, plus a ring at its own mouth
// (colored by endpoint type). Reuses the same long-press-drag rotate gesture as MapCanvas/
// PlateInspector (see rotationDrag.ts) plus click-to-select and Tab/Shift+Tab to cycle rivers.
export default function RiverInspector({
  rivers, coastlineSegments, width, height, displayWidth, displayHeight, projection, rotation,
  selectedRiverId, onSelectRiver, onRotationPreview, onRotationCommitted, interactionDisabled,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);

  // Projects one flow edge (two *true* world points), unwrapping the second point's
  // longitude relative to the first's own (rotated) longitude -- the two-point analogue of
  // PlateInspector's projectLoop (which unwraps a whole loop relative to its own center),
  // and the same technique backend render_image.py's _draw_rivers uses server-side
  // (_project_offset) for exactly this reason: a naive independent projection of both
  // endpoints can bow a short edge all the way across the map at the antimeridian seam.
  const projectSegment = (
    a: Vec3, b: Vec3, previewRotation: Mat3, transform: RenderTransform,
  ): [[number, number], [number, number]] => {
    const ra = matApply(previewRotation, a);
    const latA = Math.asin(Math.min(1, Math.max(-1, ra[2])));
    const lonA = Math.atan2(ra[1], ra[0]);
    const rb = matApply(previewRotation, b);
    const latB = Math.asin(Math.min(1, Math.max(-1, rb[2])));
    const lonB = wrapLongitudeNear(Math.atan2(rb[1], rb[0]), lonA);
    const [xA, yA] = project(projection, latA, lonA);
    const [xB, yB] = project(projection, latB, lonB);
    return [toPixels(transform, xA, yA), toPixels(transform, xB, yB)];
  };

  const drawRivers = (previewRotation: Mat3) => {
    const canvas = canvasRef.current;
    const ctx = canvas?.getContext("2d");
    if (!canvas || !ctx) return;
    ctx.fillStyle = BACKGROUND;
    ctx.fillRect(0, 0, width, height);

    const transform = getRenderTransform(projection, width, height);
    const pixelScale = width / 1100;
    const lineWidth = Math.max(1, pixelScale);

    const strokeSegments = (segments: [Vec3, Vec3][], color: string, width: number) => {
      ctx.strokeStyle = color;
      ctx.lineWidth = width;
      for (const [a, b] of segments) {
        const [[x1, y1], [x2, y2]] = projectSegment(a, b, previewRotation, transform);
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
      }
    };
    // Halo first, then a lighter line on top -- see COASTLINE_HALO_RGB's own comment.
    strokeSegments(coastlineSegments, `rgba(${COASTLINE_HALO_RGB}, 1.0)`, lineWidth * 2.6);
    strokeSegments(coastlineSegments, `rgba(${COASTLINE_RGB}, 1.0)`, lineWidth * 1.1);

    const drawOne = (river: RiverSummary, selected: boolean) => {
      const alpha = selected ? SELECTED_RIVER_ALPHA : RIVER_ALPHA;
      ctx.strokeStyle = `rgba(${RIVER_RGB}, ${alpha})`;
      ctx.lineWidth = selected ? lineWidth * 2.5 : lineWidth;
      for (const [a, b] of river.segments) {
        const [[x1, y1], [x2, y2]] = projectSegment(a, b, previewRotation, transform);
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
      }

      if (!selected) return;
      const r = matApply(previewRotation, river.mouth_xyz);
      const lat = Math.asin(Math.min(1, Math.max(-1, r[2])));
      const lon = Math.atan2(r[1], r[0]);
      const [x, y] = project(projection, lat, lon);
      const [px, py] = toPixels(transform, x, y);
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

    for (const river of rivers) {
      if (river.river_id === selectedRiverId) continue;
      drawOne(river, false);
    }
    const selected = rivers.find((r) => r.river_id === selectedRiverId);
    if (selected) drawOne(selected, true);
  };

  useEffect(() => {
    drawRivers(rotation);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [rivers, coastlineSegments, selectedRiverId, projection, rotation, width, height]);

  const handleClick = (backingX: number, backingY: number) => {
    const transform = getRenderTransform(projection, width, height);
    const x = (backingX - transform.offsetX) / transform.scale;
    const y = -(backingY - transform.offsetY) / transform.scale;
    const latLon = unproject(projection, x, y);
    if (!latLon) return; // clicked outside the sphere's projected silhouette
    const displayXyz = latLonToXyz(latLon[0], latLon[1]); // in the *display* (rotated) frame
    const trueXyz = matApply(matTranspose(rotation), displayXyz); // transpose == inverse for a rotation matrix
    const [trueLat, trueLon] = xyzToLatLon(trueXyz);
    fetchRiverAt((trueLat * 180) / Math.PI, (trueLon * 180) / Math.PI)
      .then(({ river_id }) => {
        onSelectRiver(river_id);
        containerRef.current?.focus();
      })
      .catch(() => {
        // A click that races a generate/step (no world, or a stale river list) -- ignored,
        // matching how PlateInspector's own click handler fails silently for the same reason.
      });
  };

  useRotationDrag({
    elementRef: canvasRef,
    width, height, displayWidth, displayHeight, projection, rotation,
    onFrame: drawRivers,
    onRotationPreview,
    onRotationCommitted,
    onClick: handleClick,
    disabled: interactionDisabled,
  });

  useEffect(() => {
    containerRef.current?.focus();
  }, []);

  const handleKeyDown = (e: KeyboardEvent<HTMLDivElement>) => {
    if (e.key !== "Tab") return;
    e.preventDefault();
    if (rivers.length === 0) return;
    const ids = rivers.map((r) => r.river_id).sort((a, b) => a - b);
    const currentIndex = selectedRiverId == null ? -1 : ids.indexOf(selectedRiverId);
    const delta = e.shiftKey ? -1 : 1;
    const nextIndex = currentIndex === -1 ? (e.shiftKey ? ids.length - 1 : 0) : (currentIndex + delta + ids.length) % ids.length;
    onSelectRiver(ids[nextIndex]);
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
