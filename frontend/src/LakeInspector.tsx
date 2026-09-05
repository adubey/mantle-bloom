import { useEffect, useRef } from "react";
import type { KeyboardEvent } from "react";
import type { LakeAtResponse, LakeSummary, Projection, Segment } from "./api";
import { fetchLakeAt } from "./api";
import type { Mat3, Vec3 } from "./rotation";
import { getRenderTransform, latLonToXyz, matApply, matTranspose, project, toPixels, unproject, wrapLongitudeNear, xyzToLatLon } from "./rotation";
import { useRotationDrag } from "./rotationDrag";

interface Props {
  lakes: LakeSummary[];
  coastlineSegments: Segment[];
  width: number;
  height: number;
  displayWidth: number;
  displayHeight: number;
  projection: Projection;
  // The currently *committed* view rotation -- same rotation state MapCanvas/PlateInspector/
  // RiverInspector share, so switching between views preserves orientation.
  rotation: Mat3;
  // Whatever's currently displayed -- either one of `lakes` (selectedBasin.is_lake) or a dry
  // basin returned by a land click (not itself a member of `lakes` -- see api.ts's
  // LakeAtResponse). `null` means nothing selected yet.
  selectedBasin: LakeSummary | null;
  onSelect: (kind: LakeAtResponse["kind"], basin: LakeSummary | null) => void;
  onRotationPreview: (latDeg: number, lonDeg: number) => void;
  onRotationCommitted: (rotation: Mat3) => void;
  // Inert rotate-drag + click while a background animation holds the world lock (see App.tsx).
  interactionDisabled?: boolean;
}

const BACKGROUND = "#0b1020";
// Matches backend render_image.py's COASTLINE_COLOR_RGB/COASTLINE_HALO_RGB -- same
// dark-halo-plus-light-line stroke RiverInspector.tsx uses for the same reason: this view has
// no filled backdrop, so without a coastline there's no land/ocean cue at all.
const COASTLINE_RGB = "235, 235, 235";
const COASTLINE_HALO_RGB = "15, 15, 15";

// Lake water -- close to render_image.py's own LAKE_COLOR_RGB in hue, brightened for
// visibility as scattered dots against BACKGROUND rather than a filled hypsometric wash.
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

// Renders every currently-visible lake (see backend app/main.py's GET /world/lakes) as an
// interactive point-cloud display -- no server-baked PNG, same "raw JSON, client draws it"
// philosophy as PlateInspector/RiverInspector. Every lake is always drawn, dim; the selected
// basin (a lake from `lakes`, or a dry basin returned by a land click -- see LakeAtResponse) is
// drawn last, bright, plus its own floor and outlet markers, a ring at every inflowing
// river's mouth, and (when it's actively spilling and that spill has grown into a real river
// network -- see api.ts's OutflowRiver) its own live outflow drawn as a full river path, not
// just another ring. Reuses the same long-press-drag rotate gesture as the other inspectors, plus
// click-to-select and Tab/Shift+Tab to cycle through `lakes` specifically (a dry basin has no
// enumerable identity to cycle through -- see this module's own Tab handler).
export default function LakeInspector({
  lakes, coastlineSegments, width, height, displayWidth, displayHeight, projection, rotation,
  selectedBasin, onSelect, onRotationPreview, onRotationCommitted, interactionDisabled,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);

  const drawLakes = (previewRotation: Mat3) => {
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

    const toScreen = (p: Vec3): [number, number] => {
      const r = matApply(previewRotation, p);
      const lat = Math.asin(Math.min(1, Math.max(-1, r[2])));
      const lon = Math.atan2(r[1], r[0]);
      const [x, y] = project(projection, lat, lon);
      return toPixels(transform, x, y);
    };

    // Projects one coastline edge (two *true* world points a real, short 3D hop apart),
    // unwrapping the second point's longitude relative to the first's own (rotated) longitude
    // -- same technique RiverInspector.tsx's own projectSegment uses, and for the same reason:
    // independently projecting each endpoint with plain toScreen (as this used to) can bow a
    // short edge all the way across the map at the antimeridian seam once the view is rotated,
    // drawing a spurious line clear across the whole canvas.
    const projectSegment = (a: Vec3, b: Vec3): [[number, number], [number, number]] => {
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

    const strokeSegments = (segments: [Vec3, Vec3][], color: string, width: number) => {
      ctx.strokeStyle = color;
      ctx.lineWidth = width;
      for (const [a, b] of segments) {
        const [[x1, y1], [x2, y2]] = projectSegment(a, b);
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
      }
    };
    strokeSegments(coastlineSegments, `rgba(${COASTLINE_HALO_RGB}, 1.0)`, lineWidth * 2.6);
    strokeSegments(coastlineSegments, `rgba(${COASTLINE_RGB}, 1.0)`, lineWidth * 1.1);

    const drawPoints = (points: Vec3[], color: string) => {
      ctx.fillStyle = color;
      for (const p of points) {
        const [px, py] = toScreen(p);
        ctx.fillRect(px - pointRadius, py - pointRadius, pointSize, pointSize);
      }
    };

    const drawRing = (p: Vec3, color: string) => {
      const [px, py] = toScreen(p);
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
        strokeSegments(selectedBasin.outflow_river.segments, `rgba(${OUTFLOW_RIVER_RGB}, 0.95)`, lineWidth * 2.5);
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

  const handleClick = (backingX: number, backingY: number) => {
    const transform = getRenderTransform(projection, width, height);
    const x = (backingX - transform.offsetX) / transform.scale;
    const y = -(backingY - transform.offsetY) / transform.scale;
    const latLon = unproject(projection, x, y);
    if (!latLon) return; // clicked outside the sphere's projected silhouette
    const displayXyz = latLonToXyz(latLon[0], latLon[1]); // in the *display* (rotated) frame
    const trueXyz = matApply(matTranspose(rotation), displayXyz); // transpose == inverse for a rotation matrix
    const [trueLat, trueLon] = xyzToLatLon(trueXyz);
    fetchLakeAt((trueLat * 180) / Math.PI, (trueLon * 180) / Math.PI)
      .then(({ kind, basin }) => {
        onSelect(kind, basin);
        containerRef.current?.focus();
      })
      .catch(() => {
        // A click that races a generate/step (no world, or a stale lake list) -- ignored,
        // matching how PlateInspector/RiverInspector's own click handlers fail silently for
        // the same reason.
      });
  };

  useRotationDrag({
    elementRef: canvasRef,
    width, height, displayWidth, displayHeight, projection, rotation,
    onFrame: drawLakes,
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
    if (lakes.length === 0) return;
    const ids = lakes.map((l) => l.lake_id as number).sort((a, b) => a - b);
    const currentId = selectedBasin && selectedBasin.is_lake ? selectedBasin.lake_id : null;
    const currentIndex = currentId == null ? -1 : ids.indexOf(currentId);
    const delta = e.shiftKey ? -1 : 1;
    const nextIndex = currentIndex === -1 ? (e.shiftKey ? ids.length - 1 : 0) : (currentIndex + delta + ids.length) % ids.length;
    const nextLake = lakes.find((l) => l.lake_id === ids[nextIndex]) ?? null;
    onSelect("lake", nextLake);
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
