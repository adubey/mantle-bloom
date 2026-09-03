import { useEffect, useRef } from "react";
import type { KeyboardEvent } from "react";
import type { FaultSummary, Projection, Segment } from "./api";
import { fetchFaultAt } from "./api";
import type { Mat3, RenderTransform, Vec3 } from "./rotation";
import { getRenderTransform, latLonToXyz, matApply, matTranspose, project, toPixels, unproject, wrapLongitudeNear, xyzToLatLon } from "./rotation";
import { useRotationDrag } from "./rotationDrag";

interface Props {
  faults: FaultSummary[];
  coastlineSegments: Segment[];
  width: number;
  height: number;
  displayWidth: number;
  displayHeight: number;
  projection: Projection;
  // The currently *committed* view rotation -- same rotation state MapCanvas/PlateInspector
  // share, so switching between views preserves orientation.
  rotation: Mat3;
  selectedFaultId: number | null;
  onSelectFault: (id: number | null) => void;
  onRotationPreview: (latDeg: number, lonDeg: number) => void;
  onRotationCommitted: (rotation: Mat3) => void;
}

const BACKGROUND = "#0b1020";

// Matches backend render_image.py's _ELEV_REASON_RGB fault stops (codes 15/16/17) and the
// legendData.ts entries, so a fault reads the same color here as baked into the "Last
// elevation change" view.
const KIND_RGB: Record<FaultSummary["kind"], string> = {
  normal: "120, 190, 90",
  reverse: "176, 60, 90",
  strike_slip: "230, 190, 70",
};
const KIND_LABEL: Record<FaultSummary["kind"], string> = {
  normal: "normal (graben)",
  reverse: "reverse (thrust)",
  strike_slip: "strike-slip",
};

const ACTIVE_ALPHA = 0.85;
const SCAR_ALPHA = 0.3; // locked-up faults -- present but recessive
const SELECTED_ALPHA = 1.0;

// Matches backend render_image.py's COASTLINE_COLOR_RGB/COASTLINE_HALO_RGB -- same reason
// RiverInspector needs it: this view draws no filled backdrop, so without a coastline there's
// no land/ocean cue at all.
const COASTLINE_RGB = "235, 235, 235";
const COASTLINE_HALO_RGB = "15, 15, 15";

// Renders every intraplate fault trace (see backend app/faults.py and GET /world/faults) as
// an interactive display -- no server-baked PNG, same "raw JSON, client draws it" philosophy
// as RiverInspector/PlateInspector. Active faults are drawn solid and colored by regime
// (normal/reverse/strike-slip); locked-up faults ("scars") are drawn dashed and dim. The
// selected fault is drawn last, bright and thick, with its endpoints ringed. Reuses the same
// long-press-drag rotate gesture as the other inspectors (rotationDrag.ts) plus
// click-to-select and Tab/Shift+Tab to cycle faults.
export default function FaultInspector({
  faults, coastlineSegments, width, height, displayWidth, displayHeight, projection, rotation,
  selectedFaultId, onSelectFault, onRotationPreview, onRotationCommitted,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);

  // Projects one edge (two *true* world points), unwrapping the second point's longitude
  // relative to the first's -- see RiverInspector.projectSegment for why a naive independent
  // projection bows an edge across the antimeridian seam.
  const projectEdge = (
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

  const draw = (previewRotation: Mat3) => {
    const canvas = canvasRef.current;
    const ctx = canvas?.getContext("2d");
    if (!canvas || !ctx) return;
    ctx.fillStyle = BACKGROUND;
    ctx.fillRect(0, 0, width, height);

    const transform = getRenderTransform(projection, width, height);
    const pixelScale = width / 1100;
    const lineWidth = Math.max(1, pixelScale);

    const strokeEdges = (edges: [Vec3, Vec3][], color: string, w: number, dash: number[] = []) => {
      ctx.strokeStyle = color;
      ctx.lineWidth = w;
      ctx.setLineDash(dash);
      for (const [a, b] of edges) {
        const [[x1, y1], [x2, y2]] = projectEdge(a, b, previewRotation, transform);
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
      }
      ctx.setLineDash([]);
    };

    // Halo first, then a lighter line on top -- see COASTLINE_HALO_RGB's own comment.
    strokeEdges(coastlineSegments, `rgba(${COASTLINE_HALO_RGB}, 1.0)`, lineWidth * 2.6);
    strokeEdges(coastlineSegments, `rgba(${COASTLINE_RGB}, 1.0)`, lineWidth * 1.1);

    const traceEdges = (fault: FaultSummary): [Vec3, Vec3][] => {
      const edges: [Vec3, Vec3][] = [];
      for (let i = 0; i + 1 < fault.trace.length; i++) {
        edges.push([fault.trace[i] as Vec3, fault.trace[i + 1] as Vec3]);
      }
      return edges;
    };

    const projectPoint = (v: Vec3): [number, number] => {
      const r = matApply(previewRotation, v);
      const lat = Math.asin(Math.min(1, Math.max(-1, r[2])));
      const lon = Math.atan2(r[1], r[0]);
      const [x, y] = project(projection, lat, lon);
      return toPixels(transform, x, y);
    };

    const drawOne = (fault: FaultSummary, selected: boolean) => {
      const alpha = selected ? SELECTED_ALPHA : fault.active ? ACTIVE_ALPHA : SCAR_ALPHA;
      const w = selected ? lineWidth * 3.2 : fault.active ? lineWidth * 2 : lineWidth * 1.3;
      const dash = fault.active ? [] : [lineWidth * 3, lineWidth * 3];
      const rgb = KIND_RGB[fault.kind];
      strokeEdges(traceEdges(fault), `rgba(${rgb}, ${alpha})`, w, dash);

      // A midpoint mark, so a fault whose whole trace projects to a few pixels at world
      // scale still registers -- filled for an active fault, a hollow ring for a scar.
      const [mx, my] = projectPoint(fault.trace[Math.floor(fault.trace.length / 2)] as Vec3);
      const r = (selected ? 4.5 : 2.6) * pixelScale;
      ctx.lineWidth = lineWidth * 1.2;
      ctx.strokeStyle = `rgba(${rgb}, ${Math.min(1, alpha + 0.15)})`;
      ctx.fillStyle = `rgba(${rgb}, ${alpha})`;
      ctx.beginPath();
      ctx.arc(mx, my, r, 0, 2 * Math.PI);
      if (fault.active || selected) ctx.fill();
      ctx.stroke();

      if (!selected) return;
      for (const end of [fault.trace[0], fault.trace[fault.trace.length - 1]]) {
        const [px, py] = projectPoint(end as Vec3);
        ctx.strokeStyle = `rgba(${rgb}, 1.0)`;
        ctx.lineWidth = lineWidth * 2;
        ctx.beginPath();
        ctx.arc(px, py, 6 * pixelScale, 0, 2 * Math.PI);
        ctx.stroke();
      }
    };

    // Scars first (recessive), then active faults, then the selected one on top.
    for (const f of faults) {
      if (f.fault_id === selectedFaultId || f.active) continue;
      drawOne(f, false);
    }
    for (const f of faults) {
      if (f.fault_id === selectedFaultId || !f.active) continue;
      drawOne(f, false);
    }
    const selected = faults.find((f) => f.fault_id === selectedFaultId);
    if (selected) drawOne(selected, true);
  };

  useEffect(() => {
    draw(rotation);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [faults, coastlineSegments, selectedFaultId, projection, rotation, width, height]);

  const handleClick = (backingX: number, backingY: number) => {
    const transform = getRenderTransform(projection, width, height);
    const x = (backingX - transform.offsetX) / transform.scale;
    const y = -(backingY - transform.offsetY) / transform.scale;
    const latLon = unproject(projection, x, y);
    if (!latLon) return; // clicked outside the sphere's projected silhouette
    const displayXyz = latLonToXyz(latLon[0], latLon[1]); // in the *display* (rotated) frame
    const trueXyz = matApply(matTranspose(rotation), displayXyz); // transpose == inverse for a rotation matrix
    const [trueLat, trueLon] = xyzToLatLon(trueXyz);
    fetchFaultAt((trueLat * 180) / Math.PI, (trueLon * 180) / Math.PI)
      .then(({ fault_id }) => {
        onSelectFault(fault_id);
        containerRef.current?.focus();
      })
      .catch(() => {
        // A click that races a generate/step -- ignored, matching the other inspectors.
      });
  };

  useRotationDrag({
    elementRef: canvasRef,
    width, height, displayWidth, displayHeight, projection, rotation,
    onFrame: draw,
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
    if (faults.length === 0) return;
    const ids = faults.map((f) => f.fault_id).sort((a, b) => a - b);
    const currentIndex = selectedFaultId == null ? -1 : ids.indexOf(selectedFaultId);
    const delta = e.shiftKey ? -1 : 1;
    const nextIndex = currentIndex === -1 ? (e.shiftKey ? ids.length - 1 : 0) : (currentIndex + delta + ids.length) % ids.length;
    onSelectFault(ids[nextIndex]);
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

export { KIND_LABEL as FAULT_KIND_LABEL };
