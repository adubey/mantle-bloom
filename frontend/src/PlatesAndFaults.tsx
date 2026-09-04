import { useEffect, useRef } from "react";
import type { KeyboardEvent } from "react";
import type {
  EarthquakeSummary, FaultKind, FaultSummary, FaultSystemSummary, PlateSummary, Projection, Segment, VolcanoSummary,
} from "./api";
import { fetchPlateAt } from "./api";
import type { Mat3, RenderTransform, Vec3 } from "./rotation";
import {
  getRenderTransform, latLonToXyz, matApply, matTranspose, project, rotationMatrix, toPixels, unproject,
  wrapLongitudeNear, xyzToLatLon,
} from "./rotation";
import { useRotationDrag } from "./rotationDrag";
import { plateColor } from "./platePalette";

interface Props {
  plates: PlateSummary[];
  faults: FaultSummary[];
  faultSystems: FaultSystemSummary[];
  earthquakes: EarthquakeSummary[];
  volcanoes: VolcanoSummary[];
  coastlineSegments: Segment[];
  // The "Earthquakes & volcanoes" checkbox in App's sidebar -- gates the two transient
  // overlays together (both are seismic/magmatic activity markers, distinct from the standing
  // plate + fault geometry the view always shows).
  showQuakesVolcanoes: boolean;
  // Set when a fault-type row in the legend is clicked (see Legend.tsx / legendData.ts's
  // faultKindForLegendLabel): that kind's strands + systems are haloed and everything else on
  // the map is dimmed right back, so you can see where that regime sits. null = no isolation.
  highlightedFaultKind: FaultKind | null;
  width: number;
  height: number;
  displayWidth: number;
  displayHeight: number;
  projection: Projection;
  // The currently *committed* view rotation -- same rotation state MapCanvas/the other
  // inspectors share, so switching between views preserves orientation.
  rotation: Mat3;
  selectedPlateId: number | null;
  onSelectPlate: (id: number | null) => void;
  onRotationPreview: (latDeg: number, lonDeg: number) => void;
  onRotationCommitted: (rotation: Mat3) => void;
  // Inert rotate-drag + click while a background animation holds the world lock (see App.tsx).
  interactionDisabled?: boolean;
}

const BACKGROUND = "#0b1020";

// --- Plate drawing --- just the outline loop now (no bounding ellipse, no per-node dots): the
// view is the plate + fault geometry over a coastline, nothing else. Non-selected outlines sit
// at OUTLINE_ALPHA in each plate's own palette colour, the selected one opaque and thicker.
const OUTLINE_ALPHA = 0.6;
const OUTLINE_ALPHA_SELECTED = 1.0;
const SEGMENT_BREAK_FACTOR = 6;

// The selected plate's Euler pole + motion arc, and a bright cyan distinct from every plate
// palette hue and from the three fault regime colours.
const MOTION_RGB = "125, 225, 255";
// Δθ swept about the Euler pole per (cm/yr) of plate speed -- scaled so a typical few-cm/yr
// plate draws a clearly readable arc, clamped so a MAX-railed plate can't wrap the globe.
const MOTION_ARC_DEG_PER_CM_YR = 8;
const MOTION_ARC_MIN_DEG = 6;
const MOTION_ARC_MAX_DEG = 95;

// --- Fault drawing (same regime colours / earthquake overlay the old Fault Lines view used) ---
const KIND_RGB: Record<FaultSummary["kind"], string> = {
  normal: "120, 190, 90",
  reverse: "176, 60, 90",
  strike_slip: "230, 190, 70",
};
const ACTIVE_ALPHA = 0.85;
const SCAR_ALPHA = 0.3;

const EARTHQUAKE_RGB = "255, 180, 70";
const EARTHQUAKE_RETAIN_MYR = 5.0;
// Warm magmatic orange, distinct from the earthquake amber -- an upward triangle (a cinder
// cone), filled when the vent still has eruption potential, hollow once dormant.
const VOLCANO_RGB = "255, 120, 40";

// A cool blue-grey -- deliberately not white, so it reads as clearly distinct from the plate
// outlines drawn over it (which carry each plate's own saturated palette hue).
const COASTLINE_RGB = "150, 170, 200";
const COASTLINE_HALO_RGB = "15, 15, 15";

// Unit-vector centroid of a loop of world points -- the anchor for the motion arc and the
// longitude-unwrap centre for projecting the outline.
const loopCentroid = (pts: Vec3[]): Vec3 => {
  let x = 0, y = 0, z = 0;
  for (const p of pts) {
    x += p[0];
    y += p[1];
    z += p[2];
  }
  const n = Math.hypot(x, y, z) || 1;
  return [x / n, y / n, z / n];
};

// The "Plates & Faults" view (replaces the old PNG "Plates" view + "Fault lines" inspector):
// plate outline loops with click-to-select + Tab cycling and a metadata panel (see App.tsx),
// the selected plate's Euler pole + a speed-scaled motion arc, the intraplate fault systems +
// strands drawn for context (not individually selectable here -- selecting a plate emphasises
// its own faults instead; clicking a fault type in the legend isolates that regime), and an
// earthquake + volcano activity overlay toggled together from the sidebar. All client-drawn,
// same "raw JSON, the client renders it" approach as the other inspectors.
export default function PlatesAndFaults({
  plates, faults, faultSystems, earthquakes, volcanoes, coastlineSegments, showQuakesVolcanoes,
  highlightedFaultKind,
  width, height, displayWidth, displayHeight, projection, rotation,
  selectedPlateId, onSelectPlate, onRotationPreview, onRotationCommitted, interactionDisabled,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const containerRef = useRef<HTMLDivElement>(null);

  const selectedSystemIds = new Set<number>(
    faults
      .filter((f) => f.plate_id === selectedPlateId && f.system_id != null)
      .map((f) => f.system_id as number),
  );

  // Projects a closed loop of *true* world points, unwrapping every point's longitude
  // relative to the loop's own (rotated) center -- see PlateInspector.projectLoop.
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

  // Projects one edge (two *true* world points), unwrapping the second relative to the first
  // -- see RiverInspector.projectSegment for why a naive independent projection bows an edge
  // across the antimeridian seam.
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
    // A fault type is isolated from the legend -- dim the standing geometry (coastline, plate
    // outlines, other regimes, the activity overlay) so that regime reads clearly on top.
    const hk = highlightedFaultKind;
    const contextDim = hk ? 0.32 : 1;
    // The transient activity overlay isn't a fault regime -- push it further back than the
    // standing geometry while a regime is isolated so it doesn't compete.
    const activityDim = hk ? 0.12 : 1;
    // A pixel jump longer than this between consecutive projected points is the antimeridian
    // seam wrapping, not a real edge -- skip it (see strokeRobustLoop's own break logic).
    const seamJump = Math.max(40 * pixelScale, width * 0.25);

    const projectPoint = (v: Vec3): [number, number] => {
      const r = matApply(previewRotation, v);
      const lat = Math.asin(Math.min(1, Math.max(-1, r[2])));
      const lon = Math.atan2(r[1], r[0]);
      const [x, y] = project(projection, lat, lon);
      return toPixels(transform, x, y);
    };

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

    const strokePolyline = (trace: [number, number, number][], color: string, w: number, dash: number[] = []) => {
      if (trace.length < 2) return;
      ctx.strokeStyle = color;
      ctx.lineWidth = w;
      ctx.lineJoin = "round";
      ctx.lineCap = "round";
      ctx.setLineDash(dash);
      ctx.beginPath();
      for (let i = 0; i + 1 < trace.length; i++) {
        const [[x1, y1], [x2, y2]] = projectEdge(trace[i] as Vec3, trace[i + 1] as Vec3, previewRotation, transform);
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
      }
      ctx.stroke();
      ctx.setLineDash([]);
      ctx.lineJoin = "miter";
      ctx.lineCap = "butt";
    };

    const traceEdges = (trace: [number, number, number][]): [Vec3, Vec3][] => {
      const edges: [Vec3, Vec3][] = [];
      for (let i = 0; i + 1 < trace.length; i++) edges.push([trace[i] as Vec3, trace[i + 1] as Vec3]);
      return edges;
    };

    // 1. Coastline for orientation (halo then a lighter line on top).
    strokeEdges(coastlineSegments, `rgba(${COASTLINE_HALO_RGB}, ${0.9 * activityDim})`, lineWidth * 2.6);
    strokeEdges(coastlineSegments, `rgba(${COASTLINE_RGB}, ${contextDim})`, lineWidth * 1.1);

    // 2. Plate outlines -- non-selected first, then the selected one opaque + thicker on top.
    const drawPlate = (plate: PlateSummary, selected: boolean) => {
      if (plate.outline.length === 0) return;
      const [r, g, b] = plateColor(plate.plate_id);
      const center = loopCentroid(plate.outline as Vec3[]);
      const pixels = projectLoop(plate.outline, center, previewRotation, transform);
      const alpha = (selected ? OUTLINE_ALPHA_SELECTED : OUTLINE_ALPHA) * contextDim;
      strokeRobustLoop(ctx, pixels, `rgba(${r}, ${g}, ${b}, ${alpha})`, selected ? lineWidth * 2 : lineWidth);
    };

    for (const plate of plates) {
      if (plate.plate_id === selectedPlateId) continue;
      drawPlate(plate, false);
    }
    const selectedPlate = plates.find((p) => p.plate_id === selectedPlateId);
    if (selectedPlate) drawPlate(selectedPlate, true);

    // 3. Fault-system master lineaments -- broad translucent belt + thin dashed centerline.
    for (const sys of faultSystems) {
      const rgb = KIND_RGB[sys.kind];
      const sel = selectedSystemIds.has(sys.system_id);
      const muted = hk != null && sys.kind !== hk;
      const beltBase = sys.active ? (sel ? 0.26 : 0.12) : 0.06;
      const beltAlpha = muted ? beltBase * 0.1 : hk != null ? Math.min(0.4, beltBase * 1.7) : beltBase;
      strokePolyline(sys.trace, `rgba(${rgb}, ${beltAlpha})`, lineWidth * (sel ? 24 : 16));
      const lineBase = sys.active ? (sel ? 0.85 : 0.45) : 0.22;
      strokePolyline(
        sys.trace, `rgba(${rgb}, ${muted ? lineBase * 0.12 : lineBase})`,
        lineWidth * (sel ? 1.8 : 1.1), [lineWidth * 6, lineWidth * 5],
      );
    }

    // 4. Fault strands -- scars (recessive) first, then active, then the isolated regime on top.
    const drawFault = (fault: FaultSummary) => {
      const onSelectedPlate = fault.plate_id === selectedPlateId;
      const muted = hk != null && fault.kind !== hk;
      const isolated = hk != null && fault.kind === hk;
      const alpha = (fault.active ? (onSelectedPlate ? 1.0 : ACTIVE_ALPHA) : SCAR_ALPHA) * (muted ? 0.1 : 1);
      const w = fault.active ? lineWidth * (onSelectedPlate ? 2.6 : 2) : lineWidth * 1.3;
      const dash = fault.active ? [] : [lineWidth * 3, lineWidth * 3];
      const edges = traceEdges(fault.trace);
      if (isolated) {
        strokeEdges(edges, `rgba(255, 255, 255, ${fault.active ? 0.5 : 0.25})`, w + lineWidth * 2.6, dash);
      }
      strokeEdges(edges, `rgba(${KIND_RGB[fault.kind]}, ${alpha})`, muted ? Math.max(1, w * 0.7) : w, dash);
      const [mx, my] = projectPoint(fault.trace[Math.floor(fault.trace.length / 2)] as Vec3);
      const r = 2.6 * pixelScale;
      ctx.fillStyle = `rgba(${KIND_RGB[fault.kind]}, ${alpha})`;
      ctx.strokeStyle = `rgba(${KIND_RGB[fault.kind]}, ${Math.min(1, alpha + 0.15)})`;
      ctx.lineWidth = lineWidth * 1.2;
      ctx.beginPath();
      ctx.arc(mx, my, r, 0, 2 * Math.PI);
      if (fault.active) ctx.fill();
      ctx.stroke();
    };
    const faultRank = (f: FaultSummary) => (hk != null && f.kind === hk ? 2 : 0) + (f.active ? 1 : 0);
    for (const f of [...faults].sort((a, b) => faultRank(a) - faultRank(b))) drawFault(f);

    // 4b. Selected plate's Euler pole + a motion arc whose ground length tracks plate speed.
    if (selectedPlate && selectedPlate.euler_pole && selectedPlate.outline.length > 0) {
      const ma = hk ? 0.3 : 0.95;
      const poleAxis = latLonToXyz(
        (selectedPlate.euler_pole.lat_deg * Math.PI) / 180,
        (selectedPlate.euler_pole.lon_deg * Math.PI) / 180,
      );
      // The centroid's small-circle path about the pole, swept over a fixed interval: Δθ ∝
      // speed_cm_per_yr, and the small-circle radius sin(φ) scales it back to a true ground
      // distance -- so the drawn arc length is proportional to the plate's surface speed here.
      const centroid = loopCentroid(selectedPlate.outline as Vec3[]);
      const D2R = Math.PI / 180;
      const sweep = Math.min(
        MOTION_ARC_MAX_DEG * D2R,
        Math.max(MOTION_ARC_MIN_DEG * D2R, selectedPlate.speed_cm_per_yr * MOTION_ARC_DEG_PER_CM_YR * D2R),
      );
      const STEPS = 48;
      const arcPts: Vec3[] = [];
      for (let i = 0; i <= STEPS; i++) {
        arcPts.push(matApply(rotationMatrix(poleAxis, (sweep * i) / STEPS), centroid) as Vec3);
      }
      ctx.strokeStyle = `rgba(${MOTION_RGB}, ${ma})`;
      ctx.lineWidth = lineWidth * 2.2;
      ctx.lineCap = "round";
      for (let i = 0; i + 1 < arcPts.length; i++) {
        const [[x1, y1], [x2, y2]] = projectEdge(arcPts[i], arcPts[i + 1], previewRotation, transform);
        if (Math.hypot(x2 - x1, y2 - y1) > seamJump) continue;
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
      }
      ctx.lineCap = "butt";
      // Arrowhead at the leading end (direction = right-hand rotation about the Euler pole).
      const [tailPx, tipPx] = projectEdge(arcPts[STEPS - 1], arcPts[STEPS], previewRotation, transform);
      if (Math.hypot(tipPx[0] - tailPx[0], tipPx[1] - tailPx[1]) < seamJump) {
        const ang = Math.atan2(tipPx[1] - tailPx[1], tipPx[0] - tailPx[0]);
        const ah = 9 * pixelScale;
        ctx.fillStyle = `rgba(${MOTION_RGB}, ${ma})`;
        ctx.beginPath();
        ctx.moveTo(tipPx[0], tipPx[1]);
        ctx.lineTo(tipPx[0] - ah * Math.cos(ang - 0.42), tipPx[1] - ah * Math.sin(ang - 0.42));
        ctx.lineTo(tipPx[0] - ah * Math.cos(ang + 0.42), tipPx[1] - ah * Math.sin(ang + 0.42));
        ctx.closePath();
        ctx.fill();
      }
      // Euler pole glyph: a ringed crosshair at the pole, a hollow ring at its antipode.
      const drawPole = (axis: Vec3, filled: boolean) => {
        const [px, py] = projectPoint(axis);
        const rr = 5 * pixelScale;
        ctx.strokeStyle = `rgba(${MOTION_RGB}, ${ma})`;
        ctx.fillStyle = `rgba(${MOTION_RGB}, ${ma})`;
        ctx.lineWidth = lineWidth * 1.4;
        ctx.beginPath();
        ctx.arc(px, py, rr, 0, 2 * Math.PI);
        ctx.stroke();
        if (filled) {
          ctx.beginPath();
          ctx.arc(px, py, rr * 0.42, 0, 2 * Math.PI);
          ctx.fill();
        }
        ctx.beginPath();
        ctx.moveTo(px - rr * 1.7, py);
        ctx.lineTo(px + rr * 1.7, py);
        ctx.moveTo(px, py - rr * 1.7);
        ctx.lineTo(px, py + rr * 1.7);
        ctx.stroke();
      };
      drawPole(poleAxis, true);
      drawPole([-poleAxis[0], -poleAxis[1], -poleAxis[2]], false);
    }

    if (!showQuakesVolcanoes) return;

    // 5. Earthquake epicentres -- filled dot + ring, sized by magnitude, fading with age.
    for (const q of earthquakes) {
      const recency = Math.max(0, 1 - q.age_myr / EARTHQUAKE_RETAIN_MYR);
      if (recency <= 0) continue;
      const [ex, ey] = projectPoint(q.epicenter as Vec3);
      const r = (1.5 + 1.6 * Math.max(0, q.magnitude - 4)) * pixelScale;
      ctx.beginPath();
      ctx.arc(ex, ey, r, 0, 2 * Math.PI);
      ctx.fillStyle = `rgba(${EARTHQUAKE_RGB}, ${(0.15 + 0.5 * recency) * activityDim})`;
      ctx.fill();
      ctx.lineWidth = lineWidth * 1.4;
      ctx.strokeStyle = `rgba(${EARTHQUAKE_RGB}, ${(0.4 + 0.6 * recency) * activityDim})`;
      ctx.stroke();
    }

    // 6. Volcanoes -- an upward triangle, filled when still active, hollow when dormant.
    for (const v of volcanoes) {
      const [vx, vy] = projectPoint(v.position as Vec3);
      const s = 3.4 * pixelScale;
      ctx.beginPath();
      ctx.moveTo(vx, vy - s);
      ctx.lineTo(vx - s * 0.9, vy + s * 0.7);
      ctx.lineTo(vx + s * 0.9, vy + s * 0.7);
      ctx.closePath();
      if (v.active) {
        ctx.fillStyle = `rgba(${VOLCANO_RGB}, ${0.9 * activityDim})`;
        ctx.fill();
      }
      ctx.lineWidth = lineWidth * 1.3;
      ctx.strokeStyle = `rgba(${VOLCANO_RGB}, ${(v.active ? 1 : 0.6) * activityDim})`;
      ctx.stroke();
    }
  };

  useEffect(() => {
    draw(rotation);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [
    plates, faults, faultSystems, earthquakes, volcanoes, coastlineSegments, showQuakesVolcanoes,
    highlightedFaultKind, selectedPlateId, projection, rotation, width, height,
  ]);

  const handleClick = (backingX: number, backingY: number) => {
    const transform = getRenderTransform(projection, width, height);
    const x = (backingX - transform.offsetX) / transform.scale;
    const y = -(backingY - transform.offsetY) / transform.scale;
    const latLon = unproject(projection, x, y);
    if (!latLon) return;
    const displayXyz = latLonToXyz(latLon[0], latLon[1]);
    const trueXyz = matApply(matTranspose(rotation), displayXyz);
    const [trueLat, trueLon] = xyzToLatLon(trueXyz);
    fetchPlateAt((trueLat * 180) / Math.PI, (trueLon * 180) / Math.PI)
      .then(({ plate_id }) => {
        onSelectPlate(plate_id);
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
    disabled: interactionDisabled,
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
