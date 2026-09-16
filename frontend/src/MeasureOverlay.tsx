import { useEffect, useRef, useState } from "react";
import type { Projection } from "./api";
import { backingPixelsToDisplayLatLon, greatCircleDistanceKm } from "./rotation";

// Shift+click to place the first point, then a second click (Shift no longer needed) to place
// the second and measure the great-circle distance between them (issue #158). Shift+click
// specifically (rather than a plain click) is what arms this, since a plain long-press-and-drag
// on the bare map already means "rotate the globe" (see rotationDrag.ts) -- this overlay stays
// pointer-events:none (letting mouse events fall through to the normal rotate/probe handlers
// underneath) except while Shift is held or a measurement is still pending its second click.

interface Props {
  projection: Projection;
  width: number; // backing-store px
  height: number;
  displayWidth: number; // CSS px
  displayHeight: number;
  disabled?: boolean;
}

interface MeasureState {
  anchorDisplay: { x: number; y: number };
  anchorLatLon: [number, number];
  currentDisplay: { x: number; y: number };
  distanceKm: number;
  active: boolean; // still waiting for the second click vs. finalized and dismissable
}

export default function MeasureOverlay({ projection, width, height, displayWidth, displayHeight, disabled }: Props) {
  const [shiftHeld, setShiftHeld] = useState(false);
  const [measure, setMeasure] = useState<MeasureState | null>(null);
  const overlayRef = useRef<HTMLDivElement | null>(null);

  useEffect(() => {
    const onKeyDown = (e: KeyboardEvent) => {
      if (e.key === "Shift") setShiftHeld(true);
      if (e.key === "Escape") setMeasure(null);
    };
    const onKeyUp = (e: KeyboardEvent) => {
      if (e.key === "Shift") setShiftHeld(false);
    };
    // Window can lose focus (alt-tab, devtools) while Shift is physically still down but no
    // keyup ever arrives -- reset on blur so the overlay doesn't get stuck capturing clicks.
    const onBlur = () => setShiftHeld(false);
    window.addEventListener("keydown", onKeyDown);
    window.addEventListener("keyup", onKeyUp);
    window.addEventListener("blur", onBlur);
    return () => {
      window.removeEventListener("keydown", onKeyDown);
      window.removeEventListener("keyup", onKeyUp);
      window.removeEventListener("blur", onBlur);
    };
  }, []);

  const toDisplayLatLon = (clientX: number, clientY: number) => {
    const rect = overlayRef.current?.getBoundingClientRect();
    const displayX = clientX - (rect?.left ?? 0);
    const displayY = clientY - (rect?.top ?? 0);
    const backingX = displayX * (width / displayWidth);
    const backingY = displayY * (height / displayHeight);
    return { display: { x: displayX, y: displayY }, latLon: backingPixelsToDisplayLatLon(projection, width, height, backingX, backingY) };
  };

  // While a measurement is pending its second click, the line follows the free-floating mouse
  // (no button held, unlike a drag) so the distance updates live before it's placed.
  useEffect(() => {
    if (!measure?.active) return;
    const onMove = (e: MouseEvent) => {
      const { display, latLon } = toDisplayLatLon(e.clientX, e.clientY);
      setMeasure((prev) => {
        if (!prev || !prev.active) return prev;
        const distanceKm = latLon ? greatCircleDistanceKm(prev.anchorLatLon[0], prev.anchorLatLon[1], latLon[0], latLon[1]) : prev.distanceKm;
        return { ...prev, currentDisplay: display, distanceKm };
      });
    };
    window.addEventListener("mousemove", onMove);
    return () => window.removeEventListener("mousemove", onMove);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [measure?.active, projection, width, height, displayWidth, displayHeight]);

  if (disabled) return null;

  const onClick = (e: React.MouseEvent) => {
    if (e.button !== 0) return;
    const { display, latLon } = toDisplayLatLon(e.clientX, e.clientY);

    if (measure?.active) {
      // Second click: place the endpoint and finalize -- doesn't need Shift.
      if (!latLon) return; // clicked off the globe -- keep waiting for a valid endpoint
      const distanceKm = greatCircleDistanceKm(measure.anchorLatLon[0], measure.anchorLatLon[1], latLon[0], latLon[1]);
      setMeasure({ ...measure, currentDisplay: display, distanceKm, active: false });
      return;
    }

    // First click: only arms a new measurement while Shift is held.
    if (!e.shiftKey || !latLon) return;
    setMeasure({ anchorDisplay: display, anchorLatLon: latLon, currentDisplay: display, distanceKm: 0, active: true });
  };

  const midX = measure ? (measure.anchorDisplay.x + measure.currentDisplay.x) / 2 : 0;
  const midY = measure ? (measure.anchorDisplay.y + measure.currentDisplay.y) / 2 : 0;

  return (
    <div
      ref={overlayRef}
      onClick={onClick}
      style={{
        position: "absolute",
        left: 0,
        top: 0,
        width: displayWidth,
        height: displayHeight,
        cursor: shiftHeld || measure?.active ? "crosshair" : "default",
        // Only captures clicks while arming/mid-measurement -- once finalized (or never
        // started), this must go back to letting every mouse event fall through to the normal
        // rotate/probe handlers on the map underneath, or a finished-but-not-yet-dismissed
        // measurement would permanently block rotation. The result label's own close button
        // separately re-enables pointer-events on just itself (a child CAN opt back into
        // events even while its ancestor is none), so it stays clickable to dismiss without
        // Shift.
        pointerEvents: shiftHeld || measure?.active ? "auto" : "none",
      }}
    >
      {measure && (
        <svg width={displayWidth} height={displayHeight} style={{ position: "absolute", left: 0, top: 0, pointerEvents: "none" }}>
          <line
            x1={measure.anchorDisplay.x}
            y1={measure.anchorDisplay.y}
            x2={measure.currentDisplay.x}
            y2={measure.currentDisplay.y}
            stroke="#ffcf5c"
            strokeWidth={1.5}
            strokeDasharray="4 3"
          />
          <circle cx={measure.anchorDisplay.x} cy={measure.anchorDisplay.y} r={3} fill="#ffcf5c" />
          <circle cx={measure.currentDisplay.x} cy={measure.currentDisplay.y} r={3} fill="#ffcf5c" />
        </svg>
      )}
      {measure && (
        <div
          style={{
            position: "absolute",
            left: Math.max(4, Math.min(midX + 8, displayWidth - 90)),
            top: Math.max(4, Math.min(midY - 22, displayHeight - 30)),
            background: "#151a2e",
            border: "1px solid #4b5060",
            borderRadius: 4,
            padding: "3px 6px",
            fontSize: 11,
            color: "#dee2eb",
            whiteSpace: "nowrap",
            pointerEvents: measure.active ? "none" : "auto",
            display: "flex",
            alignItems: "center",
            gap: 6,
          }}
        >
          <span>{Math.round(measure.distanceKm).toLocaleString()} km</span>
          {!measure.active && (
            <button
              type="button"
              title="Clear"
              onClick={(e) => {
                // Stop this from also bubbling to the overlay's own onClick, which -- reading
                // stale state from the same event -- could otherwise misfire as a fresh
                // measurement's first click (e.g. if Shift happened to be held while clicking ×).
                e.stopPropagation();
                setMeasure(null);
              }}
              style={{
                border: "none", background: "transparent", color: "#999", cursor: "pointer",
                fontSize: 12, lineHeight: 1, padding: 0,
              }}
            >
              ×
            </button>
          )}
        </div>
      )}
    </div>
  );
}
