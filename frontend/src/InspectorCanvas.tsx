import { useEffect, useRef } from "react";
import type { KeyboardEvent, RefObject } from "react";
import type { Projection } from "./api";
import type { Mat3 } from "./rotation";
import { useRotationDrag } from "./rotationDrag";
import { trueLatLonDegAt } from "./screenProjector";

// The props every client-drawn inspector view takes from App.
export interface InspectorViewProps {
  width: number;
  height: number;
  displayWidth: number;
  displayHeight: number;
  projection: Projection;
  // The currently *committed* view rotation -- the same rotation state MapCanvas and every
  // inspector share, so switching between views preserves orientation.
  rotation: Mat3;
  onRotationPreview: (latDeg: number, lonDeg: number) => void;
  onRotationCommitted: (rotation: Mat3) => void;
  // Inert rotate-drag + click while a background animation holds the world lock (see App.tsx).
  interactionDisabled?: boolean;
}

interface Props extends InspectorViewProps {
  canvasRef: RefObject<HTMLCanvasElement | null>;
  // Redraws the view at a (possibly in-progress) rotation -- see useRotationDrag's onFrame.
  onFrame: (previewRotation: Mat3) => void;
  // A plain click at a true lat/lon on the sphere. Resolves once the selection is made, and
  // the view then takes keyboard focus so Tab cycles from there.
  onPick: (latDeg: number, lonDeg: number) => Promise<void>;
  // Tab (forward) or Shift+Tab (backwards).
  onCycle: (backwards: boolean) => void;
}

// The interactive shell shared by PlateInspector, PlatesAndFaults, LakeInspector and
// RiverInspector: the focusable canvas, the long-press-drag rotate gesture (see
// rotationDrag.ts), click-to-select and Tab/Shift+Tab cycling. The view itself draws into
// `canvasRef`.
export default function InspectorCanvas({
  canvasRef, width, height, displayWidth, displayHeight, projection, rotation,
  onFrame, onPick, onCycle, onRotationPreview, onRotationCommitted, interactionDisabled,
}: Props) {
  const containerRef = useRef<HTMLDivElement>(null);

  const handleClick = (backingX: number, backingY: number) => {
    const latLon = trueLatLonDegAt(projection, rotation, width, height, backingX, backingY);
    if (!latLon) return; // clicked outside the sphere's projected silhouette
    onPick(latLon[0], latLon[1])
      .then(() => containerRef.current?.focus())
      .catch(() => {
        // A click that races a generate/step (no world, or a stale list) -- ignored, since a
        // purely exploratory gesture shouldn't surface a toast.
      });
  };

  useRotationDrag({
    elementRef: canvasRef,
    width, height, displayWidth, displayHeight, projection, rotation,
    onFrame,
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
    onCycle(e.shiftKey);
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
