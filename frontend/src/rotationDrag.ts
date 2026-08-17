import { useEffect, useRef } from "react";
import type { RefObject } from "react";
import type { Mat3 } from "./rotation";
import { IDENTITY_ROTATION, matApply, matMultiply, rotationBetween, screenToArcballVector, xyzToLatLon } from "./rotation";

// The long-press-then-arcball-drag "rotate the planet" gesture (see
// docs/simulation-model.md#rotating-the-view), extracted out of MapCanvas.tsx so
// PlateInspector.tsx can reuse the exact same interaction instead of duplicating this
// state machine -- both need "drag to rotate the view," they just draw completely different
// content while doing it (a PNG + graticule wireframe vs. translucent plate ellipses).

const LONG_PRESS_MS = 350;
// In *display* (CSS) pixels -- movement past this before the long-press timer fires cancels
// it, so an ordinary click/drag-to-select-text gesture elsewhere on the page never
// accidentally starts a rotation.
const MOVE_CANCEL_PX = 6;
const ARCBALL_RADIUS_FRACTION = 0.45; // of min(width, height), in backing-store pixels

export interface UseRotationDragOptions {
  elementRef: RefObject<HTMLElement | null>;
  width: number;
  height: number;
  displayWidth: number;
  displayHeight: number;
  // The currently *committed* view rotation -- each new drag starts from this and composes
  // its own incremental rotation on top.
  rotation: Mat3;
  // Called once when a drag begins (with the committed rotation, before any movement) and
  // again on every mouse move while dragging, so the caller can redraw its own content at
  // the in-progress rotation.
  onFrame: (previewRotation: Mat3) => void;
  onRotationCommitted: (rotation: Mat3) => void;
  // The map center's live lat/lon under the in-progress rotation -- optional since not every
  // caller needs a live readout.
  onRotationPreview?: (latDeg: number, lonDeg: number) => void;
  // An ordinary press+release that never turned into a drag (mouseup arrived before the
  // long-press timer armed, with movement under MOVE_CANCEL_PX) -- reported in *backing
  // store* pixel coordinates. Routed through this same state machine rather than a separate
  // click listener specifically so a completed drag's terminating mouseup can never also be
  // mistaken for a click on whatever happens to be under the release point.
  onClick?: (backingX: number, backingY: number) => void;
}

export function useRotationDrag(options: UseRotationDragOptions): void {
  const { elementRef, width, height, displayWidth, displayHeight, rotation, onFrame, onRotationCommitted, onRotationPreview, onClick } =
    options;

  // Interaction state lives in refs, not React state -- none of it should trigger a
  // re-render; content is painted directly, imperatively, by the caller's onFrame.
  const longPressTimer = useRef<ReturnType<typeof setTimeout> | null>(null);
  const dragStartClientPos = useRef<{ clientX: number; clientY: number } | null>(null);
  const dragStartArcballVec = useRef<[number, number, number] | null>(null);
  const dragBaseRotation = useRef<Mat3>(IDENTITY_ROTATION);
  const isRotating = useRef(false);

  // Callbacks are read through refs updated on every render (not via the effect's own
  // dependency array) so they can close over whatever data the caller's latest render has
  // (e.g. PlateInspector's current `plates`/`selectedPlateId`) without going stale, *and*
  // without tearing down and rebuilding the drag's window-level mousemove/mouseup listeners
  // on every render -- critical here specifically because onRotationPreview typically
  // triggers a parent state update (a live lat/lon readout) on every single mousemove, which
  // would otherwise mean a teardown+rebuild mid-drag on nearly every frame.
  const onFrameRef = useRef(onFrame);
  onFrameRef.current = onFrame;
  const onRotationCommittedRef = useRef(onRotationCommitted);
  onRotationCommittedRef.current = onRotationCommitted;
  const onRotationPreviewRef = useRef(onRotationPreview);
  onRotationPreviewRef.current = onRotationPreview;
  const onClickRef = useRef(onClick);
  onClickRef.current = onClick;

  useEffect(() => {
    const element = elementRef.current;
    if (!element) return;

    const toBackingPixels = (clientX: number, clientY: number) => {
      const rect = element.getBoundingClientRect();
      const displayX = clientX - rect.left;
      const displayY = clientY - rect.top;
      return { x: displayX * (width / displayWidth), y: displayY * (height / displayHeight) };
    };

    const previewRotationFor = (backingX: number, backingY: number): Mat3 => {
      const start = dragStartArcballVec.current;
      if (!start) return dragBaseRotation.current;
      const radius = ARCBALL_RADIUS_FRACTION * Math.min(width, height);
      const current = screenToArcballVector(backingX, backingY, width / 2, height / 2, radius);
      const delta = rotationBetween(start, current);
      return matMultiply(delta, dragBaseRotation.current);
    };

    const reportPreview = (previewRotation: Mat3) => {
      if (!onRotationPreviewRef.current) return;
      const [lat, lon] = xyzToLatLon(matApply(previewRotation, [1, 0, 0]));
      onRotationPreviewRef.current((lat * 180) / Math.PI, (lon * 180) / Math.PI);
    };

    const onWindowMouseMove = (e: MouseEvent) => {
      if (!isRotating.current) return;
      const { x, y } = toBackingPixels(e.clientX, e.clientY);
      const preview = previewRotationFor(x, y);
      onFrameRef.current(preview);
      reportPreview(preview);
    };

    const onWindowMouseUp = (e: MouseEvent) => {
      window.removeEventListener("mousemove", onWindowMouseMove);
      window.removeEventListener("mouseup", onWindowMouseUp);
      if (!isRotating.current) return;
      isRotating.current = false;
      const { x, y } = toBackingPixels(e.clientX, e.clientY);
      const finalRotation = previewRotationFor(x, y);
      dragStartArcballVec.current = null;
      onRotationCommittedRef.current(finalRotation);
    };

    const beginRotating = () => {
      isRotating.current = true;
      dragBaseRotation.current = rotation;
      const pos = dragStartClientPos.current;
      if (pos) {
        const { x, y } = toBackingPixels(pos.clientX, pos.clientY);
        const radius = ARCBALL_RADIUS_FRACTION * Math.min(width, height);
        dragStartArcballVec.current = screenToArcballVector(x, y, width / 2, height / 2, radius);
      }
      const preview = dragBaseRotation.current;
      onFrameRef.current(preview);
      reportPreview(preview);
      window.addEventListener("mousemove", onWindowMouseMove);
      window.addEventListener("mouseup", onWindowMouseUp);
    };

    const onMouseDown = (e: MouseEvent) => {
      if (e.button !== 0) return;
      dragStartClientPos.current = { clientX: e.clientX, clientY: e.clientY };
      longPressTimer.current = setTimeout(beginRotating, LONG_PRESS_MS);

      const cancelIfMoved = (moveEvent: MouseEvent) => {
        const dx = moveEvent.clientX - e.clientX;
        const dy = moveEvent.clientY - e.clientY;
        if (Math.hypot(dx, dy) > MOVE_CANCEL_PX && longPressTimer.current) {
          clearTimeout(longPressTimer.current);
          longPressTimer.current = null;
          window.removeEventListener("mousemove", cancelIfMoved);
          window.removeEventListener("mouseup", cancelOnEarlyRelease);
        }
      };
      const cancelOnEarlyRelease = (upEvent: MouseEvent) => {
        // Only reachable when no movement past MOVE_CANCEL_PX has happened yet (otherwise
        // cancelIfMoved already unregistered this listener) -- i.e. this genuinely is a
        // plain click, not a drag that happened to end quickly.
        if (longPressTimer.current) {
          clearTimeout(longPressTimer.current);
          longPressTimer.current = null;
        }
        window.removeEventListener("mousemove", cancelIfMoved);
        window.removeEventListener("mouseup", cancelOnEarlyRelease);
        if (onClickRef.current) {
          const { x, y } = toBackingPixels(upEvent.clientX, upEvent.clientY);
          onClickRef.current(x, y);
        }
      };
      window.addEventListener("mousemove", cancelIfMoved);
      window.addEventListener("mouseup", cancelOnEarlyRelease);
    };

    element.addEventListener("mousedown", onMouseDown);
    return () => {
      element.removeEventListener("mousedown", onMouseDown);
      window.removeEventListener("mousemove", onWindowMouseMove);
      window.removeEventListener("mouseup", onWindowMouseUp);
      if (longPressTimer.current) clearTimeout(longPressTimer.current);
    };
    // onFrame/onRotationCommitted/onRotationPreview/onClick are deliberately excluded --
    // read through the *Ref indirection above instead, so a new callback identity on every
    // render (very likely: onRotationPreview typically drives a parent state update on every
    // mousemove) never tears down and rebuilds the drag's window-level listeners mid-drag.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [elementRef, width, height, displayWidth, displayHeight, rotation]);
}
