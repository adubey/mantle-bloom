import { useEffect, useRef, useState } from "react";
import type { GlobeBasis, MapPixels } from "./globeTexture";
import { loadWorldMapPixels, paintGlobe } from "./globeTexture";
import type { Mat3 } from "./rotation";
import { centerOfRotation, matApply, matTranspose, rotationForCenter, wrapPanLatLon } from "./rotation";

// The sidebar's globe selector (#497): a small globe of the world, centred on the map's view
// centre, that pans the map when dragged. Grabbing it moves the surface under the pointer, so
// dragging by the globe's radius turns it by one radian -- the same linear panning as the
// map's own long-press drag (see rotationDrag.ts). The map re-renders only when the drag ends;
// while dragging, the globe itself and the "Center" readout follow the pointer.

const SIZE_PX = 140;
const BACKING_SCALE = 2;
const RADIUS_PX = SIZE_PX / 2 - 3;
const KEYBOARD_STEP_RAD = (5 * Math.PI) / 180;

interface Props {
  worldRevision: number;
  rotation: Mat3;
  onRotationPreview: (latDeg: number, lonDeg: number) => void;
  onRotationCommitted: (rotation: Mat3) => void;
  disabled?: boolean;
}

// The display frame shows its own lat=0/lon=0 at the centre, east to the right, north up;
// the transpose carries those back to the true frame the texture is in.
function basisFor(rotation: Mat3): GlobeBasis {
  const inverse = matTranspose(rotation);
  return [matApply(inverse, [1, 0, 0]), matApply(inverse, [0, 1, 0]), matApply(inverse, [0, 0, 1])];
}

export default function GlobeSelector({ worldRevision, rotation, onRotationPreview, onRotationCommitted, disabled }: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [map, setMap] = useState<MapPixels | null>(null);
  const [preview, setPreview] = useState<Mat3 | null>(null);
  const drag = useRef<{ pointerId: number; startX: number; startY: number; base: [number, number] } | null>(null);
  const shown = preview ?? rotation;

  useEffect(() => {
    let cancelled = false;
    loadWorldMapPixels(worldRevision)
      .then((pixels) => {
        if (!cancelled) setMap(pixels);
      })
      .catch(() => {
        if (!cancelled) setMap(null);
      });
    return () => {
      cancelled = true;
    };
  }, [worldRevision]);

  useEffect(() => {
    const context = canvasRef.current?.getContext("2d");
    if (!context) return;
    const size = SIZE_PX * BACKING_SCALE;
    const center = size / 2;
    const radius = RADIUS_PX * BACKING_SCALE;
    context.clearRect(0, 0, size, size);
    context.fillStyle = "#122239";
    context.beginPath();
    context.arc(center, center, radius, 0, 2 * Math.PI);
    context.fill();
    if (map) {
      const pixels = context.createImageData(size, size);
      paintGlobe(map, basisFor(shown), pixels, center, center, radius);
      const layer = document.createElement("canvas");
      layer.width = size;
      layer.height = size;
      layer.getContext("2d")?.putImageData(pixels, 0, 0);
      context.drawImage(layer, 0, 0);
    }
    context.strokeStyle = "#8796ad";
    context.lineWidth = 1.5 * BACKING_SCALE;
    context.beginPath();
    context.arc(center, center, radius, 0, 2 * Math.PI);
    context.stroke();
    // The map's centre.
    context.strokeStyle = "rgba(255, 255, 255, 0.85)";
    context.lineWidth = BACKING_SCALE;
    context.beginPath();
    context.moveTo(center - 6 * BACKING_SCALE, center);
    context.lineTo(center + 6 * BACKING_SCALE, center);
    context.moveTo(center, center - 6 * BACKING_SCALE);
    context.lineTo(center, center + 6 * BACKING_SCALE);
    context.stroke();
  }, [map, shown]);

  const panned = (base: [number, number], dLat: number, dLon: number): Mat3 => {
    const [lat, lon] = wrapPanLatLon(base[0] + dLat, base[1] + dLon);
    onRotationPreview((lat * 180) / Math.PI, (lon * 180) / Math.PI);
    return rotationForCenter(lat, lon);
  };

  const onPointerDown = (event: React.PointerEvent<HTMLCanvasElement>) => {
    if (disabled || event.button !== 0) return;
    event.currentTarget.setPointerCapture(event.pointerId);
    drag.current = { pointerId: event.pointerId, startX: event.clientX, startY: event.clientY, base: centerOfRotation(rotation) };
  };

  const onPointerMove = (event: React.PointerEvent<HTMLCanvasElement>) => {
    const current = drag.current;
    if (!current || current.pointerId !== event.pointerId) return;
    // Dragging down brings the north into view, as when grabbing a real globe.
    const dLat = (event.clientY - current.startY) / RADIUS_PX;
    const dLon = -(event.clientX - current.startX) / RADIUS_PX;
    setPreview(panned(current.base, dLat, dLon));
  };

  const endDrag = (event: React.PointerEvent<HTMLCanvasElement>) => {
    const current = drag.current;
    if (!current || current.pointerId !== event.pointerId) return;
    drag.current = null;
    const moved = event.clientX !== current.startX || event.clientY !== current.startY;
    setPreview(null);
    if (moved && event.type === "pointerup") {
      onRotationCommitted(panned(current.base, (event.clientY - current.startY) / RADIUS_PX, -(event.clientX - current.startX) / RADIUS_PX));
    } else {
      const [lat, lon] = current.base;
      onRotationPreview((lat * 180) / Math.PI, (lon * 180) / Math.PI);
    }
  };

  const onKeyDown = (event: React.KeyboardEvent<HTMLCanvasElement>) => {
    const step: Record<string, [number, number]> = {
      ArrowUp: [KEYBOARD_STEP_RAD, 0],
      ArrowDown: [-KEYBOARD_STEP_RAD, 0],
      ArrowLeft: [0, -KEYBOARD_STEP_RAD],
      ArrowRight: [0, KEYBOARD_STEP_RAD],
    };
    const delta = step[event.key];
    if (disabled || !delta) return;
    event.preventDefault();
    onRotationCommitted(panned(centerOfRotation(rotation), delta[0], delta[1]));
  };

  return (
    <canvas
      ref={canvasRef}
      width={SIZE_PX * BACKING_SCALE}
      height={SIZE_PX * BACKING_SCALE}
      tabIndex={disabled ? -1 : 0}
      role="application"
      className="globe-control"
      aria-label="Globe selector: drag, or use the arrow keys, to move the map's centre"
      title="Drag the globe to move the map"
      onPointerDown={onPointerDown}
      onPointerMove={onPointerMove}
      onPointerUp={endDrag}
      onPointerCancel={endDrag}
      onKeyDown={onKeyDown}
      style={{
        display: "block",
        width: SIZE_PX,
        height: SIZE_PX,
        margin: "6px auto 0",
        cursor: disabled ? "default" : "grab",
        touchAction: "none",
        opacity: disabled ? 0.5 : 1,
      }}
    />
  );
}
