import { useEffect, useRef } from "react";
import type { Projection } from "./api";
import type { Mat3 } from "./rotation";
import { getGraticule, getRenderTransform, matApply, project, toPixels } from "./rotation";
import { useRotationDrag } from "./rotationDrag";

interface Props {
  imageBase64: string | null;
  // width/height: the canvas's backing-store resolution, matching the image's own pixel
  // dimensions (see api.ts's renderWorld). displayWidth/displayHeight: its CSS size --
  // requesting a higher backing-store resolution than the CSS size is what makes the map
  // render sharper (retina-style) without taking up more room on the page.
  width: number;
  height: number;
  displayWidth: number;
  displayHeight: number;
  projection: Projection;
  // The currently *committed* view rotation (see rotation.ts) -- each new drag starts from
  // this and composes its own incremental rotation on top; a completed drag is reported back
  // via onRotationCommitted, which the parent stores and re-renders with.
  rotation: Mat3;
  // Fired continuously while the user drags (long-press then move) with the map center's
  // live lat/lon under the *in-progress* rotation -- needed because the real legend is baked
  // server-side into the PNG and can't update mid-drag (see docs/simulation-model.md#rotating-the-view).
  onRotationPreview: (latDeg: number, lonDeg: number) => void;
  onRotationCommitted: (rotation: Mat3) => void;
}

const BACKGROUND = "#0b1020";
const GRATICULE_COLOR = "rgba(235, 238, 245, 0.85)";
const GRATICULE_LINE_WIDTH = 1;
// A graticule line segment longer than this multiple of the line's own median segment length
// is skipped -- the same "don't draw across a projection discontinuity" technique
// render_image.py's _stroke_robust_loop uses for plate boundaries.
const SEGMENT_BREAK_FACTOR = 6;

// All map drawing (fill colors, plate boundaries/poles/rotation arcs, per-plate node
// dots) happens server-side per requested view -- see backend app/render_image.py. This
// component decodes the returned PNG and paints it onto the canvas, and additionally handles
// the long-press-and-drag "rotate the planet" gesture (see rotationDrag.ts): while dragging,
// it draws a cheap wireframe graticule preview client-side (see rotation.ts) instead of
// re-requesting the real, much more expensive, detailed render on every mouse move.
export default function MapCanvas({
  imageBase64, width, height, displayWidth, displayHeight, projection, rotation,
  onRotationPreview, onRotationCommitted,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  // One Image element, reused for the component's whole lifetime rather than a fresh
  // `new Image()` per update -- lazily created here (the standard ref-during-render pattern
  // for a one-time, non-render-affecting object) since it doesn't depend on any prop.
  const imgRef = useRef<HTMLImageElement | null>(null);
  if (imgRef.current === null) {
    imgRef.current = new Image();
  }

  useEffect(() => {
    const canvas = canvasRef.current;
    const img = imgRef.current;
    if (!canvas || !img) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    if (!imageBase64) {
      ctx.fillStyle = BACKGROUND;
      ctx.fillRect(0, 0, width, height);
      return;
    }

    // The flicker fix: don't clear the canvas up front. The previous frame stays on screen,
    // fully intact, right up until the new one has actually finished decoding -- then
    // drawImage swaps it in with a single paint covering the whole canvas, so there's never
    // a moment where the canvas shows neither frame. Reassigning .src on the same element
    // (rather than creating a new Image) also means a still-decoding previous frame is
    // simply superseded -- the browser never fires onload for an aborted load, so stepping
    // faster than a decode completes can't paint a stale frame after a newer one.
    img.onload = () => {
      ctx.drawImage(img, 0, 0, width, height);
    };
    img.src = `data:image/png;base64,${imageBase64}`;
  }, [imageBase64, width, height]);

  const drawGraticule = (previewRotation: Mat3) => {
    const canvas = canvasRef.current;
    const ctx = canvas?.getContext("2d");
    const img = imgRef.current;
    if (!canvas || !ctx || !img) return;
    if (imageBase64) ctx.drawImage(img, 0, 0, width, height);
    else {
      ctx.fillStyle = BACKGROUND;
      ctx.fillRect(0, 0, width, height);
    }

    const transform = getRenderTransform(projection, width, height);
    ctx.strokeStyle = GRATICULE_COLOR;
    ctx.lineWidth = GRATICULE_LINE_WIDTH * (width / 1100);

    for (const line of getGraticule()) {
      const pixels = line.points.map((p) => {
        const rotated = matApply(previewRotation, p);
        const lon = Math.atan2(rotated[1], rotated[0]);
        const lat = Math.asin(Math.min(1, Math.max(-1, rotated[2])));
        const [x, y] = project(projection, lat, lon);
        return toPixels(transform, x, y);
      });
      const segLengths: number[] = [];
      for (let i = 0; i < pixels.length - 1; i++) {
        segLengths.push(Math.hypot(pixels[i + 1][0] - pixels[i][0], pixels[i + 1][1] - pixels[i][1]));
      }
      const sorted = [...segLengths].sort((a, b) => a - b);
      const median = sorted[Math.floor(sorted.length / 2)] || 0;
      const breakThreshold = Math.max(median * SEGMENT_BREAK_FACTOR, 20 * (width / 1100));

      ctx.beginPath();
      let penDown = false;
      for (let i = 0; i < pixels.length; i++) {
        if (i > 0 && segLengths[i - 1] > breakThreshold) penDown = false;
        if (!penDown) {
          ctx.moveTo(pixels[i][0], pixels[i][1]);
          penDown = true;
        } else {
          ctx.lineTo(pixels[i][0], pixels[i][1]);
        }
      }
      ctx.stroke();
    }
  };

  useRotationDrag({
    elementRef: canvasRef,
    width, height, displayWidth, displayHeight, rotation,
    onFrame: drawGraticule,
    onRotationPreview,
    onRotationCommitted,
  });

  return (
    <canvas
      ref={canvasRef}
      width={width}
      height={height}
      style={{ borderRadius: 8, width: displayWidth, height: displayHeight }}
    />
  );
}
