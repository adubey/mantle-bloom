import { useEffect, useRef } from "react";

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
}

const BACKGROUND = "#0b1020";

// All map drawing (fill colors, plate boundaries/poles/rotation arcs, per-plate node
// dots) now happens server-side per requested view -- see backend app/render_image.py.
// This component's only job is to decode the returned PNG and paint it onto the canvas.
export default function MapCanvas({ imageBase64, width, height, displayWidth, displayHeight }: Props) {
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

  return (
    <canvas
      ref={canvasRef}
      width={width}
      height={height}
      style={{ borderRadius: 8, width: displayWidth, height: displayHeight }}
    />
  );
}
