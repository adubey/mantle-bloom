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

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;

    ctx.fillStyle = BACKGROUND;
    ctx.fillRect(0, 0, width, height);
    if (!imageBase64) return;

    // Image decoding is async -- if a newer render arrives (or this component unmounts)
    // before this one finishes loading, drop it rather than painting a stale frame.
    let cancelled = false;
    const img = new Image();
    img.onload = () => {
      if (cancelled) return;
      ctx.drawImage(img, 0, 0, width, height);
    };
    img.src = `data:image/png;base64,${imageBase64}`;
    return () => {
      cancelled = true;
    };
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
