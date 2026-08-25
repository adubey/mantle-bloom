import { useCallback, useEffect, useRef } from "react";
import type { Projection } from "./api";
import type { HighlightTarget } from "./legendData";
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
  // Legend-click-to-highlight (Biome and Combined views -- see Legend.tsx/App.tsx): when set,
  // every decoded pixel whose nearest color in `palette` isn't `selected` is faded toward gray
  // so the selected swatch's cells visibly pop against the rest of the map (see applyHighlight
  // below). `null`/omitted paints the decoded frame as-is, same as before this feature existed.
  highlightTarget?: HighlightTarget | null;
}

const BACKGROUND = "#0b1020";
const GRATICULE_COLOR = "rgba(235, 238, 245, 0.85)";
const GRATICULE_LINE_WIDTH = 1;
// A graticule line segment longer than this multiple of the line's own median segment length
// is skipped -- the same "don't draw across a projection discontinuity" technique
// render_image.py's _stroke_robust_loop uses for plate boundaries.
const SEGMENT_BREAK_FACTOR = 6;
// How much a non-matching pixel's grayscale brightness is scaled down by when a legend
// highlight is active -- low enough that landmass/ocean shapes are still legible as context,
// enough contrast that the exact-match biome cells still read as clearly "lit up" next to it.
const HIGHLIGHT_DIM_FACTOR = 0.35;

// The Biome view (see backend app/render_image.py's _render_biome_view) draws every pixel as
// exactly one of biomes.BIOME_COLORS' fixed palette -- no coastline/graticule overlay on top
// -- so an exact RGB match (tolerance 0) against the clicked legend swatch's color is enough
// to pick out that biome's cells, entirely client-side, with no new server render mode
// needed. The Combined view instead matches within `tolerance` against several candidate
// colors (see legendData.ts's highlightTargetFor) since a biome's flat color there gets
// shaded by elevation (see app/biomes.py's biome_relative_shade_factor) and, at real peaks,
// further blended toward the elevation gradient on top of that.
//
// Classification is nearest-neighbor across the *entire* palette, not just a within-tolerance
// check against the selected label's own colors: two land biomes' shaded variants can land
// close together in RGB space (e.g. two similarly-hued forests at different elevations), and
// checking the selected label in isolation would highlight both whenever a pixel merely came
// within `tolerance` of it, even if some other label was actually the closer match. Picking
// the closest label first, and only then comparing it to the selection, keeps one biome's
// highlight from bleeding into another's cells. `tolerance` still bounds how far the *closest*
// match can be before a pixel counts as belonging to no known label at all (plain ocean,
// rivers, coastline, graticule, etc.), so those never get swept into whichever label happens
// to be nearest.
function classifyPixel(r: number, g: number, b: number, palette: HighlightTarget["palette"], tolerance2: number): string | null {
  let bestLabel: string | null = null;
  let bestDist2 = Infinity;
  for (const entry of palette) {
    for (const [cr, cg, cb] of entry.colors) {
      const dr = r - cr, dg = g - cg, db = b - cb;
      const dist2 = dr * dr + dg * dg + db * db;
      if (dist2 < bestDist2) {
        bestDist2 = dist2;
        bestLabel = entry.label;
      }
    }
  }
  return bestDist2 <= tolerance2 ? bestLabel : null;
}

function applyHighlight(ctx: CanvasRenderingContext2D, width: number, height: number, target: HighlightTarget): void {
  const { selected, palette, tolerance } = target;
  const tolerance2 = tolerance * tolerance;
  const imageData = ctx.getImageData(0, 0, width, height);
  const data = imageData.data;
  for (let i = 0; i < data.length; i += 4) {
    const r = data[i], g = data[i + 1], b = data[i + 2];
    if (classifyPixel(r, g, b, palette, tolerance2) === selected) continue;
    const gray = (0.3 * r + 0.59 * g + 0.11 * b) * HIGHLIGHT_DIM_FACTOR;
    data[i] = gray;
    data[i + 1] = gray;
    data[i + 2] = gray;
  }
  ctx.putImageData(imageData, 0, 0);
}

// All map drawing (fill colors, plate boundaries/poles/rotation arcs, per-plate node
// dots) happens server-side per requested view -- see backend app/render_image.py. This
// component decodes the returned PNG and paints it onto the canvas, and additionally handles
// the long-press-and-drag "rotate the planet" gesture (see rotationDrag.ts): while dragging,
// it draws a cheap wireframe graticule preview client-side (see rotation.ts) instead of
// re-requesting the real, much more expensive, detailed render on every mouse move.
export default function MapCanvas({
  imageBase64, width, height, displayWidth, displayHeight, projection, rotation,
  onRotationPreview, onRotationCommitted, highlightTarget,
}: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  // One Image element, reused for the component's whole lifetime rather than a fresh
  // `new Image()` per update -- lazily created here (the standard ref-during-render pattern
  // for a one-time, non-render-affecting object) since it doesn't depend on any prop.
  const imgRef = useRef<HTMLImageElement | null>(null);
  if (imgRef.current === null) {
    imgRef.current = new Image();
  }
  // Read from the highlight-toggle effect below without also making the (decode-driven) base
  // paint effect re-run just because the highlight selection changed.
  const highlightTargetRef = useRef(highlightTarget);
  highlightTargetRef.current = highlightTarget;

  // Draws the already-decoded frame plus, if a legend highlight is active, the filter on top
  // of it -- shared by both the initial decode (below) and the highlight-toggle effect, so
  // toggling a legend swatch doesn't need a fresh server render to update the map. Reads
  // highlightTarget via a ref (not a direct closure) so its identity only changes with
  // width/height, not with the highlight selection -- see that effect's own comment for why.
  const paintDecodedFrame = useCallback(() => {
    const canvas = canvasRef.current;
    const img = imgRef.current;
    const ctx = canvas?.getContext("2d");
    if (!canvas || !img || !ctx) return;
    ctx.drawImage(img, 0, 0, width, height);
    if (highlightTargetRef.current) applyHighlight(ctx, width, height, highlightTargetRef.current);
  }, [width, height]);

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
    img.onload = () => paintDecodedFrame();
    img.src = `data:image/png;base64,${imageBase64}`;
  }, [imageBase64, width, height, paintDecodedFrame]);

  // Re-applies (or clears) the highlight filter the instant the legend selection changes,
  // without waiting for a fresh render -- the image element already holds the fully decoded
  // frame at this point, so it can be redrawn synchronously. Deliberately keyed on
  // highlightTarget alone (imageBase64/paintDecodedFrame omitted from deps on purpose): this
  // effect exists to react to the highlight selection specifically, not to re-run
  // redundantly, one render after the effect above, on every image change too.
  useEffect(() => {
    if (imageBase64 && imgRef.current?.complete) paintDecodedFrame();
    // oxlint-disable-next-line react-hooks/exhaustive-deps
  }, [highlightTarget]);

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
    width, height, displayWidth, displayHeight, projection, rotation,
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
