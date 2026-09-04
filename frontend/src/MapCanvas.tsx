import { useCallback, useEffect, useRef } from "react";
import type { Projection } from "./api";
import type { HighlightTarget } from "./legendData";
import type { Mat3 } from "./rotation";
import {
  getGraticule, getRenderTransform, latLonToXyz, matApply, matTranspose, project, toPixels, unproject, xyzToLatLon,
} from "./rotation";
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
  // every decoded pixel that isn't classified as `selected` is faded toward gray so the
  // selected swatch's cells visibly pop against the rest of the map (see applyHighlight
  // below). `null`/omitted paints the decoded frame as-is, same as before this feature existed.
  highlightTarget?: HighlightTarget | null;
  // Elevation & Biome / Elevation / Biome views only (App.tsx passes it just for those): a
  // plain click -- not a drag -- reports the clicked point's true-frame lat/lon plus its
  // position in display (CSS) pixels, for the click-to-inspect popup App renders over the
  // map. A click outside the projected globe silhouette reports null (dismiss the popup).
  // Omitted on every other view, which leaves a click doing nothing, exactly as before.
  onProbe?: (probe: { displayX: number; displayY: number; latDeg: number; lonDeg: number } | null) => void;
  // Biome and Combined mode both encode a per-pixel *dominant* biome id in the PNG's alpha
  // channel (see backend render_image.py's COMBINED_LAKE_ID_CODE comment / legendData.ts) --
  // both views now blend a boundary cell's RGB toward its runner-up class, so RGB alone no
  // longer identifies a pixel's class. When true, every painted frame gets a full-canvas pass
  // that reads those ids and then resets alpha to fully opaque -- the 237..255 alpha spread is
  // a data channel, not meant to actually composite.
  alphaEncodedIds?: boolean;
  // While a background animation holds the world lock (see App.tsx), the rotate-drag + probe
  // click are inert -- a rotation/render there would only 503.
  interactionDisabled?: boolean;
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

// Neither the Biome nor the Combined view matches on color: both blend a boundary cell's RGB
// toward its runner-up class (see backend app/render_image.py's _biome_blend_rgb), so a pixel's
// RGB alone no longer reliably identifies its class. Each carries its per-pixel *dominant*
// (>50% share) biome id in the alpha channel instead (see legendData.ts's highlightTargetFor /
// the `idCodes` branch below), read straight off each pixel with no color-distance ambiguity
// and no risk of two biomes' colors colliding.
//
// The elevReason debug view still matches on color (see legendData.ts) -- every pixel there is
// exactly one of a small fixed palette, no blending -- so `classifyPixel` stays for that case:
// nearest-neighbor across the *entire* palette, not just a within-tolerance check against the
// selected label's own colors, so one label's highlight can't bleed into a neighbour whose
// color lands nearby. `tolerance` bounds how far the *closest* match can be before a pixel
// counts as belonging to no known label at all (coastline overlay, etc).
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

// Paints the highlight dim-filter and/or (for alpha-encoded-id frames) resets alpha to
// opaque, in a single full-canvas pixel pass. Called on every frame when `alphaEncodedIds`
// so the encoded 237..255 alpha never actually composites; called only while a highlight is
// active otherwise. A no-op combination (no highlight, no alpha reset) skips the pass.
function applyHighlight(
  ctx: CanvasRenderingContext2D,
  width: number,
  height: number,
  target: HighlightTarget | null,
  resetAlpha: boolean,
): void {
  if (!target && !resetAlpha) return;
  const tolerance2 = target ? target.tolerance * target.tolerance : 0;
  const idCodes = target?.idCodes ?? null;
  const imageData = ctx.getImageData(0, 0, width, height);
  const data = imageData.data;
  for (let i = 0; i < data.length; i += 4) {
    const r = data[i], g = data[i + 1], b = data[i + 2];
    const matches = target
      ? idCodes
        ? idCodes.includes(255 - data[i + 3])
        : classifyPixel(r, g, b, target.palette, tolerance2) === target.selected
      : true;
    if (!matches) {
      const gray = (0.3 * r + 0.59 * g + 0.11 * b) * HIGHLIGHT_DIM_FACTOR;
      data[i] = gray;
      data[i + 1] = gray;
      data[i + 2] = gray;
    }
    if (resetAlpha) data[i + 3] = 255;
  }
  ctx.putImageData(imageData, 0, 0);
}

// All map drawing (fill colors, plate boundaries/poles/rotation arcs, per-plate node
// dots) happens server-side per requested view -- see backend app/render_image.py. This
// component decodes the returned PNG and paints it onto the canvas, and additionally handles
// the long-press-and-drag "rotate the planet" gesture (see rotationDrag.ts): while dragging,
// it draws a cheap wireframe graticule preview client-side (see rotation.ts) instead of
// re-requesting the real, much more expensive, detailed render on every mouse move.
// Every getContext("2d") call must pass the same options or the browser returns null on a
// mismatched re-request. willReadFrequently: the highlight / alpha-reset passes call
// getImageData on most frames.
const CTX_OPTIONS: CanvasRenderingContext2DSettings = { willReadFrequently: true };

export default function MapCanvas({
  imageBase64, width, height, displayWidth, displayHeight, projection, rotation,
  onRotationPreview, onRotationCommitted, highlightTarget, onProbe, alphaEncodedIds, interactionDisabled,
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
  const alphaEncodedIdsRef = useRef(alphaEncodedIds);
  alphaEncodedIdsRef.current = alphaEncodedIds;

  // Draws the already-decoded frame plus, if a legend highlight is active, the dim-filter on
  // top of it -- shared by both the initial decode (below) and the highlight-toggle effect, so
  // toggling a legend swatch doesn't need a fresh server render to update the map. Reads
  // highlightTarget via a ref (not a direct closure) so its identity only changes with
  // width/height, not with the highlight selection -- see that effect's own comment for why.
  const paintDecodedFrame = useCallback(() => {
    const canvas = canvasRef.current;
    const img = imgRef.current;
    const ctx = canvas?.getContext("2d", CTX_OPTIONS);
    if (!canvas || !img || !ctx) return;
    // "copy" (not the default "source-over") so the source's alpha replaces the destination
    // rather than compositing over the retained previous frame -- keeps the flicker fix (still
    // one full-canvas paint) while leaving Combined's alpha-encoded ids intact to read back.
    ctx.globalCompositeOperation = "copy";
    ctx.drawImage(img, 0, 0, width, height);
    ctx.globalCompositeOperation = "source-over";
    applyHighlight(ctx, width, height, highlightTargetRef.current ?? null, alphaEncodedIdsRef.current ?? false);
  }, [width, height]);

  useEffect(() => {
    const canvas = canvasRef.current;
    const img = imgRef.current;
    if (!canvas || !img) return;
    const ctx = canvas.getContext("2d", CTX_OPTIONS);
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
    const ctx = canvas?.getContext("2d", CTX_OPTIONS);
    const img = imgRef.current;
    if (!canvas || !ctx || !img) return;
    if (imageBase64) {
      // Fill first so a Combined frame's sub-255 alpha (see alphaEncodedIds) blends onto solid
      // background rather than showing the page through, then draws effectively opaque.
      if (alphaEncodedIds) {
        ctx.fillStyle = BACKGROUND;
        ctx.fillRect(0, 0, width, height);
      }
      ctx.drawImage(img, 0, 0, width, height);
    } else {
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

  // A plain click (routed through useRotationDrag so a completed rotate-drag's terminating
  // mouseup is never mistaken for one -- same as the inspector views): unproject it through
  // the current view rotation to a true lat/lon and hand it up for the inspect popup. Same
  // math as LakeInspector.handleClick.
  const handleProbeClick = useCallback((backingX: number, backingY: number) => {
    if (!onProbe) return;
    const transform = getRenderTransform(projection, width, height);
    const x = (backingX - transform.offsetX) / transform.scale;
    const y = -(backingY - transform.offsetY) / transform.scale;
    const latLon = unproject(projection, x, y);
    if (!latLon) {
      onProbe(null); // clicked off the globe -- dismiss any open popup
      return;
    }
    const displayXyz = latLonToXyz(latLon[0], latLon[1]); // display (rotated) frame
    const trueXyz = matApply(matTranspose(rotation), displayXyz); // transpose == inverse for a rotation matrix
    const [trueLat, trueLon] = xyzToLatLon(trueXyz);
    onProbe({
      displayX: (backingX / width) * displayWidth,
      displayY: (backingY / height) * displayHeight,
      latDeg: (trueLat * 180) / Math.PI,
      lonDeg: (trueLon * 180) / Math.PI,
    });
  }, [onProbe, projection, width, height, displayWidth, displayHeight, rotation]);

  useRotationDrag({
    elementRef: canvasRef,
    width, height, displayWidth, displayHeight, projection, rotation,
    onFrame: drawGraticule,
    onRotationPreview,
    onRotationCommitted,
    onClick: onProbe ? handleProbeClick : undefined,
    disabled: interactionDisabled,
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
