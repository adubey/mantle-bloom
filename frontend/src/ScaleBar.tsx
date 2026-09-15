import { useRef, useState } from "react";
import type { RefObject } from "react";
import type { Projection } from "./api";
import { PLANET_RADIUS_KM, backingPixelsToDisplayLatLon, getLocalPixelScale, type LocalPixelScale } from "./rotation";

// The map's scale bar (issue #158): normally a plain inline item in the Legend row below the
// map, exact only at the equator (like any small-scale world map's static scale bar -- see the
// comment on KM_PER_PIXEL_AT_EQUATOR below). Grabbing and dragging it onto the map instead
// *detaches* it into a `position: fixed` overlay pinned to that screen point, at which point it
// switches to a live per-position calculation (rotation.ts's getLocalPixelScale) so its labeled
// length stays accurate wherever it's been dropped -- dragging it back off the map re-docks it,
// reverting to the original always-equator reading unchanged. Space bar (while it's focused,
// only meaningful once detached) flips it between measuring horizontally (east-west) and
// vertically (north-south).

interface Props {
  projection: Projection;
  width: number; // backing-store px, matching MapCanvas's own width/height props
  height: number;
  displayWidth: number; // CSS px
  displayHeight: number;
  mapWrapperRef: RefObject<HTMLDivElement | null>;
}

type Orientation = "horizontal" | "vertical";
interface Tick {
  offsetPx: number;
  major: boolean;
}

// The map's own km-per-pixel scale while docked in the Legend -- unchanged from this bar's
// original, non-draggable version. The world is always Earth-sized (see backend
// app/elevation_lines.py's PLANET_RADIUS_KM) and both current projections (behrmann/eckert4)
// fit the whole sphere to the rendered width/height with their widest point at the equator, so
// this reading is just circumference / available map width -- exact only at the equator, same
// caveat any small-scale world map's static scale bar carries. Once dropped onto the map,
// getLocalPixelScale below replaces this with an actual per-position reading instead.
const MAP_DISPLAY_WIDTH_PX = 1100; // matches App.tsx's DISPLAY_WIDTH
const MAP_PADDING_PX = 20; // matches backend app/render_image.py's PADDING_PX
const KM_PER_PIXEL_AT_EQUATOR = (2 * Math.PI * PLANET_RADIUS_KM) / (MAP_DISPLAY_WIDTH_PX - 2 * MAP_PADDING_PX);
const PX_PER_KM_AT_EQUATOR = 1 / KM_PER_PIXEL_AT_EQUATOR;
const DOCKED_TARGET_PX = 150;
const DOCKED_MAJOR_STEP_KM = 1000;
const DOCKED_MINOR_STEP_KM = 100;

// Once dropped, the "nice" round distance varies wildly by position/zoom (could be 5 km or
// 5,000 km), so ticks are drawn as even fractions of whatever that distance turns out to be,
// rather than the docked bar's fixed 1,000/100 km major/minor scheme.
const DROPPED_TARGET_PX = 150;
const DROPPED_TICK_SUBDIVISIONS = 5;

const BASELINE_INSET = 4;
const MAJOR_TICK_LEN = 7;
const MINOR_TICK_LEN = 3.5;
const STROKE = "#dee2eb";
// A plain click (mousedown+mouseup with negligible movement) never lifts the bar out of the
// Legend -- only an actual drag past this many viewport px does.
const LIFT_THRESHOLD_PX = 4;

// Standard cartographic 1-2-5 rounding, so the labeled distance is always a clean number
// instead of whatever `DROPPED_TARGET_PX / pxPerKm` happens to work out to.
function niceScaleKm(roughKm: number): number {
  if (!(roughKm > 0)) return 1;
  const magnitude = Math.pow(10, Math.floor(Math.log10(roughKm)));
  const residual = roughKm / magnitude;
  const niceResidual = residual < 1.5 ? 1 : residual < 3.5 ? 2 : residual < 7.5 ? 5 : 10;
  return niceResidual * magnitude;
}

export default function ScaleBar({ projection, width, height, displayWidth, displayHeight, mapWrapperRef }: Props) {
  const [dropped, setDropped] = useState<{ left: number; top: number } | null>(null); // viewport px
  const [orientation, setOrientation] = useState<Orientation>("horizontal");
  const barRef = useRef<HTMLDivElement | null>(null);
  // Kept from the last position that actually landed on the globe, so a drag that briefly
  // crosses into the map's background padding doesn't make the bar disappear or snap to zero
  // mid-gesture -- it just holds its last reading until the pointer is back over the sphere.
  const lastScale = useRef<LocalPixelScale>({ pxPerKmHorizontal: 0, pxPerKmVertical: 0 });

  const isDropped = dropped !== null;
  const vertical = isDropped && orientation === "vertical";

  const onMouseDown = (e: React.MouseEvent) => {
    if (e.button !== 0) return;
    e.preventDefault();
    const rect = barRef.current?.getBoundingClientRect();
    const grabDx = e.clientX - (rect?.left ?? e.clientX);
    const grabDy = e.clientY - (rect?.top ?? e.clientY);
    const startClientX = e.clientX;
    const startClientY = e.clientY;
    // Repositioning an already-dropped bar starts moving immediately; picking it up fresh out
    // of the Legend needs an actual drag past LIFT_THRESHOLD_PX first, so a plain click doesn't
    // yank it out of its docked spot.
    let lifted = isDropped;

    const onMove = (moveEvent: MouseEvent) => {
      if (!lifted) {
        const dx = moveEvent.clientX - startClientX;
        const dy = moveEvent.clientY - startClientY;
        if (Math.hypot(dx, dy) < LIFT_THRESHOLD_PX) return;
        lifted = true;
      }
      setDropped({ left: moveEvent.clientX - grabDx, top: moveEvent.clientY - grabDy });
    };
    const onUp = (upEvent: MouseEvent) => {
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
      if (!lifted) return; // a plain click -- stays docked, nothing to resolve

      const mapRect = mapWrapperRef.current?.getBoundingClientRect();
      const onMap =
        !!mapRect &&
        upEvent.clientX >= mapRect.left &&
        upEvent.clientX <= mapRect.right &&
        upEvent.clientY >= mapRect.top &&
        upEvent.clientY <= mapRect.bottom;
      if (onMap) {
        setDropped({ left: upEvent.clientX - grabDx, top: upEvent.clientY - grabDy });
      } else {
        // Released back over the Legend/sidebar, not the map -- snap back to docked.
        setDropped(null);
        setOrientation("horizontal");
      }
    };
    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);
    barRef.current?.focus();
  };

  const onKeyDown = (e: React.KeyboardEvent) => {
    if (!isDropped) return; // orientation only means anything once it's measuring on the map
    if (e.key === " " || e.code === "Space") {
      // Prevents the browser's default page-scroll/button-activation behavior for space.
      e.preventDefault();
      setOrientation((o) => (o === "horizontal" ? "vertical" : "horizontal"));
    }
  };

  let label: number;
  let lengthPx: number;
  let ticks: Tick[];
  if (isDropped && dropped) {
    const mapRect = mapWrapperRef.current?.getBoundingClientRect();
    if (mapRect) {
      const backingX = ((dropped.left - mapRect.left) / displayWidth) * width;
      const backingY = ((dropped.top - mapRect.top) / displayHeight) * height;
      const displayLatLon = backingPixelsToDisplayLatLon(projection, width, height, backingX, backingY);
      if (displayLatLon) {
        lastScale.current = getLocalPixelScale(projection, width, height, displayLatLon[0], displayLatLon[1]);
      }
    }
    // getLocalPixelScale reports backing-store px/km (like rotation.ts's own
    // getPixelsPerRadian); the bar itself is drawn in CSS/viewport px, which can be a different
    // resolution (see App.tsx's RENDER_SCALE), so convert before using it for anything on-screen.
    const backingToDisplay = displayWidth / width;
    const pxPerKmBacking = orientation === "horizontal" ? lastScale.current.pxPerKmHorizontal : lastScale.current.pxPerKmVertical;
    const pxPerKm = pxPerKmBacking * backingToDisplay;
    const roughKm = pxPerKm > 1e-9 && Number.isFinite(pxPerKm) ? DROPPED_TARGET_PX / pxPerKm : 0;
    label = niceScaleKm(roughKm);
    lengthPx = Number.isFinite(pxPerKm) && pxPerKm > 0 ? Math.max(2, label * pxPerKm) : DROPPED_TARGET_PX;
    const stepPx = lengthPx / DROPPED_TICK_SUBDIVISIONS;
    ticks = Array.from({ length: DROPPED_TICK_SUBDIVISIONS + 1 }, (_, i) => ({
      offsetPx: i * stepPx,
      major: i === 0 || i === DROPPED_TICK_SUBDIVISIONS,
    }));
  } else {
    // Docked: the original always-equator reading and its 1,000/100 km tick scheme, unchanged
    // from before this bar was draggable.
    const roughKm = KM_PER_PIXEL_AT_EQUATOR * DOCKED_TARGET_PX;
    label = Math.max(DOCKED_MAJOR_STEP_KM, Math.floor(roughKm / DOCKED_MAJOR_STEP_KM) * DOCKED_MAJOR_STEP_KM);
    lengthPx = label * PX_PER_KM_AT_EQUATOR;
    ticks = [];
    for (let km = 0; km <= label; km += DOCKED_MINOR_STEP_KM) {
      ticks.push({ offsetPx: km * PX_PER_KM_AT_EQUATOR, major: km % DOCKED_MAJOR_STEP_KM === 0 });
    }
  }

  const svgWidth = vertical ? MAJOR_TICK_LEN + BASELINE_INSET + 1 : lengthPx + 1;
  const svgHeight = vertical ? lengthPx + 1 : MAJOR_TICK_LEN + BASELINE_INSET + 1;

  return (
    <>
      {/* Fills the Legend's own slot once the real bar has been dragged off onto the map, so
          the row doesn't just collapse -- also carries the space-to-rotate reminder, since the
          bar itself is no longer sitting here to show it. */}
      {isDropped && (
        <div style={{ flexShrink: 0, marginLeft: "auto", fontSize: 10, opacity: 0.6, textAlign: "right", alignSelf: "center" }}>
          <div>Scale moved to the map --</div>
          <div>press &apos;space&apos; to rotate</div>
        </div>
      )}
      <div
        ref={barRef}
        tabIndex={0}
        onMouseDown={onMouseDown}
        onKeyDown={onKeyDown}
        title={
          isDropped
            ? "Drag to measure elsewhere on the map. Press space to rotate horizontal/vertical."
            : `Ticks every ${DOCKED_MAJOR_STEP_KM.toLocaleString()}/${DOCKED_MINOR_STEP_KM} km. Exact at the equator -- drag onto the map to measure accurately elsewhere.`
        }
        style={
          isDropped
            ? {
                position: "fixed",
                left: dropped!.left,
                top: dropped!.top,
                zIndex: 20,
                cursor: "grab",
                outline: "none",
                userSelect: "none",
                display: "flex",
                flexDirection: vertical ? "row" : "column",
                alignItems: vertical ? "center" : "flex-start",
                gap: 2,
              }
            : {
                flexShrink: 0,
                marginLeft: "auto",
                fontSize: 10,
                opacity: 0.85,
                textAlign: "right",
                cursor: "grab",
                outline: "none",
                userSelect: "none",
              }
        }
      >
        <svg width={svgWidth} height={svgHeight} aria-hidden>
          {!vertical ? (
            <>
              <line x1={0} y1={BASELINE_INSET} x2={lengthPx} y2={BASELINE_INSET} stroke={STROKE} strokeWidth={1} />
              {ticks.map(({ offsetPx, major }) => (
                <line
                  key={offsetPx}
                  x1={offsetPx}
                  y1={BASELINE_INSET}
                  x2={offsetPx}
                  y2={BASELINE_INSET + (major ? MAJOR_TICK_LEN : MINOR_TICK_LEN)}
                  stroke={STROKE}
                  strokeWidth={major ? 1.5 : 1}
                />
              ))}
            </>
          ) : (
            <>
              <line x1={BASELINE_INSET} y1={0} x2={BASELINE_INSET} y2={lengthPx} stroke={STROKE} strokeWidth={1} />
              {ticks.map(({ offsetPx, major }) => (
                <line
                  key={offsetPx}
                  x1={BASELINE_INSET}
                  y1={offsetPx}
                  x2={BASELINE_INSET + (major ? MAJOR_TICK_LEN : MINOR_TICK_LEN)}
                  y2={offsetPx}
                  stroke={STROKE}
                  strokeWidth={major ? 1.5 : 1}
                />
              ))}
            </>
          )}
        </svg>
        <div style={{ marginTop: isDropped ? 0 : 2, fontSize: 10, color: STROKE, opacity: 0.9, whiteSpace: "nowrap" }}>
          {label.toLocaleString()} km{isDropped ? "" : " at equator"}
        </div>
      </div>
    </>
  );
}
