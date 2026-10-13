import { useEffect, useRef, useState } from "react";
import { fetchCrossSection } from "./api";
import { nearbyCrossSectionPositions } from "./crossSectionFrames";

const PLANET_RADIUS_KM = 6_371;
const ARC_LENGTH_KM = 3_000;
const ARC_POSITION_COUNT = 120;
const SECTION_DEPTH_KM = 350;
const VERTICAL_EXAGGERATION = 2;
const SECTION_STYLE_VERSION = 1;
const MAX_CACHED_SECTION_FRAMES = 9;
const PREFETCH_CONCURRENCY = 3;
const PREFETCHED_FRAME_COUNT = MAX_CACHED_SECTION_FRAMES;

interface FrameKey {
  worldRevision: number;
  planeAngleDeg: number;
  arcPosition: number;
  width: number;
  height: number;
}

interface CachedFrame {
  key: FrameKey;
  dataUrl: string;
  image: HTMLImageElement;
  lastUsed: number;
}

const frameCache = new Map<string, CachedFrame>();
const inFlightFrames = new Map<string, { signal: AbortSignal; promise: Promise<string> }>();
let activeFrameKey: FrameKey | null = null;
let cacheClock = 0;

function frameCacheKey(key: FrameKey): string {
  return [
    key.worldRevision, key.planeAngleDeg, key.arcPosition, key.width, key.height, SECTION_STYLE_VERSION,
  ].join(":");
}

function sameFrameSet(left: FrameKey, right: FrameKey): boolean {
  return left.worldRevision === right.worldRevision
    && left.planeAngleDeg === right.planeAngleDeg
    && left.width === right.width
    && left.height === right.height;
}

function rememberActiveFrame(key: FrameKey): void {
  activeFrameKey = key;
  pruneFrameCache();
}

function pruneFrameCache(): void {
  if (frameCache.size <= MAX_CACHED_SECTION_FRAMES) return;
  const distance = (entry: CachedFrame) => {
    if (!activeFrameKey || !sameFrameSet(entry.key, activeFrameKey)) return Number.POSITIVE_INFINITY;
    const direct = Math.abs(entry.key.arcPosition - activeFrameKey.arcPosition);
    return Math.min(direct, ARC_POSITION_COUNT - direct);
  };
  const keep = [...frameCache.entries()]
    .sort((left, right) => distance(left[1]) - distance(right[1]) || right[1].lastUsed - left[1].lastUsed)
    .slice(0, MAX_CACHED_SECTION_FRAMES);
  const retained = new Set(keep.map(([key]) => key));
  for (const key of frameCache.keys()) {
    if (!retained.has(key)) frameCache.delete(key);
  }
}

function cachedFrame(key: FrameKey): string | null {
  const entry = frameCache.get(frameCacheKey(key));
  if (!entry) return null;
  entry.lastUsed = ++cacheClock;
  return entry.dataUrl;
}

async function loadFrame(key: FrameKey, signal: AbortSignal): Promise<string> {
  const cacheKey = frameCacheKey(key);
  const existing = cachedFrame(key);
  if (existing) return existing;
  const pending = inFlightFrames.get(cacheKey);
  if (pending && !pending.signal.aborted) return pending.promise;

  const promise = (async () => {
    const response = await fetchCrossSection(key.planeAngleDeg, key.arcPosition, key.width, key.height, signal);
    const dataUrl = `data:image/png;base64,${response.image_base64}`;
    // Decode in the background before a frame is eligible for reuse so slider changes can
    // swap to a nearby cached frame without waiting for PNG decoding on the paint path.
    const image = new Image();
    if (typeof image.decode === "function") {
      image.src = dataUrl;
      await image.decode();
    } else {
      await new Promise<void>((resolve, reject) => {
        image.onload = () => resolve();
        image.onerror = () => reject(new Error("Could not decode cross-section image"));
        image.src = dataUrl;
      });
    }
    if (signal.aborted) throw new DOMException("Cross-section request aborted", "AbortError");

    frameCache.set(cacheKey, { key, dataUrl, image, lastUsed: ++cacheClock });
    pruneFrameCache();
    return dataUrl;
  })();
  const entry = { signal, promise };
  inFlightFrames.set(cacheKey, entry);
  try {
    return await promise;
  } finally {
    if (inFlightFrames.get(cacheKey) === entry) inFlightFrames.delete(cacheKey);
  }
}

function isNotImplementedError(error: unknown): boolean {
  return String(error).includes("501");
}

interface Props {
  width: number;
  height: number;
  worldRevision: number;
  planeAngleDeg: number;
  arcPosition: number;
  onPlaneAngleChange: (value: number) => void;
  onArcPositionChange: (value: number) => void;
}

interface ProjectedPoint {
  x: number;
  y: number;
  depth: number;
}

interface CurveSegment {
  front: boolean;
  points: ProjectedPoint[];
}

function globePoint(planeAngleDeg: number, arcAngle: number, cx: number, cy: number, radius: number): ProjectedPoint {
  const planeAngle = (planeAngleDeg * Math.PI) / 180;
  return {
    x: cx + radius * Math.cos(planeAngle) * Math.cos(arcAngle),
    y: cy - radius * Math.sin(arcAngle),
    depth: Math.sin(planeAngle) * Math.cos(arcAngle),
  };
}

function projectedCurveSegments(
  planeAngleDeg: number,
  startAngle: number,
  endAngle: number,
  sampleCount: number,
  cx: number,
  cy: number,
  radius: number,
): CurveSegment[] {
  const points: ProjectedPoint[] = Array.from({ length: sampleCount + 1 }, (_, index) => {
    const angle = startAngle + ((endAngle - startAngle) * index) / sampleCount;
    return globePoint(planeAngleDeg, angle, cx, cy, radius);
  });
  const segments: CurveSegment[] = [];

  for (const point of points) {
    const front = point.depth >= -0.0001;
    const current = segments[segments.length - 1];
    if (!current || current.front !== front) {
      segments.push({ front, points: [point] });
    } else {
      current.points.push(point);
    }
  }

  return segments.filter((segment) => segment.points.length > 1);
}

function pathForPoints(points: ProjectedPoint[]): string {
  return points.map((point, index) => `${index === 0 ? "M" : "L"} ${point.x.toFixed(2)} ${point.y.toFixed(2)}`).join(" ");
}

function sectionPoint(
  radiusKm: number,
  distanceKm: number,
  centerX: number,
  centerY: number,
  pixelsPerKm: number,
): ProjectedPoint {
  const arcAngle = distanceKm / PLANET_RADIUS_KM;
  return {
    x: centerX + PLANET_RADIUS_KM * Math.sin(arcAngle) * pixelsPerKm,
    y: centerY + VERTICAL_EXAGGERATION * (PLANET_RADIUS_KM - radiusKm * Math.cos(arcAngle)) * pixelsPerKm,
    depth: 0,
  };
}

function shellPath(radiusKm: number, halfArcKm: number, centerX: number, centerY: number, pixelsPerKm: number): string {
  const sampleCount = 100;
  return Array.from({ length: sampleCount + 1 }, (_, index) => {
    const distanceKm = -halfArcKm + (2 * halfArcKm * index) / sampleCount;
    const point = sectionPoint(radiusKm, distanceKm, centerX, centerY, pixelsPerKm);
    return `${index === 0 ? "M" : "L"} ${point.x.toFixed(2)} ${point.y.toFixed(2)}`;
  }).join(" ");
}

function sectionAreaPath(
  halfArcKm: number,
  centerX: number,
  centerY: number,
  pixelsPerKm: number,
): string {
  const top = shellPath(PLANET_RADIUS_KM, halfArcKm, centerX, centerY, pixelsPerKm);
  const sampleCount = 100;
  const bottom = Array.from({ length: sampleCount + 1 }, (_, index) => {
    const distanceKm = halfArcKm - (2 * halfArcKm * index) / sampleCount;
    const point = sectionPoint(PLANET_RADIUS_KM - SECTION_DEPTH_KM, distanceKm, centerX, centerY, pixelsPerKm);
    return `L ${point.x.toFixed(2)} ${point.y.toFixed(2)}`;
  }).join(" ");
  return `${top} ${bottom} Z`;
}

function SectionPlaceholder({ width, height }: { width: number; height: number }) {
  const viewWidth = width;
  const viewHeight = height;
  const marginX = Math.max(16, Math.round(viewWidth * 0.025));
  const marginTop = Math.max(30, Math.round(viewHeight * 0.075));
  const centerX = viewWidth / 2;
  const centerY = marginTop + Math.round(viewHeight * 0.08);
  const halfArcKm = ARC_LENGTH_KM / 2;
  const arcHalfAngle = halfArcKm / PLANET_RADIUS_KM;
  // Mirror cross_section.py's frame geometry so the schematic stays in place as its PNG fades in.
  const arcWidth = (viewWidth - 2 * marginX) * (690 / 820);
  const pixelsPerKm = arcWidth / (2 * PLANET_RADIUS_KM * Math.sin(arcHalfAngle));
  const surfacePath = shellPath(PLANET_RADIUS_KM, halfArcKm, centerX, centerY, pixelsPerKm);
  const areaPath = sectionAreaPath(halfArcKm, centerX, centerY, pixelsPerKm);
  return (
    <svg
      role="img"
      aria-label="Schematic 3,000 kilometre arc section to 350 kilometres below the local surface; world terrain and geology are not connected yet"
      viewBox={`0 0 ${viewWidth} ${viewHeight}`}
      width="100%"
      height="100%"
      preserveAspectRatio="xMidYMid meet"
      style={{ display: "block", minHeight: 0, flex: 1, width: "100%" }}
    >
      <rect x="0" y="0" width={viewWidth} height={viewHeight} fill="#0e1422" />
      <text x={marginX} y="17" fill="#ebeff6" fontSize="11">Surface-sampled lithosphere section · 3,000 km arc</text>
      <text x={marginX} y="31" fill="#a8b4c5" fontSize="10">350 km below local surface · vertical relief 2× · geotherm illustrative</text>
      <path d={areaPath} fill="#756C94" />
      <path d={surfacePath} fill="none" stroke="#e5e7ca" strokeWidth="1" />

      {[100, 200, 300, SECTION_DEPTH_KM].map((depthKm) => {
        const point = sectionPoint(PLANET_RADIUS_KM - depthKm, 0, centerX, centerY, pixelsPerKm);
        return (
          <g key={depthKm}>
            <line x1={centerX - 5} y1={point.y} x2={centerX + 5} y2={point.y} stroke="#e6e7ef" strokeWidth="1" />
            <text x={centerX + 8} y={point.y - 2} fill="#dee1eb" fontSize="9">{depthKm} km</text>
          </g>
        );
      })}

    </svg>
  );
}

function GlobeLocator({ planeAngleDeg, arcPosition }: { planeAngleDeg: number; arcPosition: number }) {
  const centerX = 130;
  const centerY = 130;
  const radius = 106;
  const centerAngle = (2 * Math.PI * arcPosition) / ARC_POSITION_COUNT;
  const halfArcAngle = ARC_LENGTH_KM / (2 * PLANET_RADIUS_KM);
  const wholeCircle = projectedCurveSegments(planeAngleDeg, 0, 2 * Math.PI, 720, centerX, centerY, radius);
  const selectedArc = projectedCurveSegments(
    planeAngleDeg,
    centerAngle - halfArcAngle,
    centerAngle + halfArcAngle,
    96,
    centerX,
    centerY,
    radius,
  );
  const arcCenter = globePoint(planeAngleDeg, centerAngle, centerX, centerY, radius);

  return (
    <svg
      role="img"
      aria-label={`Globe locator, plane orientation ${planeAngleDeg} degrees, arc center at position ${arcPosition + 1} of ${ARC_POSITION_COUNT}`}
      viewBox="0 0 260 260"
      style={{ width: "100%", maxWidth: 260, height: "auto", display: "block", alignSelf: "center" }}
    >
      <defs>
        <clipPath id="cross-section-globe-clip">
          <circle cx={centerX} cy={centerY} r={radius} />
        </clipPath>
      </defs>
      <circle cx={centerX} cy={centerY} r={radius} fill="#122239" fillOpacity="0.82" stroke="#8796ad" strokeWidth="1.5" />
      <g clipPath="url(#cross-section-globe-clip)" fill="none" stroke="#435875" strokeWidth="1" opacity="0.55">
        <ellipse cx={centerX} cy={centerY} rx={radius} ry={radius * 0.36} />
        <ellipse cx={centerX} cy={centerY} rx={radius * 0.36} ry={radius} />
        <ellipse cx={centerX} cy={centerY} rx={radius * 0.72} ry={radius} />
      </g>
      <g clipPath="url(#cross-section-globe-clip)" fill="none" strokeLinecap="round">
        {wholeCircle.map((segment, index) => (
          <path
            key={`circle-${index}`}
            d={pathForPoints(segment.points)}
            stroke={segment.front ? "#8797af" : "#718099"}
            strokeWidth="2"
            strokeDasharray={segment.front ? undefined : "5 5"}
            opacity={segment.front ? 0.9 : 0.66}
          />
        ))}
        {selectedArc.map((segment, index) => (
          <path
            key={`arc-${index}`}
            d={pathForPoints(segment.points)}
            stroke={segment.front ? "#61ded0" : "#d7a95c"}
            strokeWidth="5"
            strokeDasharray={segment.front ? undefined : "4 4"}
            opacity={segment.front ? 1 : 0.86}
          />
        ))}
      </g>
      {arcCenter.depth >= -0.0001 && (
        <circle cx={arcCenter.x} cy={arcCenter.y} r="4" fill="#fff4d2" stroke="#253147" strokeWidth="1.5" />
      )}
      <circle cx={centerX} cy={centerY} r={radius} fill="none" stroke="#8796ad" strokeWidth="1.5" />
    </svg>
  );
}

function LegendKey({ color, label, line = false }: { color: string; label: string; line?: boolean }) {
  return (
    <span style={{ display: "inline-flex", alignItems: "center", gap: 6, marginRight: 14, whiteSpace: "nowrap" }}>
      {line ? (
        <span aria-hidden style={{ width: 16, borderTop: `2px solid ${color}` }} />
      ) : (
        <span aria-hidden style={{ width: 12, height: 12, border: "1px solid #9aa4b5", borderRadius: 2, background: color }} />
      )}
      <span>{label}</span>
    </span>
  );
}

export default function MantleCrossSection({
  width,
  height,
  worldRevision,
  planeAngleDeg,
  arcPosition,
  onPlaneAngleChange,
  onArcPositionChange,
}: Props) {
  const centerDegrees = (360 * arcPosition) / ARC_POSITION_COUNT;
  const [displayedFrame, setDisplayedFrame] = useState<{ key: string; dataUrl: string } | null>(null);
  const [loading, setLoading] = useState(true);
  const [prefetchProgress, setPrefetchProgress] = useState<number | null>(null);
  const [loadMessage, setLoadMessage] = useState<string | null>(null);
  const scopeControllerRef = useRef<AbortController | null>(null);
  const arcPositionRef = useRef(arcPosition);
  arcPositionRef.current = arcPosition;
  const currentFrameKey: FrameKey = {
    worldRevision, planeAngleDeg, arcPosition, width: Math.round(width), height: Math.round(height),
  };
  const currentCacheKey = frameCacheKey(currentFrameKey);
  const cachedDataUrl = frameCache.get(currentCacheKey)?.dataUrl ?? null;
  // A prop change can render before its effect runs; only show a frame whose full cache key
  // still matches, so old geology never flashes after a world/plane/size change. Read the
  // decoded frame cache during render so a cached arc replaces the placeholder immediately.
  const imageDataUrl = displayedFrame?.key === currentCacheKey
    ? displayedFrame.dataUrl
    : cachedDataUrl;
  const statusMessage = loadMessage ?? (
    prefetchProgress !== null && prefetchProgress < PREFETCHED_FRAME_COUNT
      ? `Preparing nearby arc positions ${prefetchProgress}/${PREFETCHED_FRAME_COUNT}`
      : ""
  );

  const frameScopeKey = [worldRevision, planeAngleDeg, Math.round(width), Math.round(height), SECTION_STYLE_VERSION].join(":");

  useEffect(() => {
    const controller = new AbortController();
    scopeControllerRef.current = controller;
    const initialPosition = arcPositionRef.current;
    const activeKey: FrameKey = {
      worldRevision, planeAngleDeg, arcPosition: initialPosition, width: Math.round(width), height: Math.round(height),
    };
    rememberActiveFrame(activeKey);
    setPrefetchProgress(null);

    void (async () => {
      try {
        await loadFrame(activeKey, controller.signal);
        if (controller.signal.aborted) return;
        setPrefetchProgress(1);

        const positions = nearbyCrossSectionPositions(
          initialPosition, ARC_POSITION_COUNT, PREFETCHED_FRAME_COUNT - 1,
        );
        let nextPosition = 0;
        let completed = 1;
        const worker = async () => {
          while (!controller.signal.aborted) {
            const queueIndex = nextPosition++;
            if (queueIndex >= positions.length) return;
            const position = positions[queueIndex];
            const frameKey = { ...activeKey, arcPosition: position };
            try {
              await loadFrame(frameKey, controller.signal);
              if (controller.signal.aborted) return;
              completed++;
              setPrefetchProgress(completed);
            } catch (error) {
              if (controller.signal.aborted) return;
              if (isNotImplementedError(error)) {
                setLoadMessage("The cross-section route is reserved, but this backend does not render section images yet.");
                controller.abort();
                return;
              }
              setLoadMessage(`Could not cache arc positions: ${String(error)}`);
              controller.abort();
              return;
            }
          }
        };
        await Promise.all(Array.from({ length: PREFETCH_CONCURRENCY }, () => worker()));
      } catch (error) {
        if (controller.signal.aborted) return;
        setLoading(false);
        if (isNotImplementedError(error)) {
          setLoadMessage("The cross-section route is reserved, but this backend does not render section images yet.");
        } else {
          setLoadMessage(`Could not load this cross-section: ${String(error)}`);
        }
      }
    })();

    return () => {
      controller.abort();
      if (scopeControllerRef.current === controller) scopeControllerRef.current = null;
    };
  }, [worldRevision, planeAngleDeg, width, height, frameScopeKey]);

  useEffect(() => {
    const activeKey: FrameKey = {
      worldRevision, planeAngleDeg, arcPosition, width: Math.round(width), height: Math.round(height),
    };
    rememberActiveFrame(activeKey);
    const cached = cachedFrame(activeKey);
    setDisplayedFrame(cached ? { key: currentCacheKey, dataUrl: cached } : null);
    setLoading(cached === null);
    setLoadMessage(null);
    let isCurrentSelection = true;
    const controller = scopeControllerRef.current;

    if (!cached && controller) {
      void loadFrame(activeKey, controller.signal)
        .then((image) => {
          if (!isCurrentSelection || controller.signal.aborted) return;
          setDisplayedFrame({ key: currentCacheKey, dataUrl: image });
          setLoading(false);
        })
        .catch((error) => {
          if (!isCurrentSelection || controller.signal.aborted) return;
          setLoading(false);
          setLoadMessage(isNotImplementedError(error)
            ? "The cross-section route is reserved, but this backend does not render section images yet."
            : `Could not load this cross-section: ${String(error)}`);
        });
    }

    return () => {
      isCurrentSelection = false;
    };
  }, [worldRevision, planeAngleDeg, arcPosition, width, height, currentCacheKey]);

  return (
    <div
      style={{
        boxSizing: "border-box",
        display: "grid",
        gridTemplateColumns: "minmax(0, 1fr) 286px",
        gap: 12,
        width,
        height,
        padding: 12,
        overflow: "hidden",
        color: "#e5eaf2",
        background: "#0e1422",
        border: "1px solid #2c3447",
        borderRadius: 6,
      }}
    >
      <section style={{ minWidth: 0, minHeight: 0, display: "flex", flexDirection: "column", gap: 5 }}>
        <div style={{ display: "flex", alignItems: "baseline", justifyContent: "space-between", gap: 8 }}>
          <h2 style={{ margin: 0, fontSize: 17, fontWeight: 600 }}>Shallow mantle cross-section</h2>
          <span style={{ color: imageDataUrl ? "#8fe0b0" : "#d6bb85", fontSize: 9, letterSpacing: "0.09em", whiteSpace: "nowrap" }}>
            {imageDataUrl ? "SAMPLED SECTION" : loading ? "PREPARING CROSS-SECTION" : "SCHEMATIC PLACEHOLDER"}
          </span>
        </div>
        <p style={{ margin: "0 0 2px", color: "#a9b3c4", fontSize: 11 }}>
          A stable local frame follows the selected surface arc as its center moves around the planet.
        </p>
        <div style={{ display: "flex", position: "relative", minHeight: 0, flex: 1, width: "100%" }}>
          <div aria-hidden="true" style={{ display: "flex", minHeight: 0, flex: 1, width: "100%" }}>
            {!imageDataUrl && <SectionPlaceholder width={width} height={height} />}
          </div>
          {imageDataUrl && (
            <img
              src={imageDataUrl}
              alt={`Rendered shallow mantle section at arc position ${arcPosition + 1} of ${ARC_POSITION_COUNT}`}
              style={{
                position: "absolute",
                inset: 0,
                display: "block",
                minHeight: 0,
                width: "100%",
                height: "100%",
                objectFit: "contain",
                pointerEvents: "none",
              }}
            />
          )}
        </div>
        <div
          role="status"
          style={{
            boxSizing: "border-box",
            height: 12,
            flexShrink: 0,
            overflow: "hidden",
            color: loadMessage ? "#d6bb85" : "#9daac0",
            fontSize: 10,
            lineHeight: "12px",
            whiteSpace: "nowrap",
            textOverflow: "ellipsis",
            visibility: statusMessage ? "visible" : "hidden",
          }}
        >
          {statusMessage}
        </div>
        <div style={{ color: "#8f9bb0", fontSize: 9, display: "flex", flexWrap: "wrap", gap: 4 }}>
          <LegendKey color="#539FCB" label="Water" />
          <LegendKey color="#E0EDF2" label="Stored glacier" />
          <LegendKey color="#C98F82" label="Continental crust" />
          <LegendKey color="#526779" label="Oceanic crust" />
          <LegendKey color="#425A5A" label="Mantle lithosphere" />
          <LegendKey color="#756C94" label="Schematic asthenosphere" />
          <LegendKey color="#F47A3C" label="Stored melt source" />
          <LegendKey color="#FFD278" label="Melt-eligible estimate" />
          <LegendKey color="#47ECD2" label="Mantle flow (schematic depth)" line />
          <LegendKey color="#ED7E3C" label="Active vent" />
          <LegendKey color="#B7BEC1" label="Dormant vent" />
        </div>
      </section>

      <aside style={{ minWidth: 0, minHeight: 0, display: "flex", flexDirection: "column", gap: 7, borderLeft: "1px solid #303a4d", paddingLeft: 12 }}>
        <h3 style={{ margin: "1px 0 0", fontSize: 13, fontWeight: 600 }}>Globe locator</h3>
        <GlobeLocator planeAngleDeg={planeAngleDeg} arcPosition={arcPosition} />
        <div style={{ color: "#aeb8ca", fontSize: 10, textAlign: "center", marginTop: -3 }}>
          <span style={{ color: "#91a0b8" }}>Selected great circle</span>
          <span aria-hidden> · </span>
          <span style={{ color: "#61ded0" }}>3,000 km arc</span>
        </div>

        <label style={{ display: "block", marginTop: 2, fontSize: 11 }}>
          <span style={{ display: "flex", justifyContent: "space-between", gap: 6 }}>
            <span>Plane orientation</span>
            <span style={{ color: "#b8c3d3", fontVariantNumeric: "tabular-nums" }}>{planeAngleDeg}°</span>
          </span>
          <input
            type="range"
            min={0}
            max={180}
            step={1}
            value={planeAngleDeg}
            onChange={(event) => onPlaneAngleChange(Number(event.target.value))}
            aria-label="Great-circle plane orientation"
            style={{ width: "100%", margin: "3px 0 0", accentColor: "#61ded0" }}
          />
        </label>
        <div style={{ color: "#7f8ba0", fontSize: 9, marginTop: -7 }}>Rotates the great-circle plane through the globe.</div>

        <label style={{ display: "block", marginTop: 1, fontSize: 11 }}>
          <span style={{ display: "flex", justifyContent: "space-between", gap: 6 }}>
            <span>Arc center position</span>
            <span style={{ color: "#b8c3d3", fontVariantNumeric: "tabular-nums" }}>{arcPosition + 1} / {ARC_POSITION_COUNT}</span>
          </span>
          <input
            type="range"
            min={0}
            max={ARC_POSITION_COUNT - 1}
            step={1}
            value={arcPosition}
            onChange={(event) => onArcPositionChange(Number(event.target.value))}
            onKeyDown={(event) => {
              if (event.key === "ArrowRight" && arcPosition === ARC_POSITION_COUNT - 1) {
                event.preventDefault();
                onArcPositionChange(0);
              } else if (event.key === "ArrowLeft" && arcPosition === 0) {
                event.preventDefault();
                onArcPositionChange(ARC_POSITION_COUNT - 1);
              }
            }}
            aria-label="Position of the 3,000 kilometre arc around the great circle"
            style={{ width: "100%", margin: "3px 0 0", accentColor: "#61ded0" }}
          />
        </label>
        <div style={{ color: "#7f8ba0", fontSize: 9, marginTop: -7 }}>
          Center {centerDegrees.toFixed(0)}° of 360° · {ARC_LENGTH_KM.toLocaleString()} km fixed length
        </div>
      </aside>
    </div>
  );
}
