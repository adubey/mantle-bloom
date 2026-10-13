import { useEffect, useState } from "react";
import { fetchCrossSection } from "./api";

const PLANET_RADIUS_KM = 6_371;
const ARC_LENGTH_KM = 3_000;
const ARC_POSITION_COUNT = 120;
const SECTION_DEPTH_KM = 350;
const VERTICAL_EXAGGERATION = 2;
const SECTION_STYLE_VERSION = 1;
const MAX_CACHED_SECTION_FRAMES = 9;
const PREFETCH_CONCURRENCY = 3;

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

function SectionPlaceholder() {
  const viewWidth = 820;
  const viewHeight = 420;
  const centerX = 425;
  const centerY = 80;
  const halfArcKm = ARC_LENGTH_KM / 2;
  const arcHalfAngle = halfArcKm / PLANET_RADIUS_KM;
  const halfSurfaceWidthKm = PLANET_RADIUS_KM * Math.sin(arcHalfAngle);
  const pixelsPerKm = 690 / (2 * halfSurfaceWidthKm);
  const surfacePath = shellPath(PLANET_RADIUS_KM, halfArcKm, centerX, centerY, pixelsPerKm);
  const areaPath = sectionAreaPath(halfArcKm, centerX, centerY, pixelsPerKm);
  const centerSurface = sectionPoint(PLANET_RADIUS_KM, 0, centerX, centerY, pixelsPerKm);

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
      <defs>
        <linearGradient id="cross-section-placeholder-fill" x1="0" y1="0" x2="0" y2="1">
          <stop offset="0%" stopColor="#756C94" stopOpacity="0.48" />
          <stop offset="100%" stopColor="#423D5A" stopOpacity="0.58" />
        </linearGradient>
        <pattern id="cross-section-pending-hatch" width="8" height="8" patternTransform="rotate(45)" patternUnits="userSpaceOnUse">
          <line x1="0" y1="0" x2="0" y2="8" stroke="#c9b9e8" strokeOpacity="0.12" strokeWidth="2" />
        </pattern>
      </defs>
      <rect x="0.5" y="0.5" width={viewWidth - 1} height={viewHeight - 1} rx="9" fill="#101827" stroke="#2d394d" />

      <text x="24" y="28" fill="#e7ebf3" fontSize="15" fontWeight="600">Schematic section · 3,000 km surface arc</text>
      <text x="796" y="28" textAnchor="end" fill="#a6b0c2" fontSize="11">local tangent left → right · radial up</text>

      <path d={areaPath} fill="url(#cross-section-placeholder-fill)" />
      <path d={areaPath} fill="url(#cross-section-pending-hatch)" />

      {[-1500, -750, 0, 750, 1500].map((distanceKm) => {
        const point = sectionPoint(PLANET_RADIUS_KM, distanceKm, centerX, centerY, pixelsPerKm);
        return (
          <g key={distanceKm}>
            <line x1={point.x} y1={point.y - 20} x2={point.x} y2={point.y + 5} stroke="#64738a" strokeWidth="1" />
            <text x={point.x} y={point.y - 27} textAnchor="middle" fill="#9aa8bd" fontSize="10">
              {distanceKm > 0 ? `+${distanceKm}` : distanceKm}
            </text>
          </g>
        );
      })}

      {[100, 200, SECTION_DEPTH_KM].map((depthKm) => {
        const path = shellPath(PLANET_RADIUS_KM - depthKm, halfArcKm, centerX, centerY, pixelsPerKm);
        const point = sectionPoint(PLANET_RADIUS_KM - depthKm, 0, centerX, centerY, pixelsPerKm);
        return (
          <g key={depthKm}>
            <path
              d={path}
              fill="none"
              stroke={depthKm === SECTION_DEPTH_KM ? "#d5d8e6" : "#9f9db5"}
              strokeWidth={depthKm === SECTION_DEPTH_KM ? 1.5 : 1}
              strokeDasharray={depthKm === SECTION_DEPTH_KM ? "5 4" : "2 5"}
              opacity={depthKm === SECTION_DEPTH_KM ? 0.95 : 0.48}
            />
            <line x1={centerX - 7} y1={point.y} x2={centerX + 7} y2={point.y} stroke="#d5d8e6" strokeWidth="1" opacity="0.8" />
            <text x={centerX + 12} y={point.y + 4} fill="#c3c7d6" fontSize="10">
              {depthKm} km
            </text>
          </g>
        );
      })}

      <path d={surfacePath} fill="none" stroke="#a8d19e" strokeWidth="3" />
      <text x={centerX - 14} y={centerSurface.y - 8} textAnchor="end" fill="#c8e1c0" fontSize="10">idealized surface</text>

      <g>
        <rect x="231" y="286" width="389" height="64" rx="7" fill="#111726" fillOpacity="0.94" stroke="#3a4355" />
        <text x="425.5" y="310" textAnchor="middle" fill="#e0e4ed" fontSize="12" fontWeight="600">World section data placeholder</text>
        <text x="425.5" y="331" textAnchor="middle" fill="#a7b1c2" fontSize="11">Terrain and geological layers will appear when sampled fields are connected.</text>
      </g>

      <text x="24" y="390" fill="#a7b1c2" fontSize="10">Vertical scale exaggerated {VERTICAL_EXAGGERATION}× (fixed); the section window remains 350 km below each local surface point.</text>
      <text x="796" y="390" textAnchor="end" fill="#a7b1c2" fontSize="10">Depth labels are below the local surface</text>
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
  const currentFrameKey: FrameKey = {
    worldRevision, planeAngleDeg, arcPosition, width: Math.round(width), height: Math.round(height),
  };
  const currentCacheKey = frameCacheKey(currentFrameKey);
  // A prop change can render before its effect runs; only show a frame whose full cache key
  // still matches, so old geology never flashes after a world/plane/size change.
  const imageDataUrl = displayedFrame?.key === currentCacheKey ? displayedFrame.dataUrl : null;

  useEffect(() => {
    const controller = new AbortController();
    const activeKey: FrameKey = {
      worldRevision, planeAngleDeg, arcPosition, width: Math.round(width), height: Math.round(height),
    };
    rememberActiveFrame(activeKey);
    const cached = cachedFrame(activeKey);
    setDisplayedFrame(cached ? { key: currentCacheKey, dataUrl: cached } : null);
    setLoading(true);
    setPrefetchProgress(null);
    setLoadMessage(null);

    void (async () => {
      try {
        const selectedImage = await loadFrame(activeKey, controller.signal);
        if (controller.signal.aborted) return;
        setDisplayedFrame({ key: currentCacheKey, dataUrl: selectedImage });
        setLoading(false);
        setPrefetchProgress(1);

        const positions = Array.from({ length: ARC_POSITION_COUNT }, (_, index) => index)
          .filter((index) => index !== arcPosition);
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

    return () => controller.abort();
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
          {imageDataUrl ? (
            <img
              src={imageDataUrl}
              alt={`Rendered shallow mantle section at arc position ${arcPosition + 1} of ${ARC_POSITION_COUNT}`}
              style={{ display: "block", minHeight: 0, width: "100%", height: "100%", objectFit: "contain" }}
            />
          ) : (
            <SectionPlaceholder />
          )}
        </div>
        {loadMessage && <div role="status" style={{ color: "#d6bb85", fontSize: 10 }}>{loadMessage}</div>}
        {!loadMessage && prefetchProgress !== null && prefetchProgress < ARC_POSITION_COUNT && (
          <div role="status" style={{ color: "#9daac0", fontSize: 10 }}>
            Caching arc positions {prefetchProgress}/{ARC_POSITION_COUNT}
          </div>
        )}
        <div style={{ color: "#8f9bb0", fontSize: 10, display: "flex", flexWrap: "wrap", gap: 4 }}>
          {imageDataUrl ? (
            <>
              <LegendKey color="#C98F82" label="Continental crust" />
              <LegendKey color="#526779" label="Oceanic crust" />
              <LegendKey color="#425A5A" label="Mantle lithosphere" />
              <LegendKey color="#756C94" label="Schematic asthenosphere" />
              <LegendKey color="#C93D3D" label="Stored magma / melt" />
              <LegendKey color="#d5d8e6" label="Melt eligibility estimate" line />
            </>
          ) : (
            <>
              <LegendKey color="#a8d19e" label="Idealized sphere surface" line />
              <LegendKey color="#756c94" label="Subsurface placeholder" />
              <LegendKey color="#d5d8e6" label="350 km depth limit" line />
            </>
          )}
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
