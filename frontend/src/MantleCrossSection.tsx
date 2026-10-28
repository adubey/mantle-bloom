import { useEffect, useRef, useState } from "react";
import { fetchCrossSection } from "./api";
import { advanceCrossSectionPosition, nearbyCrossSectionPositions } from "./crossSectionFrames";
import { loadWorldMapPixels, paintGlobe } from "./globeTexture";

const PLANET_RADIUS_KM = 6_371;
const ARC_LENGTH_KM = 3_000;
const ARC_POSITION_COUNT = 120;
const SECTION_DEPTH_KM = 350;
const VERTICAL_EXAGGERATION = 2;
const SECTION_STYLE_VERSION = 2;
const MAX_CACHED_SECTION_FRAMES = ARC_POSITION_COUNT;
const PREFETCH_CONCURRENCY = 3;
const PREFETCHED_FRAME_COUNT = MAX_CACHED_SECTION_FRAMES;
const ARC_KEYBOARD_STEP = 2;
const TILT_LIMIT_DEG = 90;
const TILT_STEP_DEG = 0.25;
const TILT_KEYBOARD_STEP_DEG = 1;
const LOCATOR_SIZE = 260;

interface FrameKey {
  worldId: string;
  worldRevision: number;
  planeAngleDeg: number;
  arcPosition: number;
  width: number;
  height: number;
  showBoundaries: boolean;
  showFaults: boolean;
  showPlumes: boolean;
}

interface CachedFrame {
  key: FrameKey;
  dataUrl: string;
  lastUsed: number;
}

const frameCache = new Map<string, CachedFrame>();
const inFlightFrames = new Map<string, { signal: AbortSignal; promise: Promise<string> }>();
let activeFrameKey: FrameKey | null = null;
let cacheClock = 0;

function frameCacheKey(key: FrameKey): string {
  return [
    key.worldId, key.worldRevision, key.planeAngleDeg, key.arcPosition, key.width, key.height,
    Number(key.showBoundaries), Number(key.showFaults), Number(key.showPlumes), SECTION_STYLE_VERSION,
  ].join(":");
}

function sameFrameSet(left: FrameKey, right: FrameKey): boolean {
  return left.worldId === right.worldId
    && left.worldRevision === right.worldRevision
    && left.planeAngleDeg === right.planeAngleDeg
    && left.width === right.width
    && left.height === right.height
    && left.showBoundaries === right.showBoundaries
    && left.showFaults === right.showFaults
    && left.showPlumes === right.showPlumes;
}

function rememberActiveFrame(key: FrameKey): void {
  activeFrameKey = key;
  pruneFrameCache();
}

function pruneFrameCache(): void {
  for (const [key, entry] of frameCache) {
    if (!activeFrameKey || !sameFrameSet(entry.key, activeFrameKey)) frameCache.delete(key);
  }
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

function clearMantleCrossSectionCache(): void {
  frameCache.clear();
  inFlightFrames.clear();
  activeFrameKey = null;
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
    const response = await fetchCrossSection(
      key.worldId, key.planeAngleDeg, key.arcPosition, key.width, key.height, signal, key.showBoundaries, key.showFaults, key.showPlumes,
    );
    const dataUrl = `data:image/png;base64,${response.image_base64}`;
    if (signal.aborted) throw new DOMException("Cross-section request aborted", "AbortError");

    // Keep the compressed PNG string for every arc; the browser decodes only the selected
    // frame when it is attached to the visible image.
    frameCache.set(cacheKey, { key, dataUrl, lastUsed: ++cacheClock });
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
  worldId: string;
  width: number;
  height: number;
  worldRevision: number;
  planeAngleDeg: number;
  arcPosition: number;
  showBoundaries: boolean;
  showFaults: boolean;
  showPlumes: boolean;
  onLayerToggle: (layer: "boundaries" | "faults" | "plumes", enabled: boolean) => void;
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

type Vec3 = [number, number, number];

function locatorVector(planeAngleDeg: number, arcAngle: number): Vec3 {
  const plane = (planeAngleDeg * Math.PI) / 180;
  return [Math.cos(plane) * Math.cos(arcAngle), Math.sin(arcAngle), Math.sin(plane) * Math.cos(arcAngle)];
}

function locatorBasis(planeAngleDeg: number, centerAngle: number): [Vec3, Vec3, Vec3] {
  const center = locatorVector(planeAngleDeg, centerAngle);
  const plane = (planeAngleDeg * Math.PI) / 180;
  const right: Vec3 = [-Math.sin(centerAngle) * Math.cos(plane), Math.cos(centerAngle), -Math.sin(centerAngle) * Math.sin(plane)];
  const up: Vec3 = [
    center[1] * right[2] - center[2] * right[1],
    center[2] * right[0] - center[0] * right[2],
    center[0] * right[1] - center[1] * right[0],
  ];
  return [center, right, up];
}

function locatorPoint(point: Vec3, center: Vec3, right: Vec3, up: Vec3, cx: number, cy: number, radius: number): ProjectedPoint {
  const dot = (a: Vec3, b: Vec3) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  return { x: cx + radius * dot(point, right), y: cy - radius * dot(point, up), depth: dot(point, center) };
}

function projectedCurveSegments(
  planeAngleDeg: number,
  centerAngle: number,
  startAngle: number,
  endAngle: number,
  sampleCount: number,
  cx: number,
  cy: number,
  radius: number,
): CurveSegment[] {
  const [center, right, up] = locatorBasis(planeAngleDeg, centerAngle);
  const points: ProjectedPoint[] = Array.from({ length: sampleCount + 1 }, (_, index) => {
    const angle = startAngle + ((endAngle - startAngle) * index) / sampleCount;
    return locatorPoint(locatorVector(planeAngleDeg, angle), center, right, up, cx, cy, radius);
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

function clampTilt(degrees: number): number {
  const stepped = Math.round(degrees / TILT_STEP_DEG) * TILT_STEP_DEG;
  return Math.max(-TILT_LIMIT_DEG, Math.min(TILT_LIMIT_DEG, stepped));
}

interface GlobeLocatorProps {
  worldId: string;
  planeAngleDeg: number;
  arcPosition: number;
  worldRevision: number;
  onPlaneAngleChange: (value: number) => void;
  // The tilt under an in-progress drag (null when idle), held by the parent so its readout
  // follows the drag.
  previewTilt: number | null;
  onPreviewTiltChange: (value: number | null) => void;
}

// Dragging the globe up or down tilts the great circle (#497, replacing the tilt slider):
// moving the pointer by the globe's radius tilts it by one radian. The globe follows the
// pointer, and the new tilt is committed (and its section frames fetched) on release.
function GlobeLocator({
  worldId, planeAngleDeg: committedTilt, arcPosition, worldRevision, onPlaneAngleChange, previewTilt, onPreviewTiltChange: setPreviewTilt,
}: GlobeLocatorProps) {
  const [textureUrl, setTextureUrl] = useState<string | null>(null);
  const drag = useRef<{ pointerId: number; startY: number; unitsPerPixel: number } | null>(null);
  const planeAngleDeg = previewTilt ?? committedTilt;
  const centerX = LOCATOR_SIZE / 2;
  const centerY = LOCATOR_SIZE / 2;
  const radius = 106;
  const centerAngle = (2 * Math.PI * arcPosition) / ARC_POSITION_COUNT;
  const halfArcAngle = ARC_LENGTH_KM / (2 * PLANET_RADIUS_KM);
  const wholeCircle = projectedCurveSegments(planeAngleDeg, centerAngle, 0, 2 * Math.PI, 720, centerX, centerY, radius);
  const selectedArc = projectedCurveSegments(
    planeAngleDeg,
    centerAngle,
    centerAngle - halfArcAngle,
    centerAngle + halfArcAngle,
    96,
    centerX,
    centerY,
    radius,
  );
  const [center, right, up] = locatorBasis(planeAngleDeg, centerAngle);
  const arcCenter = locatorPoint(center, center, right, up, centerX, centerY, radius);

  useEffect(() => {
    let cancelled = false;
    loadWorldMapPixels(worldId, worldRevision)
      .then((map) => {
        if (cancelled) return;
        const canvas = document.createElement("canvas");
        canvas.width = LOCATOR_SIZE;
        canvas.height = LOCATOR_SIZE;
        const context = canvas.getContext("2d");
        if (!context) return;
        const pixels = context.createImageData(canvas.width, canvas.height);
        paintGlobe(map, locatorBasis(planeAngleDeg, centerAngle), pixels, centerX, centerY, radius);
        context.putImageData(pixels, 0, 0);
        setTextureUrl(canvas.toDataURL("image/png"));
      })
      .catch(() => {
        if (!cancelled) setTextureUrl(null);
      });
    return () => { cancelled = true; };
  }, [worldId, worldRevision, planeAngleDeg, centerAngle, centerX, centerY]);

  const tiltAt = (clientY: number, startY: number, unitsPerPixel: number) =>
    clampTilt(committedTilt + ((clientY - startY) * unitsPerPixel * 180) / (Math.PI * radius));

  const onPointerDown = (event: React.PointerEvent<SVGSVGElement>) => {
    if (event.button !== 0) return;
    event.currentTarget.setPointerCapture(event.pointerId);
    const rect = event.currentTarget.getBoundingClientRect();
    drag.current = { pointerId: event.pointerId, startY: event.clientY, unitsPerPixel: LOCATOR_SIZE / Math.max(1, rect.height) };
  };

  // Dragging down brings the north into view, as when grabbing a real globe.
  const onPointerMove = (event: React.PointerEvent<SVGSVGElement>) => {
    const current = drag.current;
    if (!current || current.pointerId !== event.pointerId) return;
    setPreviewTilt(tiltAt(event.clientY, current.startY, current.unitsPerPixel));
  };

  const endDrag = (event: React.PointerEvent<SVGSVGElement>) => {
    const current = drag.current;
    if (!current || current.pointerId !== event.pointerId) return;
    drag.current = null;
    setPreviewTilt(null);
    const tilt = tiltAt(event.clientY, current.startY, current.unitsPerPixel);
    if (event.type === "pointerup" && tilt !== committedTilt) onPlaneAngleChange(tilt);
  };

  const onKeyDown = (event: React.KeyboardEvent<SVGSVGElement>) => {
    // Up looks further north, as on the sidebar's globe selector (and as dragging down does).
    const delta = { ArrowUp: TILT_KEYBOARD_STEP_DEG, ArrowDown: -TILT_KEYBOARD_STEP_DEG }[event.key];
    if (delta === undefined) return;
    event.preventDefault();
    const tilt = clampTilt(committedTilt + delta);
    if (tilt !== committedTilt) onPlaneAngleChange(tilt);
  };

  return (
    <svg
      role="slider"
      tabIndex={0}
      className="globe-control"
      aria-label={`Globe locator, arc center at position ${arcPosition + 1} of ${ARC_POSITION_COUNT}. Drag up or down, or use the arrow keys, to tilt the great circle`}
      aria-orientation="vertical"
      aria-valuemin={-TILT_LIMIT_DEG}
      aria-valuemax={TILT_LIMIT_DEG}
      aria-valuenow={planeAngleDeg}
      aria-valuetext={`Great-circle tilt ${planeAngleDeg} degrees`}
      viewBox={`0 0 ${LOCATOR_SIZE} ${LOCATOR_SIZE}`}
      onPointerDown={onPointerDown}
      onPointerMove={onPointerMove}
      onPointerUp={endDrag}
      onPointerCancel={endDrag}
      onKeyDown={onKeyDown}
      style={{
        width: "100%", maxWidth: LOCATOR_SIZE, height: "auto", display: "block", alignSelf: "center",
        cursor: previewTilt === null ? "grab" : "grabbing", touchAction: "none",
      }}
    >
      <defs>
        <clipPath id="cross-section-globe-clip">
          <circle cx={centerX} cy={centerY} r={radius} />
        </clipPath>
      </defs>
      <circle cx={centerX} cy={centerY} r={radius} fill="#122239" fillOpacity="0.82" />
      {textureUrl && <image href={textureUrl} x="0" y="0" width={LOCATOR_SIZE} height={LOCATOR_SIZE} clipPath="url(#cross-section-globe-clip)" />}
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
        <circle cx={arcCenter.x} cy={arcCenter.y} r="4" fill="#61ded0" stroke="#253147" strokeWidth="1.5" />
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
  worldId,
  width,
  height,
  worldRevision,
  planeAngleDeg,
  arcPosition,
  showBoundaries,
  showFaults,
  showPlumes,
  onLayerToggle,
  onPlaneAngleChange,
  onArcPositionChange,
}: Props) {
  const centerDegrees = (360 * arcPosition) / ARC_POSITION_COUNT;
  const [zoom, setZoom] = useState(1);
  const [previewTilt, setPreviewTilt] = useState<number | null>(null);
  const [displayedFrame, setDisplayedFrame] = useState<{ key: string; dataUrl: string } | null>(null);
  const [loading, setLoading] = useState(true);
  const [prefetchProgress, setPrefetchProgress] = useState<number | null>(null);
  const [loadMessage, setLoadMessage] = useState<string | null>(null);
  const scopeControllerRef = useRef<AbortController | null>(null);
  const arcPositionRef = useRef(arcPosition);
  arcPositionRef.current = arcPosition;
  const currentFrameKey: FrameKey = {
    worldId, worldRevision, planeAngleDeg, arcPosition, width: Math.round(width), height: Math.round(height), showBoundaries, showFaults, showPlumes,
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

  const frameScopeKey = [worldId, worldRevision, planeAngleDeg, Math.round(width), Math.round(height), Number(showBoundaries), Number(showFaults), Number(showPlumes), SECTION_STYLE_VERSION].join(":");

  useEffect(() => {
    const controller = new AbortController();
    scopeControllerRef.current = controller;
    const initialPosition = arcPositionRef.current;
    const activeKey: FrameKey = {
      worldId, worldRevision, planeAngleDeg, arcPosition: initialPosition, width: Math.round(width), height: Math.round(height), showBoundaries, showFaults, showPlumes,
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
              // Aborting also ends the selected arc's own load if it shares this controller, so
              // clear its hourglass here; the message says why (#497).
              setLoading(false);
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
  }, [worldId, worldRevision, planeAngleDeg, width, height, showBoundaries, showFaults, showPlumes, frameScopeKey]);

  useEffect(() => {
    const activeKey: FrameKey = {
      worldId, worldRevision, planeAngleDeg, arcPosition, width: Math.round(width), height: Math.round(height), showBoundaries, showFaults, showPlumes,
    };
    rememberActiveFrame(activeKey);
    const cached = cachedFrame(activeKey);
    setDisplayedFrame(cached ? { key: currentCacheKey, dataUrl: cached } : null);
    setLoading(cached === null);
    setLoadMessage(null);
    let isCurrentSelection = true;
    // Once a prefetch error has aborted the scope's controller, loads on it end silently: give
    // this arc its own, or its hourglass would never clear (#497).
    const scopeController = scopeControllerRef.current;
    const controller = scopeController && !scopeController.signal.aborted ? scopeController : new AbortController();
    const ownsController = controller !== scopeController;

    if (!cached) {
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
      if (ownsController) controller.abort();
    };
  }, [worldId, worldRevision, planeAngleDeg, arcPosition, width, height, showBoundaries, showFaults, showPlumes, currentCacheKey]);

  useEffect(() => () => clearMantleCrossSectionCache(), []);

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
        <label style={{ display: "flex", alignItems: "center", gap: 8, flexShrink: 0, color: "#aeb8ca", fontSize: 10 }}>
          <span>Surface zoom · 10 km relief</span>
          <input type="range" min={1} max={10} step={0.1} value={zoom} onChange={(event) => setZoom(Number(event.target.value))} aria-label="Cross-section zoom for mountain relief up to 10 kilometres above sea level" style={{ flex: 1, accentColor: "#61ded0" }} />
          <span style={{ width: 30, textAlign: "right", fontVariantNumeric: "tabular-nums" }}>{zoom.toFixed(1)}×</span>
        </label>
        <div style={{ display: "flex", position: "relative", minHeight: 0, flex: 1, width: "100%", overflow: "hidden" }}>
          <div aria-hidden="true" style={{ display: "flex", minHeight: 0, flex: 1, width: "100%" }}>
            {!imageDataUrl && <SectionPlaceholder width={width} height={height} />}
          </div>
          {/* The selected arc's section is still being computed (#497). */}
          {!imageDataUrl && loading && (
            <div
              role="img"
              aria-label="Computing the lithosphere section"
              title="Computing the lithosphere section"
              style={{
                position: "absolute",
                inset: 0,
                display: "grid",
                placeItems: "center",
                background: "rgba(14, 20, 34, 0.45)",
                fontSize: 44,
                pointerEvents: "none",
              }}
            >
              <span className="hourglass-turning" aria-hidden="true">⏳</span>
            </div>
          )}
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
                transform: `scale(${zoom})`,
                transformOrigin: "50% 8%",
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
          {showBoundaries && <LegendKey color="#58BEFF" label="Plate boundary" line />}
          {showFaults && <LegendKey color="#FF7D57" label="Fault line" line />}
          {showPlumes && <LegendKey color="#F67A3C" label="Mantle plume" />}
          <LegendKey color="#ED7E3C" label="Active vent" />
          <LegendKey color="#B7BEC1" label="Dormant vent" />
        </div>
      </section>

      <aside style={{ minWidth: 0, minHeight: 0, display: "flex", flexDirection: "column", gap: 7, borderLeft: "1px solid #303a4d", paddingLeft: 12 }}>
        <h3 style={{ margin: "1px 0 0", fontSize: 13, fontWeight: 600 }}>Globe locator</h3>
        <GlobeLocator
          worldId={worldId}
          planeAngleDeg={planeAngleDeg}
          arcPosition={arcPosition}
          worldRevision={worldRevision}
          onPlaneAngleChange={onPlaneAngleChange}
          previewTilt={previewTilt}
          onPreviewTiltChange={setPreviewTilt}
        />
        <div style={{ color: "#aeb8ca", fontSize: 10, textAlign: "center", marginTop: -3 }}>
          <div style={{ color: "#91a0b8" }}>Drag the globe up or down to tilt the great circle</div>
          <div>
            Great-circle tilt{" "}
            <span style={{ color: "#b8c3d3", fontVariantNumeric: "tabular-nums" }}>{previewTilt ?? planeAngleDeg}°</span>
            {" · "}
            <span style={{ color: "#61ded0" }}>3,000 km cut</span>
          </div>
        </div>

        <div style={{ marginTop: 2 }}>
          <div>
            <label style={{ display: "block", fontSize: 10 }}>
              <span style={{ display: "flex", justifyContent: "space-between", gap: 4 }}>
                <span>Arc center</span>
                <span style={{ color: "#b8c3d3", fontVariantNumeric: "tabular-nums" }}>{arcPosition + 1}/{ARC_POSITION_COUNT}</span>
              </span>
              <input
                type="range"
                min={0}
                max={ARC_POSITION_COUNT - 1}
                step={1}
                value={arcPosition}
                onChange={(event) => onArcPositionChange(Number(event.target.value))}
                onKeyDown={(event) => {
                  if (event.key === "ArrowRight" || event.key === "ArrowLeft") {
                    event.preventDefault();
                    onArcPositionChange(advanceCrossSectionPosition(
                      arcPosition,
                      ARC_POSITION_COUNT,
                      event.key === "ArrowRight" ? 1 : -1,
                      ARC_KEYBOARD_STEP,
                    ));
                  }
                }}
                aria-label="Position of the 3,000 kilometre arc around the great circle"
                style={{ width: "100%", margin: "3px 0 0", accentColor: "#61ded0" }}
              />
            </label>
            <div style={{ color: "#7f8ba0", fontSize: 8, textAlign: "right" }}>
              {centerDegrees.toFixed(0)}° of 360°
            </div>
          </div>
        </div>

        <fieldset style={{ border: 0, padding: 0, margin: "1px 0 0", display: "grid", gap: 4, fontSize: 10 }}>
          <legend style={{ color: "#aeb8ca", padding: 0, marginBottom: 2 }}>Cross-section overlays</legend>
          <label style={{ display: "flex", gap: 6, alignItems: "center" }}><input type="checkbox" checked={showBoundaries} onChange={(event) => onLayerToggle("boundaries", event.target.checked)} /> Plate boundaries</label>
          <label style={{ display: "flex", gap: 6, alignItems: "center" }}><input type="checkbox" checked={showFaults} onChange={(event) => onLayerToggle("faults", event.target.checked)} /> Fault lines</label>
          <label style={{ display: "flex", gap: 6, alignItems: "center" }}><input type="checkbox" checked={showPlumes} onChange={(event) => onLayerToggle("plumes", event.target.checked)} /> Mantle plumes</label>
        </fieldset>

      </aside>
    </div>
  );
}
