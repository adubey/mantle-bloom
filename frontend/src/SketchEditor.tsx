import { useEffect, useRef, useState } from "react";

// The full-screen "Draw a map" editor opened from the Generate World dialog's "Human-made"
// tab (see App.tsx). Draws directly onto one <canvas> at a fixed 2:1 (equirectangular) aspect
// ratio -- the backend (app/worldsketch.py) resamples onto its own working grid regardless of
// exact pixel size, but keeping this fixed means the drawing always covers the whole globe the
// same way "Load an image" is documented to assume. Internal resolution is higher than the
// on-screen display size purely for smoother strokes; CSS scales it down.
const INTERNAL_WIDTH = 1440;
const INTERNAL_HEIGHT = 720;
const DISPLAY_WIDTH = 1000;
const DISPLAY_HEIGHT = 500;

// Must match backend app/worldsketch.py's COAST_COLOR/RIVER_COLOR/MOUNTAIN_COLOR exactly --
// those constants exist there so a script/test can draw a synthetic sketch guaranteed to
// classify correctly; this is the other half of that same convention. "erase" isn't really an
// ink color, it just repaints the background over whatever's there.
const TOOL_COLORS: Record<Tool, string> = {
  coast: "#000000",
  river: "#1e6fd9",
  mountain: "#8b4513",
  erase: "#ffffff",
};
const BACKGROUND_COLOR = "#ffffff";
const GRATICULE_COLOR = "#d8d8d8";

type Tool = "coast" | "erase" | "river" | "mountain";

const TOOL_CHOICES: { value: Tool; label: string; hint: string }[] = [
  { value: "coast", label: "✏️ Coast", hint: "Draw a rough coastline" },
  { value: "river", label: "🏞️ River", hint: "Optional -- a hint for where water should flow" },
  { value: "mountain", label: "⛰️ Mountain", hint: "Optional -- adds extra uplift there" },
  { value: "erase", label: "🧽 Erase", hint: "Erase any stroke" },
];

const DEFAULT_BRUSH_SIZE = 6;

interface Props {
  // A previous drawing to resume editing (re-opening "Draw a map" after already drawing one),
  // as a full data URL (`data:image/png;base64,...`) -- the same shape state gets stored in.
  initialImageDataUrl?: string | null;
  // Called with the finished drawing as a full data URL when "Done" is clicked.
  onDone: (dataUrl: string) => void;
  onCancel: () => void;
}

function drawGraticule(ctx: CanvasRenderingContext2D) {
  ctx.fillStyle = BACKGROUND_COLOR;
  ctx.fillRect(0, 0, INTERNAL_WIDTH, INTERNAL_HEIGHT);
  ctx.strokeStyle = GRATICULE_COLOR;
  ctx.lineWidth = 1;
  // Equator + prime meridian only -- just enough to orient a rough sketch, not a real grid.
  ctx.beginPath();
  ctx.moveTo(0, INTERNAL_HEIGHT / 2);
  ctx.lineTo(INTERNAL_WIDTH, INTERNAL_HEIGHT / 2);
  ctx.moveTo(INTERNAL_WIDTH / 2, 0);
  ctx.lineTo(INTERNAL_WIDTH / 2, INTERNAL_HEIGHT);
  ctx.stroke();
}

export default function SketchEditor({ initialImageDataUrl, onDone, onCancel }: Props) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [tool, setTool] = useState<Tool>("coast");
  const [brushSize, setBrushSize] = useState(DEFAULT_BRUSH_SIZE);
  // Whether anything has actually been drawn yet (beyond the blank graticule) -- Done is
  // disabled until there's at least one coastline stroke, since an empty sketch would
  // generate as an all-ocean world with no useful "human-made" intent behind it.
  const [hasCoastStroke, setHasCoastStroke] = useState(false);
  const drawingRef = useRef(false);
  const lastPointRef = useRef<{ x: number; y: number } | null>(null);

  // Initial paint -- the graticule, or the previous drawing being resumed, scaled to fill the
  // canvas (matches how the backend treats any input image: stretched onto its own working
  // grid regardless of source aspect ratio).
  useEffect(() => {
    const canvas = canvasRef.current;
    const ctx = canvas?.getContext("2d");
    if (!canvas || !ctx) return;
    drawGraticule(ctx);
    if (initialImageDataUrl) {
      const img = new Image();
      img.onload = () => ctx.drawImage(img, 0, 0, INTERNAL_WIDTH, INTERNAL_HEIGHT);
      img.src = initialImageDataUrl;
      setHasCoastStroke(true);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const canvasPoint = (e: React.PointerEvent<HTMLCanvasElement>): { x: number; y: number } => {
    const rect = e.currentTarget.getBoundingClientRect();
    return {
      x: ((e.clientX - rect.left) / rect.width) * INTERNAL_WIDTH,
      y: ((e.clientY - rect.top) / rect.height) * INTERNAL_HEIGHT,
    };
  };

  const handlePointerDown = (e: React.PointerEvent<HTMLCanvasElement>) => {
    e.currentTarget.setPointerCapture(e.pointerId);
    drawingRef.current = true;
    const point = canvasPoint(e);
    lastPointRef.current = point;
    // A single click/tap with no drag should still leave a dot, not nothing.
    const ctx = canvasRef.current?.getContext("2d");
    if (!ctx) return;
    ctx.fillStyle = TOOL_COLORS[tool];
    ctx.beginPath();
    ctx.arc(point.x, point.y, brushSize / 2, 0, Math.PI * 2);
    ctx.fill();
    if (tool === "coast") setHasCoastStroke(true);
  };

  const handlePointerMove = (e: React.PointerEvent<HTMLCanvasElement>) => {
    if (!drawingRef.current) return;
    const ctx = canvasRef.current?.getContext("2d");
    const from = lastPointRef.current;
    if (!ctx || !from) return;
    const to = canvasPoint(e);
    ctx.strokeStyle = TOOL_COLORS[tool];
    ctx.lineWidth = brushSize;
    ctx.lineCap = "round";
    ctx.lineJoin = "round";
    ctx.beginPath();
    ctx.moveTo(from.x, from.y);
    ctx.lineTo(to.x, to.y);
    ctx.stroke();
    lastPointRef.current = to;
    if (tool === "coast") setHasCoastStroke(true);
  };

  const stopDrawing = () => {
    drawingRef.current = false;
    lastPointRef.current = null;
  };

  const handleClear = () => {
    const ctx = canvasRef.current?.getContext("2d");
    if (!ctx) return;
    drawGraticule(ctx);
    setHasCoastStroke(false);
  };

  const handleDone = () => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    onDone(canvas.toDataURL("image/png"));
  };

  return (
    <div
      style={{
        position: "fixed",
        inset: 0,
        background: "rgba(0, 0, 0, 0.85)",
        zIndex: 300,
        display: "flex",
        flexDirection: "column",
        alignItems: "center",
        justifyContent: "center",
        color: "#e6e8ef",
        fontFamily: "system-ui, sans-serif",
      }}
    >
      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", width: DISPLAY_WIDTH, marginBottom: 10 }}>
        <span style={{ fontSize: 16, fontWeight: 600 }}>Draw a map</span>
        <button type="button" onClick={onCancel} style={{ background: "none", border: "none", color: "#8b8fa3", cursor: "pointer", fontSize: 16 }}>
          ✕
        </button>
      </div>

      <div style={{ display: "flex", gap: 8, alignItems: "center", width: DISPLAY_WIDTH, marginBottom: 10, fontSize: 12 }}>
        {TOOL_CHOICES.map((t) => (
          <button
            key={t.value}
            type="button"
            title={t.hint}
            onClick={() => setTool(t.value)}
            style={{
              padding: "6px 10px",
              borderRadius: 6,
              border: tool === t.value ? "1px solid #5b8cff" : "1px solid #333",
              background: tool === t.value ? "#232a4a" : "#151a2e",
              color: "#e6e8ef",
              cursor: "pointer",
            }}
          >
            {t.label}
          </button>
        ))}
        <label style={{ display: "flex", alignItems: "center", gap: 6, marginLeft: 8 }}>
          Brush
          <input
            type="range"
            min={2}
            max={20}
            value={brushSize}
            onChange={(e) => setBrushSize(Number(e.target.value))}
            style={{ width: 100 }}
          />
        </label>
        <button type="button" onClick={handleClear} style={{ marginLeft: "auto", fontSize: 12 }}>
          Clear
        </button>
      </div>

      <canvas
        ref={canvasRef}
        width={INTERNAL_WIDTH}
        height={INTERNAL_HEIGHT}
        style={{
          width: DISPLAY_WIDTH,
          height: DISPLAY_HEIGHT,
          borderRadius: 6,
          border: "1px solid #333",
          touchAction: "none",
          cursor: "crosshair",
        }}
        onPointerDown={handlePointerDown}
        onPointerMove={handlePointerMove}
        onPointerUp={stopDrawing}
        onPointerLeave={stopDrawing}
        onPointerCancel={stopDrawing}
      />

      <p style={{ fontSize: 11, opacity: 0.7, width: DISPLAY_WIDTH, marginTop: 8, marginBottom: 0 }}>
        Draw a rough coastline (corners of the map are assumed ocean). Rivers and mountains are
        optional hints -- a river carves a shallow channel for water to find once the world is
        stepped; a mountain stroke adds extra uplift there.
      </p>

      <div style={{ display: "flex", gap: 8, width: DISPLAY_WIDTH, justifyContent: "flex-end", marginTop: 14 }}>
        <button type="button" onClick={onCancel}>
          Cancel
        </button>
        <button type="button" onClick={handleDone} disabled={!hasCoastStroke}>
          Done
        </button>
      </div>
    </div>
  );
}
