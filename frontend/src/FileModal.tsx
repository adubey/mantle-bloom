import { useCallback, useRef, useState } from "react";
import type { AnimateResponse, MapView, Projection, WorldSummary } from "./api";
import { animateWorld, exportHexGrid, loadWorld, saveWorld } from "./api";
import type { Mat3 } from "./rotation";

// Matching backend app/main.py's STEP_YEARS_OPTIONS-style framing for /world/step -- the
// same "a few sane presets, not a free-form input" reasoning App.tsx's own
// STEP_YEARS_OPTIONS already uses, here for the animation's own "years per frame" choice.
// 100,000 (the finest option) is the default so a fresh animation reads as a smooth
// progression rather than jumping in million-year leaps.
const YEARS_PER_FRAME_OPTIONS = [100_000, 1_000_000, 10_000_000];
const DEFAULT_YEARS_PER_FRAME = 100_000;
const DEFAULT_NUM_FRAMES = 20;
// Matching backend app/main.py's MAX_ANIMATION_FRAMES.
const MAX_NUM_FRAMES = 60;

// Matching backend app/geodesic.py's FREQUENCY_CHOICES/tile_count (10*frequency**2 + 2).
const HEX_FREQUENCY_OPTIONS: { frequency: number; label: string }[] = [
  { frequency: 8, label: "Low (~642 tiles)" },
  { frequency: 16, label: "Medium (~2,562 tiles)" },
  { frequency: 32, label: "High (~10,242 tiles)" },
];
const DEFAULT_HEX_FREQUENCY = 16;

const INSPECTOR_VIEWS: MapView[] = ["plateInspector", "riverInspector", "lakeInspector"];

interface Props {
  hasWorld: boolean;
  seed: number | null;
  elapsedYears: number | null;
  projection: Projection;
  mapView: MapView;
  rotation: Mat3;
  renderWidth: number;
  renderHeight: number;
  // The DOM node wrapping whichever map view component is currently mounted (MapCanvas /
  // PlateInspector / RiverInspector / LakeInspector all draw onto their own <canvas> inside
  // it) -- "Save Image" reads straight from that canvas's own pixels
  // (`canvas.toDataURL`) rather than threading a ref through four separate component prop
  // interfaces, so it works uniformly across every view, not just the PNG-backed ones.
  mapWrapperRef: React.RefObject<HTMLDivElement | null>;
  onClose: () => void;
  // A load or a completed animation both replace/advance the live world -- the caller runs
  // the same full refresh sequence it already runs after generateWorld/stepWorld.
  onWorldReplaced: (summary: WorldSummary) => Promise<void>;
  onWorldAdvanced: (summary: WorldSummary) => Promise<void>;
}

function downloadBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

function fmtMyr(years: number): string {
  return `${(years / 1e6).toFixed(1)} Myr`;
}

export default function FileModal({
  hasWorld, seed, elapsedYears, projection, mapView, rotation, renderWidth, renderHeight, mapWrapperRef,
  onClose, onWorldReplaced, onWorldAdvanced,
}: Props) {
  const [busy, setBusy] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [numFrames, setNumFrames] = useState(DEFAULT_NUM_FRAMES);
  const [yearsPerFrame, setYearsPerFrame] = useState(DEFAULT_YEARS_PER_FRAME);
  const [animation, setAnimation] = useState<AnimateResponse | null>(null);
  const [hexFrequency, setHexFrequency] = useState(DEFAULT_HEX_FREQUENCY);
  const fileInputRef = useRef<HTMLInputElement>(null);

  const isInspectorView = INSPECTOR_VIEWS.includes(mapView);

  const runAction = useCallback(async (name: string, action: () => Promise<void>) => {
    setBusy(name);
    setError(null);
    try {
      await action();
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(null);
    }
  }, []);

  const handleSaveWorld = useCallback(() => runAction("save", async () => {
    const blob = await saveWorld();
    downloadBlob(blob, `mantle-bloom-seed${seed}-${Math.round(elapsedYears ?? 0)}y.mbworld`);
  }), [runAction, seed, elapsedYears]);

  const handleLoadFileChosen = useCallback((file: File) => runAction("load", async () => {
    const summary = await loadWorld(file);
    setAnimation(null);
    await onWorldReplaced(summary);
  }), [runAction, onWorldReplaced]);

  const handleSaveImage = useCallback(() => runAction("image", async () => {
    const canvas = mapWrapperRef.current?.querySelector("canvas");
    if (!canvas) throw new Error("no map canvas to capture");
    const dataUrl = canvas.toDataURL("image/png");
    const blob = await (await fetch(dataUrl)).blob();
    downloadBlob(blob, `mantle-bloom-seed${seed}-${mapView}-${Math.round(elapsedYears ?? 0)}y.png`);
  }), [runAction, mapWrapperRef, seed, mapView, elapsedYears]);

  const handleMakeAnimation = useCallback(() => runAction("animate", async () => {
    const result = await animateWorld(projection, mapView, renderWidth, renderHeight, rotation, yearsPerFrame, numFrames);
    setAnimation(result);
    await onWorldAdvanced(result);
  }), [runAction, projection, mapView, renderWidth, renderHeight, rotation, yearsPerFrame, numFrames, onWorldAdvanced]);

  const handleSaveAnimation = useCallback(() => runAction("saveAnimation", async () => {
    if (!animation) return;
    const blob = await (await fetch(`data:image/gif;base64,${animation.image_base64}`)).blob();
    downloadBlob(blob, `mantle-bloom-seed${seed}-animation-${Math.round(animation.elapsed_years)}y.gif`);
  }), [runAction, animation, seed]);

  const handleExportHexGrid = useCallback(() => runAction("export", async () => {
    const result = await exportHexGrid(hexFrequency);
    const blob = new Blob([JSON.stringify(result)], { type: "application/json" });
    downloadBlob(blob, `mantle-bloom-seed${seed}-hexgrid-f${hexFrequency}.json`);
  }), [runAction, hexFrequency, seed]);

  return (
    <div
      onClick={onClose}
      style={{
        position: "fixed", inset: 0, background: "rgba(0, 0, 0, 0.6)", zIndex: 100,
        display: "flex", alignItems: "center", justifyContent: "center",
      }}
    >
      <div
        onClick={(e) => e.stopPropagation()}
        style={{
          width: 380, maxHeight: "85vh", overflowY: "auto", padding: 20, background: "#151a2e",
          border: "1px solid #333", borderRadius: 8, color: "#e6e8ef", fontSize: 13,
          boxShadow: "0 8px 32px rgba(0, 0, 0, 0.6)",
        }}
      >
        <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 14 }}>
          <span style={{ fontSize: 16, fontWeight: 600 }}>File</span>
          <button type="button" onClick={onClose} style={{ background: "none", border: "none", color: "#8b8fa3", cursor: "pointer", fontSize: 16 }}>
            ✕
          </button>
        </div>

        {error && <div style={{ color: "#ff8080", fontSize: 11, marginBottom: 10 }}>{error}</div>}

        <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, marginBottom: 12 }}>
          <legend style={{ fontSize: 11 }}>World state</legend>
          <div style={{ display: "flex", gap: 6, marginBottom: 6 }}>
            <button onClick={handleSaveWorld} disabled={!hasWorld || busy !== null} style={{ flex: 1, fontSize: 12 }}>
              {busy === "save" ? "Saving..." : "Save World"}
            </button>
            <button onClick={() => fileInputRef.current?.click()} disabled={busy !== null} style={{ flex: 1, fontSize: 12 }}>
              {busy === "load" ? "Loading..." : "Load World"}
            </button>
          </div>
          <input
            ref={fileInputRef}
            type="file"
            accept=".mbworld"
            style={{ display: "none" }}
            onChange={(e) => {
              const file = e.target.files?.[0];
              e.target.value = ""; // so choosing the same file again still fires onChange
              if (file) handleLoadFileChosen(file);
            }}
          />
          <div style={{ fontSize: 11, opacity: 0.6 }}>
            Saves the complete simulation state to a single file, with no guarantee of
            compatibility across app versions. Loading replaces the current world entirely.
          </div>
        </fieldset>

        <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, marginBottom: 12 }}>
          <legend style={{ fontSize: 11 }}>Image</legend>
          <button onClick={handleSaveImage} disabled={!hasWorld || busy !== null} style={{ width: "100%", fontSize: 12 }}>
            {busy === "image" ? "Saving..." : "Save Image of Current View"}
          </button>
        </fieldset>

        <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, marginBottom: 12 }}>
          <legend style={{ fontSize: 11 }}>Animation</legend>
          {isInspectorView ? (
            <div style={{ fontSize: 11, opacity: 0.6 }}>
              Switch to a map view (not an Inspector) to make an animation.
            </div>
          ) : (
            <>
              <label style={{ display: "block", marginBottom: 6 }}>
                Frames (incl. current): {numFrames}
                <input
                  type="range" min={2} max={MAX_NUM_FRAMES} value={numFrames}
                  onChange={(e) => setNumFrames(Number(e.target.value))}
                  style={{ width: "100%" }}
                />
              </label>
              <label style={{ display: "block", marginBottom: 6 }}>
                Years per frame
                <select
                  value={yearsPerFrame}
                  onChange={(e) => setYearsPerFrame(Number(e.target.value))}
                  style={{ width: "100%", fontSize: 12 }}
                >
                  {YEARS_PER_FRAME_OPTIONS.map((y) => (
                    <option key={y} value={y}>{y.toLocaleString()}</option>
                  ))}
                </select>
              </label>
              <div style={{ fontSize: 11, opacity: 0.6, marginBottom: 6 }}>
                Permanently advances the world by {fmtMyr((numFrames - 1) * yearsPerFrame)}, same
                as clicking Step that many times -- not a preview.
              </div>
              <button onClick={handleMakeAnimation} disabled={!hasWorld || busy !== null} style={{ width: "100%", fontSize: 12, marginBottom: 6 }}>
                {busy === "animate" ? "Simulating..." : "Make Animation"}
              </button>
              {animation && (
                <>
                  <img
                    src={`data:image/gif;base64,${animation.image_base64}`}
                    alt="animation preview"
                    style={{ width: "100%", borderRadius: 4, marginBottom: 6 }}
                  />
                  <button onClick={handleSaveAnimation} disabled={busy !== null} style={{ width: "100%", fontSize: 12 }}>
                    {busy === "saveAnimation" ? "Saving..." : "Save Animation (GIF)"}
                  </button>
                </>
              )}
            </>
          )}
        </fieldset>

        <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8 }}>
          <legend style={{ fontSize: 11 }}>Export Hex Grid</legend>
          <label style={{ display: "block", marginBottom: 6 }}>
            Detail
            <select
              value={hexFrequency}
              onChange={(e) => setHexFrequency(Number(e.target.value))}
              style={{ width: "100%", fontSize: 12 }}
            >
              {HEX_FREQUENCY_OPTIONS.map((opt) => (
                <option key={opt.frequency} value={opt.frequency}>{opt.label}</option>
              ))}
            </select>
          </label>
          <div style={{ fontSize: 11, opacity: 0.6, marginBottom: 6 }}>
            Tiles the sphere into a geodesic-dome hex/pentagon grid, mapping current
            elevation and biome onto each tile, as JSON for use in another application (see
            docs/hex-export-format.md).
          </div>
          <button onClick={handleExportHexGrid} disabled={!hasWorld || busy !== null} style={{ width: "100%", fontSize: 12 }}>
            {busy === "export" ? "Exporting..." : "Export Hex Grid"}
          </button>
        </fieldset>
      </div>
    </div>
  );
}
