import { useCallback, useRef, useState } from "react";
import type { MapView, WorldSummary } from "./api";
import { exportHexGrid, loadWorld, saveWorld } from "./api";
import ModalShell from "./ModalShell";
import type { WorldId } from "./worldRoutes";

// Matching backend app/geodesic.py's FREQUENCY_CHOICES/tile_count (10*frequency**2 + 2).
const HEX_FREQUENCY_OPTIONS: { frequency: number; label: string }[] = [
  { frequency: 8, label: "Low (~642 tiles)" },
  { frequency: 16, label: "Medium (~2,562 tiles)" },
  { frequency: 32, label: "High (~10,242 tiles)" },
];
const DEFAULT_HEX_FREQUENCY = 16;

interface Props {
  worldId: string;
  hasWorld: boolean;
  seed: number | null;
  elapsedYears: number | null;
  mapView: MapView;
  // The DOM node wrapping whichever map view component is currently mounted (MapCanvas /
  // PlateInspector / RiverInspector / LakeInspector all draw onto their own <canvas> inside
  // it) -- "Save Image" reads straight from that canvas's own pixels
  // (`canvas.toDataURL`) rather than threading a ref through four separate component prop
  // interfaces, so it works uniformly across every view, not just the PNG-backed ones.
  mapWrapperRef: React.RefObject<HTMLDivElement | null>;
  onClose: () => void;
  // Mints the world a load needs when none is selected yet (a fresh server holds
  // no worlds, and Add no longer creates empty ones, #597).
  onEnsureWorldId: () => Promise<WorldId>;
  // A load replaces the live world -- the caller runs the same full refresh sequence it
  // already runs after generateWorld.
  onWorldReplaced: (summary: WorldSummary) => Promise<void>;
}

function downloadBlob(blob: Blob, filename: string): void {
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  a.click();
  URL.revokeObjectURL(url);
}

export default function FileModal({
  worldId, hasWorld, seed, elapsedYears, mapView, mapWrapperRef,
  onClose, onWorldReplaced, onEnsureWorldId,
}: Props) {
  const [busy, setBusy] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [hexFrequency, setHexFrequency] = useState(DEFAULT_HEX_FREQUENCY);
  const fileInputRef = useRef<HTMLInputElement>(null);

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
    const blob = await saveWorld(worldId);
    downloadBlob(blob, `mantle-bloom-seed${seed}-${Math.round(elapsedYears ?? 0)}y.bloom`);
  }), [runAction, seed, elapsedYears, worldId]);

  const handleLoadFileChosen = useCallback((file: File) => runAction("load", async () => {
    const summary = await loadWorld(worldId || (await onEnsureWorldId()), file);
    await onWorldReplaced(summary);
    onClose();
  }), [runAction, onWorldReplaced, onClose, worldId, onEnsureWorldId]);

  const handleSaveImage = useCallback(() => runAction("image", async () => {
    const canvas = mapWrapperRef.current?.querySelector("canvas");
    if (!canvas) throw new Error("no map canvas to capture");
    const dataUrl = canvas.toDataURL("image/png");
    const blob = await (await fetch(dataUrl)).blob();
    downloadBlob(blob, `mantle-bloom-seed${seed}-${mapView}-${Math.round(elapsedYears ?? 0)}y.png`);
  }), [runAction, mapWrapperRef, seed, mapView, elapsedYears]);

  const handleExportHexGrid = useCallback(() => runAction("export", async () => {
    const result = await exportHexGrid(worldId, hexFrequency);
    const blob = new Blob([JSON.stringify(result)], { type: "application/json" });
    downloadBlob(blob, `mantle-bloom-seed${seed}-hexgrid-f${hexFrequency}.json`);
  }), [runAction, hexFrequency, seed, worldId]);

  return (
    <ModalShell title="File" onClose={onClose}>
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
            accept=".bloom,.mbworld"
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
    </ModalShell>
  );
}
