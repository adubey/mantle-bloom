import { useCallback, useRef, useState } from "react";
import type { MapView, WorldImportConflict, WorldRecord } from "./api";
import { exportHexGrid, importWorld, saveWorld, WorldImportConflictError } from "./api";
import ModalShell from "./ModalShell";

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
  // A load adds a brand-new world record -- the caller refreshes the list and selects it.
  onWorldImported: (record: WorldRecord) => Promise<void>;
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
  onClose, onWorldImported,
}: Props) {
  const [busy, setBusy] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  // A file whose import hit a taken name, awaiting the user's pick of the
  // dialog's three options: keep both, rename, or cancel (#652).
  const [conflict, setConflict] = useState<{ file: File; conflict: WorldImportConflict } | null>(null);
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

  const finishImport = useCallback(async (record: WorldRecord) => {
    setConflict(null);
    await onWorldImported(record);
    onClose();
  }, [onWorldImported, onClose]);

  const handleLoadFileChosen = useCallback((file: File) => runAction("load", async () => {
    // Loading adds the save as a new world -- the open planet is never touched,
    // so no selection is needed first (#652, superseding #639's disabled Load).
    try {
      await finishImport(await importWorld(file));
    } catch (e) {
      if (e instanceof WorldImportConflictError) {
        setConflict({ file, conflict: e.conflict });
        return;
      }
      throw e;
    }
  }), [runAction, finishImport]);

  // The duplicate-name dialog's options (a) and (b): keep both under one name,
  // or take the server's unused "<name> (1)" suggestion. Cancel just drops the
  // pending file and leaves this modal open.
  const handleConflictChoice = useCallback((choice: "keep" | "rename") => runAction("load", async () => {
    if (!conflict) return;
    const record = choice === "keep"
      ? await importWorld(conflict.file, { allowDuplicate: true })
      : await importWorld(conflict.file, { name: conflict.conflict.suggestedName });
    await finishImport(record);
  }), [runAction, conflict, finishImport]);

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
          {conflict && (
            <div style={{ border: "1px solid #665", borderRadius: 6, padding: 8, marginBottom: 6, fontSize: 12 }}>
              <div style={{ marginBottom: 4 }}>
                A world named &ldquo;{conflict.conflict.duplicateName}&rdquo; already exists.
                {conflict.conflict.samePlanet && (
                  <>
                    {" "}This looks like the same planet at another date
                    {conflict.conflict.existingElapsedYears !== null && (
                      <> (this file: {Math.round(conflict.conflict.uploadedElapsedYears).toLocaleString()}y, existing: {Math.round(conflict.conflict.existingElapsedYears).toLocaleString()}y)</>
                    )}
                    .
                  </>
                )}
              </div>
              <div style={{ display: "flex", gap: 6 }}>
                <button onClick={() => handleConflictChoice("keep")} disabled={busy !== null} style={{ flex: 1, fontSize: 12 }}>
                  Load anyway
                </button>
                <button onClick={() => handleConflictChoice("rename")} disabled={busy !== null} style={{ flex: 1, fontSize: 12 }}>
                  Load as &ldquo;{conflict.conflict.suggestedName}&rdquo;
                </button>
                <button onClick={() => setConflict(null)} disabled={busy !== null} style={{ flex: 1, fontSize: 12 }}>
                  Cancel
                </button>
              </div>
            </div>
          )}
          <div style={{ fontSize: 11, opacity: 0.6 }}>
            Saves the complete simulation state to a single file, with no guarantee of
            compatibility across app versions. Loading adds the file&apos;s planet as a new
            world and selects it -- the open planet is never touched.
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
