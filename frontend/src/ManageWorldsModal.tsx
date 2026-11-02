import { useCallback, useState } from "react";
import type { ViewerView, WorldRecord } from "./api";
import { worldLabel } from "./api";
import ModalShell from "./ModalShell";
import { useModalAction } from "./useModalAction";
import type { WorldId } from "./worldRoutes";

// Matching backend app/main.py's PatchWorldRequest.name bound.
const MAX_NAME_LENGTH = 100;

interface RenameFormProps {
  world: WorldRecord;
  busy: boolean;
  onSave: (name: string) => void;
  onCancel: () => void;
}

// The draft lives here, so the World list's 15 s poll can't overwrite what is
// being typed, and a failed save leaves it in place to retry.
function RenameForm({ world, busy, onSave, onCancel }: RenameFormProps) {
  const [draft, setDraft] = useState(world.name);
  const name = draft.trim();
  return (
    <form
      onSubmit={(e) => {
        e.preventDefault();
        if (name && !busy) onSave(name);
      }}
      style={{ display: "flex", flex: 1, gap: 6, minWidth: 0 }}
    >
      <input
        autoFocus
        value={draft}
        maxLength={MAX_NAME_LENGTH}
        onChange={(e) => setDraft(e.target.value)}
        onKeyDown={(e) => {
          if (e.key === "Escape" && !busy) onCancel();
        }}
        aria-label={`New name for ${world.name}`}
        style={{ flex: 1, minWidth: 0, fontSize: 12 }}
      />
      <button type="submit" disabled={!name || busy} style={{ fontSize: 12 }}>
        Save
      </button>
      <button type="button" onClick={onCancel} disabled={busy} style={{ fontSize: 12 }}>
        Cancel
      </button>
    </form>
  );
}

interface Props {
  worlds: WorldRecord[];
  selectedWorldId: WorldId | null;
  // All reject on failure, so the error shows here rather than behind the modal.
  onRename: (id: WorldId, name: string) => Promise<void>;
  onDelete: (id: WorldId) => Promise<void>;
  // The /view default and its shared presentation (#755).
  onSetDefault: (id: WorldId, isDefault: boolean) => Promise<void>;
  onViewerViewChange: (id: WorldId, view: ViewerView) => Promise<void>;
  onShareOrientation: (id: WorldId) => Promise<void>;
  onClose: () => void;
}

const VIEWER_VIEWS: { value: ViewerView; label: string }[] = [
  { value: "combined", label: "Combined" },
  { value: "elevation", label: "Elevation" },
  { value: "biome", label: "Biome" },
];

export default function ManageWorldsModal({
  worlds, selectedWorldId, onRename, onDelete, onSetDefault, onViewerViewChange, onShareOrientation, onClose,
}: Props) {
  const { busy, error, runAction } = useModalAction();
  const [editingId, setEditingId] = useState<WorldId | null>(null);
  // A world deleted elsewhere mid-edit drops out of the list without ending the edit.
  const editing = worlds.some((world) => world.world_id === editingId);

  const handleRename = useCallback((world: WorldRecord, name: string) => runAction("rename", async () => {
    if (name !== world.name) await onRename(world.world_id, name);
    setEditingId(null);
  }), [runAction, onRename]);

  const handleDelete = useCallback((world: WorldRecord) => {
    if (!window.confirm(`Delete "${world.name}" and its checkpoint? This cannot be undone.`)) return;
    void runAction(`delete:${world.world_id}`, () => onDelete(world.world_id));
  }, [runAction, onDelete]);

  const handleSetDefault = useCallback((world: WorldRecord) => {
    void runAction(`default:${world.world_id}`, () => onSetDefault(world.world_id, !world.is_default));
  }, [runAction, onSetDefault]);

  const handleViewerView = useCallback((world: WorldRecord, view: ViewerView) => {
    void runAction(`viewer-view:${world.world_id}`, () => onViewerViewChange(world.world_id, view));
  }, [runAction, onViewerViewChange]);

  const handleShareOrientation = useCallback((world: WorldRecord) => {
    void runAction(`share-orientation:${world.world_id}`, () => onShareOrientation(world.world_id));
  }, [runAction, onShareOrientation]);

  // Held open while a request is out: its failure can only be reported here.
  const handleClose = useCallback(() => {
    if (busy === null) onClose();
  }, [busy, onClose]);

  return (
    <ModalShell title="Manage Worlds" onClose={handleClose} width={460} maxWidth="90vw">
      {error && <div style={{ color: "#ff8080", fontSize: 11, marginBottom: 10 }}>{error}</div>}
      {worlds.length === 0 && <div style={{ opacity: 0.7 }}>No worlds have been created yet.</div>}
      {worlds.map((world) => {
        const running = world.run?.running ?? false;
        return (
          <div key={world.world_id} style={{ marginBottom: 10 }}>
            <div style={{ display: "flex", alignItems: "center", gap: 6 }}>
              {world.world_id === editingId ? (
                <RenameForm
                  world={world}
                  busy={busy !== null}
                  onSave={(name) => handleRename(world, name)}
                  onCancel={() => setEditingId(null)}
                />
              ) : (
                <>
                  <span
                    title={world.name}
                    style={{
                      flex: 1, minWidth: 0, overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap",
                      fontWeight: world.world_id === selectedWorldId ? 600 : undefined,
                    }}
                  >
                    {worldLabel(world)}
                  </span>
                  <button
                    onClick={() => setEditingId(world.world_id)}
                    disabled={busy !== null || editing}
                    style={{ fontSize: 12 }}
                  >
                    Rename
                  </button>
                  {/* The polled list can lag a run started elsewhere; the server refuses those too. */}
                  <button
                    onClick={() => handleDelete(world)}
                    disabled={busy !== null || running}
                    title={running ? "Stop this world's run first" : undefined}
                    style={{ fontSize: 12 }}
                  >
                    {busy === `delete:${world.world_id}` ? "Deleting..." : "Delete"}
                  </button>
                </>
              )}
            </div>
            {/* The /view default (#755): one world, shown to every visitor
                with the presentation chosen here. Setting it publishes it. */}
            <div style={{ display: "flex", alignItems: "center", gap: 6, marginTop: 4, fontSize: 12 }}>
              {world.is_default ? (
                <>
                  <span style={{ fontWeight: 600 }}>Default on /view</span>
                  <button onClick={() => handleSetDefault(world)} disabled={busy !== null} style={{ fontSize: 12 }}>
                    {busy === `default:${world.world_id}` ? "Clearing..." : "Clear default"}
                  </button>
                  <label>
                    View{" "}
                    <select
                      value={world.viewer_view}
                      onChange={(e) => handleViewerView(world, e.target.value as ViewerView)}
                      disabled={busy !== null}
                      aria-label={`Map view viewers of ${world.name} see`}
                      style={{ fontSize: 12 }}
                    >
                      {VIEWER_VIEWS.map((view) => (
                        <option key={view.value} value={view.value}>{view.label}</option>
                      ))}
                    </select>
                  </label>
                  <button
                    onClick={() => handleShareOrientation(world)}
                    disabled={busy !== null}
                    title="Show viewers this tab's current map orientation"
                    style={{ fontSize: 12 }}
                  >
                    {busy === `share-orientation:${world.world_id}` ? "Sharing..." : "Use my orientation"}
                  </button>
                </>
              ) : (
                <button
                  onClick={() => handleSetDefault(world)}
                  disabled={busy !== null}
                  title="Open /view on this world and keep its image ready; this also makes it public"
                  style={{ fontSize: 12 }}
                >
                  {busy === `default:${world.world_id}` ? "Setting..." : "Set as /view default"}
                </button>
              )}
            </div>
          </div>
        );
      })}
    </ModalShell>
  );
}
