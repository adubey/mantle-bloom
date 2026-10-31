import { useCallback, useState } from "react";
import type { WorldRecord } from "./api";
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
  // Both reject on failure, so the error shows here rather than behind the modal.
  onRename: (id: WorldId, name: string) => Promise<void>;
  onDelete: (id: WorldId) => Promise<void>;
  onClose: () => void;
}

export default function ManageWorldsModal({ worlds, selectedWorldId, onRename, onDelete, onClose }: Props) {
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
          <div key={world.world_id} style={{ display: "flex", alignItems: "center", gap: 6, marginBottom: 6 }}>
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
        );
      })}
    </ModalShell>
  );
}
