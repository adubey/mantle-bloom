import { useEffect, useRef, useState } from "react";
import type { WorldRecord } from "./api";
import type { WorldId } from "./worldRoutes";

interface Props {
  worlds: WorldRecord[];
  worldId: WorldId | null;
  disabled: boolean;
  disabledTitle?: string;
  onSelect: (id: WorldId | null) => void;
  // Called without any confirm -- the caller confirms first, same as the
  // sidebar Delete button's own window.confirm (#633).
  onDelete: (id: WorldId) => void;
}

// The World fieldset's menu: a <select> can't host per-option buttons, so this
// custom dropdown lists every world with its own delete affordance (#633). The
// empty option from the old <select> survives as the first row, so starting
// with no world selected (#597) still works.
export default function WorldMenu({ worlds, worldId, disabled, disabledTitle, onSelect, onDelete }: Props) {
  const [open, setOpen] = useState(false);
  const rootRef = useRef<HTMLDivElement | null>(null);
  const selected = worlds.find((w) => w.world_id === worldId) ?? null;

  // Close on outside click or Escape; nothing global lingers after unmount.
  useEffect(() => {
    if (!open) return;
    const onPointerDown = (e: MouseEvent) => {
      if (rootRef.current && !rootRef.current.contains(e.target as Node)) setOpen(false);
    };
    const onKeyDown = (e: KeyboardEvent) => {
      if (e.key === "Escape") setOpen(false);
    };
    document.addEventListener("mousedown", onPointerDown);
    document.addEventListener("keydown", onKeyDown);
    return () => {
      document.removeEventListener("mousedown", onPointerDown);
      document.removeEventListener("keydown", onKeyDown);
    };
  }, [open ]);

  const pick = (id: WorldId | null) => {
    setOpen(false);
    onSelect(id);
  };

  return (
    <div ref={rootRef} style={{ position: "relative", width: "100%", marginBottom: 6 }}>
      <button
        onClick={() => setOpen((o) => !o)}
        disabled={disabled}
        title={disabled ? disabledTitle : undefined}
        aria-haspopup="menu"
        aria-expanded={open}
        aria-label="World"
        style={{ width: "100%", fontSize: 12, textAlign: "left" }}
      >
        {selected ? selected.name : "Select a world…"}
      </button>
      {open && !disabled && (
        <div role="menu" aria-label="Worlds" style={{ position: "absolute", zIndex: 50, left: 0, right: 0, background: "#1a1d29", border: "1px solid #555", borderRadius: 6, padding: 4, maxHeight: 220, overflowY: "auto" }}>
          <button role="menuitem" onClick={() => pick(null)} style={{ display: "block", width: "100%", textAlign: "left", fontSize: 12, background: "transparent", border: 0, color: "inherit", padding: "4px 6px" }}>
            Select a world…
          </button>
          {worlds.map((w) => (
            <div key={w.world_id} style={{ display: "flex", alignItems: "center", gap: 4 }}>
              <button
                role="menuitem"
                onClick={() => pick(w.world_id)}
                title={w.world_id === worldId ? "Currently selected" : `Show ${w.name}`}
                style={{ flex: 1, textAlign: "left", fontSize: 12, background: w.world_id === worldId ? "#2a3050" : "transparent", border: 0, color: "inherit", padding: "4px 6px", overflow: "hidden", textOverflow: "ellipsis", whiteSpace: "nowrap" }}
              >
                {w.name}
                {w.summary ? "" : " (empty)"}
                {w.run?.running ? " ▶" : ""}
              </button>
              <button
                onClick={() => { setOpen(false); onDelete(w.world_id); }}
                title={`Delete ${w.name}`}
                aria-label={`Delete ${w.name}`}
                style={{ flex: "none", fontSize: 12, background: "transparent", border: 0, color: "#ff8080", padding: "4px 6px" }}
              >
                ✕
              </button>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
