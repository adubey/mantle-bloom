import { useState } from "react";
import type { CornerNotchLogEntry } from "./api";

interface Props {
  entries: CornerNotchLogEntry[];
  enabled: boolean;
}

// One line per entry, most-recent first -- same shape as EventConsole's own EventList, but
// each row carries the extra structured fields (outcome, nodes_added, window geometry) the
// corner-notch log needs instead of a single free-text message.
function EntryList({ entries, maxHeight }: { entries: CornerNotchLogEntry[]; maxHeight: number | string }) {
  return (
    <div style={{ maxHeight, overflowY: "auto", fontSize: 12, fontFamily: "ui-monospace, monospace" }}>
      {entries.length === 0 && <div style={{ opacity: 0.5 }}>No corner-notch activity yet.</div>}
      {entries
        .slice()
        .reverse()
        .map((e, i) => (
          <div key={i} style={{ marginBottom: 4 }}>
            <span style={{ opacity: 0.5 }}>[{(e.elapsed_years / 1e6).toFixed(1)} Myr]</span>{" "}
            <span style={{ opacity: 0.7 }}>plate {e.plate_id}</span> —{" "}
            <span style={{ color: e.outcome === "claimed" ? "#7fd88f" : "#e6e8ef" }}>{e.outcome}</span>
            {e.nodes_added > 0 && <span style={{ opacity: 0.7 }}> ({e.nodes_added} nodes)</span>}
          </div>
        ))}
    </div>
  );
}

// Debug-only, structured decision log for LithospherePlate._fill_corner_notch_frontier (see
// GET /world/corner_notch_log / World.debug_diagnostics) -- deliberately its own panel, never
// mixed into EventConsole's always-on Event Console (see docs/debugging.md: this can fire
// once per plate per step, far higher volume than that log is meant to carry). Same
// collapse/pop-out shape as EventConsole for a consistent feel.
export default function CornerNotchLogPanel({ entries, enabled }: Props) {
  const [collapsed, setCollapsed] = useState(true);
  const [poppedOut, setPoppedOut] = useState(false);

  return (
    <>
      <div style={{ border: "1px solid #333", borderRadius: 6, background: "#0d1120", opacity: enabled ? 1 : 0.6 }}>
        <div style={{ display: "flex", alignItems: "center" }}>
          <button
            type="button"
            onClick={() => setCollapsed((c) => !c)}
            aria-expanded={!collapsed}
            style={{
              flex: 1,
              textAlign: "left",
              background: "transparent",
              border: "none",
              color: "#e6e8ef",
              padding: "8px 10px",
              fontSize: 12,
              textTransform: "uppercase",
              letterSpacing: "0.06em",
              opacity: 0.8,
              cursor: "pointer",
              display: "flex",
              justifyContent: "space-between",
            }}
          >
            <span>Corner-notch log ({entries.length})</span>
            <span>{collapsed ? "▸" : "▾"}</span>
          </button>
          <button
            type="button"
            title="Pop out corner-notch log"
            aria-label="Pop out corner-notch log"
            onClick={() => setPoppedOut(true)}
            style={{
              background: "transparent",
              border: "none",
              color: "#e6e8ef",
              opacity: 0.7,
              cursor: "pointer",
              fontSize: 13,
              padding: "8px 10px",
            }}
          >
            ⤢
          </button>
        </div>
        {!collapsed && (
          <div style={{ padding: "0 10px 10px" }}>
            {!enabled && (
              <div style={{ fontSize: 11, color: "#999", marginBottom: 6 }}>
                Off — enable "Log corner-notch decisions" in Controls → Tectonics to start recording.
              </div>
            )}
            <EntryList entries={entries} maxHeight={220} />
          </div>
        )}
      </div>

      {poppedOut && (
        <div
          onClick={() => setPoppedOut(false)}
          style={{
            position: "fixed",
            inset: 0,
            background: "rgba(0, 0, 0, 0.6)",
            zIndex: 100,
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
          }}
        >
          <div
            onClick={(e) => e.stopPropagation()}
            style={{
              width: 560,
              maxHeight: "80vh",
              padding: 20,
              background: "#151a2e",
              border: "1px solid #333",
              borderRadius: 8,
              boxShadow: "0 8px 32px rgba(0, 0, 0, 0.6)",
              display: "flex",
              flexDirection: "column",
            }}
          >
            <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 12 }}>
              <span style={{ fontSize: 16, fontWeight: 600 }}>Corner-notch log ({entries.length})</span>
              <button
                type="button"
                onClick={() => setPoppedOut(false)}
                style={{ background: "none", border: "none", color: "#8b8fa3", cursor: "pointer", fontSize: 16 }}
              >
                ✕
              </button>
            </div>
            <EntryList entries={entries} maxHeight="65vh" />
          </div>
        </div>
      )}
    </>
  );
}
