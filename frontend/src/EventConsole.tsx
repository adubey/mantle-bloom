import { useState } from "react";
import type { WorldEvent } from "./api";

interface Props {
  events: WorldEvent[];
}

function EventList({ events, maxHeight }: { events: WorldEvent[]; maxHeight: number | string }) {
  return (
    <div style={{ maxHeight, overflowY: "auto", fontSize: 12, fontFamily: "ui-monospace, monospace" }}>
      {events.length === 0 && <div style={{ opacity: 0.5 }}>No events yet.</div>}
      {events
        .slice()
        .reverse()
        .map((e, i) => (
          <div key={i} style={{ marginBottom: 4 }}>
            <span style={{ opacity: 0.5 }}>[{(e.elapsed_years / 1e6).toFixed(1)} Myr]</span> {e.message}
          </div>
        ))}
    </div>
  );
}

export default function EventConsole({ events }: Props) {
  const [collapsed, setCollapsed] = useState(false);
  // The popped-out view is a separate, larger in-page overlay reusing the same fixed-scrim
  // pattern as the Generate World / Stats dialogs -- not a real separate OS window (this app
  // has never used window.open, so a true detached window would be a new, untested
  // interaction rather than a reuse of an established one).
  const [poppedOut, setPoppedOut] = useState(false);

  return (
    <>
      <div style={{ border: "1px solid #333", borderRadius: 6, background: "#0d1120" }}>
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
            <span>Console ({events.length})</span>
            <span>{collapsed ? "▸" : "▾"}</span>
          </button>
          <button
            type="button"
            title="Pop out console"
            aria-label="Pop out console"
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
            <EventList events={events} maxHeight={220} />
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
              <span style={{ fontSize: 16, fontWeight: 600 }}>Console ({events.length})</span>
              <button
                type="button"
                onClick={() => setPoppedOut(false)}
                style={{ background: "none", border: "none", color: "#8b8fa3", cursor: "pointer", fontSize: 16 }}
              >
                ✕
              </button>
            </div>
            <EventList events={events} maxHeight="65vh" />
          </div>
        </div>
      )}
    </>
  );
}
