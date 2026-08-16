import { useState } from "react";
import type { WorldEvent } from "./api";

interface Props {
  events: WorldEvent[];
}

export default function EventConsole({ events }: Props) {
  const [collapsed, setCollapsed] = useState(false);

  return (
    <div style={{ border: "1px solid #333", borderRadius: 6, background: "#0d1120" }}>
      <button
        type="button"
        onClick={() => setCollapsed((c) => !c)}
        aria-expanded={!collapsed}
        style={{
          width: "100%",
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
      {!collapsed && (
        <div
          style={{
            maxHeight: 220,
            overflowY: "auto",
            padding: "0 10px 10px",
            fontSize: 12,
            fontFamily: "ui-monospace, monospace",
          }}
        >
          {events.length === 0 && <div style={{ opacity: 0.5 }}>No events yet.</div>}
          {events
            .slice()
            .reverse()
            .map((e, i) => (
              <div key={i} style={{ marginBottom: 4 }}>
                <span style={{ opacity: 0.5 }}>[{(e.elapsed_years / 1e6).toFixed(1)} Myr]</span>{" "}
                {e.message}
              </div>
            ))}
        </div>
      )}
    </div>
  );
}
