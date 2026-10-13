import type { ReactNode } from "react";

interface Props {
  // Where the click landed, in display pixels -- the popup opens just below-right of it,
  // clamped to stay inside the map.
  displayX: number;
  displayY: number;
  title: string; // the clicked point's lat/lon
  status: "loading" | "ok" | "error";
  width: number;
  // The popup's tallest expected height, for the clamp against the map's bottom edge.
  maxHeight: number;
  mapWidth: number;
  mapHeight: number;
  onClose: () => void;
  // The sampled values, shown once `status` is "ok".
  children: ReactNode;
}

// The small inspect popup a map click opens (the point sample, the nodeAge node probe and the
// Voronoi point probe in App.tsx): a close button, the lat/lon, and the sampled values once
// they arrive.
export default function ProbePopup({
  displayX, displayY, title, status, width, maxHeight, mapWidth, mapHeight, onClose, children,
}: Props) {
  return (
    <div
      style={{
        position: "absolute",
        left: Math.max(4, Math.min(displayX + 12, mapWidth - width - 12)),
        top: Math.max(4, Math.min(displayY + 12, mapHeight - maxHeight)),
        width,
        background: "#151a2e",
        border: "1px solid #333",
        borderRadius: 6,
        padding: "8px 10px",
        fontSize: 11,
        lineHeight: 1.6,
        boxShadow: "0 2px 10px rgba(0, 0, 0, 0.5)",
        zIndex: 10,
      }}
    >
      <button
        type="button"
        title="Close"
        onClick={onClose}
        style={{
          position: "absolute", top: 3, right: 4, width: 18, height: 18, padding: 0,
          border: "none", background: "transparent", color: "#999", cursor: "pointer",
          fontSize: 14, lineHeight: "18px",
        }}
      >
        ×
      </button>
      <div style={{ opacity: 0.6, marginBottom: 4 }}>{title}</div>
      {status === "loading" && <div style={{ opacity: 0.7 }}>Sampling…</div>}
      {status === "error" && <div style={{ color: "#ff8080" }}>Couldn’t sample this point.</div>}
      {status === "ok" && children}
    </div>
  );
}
