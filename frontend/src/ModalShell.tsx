import type { ReactNode } from "react";

interface Props {
  title: string;
  onClose: () => void;
  children: ReactNode;
  width?: number;
  maxWidth?: string;
  maxHeight?: string;
  zIndex?: number;
}

export default function ModalShell({
  title, onClose, children, width = 380, maxWidth, maxHeight = "85vh", zIndex = 100,
}: Props) {
  return (
    <div
      onClick={onClose}
      style={{
        position: "fixed", inset: 0, background: "rgba(0, 0, 0, 0.6)", zIndex,
        display: "flex", alignItems: "center", justifyContent: "center",
      }}
    >
      <div
        onClick={(event) => event.stopPropagation()}
        style={{
          width, maxWidth, maxHeight, overflowY: "auto", padding: 20, background: "#151a2e",
          border: "1px solid #333", borderRadius: 8, color: "#e6e8ef", fontSize: 13,
          boxShadow: "0 8px 32px rgba(0, 0, 0, 0.6)",
        }}
      >
        <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 14 }}>
          <span style={{ fontSize: 16, fontWeight: 600 }}>{title}</span>
          <button type="button" onClick={onClose} style={{ background: "none", border: "none", color: "#8b8fa3", cursor: "pointer", fontSize: 16 }}>
            ✕
          </button>
        </div>
        {children}
      </div>
    </div>
  );
}
