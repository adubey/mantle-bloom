import ProgressBar from "./ProgressBar";
import { phaseLinks } from "./phaseInfo";

interface Props {
  // What kind of run this is: "Generating world", "Advancing world", "Recording animation".
  headline: string;
  // The backend's current phase label (or, for recordings, "Frame i of N"); null before
  // the run's first checkpoint arrives.
  phase: string | null;
  // Real percent-complete when the run reports one; undefined for indeterminate runs
  // (recordings with no fixed frame count, the debug-generate path).
  fraction?: number;
  onCollapse: () => void;
}

// The expanded run display: a floating card over the map with a spinner, the current
// phase, and learn links for it (#629). Non-modal on purpose -- it never takes a scrim,
// so Play's Stop button stays clickable underneath. The X drops back to just the
// sidebar's pulsing line, whose pop-out icon re-opens this.
export default function RunOverlay({ headline, phase, fraction, onCollapse }: Props) {
  const links = phaseLinks(phase);
  return (
    <div
      role="status"
      aria-live="polite"
      aria-label={phase ? `${headline}: ${phase}` : headline}
      style={{
        position: "fixed",
        left: "50%",
        top: "38%",
        transform: "translate(-50%, -50%)",
        zIndex: 90,
        width: 340,
        maxWidth: "calc(100vw - 32px)",
        background: "#151a2e",
        border: "1px solid #333",
        borderRadius: 8,
        boxShadow: "0 8px 32px rgba(0, 0, 0, 0.6)",
        color: "#e6e8ef",
        fontSize: 13,
        padding: 16,
      }}
    >
      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 10 }}>
        <span style={{ fontSize: 14, fontWeight: 600 }}>{headline}</span>
        <button
          type="button"
          title="Collapse to progress line"
          aria-label="Collapse to progress line"
          onClick={onCollapse}
          style={{ background: "none", border: "none", color: "#8b8fa3", cursor: "pointer", fontSize: 15 }}
        >
          ✕
        </button>
      </div>
      <div style={{ display: "flex", alignItems: "center", gap: 10, marginBottom: 10 }}>
        <span className="run-spinner" aria-hidden="true" />
        <span>{phase ?? "Working…"}</span>
      </div>
      {fraction !== undefined && (
        <div style={{ marginBottom: links.length > 0 ? 10 : 0 }}>
          <ProgressBar label={`${headline}${phase ? `: ${phase}` : ""}`} fraction={fraction} />
        </div>
      )}
      {links.length > 0 && (
        <div style={{ borderTop: "1px solid #333", paddingTop: 8, fontSize: 12 }}>
          <span style={{ opacity: 0.7 }}>Learn more: </span>
          {links.map((link, i) => (
            <span key={link.url}>
              {i > 0 && " · "}
              <a href={link.url} target="_blank" rel="noreferrer">
                {link.label}
              </a>
            </span>
          ))}
        </div>
      )}
    </div>
  );
}
