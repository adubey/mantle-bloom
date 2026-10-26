import type { CSSProperties } from "react";

interface Props {
  // Also doubles as the accessible name (aria-label) -- there's no visible text label, just
  // the animated bar itself, so screen readers need this to announce what's in progress.
  label: string;
  // When the bar fronts an expanded run overlay (see RunOverlay.tsx), the pop-out icon
  // that re-opens it after the overlay's X collapsed it (#629).
  onExpand?: () => void;
  style?: CSSProperties;
  // A real fraction (0 to 1) from the backend's NDJSON progress stream (see api.ts's
  // generateWorld `onProgress` -- world.py's generate_world_progress is where these
  // fractions actually come from). `undefined` (e.g.
  // before the first progress line of a run has arrived) falls back to the old
  // indeterminate "something's happening" animation (see index.css's .progress-track/
  // .progress-indeterminate) rather than showing a bar frozen at 0%.
  fraction?: number;
}

export default function ProgressBar({ label, onExpand, style, fraction }: Props) {
  const determinate = fraction !== undefined;
  const percent = determinate ? Math.round(fraction * 100) : undefined;
  return (
    <div style={{ display: "flex", alignItems: "center", gap: 6 }}>
      <div
        className="progress-track"
        aria-label={label}
        role="progressbar"
        aria-valuenow={percent}
        aria-valuemin={determinate ? 0 : undefined}
        aria-valuemax={determinate ? 100 : undefined}
        style={{ flex: 1, ...style }}
      >
        {determinate ? (
          <div className="progress-fill" style={{ width: `${percent}%` }} />
        ) : (
          <div className="progress-indeterminate" />
        )}
      </div>
      {onExpand && (
        <button
          type="button"
          title="Show run details"
          aria-label="Show run details"
          onClick={onExpand}
          style={{
            background: "transparent",
            border: "none",
            color: "#e6e8ef",
            opacity: 0.7,
            cursor: "pointer",
            fontSize: 13,
            padding: "0 2px",
          }}
        >
          ⤢
        </button>
      )}
    </div>
  );
}
