import type { CSSProperties } from "react";

interface Props {
  // Also doubles as the accessible name (aria-label) -- there's no visible text label, just
  // the animated bar itself, so screen readers need this to announce what's in progress.
  label: string;
  style?: CSSProperties;
  // A real fraction (0 to 1) from the backend's NDJSON progress stream (see api.ts's
  // generateWorld/stepWorld `onProgress` -- world.py's generate_world_progress/
  // step_world_progress are where these fractions actually come from). `undefined` (e.g.
  // before the first progress line of a run has arrived) falls back to the old
  // indeterminate "something's happening" animation (see index.css's .progress-track/
  // .progress-indeterminate) rather than showing a bar frozen at 0%.
  fraction?: number;
}

export default function ProgressBar({ label, style, fraction }: Props) {
  const determinate = fraction !== undefined;
  const percent = determinate ? Math.round(fraction * 100) : undefined;
  return (
    <div
      className="progress-track"
      aria-label={label}
      role="progressbar"
      aria-valuenow={percent}
      aria-valuemin={determinate ? 0 : undefined}
      aria-valuemax={determinate ? 100 : undefined}
      style={style}
    >
      {determinate ? (
        <div className="progress-fill" style={{ width: `${percent}%` }} />
      ) : (
        <div className="progress-indeterminate" />
      )}
    </div>
  );
}
