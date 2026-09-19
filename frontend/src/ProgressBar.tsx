import type { CSSProperties } from "react";

interface Props {
  // Also doubles as the accessible name (aria-label) -- there's no visible text label, just
  // the animated bar itself, so screen readers need this to announce what's in progress.
  label: string;
  style?: CSSProperties;
}

// An indeterminate "something's happening" bar (see index.css's .progress-track/
// .progress-indeterminate) -- shown while a request that reports no real incremental progress
// (world generation, a step) is in flight. Not a percentage: there's nothing real to measure
// (see App.tsx's busy/stepping usages and issue #195, tracking real streaming progress later).
export default function ProgressBar({ label, style }: Props) {
  return (
    <div className="progress-track" aria-label={label} role="progressbar" style={style}>
      <div className="progress-indeterminate" />
    </div>
  );
}
