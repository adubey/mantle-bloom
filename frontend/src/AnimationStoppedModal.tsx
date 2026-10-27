import type { AnimationRecovery } from "./api";
import ModalShell from "./ModalShell";

export interface FailedAnimation {
  // Last completed frame the client saw before the failure.
  frame: number;
  total: number;
  unbounded: boolean;
  error: string;
}

interface Props {
  failure: FailedAnimation;
  // The server's recoverable state, fetched when the failure landed (null when the
  // run never recorded a frame, already finished, or was already recovered).
  recovery: AnimationRecovery | null;
  recovering: boolean;
  onRecover: () => void;
  onClose: () => void;
}

// "Animation stopped" instead of a bare error after a failed recording (#638): the
// finished video died with the run, but the world kept every completed frame's
// simulation work, so Recover records just the remaining frames from the last
// checkpoint rather than redoing hours of compute.
export default function AnimationStoppedModal({ failure, recovery, recovering, onRecover, onClose }: Props) {
  const remaining = recovery ? recovery.total_frames - recovery.completed_frames : 0;
  const canRecover = !recovering && recovery !== null && remaining >= 1;
  return (
    <ModalShell title="Animation stopped" onClose={onClose}>
      <div style={{ fontSize: 12, marginBottom: 10 }}>
        {failure.frame} of {failure.total} frames finished before the run failed. The finished
        video was lost with it, but the world kept every completed frame&apos;s simulation
        work -- {canRecover
          ? `Recover records the remaining ${remaining} frame${remaining === 1 ? "" : "s"} from the last checkpoint.`
          : "There is nothing on the server to resume from."}
      </div>
      <details style={{ fontSize: 11, opacity: 0.85, marginBottom: 12 }}>
        <summary style={{ cursor: "pointer" }}>Error details</summary>
        <pre style={{ whiteSpace: "pre-wrap", wordBreak: "break-word" }}>{failure.error}</pre>
      </details>
      <div style={{ display: "flex", gap: 8 }}>
        <button
          onClick={onRecover}
          disabled={!canRecover}
          title={!recovery ? "The server has no interrupted run to resume" : remaining < 1 ? "No frames left to record" : undefined}
          style={{ flex: 1, fontSize: 12 }}
        >
          {recovering ? "Recovering…" : "Recover…"}
        </button>
        <button onClick={onClose} style={{ flex: 1, fontSize: 12 }}>
          Dismiss
        </button>
      </div>
    </ModalShell>
  );
}
