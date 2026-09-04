import type { AnimateResponse } from "./api";

interface Props {
  result: AnimateResponse;
  onSave: () => void;
  onDiscard: () => void;
}

// Pops up as soon as a recording ends -- whether the toolbar's Stop button cut it short or it
// ran to its full frame count (see App.tsx's handleStartAnimation) -- so saving the MP4 is a
// direct consequence of ending the recording rather than something the user has to go looking
// for in a sidebar panel. No backdrop-dismiss: the recording already permanently advanced the
// world either way (see AnimationModal.tsx), so losing the video by clicking past this needs to
// be a deliberate Discard, not an accidental click outside the box.
export default function SaveAnimationModal({ result, onSave, onDiscard }: Props) {
  const myr = (result.elapsed_years / 1e6).toFixed(1);
  return (
    <div
      style={{
        position: "fixed", inset: 0, background: "rgba(0, 0, 0, 0.6)", zIndex: 100,
        display: "flex", alignItems: "center", justifyContent: "center",
      }}
    >
      <div
        style={{
          width: 320, padding: 20, background: "#151a2e",
          border: "1px solid #333", borderRadius: 8, color: "#e6e8ef", fontSize: 13,
          boxShadow: "0 8px 32px rgba(0, 0, 0, 0.6)",
        }}
      >
        <div style={{ fontSize: 16, fontWeight: 600, marginBottom: 10 }}>⏺ Recording {result.stoppedEarly ? "stopped" : "finished"}</div>
        <div style={{ fontSize: 11, opacity: 0.75, marginBottom: 16 }}>
          The world is now at {myr} Myr -- it advanced permanently while recording, whether or
          not you keep the video. Save it as an MP4, or discard it.
        </div>
        <div style={{ display: "flex", gap: 6 }}>
          <button onClick={onSave} style={{ flex: 1, fontSize: 12 }}>
            Save
          </button>
          <button onClick={onDiscard} style={{ flex: 1, fontSize: 12 }}>
            Discard
          </button>
        </div>
      </div>
    </div>
  );
}
