import { useState } from "react";
import type { MapView } from "./api";

// How many simulation steps each animation frame advances the world -- the real years each
// frame covers is this times the app's current "Years per step" (see App.tsx's
// STEP_YEARS_OPTIONS), so the animation and Play move the world in the same unit. A few sane
// presets, not a free-form input, same reasoning App.tsx's own STEP_YEARS_OPTIONS uses. 1
// (one step per frame) is the default so a fresh animation reads as a smooth progression
// rather than jumping many steps at a time.
const STEPS_PER_FRAME_OPTIONS = [1, 10, 100];
const DEFAULT_STEPS_PER_FRAME = 1;
const DEFAULT_NUM_FRAMES = 20;
// Matching backend app/main.py's MAX_ANIMATION_FRAMES -- also the safety ceiling underneath
// "keep going until Stop is pressed" (see the `unbounded` checkbox below).
const MAX_NUM_FRAMES = 480;

const INSPECTOR_VIEWS: MapView[] = ["plateInspector", "riverInspector", "lakeInspector", "platesAndFaults"];

interface Props {
  hasWorld: boolean;
  // The app's current "Years per step" (see App.tsx's STEP_YEARS_OPTIONS) -- one animation
  // frame advances the world by `stepsPerFrame * stepYears` real years.
  stepYears: number;
  mapView: MapView;
  onClose: () => void;
  // The "Record" toolbar button's whole point -- closes this dialog and starts recording in
  // the background off the current map view (see App.tsx's handleStartAnimation).
  // `unbounded` is "keep going until Stop is pressed": `numFrames` is still a real, finite
  // safety ceiling underneath it (the caller is never expected to reach it), it's just not
  // shown to the user as the run's actual endpoint.
  onStartAnimation: (opts: { numFrames: number; yearsPerFrame: number; unbounded: boolean }) => void;
}

function fmtMyr(years: number): string {
  // Steps per frame can be as fine as one 10,000-year step, so a plain "0.0 Myr" is possible
  // -- fall back to kyr below a million years so the number stays readable.
  if (years < 1e6) return `${Math.round(years / 1e3).toLocaleString()} kyr`;
  return `${(years / 1e6).toFixed(1)} Myr`;
}

export default function AnimationModal({ hasWorld, stepYears, mapView, onClose, onStartAnimation }: Props) {
  const [numFrames, setNumFrames] = useState(DEFAULT_NUM_FRAMES);
  const [stepsPerFrame, setStepsPerFrame] = useState(DEFAULT_STEPS_PER_FRAME);
  // Defaults to on -- most recordings are "run until I've seen enough", not a pre-planned
  // frame count, so the common case shouldn't require touching the checkbox.
  const [unbounded, setUnbounded] = useState(true);
  const yearsPerFrame = stepsPerFrame * stepYears;

  const isInspectorView = INSPECTOR_VIEWS.includes(mapView);
  const canStart = hasWorld && !isInspectorView;

  return (
    <div
      onClick={onClose}
      style={{
        position: "fixed", inset: 0, background: "rgba(0, 0, 0, 0.6)", zIndex: 100,
        display: "flex", alignItems: "center", justifyContent: "center",
      }}
    >
      <div
        onClick={(e) => e.stopPropagation()}
        style={{
          width: 380, maxHeight: "85vh", overflowY: "auto", padding: 20, background: "#151a2e",
          border: "1px solid #333", borderRadius: 8, color: "#e6e8ef", fontSize: 13,
          boxShadow: "0 8px 32px rgba(0, 0, 0, 0.6)",
        }}
      >
        <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 14 }}>
          <span style={{ fontSize: 16, fontWeight: 600 }}>⏺ Record Animation</span>
          <button type="button" onClick={onClose} style={{ background: "none", border: "none", color: "#8b8fa3", cursor: "pointer", fontSize: 16 }}>
            ✕
          </button>
        </div>

        {isInspectorView ? (
          <div style={{ fontSize: 11, opacity: 0.6, marginBottom: 12 }}>
            Switch to a map view (not an Inspector) to make an animation.
          </div>
        ) : (
          <>
            <label style={{ display: "block", marginBottom: 6 }}>
              Frames (incl. current): {unbounded ? "—" : numFrames}
              <input
                type="range" min={2} max={MAX_NUM_FRAMES} value={numFrames}
                onChange={(e) => setNumFrames(Number(e.target.value))}
                disabled={unbounded}
                style={{ width: "100%" }}
              />
            </label>
            <label style={{ display: "flex", alignItems: "center", gap: 6, marginBottom: 10, cursor: "pointer" }}>
              <input type="checkbox" checked={unbounded} onChange={(e) => setUnbounded(e.target.checked)} />
              Keep going until Stop is pressed
            </label>
            <label style={{ display: "block", marginBottom: 6 }}>
              Steps per frame
              <select
                value={stepsPerFrame}
                onChange={(e) => setStepsPerFrame(Number(e.target.value))}
                style={{ width: "100%", fontSize: 12 }}
              >
                {STEPS_PER_FRAME_OPTIONS.map((s) => (
                  <option key={s} value={s}>{s.toLocaleString()}</option>
                ))}
              </select>
            </label>
            <div style={{ fontSize: 11, opacity: 0.6, marginBottom: 10 }}>
              {fmtMyr(yearsPerFrame)} per frame ({stepsPerFrame.toLocaleString()} × {stepYears.toLocaleString()} yr/step).{" "}
              {unbounded
                ? `Permanently advances the world for as long as it runs (stops on its own after ${MAX_NUM_FRAMES.toLocaleString()} frames if Stop is never pressed).`
                : `Permanently advances the world by ${fmtMyr((numFrames - 1) * yearsPerFrame)}, same as clicking Step ${((numFrames - 1) * stepsPerFrame).toLocaleString()} times -- not a preview.`}
            </div>
            <div style={{ fontSize: 11, opacity: 0.75, marginBottom: 12 }}>
              Recording renders in the background as the world steps forward, using the
              current map view and orientation. The main map doubles as the live preview while
              it runs; other controls are paused until you press Stop (⏹) or it finishes on
              its own, then a dialog lets you save or discard the MP4.
            </div>
            <button
              onClick={() => onStartAnimation({ numFrames, yearsPerFrame, unbounded })}
              disabled={!canStart}
              style={{ width: "100%", fontSize: 12 }}
            >
              ⏺ Start Recording
            </button>
          </>
        )}
      </div>
    </div>
  );
}
