import type { StopMode } from "./api";
import ModalShell from "./ModalShell";
import {
  GLACIER_STEP_YEARS_MAX,
  GLACIER_STEP_YEARS_MIN,
  LAKE_STEP_YEARS_MAX,
  LAKE_STEP_YEARS_MIN,
  PLATE_STEP_YEARS_TICKS,
  plateYearsToTickIndex,
} from "./timeSettings";

interface Props {
  stepYears: number;
  lakeStepYears: number;
  lakeStopMode: StopMode;
  lakeMaxSubsteps: number;
  glacierStepYears: number;
  glacierStopMode: StopMode;
  glacierMaxSubsteps: number;
  simulateClimateBiomes: boolean;
  simulateGlacierMarch: boolean;
  onStepYearsChange: (years: number) => void;
  onLakeStepYearsChange: (years: number) => void;
  onLakeStopModeChange: (mode: StopMode) => void;
  onLakeMaxSubstepsChange: (max: number) => void;
  onGlacierStepYearsChange: (years: number) => void;
  onGlacierStopModeChange: (mode: StopMode) => void;
  onGlacierMaxSubstepsChange: (max: number) => void;
  onSimulateGlacierMarchChange: (enabled: boolean) => void;
  onClose: () => void;
}

const STOP_MODES: { value: StopMode; label: string }[] = [
  { value: "exhaust", label: "Exhaust" },
  { value: "resolve", label: "Resolve" },
  { value: "max", label: "Max" },
];

const COLUMN_HEADING_STYLE = {
  fontSize: 11,
  textTransform: "uppercase",
  letterSpacing: 0.5,
  color: "#8b8fa3",
  marginBottom: 12,
} as const;

function clampMaxSubsteps(raw: number, maxBound: number): number {
  if (!Number.isFinite(raw)) return 1;
  return Math.min(maxBound, Math.max(1, Math.round(raw)));
}

function stopModeNote(kind: "lake" | "glacier", mode: StopMode, maxSubsteps: number): string {
  if (kind === "lake") {
    if (mode === "resolve") return "Stop slicing once lakes settle and lump the remaining time into one final slice.";
    if (mode === "max") {
      return `Stop once lakes settle, or after ${maxSubsteps} normal slices and lump the remaining time into one final slice.`;
    }
    return "Run every lake slice across the full tectonic step.";
  }
  if (mode === "resolve") return "Stop marching once ice volume settles and hold that state for the remainder of the step.";
  if (mode === "max") {
    return `Stop once ice volume settles, or after ${maxSubsteps} substeps and hold that state for the remainder of the step.`;
  }
  return "March ice sheets through the full tectonic step without early hold.";
}

function StopModeToggle({
  label,
  value,
  maxSubsteps,
  maxBound,
  disabled = false,
  note,
  onChange,
  onMaxSubstepsChange,
}: {
  label: string;
  value: StopMode;
  maxSubsteps: number;
  maxBound: number;
  disabled?: boolean;
  note: string;
  onChange: (mode: StopMode) => void;
  onMaxSubstepsChange: (max: number) => void;
}) {
  const clampedMax = clampMaxSubsteps(maxSubsteps, maxBound);
  return (
    <div style={{ marginBottom: 16, opacity: disabled ? 0.5 : 1 }}>
      <div style={{ marginBottom: 6, fontSize: 12 }}>{label}</div>
      <div style={{ display: "flex", gap: 4 }}>
        {STOP_MODES.map((opt) => {
          const active = value === opt.value;
          return (
            <button
              key={opt.value}
              type="button"
              disabled={disabled}
              aria-pressed={active}
              onClick={() => onChange(opt.value)}
              style={{
                flex: 1,
                background: active ? "#26304d" : "transparent",
                border: "1px solid",
                borderColor: active ? "#6b7cea" : "#333",
                color: active ? "#e6e8ef" : "#8b8fa3",
                cursor: disabled ? "default" : "pointer",
                fontSize: 12,
                fontWeight: active ? 600 : 400,
                borderRadius: 4,
                padding: "5px 0",
              }}
            >
              {opt.label}
            </button>
          );
        })}
      </div>
      {value === "max" && (
        <label style={{ display: "flex", alignItems: "center", justifyContent: "space-between", gap: 8, marginTop: 8, fontSize: 12 }}>
          <span>Max steps (1–{maxBound.toLocaleString()})</span>
          <input
            type="number"
            min={1}
            max={maxBound}
            step={1}
            value={clampedMax}
            disabled={disabled}
            onChange={(e) => {
              const parsed = Number(e.target.value);
              if (Number.isFinite(parsed)) onMaxSubstepsChange(clampMaxSubsteps(parsed, maxBound));
            }}
            style={{
              width: 72,
              background: "#181a24",
              border: "1px solid #333",
              borderRadius: 4,
              color: "#e6e8ef",
              fontSize: 12,
              padding: "3px 6px",
              textAlign: "right",
            }}
          />
        </label>
      )}
      <div style={{ fontSize: 11, color: "#999", marginTop: 6 }}>{note}</div>
    </div>
  );
}

export default function TimeSettingsModal({
  stepYears,
  lakeStepYears,
  lakeStopMode,
  lakeMaxSubsteps,
  glacierStepYears,
  glacierStopMode,
  glacierMaxSubsteps,
  simulateClimateBiomes,
  simulateGlacierMarch,
  onStepYearsChange,
  onLakeStepYearsChange,
  onLakeStopModeChange,
  onLakeMaxSubstepsChange,
  onGlacierStepYearsChange,
  onGlacierStopModeChange,
  onGlacierMaxSubstepsChange,
  onSimulateGlacierMarchChange,
  onClose,
}: Props) {
  const glacierActive = simulateClimateBiomes && simulateGlacierMarch;
  const lakeMaxBound = Math.max(1, Math.ceil(stepYears / lakeStepYears));
  const glacierMaxBound = Math.max(1, Math.ceil(stepYears / glacierStepYears));

  return (
    <ModalShell title="Time" onClose={onClose} width={520}>
      <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 20 }}>
        <div>
          <div style={COLUMN_HEADING_STYLE}>Step size</div>

          <label style={{ display: "block", marginBottom: 16, fontSize: 12, opacity: simulateClimateBiomes ? 1 : 0.5 }}>
            <span style={{ display: "flex", justifyContent: "space-between", marginBottom: 2 }}>
              <span>Lake</span>
              <span>{lakeStepYears.toLocaleString()} yr</span>
            </span>
            <input
              type="range"
              min={LAKE_STEP_YEARS_MIN}
              max={LAKE_STEP_YEARS_MAX}
              step={10}
              value={lakeStepYears}
              disabled={!simulateClimateBiomes}
              onChange={(e) => onLakeStepYearsChange(Number(e.target.value))}
              style={{ width: "100%" }}
            />
          </label>

          <label style={{ display: "block", marginBottom: 16, fontSize: 12, opacity: glacierActive ? 1 : 0.5 }}>
            <span style={{ display: "flex", justifyContent: "space-between", marginBottom: 2 }}>
              <span>Glacier</span>
              <span>{glacierStepYears.toLocaleString()} yr</span>
            </span>
            <input
              type="range"
              min={GLACIER_STEP_YEARS_MIN}
              max={GLACIER_STEP_YEARS_MAX}
              step={50}
              value={glacierStepYears}
              disabled={!glacierActive}
              onChange={(e) => onGlacierStepYearsChange(Number(e.target.value))}
              style={{ width: "100%" }}
            />
          </label>

          <label style={{ display: "block", marginBottom: 4, fontSize: 12 }}>
            <span style={{ display: "flex", justifyContent: "space-between", marginBottom: 2 }}>
              <span>Plates</span>
              <span>{stepYears.toLocaleString()} yr</span>
            </span>
            <input
              type="range"
              min={0}
              max={PLATE_STEP_YEARS_TICKS.length - 1}
              step={1}
              value={plateYearsToTickIndex(stepYears)}
              onChange={(e) => onStepYearsChange(PLATE_STEP_YEARS_TICKS[Number(e.target.value)])}
              style={{ width: "100%" }}
            />
          </label>
          <div style={{ display: "flex", justifyContent: "space-between", fontSize: 10, color: "#8b8fa3" }}>
            <span>1,000</span>
            <span>10,000</span>
            <span>100,000</span>
          </div>
        </div>

        <div>
          <div style={COLUMN_HEADING_STYLE}>Stopping conditions</div>

          <StopModeToggle
            label="Lake"
            value={lakeStopMode}
            maxSubsteps={lakeMaxSubsteps}
            maxBound={lakeMaxBound}
            disabled={!simulateClimateBiomes}
            note={stopModeNote("lake", lakeStopMode, clampMaxSubsteps(lakeMaxSubsteps, lakeMaxBound))}
            onChange={onLakeStopModeChange}
            onMaxSubstepsChange={onLakeMaxSubstepsChange}
          />

          <StopModeToggle
            label="Glacier"
            value={glacierStopMode}
            maxSubsteps={glacierMaxSubsteps}
            maxBound={glacierMaxBound}
            disabled={!glacierActive}
            note={stopModeNote("glacier", glacierStopMode, clampMaxSubsteps(glacierMaxSubsteps, glacierMaxBound))}
            onChange={onGlacierStopModeChange}
            onMaxSubstepsChange={onGlacierMaxSubstepsChange}
          />
        </div>
      </div>

      <div style={{ borderTop: "1px solid #333", paddingTop: 14, marginTop: 14, opacity: simulateClimateBiomes ? 1 : 0.5 }}>
        <label style={{ display: "flex", alignItems: "center", gap: 8, cursor: simulateClimateBiomes ? "pointer" : "default" }}>
          <input
            type="checkbox"
            checked={simulateGlacierMarch}
            disabled={!simulateClimateBiomes}
            onChange={(e) => onSimulateGlacierMarchChange(e.target.checked)}
          />
          March glaciers to end-of-step
        </label>
        <div style={{ fontSize: 11, color: "#999", marginTop: 6 }}>
          {simulateGlacierMarch
            ? `Each step marches ice sheets in ${glacierStepYears.toLocaleString()}-year substeps so ice, meltwater, sea level and glacial erosion fully respond.`
            : "Off: ice is solved once per step and frozen there — much faster steps. Turn on for the full ice-sheet response."}
        </div>
      </div>
    </ModalShell>
  );
}
