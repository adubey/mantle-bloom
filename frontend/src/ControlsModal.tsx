import { useState } from "react";
import type { TuningKey, TuningMultipliers } from "./api";

interface Props {
  seaLevelM: number;
  solarMultiplier: number;
  simulatePlateMovement: boolean;
  simulateClimateBiomes: boolean;
  windModel: string;
  tuning: TuningMultipliers;
  onSeaLevelChange: (v: number) => void;
  onSolarMultiplierChange: (v: number) => void;
  onSimulatePlateMovementChange: (v: boolean) => void;
  onSimulateClimateBiomesChange: (v: boolean) => void;
  onWindModelChange: (v: string) => void;
  onTuningChange: (key: TuningKey, v: number) => void;
  onTuningReset: () => void;
  onClose: () => void;
}

// The geomorphic-budget tuning knobs, grouped for the panel. Every one is a dimensionless
// multiplier on backend app/world.py's matching World.*_multiplier (1.0 == untuned) -- see
// api.ts's TUNING_MULTIPLIER_KEYS and the backend field group. Order/labels here are the
// only place the UI copy lives.
const TUNING_GROUPS: { heading: string; knobs: { key: TuningKey; label: string }[] }[] = [
  {
    heading: "Erosion",
    knobs: [
      { key: "rain_erosion_multiplier", label: "Rain / sheet erosion" },
      { key: "river_erosion_multiplier", label: "River erosion" },
      { key: "wind_erosion_multiplier", label: "Wind erosion (weathering)" },
      { key: "ocean_erosion_multiplier", label: "Ocean erosion (submarine + coastal)" },
      { key: "coastal_leveling_multiplier", label: "Coastal planation" },
      { key: "glacier_erosion_multiplier", label: "Glacier erosion" },
      { key: "seismic_erosion_multiplier", label: "Seismic erosion" },
    ],
  },
  {
    heading: "Deposition",
    knobs: [
      { key: "river_deposition_multiplier", label: "River deposition (floodplains / deltas)" },
      { key: "ocean_deposition_multiplier", label: "Ocean deposition (shelf / beaches)" },
    ],
  },
  {
    heading: "Tectonics & volcanism",
    knobs: [
      { key: "collision_uplift_multiplier", label: "Collision uplift — amount" },
      { key: "collision_uplift_reach_multiplier", label: "Collision uplift — reach" },
      { key: "volcanism_multiplier", label: "Volcanism" },
    ],
  },
];

// Real-time controls for World's live-adjustable properties (see backend app/world.py's
// World.sea_level_m/World.solar_multiplier/World.simulate_plate_movement/
// World.simulate_climate_biomes/World.wind_model and its tuning-knob field group, plus
// main.py's /world/controls) -- unlike the "Generate World" dialog's sliders, these apply to
// the *current* world immediately, no regenerate needed. Controlled from App.tsx's own state
// (not local state here) since the same values also need to reset on a fresh Generate.
export default function ControlsModal({
  seaLevelM,
  solarMultiplier,
  simulatePlateMovement,
  simulateClimateBiomes,
  windModel,
  tuning,
  onSeaLevelChange,
  onSolarMultiplierChange,
  onSimulatePlateMovementChange,
  onSimulateClimateBiomesChange,
  onWindModelChange,
  onTuningChange,
  onTuningReset,
  onClose,
}: Props) {
  const [showTuning, setShowTuning] = useState(false);
  const anyTuned = Object.values(tuning).some((v) => v !== 1);

  return (
    <div
      onClick={onClose}
      style={{
        position: "fixed",
        inset: 0,
        background: "rgba(0, 0, 0, 0.6)",
        zIndex: 100,
        display: "flex",
        alignItems: "center",
        justifyContent: "center",
      }}
    >
      <div
        onClick={(e) => e.stopPropagation()}
        style={{
          width: 340,
          maxHeight: "85vh",
          overflowY: "auto",
          padding: 20,
          background: "#151a2e",
          border: "1px solid #333",
          borderRadius: 8,
          color: "#e6e8ef",
          fontSize: 13,
          boxShadow: "0 8px 32px rgba(0, 0, 0, 0.6)",
        }}
      >
        <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 14 }}>
          <span style={{ fontSize: 16, fontWeight: 600 }}>Controls</span>
          <button
            type="button"
            onClick={onClose}
            style={{ background: "none", border: "none", color: "#8b8fa3", cursor: "pointer", fontSize: 16 }}
          >
            ✕
          </button>
        </div>

        <label style={{ display: "block", marginBottom: 16 }}>
          Sea level: {seaLevelM} m
          <input
            type="range"
            min={-3000}
            max={3000}
            step={50}
            value={seaLevelM}
            onChange={(e) => onSeaLevelChange(Number(e.target.value))}
            style={{ width: "100%" }}
          />
        </label>

        <label style={{ display: "block", marginBottom: 16 }}>
          Solar heat: {solarMultiplier.toFixed(2)}×
          <input
            type="range"
            min={0.5}
            max={1.5}
            step={0.05}
            value={solarMultiplier}
            onChange={(e) => onSolarMultiplierChange(Number(e.target.value))}
            style={{ width: "100%" }}
          />
        </label>

        <div style={{ borderTop: "1px solid #333", paddingTop: 14 }}>
          <label style={{ display: "flex", alignItems: "center", gap: 8, marginBottom: 10, cursor: "pointer" }}>
            <input
              type="checkbox"
              checked={simulatePlateMovement}
              onChange={(e) => onSimulatePlateMovementChange(e.target.checked)}
            />
            Simulate plate movement
          </label>
          <label style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer" }}>
            <input
              type="checkbox"
              checked={simulateClimateBiomes}
              onChange={(e) => onSimulateClimateBiomesChange(e.target.checked)}
            />
            Simulate climate & biomes
          </label>
          {!simulatePlateMovement && !simulateClimateBiomes && (
            <div style={{ fontSize: 11, color: "#999", marginTop: 8 }}>
              Both are off -- stepping will only advance elapsed years, nothing else changes.
            </div>
          )}
        </div>

        <div style={{ borderTop: "1px solid #333", paddingTop: 14, marginTop: 14 }}>
          <label style={{ display: "block", marginBottom: 6 }}>Wind model</label>
          <select
            value={windModel}
            onChange={(e) => onWindModelChange(e.target.value)}
            style={{
              width: "100%",
              background: "#0f1424",
              color: "#e6e8ef",
              border: "1px solid #333",
              borderRadius: 4,
              padding: "4px 6px",
              fontSize: 13,
            }}
          >
            <option value="cfd">Shallow-water CFD (accurate)</option>
            <option value="diagnostic">Diagnostic / ABL (fast)</option>
          </select>
          <div style={{ fontSize: 11, color: "#999", marginTop: 8 }}>
            {windModel === "diagnostic"
              ? "Skips the shallow-water solve -- much faster steps, ~85-90% of the CFD biome map. Wind/temperature maps show the closed-form field."
              : "Genuine time-integrated shallow-water solve -- the largest single cost of a step at high fluid-dynamics resolution."}
          </div>
        </div>

        <div style={{ borderTop: "1px solid #333", paddingTop: 14, marginTop: 14 }}>
          <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center" }}>
            <button
              type="button"
              onClick={() => setShowTuning((v) => !v)}
              style={{
                background: "none",
                border: "none",
                color: "#e6e8ef",
                cursor: "pointer",
                fontSize: 13,
                fontWeight: 600,
                padding: 0,
              }}
            >
              {showTuning ? "▾" : "▸"} Tuning {anyTuned && !showTuning ? "•" : ""}
            </button>
            {anyTuned && (
              <button
                type="button"
                onClick={onTuningReset}
                style={{
                  background: "none",
                  border: "1px solid #333",
                  color: "#8b8fa3",
                  cursor: "pointer",
                  fontSize: 11,
                  borderRadius: 4,
                  padding: "2px 6px",
                }}
              >
                Reset to 1×
              </button>
            )}
          </div>

          {showTuning && (
            <div style={{ marginTop: 10 }}>
              <div style={{ fontSize: 11, color: "#999", marginBottom: 10 }}>
                Dimensionless multipliers on individual geomorphic processes -- 1× is the
                model's default. Use these to rebalance the long-run land budget on the
                current world; changes apply immediately.
              </div>
              {TUNING_GROUPS.map((group) => (
                <div key={group.heading} style={{ marginBottom: 12 }}>
                  <div style={{ fontSize: 11, textTransform: "uppercase", letterSpacing: 0.5, color: "#8b8fa3", marginBottom: 6 }}>
                    {group.heading}
                  </div>
                  {group.knobs.map(({ key, label }) => (
                    <label key={key} style={{ display: "block", marginBottom: 10, fontSize: 12 }}>
                      <span style={{ display: "flex", justifyContent: "space-between" }}>
                        <span>{label}</span>
                        <span style={{ color: tuning[key] === 1 ? "#8b8fa3" : "#e6e8ef" }}>{tuning[key].toFixed(1)}×</span>
                      </span>
                      <input
                        type="range"
                        min={0}
                        max={3}
                        step={0.1}
                        value={tuning[key]}
                        onChange={(e) => onTuningChange(key, Number(e.target.value))}
                        style={{ width: "100%" }}
                      />
                    </label>
                  ))}
                </div>
              ))}
            </div>
          )}
        </div>
      </div>
    </div>
  );
}
