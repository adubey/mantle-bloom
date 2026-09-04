import { useState } from "react";
import type { TuningKey, TuningMultipliers } from "./api";

interface Props {
  seaLevelM: number;
  solarMultiplier: number;
  iceAgePeriodYears: number;
  simulatePlateMovement: boolean;
  simulateClimateBiomes: boolean;
  windModel: string;
  faultDeformationMode: string;
  tuning: TuningMultipliers;
  onSeaLevelChange: (v: number) => void;
  onSolarMultiplierChange: (v: number) => void;
  onIceAgePeriodChange: (v: number) => void;
  onSimulatePlateMovementChange: (v: boolean) => void;
  onSimulateClimateBiomesChange: (v: boolean) => void;
  onWindModelChange: (v: string) => void;
  onFaultDeformationModeChange: (v: string) => void;
  onTuningChange: (key: TuningKey, v: number) => void;
  onTuningReset: () => void;
  onClose: () => void;
}

type TabKey = "climate" | "erosion" | "tectonics";

// The geomorphic-budget tuning knobs, grouped and assigned to a Controls tab. Every one is a
// dimensionless multiplier on backend app/world.py's matching World.*_multiplier (1.0 ==
// untuned) -- see api.ts's TUNING_MULTIPLIER_KEYS and the backend field group. Order/labels
// here are the only place the UI copy lives. `tab` decides which tab the group renders under:
// erosion & deposition go under "Erosion", everything plate/volcano-related under "Tectonics".
const TUNING_GROUPS: { heading: string; tab: TabKey; knobs: { key: TuningKey; label: string }[] }[] = [
  {
    heading: "Erosion",
    tab: "erosion",
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
    tab: "erosion",
    knobs: [
      { key: "river_deposition_multiplier", label: "River deposition (floodplains / deltas)" },
      { key: "ocean_deposition_multiplier", label: "Ocean deposition (shelf / beaches)" },
    ],
  },
  {
    heading: "Tectonics & volcanism",
    tab: "tectonics",
    knobs: [
      { key: "collision_uplift_multiplier", label: "Collision uplift — amount" },
      { key: "collision_uplift_reach_multiplier", label: "Collision uplift — reach" },
      { key: "volcanism_multiplier", label: "Volcanism" },
    ],
  },
];

const TABS: { key: TabKey; label: string }[] = [
  { key: "climate", label: "Climate" },
  { key: "erosion", label: "Erosion" },
  { key: "tectonics", label: "Tectonics" },
];

const ICE_AGE_MAX_YEARS = 1_000_000;
const ICE_AGE_STEP_YEARS = 100_000;

const selectStyle = {
  width: "100%",
  background: "#0f1424",
  color: "#e6e8ef",
  border: "1px solid #333",
  borderRadius: 4,
  padding: "4px 6px",
  fontSize: 13,
} as const;

// Real-time controls for World's live-adjustable properties (see backend app/world.py's
// World.sea_level_m/World.solar_multiplier/World.ice_age_period_years/
// World.simulate_plate_movement/World.simulate_climate_biomes/World.wind_model and its
// tuning-knob field group, plus main.py's /world/controls) -- unlike the "Generate World"
// dialog's sliders, these apply to the *current* world immediately, no regenerate needed.
// Controlled from App.tsx's own state (not local state here) since the same values also need
// to reset on a fresh Generate. Grouped into Climate / Erosion / Tectonics tabs.
export default function ControlsModal({
  seaLevelM,
  solarMultiplier,
  iceAgePeriodYears,
  simulatePlateMovement,
  simulateClimateBiomes,
  windModel,
  faultDeformationMode,
  tuning,
  onSeaLevelChange,
  onSolarMultiplierChange,
  onIceAgePeriodChange,
  onSimulatePlateMovementChange,
  onSimulateClimateBiomesChange,
  onWindModelChange,
  onFaultDeformationModeChange,
  onTuningChange,
  onTuningReset,
  onClose,
}: Props) {
  const [tab, setTab] = useState<TabKey>("climate");
  const anyTuned = Object.values(tuning).some((v) => v !== 1);
  const bothSimsOff = !simulatePlateMovement && !simulateClimateBiomes;

  const tuningGroup = (heading: string, knobs: { key: TuningKey; label: string }[]) => (
    <div key={heading} style={{ marginBottom: 12 }}>
      <div style={{ fontSize: 11, textTransform: "uppercase", letterSpacing: 0.5, color: "#8b8fa3", marginBottom: 6 }}>
        {heading}
      </div>
      {knobs.map(({ key, label }) => (
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
  );

  const tuningResetRow = (
    <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 10 }}>
      <span style={{ fontSize: 11, color: "#999" }}>
        Dimensionless multipliers on individual processes — 1× is the model's default; changes apply immediately.
      </span>
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
            marginLeft: 8,
            flexShrink: 0,
          }}
        >
          Reset to 1×
        </button>
      )}
    </div>
  );

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

        <div style={{ display: "flex", gap: 4, marginBottom: 16 }}>
          {TABS.map(({ key, label }) => (
            <button
              key={key}
              type="button"
              onClick={() => setTab(key)}
              style={{
                flex: 1,
                background: tab === key ? "#26304d" : "transparent",
                border: "1px solid #333",
                borderColor: tab === key ? "#6b7cea" : "#333",
                color: tab === key ? "#e6e8ef" : "#8b8fa3",
                cursor: "pointer",
                fontSize: 13,
                fontWeight: tab === key ? 600 : 400,
                borderRadius: 4,
                padding: "5px 0",
              }}
            >
              {label}
            </button>
          ))}
        </div>

        {tab === "climate" && (
          <div>
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

            <label style={{ display: "block", marginBottom: 6 }}>
              Ice Age Frequency: {iceAgePeriodYears === 0 ? "Never" : `every ${iceAgePeriodYears.toLocaleString()} yr`}
              <input
                type="range"
                min={0}
                max={ICE_AGE_MAX_YEARS}
                step={ICE_AGE_STEP_YEARS}
                value={iceAgePeriodYears}
                onChange={(e) => onIceAgePeriodChange(Number(e.target.value))}
                style={{ width: "100%" }}
              />
            </label>
            <div style={{ fontSize: 11, color: "#999", marginBottom: 16 }}>
              {iceAgePeriodYears === 0
                ? "No glacial cycle — climate is driven only by tectonics, insolation and the solar-heat slider."
                : "Full period of a slow glacial↔interglacial temperature swing. Glacial maxima cool the whole world, spread polar ice over the sea, and lower sea level as ice sheets lock up ocean water."}
            </div>

            <div style={{ borderTop: "1px solid #333", paddingTop: 14 }}>
              <label style={{ display: "flex", alignItems: "center", gap: 8, cursor: "pointer" }}>
                <input
                  type="checkbox"
                  checked={simulateClimateBiomes}
                  onChange={(e) => onSimulateClimateBiomesChange(e.target.checked)}
                />
                Simulate climate &amp; biomes
              </label>
              {bothSimsOff && (
                <div style={{ fontSize: 11, color: "#999", marginTop: 8 }}>
                  Both this and "Simulate plate movement" (Tectonics tab) are off — stepping will only advance elapsed
                  years, nothing else changes.
                </div>
              )}
            </div>

            <div style={{ borderTop: "1px solid #333", paddingTop: 14, marginTop: 14 }}>
              <label style={{ display: "block", marginBottom: 6 }}>Wind model</label>
              <select value={windModel} onChange={(e) => onWindModelChange(e.target.value)} style={selectStyle}>
                <option value="cfd">Shallow-water CFD (accurate)</option>
                <option value="diagnostic">Diagnostic / ABL (fast)</option>
              </select>
              <div style={{ fontSize: 11, color: "#999", marginTop: 8 }}>
                {windModel === "diagnostic"
                  ? "Skips the shallow-water solve -- much faster steps, ~85-90% of the CFD biome map. Wind/temperature maps show the closed-form field."
                  : "Genuine time-integrated shallow-water solve -- the largest single cost of a step at high fluid-dynamics resolution."}
              </div>
            </div>
          </div>
        )}

        {tab === "erosion" && (
          <div>
            {tuningResetRow}
            {TUNING_GROUPS.filter((g) => g.tab === "erosion").map((g) => tuningGroup(g.heading, g.knobs))}
          </div>
        )}

        {tab === "tectonics" && (
          <div>
            <label style={{ display: "flex", alignItems: "center", gap: 8, marginBottom: 8, cursor: "pointer" }}>
              <input
                type="checkbox"
                checked={simulatePlateMovement}
                onChange={(e) => onSimulatePlateMovementChange(e.target.checked)}
              />
              Simulate plate movement
            </label>
            {bothSimsOff && (
              <div style={{ fontSize: 11, color: "#999", marginBottom: 8 }}>
                Both this and "Simulate climate &amp; biomes" (Climate tab) are off — stepping will only advance elapsed
                years, nothing else changes.
              </div>
            )}

            <div style={{ borderTop: "1px solid #333", paddingTop: 14, marginTop: 6 }}>
              <label style={{ display: "block", marginBottom: 6 }}>Fault deformation model</label>
              <select
                value={faultDeformationMode}
                onChange={(e) => onFaultDeformationModeChange(e.target.value)}
                style={selectStyle}
              >
                <option value="fault">Along fault lines (default)</option>
                <option value="boundary">Boundary bands (classic)</option>
                <option value="both">Both (superimposed)</option>
              </select>
              <div style={{ fontSize: 11, color: "#999", marginTop: 8 }}>
                {faultDeformationMode === "boundary"
                  ? "Uplift / rifting apply as smooth bands at the plate-polygon edge -- the pre-faults-rework behaviour."
                  : faultDeformationMode === "fault"
                    ? "Boundary uplift / rifting concentrates onto active fault traces, and the fault relief layer is scaled up to carry it. Faults spawn boundary-hugging, so the collision zone reads as fault-tracking ridges rather than one smooth swell."
                    : "Boundary bands at full strength plus the scaled-up fault relief layer on top."}
              </div>
            </div>

            <div style={{ borderTop: "1px solid #333", paddingTop: 14, marginTop: 14 }}>
              {tuningResetRow}
              {TUNING_GROUPS.filter((g) => g.tab === "tectonics").map((g) => tuningGroup(g.heading, g.knobs))}
            </div>
          </div>
        )}
      </div>
    </div>
  );
}
