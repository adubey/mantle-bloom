interface Props {
  seaLevelM: number;
  solarMultiplier: number;
  onSeaLevelChange: (v: number) => void;
  onSolarMultiplierChange: (v: number) => void;
  onClose: () => void;
}

// Real-time controls for two of World's live-adjustable properties (see backend
// app/world.py's World.sea_level_m/World.solar_multiplier and main.py's /world/controls) --
// unlike the "Generate World" dialog's sliders, these apply to the *current* world
// immediately, no regenerate needed. Sliders are controlled from App.tsx's own state (not
// local state here) since the same values also need to reset on a fresh Generate.
export default function ControlsModal({ seaLevelM, solarMultiplier, onSeaLevelChange, onSolarMultiplierChange, onClose }: Props) {
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
          width: 320,
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

        <label style={{ display: "block" }}>
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
      </div>
    </div>
  );
}
