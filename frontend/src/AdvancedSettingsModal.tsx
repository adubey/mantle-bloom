interface Props {
  landPercent: number;
  continentalPercent: number;
  autoPlates: boolean;
  numPlates: number;
  minPlates: number;
  maxPlates: number;
  axialTiltDeg: number;
  initialSoilMaturityPercent: number;
  fluidDensity: number;
  // Same {value, label} choices App.tsx's own "Detail" dropdown uses (App.tsx's
  // DETAIL_CHOICES) -- passed down rather than imported here to avoid a circular import
  // between this file and App.tsx.
  fluidDensityChoices: { value: number; label: string }[];
  onLandPercentChange: (v: number) => void;
  onContinentalPercentChange: (v: number) => void;
  onAutoPlatesChange: (v: boolean) => void;
  onNumPlatesChange: (v: number) => void;
  onAxialTiltDegChange: (v: number) => void;
  onInitialSoilMaturityPercentChange: (v: number) => void;
  onFluidDensityChange: (v: number) => void;
  onClose: () => void;
}

// The "Advanced settings" sub-window opened from the "Generate World" dialog (see App.tsx),
// holding the less commonly tweaked generation options -- unlike ControlsModal, these only
// take effect on the next Generate, not immediately. Controlled from App.tsx's own state (not
// local state here) since Generate itself reads these same values.
export default function AdvancedSettingsModal({
  landPercent,
  continentalPercent,
  autoPlates,
  numPlates,
  minPlates,
  maxPlates,
  axialTiltDeg,
  initialSoilMaturityPercent,
  fluidDensity,
  fluidDensityChoices,
  onLandPercentChange,
  onContinentalPercentChange,
  onAutoPlatesChange,
  onNumPlatesChange,
  onAxialTiltDegChange,
  onInitialSoilMaturityPercentChange,
  onFluidDensityChange,
  onClose,
}: Props) {
  return (
    <div
      onClick={onClose}
      style={{
        position: "fixed",
        inset: 0,
        background: "rgba(0, 0, 0, 0.6)",
        zIndex: 200,
        display: "flex",
        alignItems: "center",
        justifyContent: "center",
      }}
    >
      <div
        onClick={(e) => e.stopPropagation()}
        style={{
          width: 380,
          maxWidth: "66vw",
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
          <span style={{ fontSize: 16, fontWeight: 600 }}>Advanced settings</span>
          <button
            type="button"
            onClick={onClose}
            style={{ background: "none", border: "none", color: "#8b8fa3", cursor: "pointer", fontSize: 16 }}
          >
            ✕
          </button>
        </div>

        <label style={{ display: "block", marginBottom: 16 }}>
          Initial land: {landPercent}%
          <input
            type="range"
            min={0}
            max={100}
            value={landPercent}
            onChange={(e) => onLandPercentChange(Number(e.target.value))}
            style={{ width: "100%" }}
          />
        </label>

        <label style={{ display: "block", marginBottom: 16 }}>
          Continental plates: {continentalPercent}%
          <input
            type="range"
            min={0}
            max={100}
            value={continentalPercent}
            onChange={(e) => onContinentalPercentChange(Number(e.target.value))}
            style={{ width: "100%" }}
          />
        </label>

        <label style={{ display: "flex", alignItems: "center", gap: 6, marginBottom: 6, fontSize: 12 }}>
          <input type="checkbox" checked={autoPlates} onChange={(e) => onAutoPlatesChange(e.target.checked)} />
          Number of plates: Auto (seed-based)
        </label>
        {!autoPlates && (
          <label style={{ display: "block", marginBottom: 16 }}>
            Plates: {numPlates}
            <input
              type="range"
              min={minPlates}
              max={maxPlates}
              value={numPlates}
              onChange={(e) => onNumPlatesChange(Number(e.target.value))}
              style={{ width: "100%" }}
            />
          </label>
        )}

        <label style={{ display: "block", marginBottom: 16 }}>
          Axial tilt: {axialTiltDeg}°
          <input
            type="range"
            min={0}
            max={45}
            step={0.5}
            value={axialTiltDeg}
            onChange={(e) => onAxialTiltDegChange(Number(e.target.value))}
            style={{ width: "100%" }}
          />
        </label>

        <label style={{ display: "block", marginBottom: 16 }}>
          Initial soil maturity: {initialSoilMaturityPercent}%
          <input
            type="range"
            min={0}
            max={100}
            value={initialSoilMaturityPercent}
            onChange={(e) => onInitialSoilMaturityPercentChange(Number(e.target.value))}
            style={{ width: "100%" }}
          />
          <div style={{ fontSize: 11, color: "#999", marginTop: 4 }}>
            0% starts fully barren (bare rock, no soil) -- soil then forms gradually as the
            world steps forward. Higher values start with some soil already in place.
          </div>
        </label>

        <label style={{ display: "block", marginBottom: 0 }}>
          Fluid dynamics resolution
          <select
            value={fluidDensity}
            onChange={(e) => onFluidDensityChange(Number(e.target.value))}
            style={{ width: "100%", marginTop: 4 }}
          >
            {fluidDensityChoices.map((d) => (
              <option key={d.value} value={d.value}>
                {d.label}
              </option>
            ))}
          </select>
          <div style={{ fontSize: 11, color: "#999", marginTop: 4 }}>
            How finely Ocean/Atmospheric Fluid Dynamics mode resolves currents and wind --
            independent of Detail above, so you can keep sharp climate/biome maps while running
            Fluid Dynamics mode at a coarser (faster) resolution. Lower runs faster but coarser.
          </div>
        </label>
      </div>
    </div>
  );
}
