import type { MapView } from "./api";
import { faultKindForLegendLabel, legendFor } from "./legendData";
import type { LegendGradient, LegendSymbol, SwatchKind } from "./legendData";

interface Props {
  mapView: MapView;
  // False before any world has been generated or loaded (see App.tsx's `summary`) -- the real
  // legend has nothing to describe yet, so the panel shows a prompt to get started instead.
  hasWorld: boolean;
  // Legend-click-to-highlight (Biome and Combined views only -- see App.tsx/MapCanvas.tsx):
  // the currently highlighted swatch's label, and a callback fired with a swatch's label on
  // click (App.tsx toggles it off if the same label is clicked again). Both omitted on every
  // other view, where the legend stays the plain, non-interactive read-only panel it always was.
  highlightedBiome?: string | null;
  onBiomeClick?: (label: string) => void;
  // Elevation view's own "Mountains" / "Plains & Plateaus" toggle rows -- independent
  // checkboxes (either, both, or neither can be on), unlike highlightedBiome's single-select
  // click-to-highlight above, and baked server-side into the rendered image rather than a
  // client-side redraw (see App.tsx's terrainOverlayRef/refresh). Omitted on every other view.
  showMountains?: boolean;
  onToggleMountains?: () => void;
  showPlainsPlateaus?: boolean;
  onTogglePlainsPlateaus?: () => void;
}

const SWATCH_SIZE = 14;

// The map's own km-per-pixel scale, for the ScaleBar below. The world is always Earth-sized
// (see backend app/elevation_lines.py's PLANET_RADIUS_KM) and both current projections
// (behrmann/eckert4 -- see backend app/projections.py) fit the whole sphere to the rendered
// width/height with their widest point at the equator (backend render_image.py's
// _fit_and_project_sphere, fit against PADDING_PX at REFERENCE_WIDTH_PX -- both equal to this
// app's own MAP_DISPLAY_WIDTH_PX/MAP_PADDING_PX below). For any projection shaped that way, the
// equatorial scale works out to just the planet's circumference divided by the available map
// width, independent of the projection's own math -- so this doesn't need to duplicate
// projections.py's per-projection formulas, only the constants that decide how big the map is.
// Away from the equator this equal-area map's horizontal scale changes with latitude (area is
// preserved, not shape), so the bar is only exact there -- the same caveat any small-scale
// world map's scale bar carries.
const PLANET_RADIUS_KM = 6371; // matches backend app/elevation_lines.py's PLANET_RADIUS_KM
const MAP_DISPLAY_WIDTH_PX = 1100; // matches App.tsx's DISPLAY_WIDTH
const MAP_PADDING_PX = 20; // matches backend app/render_image.py's PADDING_PX
const KM_PER_PIXEL_AT_EQUATOR =
  (2 * Math.PI * PLANET_RADIUS_KM) / (MAP_DISPLAY_WIDTH_PX - 2 * MAP_PADDING_PX);
const PX_PER_KM_AT_EQUATOR = 1 / KM_PER_PIXEL_AT_EQUATOR;
// The ruler's target on-screen length -- floored down to the nearest whole major (1,000 km)
// division, never past it, so the bar always ends exactly on a major tick rather than a
// partial one. Major ticks (taller) fall every SCALE_BAR_MAJOR_STEP_KM, minor ticks (shorter)
// every SCALE_BAR_MINOR_STEP_KM in between.
const SCALE_BAR_TARGET_PX = 150;
const SCALE_BAR_MAJOR_STEP_KM = 1000;
const SCALE_BAR_MINOR_STEP_KM = 100;

function ScaleBar() {
  const roughKm = KM_PER_PIXEL_AT_EQUATOR * SCALE_BAR_TARGET_PX;
  const totalKm = Math.max(SCALE_BAR_MAJOR_STEP_KM, Math.floor(roughKm / SCALE_BAR_MAJOR_STEP_KM) * SCALE_BAR_MAJOR_STEP_KM);
  const widthPx = totalKm * PX_PER_KM_AT_EQUATOR;
  const baselineY = 3;
  const majorTickLen = 6;
  const minorTickLen = 3;
  const minorTicks: number[] = [];
  for (let km = SCALE_BAR_MINOR_STEP_KM; km < totalKm; km += SCALE_BAR_MINOR_STEP_KM) {
    if (km % SCALE_BAR_MAJOR_STEP_KM !== 0) minorTicks.push(km);
  }
  const majorTicks: number[] = [];
  for (let km = 0; km <= totalKm; km += SCALE_BAR_MAJOR_STEP_KM) majorTicks.push(km);

  return (
    <div
      style={{ flexShrink: 0, marginLeft: "auto", fontSize: 10, opacity: 0.85, textAlign: "right" }}
      title={`Ticks every ${SCALE_BAR_MAJOR_STEP_KM.toLocaleString()}/${SCALE_BAR_MINOR_STEP_KM} km. Exact at the equator -- this equal-area map's horizontal scale changes with latitude.`}
    >
      <svg width={widthPx} height={baselineY + majorTickLen + 1} aria-hidden>
        <line x1={0} y1={baselineY} x2={widthPx} y2={baselineY} stroke="#dee2eb" strokeWidth={1} />
        {minorTicks.map((km) => {
          const x = km * PX_PER_KM_AT_EQUATOR;
          return <line key={km} x1={x} y1={baselineY} x2={x} y2={baselineY + minorTickLen} stroke="#dee2eb" strokeWidth={1} />;
        })}
        {majorTicks.map((km) => {
          const x = km * PX_PER_KM_AT_EQUATOR;
          return <line key={km} x1={x} y1={baselineY} x2={x} y2={baselineY + majorTickLen} stroke="#dee2eb" strokeWidth={1.5} />;
        })}
      </svg>
      <div style={{ marginTop: 2, whiteSpace: "nowrap" }}>{totalKm.toLocaleString()} km at equator</div>
    </div>
  );
}

function Swatch({ kind, color, outline }: { kind: SwatchKind; color: string; outline?: string }) {
  const s = SWATCH_SIZE;
  switch (kind) {
    case "line":
      return (
        <svg width={s} height={s} aria-hidden>
          <line x1={1} y1={s / 2} x2={s - 1} y2={s / 2} stroke={color} strokeWidth={2} />
        </svg>
      );
    case "square":
      return (
        <svg width={s} height={s} aria-hidden>
          <rect x={0} y={2} width={s} height={s - 4} fill={color} />
        </svg>
      );
    case "circle":
      return (
        <svg width={s} height={s} aria-hidden>
          <circle cx={s / 2} cy={s / 2} r={s / 2 - 1} fill={color} stroke={outline ?? "#ffffff"} strokeWidth={1} />
        </svg>
      );
    case "ring":
      return (
        <svg width={s} height={s} aria-hidden>
          <circle cx={s / 2} cy={s / 2} r={s / 2 - 1.5} fill="none" stroke={color} strokeWidth={1.5} />
        </svg>
      );
    case "arrow":
      return (
        <svg width={s} height={s} aria-hidden>
          <line x1={1} y1={s / 2} x2={s - 4} y2={s / 2} stroke={color} strokeWidth={1.5} />
          <polygon points={`${s - 1},${s / 2} ${s - 5},${s / 2 - 3} ${s - 5},${s / 2 + 3}`} fill={color} />
        </svg>
      );
    case "arc":
      return (
        <svg width={s} height={s} viewBox={`0 0 ${s} ${s}`} aria-hidden>
          <path
            d={`M ${s * 0.15} ${s * 0.75} A ${s * 0.4} ${s * 0.4} 0 1 1 ${s * 0.85} ${s * 0.75}`}
            fill="none"
            stroke={color}
            strokeWidth={1.5}
          />
          <polygon points={`${s * 0.85},${s * 0.75} ${s * 0.68},${s * 0.66} ${s * 0.72},${s * 0.9}`} fill={color} />
        </svg>
      );
    default:
      return null;
  }
}

function GradientBar({ gradient }: { gradient: LegendGradient }) {
  const span = gradient.max - gradient.min;
  const pct = (v: number) => (span <= 0 ? 0 : ((v - gradient.min) / span) * 100);
  const css = `linear-gradient(to right, ${gradient.stops.map((s) => `${s.color} ${pct(s.value)}%`).join(", ")})`;
  return (
    <div style={{ marginBottom: 4 }}>
      <div style={{ height: 14, borderRadius: 2, background: css }} />
      <div style={{ position: "relative", height: 14 }}>
        {gradient.ticks.map((t) => (
          <span
            key={t.label}
            style={{
              position: "absolute",
              left: `${pct(t.value)}%`,
              transform: "translateX(-50%)",
              fontSize: 10,
              opacity: 0.85,
              whiteSpace: "nowrap",
            }}
          >
            {t.label}
          </span>
        ))}
      </div>
    </div>
  );
}

function SymbolRow({ symbol, onClick, selected }: { symbol: LegendSymbol; onClick?: () => void; selected?: boolean }) {
  return (
    <div
      onClick={onClick}
      style={{
        display: "flex",
        alignItems: "center",
        gap: 6,
        padding: "2px 4px",
        margin: "0 -4px",
        borderRadius: 3,
        cursor: onClick ? "pointer" : undefined,
        background: selected ? "rgba(255, 255, 255, 0.16)" : "transparent",
      }}
    >
      <Swatch kind={symbol.kind} color={symbol.color} outline={symbol.outline} />
      <span style={{ fontSize: 11, fontWeight: selected ? 700 : 400 }}>{symbol.label}</span>
    </div>
  );
}

// Rendered client-side (see legendData.ts for why) as a full-width HTML bar below the map,
// in normal document flow (see App.tsx) rather than overlaid on top of it -- unlike that old
// server-baked panel, which sat in the map's bottom-left corner and partially covered it,
// this one never obscures any part of the map. Still updates instantly on a projection/view
// change and stays legible during a live rotation-preview drag (which only ever re-draws the
// canvas above it).
export default function Legend({
  mapView,
  hasWorld,
  highlightedBiome,
  onBiomeClick,
  showMountains,
  onToggleMountains,
  showPlainsPlateaus,
  onTogglePlainsPlateaus,
}: Props) {
  if (!hasWorld) {
    return (
      <div
        style={{
          marginTop: 10,
          padding: 10,
          background: "rgba(16, 20, 34, 0.92)",
          border: "1px solid #4b5060",
          borderRadius: 4,
          color: "#dee2eb",
          userSelect: "none",
          fontSize: 13,
        }}
      >
        Click &apos;Generate World&apos; or &apos;File&apos; to begin.
      </div>
    );
  }

  const spec = legendFor(mapView);
  if (!spec) return null;

  // Views with clickable swatches: Biome / Combined / "Last elevation change" match on exact
  // pixel colour (see the Props doc comment); "Plates & Faults" instead resolves a clicked
  // fault-type row to a FaultKind the client-drawn view isolates (see faultKindForLegendLabel).
  const clickable =
    (mapView === "biome" || mapView === "combined" || mapView === "elevReason" || mapView === "platesAndFaults") &&
    !!onBiomeClick;
  // Elevation view's own toggle rows (see the Props doc comment) -- independent checkboxes,
  // handled entirely separately from the click-to-highlight `clickable` mechanism above.
  const TERRAIN_TOGGLES: Record<string, { on?: boolean; toggle?: () => void }> = {
    Mountains: { on: showMountains, toggle: onToggleMountains },
    "Plains & Plateaus": { on: showPlainsPlateaus, toggle: onTogglePlainsPlateaus },
  };

  return (
    <div
      style={{
        marginTop: 10,
        padding: 10,
        background: "rgba(16, 20, 34, 0.92)",
        border: "1px solid #4b5060",
        borderRadius: 4,
        color: "#dee2eb",
        userSelect: "none",
      }}
    >
      <div style={{ display: "flex", alignItems: "center", gap: 16, flexWrap: "wrap" }}>
        <div style={{ fontSize: 13, fontWeight: 600, flexShrink: 0 }}>{spec.title}</div>
        {spec.gradient && (
          <div style={{ flex: "1 1 260px", minWidth: 200, maxWidth: 420 }}>
            <GradientBar gradient={spec.gradient} />
          </div>
        )}
        {spec.gradients?.map(({ label, gradient }) => (
          <div key={label} style={{ flex: "1 1 260px", minWidth: 200, maxWidth: 420 }}>
            <div style={{ fontSize: 10, opacity: 0.65, marginBottom: 1 }}>{label}</div>
            <GradientBar gradient={gradient} />
          </div>
        ))}
        <ScaleBar />
      </div>
      <div style={{ display: "flex", flexWrap: "wrap", columnGap: 16, rowGap: 2, marginTop: spec.gradient || spec.gradients ? 2 : 8 }}>
        {spec.symbols.map((sym) => {
          const terrainToggle = mapView === "elevation" ? TERRAIN_TOGGLES[sym.label] : undefined;
          if (terrainToggle) {
            return (
              <SymbolRow key={sym.label} symbol={sym} onClick={terrainToggle.toggle} selected={!!terrainToggle.on} />
            );
          }
          // On "Plates & Faults" only the three fault-type rows do anything; elsewhere every
          // row but "Coastline" (a plain orientation cue) is clickable.
          const rowClickable =
            clickable &&
            (mapView === "platesAndFaults"
              ? faultKindForLegendLabel(sym.label) !== null
              : sym.label !== "Coastline");
          return (
            <SymbolRow
              key={sym.label}
              symbol={sym}
              onClick={rowClickable ? () => onBiomeClick!(sym.label) : undefined}
              selected={rowClickable && highlightedBiome === sym.label}
            />
          );
        })}
      </div>
    </div>
  );
}
