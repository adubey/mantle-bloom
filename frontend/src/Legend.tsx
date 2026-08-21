import type { MapView } from "./api";
import { legendFor } from "./legendData";
import type { LegendGradient, LegendSymbol, SwatchKind } from "./legendData";

interface Props {
  mapView: MapView;
  // Legend-click-to-highlight (Biome and Combined views only -- see App.tsx/MapCanvas.tsx):
  // the currently highlighted swatch's label, and a callback fired with a swatch's label on
  // click (App.tsx toggles it off if the same label is clicked again). Both omitted on every
  // other view, where the legend stays the plain, non-interactive read-only panel it always was.
  highlightedBiome?: string | null;
  onBiomeClick?: (label: string) => void;
}

const SWATCH_SIZE = 14;

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
export default function Legend({ mapView, highlightedBiome, onBiomeClick }: Props) {
  const spec = legendFor(mapView);
  if (!spec) return null;

  // Biome and Combined are the only views whose swatches are clickable -- see the Props doc
  // comment.
  const clickable = (mapView === "biome" || mapView === "combined") && !!onBiomeClick;

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
      </div>
      <div style={{ display: "flex", flexWrap: "wrap", columnGap: 16, rowGap: 2, marginTop: spec.gradient ? 2 : 8 }}>
        {spec.symbols.map((sym) => (
          <SymbolRow
            key={sym.label}
            symbol={sym}
            onClick={clickable ? () => onBiomeClick!(sym.label) : undefined}
            selected={clickable && highlightedBiome === sym.label}
          />
        ))}
      </div>
    </div>
  );
}
