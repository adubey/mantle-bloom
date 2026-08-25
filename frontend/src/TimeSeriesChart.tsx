import { useMemo, useRef, useState } from "react";
import type { PointerEvent as ReactPointerEvent } from "react";

export interface ChartSeries {
  key: string;
  label: string;
  color: string;
}

export interface ChartPoint {
  x: number; // elapsed years
  values: Record<string, number | null>;
}

// A shaded region between two of a point's `values` keys (e.g. mean-std/mean+std), drawn
// underneath every line series -- used for the "avg with a std-dev band" grouped-stat charts
// (see StatsModal.tsx) rather than listing min/max/avg/std-dev as separate line choices.
export interface ChartBand {
  lowKey: string;
  highKey: string;
  color: string;
  label: string;
}

interface TimeSeriesChartProps {
  series: ChartSeries[];
  data: ChartPoint[];
  yFormat: (v: number) => string;
  bands?: ChartBand[];
  height?: number;
}

const CHART_WIDTH = 460;
const MARGIN = { top: 10, right: 14, bottom: 26, left: 54 };
// Matches this app's own established dark chrome (EventConsole/App.tsx's Generate World
// dialog already use #0d1120/#151a2e surfaces, #333 hairlines, opacity-based muted text).
const GRIDLINE_COLOR = "#333";
const AXIS_TEXT_COLOR = "#8b8fa3";
const SURFACE_COLOR = "#151a2e";

// Matches EventConsole's own elapsed-years formatting ("X.X Myr") rather than inventing a
// second convention for the same unit.
function formatMyr(years: number): string {
  return `${(years / 1e6).toFixed(1)} Myr`;
}

function niceTicks(rawMin: number, rawMax: number, count = 4): number[] {
  let min = rawMin;
  let max = rawMax;
  if (min === max) {
    // A flat series (common early in a run, or for a metric that barely moves) needs
    // headroom scaled to its own magnitude, not a flat +-1 -- which would dwarf a small
    // fraction like land_fraction=0.35 into a -100%..150% axis.
    const pad = min === 0 ? 1 : Math.abs(min) * 0.1;
    min -= pad;
    max += pad;
  }
  const range = max - min;
  const magnitude = 10 ** Math.floor(Math.log10(range / count));
  const candidates = [1, 2, 2.5, 5, 10];
  let step = candidates[candidates.length - 1] * magnitude;
  for (const c of candidates) {
    if (range / (c * magnitude) <= count) {
      step = c * magnitude;
      break;
    }
  }
  const niceMin = Math.floor(min / step) * step;
  const niceMax = Math.ceil(max / step) * step;
  const ticks: number[] = [];
  for (let v = niceMin; v <= niceMax + step * 0.5; v += step) {
    ticks.push(Math.round(v / step) * step);
  }
  return ticks;
}

// A small, dependency-free multi-series line chart for a value plotted against elapsed
// years, avoiding a charting library for a relatively simple need. 2px round-cap lines,
// hairline gridlines, an always-on crosshair+tooltip, a legend for 2+ series, an 8px
// surface-ringed end-dot per series. Single-series charts skip the legend -- the caller's
// own heading/dropdown already names what's plotted.
export default function TimeSeriesChart({ series, data, yFormat, bands, height = 200 }: TimeSeriesChartProps) {
  const [hoverIndex, setHoverIndex] = useState<number | null>(null);
  const svgRef = useRef<SVGSVGElement>(null);

  const plotW = CHART_WIDTH - MARGIN.left - MARGIN.right;
  const plotH = height - MARGIN.top - MARGIN.bottom;

  const xs = data.map((d) => d.x);
  const xMin = xs.length ? Math.min(...xs) : 0;
  const xMax = xs.length ? Math.max(...xs) : 1;

  let yMin = Infinity;
  let yMax = -Infinity;
  for (const point of data) {
    for (const s of series) {
      const v = point.values[s.key];
      if (v !== null && Number.isFinite(v)) {
        if (v < yMin) yMin = v;
        if (v > yMax) yMax = v;
      }
    }
    // A band isn't necessarily bounded by its series' own values (e.g. mean +/- std can
    // exceed the plotted min/max lines for a skewed distribution), so it needs its own say
    // in the y-domain rather than assuming the lines already cover it.
    for (const b of bands ?? []) {
      for (const key of [b.lowKey, b.highKey]) {
        const v = point.values[key];
        if (v !== null && Number.isFinite(v)) {
          if (v < yMin) yMin = v;
          if (v > yMax) yMax = v;
        }
      }
    }
  }
  const hasData = yMin !== Infinity;
  const yTicks = hasData ? niceTicks(yMin, yMax) : [0, 1];
  const yLo = yTicks[0];
  const yHi = yTicks[yTicks.length - 1];

  const xScale = (x: number) => (xMax === xMin ? plotW / 2 : ((x - xMin) / (xMax - xMin)) * plotW);
  const yScale = (v: number) => (yHi === yLo ? plotH / 2 : plotH - ((v - yLo) / (yHi - yLo)) * plotH);

  const paths = useMemo(() => {
    return series.map((s) => {
      let d = "";
      let open = false;
      for (const point of data) {
        const v = point.values[s.key];
        if (v === null || !Number.isFinite(v)) {
          open = false;
          continue;
        }
        const x = xScale(point.x);
        const y = yScale(v);
        d += open ? ` L ${x.toFixed(1)} ${y.toFixed(1)}` : ` M ${x.toFixed(1)} ${y.toFixed(1)}`;
        open = true;
      }
      return { key: s.key, color: s.color, d };
    });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [series, data, xMin, xMax, yLo, yHi]);

  // Filled low/high polygons, one contiguous run at a time (same gap-handling as `paths`
  // above, just building a closed shape -- forward along the low edge, then back along the
  // high edge -- instead of an open line).
  const bandPaths = useMemo(() => {
    return (bands ?? []).map((b) => {
      const segments: string[] = [];
      let lowPts: [number, number][] = [];
      let highPts: [number, number][] = [];
      const flush = () => {
        if (lowPts.length >= 2) {
          const forward = lowPts.map(([x, y], i) => `${i === 0 ? "M" : "L"} ${x.toFixed(1)} ${y.toFixed(1)}`).join(" ");
          const backward = [...highPts].reverse().map(([x, y]) => `L ${x.toFixed(1)} ${y.toFixed(1)}`).join(" ");
          segments.push(`${forward} ${backward} Z`);
        }
        lowPts = [];
        highPts = [];
      };
      for (const point of data) {
        const lo = point.values[b.lowKey];
        const hi = point.values[b.highKey];
        if (lo === null || hi === null || !Number.isFinite(lo) || !Number.isFinite(hi)) {
          flush();
          continue;
        }
        const x = xScale(point.x);
        lowPts.push([x, yScale(lo)]);
        highPts.push([x, yScale(hi)]);
      }
      flush();
      return { key: `${b.lowKey}-${b.highKey}`, color: b.color, d: segments.join(" ") };
    });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [bands, data, xMin, xMax, yLo, yHi]);

  const lastValid = (key: string) => {
    for (let i = data.length - 1; i >= 0; i--) {
      const v = data[i].values[key];
      if (v !== null && Number.isFinite(v)) return { point: data[i], value: v };
    }
    return null;
  };

  function handlePointerMove(e: ReactPointerEvent<SVGSVGElement>) {
    if (!svgRef.current || data.length === 0) return;
    const rect = svgRef.current.getBoundingClientRect();
    const userX = ((e.clientX - rect.left) / rect.width) * CHART_WIDTH - MARGIN.left;
    let closest = 0;
    let closestDist = Infinity;
    for (let i = 0; i < data.length; i++) {
      const dist = Math.abs(xScale(data[i].x) - userX);
      if (dist < closestDist) {
        closestDist = dist;
        closest = i;
      }
    }
    setHoverIndex(closest);
  }

  const hovered = hoverIndex !== null ? data[hoverIndex] : null;
  const hoveredX = hovered ? MARGIN.left + xScale(hovered.x) : null;

  if (!hasData) {
    return <div style={{ padding: "24px 0", textAlign: "center", opacity: 0.5, fontSize: 12 }}>No data yet</div>;
  }

  return (
    <div style={{ position: "relative" }}>
      <svg
        ref={svgRef}
        viewBox={`0 0 ${CHART_WIDTH} ${height}`}
        width="100%"
        height={height}
        style={{ display: "block", touchAction: "none" }}
        onPointerMove={handlePointerMove}
        onPointerLeave={() => setHoverIndex(null)}
      >
        <g transform={`translate(${MARGIN.left}, ${MARGIN.top})`}>
          {yTicks.map((t) => (
            <g key={t}>
              <line x1={0} x2={plotW} y1={yScale(t)} y2={yScale(t)} stroke={GRIDLINE_COLOR} strokeWidth={1} />
              <text x={-8} y={yScale(t)} dy={4} textAnchor="end" fontSize={10} fill={AXIS_TEXT_COLOR}>
                {yFormat(t)}
              </text>
            </g>
          ))}

          {bandPaths.map((b) => (
            <path key={b.key} d={b.d} fill={b.color} stroke="none" />
          ))}

          {/* Drawn after the bands, and in `series` order -- the caller lists its "avg" line
              last precisely so it paints on top of both the std-dev band and the min/max
              lines, per the panel's own "avg visible atop the band" requirement. */}
          {paths.map((p) => (
            <path key={p.key} d={p.d} fill="none" stroke={p.color} strokeWidth={2} strokeLinecap="round" strokeLinejoin="round" />
          ))}

          {series.map((s) => {
            const last = lastValid(s.key);
            if (!last) return null;
            const cx = xScale(last.point.x);
            const cy = yScale(last.value);
            return (
              <g key={s.key}>
                <circle cx={cx} cy={cy} r={5} fill={SURFACE_COLOR} />
                <circle cx={cx} cy={cy} r={4} fill={s.color} />
              </g>
            );
          })}

          {hoveredX !== null && (
            <line
              x1={hoveredX - MARGIN.left}
              x2={hoveredX - MARGIN.left}
              y1={0}
              y2={plotH}
              stroke={AXIS_TEXT_COLOR}
              strokeWidth={1}
              strokeDasharray="3 3"
            />
          )}

          <text x={0} y={plotH + 18} fontSize={10} fill={AXIS_TEXT_COLOR}>
            {formatMyr(xMin)}
          </text>
          <text x={plotW} y={plotH + 18} textAnchor="end" fontSize={10} fill={AXIS_TEXT_COLOR}>
            {formatMyr(xMax)}
          </text>
        </g>
      </svg>

      {hovered && (
        <div
          style={{
            position: "absolute",
            top: 4,
            left: hoveredX !== null && hoveredX > CHART_WIDTH / 2 ? undefined : "60%",
            right: hoveredX !== null && hoveredX > CHART_WIDTH / 2 ? "4%" : undefined,
            background: "#0d1120",
            border: "1px solid #333",
            borderRadius: 4,
            padding: "6px 8px",
            fontSize: 11,
            pointerEvents: "none",
            minWidth: 120,
          }}
        >
          <div style={{ opacity: 0.6, marginBottom: 4 }}>{formatMyr(hovered.x)}</div>
          {series.map((s) => {
            const v = hovered.values[s.key];
            return (
              <div key={s.key} style={{ display: "flex", alignItems: "center", gap: 6, padding: "1px 0" }}>
                <span style={{ display: "inline-block", width: 10, height: 2, background: s.color, flexShrink: 0 }} />
                <span style={{ opacity: 0.7, flex: 1 }}>{s.label}</span>
                <span style={{ fontWeight: 600 }}>{v === null || !Number.isFinite(v) ? "--" : yFormat(v)}</span>
              </div>
            );
          })}
        </div>
      )}

      {(series.length > 1 || (bands?.length ?? 0) > 0) && (
        <div style={{ display: "flex", flexWrap: "wrap", gap: "4px 12px", marginTop: 6 }}>
          {series.map((s) => (
            <div key={s.key} style={{ display: "flex", alignItems: "center", gap: 5, fontSize: 11 }}>
              <span style={{ display: "inline-block", width: 10, height: 2, flexShrink: 0, background: s.color }} />
              <span style={{ opacity: 0.75 }}>{s.label}</span>
            </div>
          ))}
          {(bands ?? []).map((b) => (
            <div key={`${b.lowKey}-${b.highKey}`} style={{ display: "flex", alignItems: "center", gap: 5, fontSize: 11 }}>
              <span style={{ display: "inline-block", width: 10, height: 10, flexShrink: 0, background: b.color }} />
              <span style={{ opacity: 0.75 }}>{b.label}</span>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
