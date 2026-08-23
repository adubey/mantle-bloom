import { useMemo, useState } from "react";
import type { WorldStats } from "./api";
import TimeSeriesChart from "./TimeSeriesChart";
import type { ChartPoint, ChartSeries } from "./TimeSeriesChart";

interface Props {
  stats: WorldStats | null;
  history: WorldStats[];
  onClose: () => void;
}

interface Metric {
  key: string;
  label: string;
  get: (s: WorldStats) => number | null;
  tableFormat: (v: number | null) => string;
  yFormat: (v: number) => string;
}

function pctMetric(key: string, label: string, get: (s: WorldStats) => number | null): Metric {
  return {
    key,
    label,
    get,
    tableFormat: (v) => (v === null ? "--" : `${(v * 100).toFixed(1)}%`),
    yFormat: (v) => `${(v * 100).toFixed(0)}%`,
  };
}

function numMetric(key: string, label: string, get: (s: WorldStats) => number | null, digits = 1, suffix = ""): Metric {
  return {
    key,
    label,
    get,
    tableFormat: (v) => (v === null ? "--" : `${v.toFixed(digits)}${suffix}`),
    yFormat: (v) => `${v.toFixed(digits)}${suffix}`,
  };
}

type TabKey = "physical" | "temperature" | "precipitation" | "biome" | "simulation";

const TABS: { key: TabKey; label: string }[] = [
  { key: "physical", label: "Physical" },
  { key: "temperature", label: "Temperature" },
  { key: "precipitation", label: "Precipitation" },
  { key: "biome", label: "Biome" },
  { key: "simulation", label: "Simulation" },
];

// Fixed order (biomes.BIOME_NAMES minus "Ocean") rather than sorted by current fraction --
// keeps the tab's rows, and the Graph dropdown/legend, in a stable order across generate/
// step calls instead of visibly reshuffling every time the world changes.
const BIOME_NAMES = [
  "Ice", "Tundra", "Boreal Forest", "Temperate Desert", "Temperate Grassland",
  "Woodland/Shrubland", "Temperate Seasonal Forest", "Temperate Rainforest",
  "Subtropical Desert", "Savanna", "Tropical Seasonal Forest", "Tropical Rainforest",
];

const TAB_METRICS: Record<TabKey, Metric[]> = {
  physical: [
    pctMetric("land_fraction", "Land", (s) => s.land_fraction),
    pctMetric("ocean_fraction", "Water", (s) => s.ocean_fraction),
    numMetric("elevation_min_m", "Elevation min (land)", (s) => s.elevation_min_m, 0, " m"),
    numMetric("elevation_mean_m", "Elevation avg (land)", (s) => s.elevation_mean_m, 0, " m"),
    numMetric("elevation_max_m", "Elevation max (land)", (s) => s.elevation_max_m, 0, " m"),
    numMetric("ocean_depth_min_m", "Ocean depth min", (s) => s.ocean_depth_min_m, 0, " m"),
    numMetric("ocean_depth_mean_m", "Ocean depth avg", (s) => s.ocean_depth_mean_m, 0, " m"),
    numMetric("ocean_depth_max_m", "Ocean depth max", (s) => s.ocean_depth_max_m, 0, " m"),
  ],
  biome: BIOME_NAMES.map((name) =>
    pctMetric(name, name, (s) => (name in s.biome_land_fraction ? s.biome_land_fraction[name] : null)),
  ),
  temperature: [
    numMetric("land_temperature_min_c", "Land temp min", (s) => s.land_temperature_min_c, 1, "°C"),
    numMetric("land_temperature_mean_c", "Land temp avg", (s) => s.land_temperature_mean_c, 1, "°C"),
    numMetric("land_temperature_max_c", "Land temp max", (s) => s.land_temperature_max_c, 1, "°C"),
    numMetric("air_temperature_min_c", "Air temp min", (s) => s.air_temperature_min_c, 1, "°C"),
    numMetric("air_temperature_mean_c", "Air temp avg", (s) => s.air_temperature_mean_c, 1, "°C"),
    numMetric("air_temperature_max_c", "Air temp max", (s) => s.air_temperature_max_c, 1, "°C"),
    numMetric("ocean_temperature_min_c", "Ocean temp min", (s) => s.ocean_temperature_min_c, 1, "°C"),
    numMetric("ocean_temperature_mean_c", "Ocean temp avg", (s) => s.ocean_temperature_mean_c, 1, "°C"),
    numMetric("ocean_temperature_max_c", "Ocean temp max", (s) => s.ocean_temperature_max_c, 1, "°C"),
  ],
  precipitation: [
    numMetric("precipitation_min_mm", "Precipitation min", (s) => s.precipitation_min_mm, 0, " mm/yr"),
    numMetric("precipitation_mean_mm", "Precipitation avg", (s) => s.precipitation_mean_mm, 0, " mm/yr"),
    numMetric("precipitation_max_mm", "Precipitation max", (s) => s.precipitation_max_mm, 0, " mm/yr"),
  ],
};

// Unlike every other tab's metrics (a spatial min/max/mean snapshot of the *current* world,
// straight off `current`), these two are single running totals with no per-call distribution
// of their own -- see api.ts's WorldStats docstring. Reused for this tab's Graph mode (the
// raw series over time, same dropdown+chart pattern as every other tab); Table mode instead
// runs runHistoryStats over `history` to get an actual min/max/mean/std-dev out of them.
const SIMULATION_METRICS: Metric[] = [
  numMetric("elevation_point_count", "Elevation points", (s) => s.elevation_point_count, 0),
  numMetric("plate_count", "Plates", (s) => s.plate_count, 0),
];

interface HistoryStats {
  min: number;
  max: number;
  mean: number;
  stdDev: number;
}

function historyStats(history: WorldStats[], get: (s: WorldStats) => number): HistoryStats | null {
  if (history.length === 0) return null;
  const values = history.map(get);
  const min = Math.min(...values);
  const max = Math.max(...values);
  const mean = values.reduce((a, b) => a + b, 0) / values.length;
  const variance = values.reduce((a, b) => a + (b - mean) ** 2, 0) / values.length;
  return { min, max, mean, stdDev: Math.sqrt(variance) };
}

// Single accent for every series -- only one metric is ever plotted at a time (picked via
// the dropdown below), so identity comes from the dropdown/heading, not the line color.
const ACCENT_COLOR = "#4f9dff";

function Row({ label, value }: { label: string; value: string }) {
  return (
    <div style={{ display: "flex", justifyContent: "space-between", gap: 12, padding: "4px 0" }}>
      <span style={{ opacity: 0.7 }}>{label}</span>
      <span>{value}</span>
    </div>
  );
}

function TabButton({ label, active, onClick }: { label: string; active: boolean; onClick: () => void }) {
  return (
    <button
      type="button"
      onClick={onClick}
      style={{
        padding: "6px 10px",
        fontSize: 12,
        background: active ? "#1b2340" : "none",
        border: "none",
        borderBottom: active ? "2px solid #4f9dff" : "2px solid transparent",
        color: active ? "#e6e8ef" : "#8b8fa3",
        cursor: "pointer",
        whiteSpace: "nowrap",
      }}
    >
      {label}
    </button>
  );
}

function ViewModeToggle({ viewMode, onChange }: { viewMode: "table" | "graph"; onChange: (v: "table" | "graph") => void }) {
  return (
    <div style={{ display: "flex", gap: 6, marginBottom: 12 }}>
      {(["table", "graph"] as const).map((mode) => (
        <button
          key={mode}
          type="button"
          onClick={() => onChange(mode)}
          style={{
            padding: "4px 10px",
            fontSize: 11,
            background: viewMode === mode ? "#2a3050" : "none",
            border: "1px solid #333",
            borderRadius: 4,
            color: "#e6e8ef",
            cursor: "pointer",
          }}
        >
          {mode === "table" ? "Table" : "Graph"}
        </button>
      ))}
    </div>
  );
}

function MetricTab({ metrics, history, current }: { metrics: Metric[]; history: WorldStats[]; current: WorldStats }) {
  const [viewMode, setViewMode] = useState<"table" | "graph">("table");
  const [selectedKey, setSelectedKey] = useState(metrics[0].key);
  const selected = metrics.find((m) => m.key === selectedKey) ?? metrics[0];

  const chartData: ChartPoint[] = useMemo(
    () => history.map((h) => ({ x: h.elapsed_years, values: { [selected.key]: selected.get(h) } })),
    [history, selected],
  );
  const chartSeries: ChartSeries[] = useMemo(() => [{ key: selected.key, label: selected.label, color: ACCENT_COLOR }], [selected]);

  return (
    <>
      <ViewModeToggle viewMode={viewMode} onChange={setViewMode} />
      {viewMode === "table" ? (
        metrics.map((m) => <Row key={m.key} label={m.label} value={m.tableFormat(m.get(current))} />)
      ) : (
        <>
          <select
            value={selectedKey}
            onChange={(e) => setSelectedKey(e.target.value)}
            style={{ width: "100%", padding: "5px 4px", marginBottom: 10, fontSize: 12 }}
          >
            {metrics.map((m) => (
              <option key={m.key} value={m.key}>
                {m.label}
              </option>
            ))}
          </select>
          <TimeSeriesChart series={chartSeries} data={chartData} yFormat={selected.yFormat} />
        </>
      )}
    </>
  );
}

// `current` isn't enough on its own here (unlike MetricTab's tabs, "min/max/avg/std dev"
// for these two metrics means over the *run*, not a snapshot) -- Table mode is a dedicated
// historyStats summary per metric instead of MetricTab's plain current-value rows; Graph
// mode reuses MetricTab's own dropdown+chart exactly, since the raw series over time is the
// same shape for these metrics as for every other tab's.
function SimulationTab({ history }: { history: WorldStats[] }) {
  const [viewMode, setViewMode] = useState<"table" | "graph">("table");
  const [selectedKey, setSelectedKey] = useState(SIMULATION_METRICS[0].key);
  const selected = SIMULATION_METRICS.find((m) => m.key === selectedKey) ?? SIMULATION_METRICS[0];

  const chartData: ChartPoint[] = useMemo(
    () => history.map((h) => ({ x: h.elapsed_years, values: { [selected.key]: selected.get(h) } })),
    [history, selected],
  );
  const chartSeries: ChartSeries[] = useMemo(() => [{ key: selected.key, label: selected.label, color: ACCENT_COLOR }], [selected]);

  return (
    <>
      <ViewModeToggle viewMode={viewMode} onChange={setViewMode} />
      {viewMode === "table" ? (
        history.length === 0 ? (
          <div style={{ opacity: 0.6 }}>No history yet -- step the world to start tracking.</div>
        ) : (
          SIMULATION_METRICS.map((m) => {
            const s = historyStats(history, (h) => m.get(h) ?? 0);
            return (
              <div key={m.key} style={{ marginBottom: 14 }}>
                <div style={{ fontWeight: 600, marginBottom: 4 }}>{m.label}</div>
                <Row label="Min" value={m.tableFormat(s?.min ?? null)} />
                <Row label="Max" value={m.tableFormat(s?.max ?? null)} />
                <Row label="Avg" value={m.tableFormat(s?.mean ?? null)} />
                <Row label="Std dev" value={m.tableFormat(s?.stdDev ?? null)} />
              </div>
            );
          })
        )
      ) : (
        <>
          <select
            value={selectedKey}
            onChange={(e) => setSelectedKey(e.target.value)}
            style={{ width: "100%", padding: "5px 4px", marginBottom: 10, fontSize: 12 }}
          >
            {SIMULATION_METRICS.map((m) => (
              <option key={m.key} value={m.key}>
                {m.label}
              </option>
            ))}
          </select>
          <TimeSeriesChart series={chartSeries} data={chartData} yFormat={selected.yFormat} />
        </>
      )}
    </>
  );
}

export default function StatsModal({ stats, history, onClose }: Props) {
  const [activeTab, setActiveTab] = useState<TabKey>("physical");

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
          width: 520,
          maxHeight: "80vh",
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
          <span style={{ fontSize: 16, fontWeight: 600 }}>World Stats</span>
          <button
            type="button"
            onClick={onClose}
            style={{ background: "none", border: "none", color: "#8b8fa3", cursor: "pointer", fontSize: 16 }}
          >
            ✕
          </button>
        </div>

        {!stats && <div style={{ opacity: 0.6 }}>No world yet.</div>}

        {stats && (
          <>
            <div style={{ display: "flex", flexWrap: "wrap", gap: 2, borderBottom: "1px solid #333", marginBottom: 14 }}>
              {TABS.map((t) => (
                <TabButton key={t.key} label={t.label} active={activeTab === t.key} onClick={() => setActiveTab(t.key)} />
              ))}
            </div>
            {/* key={activeTab} forces a remount on tab change -- without it, MetricTab's
                internal viewMode/selectedKey state (useState has no way to know its props'
                *meaning* changed) would carry over from the previous tab, e.g. leaving
                Graph mode on and an arbitrary fallback metric selected right after
                switching to a tab whose metrics don't include the old selection. */}
            {activeTab === "simulation" ? (
              <SimulationTab history={history} />
            ) : (
              <MetricTab key={activeTab} metrics={TAB_METRICS[activeTab]} history={history} current={stats} />
            )}
          </>
        )}
      </div>
    </div>
  );
}
