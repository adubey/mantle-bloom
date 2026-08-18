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

type TabKey = "physical" | "temperature" | "precipitation";

const TABS: { key: TabKey; label: string }[] = [
  { key: "physical", label: "Physical" },
  { key: "temperature", label: "Temperature" },
  { key: "precipitation", label: "Precipitation" },
];

const TAB_METRICS: Record<TabKey, Metric[]> = {
  physical: [
    pctMetric("land_fraction", "Land", (s) => s.land_fraction),
    pctMetric("ocean_fraction", "Water", (s) => s.ocean_fraction),
    numMetric("elevation_min_m", "Elevation min", (s) => s.elevation_min_m, 0, " m"),
    numMetric("elevation_mean_m", "Elevation avg", (s) => s.elevation_mean_m, 0, " m"),
    numMetric("elevation_max_m", "Elevation max", (s) => s.elevation_max_m, 0, " m"),
  ],
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
            <MetricTab key={activeTab} metrics={TAB_METRICS[activeTab]} history={history} current={stats} />
          </>
        )}
      </div>
    </div>
  );
}
