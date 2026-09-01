import { useMemo, useState } from "react";
import type { WorldStats } from "./api";
import TimeSeriesChart from "./TimeSeriesChart";
import type { ChartBand, ChartPoint, ChartSeries } from "./TimeSeriesChart";

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

// A backend-computed min/max/mean/std-dev quadruple for one quantity (e.g. elevation), as
// opposed to `Metric`'s single current-value snapshot -- see TAB_METRICS below for which
// stats get this treatment. Table mode lists all four; Graph mode plots avg/min/max as
// separate lines plus a std-dev band, rather than offering min/max/avg/std-dev as four
// separate single-line choices the way each used to be (see MetricTab).
interface StatGroup {
  key: string;
  label: string;
  min: (s: WorldStats) => number | null;
  mean: (s: WorldStats) => number | null;
  max: (s: WorldStats) => number | null;
  std: (s: WorldStats) => number | null;
  tableFormat: (v: number | null) => string;
  yFormat: (v: number) => string;
}

function numGroup(
  key: string,
  label: string,
  min: (s: WorldStats) => number | null,
  mean: (s: WorldStats) => number | null,
  max: (s: WorldStats) => number | null,
  std: (s: WorldStats) => number | null,
  digits = 1,
  suffix = "",
): StatGroup {
  return {
    key,
    label,
    min,
    mean,
    max,
    std,
    tableFormat: (v) => (v === null ? "--" : `${v.toFixed(digits)}${suffix}`),
    yFormat: (v) => `${v.toFixed(digits)}${suffix}`,
  };
}

// A tab's row order mixes plain single-value Metrics (e.g. Land/Water fraction) with
// StatGroups (e.g. Elevation min/avg/max/std-dev) -- this tags each row with which one it is
// so MetricTab can render/plot them differently without a runtime type check on shape.
type TabEntry = { kind: "metric"; metric: Metric } | { kind: "group"; group: StatGroup };

function metricEntry(metric: Metric): TabEntry {
  return { kind: "metric", metric };
}

function groupEntry(group: StatGroup): TabEntry {
  return { kind: "group", group };
}

function entryKey(e: TabEntry): string {
  return e.kind === "metric" ? e.metric.key : e.group.key;
}

function entryLabel(e: TabEntry): string {
  return e.kind === "metric" ? e.metric.label : e.group.label;
}

type TabKey = "physical" | "temperature" | "precipitation" | "biome" | "simulation";

const TABS: { key: TabKey; label: string }[] = [
  { key: "physical", label: "Physical" },
  { key: "temperature", label: "Temperature" },
  { key: "precipitation", label: "Precipitation" },
  { key: "biome", label: "Biome" },
  { key: "simulation", label: "Simulation" },
];

// Fixed order (the Köppen land classes of biomes.BIOME_NAMES -- the pelagic ocean classes are
// excluded from biome_land_fraction) rather than sorted by current fraction, so the tab's
// rows and the Graph dropdown/legend stay stable across generate/step calls. Kept in sync by
// hand with backend biomes.KOPPEN_NAMES, same precedent as legendData.ts's BIOME_RGB_ENTRIES.
const BIOME_NAMES = [
  "Tropical Rainforest", "Tropical Monsoon", "Tropical Savanna", "Tropical Savanna (Dry Summer)",
  "Hot Desert", "Cold Desert", "Hot Semi-Arid", "Cold Semi-Arid",
  "Hot-Summer Mediterranean", "Warm-Summer Mediterranean", "Cold-Summer Mediterranean",
  "Humid Subtropical (Dry Winter)", "Subtropical Highland", "Cold Subtropical Highland",
  "Humid Subtropical", "Oceanic", "Subpolar Oceanic",
  "Mediterranean Continental (Hot Summer)", "Mediterranean Continental (Warm Summer)",
  "Mediterranean Subarctic", "Extremely Cold Mediterranean Subarctic",
  "Monsoon Continental (Hot Summer)", "Monsoon Continental (Warm Summer)",
  "Monsoon Subarctic", "Extremely Cold Monsoon Subarctic",
  "Hot-Summer Humid Continental", "Warm-Summer Humid Continental",
  "Subarctic (Boreal)", "Extremely Cold Subarctic",
  "Tundra", "Ice Cap",
];

const TAB_METRICS: Record<Exclude<TabKey, "simulation">, TabEntry[]> = {
  physical: [
    metricEntry(pctMetric("land_fraction", "Land", (s) => s.land_fraction)),
    metricEntry(pctMetric("ocean_fraction", "Water", (s) => s.ocean_fraction)),
    groupEntry(
      numGroup(
        "elevation_m", "Elevation (land)",
        (s) => s.elevation_min_m, (s) => s.elevation_mean_m, (s) => s.elevation_max_m, (s) => s.elevation_std_m,
        0, " m",
      ),
    ),
    groupEntry(
      numGroup(
        "ocean_depth_m", "Ocean depth",
        (s) => s.ocean_depth_min_m, (s) => s.ocean_depth_mean_m, (s) => s.ocean_depth_max_m, (s) => s.ocean_depth_std_m,
        0, " m",
      ),
    ),
  ],
  biome: BIOME_NAMES.map((name) =>
    metricEntry(pctMetric(name, name, (s) => (name in s.biome_land_fraction ? s.biome_land_fraction[name] : null))),
  ),
  temperature: [
    groupEntry(
      numGroup(
        "land_temperature_c", "Land temp",
        (s) => s.land_temperature_min_c, (s) => s.land_temperature_mean_c, (s) => s.land_temperature_max_c, (s) => s.land_temperature_std_c,
        1, "°C",
      ),
    ),
    groupEntry(
      numGroup(
        "air_temperature_c", "Air temp",
        (s) => s.air_temperature_min_c, (s) => s.air_temperature_mean_c, (s) => s.air_temperature_max_c, (s) => s.air_temperature_std_c,
        1, "°C",
      ),
    ),
    groupEntry(
      numGroup(
        "ocean_temperature_c", "Ocean temp",
        (s) => s.ocean_temperature_min_c, (s) => s.ocean_temperature_mean_c, (s) => s.ocean_temperature_max_c, (s) => s.ocean_temperature_std_c,
        1, "°C",
      ),
    ),
  ],
  precipitation: [
    groupEntry(
      numGroup(
        "precipitation_mm", "Precipitation",
        (s) => s.precipitation_min_mm, (s) => s.precipitation_mean_mm, (s) => s.precipitation_max_mm, (s) => s.precipitation_std_mm,
        0, " mm/yr",
      ),
    ),
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

// One quantity is plotted at a time (picked via the dropdown below), so a plain Metric's
// single line just uses the accent color -- identity comes from the dropdown/heading, not
// the line color. A StatGroup instead plots three lines at once (avg/min/max) plus a std-dev
// band, so those get their own fixed colors so the legend can tell them apart.
const ACCENT_COLOR = "#4f9dff"; // avg
const MIN_COLOR = "#4fd68c";
const MAX_COLOR = "#ff9f4f";
const BAND_COLOR = "rgba(139, 143, 163, 0.28)"; // translucent AXIS_TEXT_COLOR grey

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

function MetricTab({ entries, history, current }: { entries: TabEntry[]; history: WorldStats[]; current: WorldStats }) {
  const [viewMode, setViewMode] = useState<"table" | "graph">("graph");
  const [selectedKey, setSelectedKey] = useState(entryKey(entries[0]));
  const selected = entries.find((e) => entryKey(e) === selectedKey) ?? entries[0];

  const chartData: ChartPoint[] = useMemo(() => {
    if (selected.kind === "metric") {
      const m = selected.metric;
      return history.map((h) => ({ x: h.elapsed_years, values: { [m.key]: m.get(h) } }));
    }
    const g = selected.group;
    return history.map((h) => {
      const mean = g.mean(h);
      const std = g.std(h);
      const bandLow = mean !== null && std !== null ? mean - std : null;
      const bandHigh = mean !== null && std !== null ? mean + std : null;
      return { x: h.elapsed_years, values: { min: g.min(h), max: g.max(h), mean, bandLow, bandHigh } };
    });
  }, [history, selected]);

  // For a StatGroup, "mean" is listed last so TimeSeriesChart paints the avg line on top of
  // the min/max lines and the std-dev band (see that component's own z-order comment).
  const chartSeries: ChartSeries[] = useMemo(() => {
    if (selected.kind === "metric") return [{ key: selected.metric.key, label: selected.metric.label, color: ACCENT_COLOR }];
    return [
      { key: "min", label: "Min", color: MIN_COLOR },
      { key: "max", label: "Max", color: MAX_COLOR },
      { key: "mean", label: "Avg", color: ACCENT_COLOR },
    ];
  }, [selected]);

  const chartBands: ChartBand[] | undefined = useMemo(
    () => (selected.kind === "group" ? [{ lowKey: "bandLow", highKey: "bandHigh", color: BAND_COLOR, label: "±1 std dev" }] : undefined),
    [selected],
  );

  const yFormat = selected.kind === "metric" ? selected.metric.yFormat : selected.group.yFormat;

  return (
    <>
      <ViewModeToggle viewMode={viewMode} onChange={setViewMode} />
      {viewMode === "table" ? (
        entries.map((e) =>
          e.kind === "metric" ? (
            <Row key={e.metric.key} label={e.metric.label} value={e.metric.tableFormat(e.metric.get(current))} />
          ) : (
            <div key={e.group.key} style={{ marginBottom: 14 }}>
              <div style={{ fontWeight: 600, marginBottom: 4 }}>{e.group.label}</div>
              <Row label="Min" value={e.group.tableFormat(e.group.min(current))} />
              <Row label="Avg" value={e.group.tableFormat(e.group.mean(current))} />
              <Row label="Max" value={e.group.tableFormat(e.group.max(current))} />
              <Row label="Std dev" value={e.group.tableFormat(e.group.std(current))} />
            </div>
          ),
        )
      ) : (
        <>
          <select
            value={selectedKey}
            onChange={(e) => setSelectedKey(e.target.value)}
            style={{ width: "100%", padding: "5px 4px", marginBottom: 10, fontSize: 12 }}
          >
            {entries.map((e) => (
              <option key={entryKey(e)} value={entryKey(e)}>
                {entryLabel(e)}
              </option>
            ))}
          </select>
          <TimeSeriesChart series={chartSeries} data={chartData} yFormat={yFormat} bands={chartBands} />
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
  const [viewMode, setViewMode] = useState<"table" | "graph">("graph");
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
              <MetricTab key={activeTab} entries={TAB_METRICS[activeTab]} history={history} current={stats} />
            )}
          </>
        )}
      </div>
    </div>
  );
}
