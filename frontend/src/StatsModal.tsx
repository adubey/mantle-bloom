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
  // Set only for a metric that can go stale (currently land_fraction/ocean_fraction -- see
  // WorldStats.land_fraction_stale, GitHub issue #121): Table mode's current-value row shows
  // a marker when this returns true, so a frozen hydrology cache doesn't read as a normal
  // live number.
  stale?: (s: WorldStats) => boolean;
}

function pctMetric(
  key: string,
  label: string,
  get: (s: WorldStats) => number | null,
  stale?: (s: WorldStats) => boolean,
): Metric {
  return {
    key,
    label,
    get,
    tableFormat: (v) => (v === null ? "--" : `${(v * 100).toFixed(1)}%`),
    yFormat: (v) => `${(v * 100).toFixed(0)}%`,
    stale,
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

const TAB_METRICS: Record<Exclude<TabKey, "simulation" | "biome">, TabEntry[]> = {
  physical: [
    metricEntry(pctMetric("land_fraction", "Land", (s) => s.land_fraction, (s) => s.land_fraction_stale ?? false)),
    metricEntry(pctMetric("ocean_fraction", "Water", (s) => s.ocean_fraction, (s) => s.land_fraction_stale ?? false)),
    // Raw elevation>sea_level node fraction -- always fresh, unlike Land/Water above (see
    // WorldStats.land_fraction_node's own comment, GitHub issue #121).
    metricEntry(pctMetric("land_fraction_node", "Land (raw)", (s) => s.land_fraction_node ?? null)),
    metricEntry(pctMetric(
      "land_near_max_elevation_fraction", "Land within 5% of max elevation",
      (s) => s.land_near_max_elevation_fraction ?? null,
    )),
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

// Hc/Hm are spatial min/max/mean/std-dev snapshots of the *current* world, same shape as
// elevation/temperature/precipitation on the other tabs -- grouped the same way (one graph,
// avg/min/max lines plus a std-dev band) rather than four separately-plottable lines. The
// remaining entries here are single running totals with no per-call distribution of their
// own -- see api.ts's WorldStats docstring -- so SimulationTab still runs historyStats over
// `history` for their Table mode, unlike a StatGroup's own current-snapshot table.
const SIMULATION_ENTRIES: TabEntry[] = [
  metricEntry(numMetric(
    "hc_at_max_fraction", "% of Hc at max",
    (s) => s.hc_at_max_fraction == null ? null : 100 * s.hc_at_max_fraction, 2, "%",
  )),
  groupEntry(numGroup(
    "hc_m", "Hc (crustal thickness)",
    (s) => s.hc_min_m ?? null, (s) => s.hc_mean_m ?? null, (s) => s.hc_max_m ?? null, (s) => s.hc_std_m ?? null,
    0, " m",
  )),
  metricEntry(numMetric(
    "hm_at_max_fraction", "% of Hm at max",
    (s) => s.hm_at_max_fraction == null ? null : 100 * s.hm_at_max_fraction, 2, "%",
  )),
  groupEntry(numGroup(
    "hm_m", "Hm (mantle lithosphere thickness)",
    (s) => s.hm_min_m ?? null, (s) => s.hm_mean_m ?? null, (s) => s.hm_max_m ?? null, (s) => s.hm_std_m ?? null,
    0, " m",
  )),
  metricEntry(numMetric("elevation_point_count", "Elevation points", (s) => s.elevation_point_count, 0)),
  metricEntry(numMetric("plate_count", "Plates", (s) => s.plate_count, 0)),
  // Land/volume-conservation check (see api.ts's own comment): land area is read straight off
  // world.plates -- immune to the climate-grid land_fraction's hydrology-cache staleness --
  // and crust volume is continental crust only, so a real mass-conservation bug (a topology
  // change that drops a column's volume instead of preserving it) shows up as this number
  // actually shrinking, distinct from land area swinging while volume holds roughly steady
  // (ordinary tectonics moving existing crust above/below sea level).
  metricEntry(numMetric("total_land_area_km2", "Land area", (s) => s.total_land_area_km2, 0, " km²")),
  metricEntry(numMetric(
    "total_continental_crust_volume_km3", "Continental crust volume",
    (s) => s.total_continental_crust_volume_km3, 0, " km³",
  )),
  metricEntry(numMetric("sea_level_m", "Sea level", (s) => s.sea_level_m, 0, " m")),
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

// Categorical line colors for the multi-series Biome charts, assigned in fixed order (a
// series keeps its color regardless of how many others are shown) -- the dark-surface steps
// of the shared data-viz categorical palette, validated for CVD separation against this
// modal's #151a2e surface. No Biome view defines more than 5 series, so the ramp is never
// cycled (a cycled hue would collide); a view that ever needs a 9th series should group instead.
const CATEGORICAL_COLORS = [
  "#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#9085e9", "#008300", "#e66767",
];

// The Biome tab plots several biome fractions on one chart at once; the dropdown picks which
// set. Each `series` value at a given step is the summed fraction of its `members` (entries
// of backend biomes.BIOME_NAMES) -- over land cells for a "land" view, ocean cells for an
// "ocean" view (biome_land_fraction / biome_ocean_fraction respectively). Member name lists
// are hand-synced with backend biomes.KOPPEN_NAMES / PELAGIC_NAMES, same precedent as
// legendData.ts's BIOME_RGB_ENTRIES. The "groups" and "ocean" views also back Table mode's
// full per-class breakdown, so their `members` stay the complete class lists.
type BiomeDomain = "land" | "ocean";

interface BiomeSeriesDef {
  key: string;
  label: string;
  members: string[];
}

interface BiomeView {
  key: string;
  label: string;
  domain: BiomeDomain;
  series: BiomeSeriesDef[];
}

const KOPPEN_A = ["Tropical Rainforest", "Tropical Monsoon", "Tropical Savanna", "Tropical Savanna (Dry Summer)"];
const KOPPEN_B = ["Hot Desert", "Cold Desert", "Hot Semi-Arid", "Cold Semi-Arid"];
const KOPPEN_C = [
  "Hot-Summer Mediterranean", "Warm-Summer Mediterranean", "Cold-Summer Mediterranean",
  "Humid Subtropical (Dry Winter)", "Subtropical Highland", "Cold Subtropical Highland",
  "Humid Subtropical", "Oceanic", "Subpolar Oceanic",
];
const KOPPEN_D = [
  "Mediterranean Continental (Hot Summer)", "Mediterranean Continental (Warm Summer)",
  "Mediterranean Subarctic", "Extremely Cold Mediterranean Subarctic",
  "Monsoon Continental (Hot Summer)", "Monsoon Continental (Warm Summer)",
  "Monsoon Subarctic", "Extremely Cold Monsoon Subarctic",
  "Hot-Summer Humid Continental", "Warm-Summer Humid Continental",
  "Subarctic (Boreal)", "Extremely Cold Subarctic",
];
const KOPPEN_E = ["Tundra", "Ice Cap"];

const one = (name: string): BiomeSeriesDef => ({ key: name, label: name, members: [name] });

const BIOME_VIEWS: BiomeView[] = [
  {
    key: "groups",
    label: "Climate groups (Köppen A–E)",
    domain: "land",
    series: [
      { key: "A", label: "Tropical (A)", members: KOPPEN_A },
      { key: "B", label: "Arid (B)", members: KOPPEN_B },
      { key: "C", label: "Temperate (C)", members: KOPPEN_C },
      { key: "D", label: "Continental (D)", members: KOPPEN_D },
      { key: "E", label: "Polar (E)", members: KOPPEN_E },
    ],
  },
  {
    key: "ocean",
    label: "Oceanic provinces",
    domain: "ocean",
    series: [
      {
        key: "tropical",
        label: "Tropical seas",
        members: ["Tropical Open Ocean", "Subtropical Gyre", "Equatorial Divergence", "Tropical Coastal Waters"],
      },
      { key: "temperate", label: "Temperate seas", members: ["Temperate Open Ocean", "Temperate Shelf"] },
      { key: "coldTemperate", label: "Cold-temperate seas", members: ["Cold-Temperate Open Ocean", "Cold-Temperate Shelf"] },
      { key: "polar", label: "Polar ocean", members: ["Polar Ocean"] },
      { key: "seaIce", label: "Sea ice", members: ["Polar Sea Ice"] },
    ],
  },
  { key: "tropical", label: "Tropical zone (A) subtypes", domain: "land", series: KOPPEN_A.map(one) },
  { key: "arid", label: "Arid zone (B) subtypes", domain: "land", series: KOPPEN_B.map(one) },
  {
    key: "temperate",
    label: "Temperate zone (C) subtypes",
    domain: "land",
    series: [
      { key: "Cs", label: "Mediterranean (dry summer)", members: KOPPEN_C.slice(0, 3) },
      { key: "Cw", label: "Dry-winter temperate", members: KOPPEN_C.slice(3, 6) },
      { key: "Cf", label: "Year-round-wet temperate", members: KOPPEN_C.slice(6, 9) },
    ],
  },
  {
    key: "continental",
    label: "Continental zone (D) subtypes",
    domain: "land",
    series: [
      { key: "Ds", label: "Dry-summer continental", members: KOPPEN_D.slice(0, 4) },
      { key: "Dw", label: "Dry-winter continental", members: KOPPEN_D.slice(4, 8) },
      { key: "Df", label: "Year-round-wet continental", members: KOPPEN_D.slice(8, 12) },
    ],
  },
  { key: "polar", label: "Polar zone (E) subtypes", domain: "land", series: KOPPEN_E.map(one) },
];

function biomeDomainDict(s: WorldStats, domain: BiomeDomain): Record<string, number> {
  return domain === "land" ? s.biome_land_fraction : s.biome_ocean_fraction ?? {};
}

// SimulationTab's Table mode for a plain Metric entry: a run-level historyStats summary (see
// that component's own comment for why, as opposed to a StatGroup's current-snapshot rows).
function MetricHistoryRow({ metric, history }: { metric: Metric; history: WorldStats[] }) {
  const s = historyStats(history, (h) => metric.get(h) ?? 0);
  return (
    <div style={{ marginBottom: 14 }}>
      <div style={{ fontWeight: 600, marginBottom: 4 }}>{metric.label}</div>
      <Row label="Min" value={metric.tableFormat(s?.min ?? null)} />
      <Row label="Max" value={metric.tableFormat(s?.max ?? null)} />
      <Row label="Avg" value={metric.tableFormat(s?.mean ?? null)} />
      <Row label="Std dev" value={metric.tableFormat(s?.stdDev ?? null)} />
    </div>
  );
}

// Summed fraction of `members`, or null when the domain has no cells at all this step (an
// empty dict) so the chart draws a gap rather than a misleading 0.
function biomeSeriesValue(s: WorldStats, domain: BiomeDomain, members: string[]): number | null {
  const dict = biomeDomainDict(s, domain);
  if (Object.keys(dict).length === 0) return null;
  let sum = 0;
  for (const m of members) sum += dict[m] ?? 0;
  return sum;
}

const biomePctFormat = (v: number) => `${(v * 100).toFixed(0)}%`;
const biomePctTableFormat = (v: number | null) => (v === null ? "--" : `${(v * 100).toFixed(1)}%`);

function Row({ label, value, title }: { label: string; value: string; title?: string }) {
  return (
    <div style={{ display: "flex", justifyContent: "space-between", gap: 12, padding: "4px 0" }} title={title}>
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

// Shared by MetricTab and SimulationTab's Graph mode: both plot one selected TabEntry at a
// time, a plain Metric as a single accent-colored line and a StatGroup as avg/min/max lines
// plus a std-dev band.
function entryChartData(selected: TabEntry, history: WorldStats[]): ChartPoint[] {
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
}

// For a StatGroup, "mean" is listed last so TimeSeriesChart paints the avg line on top of
// the min/max lines and the std-dev band (see that component's own z-order comment).
function entryChartSeries(selected: TabEntry): ChartSeries[] {
  if (selected.kind === "metric") return [{ key: selected.metric.key, label: selected.metric.label, color: ACCENT_COLOR }];
  return [
    { key: "min", label: "Min", color: MIN_COLOR },
    { key: "max", label: "Max", color: MAX_COLOR },
    { key: "mean", label: "Avg", color: ACCENT_COLOR },
  ];
}

function entryChartBands(selected: TabEntry): ChartBand[] | undefined {
  return selected.kind === "group" ? [{ lowKey: "bandLow", highKey: "bandHigh", color: BAND_COLOR, label: "±1 std dev" }] : undefined;
}

function entryYFormat(selected: TabEntry): (v: number) => string {
  return selected.kind === "metric" ? selected.metric.yFormat : selected.group.yFormat;
}

function MetricTab({ entries, history, current }: { entries: TabEntry[]; history: WorldStats[]; current: WorldStats }) {
  const [viewMode, setViewMode] = useState<"table" | "graph">("graph");
  const [selectedKey, setSelectedKey] = useState(entryKey(entries[0]));
  const selected = entries.find((e) => entryKey(e) === selectedKey) ?? entries[0];

  const chartData = useMemo(() => entryChartData(selected, history), [history, selected]);
  const chartSeries = useMemo(() => entryChartSeries(selected), [selected]);
  const chartBands = useMemo(() => entryChartBands(selected), [selected]);
  const yFormat = entryYFormat(selected);

  return (
    <>
      <ViewModeToggle viewMode={viewMode} onChange={setViewMode} />
      {viewMode === "table" ? (
        entries.map((e) =>
          e.kind === "metric" ? (
            <Row
              key={e.metric.key}
              label={e.metric.label}
              value={e.metric.tableFormat(e.metric.get(current)) + (e.metric.stale?.(current) ? " *" : "")}
              title={
                e.metric.stale?.(current)
                  ? "Stale: erosion/hydrology hasn't run since Climate was turned off, so this may not match the world's current coastline. See \"Land (raw)\" for a value that's always current."
                  : undefined
              }
            />
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
          {entries.length > 1 && (
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
          )}
          <TimeSeriesChart series={chartSeries} data={chartData} yFormat={yFormat} bands={chartBands} />
        </>
      )}
    </>
  );
}

// Biomes get their own tab (not MetricTab): Graph mode plots a whole *set* of biome
// fractions on one chart -- the dropdown picks the set (top-level Köppen groups, the ocean
// provinces, or one zone's subtypes) rather than a single line at a time. Table mode ignores
// the dropdown and shows the full per-class breakdown for both domains at the current step.
function BiomeTab({ history, current }: { history: WorldStats[]; current: WorldStats }) {
  const [viewMode, setViewMode] = useState<"table" | "graph">("graph");
  const [selectedKey, setSelectedKey] = useState(BIOME_VIEWS[0].key);
  const view = BIOME_VIEWS.find((v) => v.key === selectedKey) ?? BIOME_VIEWS[0];

  const chartData: ChartPoint[] = useMemo(
    () =>
      history.map((h) => ({
        x: h.elapsed_years,
        values: Object.fromEntries(view.series.map((s) => [s.key, biomeSeriesValue(h, view.domain, s.members)])),
      })),
    [history, view],
  );
  const chartSeries: ChartSeries[] = useMemo(
    () => view.series.map((s, i) => ({ key: s.key, label: s.label, color: CATEGORICAL_COLORS[i % CATEGORICAL_COLORS.length] })),
    [view],
  );

  // Full breakdown for both domains -- the "groups" and "ocean" views already carry the
  // complete Köppen / pelagic class lists, so Table mode is just those rendered with a
  // group subtotal and its member classes indented under it.
  const tableSections: { title: string; view: BiomeView }[] = [
    { title: "Land — Köppen classes (% of land)", view: BIOME_VIEWS.find((v) => v.key === "groups")! },
    { title: "Ocean — pelagic provinces (% of ocean)", view: BIOME_VIEWS.find((v) => v.key === "ocean")! },
  ];

  return (
    <>
      <ViewModeToggle viewMode={viewMode} onChange={setViewMode} />
      {viewMode === "table" ? (
        tableSections.map((section) => {
          const dict = biomeDomainDict(current, section.view.domain);
          return (
            <div key={section.title} style={{ marginBottom: 16 }}>
              <div style={{ fontWeight: 600, marginBottom: 6 }}>{section.title}</div>
              {Object.keys(dict).length === 0 ? (
                <div style={{ opacity: 0.6 }}>No {section.view.domain} cells.</div>
              ) : (
                section.view.series.map((group) => (
                  <div key={group.key} style={{ marginBottom: 8 }}>
                    <Row label={group.label} value={biomePctTableFormat(biomeSeriesValue(current, section.view.domain, group.members))} />
                    <div style={{ paddingLeft: 12, opacity: 0.85 }}>
                      {group.members.map((m) => (
                        <Row key={m} label={m} value={biomePctTableFormat(m in dict ? dict[m] : null)} />
                      ))}
                    </div>
                  </div>
                ))
              )}
            </div>
          );
        })
      ) : (
        <>
          <select
            value={selectedKey}
            onChange={(e) => setSelectedKey(e.target.value)}
            style={{ width: "100%", padding: "5px 4px", marginBottom: 10, fontSize: 12 }}
          >
            {BIOME_VIEWS.map((v) => (
              <option key={v.key} value={v.key}>
                {v.label}
              </option>
            ))}
          </select>
          <TimeSeriesChart series={chartSeries} data={chartData} yFormat={biomePctFormat} />
        </>
      )}
    </>
  );
}

// Graph mode reuses the same dropdown+chart pattern (and the entryChart* helpers) as every
// other tab's MetricTab, since the raw series over time is the same shape for a plain Metric
// or a StatGroup here as anywhere else. Table mode differs per entry kind: a StatGroup (Hc,
// Hm) is already a spatial min/max/mean/std-dev snapshot of `current`, so it gets the same
// current-value rows MetricTab's own table would show it; a plain Metric here is a single
// running total with no per-call distribution of its own (see SIMULATION_ENTRIES' comment),
// so it still needs a dedicated historyStats summary over `history` instead.
function SimulationTab({ history, current }: { history: WorldStats[]; current: WorldStats }) {
  const [viewMode, setViewMode] = useState<"table" | "graph">("graph");
  const [selectedKey, setSelectedKey] = useState(entryKey(SIMULATION_ENTRIES[0]));
  const selected = SIMULATION_ENTRIES.find((e) => entryKey(e) === selectedKey) ?? SIMULATION_ENTRIES[0];

  const chartData = useMemo(() => entryChartData(selected, history), [history, selected]);
  const chartSeries = useMemo(() => entryChartSeries(selected), [selected]);
  const chartBands = useMemo(() => entryChartBands(selected), [selected]);
  const yFormat = entryYFormat(selected);

  return (
    <>
      <ViewModeToggle viewMode={viewMode} onChange={setViewMode} />
      {viewMode === "table" ? (
        history.length === 0 ? (
          <div style={{ opacity: 0.6 }}>No history yet -- step the world to start tracking.</div>
        ) : (
          SIMULATION_ENTRIES.map((e) =>
            e.kind === "group" ? (
              <div key={e.group.key} style={{ marginBottom: 14 }}>
                <div style={{ fontWeight: 600, marginBottom: 4 }}>{e.group.label}</div>
                <Row label="Min" value={e.group.tableFormat(e.group.min(current))} />
                <Row label="Avg" value={e.group.tableFormat(e.group.mean(current))} />
                <Row label="Max" value={e.group.tableFormat(e.group.max(current))} />
                <Row label="Std dev" value={e.group.tableFormat(e.group.std(current))} />
              </div>
            ) : (
              <MetricHistoryRow key={e.metric.key} metric={e.metric} history={history} />
            ),
          )
        )
      ) : (
        <>
          <select
            value={selectedKey}
            onChange={(e) => setSelectedKey(e.target.value)}
            style={{ width: "100%", padding: "5px 4px", marginBottom: 10, fontSize: 12 }}
          >
            {SIMULATION_ENTRIES.map((e) => (
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
              <SimulationTab history={history} current={stats} />
            ) : activeTab === "biome" ? (
              <BiomeTab key={activeTab} history={history} current={stats} />
            ) : (
              <MetricTab key={activeTab} entries={TAB_METRICS[activeTab]} history={history} current={stats} />
            )}
          </>
        )}
      </div>
    </div>
  );
}
