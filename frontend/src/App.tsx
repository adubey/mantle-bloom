import { useCallback, useEffect, useRef, useState } from "react";
import "./index.css";
import { generateWorld, renderWorld, stepWorld } from "./api";
import type { Projection, RenderResponse, WorldSummary } from "./api";
import MapCanvas from "./MapCanvas";
import type { MapView } from "./MapCanvas";
import EventConsole from "./EventConsole";

const CANVAS_WIDTH = 1100;
const CANVAS_HEIGHT = 611;
const STEP_YEARS_OPTIONS = [10_000, 100_000, 1_000_000, 10_000_000];
const PLAY_INTERVAL_MS = 400;
// Matches backend app/plates.py's MIN_CONTINENTS/MAX_CONTINENTS.
const MIN_CONTINENTS = 1;
const MAX_CONTINENTS = 12;
const DEFAULT_CONTINENTS = 7;

function randomSeed(): number {
  return Math.floor(Math.random() * 1_000_000_000);
}

export default function App() {
  const [showGenerateDialog, setShowGenerateDialog] = useState(false);
  const [seed, setSeed] = useState(randomSeed());
  const [numContinents, setNumContinents] = useState(DEFAULT_CONTINENTS);

  const [stepYears, setStepYears] = useState(STEP_YEARS_OPTIONS[1]);
  const [projection, setProjection] = useState<Projection>("behrmann");
  const [mapView, setMapView] = useState<MapView>("elevation");
  const [summary, setSummary] = useState<WorldSummary | null>(null);
  const [renderData, setRenderData] = useState<RenderResponse | null>(null);
  const [playing, setPlaying] = useState(false);
  const [busy, setBusy] = useState(false);
  const [stepping, setStepping] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // include_lines is only requested for "Plates (details)" -- the raw per-plate node data
  // it carries (~2MB of JSON) goes unused by every other view, which paint entirely from
  // the grid (see api.ts's renderWorld).
  const refresh = useCallback(async (proj: Projection, view: MapView) => {
    try {
      const data = await renderWorld(proj, view === "platesDetail");
      setRenderData(data);
    } catch (e) {
      setError(String(e));
    }
  }, []);

  const handleGenerate = useCallback(async () => {
    setBusy(true);
    setError(null);
    try {
      const s = await generateWorld(seed, numContinents);
      setSummary(s);
      setShowGenerateDialog(false);
      await refresh(projection, mapView);
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }, [seed, numContinents, projection, mapView, refresh]);

  const handleStep = useCallback(async () => {
    if (!summary) return;
    setStepping(true);
    setError(null);
    try {
      const s = await stepWorld(stepYears);
      setSummary(s);
      await refresh(projection, mapView);
    } catch (e) {
      setError(String(e));
      setPlaying(false);
    } finally {
      setStepping(false);
    }
  }, [summary, stepYears, projection, mapView, refresh]);

  // Re-render with the current world whenever the projection or map view changes (the view
  // matters because it decides whether include_lines is requested, above).
  useEffect(() => {
    if (summary) {
      refresh(projection, mapView);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [projection, mapView]);

  const stepRef = useRef(handleStep);
  stepRef.current = handleStep;

  // Self-scheduling rather than setInterval: each step must finish before the next is
  // scheduled, so a slow step (larger worlds take longer -- a step that also lands on a
  // garbage-collection/gap-fill interval, or triggers a merge/split, can take noticeably
  // longer than a routine one) can never overlap with the next tick. A fixed-cadence
  // setInterval would keep firing regardless of whether the previous request had returned.
  useEffect(() => {
    if (!playing) return;
    let cancelled = false;
    let timeoutId: ReturnType<typeof setTimeout>;
    const tick = () => {
      timeoutId = setTimeout(async () => {
        if (cancelled) return;
        await stepRef.current();
        if (!cancelled) tick();
      }, PLAY_INTERVAL_MS);
    };
    tick();
    return () => {
      cancelled = true;
      clearTimeout(timeoutId);
    };
  }, [playing]);

  return (
    <div style={{ fontFamily: "system-ui, sans-serif", color: "#e6e8ef", padding: 24 }}>
      <h1 style={{ fontSize: 20, marginBottom: 4 }}>mantle-bloom</h1>
      <p style={{ opacity: 0.7, marginTop: 0, marginBottom: 16 }}>
        sphere-native plate tectonics
      </p>

      <div style={{ display: "flex", gap: 16, alignItems: "flex-start" }}>
        <div style={{ display: "flex", flexDirection: "column", gap: 8, width: 170, flexShrink: 0 }}>
          <button onClick={() => setShowGenerateDialog(true)} disabled={busy} style={{ fontSize: 12 }}>
            Generate World
          </button>

          <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
            <legend style={{ fontSize: 11 }}>Time</legend>
            <label style={{ display: "block", marginBottom: 6 }}>
              Years per step
              <select
                value={stepYears}
                onChange={(e) => setStepYears(Number(e.target.value))}
                style={{ width: "100%", fontSize: 12 }}
              >
                {STEP_YEARS_OPTIONS.map((y) => (
                  <option key={y} value={y}>
                    {y.toLocaleString()}
                  </option>
                ))}
              </select>
            </label>
            <div style={{ display: "flex", gap: 6 }}>
              <button onClick={handleStep} disabled={busy || stepping || !summary} style={{ flex: 1, fontSize: 12 }}>
                Step
              </button>
              <button onClick={() => setPlaying((p) => !p)} disabled={busy || !summary} style={{ flex: 1, fontSize: 12 }}>
                {playing ? "Pause" : "Play"}
              </button>
            </div>
          </fieldset>

          <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
            <legend style={{ fontSize: 11 }}>Map View</legend>
            <select
              value={mapView}
              onChange={(e) => setMapView(e.target.value as MapView)}
              style={{ width: "100%", marginBottom: 6, fontSize: 12 }}
            >
              <option value="plates">Plates</option>
              <option value="platesDetail">Plates (details)</option>
              <option value="elevation">Elevation</option>
            </select>
            <select
              value={projection}
              onChange={(e) => setProjection(e.target.value as Projection)}
              style={{ width: "100%", fontSize: 12 }}
            >
              <option value="behrmann">Behrmann (cylindrical equal-area)</option>
              <option value="eckert4">Eckert IV (pseudocylindrical equal-area)</option>
            </select>
          </fieldset>

          {summary && (
            <div style={{ fontSize: 11, opacity: 0.8 }}>
              <div>seed: {summary.seed}</div>
              <div>plates: {summary.num_plates}</div>
              <div>elapsed: {(summary.elapsed_years / 1e6).toFixed(1)} Myr</div>
            </div>
          )}
          {error && <div style={{ color: "#ff8080", fontSize: 11 }}>{error}</div>}

          <EventConsole events={summary?.events ?? []} />
        </div>

        <div style={{ marginTop: -50 }}>
          <MapCanvas data={renderData} view={mapView} width={CANVAS_WIDTH} height={CANVAS_HEIGHT} />
        </div>
      </div>

      {showGenerateDialog && (
        <div
          style={{
            position: "fixed",
            inset: 0,
            background: "rgba(0, 0, 0, 0.6)",
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
          }}
        >
          <div
            style={{
              background: "#151a2e",
              border: "1px solid #333",
              borderRadius: 8,
              padding: 20,
              minWidth: 280,
            }}
          >
            <h2 style={{ fontSize: 16, marginTop: 0 }}>Generate World</h2>

            <label style={{ display: "block", marginBottom: 16 }}>
              Seed
              <div style={{ display: "flex", gap: 6 }}>
                <input
                  type="number"
                  value={seed}
                  onChange={(e) => setSeed(Number(e.target.value))}
                  style={{ flex: 1 }}
                />
                <button
                  type="button"
                  title="Randomize seed"
                  aria-label="Randomize seed"
                  onClick={() => setSeed(randomSeed())}
                >
                  🎲
                </button>
              </div>
            </label>

            <label style={{ display: "block", marginBottom: 16 }}>
              Continents: {numContinents}
              <input
                type="range"
                min={MIN_CONTINENTS}
                max={MAX_CONTINENTS}
                value={numContinents}
                onChange={(e) => setNumContinents(Number(e.target.value))}
                style={{ width: "100%" }}
              />
            </label>

            <div style={{ display: "flex", gap: 8, justifyContent: "flex-end" }}>
              <button onClick={() => setShowGenerateDialog(false)} disabled={busy}>
                Cancel
              </button>
              <button onClick={handleGenerate} disabled={busy}>
                Generate
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}
