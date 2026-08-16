import { useCallback, useEffect, useRef, useState } from "react";
import "./index.css";
import { generateWorld, renderWorld, stepWorld } from "./api";
import type { Projection, RenderResponse, WorldSummary } from "./api";
import MapCanvas from "./MapCanvas";

const CANVAS_WIDTH = 900;
const CANVAS_HEIGHT = 500;
const DEFAULT_STEP_YEARS = 2_000_000;
const PLAY_INTERVAL_MS = 400;

export default function App() {
  const [seed, setSeed] = useState(1);
  const [numPlates, setNumPlates] = useState(12);
  const [stepYears, setStepYears] = useState(DEFAULT_STEP_YEARS);
  const [projection, setProjection] = useState<Projection>("behrmann");
  const [summary, setSummary] = useState<WorldSummary | null>(null);
  const [renderData, setRenderData] = useState<RenderResponse | null>(null);
  const [playing, setPlaying] = useState(false);
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const refresh = useCallback(
    async (proj: Projection) => {
      try {
        const data = await renderWorld(proj);
        setRenderData(data);
      } catch (e) {
        setError(String(e));
      }
    },
    [],
  );

  const handleGenerate = useCallback(async () => {
    setBusy(true);
    setError(null);
    try {
      const s = await generateWorld(seed, numPlates);
      setSummary(s);
      await refresh(projection);
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }, [seed, numPlates, projection, refresh]);

  const handleStep = useCallback(async () => {
    if (!summary) return;
    setBusy(true);
    setError(null);
    try {
      const s = await stepWorld(stepYears);
      setSummary(s);
      await refresh(projection);
    } catch (e) {
      setError(String(e));
      setPlaying(false);
    } finally {
      setBusy(false);
    }
  }, [summary, stepYears, projection, refresh]);

  // Re-render with the current world whenever the projection changes.
  useEffect(() => {
    if (summary) {
      refresh(projection);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [projection]);

  const stepRef = useRef(handleStep);
  stepRef.current = handleStep;

  useEffect(() => {
    if (!playing) return;
    const id = setInterval(() => {
      stepRef.current();
    }, PLAY_INTERVAL_MS);
    return () => clearInterval(id);
  }, [playing]);

  return (
    <div style={{ fontFamily: "system-ui, sans-serif", color: "#e6e8ef", padding: 24 }}>
      <h1 style={{ fontSize: 20, marginBottom: 4 }}>mantle-bloom</h1>
      <p style={{ opacity: 0.7, marginTop: 0, marginBottom: 16 }}>
        sphere-native plate tectonics -- elevation view
      </p>

      <div style={{ display: "flex", gap: 24, alignItems: "flex-start" }}>
        <div style={{ display: "flex", flexDirection: "column", gap: 10, minWidth: 220 }}>
          <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 10 }}>
            <legend>Generate</legend>
            <label style={{ display: "block", marginBottom: 6 }}>
              Seed
              <input
                type="number"
                value={seed}
                onChange={(e) => setSeed(Number(e.target.value))}
                style={{ width: "100%" }}
              />
            </label>
            <label style={{ display: "block", marginBottom: 6 }}>
              Plates
              <input
                type="number"
                min={2}
                max={40}
                value={numPlates}
                onChange={(e) => setNumPlates(Number(e.target.value))}
                style={{ width: "100%" }}
              />
            </label>
            <button onClick={handleGenerate} disabled={busy} style={{ width: "100%" }}>
              Generate World
            </button>
          </fieldset>

          <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 10 }}>
            <legend>Time</legend>
            <label style={{ display: "block", marginBottom: 6 }}>
              Years per step
              <input
                type="number"
                min={1000}
                value={stepYears}
                onChange={(e) => setStepYears(Number(e.target.value))}
                style={{ width: "100%" }}
              />
            </label>
            <div style={{ display: "flex", gap: 8 }}>
              <button onClick={handleStep} disabled={busy || !summary} style={{ flex: 1 }}>
                Step
              </button>
              <button
                onClick={() => setPlaying((p) => !p)}
                disabled={!summary}
                style={{ flex: 1 }}
              >
                {playing ? "Pause" : "Play"}
              </button>
            </div>
          </fieldset>

          <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 10 }}>
            <legend>Projection</legend>
            <select
              value={projection}
              onChange={(e) => setProjection(e.target.value as Projection)}
              style={{ width: "100%" }}
            >
              <option value="behrmann">Behrmann (cylindrical equal-area)</option>
              <option value="eckert4">Eckert IV (pseudocylindrical equal-area)</option>
            </select>
          </fieldset>

          {summary && (
            <div style={{ fontSize: 13, opacity: 0.8 }}>
              <div>seed: {summary.seed}</div>
              <div>plates: {summary.num_plates}</div>
              <div>elapsed: {(summary.elapsed_years / 1e6).toFixed(1)} Myr</div>
            </div>
          )}
          {error && <div style={{ color: "#ff8080", fontSize: 13 }}>{error}</div>}
        </div>

        <MapCanvas data={renderData} width={CANVAS_WIDTH} height={CANVAS_HEIGHT} />
      </div>
    </div>
  );
}
