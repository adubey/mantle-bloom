import { useCallback, useEffect, useRef, useState } from "react";
import "./index.css";
import { fetchPlates, fetchRivers, fetchStats, generateWorld, renderWorld, stepWorld } from "./api";
import type { MapView, PlateSummary, Projection, RenderResponse, RiverSummary, Segment, WorldStats, WorldSummary } from "./api";
import MapCanvas from "./MapCanvas";
import PlateInspector from "./PlateInspector";
import RiverInspector from "./RiverInspector";
import EventConsole from "./EventConsole";
import StatsModal from "./StatsModal";
import { IDENTITY_ROTATION } from "./rotation";
import type { Mat3 } from "./rotation";

// The map is displayed at DISPLAY_WIDTH x DISPLAY_HEIGHT (CSS pixels, unchanged from
// before), but the image requested from the server is RENDER_SCALE times bigger -- a
// sharper, retina-style render at the same on-screen size, rather than a bigger map.
const DISPLAY_WIDTH = 1100;
const DISPLAY_HEIGHT = 611;
const RENDER_SCALE = 2;
const RENDER_WIDTH = DISPLAY_WIDTH * RENDER_SCALE;
const RENDER_HEIGHT = DISPLAY_HEIGHT * RENDER_SCALE;
const STEP_YEARS_OPTIONS = [10_000, 100_000, 1_000_000, 10_000_000];
const PLAY_INTERVAL_MS = 400;
// Percent, matching backend app/plates.py's DEFAULT_CONTINENTAL_FRACTION/DEFAULT_LAND_FRACTION.
const DEFAULT_CONTINENTAL_PERCENT = 70;
const DEFAULT_LAND_PERCENT = 29;
// Degrees, matching backend app/world.py's DEFAULT_AXIAL_TILT_DEG (Earth's real tilt).
const DEFAULT_AXIAL_TILT_DEG = 23.5;

function randomSeed(): number {
  return Math.floor(Math.random() * 1_000_000_000);
}

function formatLatLon(latDeg: number, lonDeg: number): string {
  const latDir = latDeg >= 0 ? "N" : "S";
  const lonDir = lonDeg >= 0 ? "E" : "W";
  return `${Math.abs(latDeg).toFixed(1)}°${latDir}, ${Math.abs(lonDeg).toFixed(1)}°${lonDir}`;
}

function isIdentityRotation(rotation: Mat3): boolean {
  return rotation.every((v, i) => v === IDENTITY_ROTATION[i]);
}

export default function App() {
  const [showGenerateDialog, setShowGenerateDialog] = useState(false);
  const [seed, setSeed] = useState(randomSeed());
  const [continentalPercent, setContinentalPercent] = useState(DEFAULT_CONTINENTAL_PERCENT);
  const [landPercent, setLandPercent] = useState(DEFAULT_LAND_PERCENT);
  const [axialTiltDeg, setAxialTiltDeg] = useState(DEFAULT_AXIAL_TILT_DEG);

  const [stepYears, setStepYears] = useState(STEP_YEARS_OPTIONS[1]);
  const [projection, setProjection] = useState<Projection>("eckert4");
  const [mapView, setMapView] = useState<MapView>("elevation");
  // The map's current view orientation (see rotation.ts and
  // docs/simulation-model.md#rotating-the-view) -- lives entirely here, like
  // projection/mapView, sent fresh with every render call rather than stored server-side
  // (it's client-local view state, not simulation state). centerLatLon is derived display
  // state: it tracks `rotation` normally, but during an active drag MapCanvas overrides it
  // continuously via onRotationPreview, since the real legend is baked server-side into the
  // PNG and can't update mid-drag.
  const [rotation, setRotation] = useState<Mat3>(IDENTITY_ROTATION);
  const [centerLatLon, setCenterLatLon] = useState({ lat: 0, lon: 0 });
  const [summary, setSummary] = useState<WorldSummary | null>(null);
  const [renderData, setRenderData] = useState<RenderResponse | null>(null);
  // Plate Inspector's own data (see PlateInspector.tsx) -- true-frame/rotation-independent,
  // so unlike renderData it's only refetched on an actual world-state change (generate/step),
  // never on a projection/rotation change. selectedPlateId only makes sense within one
  // world's lifetime (plate ids aren't stable across a regenerate), so it resets there too.
  const [platesData, setPlatesData] = useState<PlateSummary[]>([]);
  const [selectedPlateId, setSelectedPlateId] = useState<number | null>(null);
  // River Inspector's own data (see RiverInspector.tsx) -- same true-frame/world-state-only
  // refresh pattern as platesData. river_id is only meaningful against the most recent
  // /world/rivers response (rivers are regrouped fresh every call, no persistent identity --
  // see backend app/hydrology.py's group_rivers), so it's reset on every generate *and* step,
  // not just generate like selectedPlateId.
  const [riversData, setRiversData] = useState<RiverSummary[]>([]);
  const [selectedRiverId, setSelectedRiverId] = useState<number | null>(null);
  // The land/lake-vs-ocean boundary (see backend app/coastline.py), fetched alongside rivers
  // from the same /world/rivers response -- the River Inspector has no other land/ocean cue.
  const [coastlineSegments, setCoastlineSegments] = useState<Segment[]>([]);
  // Stats panel data (see StatsModal.tsx) -- `stats` is the latest snapshot, `statsHistory`
  // accumulates one entry per generate/step (deduped by elapsed_years) for the panel's graph
  // tabs, matching how plate-sim's own Stats feature builds its history entirely
  // client-side (the backend endpoint itself is stateless, see backend app/stats.py).
  // Recorded continuously, not just while the modal is open, so opening it later still shows
  // the full history since the world was generated.
  const [stats, setStats] = useState<WorldStats | null>(null);
  const [statsHistory, setStatsHistory] = useState<WorldStats[]>([]);
  const [showStatsModal, setShowStatsModal] = useState(false);
  const [playing, setPlaying] = useState(false);
  const [busy, setBusy] = useState(false);
  const [stepping, setStepping] = useState(false);
  const [error, setError] = useState<string | null>(null);

  const refresh = useCallback(async (proj: Projection, view: MapView, viewRotation: Mat3) => {
    if (view === "plateInspector" || view === "riverInspector") return; // neither uses renderData -- see below
    try {
      const data = await renderWorld(proj, view, RENDER_WIDTH, RENDER_HEIGHT, viewRotation);
      setRenderData(data);
    } catch (e) {
      setError(String(e));
    }
  }, []);

  const refreshPlates = useCallback(async () => {
    try {
      const data = await fetchPlates();
      setPlatesData(data.plates);
    } catch (e) {
      setError(String(e));
    }
  }, []);

  const refreshRivers = useCallback(async () => {
    try {
      const data = await fetchRivers();
      setRiversData(data.rivers);
      setCoastlineSegments(data.coastline_segments);
    } catch (e) {
      setError(String(e));
    }
  }, []);

  // Stats are a secondary/best-effort feature -- a failed fetch here (e.g. a transient
  // network blip) shouldn't surface as the main error line or block generate/step, unlike
  // refresh/refreshPlates above which are core to the map actually updating.
  const recordStats = useCallback(async () => {
    try {
      const s = await fetchStats();
      setStats(s);
      setStatsHistory((prev) =>
        prev.length > 0 && prev[prev.length - 1].elapsed_years === s.elapsed_years ? prev : [...prev, s],
      );
    } catch {
      // ignored -- see comment above
    }
  }, []);

  const handleGenerate = useCallback(async () => {
    setBusy(true);
    setError(null);
    try {
      const s = await generateWorld(seed, continentalPercent / 100, landPercent / 100, axialTiltDeg);
      setSummary(s);
      setSelectedPlateId(null);
      setSelectedRiverId(null);
      setShowGenerateDialog(false);
      setStatsHistory([]); // plate ids and elapsed_years both reset with a fresh world
      await Promise.all([refresh(projection, mapView, rotation), refreshPlates(), refreshRivers(), recordStats()]);
    } catch (e) {
      setError(String(e));
    } finally {
      setBusy(false);
    }
  }, [seed, continentalPercent, landPercent, axialTiltDeg, projection, mapView, rotation, refresh, refreshPlates, refreshRivers, recordStats]);

  const handleStep = useCallback(async () => {
    if (!summary) return;
    setStepping(true);
    setError(null);
    try {
      const s = await stepWorld(stepYears);
      setSummary(s);
      setSelectedRiverId(null); // rivers are regrouped fresh every step -- a stale id could point at an unrelated network
      await Promise.all([refresh(projection, mapView, rotation), refreshPlates(), refreshRivers(), recordStats()]);
    } catch (e) {
      setError(String(e));
      setPlaying(false);
    } finally {
      setStepping(false);
    }
  }, [summary, stepYears, projection, mapView, rotation, refresh, refreshPlates, refreshRivers, recordStats]);

  // Resets the view orientation back to the default (lat=0/lon=0, see rotation.ts) -- just
  // updates state, same as a completed drag; the effect below does the actual re-fetch.
  const handleRecenter = useCallback(() => {
    setRotation(IDENTITY_ROTATION);
    setCenterLatLon({ lat: 0, lon: 0 });
  }, []);

  const selectedPlate = platesData.find((p) => p.plate_id === selectedPlateId) ?? null;
  const selectedRiver = riversData.find((r) => r.river_id === selectedRiverId) ?? null;

  // Re-render with the current world whenever the projection, map view, or view rotation
  // changes -- all three are baked server-side into the returned image (see api.ts's
  // renderWorld). A completed drag (see MapCanvas's onRotationCommitted below) just updates
  // `rotation` state; this effect is what actually re-fetches, the same pattern
  // projection/mapView already used before rotation existed.
  useEffect(() => {
    if (summary) {
      refresh(projection, mapView, rotation);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [projection, mapView, rotation]);

  const stepRef = useRef(handleStep);
  stepRef.current = handleStep;

  // Self-scheduling rather than setInterval: each step must finish before the next is
  // scheduled, so a slow step (larger worlds take longer -- a step that also lands on a
  // regularize/gap-fill or reassignment interval, or triggers a merge/split, can take
  // noticeably longer than a routine one) can never overlap with the next tick. A
  // fixed-cadence setInterval would keep firing regardless of whether the previous request
  // had returned.
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
      <div style={{ display: "flex", gap: 16, alignItems: "flex-start" }}>
        <div style={{ display: "flex", flexDirection: "column", gap: 8, width: 170, flexShrink: 0 }}>
          <h1 style={{ fontSize: 20, marginTop: 0, marginBottom: 4, textAlign: "center" }}>Mantle Bloom</h1>
          <p style={{ opacity: 0.7, marginTop: 0, marginBottom: 16 }}>
            Physical World Builder
          </p>

          <button onClick={() => setShowGenerateDialog(true)} disabled={busy} style={{ fontSize: 12 }}>
            Generate World
          </button>

          <button onClick={() => setShowStatsModal(true)} disabled={!summary} style={{ fontSize: 12 }}>
            📊 Stats
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
              <option value="plateInspector">Plate Inspector</option>
              <option value="riverInspector">River Inspector</option>
              <option value="elevation">Elevation</option>
              <option value="temperature">Temperature</option>
              <option value="wind">Wind</option>
              <option value="oceanCurrents">Ocean currents</option>
              <option value="humidity">Humidity</option>
              <option value="precipitation">Precipitation</option>
              <option value="biome">Biome</option>
            </select>
            <select
              value={projection}
              onChange={(e) => setProjection(e.target.value as Projection)}
              style={{ width: "100%", fontSize: 12 }}
            >
              <option value="behrmann">Behrmann (cylindrical equal-area)</option>
              <option value="eckert4">Eckert IV (pseudocylindrical equal-area)</option>
            </select>
            <div style={{ marginTop: 6, opacity: 0.8 }}>
              Center: {formatLatLon(centerLatLon.lat, centerLatLon.lon)}
            </div>
            <button
              onClick={handleRecenter}
              disabled={busy || !summary || isIdentityRotation(rotation)}
              style={{ width: "100%", marginTop: 6, fontSize: 12 }}
            >
              Re-center
            </button>
          </fieldset>

          {mapView === "plateInspector" && (
            <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
              <legend style={{ fontSize: 11 }}>Selected plate</legend>
              {selectedPlate ? (
                <div style={{ opacity: 0.9 }}>
                  <div>id: {selectedPlate.plate_id}</div>
                  <div>crust: {selectedPlate.crust_type}</div>
                  <div>rows: {selectedPlate.num_rows}</div>
                  <div>points: {selectedPlate.num_points}</div>
                  {selectedPlate.bounding_ellipse && (
                    <>
                      <div>diameter A: {selectedPlate.bounding_ellipse.diameter_a_km.toFixed(0)} km</div>
                      <div>diameter B: {selectedPlate.bounding_ellipse.diameter_b_km.toFixed(0)} km</div>
                    </>
                  )}
                </div>
              ) : (
                <div style={{ opacity: 0.6 }}>Click a plate, or press Tab.</div>
              )}
            </fieldset>
          )}

          {mapView === "riverInspector" && (
            <fieldset style={{ border: "1px solid #333", borderRadius: 6, padding: 8, fontSize: 12 }}>
              <legend style={{ fontSize: 11 }}>Selected river</legend>
              {selectedRiver ? (
                <div style={{ opacity: 0.9 }}>
                  <div>id: {selectedRiver.river_id}</div>
                  <div>ends at: {selectedRiver.mouth_type}</div>
                  <div>flow rate: {selectedRiver.flow_rate.toFixed(1)}</div>
                  <div>speed: {selectedRiver.speed.toFixed(2)}</div>
                  <div>tributaries: {selectedRiver.num_tributaries}</div>
                  <div>nodes: {selectedRiver.num_nodes}</div>
                </div>
              ) : (
                <div style={{ opacity: 0.6 }}>
                  {riversData.length > 0 ? "Click a river, or press Tab." : "No rivers yet -- step the world forward."}
                </div>
              )}
            </fieldset>
          )}

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

        <div>
          {mapView === "plateInspector" ? (
            <PlateInspector
              plates={platesData}
              width={RENDER_WIDTH}
              height={RENDER_HEIGHT}
              displayWidth={DISPLAY_WIDTH}
              displayHeight={DISPLAY_HEIGHT}
              projection={projection}
              rotation={rotation}
              selectedPlateId={selectedPlateId}
              onSelectPlate={setSelectedPlateId}
              onRotationPreview={(latDeg, lonDeg) => setCenterLatLon({ lat: latDeg, lon: lonDeg })}
              onRotationCommitted={(newRotation) => setRotation(newRotation)}
            />
          ) : mapView === "riverInspector" ? (
            <RiverInspector
              rivers={riversData}
              coastlineSegments={coastlineSegments}
              width={RENDER_WIDTH}
              height={RENDER_HEIGHT}
              displayWidth={DISPLAY_WIDTH}
              displayHeight={DISPLAY_HEIGHT}
              projection={projection}
              rotation={rotation}
              selectedRiverId={selectedRiverId}
              onSelectRiver={setSelectedRiverId}
              onRotationPreview={(latDeg, lonDeg) => setCenterLatLon({ lat: latDeg, lon: lonDeg })}
              onRotationCommitted={(newRotation) => setRotation(newRotation)}
            />
          ) : (
            <MapCanvas
              imageBase64={renderData?.image_base64 ?? null}
              width={RENDER_WIDTH}
              height={RENDER_HEIGHT}
              displayWidth={DISPLAY_WIDTH}
              displayHeight={DISPLAY_HEIGHT}
              projection={projection}
              rotation={rotation}
              onRotationPreview={(latDeg, lonDeg) => setCenterLatLon({ lat: latDeg, lon: lonDeg })}
              onRotationCommitted={(newRotation) => setRotation(newRotation)}
            />
          )}
          <p style={{ fontSize: 11, opacity: 0.6, marginTop: 4 }}>
            {mapView === "plateInspector"
              ? "Click a plate to select it. Tab / Shift+Tab cycles plates. Press and hold, then drag to rotate."
              : mapView === "riverInspector"
                ? "Click a river to select it. Tab / Shift+Tab cycles rivers. Press and hold, then drag to rotate."
                : "Press and hold, then drag the map to rotate it."}
          </p>
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
              Initial land: {landPercent}%
              <input
                type="range"
                min={0}
                max={100}
                value={landPercent}
                onChange={(e) => setLandPercent(Number(e.target.value))}
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
                onChange={(e) => setContinentalPercent(Number(e.target.value))}
                style={{ width: "100%" }}
              />
            </label>

            <label style={{ display: "block", marginBottom: 16 }}>
              Axial tilt: {axialTiltDeg}°
              <input
                type="range"
                min={0}
                max={45}
                step={0.5}
                value={axialTiltDeg}
                onChange={(e) => setAxialTiltDeg(Number(e.target.value))}
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

      {showStatsModal && <StatsModal stats={stats} history={statsHistory} onClose={() => setShowStatsModal(false)} />}
    </div>
  );
}
