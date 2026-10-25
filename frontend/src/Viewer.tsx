import { useCallback, useEffect, useRef, useState } from "react";
import { API_BASE } from "./api";
import Legend from "./Legend";
import MapCanvas from "./MapCanvas";
import {
  fetchPublicState,
  fetchPublicWorlds,
  PublicStepRefused,
  renderPublicWorld,
  stepPublicWorld,
  type PublicState,
  type PublicView,
  type PublicWorldListing,
} from "./publicApi";
import { withWorldId, worldIdFromSearch, type WorldId } from "./worldRoutes";
import { IDENTITY_ROTATION } from "./rotation";
import { useSession } from "./session";

// The logged-out viewer at /view (#373): the shared world in Eckert IV, elevation and biome
// only, plus a fixed-size step limited per browser. /demo/view previews this UI locally;
// everything else is the full app at /, and the server enforces the API split.

const REPO_URL = "https://github.com/adubey/mantle-bloom";
const COFFEE_URL = "https://buymeacoffee.com/mantlebloom";
const MAX_DISPLAY_WIDTH = 1100;
// The Eckert IV map's aspect, same as the full app's 1100x611 frame.
const ASPECT = 611 / 1100;
const RENDER_SCALE = 2;
const PAGE_GUTTER = 16;
const noop = () => {};
const POLL_RETRY_MS = 5000;
const DEMO_STEP_PAUSE_MS = 1000;
const isDemo = window.location.pathname.replace(/\/+$/, "") === "/demo/view";

function GitHubIcon() {
  return (
    <svg aria-hidden="true" viewBox="0 0 16 16" width="16" height="16" fill="currentColor" style={{ verticalAlign: "-3px", marginRight: 5 }}>
      <path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38l-.01-1.49c-2.23.49-2.7-.95-2.7-.95-.36-.92-.89-1.16-.89-1.16-.73-.5.06-.49.06-.49.8.06 1.22.82 1.22.82.72 1.21 1.88.86 2.34.66.07-.52.28-.86.51-1.06-1.78-.2-3.64-.89-3.64-3.96 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82a7.65 7.65 0 0 1 4 0c1.53-1.03 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.28.82 2.15 0 3.08-1.87 3.76-3.65 3.96.29.25.54.73.54 1.48l-.01 2.19c0 .21.15.46.55.38A8.01 8.01 0 0 0 16 8c0-4.42-3.58-8-8-8Z" />
    </svg>
  );
}

function useDisplayWidth(): number {
  const measure = () => Math.min(MAX_DISPLAY_WIDTH, window.innerWidth - 2 * PAGE_GUTTER);
  const [width, setWidth] = useState(measure);
  useEffect(() => {
    // The render is re-requested only on a view or world change, not per resize: the canvas
    // just scales the existing frame.
    const onResize = () => setWidth(measure());
    window.addEventListener("resize", onResize);
    return () => window.removeEventListener("resize", onResize);
  }, []);
  return width;
}

function formatMyr(years: number): string {
  return `${(years / 1e6).toFixed(1)} Myr`;
}

function refusalMessage(err: unknown): string {
  if (err instanceof PublicStepRefused && err.status === 503) {
    return "The world is busy with another run right now. Try again in a moment.";
  }
  return `Step failed: ${err instanceof Error ? err.message : String(err)}`;
}

function AllowanceNote({ state, demo = false }: { state: PublicState; demo?: boolean }) {
  if (demo) {
    return state.steps_remaining === 0 ? <span role="status">Out of quota for stepping world.</span> : null;
  }
  if (state.steps_remaining === null) return <span style={{ opacity: 0.7 }}>Unlimited steps (signed in)</span>;
  if (state.steps_remaining > 0) {
    return (
      <span style={{ opacity: 0.7 }}>
        {state.steps_remaining} of {state.step_limit} steps left in this browser
      </span>
    );
  }
  return (
    <span role="status">
      You've used all {state.step_limit} of this browser's simulation steps.{" "}
      <a href={REPO_URL}><GitHubIcon />GitHub</a> to run your own worlds.
    </span>
  );
}

function AccountLink() {
  const session = useSession();
  if (isDemo || session?.admin) return null;
  if (session?.login_available) return <a href={`${API_BASE}/login`}>Sign in</a>;
  return null;
}

// Which public world to show (#563): `?world=<id>`, the only public one when
// there is one, or a picker when there are several.
function usePublicWorldId(): { worlds: PublicWorldListing[] | null; worldId: WorldId | null; selectWorld: (id: WorldId) => void } {
  const [worlds, setWorlds] = useState<PublicWorldListing[] | null>(null);
  const [worldId, setWorldId] = useState<WorldId | null>(() => worldIdFromSearch(window.location.search));

  useEffect(() => {
    fetchPublicWorlds().then(
      ({ worlds: list }) => {
        setWorlds(list);
        const urlId = worldIdFromSearch(window.location.search);
        if (urlId && list.some((w) => w.world_id === urlId)) return;
        if (list.length === 0) return;
        // No selection and exactly one public world: adopt it. A selection naming an
        // unknown, private, or deleted world falls back to the first listed one
        // instead of leaving a blank page (#563 S8).
        if ((!urlId && list.length === 1) || urlId) {
          const id = list[0].world_id;
          setWorldId(id);
          window.history.replaceState(null, "", withWorldId(window.location.href, id));
        }
      },
      () => {},
    );
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  const selectWorld = useCallback((id: WorldId) => {
    setWorldId(id);
    window.history.replaceState(null, "", withWorldId(window.location.href, id));
  }, []);

  return { worlds, worldId, selectWorld };
}

export default function Viewer() {
  const { worlds, worldId, selectWorld } = usePublicWorldId();

  return (
    <main
      style={{
        fontFamily: "system-ui, sans-serif",
        color: "#e6e8ef",
        maxWidth: MAX_DISPLAY_WIDTH,
        margin: "0 auto",
        padding: `24px ${PAGE_GUTTER}px`,
      }}
    >
      <header style={{ display: "flex", flexWrap: "wrap", alignItems: "baseline", gap: 12, marginBottom: 12 }}>
        <h1 style={{ fontSize: 22, margin: 0 }}>Mantle Bloom</h1>
        <span style={{ flex: 1 }} />
        <a href={REPO_URL}><GitHubIcon />GitHub</a>
        <a href={COFFEE_URL} target="_blank" rel="noreferrer">☕ Buy me a coffee</a>
        <AccountLink />
      </header>
      {worlds && worlds.length > 1 && (
        <label style={{ display: "block", marginBottom: 12 }}>
          World{" "}
          <select value={worldId ?? ""} onChange={(e) => selectWorld(e.target.value)} aria-label="World">
            {worlds.map((w) => (
              <option key={w.world_id} value={w.world_id}>
                {w.name} · {formatMyr(w.elapsed_years)}
              </option>
            ))}
          </select>
        </label>
      )}
      {!worldId && worlds && worlds.length === 0 && <p>No world has been generated yet. Check back later.</p>}
      {worldId && <WorldView key={worldId} worldId={worldId} />}
    </main>
  );
}

// One public world's page: its state, stepping, and map (#563).
function WorldView({ worldId }: { worldId: WorldId }) {
  const [state, setState] = useState<PublicState | null>(null);
  const view: PublicView = "combined";
  const [image, setImage] = useState<string | null>(null);
  const [stepping, setStepping] = useState(false);
  const [notice, setNotice] = useState<string | null>(null);
  const [error, setError] = useState<string | null>(null);
  const displayWidth = useDisplayWidth();
  const displayHeight = Math.round(displayWidth * ASPECT);
  const mapWrapperRef = useRef<HTMLDivElement>(null);
  const revision = state?.world?.revision ?? null;
  const loaded = state !== null;

  useEffect(() => {
    let cancelled = false;
    fetchPublicState(worldId, undefined, undefined, isDemo).then(
      (s) => {
        if (!cancelled) setState(s);
      },
      (err) => {
        if (!cancelled) setError(String(err));
      },
    );
    return () => {
      cancelled = true;
    };
  }, [worldId]);

  // Follow the selected public world: a long poll that answers when anyone
  // (another viewer, or the owner) changes it, so this page re-renders
  // without stepping itself (#374). Restarted on every revision change,
  // including this viewer's own steps.
  useEffect(() => {
    if (!loaded) return;
    const seen = revision ?? "";
    const controller = new AbortController();
    let timer: number | undefined;
    const poll = () => {
      fetchPublicState(worldId, seen, controller.signal, isDemo).then(
        (next) => {
          if ((next.world?.revision ?? "") === seen) poll();
          else setState(next);
        },
        () => {
          // Back off on errors (server restarting, network blip) rather than spinning.
          if (!controller.signal.aborted) timer = window.setTimeout(poll, POLL_RETRY_MS);
        },
      );
    };
    poll();
    return () => {
      controller.abort();
      window.clearTimeout(timer);
    };
  }, [loaded, revision, worldId]);

  useEffect(() => {
    if (revision === null) return;
    let cancelled = false;
    renderPublicWorld(worldId, view, MAX_DISPLAY_WIDTH * RENDER_SCALE, Math.round(MAX_DISPLAY_WIDTH * ASPECT) * RENDER_SCALE)
      .then((r) => {
        if (!cancelled) setImage(r.image_base64);
      })
      .catch((err) => {
        if (!cancelled) setError(String(err));
      });
    return () => {
      cancelled = true;
    };
  }, [view, revision, worldId]);

  const step = useCallback(async () => {
    if (revision === null) return;
    setStepping(true);
    setNotice(null);
    setError(null);
    try {
      const result = await stepPublicWorld(worldId, revision, isDemo);
      setState((s) => s && { ...s, world: result.world, steps_remaining: result.steps_remaining });
      if (!result.computed) setNotice("The world had already changed, so you're caught up to it for free.");
    } catch (err) {
      if (err instanceof PublicStepRefused && err.status === 429) {
        setState((s) => s && { ...s, steps_remaining: 0 });
      } else {
        setError(refusalMessage(err));
      }
    } finally {
      setStepping(false);
    }
  }, [revision, worldId]);

  const exhausted = state?.steps_remaining === 0;

  useEffect(() => {
    if (!isDemo || revision === null || stepping || exhausted) return;
    const timer = window.setTimeout(() => void step(), DEMO_STEP_PAUSE_MS);
    return () => window.clearTimeout(timer);
  }, [exhausted, revision, stepping, step]);

  if (state && !state.world) return <p>This world is no longer public. Check back later.</p>;
  // A failed fetch used to fall through to a blank page, swallowing the message (#563 S8).
  if (!state?.world) {
    if (error) return <p role="alert" style={{ color: "#ff9b9b" }}>{error}</p>;
    return null;
  }
  const worldState = state.world;
  return (
    <>
      <div style={{ display: "flex", flexWrap: "wrap", alignItems: "center", gap: 12, marginBottom: 12 }}>
        <span style={{ opacity: 0.7 }}>{formatMyr(worldState.elapsed_years)} elapsed</span>
        <span style={{ flex: 1 }} />
        {!isDemo && (
          <>
            <button onClick={step} disabled={stepping || exhausted}>
              {stepping ? "Stepping..." : `Step ${state.step_years / 1000} kyr`}
            </button>
            <AllowanceNote state={state} />
          </>
        )}
        {isDemo && <AllowanceNote state={state} demo />}
      </div>
      {notice && <p style={{ opacity: 0.8 }}>{notice}</p>}
      <div ref={mapWrapperRef}>
        <MapCanvas
          imageBase64={image}
          width={MAX_DISPLAY_WIDTH * RENDER_SCALE}
          height={Math.round(MAX_DISPLAY_WIDTH * ASPECT) * RENDER_SCALE}
          displayWidth={displayWidth}
          displayHeight={displayHeight}
          projection="eckert4"
          rotation={IDENTITY_ROTATION}
          onRotationPreview={noop}
          onRotationCommitted={noop}
          alphaEncodedIds
          interactionDisabled
        />
      </div>
      <Legend
        mapView={view}
        projection="eckert4"
        width={MAX_DISPLAY_WIDTH * RENDER_SCALE}
        height={Math.round(MAX_DISPLAY_WIDTH * ASPECT) * RENDER_SCALE}
        displayWidth={displayWidth}
        displayHeight={displayHeight}
        mapWrapperRef={mapWrapperRef}
        hasWorld
      />
      {error && <p role="alert" style={{ color: "#ff9b9b" }}>{error}</p>}
    </>
  );
}
