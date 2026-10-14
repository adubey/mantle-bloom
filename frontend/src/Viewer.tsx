import { useCallback, useEffect, useRef, useState } from "react";
import { API_BASE } from "./api";
import Legend from "./Legend";
import MapCanvas from "./MapCanvas";
import {
  fetchPublicState,
  PublicStepRefused,
  renderPublicWorld,
  stepPublicWorld,
  type PublicState,
  type PublicView,
} from "./publicApi";
import { IDENTITY_ROTATION } from "./rotation";
import { useSession } from "./session";

// The logged-out viewer at /view (#373): the shared world in Eckert IV, biome or elevation
// only, plus a fixed-size step limited per browser. Everything else is the full app at /,
// and the server enforces that split (backend app/public_view.py), not this page.

const REPO_URL = "https://github.com/adubey/mantle-bloom";
const MAX_DISPLAY_WIDTH = 1100;
// The Eckert IV map's aspect, same as the full app's 1100x611 frame.
const ASPECT = 611 / 1100;
const RENDER_SCALE = 2;
const PAGE_GUTTER = 16;
const VIEWS: { id: PublicView; label: string }[] = [
  { id: "elevation", label: "Elevation" },
  { id: "biome", label: "Biome" },
];
const noop = () => {};
const POLL_RETRY_MS = 5000;

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

function AllowanceNote({ state }: { state: PublicState }) {
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
      <a href={REPO_URL}>Download Mantle Bloom from GitHub</a> to run your own worlds.
    </span>
  );
}

function AccountLink() {
  const session = useSession();
  if (session?.admin) return <a href="/">Open the full app</a>;
  if (session?.login_available) return <a href={`${API_BASE}/login`}>Sign in</a>;
  return null;
}

export default function Viewer() {
  const [state, setState] = useState<PublicState | null>(null);
  const [view, setView] = useState<PublicView>("elevation");
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
    fetchPublicState().then(setState, (err) => setError(String(err)));
  }, []);

  // Follow the shared world: a long poll that answers when anyone (another viewer, or the
  // owner) changes it, so this page re-renders without stepping itself (#374). Restarted on
  // every revision change, including this viewer's own steps.
  useEffect(() => {
    if (!loaded) return;
    const seen = revision ?? "";
    const controller = new AbortController();
    let timer: number | undefined;
    const poll = () => {
      fetchPublicState(seen, controller.signal).then(
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
  }, [loaded, revision]);

  useEffect(() => {
    if (revision === null) return;
    let cancelled = false;
    renderPublicWorld(view, MAX_DISPLAY_WIDTH * RENDER_SCALE, Math.round(MAX_DISPLAY_WIDTH * ASPECT) * RENDER_SCALE)
      .then((r) => {
        if (!cancelled) setImage(r.image_base64);
      })
      .catch((err) => {
        if (!cancelled) setError(String(err));
      });
    return () => {
      cancelled = true;
    };
  }, [view, revision]);

  const step = useCallback(async () => {
    if (revision === null) return;
    setStepping(true);
    setNotice(null);
    setError(null);
    try {
      const result = await stepPublicWorld(revision);
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
  }, [revision]);

  const exhausted = state?.steps_remaining === 0;
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
        {state?.world && <span style={{ opacity: 0.7 }}>{formatMyr(state.world.elapsed_years)} elapsed</span>}
        <span style={{ flex: 1 }} />
        <a href={REPO_URL}>Source &amp; download</a>
        <AccountLink />
      </header>
      {state && !state.world && <p>No world has been generated yet. Check back later.</p>}
      {state?.world && (
        <>
          <div style={{ display: "flex", flexWrap: "wrap", alignItems: "center", gap: 12, marginBottom: 12 }}>
            <div role="radiogroup" aria-label="Map view" style={{ display: "flex", gap: 4 }}>
              {VIEWS.map((v) => (
                <button key={v.id} role="radio" aria-checked={view === v.id} onClick={() => setView(v.id)}
                  style={{ fontWeight: view === v.id ? 700 : 400 }}>
                  {v.label}
                </button>
              ))}
            </div>
            <button onClick={step} disabled={stepping || exhausted}>
              {stepping ? "Stepping..." : `Step ${state.step_years / 1000} kyr`}
            </button>
            <AllowanceNote state={state} />
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
              alphaEncodedIds={view === "biome"}
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
        </>
      )}
      {error && <p role="alert" style={{ color: "#ff9b9b" }}>{error}</p>}
    </main>
  );
}
