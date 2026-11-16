import { resetDom } from "./domSetup";
import assert from "node:assert/strict";
import { afterEach, beforeEach, describe, it, type TestContext } from "node:test";
import { cleanup, fireEvent, render, screen, waitFor } from "@testing-library/react";
import Viewer from "../src/Viewer";
import type { PublicState, PublicWorldListing } from "../src/publicApi";

const TINY_PNG_BASE64 =
  "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mNkYAAAAAYAAjCB0C8AAAAASUVORK5CYII=";

function makeListing(overrides: Partial<PublicWorldListing> & { world_id: string; name: string }): PublicWorldListing {
  return {
    revision: "rev-1",
    seed: 42,
    elapsed_years: 100_000,
    running: false,
    is_default: false,
    ...overrides,
  };
}

function makeState(listing: PublicWorldListing, overrides: Partial<PublicState> = {}): PublicState {
  return {
    name: listing.name,
    view: "combined",
    rotation: [1, 0, 0, 0, 1, 0, 0, 0, 1],
    is_default: listing.is_default,
    world: {
      revision: listing.revision,
      seed: listing.seed,
      elapsed_years: listing.elapsed_years,
    },
    running: listing.running,
    steps_remaining: 3,
    step_limit: 3,
    step_years: 50_000,
    step_progress: null,
    ...overrides,
  };
}

function jsonResponse(body: unknown, status = 200): Response {
  return new Response(JSON.stringify(body), {
    status,
    headers: { "Content-Type": "application/json" },
  });
}

function abortablePending(signal?: AbortSignal | null): Promise<Response> {
  return new Promise((_resolve, reject) => {
    if (signal?.aborted) {
      reject(new DOMException("Aborted", "AbortError"));
      return;
    }
    signal?.addEventListener("abort", () => reject(new DOMException("Aborted", "AbortError")), { once: true });
  });
}

interface FetchMockOptions {
  worlds: () => PublicWorldListing[];
  states?: Record<string, Partial<PublicState>>;
  onStep?: (worldId: string) => Response | Promise<Response>;
  admin?: boolean;
  user?: { name: string | null; provider: string } | null;
}

function installFetchMock({ worlds, states = {}, onStep, admin = false, user = null }: FetchMockOptions): { worldsCalls: () => number } {
  let publicWorldsCalls = 0;
  globalThis.fetch = (async (input: string | URL | Request, init?: RequestInit): Promise<Response> => {
    const rawUrl = typeof input === "string" ? input : input instanceof URL ? input.toString() : input.url;
    const url = new URL(rawUrl, "http://localhost");

    if (url.pathname === "/auth/session") {
      return jsonResponse({ admin, local: false, login: admin ? "owner" : null, login_available: false, user });
    }
    if (url.pathname === "/public/worlds") {
      publicWorldsCalls += 1;
      return jsonResponse({ worlds: worlds() });
    }

    const match = /^\/public\/worlds\/([^/]+)\/(world|render|step)$/.exec(url.pathname);
    if (match) {
      const worldId = decodeURIComponent(match[1]);
      const endpoint = match[2];
      const listing = worlds().find((w) => w.world_id === worldId);
      if (endpoint === "world") {
        // Long-poll requests carry ?after=<revision>; keep them pending until aborted.
        if (url.searchParams.has("after")) return abortablePending(init?.signal);
        if (!listing) return jsonResponse({ detail: "Not found" }, 404);
        return jsonResponse(makeState(listing, states[worldId]));
      }
      if (endpoint === "render") {
        return jsonResponse({
          revision: listing?.revision ?? "rev-1",
          seed: listing?.seed ?? 42,
          elapsed_years: listing?.elapsed_years ?? 100_000,
          image_base64: TINY_PNG_BASE64,
        });
      }
      if (endpoint === "step" && onStep) {
        return onStep(worldId);
      }
    }

    throw new Error(`Unexpected fetch: ${url.pathname}${url.search}`);
  }) as typeof fetch;

  return { worldsCalls: () => publicWorldsCalls };
}

describe("Viewer (/view)", () => {
  const realFetch = globalThis.fetch;

  beforeEach(() => {
    resetDom("http://localhost/view");
  });

  afterEach(() => {
    cleanup();
    globalThis.fetch = realFetch;
  });

  it("adopts a lone live world when no world is selected (#597)", async () => {
    const idle = makeListing({ world_id: "w-idle", name: "Kepler-186 f", running: false });
    const live = makeListing({ world_id: "w-live", name: "TRAPPIST-1 e", running: true, elapsed_years: 2_400_000 });
    installFetchMock({ worlds: () => [idle, live] });

    render(<Viewer />);

    await screen.findByText("TRAPPIST-1 e");
    assert.equal(window.location.search, "?world=w-live");
    const select = screen.getByLabelText("World") as HTMLSelectElement;
    assert.equal(select.value, "w-live");
    assert.ok(screen.getByText("● Live"));
  });

  it("prefers the default world over a lone live world", async () => {
    const live = makeListing({ world_id: "w-live", name: "TRAPPIST-1 e", running: true });
    const def = makeListing({ world_id: "w-def", name: "Proxima Centauri b", is_default: true });
    installFetchMock({ worlds: () => [live, def] });

    render(<Viewer />);

    await screen.findByText("Proxima Centauri b");
    assert.equal(window.location.search, "?world=w-def");
  });

  it("falls back to the first listed world when multiple live worlds exist without a default", async () => {
    const idle = makeListing({ world_id: "w-idle", name: "Kepler-186 f", running: false });
    const live1 = makeListing({ world_id: "w-live-1", name: "TRAPPIST-1 e", running: true });
    const live2 = makeListing({ world_id: "w-live-2", name: "LHS 1140 b", running: true });
    installFetchMock({ worlds: () => [idle, live1, live2] });

    render(<Viewer />);

    await screen.findByText("Kepler-186 f");
    assert.equal(window.location.search, "?world=w-idle");
  });

  it("keeps an explicit valid ?world= even when another world is marked default", async () => {
    resetDom("http://localhost/view?world=w-explicit");
    const def = makeListing({ world_id: "w-def", name: "Proxima Centauri b", is_default: true });
    const explicit = makeListing({ world_id: "w-explicit", name: "Kepler-186 f" });
    installFetchMock({ worlds: () => [def, explicit] });

    render(<Viewer />);

    await screen.findByText("Kepler-186 f");
    assert.equal(window.location.search, "?world=w-explicit");
  });

  it("keeps the selected world on 15s polling when a new live world appears (#597)", async (t: TestContext) => {
    resetDom("http://localhost/view?world=w1");
    t.mock.timers.enable({ apis: ["setInterval"] });

    const w1 = makeListing({ world_id: "w1", name: "Kepler-186 f", running: false });
    const w2 = makeListing({ world_id: "w2", name: "TRAPPIST-1 e", running: true });
    let list = [w1];
    const { worldsCalls } = installFetchMock({ worlds: () => list });

    render(<Viewer />);

    await screen.findByText("Kepler-186 f");
    assert.equal(worldsCalls(), 1);

    list = [w1, w2];
    t.mock.timers.tick(15_000);

    await waitFor(() => assert.equal(worldsCalls(), 2));
    const select = (await screen.findByLabelText("World")) as HTMLSelectElement;
    assert.equal(select.value, "w1");
    assert.equal(window.location.search, "?world=w1");
    assert.ok(screen.getByRole("option", { name: /TRAPPIST-1 e ● live/ }));
  });

  it("shows empty message when no worlds exist and falls back from unknown ?world=", async () => {
    let list: PublicWorldListing[] = [];
    installFetchMock({ worlds: () => list });

    const { unmount } = render(<Viewer />);
    await screen.findByText("No world has been generated yet. Check back later.");
    unmount();

    resetDom("http://localhost/view?world=missing");
    list = [makeListing({ world_id: "w1", name: "Kepler-186 f" })];
    render(<Viewer />);
    await screen.findByText("Kepler-186 f");
    assert.equal(window.location.search, "?world=w1");
  });

  it("advances a world, handles caught-up notices, and surfaces 429/503 refusals", async () => {
    const w1 = makeListing({ world_id: "w1", name: "Kepler-186 f", elapsed_years: 100_000 });
    let stepMode: "caught-up" | "429" | "503" = "caught-up";
    installFetchMock({
      worlds: () => [w1],
      onStep: () => {
        if (stepMode === "429") return jsonResponse({ detail: "limit" }, 429);
        if (stepMode === "503") return jsonResponse({ detail: "busy" }, 503);
        return jsonResponse({
          world: { revision: "rev-2", seed: 42, elapsed_years: 150_000 },
          computed: false,
          steps_remaining: 2,
        });
      },
    });

    render(<Viewer />);

    const advanceBtn = (await screen.findByRole("button", { name: "Advance 50,000 years" })) as HTMLButtonElement;
    fireEvent.click(advanceBtn);

    await screen.findByText("The world had already changed, so you're caught up to it for free.");
    assert.ok(screen.getByText(/2 of 3 advances left in this browser/));

    stepMode = "503";
    fireEvent.click(advanceBtn);
    await screen.findByText("The world is busy with another run right now. Try again in a moment.");

    stepMode = "429";
    fireEvent.click(advanceBtn);
    await waitFor(() => assert.equal(advanceBtn.disabled, true));
    assert.ok(screen.getByText(/You've used all 3 of this browser's advances\./));
  });

  it("offers guests a closed expando of the full version's other views (#560)", async () => {
    installFetchMock({ worlds: () => [makeListing({ world_id: "w1", name: "Kepler-186 f" })] });

    render(<Viewer />);

    await screen.findByText("Kepler-186 f");
    const details = screen.getByText("Features of the full version").closest("details");
    assert.ok(details);
    // Closed until asked for: the full version isn't live yet.
    assert.equal(details.open, false);
    const images = [...details.querySelectorAll("img")];
    assert.deepEqual(
      images.map((img) => img.getAttribute("src")),
      ["/preview/temperature.webp", "/preview/ocean-currents.webp", "/preview/cross-section.webp"],
    );
    // Lazy, so a closed expando downloads nothing.
    assert.ok(images.every((img) => img.getAttribute("alt") && img.getAttribute("loading") === "lazy"));
  });

  it("offers the expando even when no world is public yet (#560)", async () => {
    installFetchMock({ worlds: () => [] });

    render(<Viewer />);

    await screen.findByText("No world has been generated yet. Check back later.");
    assert.ok(screen.getByText("Features of the full version"));
  });

  it("shows a signed-in free user their name and a sign-out form instead of Sign in (#929)", async () => {
    installFetchMock({
      worlds: () => [makeListing({ world_id: "w1", name: "Kepler-186 f" })],
      user: { name: "Ada Lovelace", provider: "google" },
    });

    render(<Viewer />);

    const signOut = await screen.findByRole("button", { name: "Sign out" });
    assert.ok(screen.getByText(/Ada Lovelace/));
    assert.equal(screen.queryByRole("link", { name: "Sign in" }), null);
    const form = signOut.closest("form");
    assert.equal(form?.getAttribute("method"), "post");
    assert.ok(form?.getAttribute("action")?.endsWith("/logout"));
    // Signed in is still the free viewer: the browser's step allowance is unchanged.
    assert.ok(screen.getByText(/3 of 3 advances left in this browser/));
  });

  it("offers Sign in to a signed-out visitor (#929)", async () => {
    installFetchMock({ worlds: () => [makeListing({ world_id: "w1", name: "Kepler-186 f" })] });

    render(<Viewer />);

    const link = await screen.findByRole("link", { name: "Sign in" });
    assert.ok(link.getAttribute("href")?.endsWith("/login"));
    assert.equal(screen.queryByRole("button", { name: "Sign out" }), null);
  });

  it("hides the preview from the admin, who already has the full app (#560)", async () => {
    installFetchMock({ worlds: () => [makeListing({ world_id: "w1", name: "Kepler-186 f" })], admin: true });

    render(<Viewer />);

    await screen.findByText("Kepler-186 f");
    await waitFor(() => assert.equal(screen.queryByText("Features of the full version"), null));
  });
});
