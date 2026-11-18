import { resetDom } from "./domSetup";
import assert from "node:assert/strict";
import { afterEach, beforeEach, describe, it } from "node:test";
import { cleanup, fireEvent, render, screen, waitFor } from "@testing-library/react";
import AdminSession from "../src/AdminSession";
import { canAccessApp, isDebugUiEnabled, type Session } from "../src/session";

function makeSession(overrides: Partial<Session> = {}): Session {
  return {
    admin: false,
    local: false,
    login: null,
    login_available: true,
    billing_available: true,
    user: null,
    ...overrides,
  };
}

describe("AdminSession", () => {
  beforeEach(() => {
    resetDom("http://localhost/");
  });

  afterEach(() => {
    cleanup();
  });

  it("allows admin and paid sessions into the app, and excludes free/logged-out sessions", () => {
    assert.equal(canAccessApp(makeSession({ admin: true })), true);
    assert.equal(canAccessApp(makeSession({ user: { name: "Ada", provider: "google", tier: "paid" } })), true);
    assert.equal(canAccessApp(makeSession({ user: { name: "Ada", provider: "google", tier: "free" } })), false);
    assert.equal(canAccessApp(makeSession()), false);
  });

  it("gates ?deb debug UI strictly on admin sessions", () => {
    const admin = makeSession({ admin: true });
    const paid = makeSession({ user: { name: "Ada", provider: "google", tier: "paid" } });
    assert.equal(isDebugUiEnabled("?deb", admin), true);
    assert.equal(isDebugUiEnabled("", admin), false);
    assert.equal(isDebugUiEnabled("?deb", paid), false);
    assert.equal(isDebugUiEnabled("?deb", null), false);
  });

  it("posts via fetch for Manage subscription and renders a logout form for a paid user", async () => {
    const realFetch = globalThis.fetch;
    const calls: string[] = [];
    globalThis.fetch = (async (input: string | URL | Request): Promise<Response> => {
      const url = typeof input === "string" ? input : input instanceof URL ? input.toString() : input.url;
      calls.push(url);
      return new Response(JSON.stringify({ url: "https://polar.sh/portal/test" }), {
        status: 200,
        headers: { "Content-Type": "application/json" },
      });
    }) as typeof fetch;
    try {
      render(
        <AdminSession
          session={makeSession({
            billing_available: true,
            user: { name: "Ada", provider: "google", tier: "paid", manage_subscription: true },
          })}
        />,
      );
      assert.ok(screen.getByText(/Ada/));
      const manage = screen.getByRole("button", { name: "Manage subscription" });
      assert.equal(manage.getAttribute("type"), "button");
      fireEvent.click(manage);
      await waitFor(() => {
        assert.ok(calls.some((u) => u.endsWith("/billing/portal")));
      });
      const signOut = screen.getByRole("button", { name: "Sign out" });
      assert.ok(signOut.closest("form")?.getAttribute("action")?.endsWith("/logout"));
    } finally {
      globalThis.fetch = realFetch;
    }
  });

  it("renders nothing for local admin mode and login + Sign out for GitHub admin", () => {
    const { rerender } = render(<AdminSession session={makeSession({ admin: true, local: true })} />);
    assert.equal(screen.queryByRole("button"), null);

    rerender(<AdminSession session={makeSession({ admin: true, login: "adubey" })} />);
    assert.ok(screen.getByText(/adubey/));
    assert.ok(screen.getByRole("button", { name: "Sign out" }));
  });
});
