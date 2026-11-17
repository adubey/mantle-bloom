import { resetDom } from "./domSetup";
import assert from "node:assert/strict";
import { afterEach, beforeEach, describe, it } from "node:test";
import { cleanup, render, screen } from "@testing-library/react";
import AdminSession from "../src/AdminSession";
import { canAccessApp, type Session } from "../src/session";

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

  it("renders Manage subscription and Sign out for a paid user", () => {
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
    assert.ok(manage.closest("form")?.getAttribute("action")?.endsWith("/billing/portal"));
    const signOut = screen.getByRole("button", { name: "Sign out" });
    assert.ok(signOut.closest("form")?.getAttribute("action")?.endsWith("/logout"));
  });

  it("renders nothing for local admin mode and login + Sign out for GitHub admin", () => {
    const { rerender } = render(<AdminSession session={makeSession({ admin: true, local: true })} />);
    assert.equal(screen.queryByRole("button"), null);

    rerender(<AdminSession session={makeSession({ admin: true, login: "adubey" })} />);
    assert.ok(screen.getByText(/adubey/));
    assert.ok(screen.getByRole("button", { name: "Sign out" }));
  });
});
