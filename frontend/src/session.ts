import { useEffect, useState } from "react";
import { API_BASE } from "./api";
import { asJson } from "./apiResponse";

// GET /auth/session (backend/app/auth.py): who this browser is to the server. `local` means
// this machine, which skips login (#372); `login` is set only for the owner's GitHub session.
// `user` is a signed-in user (#929, #945): identity and subscription state for UI controls.
export interface Session {
  admin: boolean;
  local: boolean;
  login: string | null;
  login_available: boolean;
  billing_available?: boolean;
  user: {
    name: string | null;
    provider: string;
    tier?: "free" | "paid";
    manage_subscription?: boolean;
  } | null;
}

export function canAccessApp(session: Session): boolean {
  return session.admin || session.user?.tier === "paid";
}

export function useSession(pollForPaid = false): Session | null {
  const [session, setSession] = useState<Session | null>(null);
  useEffect(() => {
    let cancelled = false;
    let timer: ReturnType<typeof setTimeout> | undefined;
    let attempts = 0;

    function load() {
      // No `credentials: "include"`: production is same-origin, so the cookie goes anyway, and the
      // dev frontend is local (no cookie needed) and its CORS setup doesn't allow credentials.
      fetch(`${API_BASE}/auth/session`)
        .then((resp) => asJson<Session>(resp))
        .then((s) => {
          if (cancelled) return;
          setSession(s);
          if (pollForPaid && s.user && s.user.tier !== "paid" && attempts < 6) {
            attempts += 1;
            timer = setTimeout(load, 1500);
          }
        })
        .catch(() => {
          // Leave it null: the page still works, minus the sign-in/out link.
        });
    }

    load();
    return () => {
      cancelled = true;
      if (timer !== undefined) clearTimeout(timer);
    };
  }, [pollForPaid]);
  return session;
}
