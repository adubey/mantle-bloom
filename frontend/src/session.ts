import { useEffect, useState } from "react";
import { API_BASE } from "./api";
import { asJson } from "./apiResponse";

// GET /auth/session (backend/app/auth.py): who this browser is to the server. `local` means
// this machine, which skips login (#372); `login` is set only for the owner's GitHub session.
// `user` is a signed-in free account (#929): it identifies the visitor and unlocks nothing.
interface Session {
  admin: boolean;
  local: boolean;
  login: string | null;
  login_available: boolean;
  user: { name: string | null; provider: string } | null;
}

export function useSession(): Session | null {
  const [session, setSession] = useState<Session | null>(null);
  useEffect(() => {
    let cancelled = false;
    // No `credentials: "include"`: production is same-origin, so the cookie goes anyway, and the
    // dev frontend is local (no cookie needed) and its CORS setup doesn't allow credentials.
    fetch(`${API_BASE}/auth/session`)
      .then((resp) => asJson<Session>(resp))
      .then((s) => {
        if (!cancelled) setSession(s);
      })
      .catch(() => {
        // Leave it null: the page still works, minus the sign-in/out link.
      });
    return () => {
      cancelled = true;
    };
  }, []);
  return session;
}
