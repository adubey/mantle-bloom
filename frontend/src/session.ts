import { useEffect, useState } from "react";
import { API_BASE } from "./api";
import { asJson } from "./apiResponse";

// GET /auth/session (backend/app/auth.py): who this browser is to the server. `local` means
// this machine, which skips login (#372); `login` is set only for a GitHub session.
interface Session {
  admin: boolean;
  local: boolean;
  login: string | null;
  login_available: boolean;
}

export function useSession(): Session | null {
  const [session, setSession] = useState<Session | null>(null);
  useEffect(() => {
    let cancelled = false;
    fetch(`${API_BASE}/auth/session`, { credentials: "include" })
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
