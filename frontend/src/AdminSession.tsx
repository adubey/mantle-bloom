import { useEffect } from "react";
import { API_BASE } from "./api";
import { useSession } from "./session";

// The full app's account control. The server already sends logged-out visits to "/" over to
// /view; this catches a tab whose session expired or was signed out elsewhere.
export default function AdminSession() {
  const session = useSession();
  useEffect(() => {
    if (session && !session.admin) window.location.replace("/view");
  }, [session]);
  if (!session?.login) return null;
  return (
    <form method="post" action={`${API_BASE}/logout`} style={{ margin: 0, fontSize: 11, textAlign: "center" }}>
      <span style={{ opacity: 0.7 }}>{session.login} · </span>
      <button type="submit" style={{ fontSize: 11 }}>
        Sign out
      </button>
    </form>
  );
}
