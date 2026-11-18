import { useEffect } from "react";
import { API_BASE } from "./api";
import { BillingActionButton } from "./FullVersionPreview";
import { canAccessApp, type Session } from "./session";

// The full app's account control. The server already sends free/logged-out visits to "/" over to
// /view; this catches a tab whose session expired, downgraded, or was signed out elsewhere.
export default function AdminSession({ session }: { session: Session | null }) {
  useEffect(() => {
    if (session && !canAccessApp(session)) window.location.replace("/view");
  }, [session]);
  if (!session) return null;
  if (session.admin) {
    if (!session.login) return null;
    return (
      <form method="post" action={`${API_BASE}/logout`} style={{ margin: 0, fontSize: 11, textAlign: "center" }}>
        <span style={{ opacity: 0.7 }}>{session.login} · </span>
        <button type="submit" style={{ fontSize: 11 }}>
          Sign out
        </button>
      </form>
    );
  }
  if (session.user?.tier !== "paid") return null;
  return (
    <div style={{ display: "flex", flexWrap: "wrap", alignItems: "center", justifyContent: "center", gap: 4, fontSize: 11 }}>
      <span style={{ opacity: 0.7 }}>{session.user.name ?? "Signed in"} ·</span>
      {session.billing_available && <BillingActionButton apiBase={API_BASE} manage />}
      <form method="post" action={`${API_BASE}/logout`} style={{ margin: 0 }}>
        <button type="submit" style={{ fontSize: 11 }}>
          Sign out
        </button>
      </form>
    </div>
  );
}
