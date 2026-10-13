import { API_BASE } from "./api";
import { useSession } from "./session";

// The logged-out landing page at /view. For now it only points the owner at the login; the
// public map viewer itself is #373.
export default function Viewer() {
  const session = useSession();
  return (
    <main
      style={{
        fontFamily: "system-ui, sans-serif",
        color: "#e6e8ef",
        maxWidth: 420,
        margin: "0 auto",
        padding: "64px 16px",
        textAlign: "center",
      }}
    >
      <h1 style={{ fontSize: 24, marginBottom: 4 }}>Mantle Bloom</h1>
      <p style={{ opacity: 0.7, marginTop: 0 }}>Physical World Builder</p>
      <p>The public viewer isn't open yet.</p>
      {session?.admin ? (
        <a href="/">Open the full app</a>
      ) : (
        session?.login_available && <a href={`${API_BASE}/login`}>Sign in</a>
      )}
    </main>
  );
}
