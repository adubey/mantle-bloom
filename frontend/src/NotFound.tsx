export default function NotFound() {
  return (
    <main
      style={{
        minHeight: "100vh",
        boxSizing: "border-box",
        display: "grid",
        placeContent: "center",
        justifyItems: "center",
        gap: 12,
        padding: 24,
        color: "#e8ebf2",
        textAlign: "center",
      }}
    >
      <p style={{ margin: 0, color: "#61ded0", fontSize: 13, letterSpacing: ".12em" }}>404</p>
      <h1 style={{ margin: 0, fontSize: 28 }}>Page not found</h1>
      <p style={{ margin: 0, color: "#aab0c0" }}>That address doesn’t lead to a page here.</p>
      <a href="/" style={{ color: "#61ded0", marginTop: 4 }}>Go home</a>
    </main>
  );
}
