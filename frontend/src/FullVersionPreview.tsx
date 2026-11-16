// A look at the full app for /view visitors (#560): static screenshots of views the public
// routes don't serve. Images are files under public/preview/, not rendered from the shared
// world, so the preview costs no simulation or render work and adds no public API surface.

const PREVIEWS = [
  {
    src: "/preview/temperature.webp",
    width: 1200,
    height: 619,
    title: "Climate maps",
    caption: "Temperature, wind, humidity and precipitation, computed from the planet's land, sea and ice.",
    alt: "A world map coloured by surface temperature, warm across the tropics and cold over mountain belts and the poles.",
  },
  {
    src: "/preview/ocean-currents.webp",
    width: 1200,
    height: 637,
    title: "Ocean currents",
    caption: "Warm and cold surface currents steered around the continents.",
    alt: "A world map with arrows tracing ocean surface currents around dark continents.",
  },
  {
    src: "/preview/cross-section.webp",
    width: 1200,
    height: 681,
    title: "Lithosphere cross-section",
    caption: "Cut along any great circle to see crust, mantle lithosphere and mantle flow beneath the surface.",
    alt: "A curved slice through the planet's outer 350 km showing crust, mantle lithosphere and asthenosphere, beside a globe marking where the cut runs.",
  },
] as const;

export default function FullVersionPreview({ repoUrl }: { repoUrl: string }) {
  return (
    <section aria-labelledby="full-version-heading" className="full-preview">
      <h2 id="full-version-heading">More in the full app</h2>
      <p className="full-preview-lead">
        This page shows one shared world in a single map view. The full app generates your own worlds from a seed or a
        sketch, runs them for as long as you like, and looks at them on a globe and in over twenty map views, with
        plates, faults, rivers, lakes and resources. You can run it today from the <a href={repoUrl}>source on GitHub</a>.
        Accounts with your own hosted worlds are planned.
      </p>
      <ul className="full-preview-grid">
        {PREVIEWS.map((p) => (
          <li key={p.src}>
            <figure>
              <a href={p.src} target="_blank" rel="noreferrer" aria-label={`${p.title}: open the full-size image`}>
                <img src={p.src} width={p.width} height={p.height} alt={p.alt} loading="lazy" decoding="async" />
              </a>
              <figcaption>
                <strong>{p.title}</strong>
                <span>{p.caption}</span>
              </figcaption>
            </figure>
          </li>
        ))}
      </ul>
    </section>
  );
}
