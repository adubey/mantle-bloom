import { useState } from "react";

// A look at the full app for /view visitors (#560): static screenshots of views the public
// routes don't serve, in a <details> that starts open so visitors see them immediately
// (and can collapse them). Images under public/preview/ set loading="lazy" and
// fetchPriority="low" as scheduling hints to keep map rendering ahead of them, and are
// static files rather than rendered from the shared world, so the preview costs no
// simulation or render work and adds no public API surface.

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

interface BillingActionProps {
  apiBase: string;
  manage: boolean;
  label?: string;
  inline?: boolean;
}

export function BillingActionButton({ apiBase, manage, label, inline = false }: BillingActionProps) {
  const [busy, setBusy] = useState(false);
  const [error, setError] = useState<string | null>(null);

  async function handleClick() {
    setBusy(true);
    setError(null);
    try {
      const resp = await fetch(`${apiBase}${manage ? "/billing/portal" : "/billing/checkout"}`, {
        method: "POST",
      });
      if (!resp.ok) {
        const body = await resp.json().catch(() => null) as { detail?: unknown } | null;
        const message = typeof body?.detail === "string" ? body.detail : `Request failed (${resp.status})`;
        throw new Error(message);
      }
      const data = await resp.json() as { url: string };
      window.location.assign(data.url);
    } catch (err) {
      setError(err instanceof Error ? err.message : "Billing unavailable");
      setBusy(false);
    }
  }

  return (
    <span style={inline ? { display: "inline", marginLeft: 8 } : { display: "inline-flex", alignItems: "center", gap: 6 }}>
      <button type="button" onClick={handleClick} disabled={busy}>
        {label ?? (manage ? "Manage subscription" : "Upgrade")}
      </button>
      {error && (
        <span role="alert" style={{ color: "#ff9b9b", fontSize: 12, marginLeft: inline ? 8 : 0 }}>
          {error}
        </span>
      )}
    </span>
  );
}

interface FullVersionPreviewProps {
  repoUrl: string;
  apiBase?: string;
  billingAvailable?: boolean;
  signedIn?: boolean;
  manageSubscription?: boolean;
}

export default function FullVersionPreview({
  repoUrl,
  apiBase = "",
  billingAvailable = false,
  signedIn = false,
  manageSubscription = false,
}: FullVersionPreviewProps) {
  return (
    <details className="full-preview" open>
      <summary>Features of the full version</summary>
      <div className="full-preview-lead">
        This page shows one shared world in a single map view. The full app generates your own worlds from a seed or a
        sketch, runs them for as long as you like, and looks at them on a globe and in over twenty map views, with
        plates, faults, rivers, lakes and resources. You can run it today from the <a href={repoUrl}>source on GitHub</a>
        {billingAvailable ? (
          signedIn ? (
            manageSubscription ? (
              <>
                .<BillingActionButton apiBase={apiBase} manage inline />
              </>
            ) : (
              <>
                , or run it from mantle-bloom.com with a{" "}
                <BillingActionButton apiBase={apiBase} manage={false} label="Paid Upgrade" inline />.
              </>
            )
          ) : (
            <>
              . <a href={`${apiBase}/login`}>Sign in to upgrade</a>.
            </>
          )
        ) : (
          ". Accounts with your own hosted worlds are planned."
        )}
      </div>
      <ul className="full-preview-grid">
        {PREVIEWS.map((p) => (
          <li key={p.src}>
            <figure>
              <a href={p.src} target="_blank" rel="noreferrer" aria-label={`${p.title}: open the full-size image`}>
                <img
                  src={p.src}
                  width={p.width}
                  height={p.height}
                  alt={p.alt}
                  loading="lazy"
                  fetchPriority="low"
                  decoding="async"
                />
              </a>
              <figcaption>
                <strong>{p.title}</strong>
                <span>{p.caption}</span>
              </figcaption>
            </figure>
          </li>
        ))}
      </ul>
    </details>
  );
}
