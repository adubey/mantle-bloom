// Minimal document.cookie helpers -- just enough to persist the map's view state (see
// App.tsx's VIEW_COOKIE_NAME) across a browser refresh. Not a general cookie library: no
// attribute options beyond the fixed ones below, since this app only ever needs the one.

const COOKIE_MAX_AGE_SECONDS = 60 * 60 * 24 * 365; // 1 year -- long-lived, like any "remember my last view" preference

export function setCookie(name: string, value: string): void {
  document.cookie = `${encodeURIComponent(name)}=${encodeURIComponent(value)};max-age=${COOKIE_MAX_AGE_SECONDS};path=/;samesite=lax`;
}

export function getCookie(name: string): string | null {
  const prefix = `${encodeURIComponent(name)}=`;
  for (const part of document.cookie.split(";")) {
    const trimmed = part.trim();
    if (trimmed.startsWith(prefix)) return decodeURIComponent(trimmed.slice(prefix.length));
  }
  return null;
}
