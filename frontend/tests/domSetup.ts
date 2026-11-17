import { Window } from "happy-dom";

// Install a single happy-dom Window on globalThis before React or Testing
// Library are imported. Call resetDom(url) in beforeEach to isolate URL,
// history, and DOM state between tests.
const domWindow = new Window({ url: "http://localhost/view" });

const DOM_GLOBALS = [
  "window",
  "document",
  "navigator",
  "location",
  "history",
  "HTMLElement",
  "HTMLInputElement",
  "HTMLSelectElement",
  "HTMLButtonElement",
  "HTMLCanvasElement",
  "HTMLImageElement",
  "Image",
  "SVGElement",
  "Element",
  "Node",
  "Text",
  "Event",
  "CustomEvent",
  "MouseEvent",
  "KeyboardEvent",
  "MutationObserver",
  "getComputedStyle",
  "matchMedia",
  "localStorage",
  "sessionStorage",
  "requestAnimationFrame",
  "cancelAnimationFrame",
] as const;

for (const key of DOM_GLOBALS) {
  Object.defineProperty(globalThis, key, {
    value: (domWindow as unknown as Record<string, unknown>)[key],
    configurable: true,
    writable: true,
  });
}

// Forward window timer methods to globalThis so node:test's t.mock.timers
// intercepts window.setInterval / window.clearInterval calls in components.
Object.defineProperties(domWindow, {
  setInterval: { get: () => globalThis.setInterval, configurable: true },
  clearInterval: { get: () => globalThis.clearInterval, configurable: true },
});

// React checks IS_REACT_ACT_ENVIRONMENT to suppress act() warnings when
// Testing Library drives updates.
(globalThis as Record<string, unknown>).IS_REACT_ACT_ENVIRONMENT = true;

export function resetDom(url = "http://localhost/view"): void {
  domWindow.happyDOM.setURL(url);
  domWindow.document.head.innerHTML = "";
  domWindow.document.body.innerHTML = "";
  domWindow.localStorage.clear();
  domWindow.sessionStorage.clear();
  for (const cookie of domWindow.document.cookie.split(";")) {
    const name = cookie.split("=")[0]?.trim();
    if (name) domWindow.document.cookie = `${name}=; expires=Thu, 01 Jan 1970 00:00:00 GMT; path=/`;
  }
}
