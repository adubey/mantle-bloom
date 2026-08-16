type Stop = [elevation: number, rgb: [number, number, number]];

// A simple hypsometric tint: deep blue -> shallow blue -> coastal green -> tan -> brown ->
// white peaks. Matches the terrain-style colormap used while developing/verifying the
// simulation (see the throwaway matplotlib previews in docs/simulation-model.md).
const STOPS: Stop[] = [
  [-11000, [10, 10, 40]],
  [-4000, [15, 40, 110]],
  [-1500, [40, 110, 190]],
  [-200, [110, 170, 210]],
  [0, [200, 210, 150]],
  [200, [90, 150, 60]],
  [1200, [170, 160, 90]],
  [3000, [120, 90, 60]],
  [6000, [230, 230, 230]],
  [9000, [255, 255, 255]],
];

export function elevationColor(elevation: number): string {
  const clamped = Math.min(Math.max(elevation, STOPS[0][0]), STOPS[STOPS.length - 1][0]);

  let lower = STOPS[0];
  let upper = STOPS[STOPS.length - 1];
  for (let i = 0; i < STOPS.length - 1; i++) {
    if (clamped >= STOPS[i][0] && clamped <= STOPS[i + 1][0]) {
      lower = STOPS[i];
      upper = STOPS[i + 1];
      break;
    }
  }

  const [e0, [r0, g0, b0]] = lower;
  const [e1, [r1, g1, b1]] = upper;
  const t = e1 === e0 ? 0 : (clamped - e0) / (e1 - e0);
  const r = Math.round(r0 + (r1 - r0) * t);
  const g = Math.round(g0 + (g1 - g0) * t);
  const b = Math.round(b0 + (b1 - b0) * t);
  return `rgb(${r}, ${g}, ${b})`;
}
