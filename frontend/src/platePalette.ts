// A direct, hand-synced port of backend app/render_image.py's PLATE_PALETTE -- same 20 RGB
// tuples, in the same order, so a plate's color in the Plate Inspector visually matches its
// color in the "Plates" view. Keep in sync by hand if the backend palette ever changes (same
// convention already established for the old frontend/src/elevationColor.ts before rendering
// moved server-side).
export const PLATE_PALETTE: [number, number, number][] = [
  [230, 25, 75], [60, 180, 75], [255, 225, 25], [67, 99, 216], [245, 130, 49],
  [66, 212, 244], [240, 50, 230], [188, 246, 12], [250, 190, 190], [70, 153, 144],
  [230, 190, 255], [154, 99, 36], [255, 250, 200], [128, 0, 0], [170, 255, 195],
  [128, 128, 0], [255, 216, 177], [0, 0, 117], [169, 169, 169], [255, 255, 255],
];

export function plateColor(plateId: number): [number, number, number] {
  const n = PLATE_PALETTE.length;
  return PLATE_PALETTE[((plateId % n) + n) % n];
}
