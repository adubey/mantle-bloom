export function nearbyCrossSectionPositions(centerIndex: number, totalPositions: number, limit: number): number[] {
  if (!Number.isInteger(totalPositions) || totalPositions < 1) {
    throw new RangeError("totalPositions must be a positive integer");
  }
  const center = ((centerIndex % totalPositions) + totalPositions) % totalPositions;
  const neighborCount = Math.max(0, Math.min(Math.floor(limit), totalPositions - 1));
  const positions: number[] = [];

  for (let distance = 1; positions.length < neighborCount; distance++) {
    positions.push((center - distance + totalPositions) % totalPositions);
    if (positions.length < neighborCount) positions.push((center + distance) % totalPositions);
  }

  return positions;
}
