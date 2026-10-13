const LOCATOR_EQUATOR_FLATTENING = 0.36;
const LOCATOR_CAMERA_TILT = Math.asin(LOCATOR_EQUATOR_FLATTENING);

export interface ProjectedGreatCirclePoint {
  x: number;
  y: number;
  depth: number;
}

/** Project backend great-circle coordinates onto a familiar tilted globe view. */
export function projectGreatCirclePoint(
  planeAngleDeg: number,
  arcAngle: number,
  cx: number,
  cy: number,
  radius: number,
): ProjectedGreatCirclePoint {
  const planeAngle = (planeAngleDeg * Math.PI) / 180;
  const planeCosine = Math.cos(planeAngle);
  const planeSine = Math.sin(planeAngle);
  const arcCosine = Math.cos(arcAngle);
  const arcSine = Math.sin(arcAngle);
  const equatorialAxis = planeSine * arcCosine;
  const alongAxis = arcSine;
  const screenVertical = equatorialAxis * Math.cos(LOCATOR_CAMERA_TILT)
    - alongAxis * Math.sin(LOCATOR_CAMERA_TILT);
  const depth = alongAxis * Math.cos(LOCATOR_CAMERA_TILT)
    + equatorialAxis * Math.sin(LOCATOR_CAMERA_TILT);

  return {
    x: cx + radius * planeCosine * arcCosine,
    y: cy - radius * screenVertical,
    depth,
  };
}
