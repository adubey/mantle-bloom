# Reviewer handoff: RFC 0010 Gate 3b evolving slab ribbons

**Requested review:** Independently review the issue #699 working changes and the preregistered
evolving-slab mini-simulation. Report correctness or interpretation findings by severity, identify
which RFC 0010 Gate 3b requirements remain unmet, and recommend the next diagnostic. This is a
review handoff; do not treat the experiment as a Gate 3b acceptance result.

## Workstream state

- RFC 0010, “Rifts follow weak lithosphere,” is **Accepted**. Gate 3a local crack mechanics is
  closed. Gate 3b and its tracking issue #699 are **open**.
- The current RFC revision requires an evolving 3-D slab body coupled to deforming plate geometry
  before slab pull can count as a rift source. It then requires a broad sustained extension event,
  local crack linkage into a viable rift, a 10–30 Myr duration bracket, timestep sensitivity, and
  plate-scale breakup histories before Gate 3b is complete. See the [RFC gate criteria](../../rfc/0010-rifts-follow-weak-lithosphere.md).
- The latest probe is explicitly a reduced-order **implementation/mechanics screen**. It did not
  produce a viable split and does not close Gate 3b. The parent issue remains open.
- Review the complete current diff in branch `codex/issue-699-plate-scale-rift`, including its nine
  commits ahead of `origin/main` and the uncommitted working changes. The probe, output, image and
  report have not all been committed. Preserve the dirty worktree while reviewing.

## Hypothesis and frozen experiment

The motivating hypothesis is that magma/mantle flow and slab pull load broad parts of a plate;
opposed slab pull may stretch its interior, and Gate 3a's local crack/link rules should accumulate
damage into a rift before the topology splits. The owner requested material-bearing 3-D slabs and
endogenous plate geometry because rigid-plate trench forces cannot establish this mechanism.

The experiment was preregistered in [issue comment 6540115898](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6540115898),
with numerical assumptions frozen in [comment 6540154180](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6540154180):

- Input: saved world seed 447276134, plate 14; embedded age 116.4 Myr. The input SHA-256 is
  `1ccae7159f4a6326831ba032078a408f1a692ca4b59ad315e798d7103c71c513`. It was selected before
  observing this response because its two active fronts pass the bilateral lighter-neighbor screen
  (75 and 41 source cells; baseline resultants 152.2° apart, weaker/stronger ratio 0.598).
- Transfer front lithosphere once into two 3-D ribbons; run matched 20 Myr cases at 2 Myr steps.
  The native-strength geometric midline mantle-current source is common to all arms. Arm A applies
  no slab force to the plate, B applies the stronger flank only, C applies both flanks. Run each at
  p10, median and p90 Gate 2 material-capacity cases. Keep the existing local damage, failed-link,
  and daughter-size rules fixed; a viable split needs exactly two daughters, each at least 10% of
  initial plate area.
- Frozen body/shell assumptions: 45° initial dip, 660 km depth parameter, 8 stations, density
  anomaly 80 kg/m³, bending modulus 70 GPa, 273 K surface, 3e-4 K/m geothermal gradient, mantle
  viscosity 1e21 Pa·s, 20 km shear zone, shell extensional viscosity 1e24 Pa·s and 15° load spread.
  These are prototype values, not calibrated Earth parameters.
- Registered validity limits: material relative residual ≤1e-10; shell force, torque, equilibrium
  and elastic-work residuals each ≤1e-3; slab bending/relaxation work residual ≤5%. The first
  loaded-arm attempt failed force-moment accounting and is explicitly excluded in
  [issue comment 6540204654](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6540204654).
  The correction preserved the same 15° footprint and source magnitudes while constraining net
  force and moment, and separated kinematic convergence advection from beam-relaxation work.

## Results to verify

All nine corrected runs completed and passed the probe's implemented checks; none produced a
viable split. Independent review found that those checks do not yet establish all preregistered
coupled-mechanics accounts; see the review findings below.
The p10 results were:

| Slab-force arm | Failed cells / links | Components | Two largest area fractions | Max shell strain |
|---|---:|---:|---:|---:|
| No slab force | 0 / 0 | 2 | 99.9580% / 0.0420% | 1.20e-4 |
| Stronger flank only | 13 / 8 | 2 | 99.9583% / 0.0417% | 9.76e-4 |
| Both flanks | 1 / 0 | 2 | 99.9589% / 0.0411% | 5.30e-3 |

Median and p90 produced no failed cells or links in any arm. The second component is a single-cell
island present before damage, so the reported component count is not evidence of breakup. The
both-flank case had the greatest strain but no connected crack path. Plate 14 has zero craton cells;
other plates in the input world do contain cratons.

At 20 Myr, reported maximum slab length was 1,599–1,607 km and maximum depth was 1,176–1,182 km.
The largest material-volume residual was 1.61e-16. Maximum shell force, torque, equilibrium and
elastic-work residuals were 2.23e-14, 5.68e-14, 3.32e-4 and 3.49e-9; maximum slab relaxation-work
residual was 3.56e-13. The screen reports a deep-reservoir volume while conserving total ribbon
volume.

**Specific preregistration question:** the 660 km parameter is exceeded by the reported maximum
body depth. The implementation treats it as the threshold for moving station material into the
deep-reservoir account, while continuing to advect and report the ribbon centerline geometry below
that depth. Review whether that interpretation satisfies the frozen “maximum depth” setting and
whether the result should instead be marked as a preregistration deviation or reported with separate
active-body and reservoir depths. Do not silently describe the ribbon as capped at 660 km.

The corrected results were posted in [issue comment 6540311299](https://github.com/adubey/mantle-bloom/issues/699#issuecomment-6540311299).
The report, compact output and montage are in `analysis/issue699/README.md`,
`evolving-slabs-seed447276134-14.json`, and `evolving-slabs-seed447276134-14.png`. The output
SHA-256 is `4939cd0d7f677198980cce0a01c3643423d96fb0cc1f950f900983712b3e822c`; the montage SHA-256
is `8e19c9553270d9efc520d88ae0f15b7e5eebdefcefe49b01e786adbe6ba910e1`.

## Files and review focus

- `backend/app/slab_body.py`: slab state, geometry, advancement, force transfer and shell/plate
  deformation coupling.
- `backend/app/sparse_quad_patch.py` and `backend/app/thin_shell.py`: current worktree changes to
  geometry/solver support.
- `backend/unit_tests/test_slab_body.py`, `test_sparse_quad_patch.py`, `test_thin_shell.py`: added or
  modified tests. They were **not run** for this handoff.
- `bin/debug/probe_gate3b_evolving_slabs.py`: preregistered experiment implementation and montage.
- `bin/debug/probe_gate3b_component_shell.py`: prior shell diagnostic changes used by this work.
- `analysis/issue699/README.md`, JSON output and image: evidence and visual results.
- `docs/architecture.md` and `docs/simulation-model/tectonics.md`: model documentation changes.

Please specifically assess:

1. Whether material transfer, station volumes, active/deep-reservoir ownership and conservation are
   internally consistent; whether old saves and state lifecycle are affected.
2. Whether slab bending, imposed convergence advection, mantle drag, buoyancy, attachment action /
   reaction, and shell load distribution have defensible units and force, moment and work accounts.
3. Whether the measured convergence speed and evolving ribbon load are coupled to plate geometry
   without double-counting existing trench slab pull or midline mantle flow.
4. Whether the applied shell strain is genuinely spatial/endogenous enough to support the stated
   local cracking inference, and whether the damage/capacity pairing and component accounting are
   correct. Check the pre-existing singleton baseline and area-weighted daughter rule.
5. Whether reported validity checks are independent and meaningful. In particular, inspect the
   surface inventory check and the 660 km depth/deep-reservoir semantics rather than relying only
   on small residuals.
6. Whether the experiment followed its preregistration and whether the conclusion is properly
   limited. The current run has no detached-slab control, no timestep partition comparison, only
   one selected plate and one 20 Myr duration, and no viable rift. Those remain explicit gaps against
   RFC 0010 Gate 3b; this negative result cannot satisfy the event/source or breakup criteria.

## First independent review (2026-10-31)

A read-only review characterized the run as a limited shell-load-transfer smoke screen, not a
validated coupled slab-mechanics result. The reviewer prioritized these findings:

1. **P1 — 660 km depth is not an enforced body-depth cap.** The active mask disables buoyancy for
   stations beyond that depth, but advection moves stations deeper and the beam/drag solve still
   includes them. The 1,176–1,182 km reported depths therefore expose a deviation from the frozen
   parameter. Resolve whether to cap/transfer the geometry, report active-body and reservoir depth
   separately, or classify the current result as exploratory and preregister a corrected run.
2. **P1 — surface inventory residual is tautological; tracer material is not transferred.** The
   probe calculates before/after surface inventories consecutively without an intervening mutation,
   so the residual cannot detect loss during the strain update. Source cells also have continental
   and craton tracers zeroed without those extensive tracer quantities entering the ribbons or a
   separate ledger. Close each tracer and compute the pre/post inventory around the actual geometry
   update.
3. **P1 — slab/plate hinge equilibrium is not closed.** The shell receives a prescribed traction
   derived from active excess-density material, while the beam's attachment reaction is not solved
   and applied back to slab dynamics. The reported bending/drag resultants do not provide hinge
   force, moment or work closure. The shell's small transfer residuals validate redistribution of
   that prescribed traction only.
4. **P1/P2 — ribbon dimensions and slice inventories need review.** Initial source volume spread
   over the centerline yields 3.1–3.5 km thickness; length later grows to about 1,600 km while
   thickness, rigidity and station volumes stay fixed. Eight slices receive equal mass despite
   endpoint volumes being about half the interior volumes. Either make the body cross-section and
   per-slice mass material-consistent as geometry changes, or label these values as geometric
   diagnostics only.
5. **P2 — plate deformation is spatially affine.** The shell computes local displacement/stress,
   but only a fitted homogeneous strain tensor updates the plate; local displacement and slab
   attachment coordinates do not deform with that field. This can screen plate-wide strain and
   heterogeneous demand, but cannot yet establish localized deformation around the slab hinges.
6. **P2 — event evidence and controls are incomplete.** The output lacks a broad-extension footprint,
   direction/capacity field and sustained-event qualification. It has one 20 Myr duration at one
   timestep on one plate, without detached-slab or current-free/tectonic-only controls. The
   pre-existing singleton also means split qualification should compare against baseline topology
   while preserving the registered two-daughter rule.

Recommended next probe: a minimal corrected beam–plate coupling with independent hinge force,
moment and work residuals; explicit 660 km/deep-reservoir semantics; material/tracer closure; and
attached, detached and no-body controls. Only after those checks should the duration/timestep screen
and later production-resolution breakup histories resume. The reviewer found no viable rift in this
run and agreed it does not close Gate 3b.

## Reproduction

Run the single memory-intensive probe alone (do not run a replay or test suite beside it), with the
registered save available locally:

```bash
PYTHONHASHSEED=0 PYTHONPATH=backend:bin/debug NUMBA_CACHE_DIR=<private-cache> \
  backend/.venv/bin/python bin/debug/probe_gate3b_evolving_slabs.py \
  --save <registered-seed-447276134.mbworld> \
  --out analysis/issue699/evolving-slabs-seed447276134-14.json \
  --image analysis/issue699/evolving-slabs-seed447276134-14.png
```

The save remains a local input and is not included in the repository. The private file path should
not be copied into review notes or commits.
