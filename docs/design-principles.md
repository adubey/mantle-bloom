# Simulation design principles

These principles guide changes to the planet model. They describe the behavior the
simulation is trying to preserve, including where an approximation is acceptable. When a
change exposes a conflict between principles, document the tradeoff and measure its effect.

## 1. State the question and the model boundary

Before choosing an equation or adding state, state what question the simulation change should
answer. Identify the physical mechanisms needed to answer it, the conditions and scales where
the model applies, and the important effects it leaves out. This makes scope testable and
keeps a local improvement from quietly claiming to model a larger system.

Represent the core mechanisms that control the outcome at the simulation's spatial and time
scales. Simplify secondary or much smaller-scale effects when that preserves the large-scale
behavior the model is meant to capture and keeps updates affordable. Record the approximation
and its limits; add detail when evidence shows the simplified mechanism changes the result of
interest.

## 2. Keep the planet moving at a useful pace

Users should see meaningful change after a reasonably quick update. Prefer bounded,
vectorized or compiled work over per-node searches and repeated passes over the same graph.
Share expensive calculations when their inputs have not changed. Measure runtime when a
change adds a solve or pass, and keep the cost in view alongside the geological result.

This is a product constraint on model design, not a reason to skip important feedback. Issue
[#297](https://github.com/adubey/mantle-bloom/issues/297) asks for a timestep-integrated
landscape solve while avoiding a path search for every sediment source; issue
[#355](https://github.com/adubey/mantle-bloom/issues/355) identifies duplicate routing as
both a performance cost and a source of inconsistent water behavior.

## 3. Model fast processes in aggregate over long steps

Tectonic steps span thousands of years to millions of years. A step cannot resolve a daily or
seasonal cycle, but those cycles can still affect its integrated result. Represent them with
an explicit average or closed-form approximation where they matter. For example, annual mean
temperature below freezing does not imply that a glacier never melts: seasonal warm periods
can produce melt even when the mean is below 0 °C.

Do not let a long step skip intermediate physical thresholds or routes. Where practical,
compare one long step with several shorter steps and state the tolerance. Record which
processes use annual means, imputed seasons, or step-integrated rates so coupled systems use
the same interpretation. See [#275](https://github.com/adubey/mantle-bloom/issues/275) and
[#297](https://github.com/adubey/mantle-bloom/issues/297).

## 4. Treat each terrain node as an area, not a point

Surface nodes represent finite control volumes. Weight water, rock, ice, sediment, and area
totals by the represented cell area. At coarse resolution, the cell center does not capture
every valley, pass, shoreline, or partially flooded margin; use shared-edge passes or compact
subgrid geometry when the process depends on those features. Avoid rules whose behavior
changes just because the same catchment has more nodes.

Check important outcomes at more than one node density when resolution sensitivity is
plausible. Keep algorithms compatible with the current sparse surface rather than solving a
fine global mesh by default. This addresses the coarse-cell and node-density failures described
in [#297](https://github.com/adubey/mantle-bloom/issues/297) and
[#355](https://github.com/adubey/mantle-bloom/issues/355).

## 5. Give quantities explicit, consistent units

Use physical units where the model can support them: metres, square metres, cubic metres,
years, and rates such as m³/s or m³/yr. Name units in state and API fields where practical,
and document units for every rate, threshold, and conversion.

If a process still uses a calibrated or dimensionless model quantity, label it as such and
keep it within that process. Do not add quantities with different dimensions or silently
reinterpret a value at a subsystem boundary. Put conversions in one named place and verify
that changing timestep or node density does not change a rate that should be independent of
either. Caps and empirical coefficients need units and a reason for their scale. Issues
[#355](https://github.com/adubey/mantle-bloom/issues/355),
[#363](https://github.com/adubey/mantle-bloom/issues/363), and
[#368](https://github.com/adubey/mantle-bloom/issues/368) show how one unit mismatch can
propagate into river flow, erosion, lake silt, and glacier response.

## 6. Account for material through every process

Erosion and redistribution move material; they do not make it disappear. Track volume and
provenance through erosion, transport, deposition, remapping, collisions, rifting, and
topology changes. Every change to a tracked budget must end in a named reservoir, transfer,
or declared sink, with an audit that can expose clipping and residuals. Keep continental
material identifiable when it sits on oceanic crust or moves between plates.

Test budget closure directly. A cap may protect numerical stability, but material clipped by
the cap still needs an explicit destination and diagnostic. See
[#275](https://github.com/adubey/mantle-bloom/issues/275),
[#340](https://github.com/adubey/mantle-bloom/issues/340), and
[#378](https://github.com/adubey/mantle-bloom/issues/378).

## 7. Preserve cratonic crust; make its loss path explicit

Cratons are long-lived, strong crustal blocks. They should resist internal ductile flow and
help determine where deformation is concentrated; ordinary shortening should preferentially
crumple weaker surrounding crust. Large-scale loss of cratonic crust should occur through
explicit tectonic attrition at margins, not as a side effect of generic relaxation, remapping,
or bookkeeping. Cratonic mantle keels can thin or delaminate without implying that the
overlying cratonic crust was destroyed.

“Resistant” does not mean exempt from every process: underthrusting, margin erosion, and
other explicitly modeled recycling paths can act on cratonic material. Separate crust from
keel and state which reservoir a process changes. This principle follows the measured loss
paths and geology in [#340](https://github.com/adubey/mantle-bloom/issues/340) and the
craton-survival checks in [#275](https://github.com/adubey/mantle-bloom/issues/275).

## 8. Let each process follow its own physics

Processes may share terrain, channels, routing infrastructure, and ledgers, but they should
not inherit one another's driving law without a physical reason. Rivers respond to water
discharge and basin filling; glacier ice moves under gravity from ice thickness and ice-surface
slope; landslides fail under shaking and slope, then run out under gravity. Keep failure,
transport, and deposition distinct when they have different controls. Couple processes at
clear interfaces and pass quantities with explicit units and provenance.

Shared infrastructure should make each process cheaper and more consistent, not make them
physically interchangeable. See the process-specific routing requirements in
[#351](https://github.com/adubey/mantle-bloom/issues/351) and the real-unit coupling questions
in [#355](https://github.com/adubey/mantle-bloom/issues/355).

## 9. Verify the implementation, validate the model, and show uncertainty

When a budget or geological result is wrong, split diagnostics by the code path that books or
changes it before choosing a fix. Report starting state, interval, seed, units, and the share
from each path. A catch-all ledger may close while hiding an unintended physical loss.

**Verification** checks that code implements its intended equations and invariants. Use
equation-level checks and known or limiting cases as evidence. **Validation** checks whether
the chosen equations and assumptions reproduce observations or accepted planetary
benchmarks. Passing software checks does not by itself show that a model is physically
adequate.

Quantify uncertainty where results depend on uncertain parameters, random seeds, resolution,
or boundary assumptions. Compare distributions and trajectories across representative runs,
not only one seed or a pooled mean. Record runtime with the geological measures. Do not tune a
parameter to compensate for unmeasured numerical loss or a unit mismatch. Issues
[#340](https://github.com/adubey/mantle-bloom/issues/340) and
[#378](https://github.com/adubey/mantle-bloom/issues/378) illustrate why source attribution
should precede a broad intervention.

## 10. Design for a planetary history, not one favorable interval

Evaluate behavior over the timescales the simulation represents: short steps for local
stability, hundreds of millions of years for tectonic and climate evolution, and multi-billion
year runs for long-term planetary persistence. As a working horizon, allow for Earth to remain
broadly present-like for only about 3–5 billion years; make the assumptions behind that span
explicit rather than treating today's state as an unlimited equilibrium.

Use multiple seeds and inspect trajectories, not only endpoints. Track whether land,
cratons, crustal thickness, water, ice, and plate topology remain plausible over time, and
distinguish intended planetary evolution from numerical drift.

## Issue and change descriptions

For a simulation bug or proposed model change, name the governing principle and the current
behavior that violates it. State the observable invariant or metric that will show the issue
is fixed. For a deliberate exception, explain the approximation and its limits. This lets the
next change address the cause instead of repeating the same failure in another subsystem.
