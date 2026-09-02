"""`LithospherePlate`: the plate representation. Subclasses `PlateWithLines` (see
plates.py) rather than reinventing its plumbing -- everything not about *how a plate moves or
deforms* (outline tracing, containment, neighbour search, node iteration, the whole render/
erosion/hydrology/stats/persistence surface) is inherited unchanged, since every one of those
consumers only ever reads `Plate`'s abstract interface (`all_points_and_elevation`, `collect`,
`contains_batch`, ...), never `shift`/`deform`'s own internals. See the plan's "Why subclass
PlateWithLines" section.

Overridden here: `shift` (torque integration, torque.py), `deform` (Mohr-Coulomb/isostasy,
rheology.py + lithosphere.py), `merge_with`/`_merge_nodes_with`/`split`/`grow_into` (carrying
Hc/Hm through the representation's own resample/partition operations, which
`PlateWithLines`'s versions don't know about).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import cKDTree

from . import geometry
from .elevation_lines import (
    ELEV_CHANGE_COLLISION,
    ELEV_CHANGE_MIN_DELTA_M,
    ELEV_CHANGE_NEW_CRUST,
    ELEV_CHANGE_RIFT,
    ELEV_CHANGE_SUBDUCTION_ARC,
    ELEV_CHANGE_TRENCH,
    ELEV_CHANGE_VOLCANO,
    ElevationLine,
    build_lines_from_lattice,
    line_spacing_rad,
    needs_regularizing,
    regularize_line,
    split_into_contiguous_runs,
)
from .noise import SphereNoise
from .plates import (
    CONTINENTAL_FRACTION,
    MIN_AUTO_PLATES,
    MAX_AUTO_PLATES,
    MIN_OCEANIC_PLATES,
    POLE_CAP_MARGIN_MULT,
    PlateWithLines,
    _INTERIOR_SUBDUCTION_MIN_RUN,
    _ROW_FULL_REVOLUTION_SLACK,
    _contested_by_any,
    _land_noise_threshold,
    _row_median_step,
)
from . import bathymetry, lithosphere, rheology, terrain_noise, torque

EXTEND_THRESHOLD_MULTIPLIER = 1.3  # same shape as v1's plates.EXTEND_THRESHOLD_RAD
MAX_EXTEND_NODES_PER_STEP = 400

# A continental line's *contested* end is allowed to retreat -- one node per step -- whether
# the overriding neighbour is oceanic (a passive margin / accretion front: the ocean slab
# descends under it and the buried continental node cedes nothing the model should keep) or
# continental (a suture whose overlapping crust is being consumed into the orogen -- see
# rheology.CONTINENTAL_COLLISION_SHORTENING_BOOST, which channels that shortening into extra
# thickening so the belt still builds real relief). Left un-retreatable, a contested end
# still grows at its *other* (divergent) side every step and never back -- the continental
# node ratchet that drives the unbounded node-count creep and the slow land-fraction decline
# (docs/TODO.md "Node-count creep") -- and, for a continent-continent pile-up, a deep
# territory overlap that just sat there for tens of Myr until the forced-merge timer fused
# the pair (the `overlapAge` view's stalled multi-plate collisions). Retreat is gated:
#   - one node per step (the existing `n_distance_cap` / `max_extend_nodes` caps already do
#     this at continental drift rates), and
#   - only where the contested node is part of a run of at least this many consecutive
#     contested nodes, so a single stray node from bounding-polygon envelope fuzz can't
#     nibble a stable coastline or, worse, sever a lobe into a spurious defragmentation plate
#     (the failure the naive "retreat every continental contested node" experiment hit -- see
#     that TODO section; the interior-subduction carve below also stays oceanic-only for the
#     same lobe-severing reason).
CONTINENTAL_CONTESTED_RETREAT_MIN_RUN = 3


def growth_seed_thickness() -> tuple[float, float]:
    """(Hc, Hm) a plate seeds *brand-new areal* nodes with -- when a line grows an end into
    open water (`_grow_or_shrink_line_for_deform`) or claims a whole new phi row
    (`_claim_adjacent_territory`).

    Always the *oceanic* reference column, regardless of the growing plate's own
    `crust_type`: any gap that opens on the sphere is floored by sea-floor spreading, not by
    the neighbouring plate's crust. Seeding a continental plate's own reference column here
    (Hc 35 km / Hm 100 km -> isostatic_elevation = +200 m) was a real land-area runaway: a
    continental plate continuously grows into the space a subducting oceanic plate vacates,
    and continental crust never subducts back, so every such step converted ocean floor into
    +200 m dry land permanently -- measured land fraction climbed 0.27 -> 0.48 and mean
    planet elevation rose ~1.7 km over 180 Myr on seed 559394024. New oceanic crust on a
    continental plate lands ~-3.5 km (a drowned passive margin / accreted terrane); genuine
    continental rifting is untouched, since that thins *existing* crust
    (`rheology.apply_divergent_deformation`) rather than growing new nodes here."""
    return lithosphere.REFERENCE_HC_OCEANIC_M, lithosphere.YOUNG_RIDGE_HM_M


def _runs_of_at_least(mask: np.ndarray, min_run: int) -> np.ndarray:
    """`mask`, with every True-run shorter than `min_run` cleared to False. Used to gate
    continental-edge retreat on a genuine multi-node contested stretch rather than a
    single stray envelope-fuzz node (see CONTINENTAL_CONTESTED_RETREAT_MIN_RUN). Runs are taken
    in the plate's concatenated node order -- a run that happens to bridge two lines' worth of
    nodes is astronomically rare (line breaks sit at a plate's theta extremes) and harmless
    if it ever happens, since `_grow_or_shrink_line_for_deform` re-checks per line anyway."""
    if min_run <= 1 or not mask.any():
        return mask.copy()
    edges = np.diff(np.concatenate(([0], mask.astype(np.int8), [0])))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1)
    long_enough = (ends - starts) >= min_run
    out = np.zeros_like(mask)
    for start, end in zip(starts[long_enough], ends[long_enough]):
        out[start:end] = True
    return out


class LithospherePlate(PlateWithLines):
    """A `PlateWithLines` whose per-node state is a lithospheric column (Hc/Hm) rather than
    an independently-set elevation -- see elevation_lines.py's own note on the two new
    OPTIONAL_FIELDS this relies on."""

    def crust_density(self) -> float:
        return lithosphere.crust_density(self.crust_type)

    # -- Motion: torque.py's real implementation -----------------------------------------

    def shift(self, world: "World", years: float) -> float:  # noqa: F821 (World only for typing)
        other_plates = [p for p in world.plates if p.plate_id != self.plate_id]
        return torque.shift_plate(self, world, other_plates, years)

    # -- Deformation: rheology.py's Mohr-Coulomb/isostasy update ---------------------------

    def deform(self, world: "World", other_plates: list, years: float, max_distance: float) -> None:  # noqa: F821
        own_points, _ = self.all_points_and_elevation()
        if not self.lines or len(own_points) == 0:
            return

        spacing_rad = line_spacing_rad(world.node_density)
        reach_rad = torque.BOUNDARY_FORCE_REACH_MULTIPLIER * spacing_rad
        extend_threshold_rad = EXTEND_THRESHOLD_MULTIPLIER * spacing_rad
        max_extend_nodes = max(1, round(MAX_EXTEND_NODES_PER_STEP * np.sqrt(world.node_density)))

        neighbours = self.get_neighbours(other_plates, threshold_rad=reach_rad)
        inputs = torque.gather_boundary_force_inputs(self, neighbours, spacing_rad, reach_rad)
        contested_all, divergent_all = torque.classify_boundary_nodes(self, neighbours, inputs, reach_rad)
        transform_all = ~contested_all & ~divergent_all & np.isfinite(inputs.dist_to_neighbor)

        neighbor_omega_all = inputs.neighbor_omega
        closing_rate_all = rheology.normal_closing_rate_m_per_s(self.omega, neighbor_omega_all, own_points, inputs.direction_to_neighbor)

        # What may retreat this step. Oceanic crust: any contested node subducts. Continental
        # crust: any contested end-node in a run of >= CONTINENTAL_CONTESTED_RETREAT_MIN_RUN
        # consecutive contested nodes -- whether the overriding neighbour is oceanic (passive
        # margin) or continental (a suture whose overlap is consumed into the orogen, with the
        # shortening channelled into extra thickening -- see the fault_factor boost below and
        # rheology.CONTINENTAL_COLLISION_SHORTENING_BOOST). Envelope fuzz (a lone contested
        # node) still can't nibble a stable margin, and the interior carve below stays
        # oceanic-only so a continental row is never severed mid-line. See
        # CONTINENTAL_CONTESTED_RETREAT_MIN_RUN for the ratchet / frozen-overlap this breaks.
        if self.crust_type != "continental":
            shrinkable_all = contested_all
        else:
            shrinkable_all = _runs_of_at_least(contested_all, CONTINENTAL_CONTESTED_RETREAT_MIN_RUN)

        years_myr = years / 1_000_000.0
        rho_c = self.crust_density()

        fault_noise = (
            SphereNoise(np.random.default_rng((world.seed, self.plate_id, 9001)), octaves=3, base_freq=9.0)
            if self.crust_type == "continental"
            else None
        )

        new_lines: list[ElevationLine] = []
        offset = 0
        for line_index, line in enumerate(self.lines):
            n = len(line)
            sl = slice(offset, offset + n)
            offset += n

            contested = contested_all[sl]
            divergent = divergent_all[sl]
            shrinkable = shrinkable_all[sl]
            closing_rate = closing_rate_all[sl]

            hc = line.crustal_thickness_m.copy()
            hm = line.mantle_lithosphere_thickness_m.copy()
            # Isostasy-driven elevation change is applied as a *delta* on top of whatever
            # elevation already holds (elevation_before -> below), not a wholesale overwrite
            # -- erosion.py (run later this same step_world call, and every step
            # thereafter until the next deform()) mutates `elevation` directly, with no
            # notion of Hc/Hm at all. An unconditional overwrite here would silently erase
            # every step's worth of erosion the instant the *next* deform() call ran,
            # confirmed directly as a real bug (a 3-step run's own elevation stopped
            # matching isostasy(Hc, Hm) exactly the way an unconditional-overwrite design
            # would have predicted, because erosion's own contribution was still baked into
            # the *un-clipped* portion of `elevation` between tectonic uplift events -- the
            # fix is this delta, not forcing elevation back to a bare isostasy readout).
            rho_c = self.crust_density()
            elevation_before = lithosphere.isostatic_elevation(hc, hm, rho_c)

            if np.any(contested):
                fault_factor = (
                    np.where(
                        fault_noise.sample(geometry.local_xyz(np.full(n, line.phi), line.theta)) < -0.15,
                        rheology.REVERSE_FAULT_VALLEY_UPLIFT_FACTOR,
                        1.0,
                    )
                    if fault_noise is not None
                    else np.ones(n)
                )
                # Continent-continent contested nodes: the overlapping crust the suture is
                # now retreating over (shrinkable, above) is thrust into the belt, not lost --
                # channel that shortening into extra plastic thickening rather than plumbing
                # the retreated column's volume through the grow/shrink pass. `fault_factor`
                # is normally <= 1 (a downthrown block accumulating less strain); here it is
                # deliberately pushed past 1 on these nodes, since it is exactly the
                # strain-accumulation multiplier apply_convergent_deformation applies and more
                # shortening is the intent. See rheology.CONTINENTAL_COLLISION_SHORTENING_BOOST.
                if self.crust_type == "continental":
                    cc_contested = contested & ~inputs.neighbor_is_oceanic[sl]
                    fault_factor = np.where(
                        cc_contested, fault_factor * rheology.CONTINENTAL_COLLISION_SHORTENING_BOOST, fault_factor
                    )
                new_hc, new_hm = rheology.apply_convergent_deformation(hc[contested], hm[contested], closing_rate[contested], years_myr, fault_factor[contested])
                hc[contested] = new_hc
                hm[contested] = new_hm

            melting = np.zeros(n, dtype=bool)
            if np.any(divergent):
                new_hc, new_hm, melt = rheology.apply_divergent_deformation(hc[divergent], hm[divergent], closing_rate[divergent], years_myr)
                hc[divergent] = new_hc
                hm[divergent] = new_hm
                melting[divergent] = melt

            prior_age = line.divergent_age_myr
            new_age = np.where(divergent, prior_age + years_myr, 0.0)
            if self.crust_type == "oceanic":
                hm = rheology.relax_young_oceanic_mantle_lithosphere(hm, new_age, years_myr)

            is_volcano = line.is_volcano.copy()
            volcano_remaining = line.volcano_active_years_remaining.copy()
            if np.any(melting):
                # Decompression melting (spec 2.3): a rift that just thinned past the
                # critical threshold erupts fresh oceanic crust in place -- same one-
                # guaranteed-eruption convention v1's stretch-volcano growth used.
                from .elevation_lines import ERUPTION_ELEVATION_M, VOLCANO_ACTIVE_MAX_YEARS, VOLCANO_ACTIVE_MIN_YEARS

                hc[melting] = lithosphere.REFERENCE_HC_OCEANIC_M
                hm[melting] = lithosphere.YOUNG_RIDGE_HM_M
                is_volcano[melting] = True
                rng = np.random.default_rng((world.seed, round(world.elapsed_years), self.plate_id, line_index))
                volcano_remaining[melting] = rng.uniform(VOLCANO_ACTIVE_MIN_YEARS, VOLCANO_ACTIVE_MAX_YEARS, size=int(melting.sum()))

            elevation_after = lithosphere.isostatic_elevation(hc, hm, rho_c)
            new_elevation = rheology.clip_elevation_bounds(line.elevation + (elevation_after - elevation_before))

            # Elevation-change provenance (diagnostic only -- see elevation_lines.ELEV_CHANGE_*
            # and render_image's "elevReason" view). Stamp whichever tectonic process moved a
            # node this step, gated on ELEV_CHANGE_MIN_DELTA_M so a node barely grazed by a
            # fading boundary force keeps its older provenance. This engine (unlike v1's
            # plates.deform) has no separate transform-uplift term, so contested convergence
            # and divergent thinning/melting are the only structural codes it emits.
            neighbor_oceanic = inputs.neighbor_is_oceanic[sl]
            reason = line.elev_change_reason.copy()
            moved = np.abs(new_elevation - line.elevation) >= ELEV_CHANGE_MIN_DELTA_M
            if self.crust_type == "continental":
                reason[contested & moved & ~neighbor_oceanic] = ELEV_CHANGE_COLLISION
                reason[contested & moved & neighbor_oceanic] = ELEV_CHANGE_SUBDUCTION_ARC
            else:
                reason[contested & moved] = ELEV_CHANGE_TRENCH
            reason[divergent & moved] = ELEV_CHANGE_RIFT
            reason[melting] = ELEV_CHANGE_VOLCANO

            updated_line = line.replace(
                elevation=new_elevation,
                crustal_thickness_m=hc,
                mantle_lithosphere_thickness_m=hm,
                divergent_age_myr=new_age,
                is_volcano=is_volcano,
                volcano_active_years_remaining=volcano_remaining,
                elev_change_reason=reason,
            )
            grown_lines = self._grow_or_shrink_line_for_deform(
                updated_line,
                inputs.dist_to_neighbor[sl],
                contested,
                shrinkable,
                spacing_rad,
                extend_threshold_rad,
                max_extend_nodes,
                max_distance,
                world,
                line_index,
                neighbours,
            )
            new_lines.extend(gl for gl in grown_lines if len(gl) > 0)

        self.set_lines(new_lines)
        self._claim_adjacent_territory(world, neighbours, spacing_rad)

        for line_index, line in enumerate(self.lines):
            if needs_regularizing(line, spacing_rad):
                self.replace_line(line_index, regularize_line(line, spacing_rad))

    def _count_open_prefix(self, theta_candidates: np.ndarray, phi: float, neighbours: list) -> int:
        if len(theta_candidates) == 0 or not neighbours:
            return len(theta_candidates)
        world_pts = geometry.to_world(self.frame, geometry.local_xyz(np.full_like(theta_candidates, phi), theta_candidates))
        contested = _contested_by_any(world_pts, neighbours)
        first_contested = np.argmax(contested) if np.any(contested) else len(contested)
        return int(first_contested)

    def _grow_or_shrink_line_for_deform(
        self,
        line: ElevationLine,
        dist: np.ndarray,
        contested: np.ndarray,
        shrinkable: np.ndarray,
        spacing_rad: float,
        extend_threshold_rad: float,
        max_extend_nodes: int,
        max_distance: float,
        world: "World",  # noqa: F821
        line_index: int,
        neighbours: list,
    ) -> list[ElevationLine]:
        """Same grow/shrink shape as `PlateWithLines._grow_or_shrink_line_for_deform` (see
        that method's own docstring -- end-only growth/shrink, plus the oceanic-only
        interior-subduction carve-out that can return a row as two contiguous
        `ElevationLine`s) -- reimplemented rather than inherited only because `grow_end`
        below needs to seed fresh Hc/Hm columns instead of a flat elevation target.
        Shrinking (end and interior) is generic over every `ElevationLine.OPTIONAL_FIELDS`
        name already (Hc/Hm included, since they're threaded through `OPTIONAL_FIELDS` --
        see elevation_lines.py), so only growth needed a new Hc/Hm-aware body."""
        theta = line.theta.copy()
        elevation = line.elevation.copy()
        contested = contested.copy()
        shrinkable = shrinkable.copy()
        dist = dist.copy()
        persistent_fields = {name: getattr(line, name).copy() for name in ElevationLine.OPTIONAL_FIELDS}
        if len(theta) == 0:
            return [ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)]

        dtheta = spacing_rad / max(np.cos(line.phi), 1e-3)
        n_distance_cap = max(1, int(max_distance / spacing_rad))

        # A row is a circle of local latitude -- its theta extent physically cannot exceed a
        # full revolution. Nothing here treats theta as periodic, so once end-growth has
        # closed the loop the "gap to nearest neighbour is wide open" test stays true forever
        # near a plate's own local pole (the pole cap belongs to nobody) and the row just
        # keeps winding. `ring_room()` is how many more `dtheta` nodes an end can take before
        # the row spans 2*pi; growth is capped by it, and at zero the end stops. Ported from
        # `PlateWithLines._grow_or_shrink_line_for_deform` (the v1 pole-winding fix) -- this
        # v2 override predates that fix and, without this, relied entirely on
        # `regularize_line`'s after-the-fact unwind, so rows still over-wound by up to a
        # revolution every step and were unwound the next (continuous churn, and near-pole
        # rings feeding overlap / node count).
        full_revolution_span = 2.0 * np.pi - _ROW_FULL_REVOLUTION_SLACK * dtheta

        def ring_room() -> int:
            if len(theta) < 2:
                return n_distance_cap
            return int(np.floor((full_revolution_span - (theta[-1] - theta[0])) / dtheta))

        def contested_run_from_end(mask: np.ndarray, from_high: bool) -> int:
            ordered = mask[::-1] if from_high else mask
            run = 0
            for value in ordered:
                if not value:
                    break
                run += 1
            return run

        if len(shrinkable) > 0 and shrinkable[-1]:
            n_remove = min(contested_run_from_end(shrinkable, from_high=True), n_distance_cap, max_extend_nodes, len(theta) - 1)
            if n_remove > 0:
                theta, elevation = theta[:-n_remove], elevation[:-n_remove]
                contested, shrinkable, dist = contested[:-n_remove], shrinkable[:-n_remove], dist[:-n_remove]
                persistent_fields = {name: values[:-n_remove] for name, values in persistent_fields.items()}

        if len(theta) == 0:
            return [ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)]

        if shrinkable[0]:
            n_remove = min(contested_run_from_end(shrinkable, from_high=False), n_distance_cap, max_extend_nodes, len(theta) - 1)
            if n_remove > 0:
                theta, elevation = theta[n_remove:], elevation[n_remove:]
                contested, shrinkable, dist = contested[n_remove:], shrinkable[n_remove:], dist[n_remove:]
                persistent_fields = {name: values[n_remove:] for name, values in persistent_fields.items()}

        if len(theta) == 0:
            return [ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)]

        # Interior subduction: carve out each substantial mid-row `shrinkable` run the two
        # end-shrinks can't reach, splitting the row into separate arcs -- see
        # `PlateWithLines._grow_or_shrink_line_for_deform` for the full rationale. Oceanic
        # self-plate only: a continental row's oceanic-contested nodes (now `shrinkable`, see
        # deform) must only ever retreat from the ends -- carving a continental row's middle
        # would sever the landmass into a spurious defragmentation plate.
        if self.crust_type == "oceanic" and len(shrinkable) >= _INTERIOR_SUBDUCTION_MIN_RUN + 2 and shrinkable[1:-1].any():
            prev_shrink = np.concatenate([[False], shrinkable[:-1]])
            run_starts = np.nonzero(shrinkable & ~prev_shrink)[0]
            keep = np.ones(len(theta), dtype=bool)
            budget = max_extend_nodes
            for start in run_starts:
                end = start
                while end < len(shrinkable) and shrinkable[end]:
                    end += 1
                if start == 0 or end >= len(shrinkable):
                    continue
                if end - start < _INTERIOR_SUBDUCTION_MIN_RUN or budget <= 0:
                    continue
                take = min(end - start, budget)
                keep[start : start + take] = False
                budget -= take
            if not keep.all():
                theta, elevation = theta[keep], elevation[keep]
                contested, shrinkable, dist = contested[keep], shrinkable[keep], dist[keep]
                persistent_fields = {name: values[keep] for name, values in persistent_fields.items()}

        # Brand-new areal crust at a growing end is oceanic regardless of this plate's own
        # type -- see growth_seed_thickness() for the land-area runaway this fixes.
        hc0, hm0 = growth_seed_thickness()
        rho_c = self.crust_density()
        new_node_elevation = float(lithosphere.isostatic_elevation(np.array([hc0]), np.array([hm0]), rho_c)[0])

        def grow_end(n_new: int) -> tuple[np.ndarray, np.ndarray]:
            return np.full(n_new, hc0), np.full(n_new, hm0)

        if not contested[-1] and dist[-1] > extend_threshold_rad and ring_room() > 0:
            gap_estimate = min(dist[-1], (n_distance_cap + 1) * spacing_rad)
            n_candidates = min(max(int(gap_estimate / spacing_rad), 1), n_distance_cap, max_extend_nodes, ring_room())
            candidate_theta = theta[-1] + dtheta * np.arange(1, n_candidates + 1)
            n_new = self._count_open_prefix(candidate_theta, line.phi, neighbours)
            if n_new > 0:
                new_theta = candidate_theta[:n_new]
                new_hc, new_hm = grow_end(n_new)
                theta = np.append(theta, new_theta)
                elevation = np.append(elevation, np.full(n_new, new_node_elevation))
                for name, values in persistent_fields.items():
                    if name == "crustal_thickness_m":
                        fill = new_hc
                    elif name == "mantle_lithosphere_thickness_m":
                        fill = new_hm
                    elif name == "elev_change_reason":
                        fill = np.full(n_new, ELEV_CHANGE_NEW_CRUST, dtype=values.dtype)
                    else:
                        fill = np.zeros(n_new, dtype=values.dtype)
                    persistent_fields[name] = np.append(values, fill)

        if not contested[0] and dist[0] > extend_threshold_rad and ring_room() > 0:
            gap_estimate = min(dist[0], (n_distance_cap + 1) * spacing_rad)
            n_candidates = min(max(int(gap_estimate / spacing_rad), 1), n_distance_cap, max_extend_nodes, ring_room())
            candidate_theta = theta[0] - dtheta * np.arange(1, n_candidates + 1)
            n_new = self._count_open_prefix(candidate_theta, line.phi, neighbours)
            if n_new > 0:
                new_theta = candidate_theta[:n_new][::-1]
                new_hc, new_hm = grow_end(n_new)
                theta = np.insert(theta, 0, new_theta)
                elevation = np.insert(elevation, 0, np.full(n_new, new_node_elevation))
                for name, values in persistent_fields.items():
                    if name == "crustal_thickness_m":
                        fill = new_hc
                    elif name == "mantle_lithosphere_thickness_m":
                        fill = new_hm
                    elif name == "elev_change_reason":
                        fill = np.full(n_new, ELEV_CHANGE_NEW_CRUST, dtype=values.dtype)
                    else:
                        fill = np.zeros(n_new, dtype=values.dtype)
                    persistent_fields[name] = np.insert(values, 0, fill)

        result = ElevationLine(phi=line.phi, theta=theta, elevation=elevation, **persistent_fields)
        return split_into_contiguous_runs(result, dtheta)

    def _claim_adjacent_territory(self, world: "World", neighbours: list, spacing_rad: float) -> None:  # noqa: F821
        """Same shape as `PlateWithLines._claim_adjacent_territory` -- a brand-new phi row
        just past this plate's own phi extremes, where open -- seeded with fresh Hc/Hm
        (oceanic reference, see `growth_seed_thickness`) plus `terrain_noise.FractalTexture`
        on Hc (an extension of an already-shaped plate, so texture rather than a fresh
        orogen), rather than a flat elevation baseline. Keyed off `(world.seed, plate_id,
        _TERRAIN_SEED_TAG)` so the texture stays attached to this plate as it grows."""
        lines_with_nodes = [line for line in self.lines if len(line) > 0]
        if not lines_with_nodes:
            return
        ordered = sorted(lines_with_nodes, key=lambda line: line.phi)
        # Keep POLE_CAP_MARGIN_MULT target spacings clear of the local pole -- see the v1
        # POLE_CAP_MARGIN_MULT comment. Right at +-pi/2 a row's theta step (spacing / cos phi)
        # blows up and the row degenerates into a handful of sub-spacing rings that read as
        # concentric circles / holes and feed the theta-winding pathology the ring_room cap in
        # `_grow_or_shrink_line_for_deform` now guards against. This v2 override predated the
        # v1 fix and still marched a plate right onto its pole (spacing_rad / 2).
        max_phi_limit = np.pi / 2 - POLE_CAP_MARGIN_MULT * spacing_rad
        # Oceanic regardless of self.crust_type -- new sphere area is floored by sea-floor
        # spreading, not by this plate's own crust (see growth_seed_thickness()).
        hc0, hm0 = growth_seed_thickness()
        amp = hc0 * 0.1  # texture on fresh Hc, same spirit as v1's noise-on-elevation
        texture = terrain_noise.FractalTexture(
            np.random.default_rng((world.seed, self.plate_id, _TERRAIN_SEED_TAG))
        )
        new_lines: list[ElevationLine] = []

        for reference, direction in ((ordered[0], -1), (ordered[-1], 1)):
            new_phi = reference.phi + direction * spacing_rad
            if abs(new_phi) > max_phi_limit:
                continue
            dtheta = spacing_rad / max(np.cos(new_phi), 1e-3)
            span = reference.theta[-1] - reference.theta[0]
            n_cols = max(int(round(span / dtheta)) + 1, 1)
            theta_candidates = reference.theta[0] + dtheta * np.arange(n_cols)
            world_pts = geometry.to_world(self.frame, geometry.local_xyz(np.full(n_cols, new_phi), theta_candidates))

            contested = _contested_by_any(world_pts, neighbours)
            open_mask = ~contested
            if not np.any(open_mask):
                continue

            theta_open = theta_candidates[open_mask]
            n_open = int(open_mask.sum())
            hc_open = np.full(n_open, hc0) + amp * texture.sample(world_pts[open_mask])
            hm_open = np.full(n_open, hm0)
            elevation_open = lithosphere.isostatic_elevation(hc_open, hm_open, self.crust_density())
            new_lines.append(
                ElevationLine(
                    phi=new_phi,
                    theta=theta_open,
                    elevation=elevation_open,
                    crustal_thickness_m=hc_open,
                    mantle_lithosphere_thickness_m=hm_open,
                    elev_change_reason=np.full(n_open, ELEV_CHANGE_NEW_CRUST, dtype=float),
                )
            )

        if new_lines:
            self.set_lines(list(self.lines) + new_lines)

    # -- Merge/split: carry Hc/Hm through, not just elevation -------------------------------

    def merge_with(self, other: "LithospherePlate", spacing_rad: float, coverage_radius_rad: float, other_points_xyz: np.ndarray) -> None:
        own_points, _ = self.all_points_and_elevation()
        other_points, _ = other.all_points_and_elevation()
        inertia_self = lithosphere.moment_of_inertia_tensor(
            own_points, self.collect("crustal_thickness_m"), self.collect("mantle_lithosphere_thickness_m"), self.crust_density(), spacing_rad
        )
        inertia_other = lithosphere.moment_of_inertia_tensor(
            other_points, other.collect("crustal_thickness_m"), other.collect("mantle_lithosphere_thickness_m"), other.crust_density(), spacing_rad
        )
        self._merge_nodes_with(other, spacing_rad, coverage_radius_rad, other_points_xyz)
        self.set_omega(torque.merge_omega(self, inertia_self, other, inertia_other))
        self.reset_age()

    def _merge_nodes_with(self, other: "LithospherePlate", spacing_rad: float, coverage_radius_rad: float, other_points_xyz: np.ndarray) -> None:
        keep_pts, _ = self.all_points_and_elevation()
        absorb_pts, _ = other.all_points_and_elevation()
        old_points = np.concatenate([keep_pts, absorb_pts], axis=0)
        old_hc = np.concatenate([self.collect("crustal_thickness_m"), other.collect("crustal_thickness_m")])
        old_hm = np.concatenate([self.collect("mantle_lithosphere_thickness_m"), other.collect("mantle_lithosphere_thickness_m")])
        exclude_tree = cKDTree(other_points_xyz) if len(other_points_xyz) else None
        self.set_lines(_lines_from_resample(self.frame, old_points, old_hc, old_hm, coverage_radius_rad, spacing_rad, exclude_tree))
        lithosphere.sync_plate_elevation(self)

    def grow_into(self, new_points_xyz: np.ndarray, new_elevation: np.ndarray, coverage_radius_rad: float, spacing_rad: float) -> None:
        old_points, _ = self.all_points_and_elevation()
        old_hc = self.collect("crustal_thickness_m")
        old_hm = self.collect("mantle_lithosphere_thickness_m")
        hc0, hm0 = lithosphere.reference_thickness(self.crust_type)
        combined_points = np.concatenate([old_points, new_points_xyz], axis=0)
        combined_hc = np.concatenate([old_hc, np.full(len(new_points_xyz), hc0)])
        combined_hm = np.concatenate([old_hm, np.full(len(new_points_xyz), hm0)])
        self.set_lines(_lines_from_resample(self.frame, combined_points, combined_hc, combined_hm, coverage_radius_rad, spacing_rad))
        lithosphere.sync_plate_elevation(self)

    def split(self, new_id: int, cut_normal: np.ndarray, min_nodes: int) -> tuple["LithospherePlate", "LithospherePlate"] | None:
        lines_a: list[ElevationLine] = []
        lines_b: list[ElevationLine] = []
        for line in self.lines:
            world_pts = line.world_xyz(self.frame)
            side = np.sum(world_pts * cut_normal, axis=-1) > 0
            # One contiguous ElevationLine per arc -- see PlateWithLines.split's own docstring
            # for why a great-circle cut can otherwise strand a row as two arcs, and why
            # carrying that whole makes the two daughters' envelopes overlap.
            ref = _row_median_step(line)
            if np.any(side):
                lines_a.extend(split_into_contiguous_runs(line.masked(side), ref))
            if np.any(~side):
                lines_b.extend(split_into_contiguous_runs(line.masked(~side), ref))

        if sum(len(l) for l in lines_a) < min_nodes or sum(len(l) for l in lines_b) < min_nodes:
            return None

        plate_a = LithospherePlate(plate_id=self.plate_id, frame=self.frame.copy(), crust_type=self.crust_type, lines=lines_a)
        plate_b = LithospherePlate(plate_id=new_id, frame=self.frame.copy(), crust_type=self.crust_type, lines=lines_b)
        return plate_a, plate_b


def _lines_from_resample(
    frame: np.ndarray,
    points: np.ndarray,
    hc: np.ndarray,
    hm: np.ndarray,
    coverage_radius_rad: float,
    spacing_rad: float,
    exclude_tree: cKDTree | None = None,
) -> list[ElevationLine]:
    """`plates._lines_from_resample`'s own algorithm, carrying Hc/Hm (nearest-point lookup)
    instead of a bare scalar elevation -- see that function's docstring for the exclusivity
    logic. `elevation` on the returned lines is a placeholder (zeros); the caller must run
    `lithosphere.sync_plate_elevation` right after to derive the real isostatic value."""
    tree = cKDTree(points)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        own_dist, _ = tree.query(world_pts)
        if exclude_tree is None:
            return own_dist < coverage_radius_rad
        other_dist, _ = exclude_tree.query(world_pts)
        return (own_dist < coverage_radius_rad) & (own_dist < other_dist)

    lines: list[ElevationLine] = []
    from .elevation_lines import iter_local_lattice

    for phi, theta_candidates, world_pts in iter_local_lattice(frame, spacing_rad=spacing_rad):
        owned = is_owned(world_pts)
        if not np.any(owned):
            continue
        theta_owned = theta_candidates[owned]
        _, idx = tree.query(world_pts[owned])
        lines.append(
            ElevationLine(
                phi=phi,
                theta=theta_owned,
                elevation=np.zeros(len(theta_owned)),
                crustal_thickness_m=hc[idx],
                mantle_lithosphere_thickness_m=hm[idx],
            )
        )
    return lines


# Hc-noise amplitude that reproduces v1's own CONTINENTAL/OCEANIC_NOISE_AMPLITUDE_M elevation
# swing through the isostasy formula rather than through elevation directly: amplitude_Hc =
# amplitude_elevation / (dz/dHc) at each crust type's own reference thickness. Continental
# dz/dHc = 1 - rho_c/rho_a =~ 0.169 (dry); oceanic uses the water-loaded factor (rho_a >
# rho_w always applies at the oceanic reference depth) =~ 0.156 -- see lithosphere.py's own
# isostatic_elevation for both branches.
_HC_NOISE_AMPLITUDE_CONTINENTAL_M = 2000.0 / (1.0 - lithosphere.RHO_CONTINENTAL_CRUST / lithosphere.RHO_ASTHENOSPHERE)
_HC_NOISE_AMPLITUDE_OCEANIC_M = 900.0 / (
    (1.0 - lithosphere.RHO_OCEANIC_CRUST / lithosphere.RHO_ASTHENOSPHERE) * (lithosphere.RHO_ASTHENOSPHERE / (lithosphere.RHO_ASTHENOSPHERE - lithosphere.RHO_WATER))
)

# Amplitudes for terrain_noise.ContinentalRelief.uplift() (orogenic belts + plateaus).
# Expressed in metres of the elevation swing each contribution should produce, then divided
# by the same 2000 m the CONTINENTAL amplitude above bakes in -- so multiplying the relief
# field (in those units) by `_HC_NOISE_AMPLITUDE_CONTINENTAL_M` lands the contribution at
# the intended elevation through the isostasy formula, and continental crust keeps one
# single Hc-noise amplitude. `_OROGENIC_RELIEF_M` is a belt crest's lift over its own
# `sample()` baseline; a between-ridge basin gets ~none of it, so that is also the depth of
# the intermontane valley below the crests. Sized (with plateaus, and the occasional
# overlap of the two) to land routine peaks near 6 km and the tallest near 7-8 km, under
# MAX_ELEVATION_M = 9000 with only occasional clipping.
_OROGENIC_RELIEF_M = 5200.0
_PLATEAU_BASE_UPLIFT_M = 2800.0
_PLATEAU_INTERNAL_RELIEF_M = 900.0
_OROGENIC_RELIEF_UNITS = _OROGENIC_RELIEF_M / 2000.0
_PLATEAU_UPLIFT_UNITS = _PLATEAU_BASE_UPLIFT_M / 2000.0
_PLATEAU_RELIEF_UNITS = _PLATEAU_INTERNAL_RELIEF_M / 2000.0

# Land gate for the uplift term: it ramps from 0 to full over `sample()` values from
# `land_threshold + _UPLIFT_SEA_MARGIN` to `+ _UPLIFT_SEA_MARGIN + _UPLIFT_SEA_RAMP` (units
# of `sample()`, ~2000 m each). So orogeny/plateaus only touch crust already well above sea
# level -- they can never lift a marine node into land or (being non-negative anyway) drop a
# coastal node into the sea, keeping the land set identical to `sample()` alone.
_UPLIFT_SEA_MARGIN = 0.05
_UPLIFT_SEA_RAMP = 0.35

# RNG tag distinguishing the terrain-noise stream from every other `(seed, plate_id, ...)`
# stream a plate draws (fault noise uses 9001, etc). numpy's SeedSequence only accepts
# integers, so this is an int, not the string "terrain".
_TERRAIN_SEED_TAG = 0x7E44A1


def _continental_sealevel_noise_offset() -> float:
    """`(Hc_at_sealevel - Hc0) / amplitude` for the continental column -- the amount by
    which a node's `relief.sample()` value can sit *below* `land_threshold` and still be
    land, because the reference continental column already floats ~+200 m above sea level.
    Passed to `_land_noise_threshold` so its quantile lands on the true land/sea crossing
    (see that function). Small and negative (~-0.10)."""
    hc0, hm0 = lithosphere.reference_thickness("continental")
    hc_sealevel = float(
        lithosphere.crustal_thickness_for_submerged_elevation(
            np.array([0.0]), np.array([float(hm0)]), lithosphere.RHO_CONTINENTAL_CRUST
        )[0]
    )
    return (hc_sealevel - hc0) / _HC_NOISE_AMPLITUDE_CONTINENTAL_M


# Each plate is seeded with one "primary" site plus this many "extra" sites, and the plate's
# territory is the *union* of its own sites' Voronoi cells rather than a single cell. Merging
# a handful of adjacent cells per plate is what turns the old one-cell-per-plate tiling (every
# plate a convex-ish blob) into lumpier, more continent-like outlines. 0 recovers the exact
# old behaviour. Kept modest on purpose: the more cells a plate fuses, the more concave its
# outline can get, and `PlateWithLines`' per-row outline is only an envelope for a genuinely
# non-convex shape (see `PlateWithLines.outline_world`).
EXTRA_SITES_PER_PLATE = 2


@dataclass
class PlateTiling:
    """A merged-Voronoi partition of the sphere into `num_plates` plates. `site_xyz` is every
    Voronoi site (unit vectors); `site_plate[s]` is which plate owns site `s`'s cell. The
    first `num_plates` sites are the per-plate "primary" sites (`site_plate[:num_plates] ==
    arange(num_plates)`), the rest are extras merged into an existing plate. A world point
    belongs to whichever plate owns its nearest site."""

    site_xyz: np.ndarray
    site_plate: np.ndarray
    num_plates: int

    def primary_site(self, plate_id: int) -> np.ndarray:
        return self.site_xyz[plate_id]


def build_plate_tiling(rng: np.random.Generator, num_plates: int, extra_sites_per_plate: int = EXTRA_SITES_PER_PLATE) -> PlateTiling:
    """Place `num_plates` primary sites plus `num_plates * extra_sites_per_plate` extra sites
    uniformly on the sphere, then hand every extra site to a plate by region-growing: each
    round, the still-unassigned site closest (angularly) to any already-assigned site joins
    that site's plate. Prim-style growth keeps each plate's set of sites a compact cluster,
    so the union of their Voronoi cells stays a single lumpy blob rather than scattering
    disconnected islands across the sphere. Deterministic in `rng`."""
    num_extra = max(0, num_plates * extra_sites_per_plate)
    site_xyz = rng.normal(size=(num_plates + num_extra, 3))
    site_xyz /= np.linalg.norm(site_xyz, axis=-1, keepdims=True)

    site_plate = np.full(len(site_xyz), -1, dtype=int)
    site_plate[:num_plates] = np.arange(num_plates)

    if num_extra > 0:
        angular = np.arccos(np.clip(site_xyz @ site_xyz.T, -1.0, 1.0))
        for _ in range(num_extra):
            unassigned = np.flatnonzero(site_plate < 0)
            assigned = np.flatnonzero(site_plate >= 0)
            block = angular[np.ix_(unassigned, assigned)]
            u, a = np.unravel_index(int(np.argmin(block)), block.shape)
            site_plate[unassigned[u]] = site_plate[assigned[a]]

    return PlateTiling(site_xyz=site_xyz, site_plate=site_plate, num_plates=num_plates)


def generate_plates(
    seed: int,
    num_plates: int | None = None,
    continental_fraction: float | None = None,
    land_fraction: float | None = None,
    node_density: float = 1.0,
    extra_sites_per_plate: int = EXTRA_SITES_PER_PLATE,
) -> list[LithospherePlate]:
    """`plates.generate_plates`'s own seed-placement/Voronoi-tiling algorithm, extended so
    each plate owns the union of several adjacent Voronoi cells (see `build_plate_tiling` and
    `EXTRA_SITES_PER_PLATE`) rather than a single cell -- still deterministic per `seed`, still
    nearest-site-owns-the-node so the tiling has no gaps/overlaps by construction, just with
    lumpier, less convex plate outlines. Only the per-plate line-building step differs from
    v1: each node gets a reference Hc/Hm plus a composite relief field on Hc (see
    `terrain_noise.py` -- a low-frequency `sample()` that decides land/sea exactly as v1's
    single noise did, plus a land-gated non-negative `uplift()` carrying orogenic belts and
    plateaus; `_HC_NOISE_AMPLITUDE_*`/`_OROGENIC_*`/`_PLATEAU_*` above set the amplitudes),
    with `elevation` itself computed once via isostasy at the end."""
    rng = np.random.default_rng(seed)
    if num_plates is None:
        num_plates = int(rng.integers(MIN_AUTO_PLATES, MAX_AUTO_PLATES + 1))

    num_continents: int | None = None
    if continental_fraction is not None:
        continental_fraction = max(0.0, min(continental_fraction, 1.0))
        num_continents = round(continental_fraction * num_plates)
        num_plates = max(num_plates, num_continents + MIN_OCEANIC_PLATES)

    tiling = build_plate_tiling(rng, num_plates, extra_sites_per_plate)
    seed_xyz = tiling.site_xyz

    if num_continents is None:
        crust_types = ["continental" if rng.random() < CONTINENTAL_FRACTION else "oceanic" for _ in range(num_plates)]
    else:
        continental_indices = set(rng.choice(num_plates, size=num_continents, replace=False).tolist())
        crust_types = ["continental" if i in continental_indices else "oceanic" for i in range(num_plates)]

    # Per-site crust type (each site inherits its owning plate's) -- `_land_noise_threshold`
    # and the `is_owned` test below both index by nearest *site*, not nearest plate.
    site_crust_types = [crust_types[tiling.site_plate[s]] for s in range(len(seed_xyz))]

    owner_tree = cKDTree(seed_xyz)
    # Composite relief fields (see terrain_noise.py) -- the last consumers of `rng`, drawn in
    # a fixed order so a given seed reproduces the same terrain. `relief.sample()` stands in
    # for the old single `SphereNoise` (same std, same land/sea decision); `relief.uplift()`
    # adds the orogenic belts and plateaus, land-gated in `hc_at` below.
    relief = terrain_noise.ContinentalRelief(
        rng,
        orogenic_units=_OROGENIC_RELIEF_UNITS,
        plateau_units=_PLATEAU_UPLIFT_UNITS,
        plateau_relief_units=_PLATEAU_RELIEF_UNITS,
    )
    ocean_relief = terrain_noise.OceanicRelief(rng)

    land_threshold = None
    if land_fraction is not None:
        land_fraction = max(0.0, min(land_fraction, 1.0))
        land_threshold = _land_noise_threshold(
            owner_tree, site_crust_types, relief, land_fraction, _continental_sealevel_noise_offset()
        )

    spacing_rad = line_spacing_rad(node_density)
    plates: list[LithospherePlate] = []
    for i in range(num_plates):
        frame = geometry.plate_frame_from_seed(tiling.primary_site(i))
        crust_type = crust_types[i]
        hc0, hm0 = lithosphere.reference_thickness(crust_type)
        hc_amp = _HC_NOISE_AMPLITUDE_CONTINENTAL_M if crust_type == "continental" else _HC_NOISE_AMPLITUDE_OCEANIC_M

        def is_owned(world_pts: np.ndarray, _i: int = i) -> np.ndarray:
            _, nearest_idx = owner_tree.query(world_pts)
            return tiling.site_plate[nearest_idx] == _i

        if crust_type == "continental":
            _lt = 0.0 if land_threshold is None else land_threshold

            def hc_at(world_pts: np.ndarray, _hc0: float = hc0, _amp: float = hc_amp, _lt: float = _lt) -> np.ndarray:
                s = relief.sample(world_pts)
                gate = np.clip((s - _lt - _UPLIFT_SEA_MARGIN) / _UPLIFT_SEA_RAMP, 0.0, 1.0)
                return _hc0 + _amp * (s - _lt) + _amp * gate * relief.uplift(world_pts)
        else:

            def hc_at(world_pts: np.ndarray, _hc0: float = hc0, _amp: float = hc_amp) -> np.ndarray:
                return _hc0 + _amp * ocean_relief.sample(world_pts)

        def elevation_at(world_pts: np.ndarray) -> np.ndarray:
            return np.zeros(len(world_pts))  # placeholder; synced from Hc/Hm below

        lines = build_lines_from_lattice(frame, is_owned, elevation_at, spacing_rad=spacing_rad)
        hc_lines = []
        for line in lines:
            world_pts = line.world_xyz(frame)
            hc = np.clip(hc_at(world_pts), lithosphere.MIN_CRUSTAL_THICKNESS_M, None)
            hm = np.full(len(line), hm0)
            hc_lines.append(line.replace(crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm))

        plate = LithospherePlate(plate_id=i, frame=frame, crust_type=crust_type, lines=hc_lines)
        lithosphere.sync_plate_elevation(plate)
        plates.append(plate)

    # Each plate above was seeded purely from its own crust type: submerged continental crust
    # is a uniform bright shelf however far from land it sits, and every continent/ocean
    # boundary is a vertical cliff. Drown the offshore continental interiors and grade the
    # boundary steps into slopes -- see bathymetry.shape_initial_bathymetry.
    bathymetry.shape_initial_bathymetry(plates)
    return plates


def new_plate(plate_id: int, frame: np.ndarray, crust_type: str, spacing_rad: float, seed: int) -> LithospherePlate:
    """A brand-new `LithospherePlate` covering `frame`'s entire local lattice at
    `spacing_rad`, seeded with reference Hc/Hm plus the same composite relief field
    `generate_plates` uses (see `terrain_noise.py`) -- the v2 analogue of
    `plates.generate_plates`' own per-plate initial-line construction. Keyed off
    `(seed, plate_id, _TERRAIN_SEED_TAG)` so the crust stays attached to this plate."""
    hc0, hm0 = lithosphere.reference_thickness(crust_type)
    rng = np.random.default_rng((seed, plate_id, _TERRAIN_SEED_TAG))
    if crust_type == "continental":
        relief = terrain_noise.ContinentalRelief(
            rng,
            orogenic_units=_OROGENIC_RELIEF_UNITS,
            plateau_units=_PLATEAU_UPLIFT_UNITS,
            plateau_relief_units=_PLATEAU_RELIEF_UNITS,
        )
        amp = _HC_NOISE_AMPLITUDE_CONTINENTAL_M

        def hc_at(world_pts: np.ndarray) -> np.ndarray:
            s = relief.sample(world_pts)
            gate = np.clip((s - _UPLIFT_SEA_MARGIN) / _UPLIFT_SEA_RAMP, 0.0, 1.0)
            return hc0 + amp * s + amp * gate * relief.uplift(world_pts)
    else:
        ocean_relief = terrain_noise.OceanicRelief(rng)
        amp = _HC_NOISE_AMPLITUDE_OCEANIC_M

        def hc_at(world_pts: np.ndarray) -> np.ndarray:
            return hc0 + amp * ocean_relief.sample(world_pts)

    def is_owned(world_pts: np.ndarray) -> np.ndarray:
        return np.ones(len(world_pts), dtype=bool)

    def elevation_at(world_pts: np.ndarray) -> np.ndarray:
        return np.zeros(len(world_pts))

    lines = build_lines_from_lattice(frame, is_owned, elevation_at, spacing_rad=spacing_rad)
    hc_lines = []
    for line in lines:
        world_pts = line.world_xyz(frame)
        hc = np.clip(hc_at(world_pts), lithosphere.MIN_CRUSTAL_THICKNESS_M, None)
        hm = np.full(len(line), hm0)
        hc_lines.append(line.replace(crustal_thickness_m=hc, mantle_lithosphere_thickness_m=hm))
    plate = LithospherePlate(plate_id=plate_id, frame=frame, crust_type=crust_type, lines=hc_lines)
    lithosphere.sync_plate_elevation(plate)
    return plate
