"""Shared command-line selection of the glacier march for the #822/#842 probes and replays.

`add_march_args(parser)` adds `--preset` (a name in `ice_sheet.MARCH_COST_PRESETS`),
repeatable `--set field=value` overrides of one `MarchCost` field, `--substep-years`
(the base substep, `MARCH_SUBSTEP_YEARS`, for a D2 convergence control) and
`--ice-age-period-myr` (sets `World.ice_age_period_years`, which saves store as 0).
`apply_march_args(ns, ice_sheet)` installs the selection and returns a JSON-ready record of it.
"""

from __future__ import annotations

import dataclasses
import json


def add_march_args(parser, default_preset: str = "d2") -> None:
    parser.add_argument("--preset", default=default_preset, help="ice_sheet.MARCH_COST_PRESETS name")
    parser.add_argument("--set", action="append", default=[], metavar="FIELD=VALUE",
                        help="override one MarchCost field (repeatable), e.g. trend_max_substep_years=1000")
    parser.add_argument("--substep-years", type=float, help="override MARCH_SUBSTEP_YEARS")
    parser.add_argument("--ice-age-period-myr", type=float, help="set World.ice_age_period_years on the loaded world")


def _parse_override(ice_sheet, text: str) -> tuple[str, object]:
    key, _, raw = text.partition("=")
    fields = {f.name for f in dataclasses.fields(ice_sheet.MarchCost)}
    if key not in fields:
        raise SystemExit(f"unknown MarchCost field {key!r}; known: {sorted(fields)}")
    return key, json.loads(raw) if raw not in ("hold", "endpoint") else raw


def apply_march_args(ns, ice_sheet) -> dict:
    """Install the chosen march on `ice_sheet` and return what was chosen."""
    base = getattr(ice_sheet, "march_cost", None)
    if base is None:  # a tree without MarchCost (origin/main before #839): D2 only
        if ns.preset != "d2" or ns.set:
            raise SystemExit("this tree has no MarchCost; only --preset d2 runs here")
        cost_record = None
    else:
        cost = dataclasses.replace(base(ns.preset), **dict(_parse_override(ice_sheet, t) for t in ns.set))
        ice_sheet.MARCH_COST = cost
        cost_record = dataclasses.asdict(cost)
    if ns.substep_years:
        ice_sheet.MARCH_SUBSTEP_YEARS = ns.substep_years
    return {
        "preset": ns.preset,
        "overrides": list(ns.set),
        "cost": cost_record,
        "substep_years": ice_sheet.MARCH_SUBSTEP_YEARS,
        "ice_age_period_myr": ns.ice_age_period_myr,
    }


def apply_world_args(ns, world) -> None:
    if ns.ice_age_period_myr is not None:
        world.ice_age_period_years = ns.ice_age_period_myr * 1e6
