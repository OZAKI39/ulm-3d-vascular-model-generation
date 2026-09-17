"""Configuration objects for staged vascular growth."""

from dataclasses import dataclass


VALID_MODES = {
    "empty",
    "global",
    "previous",
    "stages",
    "ids",
    "terminals",
    "query",
}


@dataclass(frozen=True)
class TInitSpec:
    """Describe the existing vessels used to initialize one growth stage."""

    mode: str
    stage_names: tuple = ()
    vessel_ids: tuple = ()
    terminal_ids: tuple = ()
    roles: tuple = ()
    branching_modes: tuple = ()
    include_ancestors: bool = True
    include_descendants: bool = True
    branch_from_all: bool = False


@dataclass(frozen=True)
class StageConfig:
    """Bind a stage name to its configuration and initial tree."""

    name: str
    cfg: object
    tinit: TInitSpec


def _tuple_value(value):
    if value is None:
        return ()
    if isinstance(value, (list, tuple)):
        return tuple(value)
    if isinstance(value, (str, int)) and type(value) is not bool:
        return (value,)
    return None


def parse_tinit_spec(raw, stage_index):
    """
    Convert supported YAML shorthand into a ``TInitSpec``.
    """
    if isinstance(raw, TInitSpec):
        return raw

    if isinstance(raw, str):
        token = raw.strip()
        if token in {"empty", "global", "previous"}:
            return TInitSpec(
                mode=token,
                branch_from_all=(token == "global"),
            )
        if token.startswith("stage:") and token.split(":", 1)[1].strip():
            return TInitSpec(
                mode="stages",
                stage_names=(token.split(":", 1)[1].strip(),),
            )

    elif isinstance(raw, (list, tuple)) and raw:
        if all(type(item) is int and item >= 0 for item in raw):
            return TInitSpec(mode="ids", vessel_ids=tuple(raw))
        if all(isinstance(item, str) and item.strip() for item in raw):
            return TInitSpec(
                mode="stages",
                stage_names=tuple(item.strip() for item in raw),
            )

    elif isinstance(raw, dict):
        allowed = {
            "mode",
            "stages",
            "vessel_ids",
            "terminal_ids",
            "roles",
            "branching_modes",
            "include_ancestors",
            "include_descendants",
            "branch_from_all",
        }
        mode = raw.get("mode")
        stage_names = _tuple_value(raw.get("stages"))
        vessel_ids = _tuple_value(raw.get("vessel_ids"))
        terminal_ids = _tuple_value(raw.get("terminal_ids"))
        roles = _tuple_value(raw.get("roles"))
        branching_modes = _tuple_value(raw.get("branching_modes"))
        include_ancestors = raw.get("include_ancestors", True)
        include_descendants = raw.get("include_descendants", True)
        branch_from_all = raw.get("branch_from_all", mode == "global")

        valid = (
            set(raw) <= allowed
            and mode in VALID_MODES
            and stage_names is not None
            and vessel_ids is not None
            and terminal_ids is not None
            and roles is not None
            and branching_modes is not None
            and all(isinstance(item, str) and item.strip() for item in stage_names)
            and all(type(item) is int and item >= 0 for item in vessel_ids)
            and all(type(item) is int and item >= 0 for item in terminal_ids)
            and all(isinstance(item, str) and item.strip() for item in roles)
            and all(
                isinstance(item, str) and item.strip()
                for item in branching_modes
            )
            and type(include_ancestors) is bool
            and type(include_descendants) is bool
            and type(branch_from_all) is bool
            and (mode != "stages" or bool(stage_names))
            and (mode != "ids" or bool(vessel_ids))
            and (mode != "terminals" or bool(terminal_ids))
            and (mode != "query" or bool(roles or branching_modes))
        )
        if valid:
            return TInitSpec(
                mode=mode,
                stage_names=tuple(item.strip() for item in stage_names),
                vessel_ids=vessel_ids,
                terminal_ids=terminal_ids,
                roles=tuple(item.strip() for item in roles),
                branching_modes=tuple(
                    item.strip() for item in branching_modes
                ),
                include_ancestors=include_ancestors,
                include_descendants=include_descendants,
                branch_from_all=branch_from_all,
            )

    raise ValueError(f"Invalid tinit configuration for stage {stage_index}.")
