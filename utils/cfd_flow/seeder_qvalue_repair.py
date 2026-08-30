"""Source-contract helpers for the research-only Seeder qVal repair.

This module does not launch Seeder and never writes a production mesh.  It
captures the configuration and source invariants that explain why the pinned
ray/triangle test fails when otherwise identical geometry is expressed in SI
metres at micrometre scale.
"""

from __future__ import annotations

import math
import re
from typing import Any, Mapping

import numpy as np


PINNED_SEEDER_SHA = "667109df6fafdcb39f4409e3f5d90f04d75cd33c"
PINNED_TREELM_SHA = "53f273dbb8e9dcbe7feeb3d9831a35f5ae3cd72c"
CURRENT_SEEDER_SHA = "667109df6fafdcb39f4409e3f5d90f04d75cd33c"
CURRENT_TREELM_SHA = "53f273dbb8e9dcbe7feeb3d9831a35f5ae3cd72c"
MUSUBI_TREELM_SHA = "9899d1376992c4fafc8a343d2b4ccef81de670d1"
ROOT_CAUSE_CATEGORY = "E_GEOMETRY_CANDIDATES_EXIST_INTERSECTION_ROUTINE_MISSES"
TRACE_KEYS = (
    "wall_bcid_count",
    "need_calc_qval_true_count",
    "sdr_qval_by_node_call_count",
    "candidate_geometry_count",
    "intersection_found_count",
    "q_raw_0_to_1_count",
    "q_minus_one_count",
    "q_greater_than_one_count",
    "unknown_boundary_halfway_count",
    "truncate_to_one_count",
    "missing_intersected_object_count",
)


def parallel_threshold_diagnostic(
    edge1: np.ndarray,
    edge2: np.ndarray,
    direction: np.ndarray,
    *,
    epsilon: float = np.finfo(np.float64).eps,
) -> dict[str, float | bool]:
    """Compare the pinned and repaired ray/triangle parallel predicates."""

    first = np.asarray(edge1, dtype=np.float64).reshape(3)
    second = np.asarray(edge2, dtype=np.float64).reshape(3)
    ray = np.asarray(direction, dtype=np.float64).reshape(3)
    normal = np.cross(first, second)
    normal_length = float(np.linalg.norm(normal))
    direction_length = float(np.linalg.norm(ray))
    scale = normal_length * direction_length
    determinant = float(np.dot(normal, ray))
    normalized = abs(determinant) / scale if scale > 0.0 else math.nan
    return {
        "determinant_b": determinant,
        "normal_length": normal_length,
        "direction_length": direction_length,
        "dimensional_scale": scale,
        "normalized_abs_dot": normalized,
        "machine_epsilon": float(epsilon),
        "pinned_absolute_predicate_parallel": bool(abs(determinant) < epsilon),
        "repaired_scale_invariant_predicate_parallel": bool(
            scale == 0.0 or normalized <= epsilon
        ),
    }


def _table_end(text: str, opening: int) -> int:
    depth = 0
    quote: str | None = None
    escaped = False
    for index in range(opening, len(text)):
        char = text[index]
        if quote is not None:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == quote:
                quote = None
            continue
        if char in "'\"":
            quote = char
        elif char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return index + 1
    raise ValueError("unterminated Lua table")


def _assignment(table: str, key: str) -> str | None:
    match = re.search(
        rf"\b{re.escape(key)}\s*=\s*([^,}}\n]+)", table, re.IGNORECASE
    )
    return match.group(1).strip() if match else None


def _unquote(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    if len(stripped) >= 2 and stripped[0] == stripped[-1] and stripped[0] in "'\"":
        return stripped[1:-1]
    return stripped


def extract_boundary_contracts(lua_text: str) -> list[dict[str, Any]]:
    """Extract boundary/geometry fields without executing untrusted Lua."""

    starts = list(re.finditer(r"\battribute\s*=\s*\{", lua_text, re.IGNORECASE))
    records: list[dict[str, Any]] = []
    for position, match in enumerate(starts):
        attr_open = lua_text.find("{", match.start())
        attr_end = _table_end(lua_text, attr_open)
        attr = lua_text[attr_open:attr_end]
        if _unquote(_assignment(attr, "kind")) != "boundary":
            continue
        limit = starts[position + 1].start() if position + 1 < len(starts) else len(lua_text)
        tail = lua_text[attr_end:limit]
        geometry_match = re.search(r"\bgeometry\s*=\s*\{", tail, re.IGNORECASE)
        if geometry_match is None:
            raise ValueError("boundary attribute has no following geometry table")
        geometry_open = attr_end + geometry_match.start() + geometry_match.group().rfind("{")
        geometry_end = _table_end(lua_text, geometry_open)
        geometry = lua_text[geometry_open:geometry_end]
        after_geometry = lua_text[geometry_end:limit]
        transform_match = re.search(
            r"\btransformation\s*=\s*\{", after_geometry, re.IGNORECASE
        )
        transformation = None
        if transform_match is not None:
            transform_open = (
                geometry_end
                + transform_match.start()
                + transform_match.group().rfind("{")
            )
            transformation = lua_text[transform_open:_table_end(lua_text, transform_open)]
        filename = re.search(
            r"\bfilename\s*=\s*([^,}\n]+)", geometry, re.IGNORECASE
        )
        records.append(
            {
                "kind": "boundary",
                "label": _unquote(_assignment(attr, "label")),
                "level": _assignment(attr, "level"),
                "calc_dist": _assignment(attr, "calc_dist"),
                "geometry_kind": _unquote(_assignment(geometry, "kind")),
                "stl_object": _unquote(filename.group(1)) if filename else None,
                "transformation": {
                    "present": transformation is not None,
                    "deformation": _assignment(transformation or "", "deformation"),
                    "translation": _assignment(transformation or "", "translation"),
                },
                "flood_related": {
                    "flood_diagonal": _assignment(attr, "flood_diagonal"),
                    "color": _assignment(attr, "color"),
                    "inverted": _assignment(attr, "inverted"),
                },
            }
        )
    return records


def calc_dist_config_contract(
    base_lua: str, pipe_lua: str, official_lua: str
) -> dict[str, Any]:
    """Compare generated wall attributes with the official obstacle example."""

    parsed = {
        "vascular_base": extract_boundary_contracts(base_lua),
        "pipe_axis_n27": extract_boundary_contracts(pipe_lua),
        "official_tutorial_channelGeneric": extract_boundary_contracts(official_lua),
    }

    def select(records: list[dict[str, Any]], *, official: bool = False) -> dict[str, Any]:
        if official:
            candidates = [record for record in records if record["calc_dist"] is not None]
        else:
            candidates = [record for record in records if record["label"] == "wall"]
        if len(candidates) != 1:
            raise ValueError(f"expected one calc_dist boundary, found {len(candidates)}")
        return candidates[0]

    selected = {
        "vascular_base": select(parsed["vascular_base"]),
        "pipe_axis_n27": select(parsed["pipe_axis_n27"]),
        "official_tutorial_channelGeneric": select(
            parsed["official_tutorial_channelGeneric"], official=True
        ),
    }
    core_pass = all(
        record["kind"] == "boundary"
        and record["geometry_kind"] == "stl"
        and record["calc_dist"] not in {None, "false", ".false."}
        and record["stl_object"] is not None
        for record in selected.values()
    )
    return {
        "status": "PASS" if core_pass else "FAIL",
        "selected": selected,
        "global_options": {
            name: {
                key: _assignment(text, key)
                for key in ("smoothbounds", "smoothlevels", "useObstacle", "qValues")
            }
            for name, text in (
                ("vascular_base", base_lua),
                ("pipe_axis_n27", pipe_lua),
                ("official_tutorial_channelGeneric", official_lua),
            )
        },
        "semantic_conclusion": (
            "Generated wall attributes propagate calc_dist to an STL boundary "
            "with the same core syntax as the official obstacle example."
            if core_pass
            else "At least one generated calc_dist/STL attribute is inconsistent."
        ),
    }


def source_contract(
    *,
    pinned_source: str,
    current_source: str,
    pinned_seeder_sha: str,
    current_seeder_sha: str,
    pinned_treelm_sha: str,
    current_treelm_sha: str,
) -> dict[str, Any]:
    """Prove the precise pinned/current source branch used by the repair."""

    required = (
        "b =  dot_product( n, dir)",
        "if (abs(b) < eps .and. abs(a) > tiny(a)) return",
        "if (abs(b) < eps) then",
    )
    tokens_present = all(token in pinned_source for token in required)
    revisions_match = (
        pinned_seeder_sha == PINNED_SEEDER_SHA
        and current_seeder_sha == CURRENT_SEEDER_SHA
        and pinned_treelm_sha == PINNED_TREELM_SHA
        and current_treelm_sha == CURRENT_TREELM_SHA
    )
    unchanged = pinned_source == current_source
    return {
        "status": "PASS" if tokens_present and revisions_match and unchanged else "FAIL",
        "revisions": {
            "pinned_seeder": pinned_seeder_sha,
            "current_seeder": current_seeder_sha,
            "pinned_treelm": pinned_treelm_sha,
            "current_treelm": current_treelm_sha,
            "musubi_treelm_comparison": MUSUBI_TREELM_SHA,
        },
        "required_tokens_present": tokens_present,
        "pinned_equals_current": unchanged,
        "upstream_fix_available": False,
        "exact_failing_source_branch": (
            "tem_line_module.fpp::intersect_RayTriangle classifies the ray as "
            "parallel when abs(dot(cross(u,v),dir)) < epsilon(1.0)."
        ),
        "root_cause_category": ROOT_CAUSE_CATEGORY,
    }


def parse_runtime_trace(text: str) -> dict[str, int]:
    """Parse aggregate-only Seeder instrumentation output."""

    result: dict[str, int] = {}
    for key in TRACE_KEYS:
        match = re.search(rf"(?m)^\s*{re.escape(key)}\s*[:=]\s*(\d+)\s*$", text)
        if match is None:
            raise ValueError(f"missing runtime trace counter: {key}")
        result[key] = int(match.group(1))
    return result


def forbidden_production_paths_modified(paths: Mapping[str, Any] | list[str]) -> list[str]:
    """Return forbidden production files present in a change list."""

    forbidden = {"cfd_flow.py", "configs/cfd_flow.yaml", "utils/cfd_flow/pipeline.py"}
    return sorted(forbidden.intersection(str(path).replace("\\", "/") for path in paths))
